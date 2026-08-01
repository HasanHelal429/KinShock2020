"""Render ``runs/<ID>/README.md`` — one page per run, with provenance for every number.

The point is that a reader (or a future session) can answer "where did this value come
from?" without re-deriving anything. Every row carries a source:

  * ``config.yaml:<key>``  — a PRIMARY. Editing it is the only way to change the run.
  * ``derived: <formula>`` — computed by ``units.derive`` from the primaries.
  * ``Table I``            — Schaeffer 2020's value, with the delta.

Prose between the ``<!-- prose:begin -->`` / ``<!-- prose:end -->`` markers is
hand-written and is **preserved verbatim** across regeneration; everything else is
generated. See ``scripts/make_run_readme.py``.
"""
from __future__ import annotations

import math
import os
import re

from . import units

PROSE_BEGIN = "<!-- prose:begin -->"
PROSE_END = "<!-- prose:end -->"

_PROSE_PLACEHOLDER = """_No hand-written notes yet. Anything written between the
`prose:begin` / `prose:end` markers survives regeneration — put the run's story here
(why it exists, what it showed, what to distrust)._"""

# Schaeffer 2020 Table I, in the units this repo reports. ``None`` = the paper does not
# quote it. Kept here rather than in a config so every run is scored against one copy.
TABLE_I = {
    "mass_ratio":        (100.0,   "mass ratio m_i/m_e"),
    "Cs_ab_over_c":      (0.030,   "C_s,ab / c"),
    "vp_model_over_c":   (0.104,   "piston speed v_p / c"),
    "MA":                (14.0,    "Alfven Mach v_sh/v_A"),
    "Mms":               (13.0,    "magnetosonic Mach"),
    # units.py reports beta = 2*mu0*n*T/B^2; Table I tabulates mu0*n*T/B^2 (its own
    # text says 2x, but 1150/0.2 only reproduce without the 2). Compare like with
    # like by doubling the paper value rather than silently flagging every run.
    "beta_ab":           (2300.0,  "ablation beta (Table I 1150, x2 convention)"),
    "beta_0":            (0.4,     "upstream beta (Table I 0.2, x2 convention)"),
    "di0_over_di_ab":    (11.18,   "d_i0 / d_i,ab"),
    "wci0_inv_over_tab": (33.9,    "gyroperiod in ablation times"),
    "lambda_ab":         (20.0,    "collisionality mfp/d_e,ab"),
}

# Rows that CLAUDE.md / RESULTS flag as knowingly off the paper. Keyed by report key.
# Shown only for rows that are ACTUALLY off (>5%), so a faithful run's table stays
# clean instead of carrying stale excuses.
KNOWN_DEVIATIONS = {
    "beta_ab": "n_amb 0.01 vs Table I's 0.008, and/or theta_e off 0.092",
    "wci0_inv_over_tab": "n_amb 25% high (1.118x) and/or theta_e recal (1.086x)",
    "di0_over_di_ab": "n_amb is 0.01 n0; Table I is 0.008 n0",
    "Cs_ab_over_c": "theta_e_heat recalibrated off the paper's 0.092",
    "MA": "model M_A from model.vsh_over_Csab; the by-eye settled value is in shock_fit.yaml",
}


def _fmt(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v == 0:
            return "0"
        if math.isinf(v):
            return "inf"
        if math.isnan(v):
            return "nan"
        return f"{v:.6g}"
    if isinstance(v, dict):     # boundary {lo,hi}, ppc {piston,ambient}, ...
        return ", ".join(f"{k} {_fmt(x)}" for k, x in v.items())
    return str(v)


def _delta(got, want) -> str:
    if got in (None, 0) or want in (None, 0):
        return "—"
    try:
        r = float(got) / float(want)
    except (TypeError, ValueError, ZeroDivisionError):
        return "—"
    pct = (r - 1.0) * 100.0
    mark = "ok" if abs(pct) < 5 else ("~" if abs(pct) < 20 else "**OFF**")
    return f"{r:.3f}x ({pct:+.1f}%) {mark}"


def _primaries(cfg: dict) -> list[tuple[str, str, str]]:
    """(label, value, config key) for every primary the deck is rendered from."""
    rows: list[tuple[str, str, str]] = []

    def add(label, path):
        node = cfg
        for k in path:
            if not isinstance(node, dict) or k not in node:
                return
            node = node[k]
        rows.append((label, _fmt(node), "`" + ".".join(path) + "`"))

    add("reference density n0 [m^-3]", ["reference", "n0"])
    add("mass ratio m_i/m_e", ["reference", "mass_ratio"])
    add("charge state Z", ["reference", "charge_state"])
    add("target density / n0", ["plasma", "piston", "density_over_n0"])
    add("heater theta_e,ab", ["plasma", "piston", "theta_e_heat"])
    add("piston init theta_e", ["plasma", "piston", "theta_e_init"])
    add("piston init theta_i", ["plasma", "piston", "theta_i_init"])
    add("ambient density / n0", ["plasma", "ambient", "density_over_n0"])
    add("ambient theta_0", ["plasma", "ambient", "theta_0"])
    add("field orientation", ["field", "orientation"])
    add("**B0 [T]**", ["field", "B0_tesla"])
    add("dims", ["geometry", "dims"])
    add("layout", ["geometry", "layout"])
    add("slab halfwidth [d_i]", ["geometry", "slab_halfwidth_di"])
    add("domain halfwidth [d_e]", ["geometry", "domain_halfwidth_de"])
    add("dz [d_e]", ["geometry", "dz_over_de"])
    add("boundary lo / hi", ["geometry", "boundary"])
    add("CFL", ["numerics", "cfl"])
    add("particle shape", ["numerics", "particle_shape"])
    add("max_step", ["numerics", "max_step"])
    add("ppc", ["numerics", "ppc"])
    add("heater intervals", ["operators", "heater", "intervals"])
    add("injector intervals", ["operators", "injector", "intervals"])
    add("injector tau [1/wpe]", ["operators", "injector", "tau_over_wpe_inv"])
    add("collisions target", ["collisions", "target"])
    return rows


def _derived(sc: units.Scales, rep: dict, max_step: float) -> list[tuple[str, str, str]]:
    """(label, value, formula) — the chain from primaries to everything else."""
    rows = [
        ("omega_pe [rad/s]", _fmt(sc.wpe), "sqrt(n0 q^2/(eps0 m_e))"),
        ("d_e [m]", _fmt(sc.de), "c/omega_pe"),
        ("d_i,ab [m]", _fmt(sc.di), "d_e*sqrt(m_i/m_e)"),
        ("C_s,ab / c", _fmt(rep["Cs_ab_over_c"]), "sqrt(theta_e,ab/mass_ratio)"),
        ("t_ab [s]", _fmt(sc.t_ab), "d_i,ab / C_s,ab"),
        ("n_amb [m^-3]", _fmt(sc.namb), "density_over_n0 * n0"),
        ("d_i0 [m]", _fmt(sc.di0), "c/omega_pi(n_amb)"),
        ("d_i0 / d_e", _fmt(sc.di0 / sc.de), "derived"),
        ("**omega_ci0 [rad/s]**", _fmt(sc.wci0), "**q_e*B0/m_i — B0 only, no n_amb**"),
        ("1/omega_ci0 [s]", _fmt(sc.wci0_inv), "1/omega_ci0"),
        ("**v_A / c**", _fmt(rep["vA_over_c"]), "**B0/sqrt(mu0*n_amb*m_i) — DERIVED**"),
        ("rho_i0 / d_e", _fmt(rep["rho_i0_over_de"]), "v_p/omega_ci0 / d_e"),
        ("v_p / c (model)", _fmt(rep["vp_model_over_c"]), "config model.vp_over_c"),
        ("v_sh / c (model)", _fmt(rep["vsh_model_over_c"]), "model.vsh_over_Csab * C_s,ab"),
        ("M_A", _fmt(sc.MA), "v_sh/v_A"),
        ("M_ms", _fmt(sc.Mms), "v_sh/sqrt(v_A^2+C_s0^2)"),
        ("beta_ab", _fmt(sc.beta_ab), "2*mu0*n0*T_e,ab/B0^2"),
        ("beta_0", _fmt(sc.beta_0), "2*mu0*n_amb*T_0/B0^2"),
        ("dz [m]", _fmt(sc.dz), "dz_over_de * d_e"),
        ("dt [s]", _fmt(sc.dt), "CFL-limited"),
        ("dt*omega_pe", _fmt(sc.dt_wpe), "dt * omega_pe"),
        ("n_cell", _fmt(sc.n_cell), "domain / dz (halved when one_sided)"),
        ("steps per 1/omega_ci0", _fmt(sc.steps_per_wci0), "1/(omega_ci0*dt)"),
        ("run length [1/omega_ci0]", _fmt(sc.wci0 * sc.dt * max_step),
         "max_step * dt * omega_ci0"),
        ("T_e,ab [eV]", _fmt(sc.Te_ab_eV), "theta_e,ab * m_e c^2"),
    ]
    if sc.coulomb_log is not None:
        rows += [
            ("lnLambda (used)", _fmt(sc.coulomb_log), "units.coulomb_log_for(collisions.target)"),
            ("lnLambda (physical)", _fmt(sc.coulomb_log_nrl), "NRL at (n0, T_e,ab)"),
            ("nu_ei,ab [1/s]", _fmt(sc.nu_ei_ab), "NRL electron-ion"),
            ("nu_ei*dt", _fmt(sc.nu_ei_dt), "must be << 1"),
            ("mfp_ei,ab / d_e", _fmt(rep.get("mfp_ei_ab_over_de")), "v_te/nu_ei / d_e"),
            ("lambda_ab", _fmt(sc.lambda_ab), "omega_ce,ab/nu_ei,ab = mfp/d_e,ab"),
            ("mfp_ii,amb / d_i0", _fmt(rep.get("mfp_ii_amb_over_di0")),
             "upstream ion-ion — what criterion 2 tests"),
        ]
    return rows


def _paper_rows(sc: units.Scales, rep: dict) -> list[tuple[str, str, str, str, str]]:
    got = dict(rep)
    got["di0_over_di_ab"] = sc.di0 / sc.di
    got["wci0_inv_over_tab"] = sc.wci0_inv / sc.t_ab
    rows = []
    for key, (want, label) in TABLE_I.items():
        if key not in got or got[key] is None:
            continue
        off = False
        try:
            off = abs(float(got[key]) / float(want) - 1.0) > 0.05
        except (TypeError, ValueError, ZeroDivisionError):
            off = True
        rows.append((label, _fmt(got[key]), _fmt(want), _delta(got[key], want),
                     KNOWN_DEVIATIONS.get(key, "") if off else ""))
    return rows


def _table(header: list[str], rows) -> list[str]:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def extract_prose(path: str) -> str:
    """Pull the hand-written block out of an existing README (empty if none)."""
    if not os.path.exists(path):
        return ""
    with open(path) as fh:
        text = fh.read()
    m = re.search(re.escape(PROSE_BEGIN) + r"\n?(.*?)\n?" + re.escape(PROSE_END),
                  text, re.S)
    return m.group(1).strip() if m else ""


def render(cfg: dict, run_dir: str) -> str:
    sc = units.derive(cfg)
    rep = sc.report()
    max_step = float(cfg.get("numerics", {}).get("max_step", 0) or 0)
    meta = cfg.get("meta", {}) or {}
    rid = meta.get("run_id", os.path.basename(os.path.normpath(run_dir)))
    prose = extract_prose(os.path.join(run_dir, "README.md")) or _PROSE_PLACEHOLDER

    L: list[str] = []
    L.append(f"# Run `{rid}`")
    L.append("")
    desc = " ".join(str(meta.get("description", "")).split())
    if desc:
        L.append(f"> {desc}")
        L.append("")
    bits = []
    if meta.get("tier"):
        bits.append(f"**tier** {meta['tier']}")
    if meta.get("deck"):
        bits.append(f"**deck** `{meta['deck']}`")
    if meta.get("reference"):
        bits.append(f"**paper** {meta['reference']}")
    if bits:
        L.append(" · ".join(bits))
        L.append("")

    L.append("## Notes")
    L.append("")
    L.append(PROSE_BEGIN)
    L.append(prose)
    L.append(PROSE_END)
    L.append("")

    L.append("## Primaries — `config.yaml` is the only source of truth")
    L.append("")
    L.append("Edit these; never the deck. Regenerate with "
             f"`python scripts/make_inputs.py runs/{rid}`.")
    L.append("")
    L += _table(["Quantity", "Value", "config key"], _primaries(cfg))
    L.append("")

    L.append("## Derived — computed by `units.derive`, not stored")
    L.append("")
    L += _table(["Quantity", "Value", "From"], _derived(sc, rep, max_step))
    L.append("")

    prows = _paper_rows(sc, rep)
    if prows:
        L.append("## vs Schaeffer 2020 Table I")
        L.append("")
        L += _table(["Quantity", "This run", "Table I", "ratio", "known cause if off"], prows)
        L.append("")

    L.append("## Files")
    L.append("")
    L.append("| Path | What |")
    L.append("|---|---|")
    L.append("| `config.yaml` | the primaries (tracked) |")
    if meta.get("deck"):
        L.append(f"| `{meta['deck']}` | generated deck (tracked, **never hand-edit**) |")
    L.append("| `warpx_used_inputs` | what WarpX actually ran (tracked) |")
    L.append("| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |")
    L.append("| `diags/`, `*.log` | output (gitignored, regenerable) |")
    L.append("")
    L.append("---")
    L.append("")
    L.append("_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, "
             "not the tables. Prose between the `prose:` markers is preserved._")
    L.append("")
    return "\n".join(L)
