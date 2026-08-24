#!/usr/bin/env python3
"""Schaeffer 2020 Table I side by side with the E_phase `eps` ladder.

WHAT THIS IS FOR. The ladder was *built* so that every Table I group is held fixed
and only `eps = v_te,ab/c = sqrt(theta_e,ab)` moves (see runs/E_phase/README.md).
This renders that claim as a table you can read: the paper's two printed columns,
then one column per run, every entry in the row's stated unit, with a `spread` tag
that says whether the run columns actually agree.

    [=]  max/min over the run columns < 1%   -- held fixed, as designed
    [x]  they differ                         -- this is a knob that moved

THE TWO PAPER COLUMNS ARE NOT THE SAME PLASMA.
  * `Table I PSC` is the code column: PSC's normalization, with its *reduced* speed
    of light. Sec. II sets c = sqrt(mu_p/T_e,ab) C_s,ab, so PSC's theta_e,ab = 0.092
    IS its speed of light, and every "x c" in that column is a fraction of a c that
    is 0.1 of the real one.
  * `Table I phys` is the caption's "one possible set of experimentally-relevant
    physical values" -- a real *hydrogen* (mu = 1836) plasma at T_e,ab = 470 eV. It
    is an illustration, not a unit map, so its ion rows do NOT follow from the code
    column by any single factor (see scripts/table1.py's header, and REPLICATION_PLAN
    Sec. 1: velocities x0.0234, temperatures /100, lengths unscaled).

WHERE THE RUNS SIT BETWEEN THEM. WarpX has no reduced-c option, so a rung's
theta_e,ab is a *physical* electron temperature at the real c. Rung 5 (`es_47keV`)
picks theta_e,ab = 0.092 and therefore reproduces the PSC code column numerically --
same dimensionless problem, real c, 47 keV electrons. Every lower rung is the same
dimensionless problem at a smaller eps, i.e. a colder, more magnetized-relative-to-
Debye electron population, which is the one axis this phase is scanning.

    python scripts/table1_ladder.py                     # the four E_phase rungs
    python scripts/table1_ladder.py --all               # + the 470 eV anchor and ek_47keV_t2
    python scripts/table1_ladder.py --markdown          # GitHub table
    python scripts/table1_ladder.py runs/A runs/B ...   # any runs
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import yaml                                                    # noqa: E402
from kinshock import units as U                                # noqa: E402

# The four ladder rungs live in E_phase; rung 1 is the 470 eV anchor in S_phase and
# `ek_47keV_t2` is the later 47 keV precursor run, both --all extras.
LADDER = ["runs/E_phase/es_1p5keV", "runs/E_phase/es_4p7keV",
          "runs/E_phase/es_15keV", "runs/E_phase/es_47keV"]
EXTRAS = ["runs/S_phase/ss_dz16_ppc100", "runs/E_phase/ek_47keV_t2"]


def scales(run_dir: Path):
    """Every number this table needs, from the config primaries only."""
    cfg = yaml.safe_load((run_dir / "config.yaml").read_text())
    sc = U.derive(cfg)
    mu, mi, c, qe = sc.mass_ratio, sc.mi, U.C, U.QE

    th_ab = float(cfg["plasma"]["piston"]["theta_e_heat"])
    th_0 = float(cfg["plasma"]["ambient"]["theta_0"])
    vte_ab = math.sqrt(th_ab) * c
    vte_0 = math.sqrt(th_0) * c

    lD_ab = vte_ab / sc.wpe
    wpe0 = math.sqrt(sc.namb * qe * qe / (U.EPS0 * U.ME))
    lD_0 = vte_0 / wpe0
    wce, wci = qe * sc.B0 / U.ME, qe * sc.B0 / mi
    ppc = float(cfg["numerics"]["ppc"]["ambient"])
    steps = float(cfg["numerics"]["max_step"])

    d = dict(sc.__dict__)
    d.update(
        cfg=cfg, sc=sc, theta_e_ab=th_ab, theta_0=th_0, eps=math.sqrt(th_ab),
        vte_ab=vte_ab, lD_ab=lD_ab, lD_0=lD_0, wpe0=wpe0, wce=wce, wci=wci,
        Te_ab_eV=th_ab * U.ME_C2_EV, T_0_eV=th_0 * U.ME_C2_EV,
        n_amb_frac=sc.namb / sc.n0, n_t_frac=sc.nt / sc.n0,
        Cs_ab_over_c=sc.Cs_ab / c, vA_over_c=sc.vA / c, vp_over_c=sc.vp_model / c,
        vsh_over_c=sc.vsh_model / c, c_over_Csab=c / sc.Cs_ab,
        vsh_over_Csab=sc.vsh_model / sc.Cs_ab, vp_over_Csab=sc.vp_model / sc.Cs_ab,
        vp_over_vA=sc.vp_model / sc.vA,
        di_over_de=sc.di / sc.de, di0_over_diab=sc.di0 / sc.di,
        lD_ab_over_de=lD_ab / sc.de, lD_0_over_di0=lD_0 / sc.di0,
        rho_i0_over_di0=sc.rho_i0 / sc.di0, rho_sh_over_di0=(sc.vsh_model / wci) / sc.di0,
        rho_e_ab_over_de=(vte_ab / wce) / sc.de,
        wci0_inv_over_tab=sc.wci0_inv / sc.t_ab, wpe_over_wce=sc.wpe / wce,
        wpe0_over_wce=wpe0 / wce, wpi0_over_wci=wpe0 / math.sqrt(mu) / wci,
        sigma_i=(sc.vA / c) ** 2, sigma_e=(sc.vA / c) ** 2 * mu,
        gamma_te=1.0 + 1.5 * th_ab,
        nu_ei_over_wce=(sc.nu_ei_ab / wce) if sc.nu_ei_ab else float("nan"),
        tau_ei_over_tab=(1.0 / sc.nu_ei_ab / sc.t_ab) if sc.nu_ei_ab else float("nan"),
        dz_over_de=sc.dz / sc.de, dz_over_lD_ab=sc.dz / lD_ab, dz_over_lD_0=sc.dz / lD_0,
        ppc=ppc, N_D_pic=ppc * lD_0 / sc.dz, N_D_phys=sc.namb * lD_0 ** 3,
        max_step=steps, t_run_wci0=sc.dt * steps * sc.wci0,
        cfl=float(cfg["numerics"]["cfl"]),
        particle_shape=float(cfg["numerics"]["particle_shape"]),
        filter_npass=float(cfg["numerics"].get("filter_npass", 0)),
        domain_over_di0=sc.domain_halfwidth / sc.di0,
        L_target_over_di0=2.0 * sc.di / sc.di0,
        tau_sim_over_tab=sc.dt * steps / sc.t_ab,
    )
    return d


# --------------------------------------------------------------------------- #
# the table.  (section, label, unit, key, PSC-as-printed, physical-as-printed)
# Paper strings are Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I, p.4;
# "" means the paper does not print that row and the entry is left blank on purpose.
# --------------------------------------------------------------------------- #
ROWS = [
    ("ABLATION (piston) -- the plasma the heater maintains", None, None, None, None, None),
    ("n_e,ab",           "m^-3",      "n0",            "1.25 (code)", "6e20 cm^-3"),
    ("n_t/n_e,ab",       "",          "n_t_frac",      "2.0",         ""),
    ("T_e,ab",           "eV",        "Te_ab_eV",      "0.092 m_ec^2", "470 eV"),
    ("theta_e,ab = kT/m_ec^2", "",    "theta_e_ab",    "0.092",       "9.2e-4 (real c)"),
    ("eps = v_te,ab/c = sqrt(theta)", "", "eps",       "0.3033",      "0.03033 (real c)"),
    ("v_te,ab/c",        "",          "eps",           "0.3033",      ""),
    ("C_s,ab/c",         "",          "Cs_ab_over_c",  "0.030",       "210 km/s"),
    ("d_e,ab",           "m",         "de",            "1 d_e,ab",    ""),
    ("d_i,ab/d_e,ab",    "",          "di_over_de",    "10",          ""),
    ("lambda_D,ab/d_e,ab", "",        "lD_ab_over_de", "0.3033",      ""),
    ("t_ab = d_i,ab/C_s,ab", "s",     "t_ab",          "1 t_ab",      "50 ps"),
    ("v_p/c",            "",          "vp_over_c",     "0.104",       "730 km/s"),
    ("v_p/C_s,ab",       "",          "vp_over_Csab",  "3.43",        ""),
    ("v_sh/C_s,ab",      "",          "vsh_over_Csab", "4.6",         "980 km/s"),

    ("UPSTREAM (ambient) -- what the shock runs into", None, None, None, None, None),
    ("B_0",              "T",         "B0",            "0.01 sqrt(m_ec^2)", "7 T"),
    ("n_e0/n_e,ab",      "",          "n_amb_frac",    "0.008",       ""),
    ("n_e0",             "m^-3",      "namb",          "0.01 (code)", "4.8e18 cm^-3"),
    ("T_0",              "eV",        "T_0_eV",        "0.002 m_ec^2", "10 eV"),
    ("d_i0/d_i,ab",      "",          "di0_over_diab", "11.2",        ""),
    ("d_i0",             "m",         "di0",           "",            "104 um"),
    ("v_A/c",            "",          "vA_over_c",     "0.0100",      ""),
    ("1/w_ci0",          "t_ab",      "wci0_inv_over_tab", "33.9",    "1.5 ns"),
    ("1/w_ci0",          "s",         "wci0_inv",      "",            "1.5 ns"),
    ("rho_i0/d_i0 (=v_p/v_A)", "",    "rho_i0_over_di0", "10.4",      ""),
    ("rho_sh/d_i0 (=v_sh/v_A)", "",   "rho_sh_over_di0", "14",        ""),

    ("DIMENSIONLESS -- Table I's invariants; these define the problem", None, None, None, None, None),
    ("m_i/m_e",          "",          "mass_ratio",    "100",         "1836 (H+)"),
    ("beta_ab = mu0 nkT/B^2", "",     "beta_ab",       "1150",        ""),
    ("beta_0",           "",          "beta_0",        "0.2",         ""),
    ("M_A = v_sh/v_A",   "",          "MA",            "14",          ""),
    ("M_ms",             "",          "Mms",           "13",          ""),
    ("T_0/T_e,ab",       "",          "_T_ratio",      "0.0217",      ""),

    ("THE eps FAMILY -- everything that IS the speed of light", None, None, None, None, None),
    ("c/C_s,ab  <- THE reduced c", "", "c_over_Csab",  "33.0",        "1429 (real c)"),
    ("w_pe,ab/w_ce",     "",          "wpe_over_wce",  "111.8",       ""),
    ("w_pe,0/w_ce (= rho_e/lambda_D)", "", "wpe0_over_wce", "10.0",   ""),
    ("w_pi,0/w_ci",      "",          "wpi0_over_wci", "100",         ""),
    ("lambda_D,0/d_i0",  "",          "lD_0_over_di0", "4.47e-3",     ""),
    ("rho_e,ab/d_e,ab",  "",          "rho_e_ab_over_de", "33.9",     ""),
    ("sigma_i = B^2/mu0 n_i m_ic^2", "", "sigma_i",    "1e-4",        ""),
    ("sigma_e",          "",          "sigma_e",       "1e-2",        ""),
    ("gamma(T_e,ab) = 1+1.5theta", "", "gamma_te",     "1.138",       ""),

    ("COLLISIONS", None, None, None, None, None),
    ("lambda_ab = mfp_e/d_e,ab", "",  "lambda_ab",     "20",          ""),
    ("ln Lambda (used in deck)", "",  "coulomb_log",   "n/a (Takizuka-Abe)", ""),
    ("ln Lambda (NRL physical)", "",  "coulomb_log_nrl", "",          ""),
    ("nu_ei,ab/w_ce",    "",          "nu_ei_over_wce", "",           ""),
    ("tau_ei,ab",        "t_ab",      "tau_ei_over_tab", "0.009",     "0.43 ps"),

    ("NUMERICS -- not physics; listed to separate the confound", None, None, None, None, None),
    ("dz/d_e,ab",        "",          "dz_over_de",    "0.3",         ""),
    ("dz/lambda_D,0",    "",          "dz_over_lD_0",  "",            ""),
    ("dz",               "m",         "dz",            "",            ""),
    ("dt",               "s",         "dt",            "",            ""),
    ("n_cell",           "",          "n_cell",        "30000",       ""),
    ("domain half-width", "d_i0",     "domain_over_di0", "900 d_i,ab", "8.4 mm"),
    ("L_target",         "d_i0",      "L_target_over_di0", "2 d_i,ab", ""),
    ("ppc (ambient)",    "",          "ppc",           "",            ""),
    ("CFL",              "",          "cfl",           "",            ""),
    ("particle_shape",   "",          "particle_shape", "",           ""),
    ("filter_npass",     "",          "filter_npass",  "",            ""),
    ("N_D,PIC = ppc*lambda_D0/dz", "", "N_D_pic",      "",            ""),
    ("max_step",         "",          "max_step",      "400000",      ""),
    ("t_run",            "1/w_ci0",   "t_run_wci0",    "6.5",         ""),
    ("t_run",            "t_ab",      "tau_sim_over_tab", "220",      "10.9 ns"),
]


# Rows that are counts, not measurements -- printed as integers.
INTEGRAL = {"n_cell", "max_step", "ppc", "mass_ratio",
            "particle_shape", "filter_npass"}


def value(d, key):
    if key == "_T_ratio":
        return d["theta_0"] / d["theta_e_ab"]
    v = d.get(key)
    return float("nan") if v is None else float(v)


def fmt(v, integral=False):
    if v != v:
        return "--"
    if v == 0:
        return "0"
    if integral:
        return f"{int(round(v)):,d}"
    return f"{v:.4g}"


def spread_tag(vals):
    """[=] when the run columns agree to <1%, [x] when they do not."""
    good = [v for v in vals if v == v]
    if len(good) < 2:
        return "", ""
    if min(good) == 0:
        return ("[=]", "") if max(good) == 0 else ("[x]", "0 vs nonzero")
    r = max(good) / min(good)
    return ("[=]" if r < 1.01 else "[x]"), ("" if r < 1.01 else f"x{r:.4g}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*", default=None,
                    help="run dirs (default: the four E_phase eps rungs)")
    ap.add_argument("--all", action="store_true",
                    help="also the 470 eV anchor (S_phase/ss_dz16_ppc100) and ek_47keV_t2")
    ap.add_argument("--markdown", action="store_true", help="emit a GitHub table")
    ap.add_argument("--json", action="store_true",
                    help="emit the same table as JSON, for a renderer downstream")
    ap.add_argument("--only-differing", action="store_true",
                    help="drop the rows the runs hold fixed")
    args = ap.parse_args(argv)

    paths = [Path(p) for p in (args.runs or LADDER)]
    if args.all and not args.runs:
        paths = [Path(EXTRAS[0])] + paths + [Path(EXTRAS[1])]
    paths = [p if p.is_absolute() else ROOT / p for p in paths]
    data = [(p.name, scales(p)) for p in paths]
    names = [n for n, _ in data]

    hdr = ["Parameter", "Unit", "Table I PSC", "Table I phys"] + names + ["", "spread"]
    body = []
    for row in ROWS:
        if row[1] is None:                                   # section heading
            body.append(("SECTION", row[0]))
            continue
        label, unit, key, psc, phys = row
        vals = [value(d, key) for _, d in data]
        tag, ratio = spread_tag(vals)
        if args.only_differing and tag != "[x]":
            continue
        body.append(("ROW", [label, unit, psc, phys]
                    + [fmt(v, key in INTEGRAL) for v in vals]
                     + [tag, ratio]))

    pruned, i = [], 0
    while i < len(body):
        if body[i][0] == "SECTION" and (i + 1 == len(body) or body[i + 1][0] == "SECTION"):
            i += 1
            continue
        pruned.append(body[i])
        i += 1
    body = pruned

    if args.json:
        emit_json(hdr, body, data)
    elif args.markdown:
        emit_markdown(hdr, body, data)
    else:
        emit_text(hdr, body, data)
    return 0


def emit_text(hdr, body, data):
    rows = [r for k, r in body if k == "ROW"]
    w = [max(len(hdr[i]), *(len(r[i]) for r in rows)) for i in range(len(hdr))]
    total = sum(w) + 2 * len(w)
    print()
    print(f"Schaeffer 2020 Table I  vs  {len(data)} WarpX runs "
          f"({', '.join(n for n, _ in data)})")
    print("Table I PSC = code units at PSC's REDUCED c;  Table I phys = the paper's "
          "mu=1836 hydrogen illustration.")
    print("Run columns are at the REAL c, in the row's stated unit.  "
          "[=] runs agree <1%   [x] runs differ")
    print("=" * total)
    print("  ".join(h.ljust(w[i]) if i < 4 else h.rjust(w[i])
                    for i, h in enumerate(hdr)))
    print("=" * total)
    for kind, r in body:
        if kind == "SECTION":
            print("-" * total)
            print(f"  == {r} ==")
            print("-" * total)
        else:
            print("  ".join(c.ljust(w[i]) if i < 4 else c.rjust(w[i])
                            for i, c in enumerate(r)))
    print("=" * total)
    diff = sum(1 for k, r in body if k == "ROW" and r[-2] == "[x]")
    same = sum(1 for k, r in body if k == "ROW" and r[-2] == "[=]")
    print(f"\n{same} rows held fixed across the runs, {diff} rows differ.")
    eps_powers(body, data)
    print(NOTES)


def eps_powers(body, data):
    """Is the difference one-dimensional?  Fit every differing row as eps^p.

    The ladder's whole claim is that eps is the ONLY knob that moved, so each
    differing row should be an exact power of it.  A least-squares fit of
    log(value) against log(eps) reports that power and the residual; anything with
    an rms above ~1e-3 is a row that is NOT a pure function of eps, i.e. a second
    knob hiding in the configs.
    """
    # Two rows are NOT homogeneous in eps and must not fire the alarm: gamma is
    # 1 + 1.5 eps^2 (an affine function of a power, so its log-log slope drifts) and
    # the NRL lnLambda is 24 - ln(sqrt(n)/T), i.e. logarithmic in eps. Both are
    # exempt by construction, not by fit quality.
    exempt = {"gamma(T_e,ab) = 1+1.5theta": "affine 1 + 1.5 eps^2",
              "ln Lambda (NRL physical)": "logarithmic, 24 - ln(sqrt n / T)"}
    eps = [d["eps"] for _, d in data]
    if len(eps) < 3 or max(eps) / min(eps) < 1.01:
        return
    le = [math.log(e) for e in eps]
    n = len(le)
    mx = sum(le) / n
    sxx = sum((x - mx) ** 2 for x in le)
    print("\nEvery differing row as a power of eps = v_te,ab/c   "
          "(least squares on log-log; rms is the residual)")
    print(f"  {'row':38}{'p in eps^p':>12}{'rms':>10}")
    worst = 0.0
    for kind, r in body:
        if kind != "ROW" or r[-2] != "[x]":
            continue
        try:
            ly = [math.log(abs(float(c.replace(",", "")))) for c in r[4:4 + n]]
        except ValueError:
            continue
        my = sum(ly) / n
        p = sum((x - mx) * (y - my) for x, y in zip(le, ly)) / sxx
        rms = math.sqrt(sum((y - my - p * (x - mx)) ** 2 for x, y in zip(le, ly)) / n)
        note = exempt.get(r[0], "")
        if not note:
            worst = max(worst, rms)
        print(f"  {r[0]:38}{p:>+12.4f}{rms:>10.1e}"
              + (f"   exempt: {note}" if note else ""))
    print(f"  -> worst residual {worst:.1e} (exempt rows excluded).  "
          + ("Every differing row IS a pure power of eps: the difference between "
             "these runs is ONE-DIMENSIONAL." if worst < 5e-3 else
             "At least one row is NOT a pure power of eps -- a second knob moved."))


NOTES = """
Reading the paper columns
  * `tau_ei,ab`.  Table I's 0.009 t_ab is the PHYSICAL column's number: 0.43 ps / 50 ps,
    both at mu = 1836.  At the code column's mu = 100, lambda_ab = 20 forces
    tau_ei/t_ab = lambda_ab/(mu/sqrt(mu) * d_i/d_e) = 20/100 = 0.2 exactly, which is what
    every run reports.  The two are not the same quantity; the runs match lambda_ab = 20,
    which is the row Sec. II actually defines.
  * `beta_0` 0.2 vs 0.1957.  beta_0 is DERIVED from n_e0, T_0 and a B0 already fixed by
    beta_ab, so it is not free; 0.1957 is what Table I's own primaries give and 0.2 is
    the printed rounding (scripts/table1.py header).
  * `M_A` 14 / `M_ms` 13 / `d_i0/d_i,ab` 11.2 / `rho_sh/d_i0` 14 are Table I roundings of
    13.95, 12.76, 11.18, 13.95.
  * Dimensional rows (`d_e,ab`, `d_i0`, `t_ab`, `1/w_ci0` [s], `dz`, `dt`) do NOT match
    the physical column, and cannot: the runs pin n_e,ab = 6e26 m^-3 at mu = 100, while
    the physical column is a mu = 1836 hydrogen illustration.  Only the ratios transfer.
  * `n_cell`, `domain`, `max_step`, `t_run`.  The ladder is the S_phase EARLY-FORMATION
    box -- 12.02 d_i0 to t*w_ci0 = 0.302 -- not the paper's 900 d_i,ab to 220 t_ab.  That
    is a deliberate truncation shared by all four rungs, so it cannot confound them.
  * `ln Lambda (used in deck)`.  A knob, not a physical value: it is whatever reproduces
    lambda_ab = 20 at that rung's T_e,ab, and nu_ei ~ T^-3/2 makes it climb as eps^-4.
    The NRL physical lnLambda at the same (n, T) is the row below it.
"""


def emit_json(hdr, body, data):
    import json
    out = dict(
        runs=[dict(name=n, eps=d["eps"], Te_ab_eV=d["Te_ab_eV"], B0=d["B0"],
                   wpe0_over_wce=d["wpe0_over_wce"], n_cell=int(d["n_cell"]),
                   max_step=int(d["max_step"]), dz_over_lD_0=d["dz_over_lD_0"])
              for n, d in data],
        header=hdr,
        rows=[dict(kind=k, section=r) if k == "SECTION" else
              dict(kind=k, label=r[0], unit=r[1], psc=r[2], phys=r[3],
                   values=r[4:-2], tag=r[-2], ratio=r[-1])
              for k, r in body],
        notes=NOTES.strip(),
    )
    print(json.dumps(out, indent=1))


def emit_markdown(hdr, body, data):
    cols = [h for h in hdr[:-2]] + ["spread"]
    print(f"\n### Schaeffer 2020 Table I vs the eps ladder\n")
    print("| " + " | ".join(cols) + " |")
    print("|" + "|".join(["---"] * len(cols)) + "|")
    for kind, r in body:
        if kind == "SECTION":
            print(f"| **{r}** |" + " |" * (len(cols) - 1))
        else:
            cells = r[:-2] + [f"{r[-2]} {r[-1]}".strip()]
            print("| " + " | ".join(c if c else "" for c in cells) + " |")


if __name__ == "__main__":
    raise SystemExit(main())
