"""Generate a WarpX input deck from a KinShock2020 ``config.yaml`` (the single
source of truth) — the forward direction of REPLICATION_PLAN.md §6.0a.

``config.yaml`` holds only the PRIMARY, physically-intuitive parameters
(densities as fractions of n0, temperatures as theta = kT/m_e c^2, lengths in
d_e / d_i, speeds as fractions of c). :func:`render` maps those onto a WarpX deck
whose ``my_constants`` are written *symbolically* — ``nt = 2.5*n0``,
``slab = 2.0*di``, ``di = de*sqrt(mass_ratio)`` — so the deck stays readable and
WarpX still records the fully-resolved values in ``warpx_used_inputs``.

Because the mapping is one-way (config -> deck), the deck never has to be hand-
edited: change ``config.yaml`` and regenerate. :func:`verify` closes the loop
after a run by parsing ``warpx_used_inputs`` and confirming the numbers WarpX
actually used match the config.
"""

from __future__ import annotations

import math

from . import units

# --- WarpX physical-constant namespace (names as used in input decks) ---
CONSTS = {
    "m_e": units.ME, "q_e": units.QE, "clight": units.C,
    "epsilon0": units.EPS0, "mu0": units.MU0, "kb": 1.380649e-23,
    "pi": math.pi, "m_p": 1.67262192369e-27,
}
_FUNCS = {"sqrt": math.sqrt, "abs": abs, "exp": math.exp, "log": math.log,
          "sin": math.sin, "cos": math.cos, "tan": math.tan, "pow": pow}


# --------------------------------------------------------------------------- #
# forward direction: config -> deck
# --------------------------------------------------------------------------- #
def _num(x: float) -> str:
    """Format a number the WarpX way: integer-valued floats get a trailing dot,
    scientific notation drops the redundant ``+`` (``1e+18`` -> ``1e18``)."""
    x = float(x)
    if x == int(x) and abs(x) < 1e6:
        return f"{int(x)}."
    return repr(x).replace("e+", "e")


# Semantic boundary name -> (field_bc, particle_bc) WarpX tokens.
#   periodic   : wrap (both faces must be periodic together)
#   reflecting : foil wall at z=0 — tangential E = 0 (pec) + specular particles.
#                See REPLICATION_PLAN.md §7 (one-sided ablation): the specular
#                approximation flips only the normal v_z, not the gyro-coupled v_y,
#                so it carries a ~5% near-wall artifact (RESULTS.md 2026-07-23).
#   symmetry   : the CORRECT one-sided perpendicular-shock wall — pec fields +
#                a pi-rotation reflection about the B0 axis (flips v_z AND the
#                transverse v_perp, preserves v_parallel), the exact discrete
#                symmetry of the two-sided problem. Emits the fork-only input
#                boundary.reflect_symmetry_axis. Use in place of 'reflecting'.
#   open       : far boundary with a background field — particles leave (absorbing)
#                while fields use pec. pec is required (not silver-mueller) because
#                the projection B-field divergence cleaner, active whenever an
#                external B is set, only accepts periodic/pec/pmc/neumann; pmc and
#                neumann would zero the tangential B0, so pec is the one div-safe
#                choice that preserves B0. Harmless here: the domain is sized so
#                nothing reaches the far boundary within the run.
#   absorbing  : true EM outflow (Silver-Mueller + particles leave). NOT compatible
#                with the B-field div cleaner, so only for field.orientation != a
#                background field (e.g. the B0=0 negative control). Use 'open' when
#                a background B is present.
_BC_MAP = {
    "periodic":   ("periodic", "periodic"),
    "reflecting": ("pec", "reflecting"),
    "symmetry":   ("pec", "reflecting"),   # + boundary.reflect_symmetry_axis (see _boundaries)
    "open":       ("pec", "absorbing"),
    "absorbing":  ("absorbing_silver_mueller", "absorbing"),
}

# Boundary names that use the rotational (pi-rotation about B0) reflection.
_SYMMETRY_BCS = {"symmetry"}


def _boundaries(geo: dict, field_axis: str = "x"):
    """((field_lo, field_hi), (particle_lo, particle_hi), symmetry_axis) from
    ``geo['boundary']``.

    ``boundary`` is either a single semantic name applied to both faces (legacy,
    e.g. ``periodic``) or a ``{lo, hi}`` mapping (e.g. ``{lo: symmetry, hi:
    open}`` for a one-sided half-domain). Each name resolves via :data:`_BC_MAP`
    to its (field, particle) WarpX token pair. ``symmetry_axis`` is ``field_axis``
    (the B0 direction) if either face is a rotational-``symmetry`` boundary, else
    ``None`` — it drives the ``boundary.reflect_symmetry_axis`` input line."""
    b = geo.get("boundary", "periodic")
    lo_name, hi_name = (b.get("lo", "periodic"), b.get("hi", "periodic")) \
        if isinstance(b, dict) else (b, b)
    try:
        flo, plo = _BC_MAP[lo_name]
        fhi, phi = _BC_MAP[hi_name]
    except KeyError as e:
        raise ValueError(f"unknown boundary {e.args[0]!r}; "
                         f"expected one of {sorted(_BC_MAP)}") from None
    sym_axis = field_axis if ({lo_name, hi_name} & _SYMMETRY_BCS) else None
    return (flo, fhi), (plo, phi), sym_axis


def _n_cell(geo: dict) -> int:
    """Number of cells: a one-sided domain [0, half] has half the cells of the
    symmetric domain [-half, +half] at the same dz."""
    span = 1.0 if geo.get("layout", "symmetric") == "one_sided" else 2.0
    return int(round(span * float(geo["domain_halfwidth_de"]) / float(geo["dz_over_de"])))


def _init_theta(role: str, kind: str) -> str:
    """my_constants expression for a species' initial thermal u_std^2 (= theta)."""
    if role == "piston":
        return "theta_e_init" if kind == "electron" else "theta_i_init"
    # ambient: electrons at theta0; ions at the same temperature -> theta0/mass_ratio
    return "theta0" if kind == "electron" else "theta0/mass_ratio"


def _alias(name: str) -> str:
    """Short tag for a species name: initials of its underscore-separated tokens
    (``piston_electrons`` -> ``pe``). Used to build collision names."""
    return "".join(tok[0] for tok in name.split("_") if tok)


def _collision_names(cfg: dict) -> list[tuple[str, str, str]]:
    """[(collision_name, species_a, species_b)] for the config's collision pairs.

    Names are ``col_<a>_<b>`` from :func:`_alias`; if two species collapse onto
    the same alias the full species names are used instead, so the mapping stays
    injective (the deck is regenerated from config, so names must be stable)."""
    pairs = units.collision_pairs(cfg)
    names = list(cfg["species"])
    tags = ({n: _alias(n) for n in names}
            if len({_alias(n) for n in names}) == len(names) else {n: n for n in names})
    return [(f"col_{tags[a]}_{tags[b]}", a, b) for a, b in pairs]


def _diag_intervals(cfg: dict):
    """(plotfile_intervals, reduced_intervals, field_intervals) from config.

    ``field_intervals`` is opt-in: if the config's diagnostics block sets it, the
    deck gains a second, field-only diagnostic (no particles) at that high cadence
    for streaks / field structure; if absent it is ``None`` and only the single
    fields+particles ``diag1`` is written (unchanged legacy behaviour)."""
    d = cfg.get("diagnostics", {}) or {}
    max_step = int(cfg["numerics"]["max_step"])
    plot = int(d.get("plotfile_intervals", max(1, round(max_step / 50))))
    red = int(d.get("reduced_intervals", max(1, round(max_step / 250))))
    field = d.get("field_intervals")
    field = int(field) if field else None
    return plot, red, field


def render(cfg: dict) -> str:
    """Build a WarpX input deck (as a string) from a loaded KinShock2020 config."""
    ref, pis = cfg["reference"], cfg["plasma"]["piston"]
    amb, geo, num = cfg["plasma"]["ambient"], cfg["geometry"], cfg["numerics"]
    ops, meta = cfg["operators"], cfg.get("meta", {})
    sc = units.derive(cfg)

    one_sided = geo.get("layout", "symmetric") == "one_sided"
    n_cell = _n_cell(geo)
    ppc_piston = int(num["ppc"]["piston"])
    ppc_ambient = int(num["ppc"]["ambient"])
    bx = "B0" if cfg["field"].get("orientation", "perpendicular") == "perpendicular" else "0."
    # B0 lies along x in the perpendicular geometry -> that is the rotational-symmetry axis.
    (field_lo, field_hi), (part_lo, part_hi), sym_axis = _boundaries(geo, field_axis="x")
    plot_int, red_int, field_int = _diag_intervals(cfg)

    L = []  # lines
    a = L.append

    # header --------------------------------------------------------------
    a("# " + "=" * 74)
    a(f"# KinShock2020 / {meta.get('run_id', '?')} -- GENERATED from config.yaml by")
    a("# scripts/make_inputs.py. DO NOT EDIT THIS FILE: edit config.yaml and")
    a("# regenerate (the config is the single source of truth). Derived scales:")
    a(f"#   d_e={sc.de:.4g} m  d_i0={sc.di0:.4g} m  B0={sc.B0:.4g} T  "
      f"dt*wpe={sc.dt_wpe:.4g}")
    domain_desc = (f"0..{geo['domain_halfwidth_de']} d_e (one-sided, wall at z=0)"
                   if one_sided else f"+-{geo['domain_halfwidth_de']} d_e")
    a(f"#   M_A={sc.MA:.3g}  M_ms={sc.Mms:.3g}  n_cell={n_cell}  "
      f"dz={geo['dz_over_de']} d_e  domain={domain_desc}")
    if meta.get("description"):
        a(f"# {' '.join(str(meta['description']).split())}")
    a("# " + "=" * 74)
    a("")

    # my_constants (symbolic; order respects dependencies) ----------------
    a("# --- primaries (mirror config.yaml) ---")
    a(f"my_constants.n0           = {_num(ref['n0'])}")
    a(f"my_constants.mass_ratio   = {_num(ref['mass_ratio'])}")
    a("my_constants.Mi           = mass_ratio*m_e")
    a("my_constants.wpe          = sqrt(n0*q_e^2/(epsilon0*m_e))")
    a("my_constants.de           = clight/wpe")
    a("my_constants.di           = de*sqrt(mass_ratio)")
    a(f"my_constants.nt           = {_num(pis['density_over_n0'])}*n0")
    a(f"my_constants.namb         = {_num(amb['density_over_n0'])}*n0")
    # B0 is a primary: emitted as a literal so the deck cannot drift from the config
    # if namb changes. v_A is derived (reported in the header) and is not needed here.
    a(f"my_constants.B0           = {_num(sc.B0)}")
    a(f"my_constants.theta_e_heat = {_num(pis['theta_e_heat'])}")
    a(f"my_constants.theta0       = {_num(amb['theta_0'])}")
    a(f"my_constants.theta_e_init = {_num(pis['theta_e_init'])}")
    a(f"my_constants.theta_i_init = {_num(pis['theta_i_init'])}")
    a(f"my_constants.slab         = {_num(geo['slab_halfwidth_di'])}*di")
    a(f"my_constants.half         = {_num(geo['domain_halfwidth_de'])}*de")
    a("")

    # time / solver -------------------------------------------------------
    a("# --- time / solver ---")
    a(f"max_step      = {int(num['max_step'])}")
    a(f"warpx.cfl     = {_num(num['cfl'])}")
    a("warpx.verbose = 1")
    a("")

    # grid ----------------------------------------------------------------
    a(f"# --- grid (dz = {geo['dz_over_de']} d_e; z = ablation / propagation axis) ---")
    a(f"amr.n_cell        = {n_cell}")
    a("amr.max_level     = 0")
    a(f"geometry.dims     = {int(geo['dims'])}")
    a("geometry.prob_lo  =  0." if one_sided else "geometry.prob_lo  = -half")
    a("geometry.prob_hi  =  half")
    a("")
    a(f"boundary.field_lo    = {field_lo}")
    a(f"boundary.field_hi    = {field_hi}")
    a(f"boundary.particle_lo = {part_lo}")
    a(f"boundary.particle_hi = {part_hi}")
    if sym_axis:
        # Rotational reflection: reflecting bounces do a pi-rotation about the B0
        # axis (flip v_z + transverse v_perp, keep v_parallel) -- the exact one-sided
        # perpendicular-shock symmetry. Fork-only input (warpx-cda).
        a(f"boundary.reflect_symmetry_axis = {sym_axis}")
    a("")
    a(f"algo.particle_shape = {int(num['particle_shape'])}")
    a("")

    # background field ----------------------------------------------------
    a(f"# --- background field ({cfg['field'].get('orientation', 'perpendicular')}: "
      "B0 along x, out of the z propagation axis) ---")
    a("warpx.B_ext_grid_init_style       = parse_B_ext_grid_function")
    a(f'warpx.Bx_external_grid_function(x,y,z) = "{bx}"')
    a('warpx.By_external_grid_function(x,y,z) = "0."')
    a('warpx.Bz_external_grid_function(x,y,z) = "0."')
    a("")

    # species -------------------------------------------------------------
    species = cfg["species"]
    a(f"particles.species_names = {' '.join(species)}")
    a("")
    for name, spec in species.items():
        kind, role = spec["kind"], spec["role"]
        ppc = ppc_piston if role == "piston" else ppc_ambient
        dens = "nt*(abs(z)<slab)" if role == "piston" else "namb*(abs(z)>=slab)"
        theta = _init_theta(role, kind)
        a(f"# --- {role} {kind}s ---")
        if kind == "electron":
            a(f"{name}.species_type = electron")
        else:
            Z = int(spec.get("charge_state", 1))
            a(f"{name}.charge = {'' if Z == 1 else f'{Z}*'}q_e")
            a(f"{name}.mass   = Mi")
        a(f'{name}.injection_style = "NUniformPerCell"')
        a(f"{name}.num_particles_per_cell_each_dim = {ppc}")
        a(f"{name}.profile = parse_density_function")
        a(f'{name}.density_function(x,y,z) = "{dens}"')
        a(f"{name}.density_min = 1.")
        a(f"{name}.momentum_distribution_type = maxwellian")
        a(f'{name}.maxwellian_u_std_distribution_type = "constant"')
        for comp in ("ux_std", "uy_std", "uz_std"):
            a(f"{name}.{comp} = sqrt({theta})")
        a("")

    # heater --------------------------------------------------------------
    h = ops["heater"]
    a("# --- laser-ablation surrogate: heat the piston slab (ParticleHeater) ---")
    a(f"particle_heater.species   = {h['species']}")
    a(f"particle_heater.intervals = {int(h['intervals'])}")
    a("particle_heater.profile   = foil")
    a("particle_heater.foil.normal      = z")
    # foil.lo/hi stay [-slab, +slab] even one-sided: the foil WIDTH sets the PSC
    # heating rate (H ~ 1/width), so keeping the full symmetric width reproduces
    # the full-domain rate; the one-sided domain then clips the heated region to
    # [0, slab]. (Using [0, slab] here would halve the width and double the rate.)
    a("particle_heater.foil.lo          = -slab")
    a("particle_heater.foil.hi          =  slab")
    a("particle_heater.foil.spot_radius = 0.")
    a("particle_heater.foil.n0          = n0")
    a("particle_heater.foil.mass_ratio  = mass_ratio")
    a(f"particle_heater.{h['species']}.theta = theta_e_heat")
    a("")

    # injector ------------------------------------------------------------
    inj = ops["injector"]
    isp, nsp = inj["species"], inj["neutralizing_species"]
    ispec, nspec = species[isp], species[nsp]
    a("# --- replenish the ablated piston slab (TargetInjector) ---")
    a(f"target_injector.species              = {isp}")
    a(f"target_injector.neutralizing_species = {nsp}")
    a(f"target_injector.intervals            = {int(inj['intervals'])}")
    a(f"target_injector.tau                  = {_num(inj['tau_over_wpe_inv'])}/wpe")
    a("target_injector.lo                   = -slab")
    a("target_injector.hi                   =  slab")
    a("target_injector.density              = nt")
    a("target_injector.reference_density    = n0")
    a(f"target_injector.ppc_reference        = {ppc_piston}")
    a(f"target_injector.{isp}.u_std = sqrt({_init_theta(ispec['role'], ispec['kind'])})")
    a(f"target_injector.{nsp}.u_std = sqrt({_init_theta(nspec['role'], nspec['kind'])})")
    a("")

    # collisions ----------------------------------------------------------
    coll = cfg.get("collisions") or {}
    if coll:
        tgt = coll.get("target") or {}
        cols = _collision_names(cfg)
        a("# --- Coulomb collisions (pairwise relativistic binary, Perez 2012 -- the")
        a("#     WarpX equivalent of the paper's Takizuka-Abe operator) ---")
        a(f"# target: {tgt.get('quantity', 'coulomb_log')} = {tgt.get('value', 'physical')}")
        a(f"#   T_e,ab = {sc.Te_ab_eV:.4g} eV, n_e,ab = {sc.n0/1e6:.4g} cm^-3")
        a(f"#   -> nu_ei,ab = {sc.nu_ei_ab:.4g} 1/s  (nu_ei*dt = {sc.nu_ei_dt:.3g})")
        a(f"#   -> mfp_ei,ab = {sc.mfp_ei_ab/sc.de:.4g} d_e,ab = "
          f"{sc.mfp_ei_ab/sc.di0:.4g} d_i0,  lambda_ab = wce/nu_ei = {sc.lambda_ab:.4g}")
        a(f"#   physical (NRL) lnLambda here would be {sc.coulomb_log_nrl:.4g}")
        a(f"my_constants.coulomb_log = {_num(sc.coulomb_log)}")
        a(f"collisions.collision_names = {' '.join(n for n, _, _ in cols)}")
        ndt = int(coll.get("ndt_supercycle", 1))
        for name, sa, sb in cols:
            a(f"{name}.type       = pairwisecoulomb")
            a(f"{name}.species    = {sa} {sb}")
            a(f"{name}.CoulombLog = coulomb_log")
            if ndt != 1:
                a(f"{name}.ndt_supercycle = {ndt}")
        a("")

    # diagnostics ---------------------------------------------------------
    a("# --- diagnostics ---")
    a("# Reduced diags: operator sanity (energy balance, piston replenishment).")
    a("warpx.reduced_diags_names = EP PN")
    a("EP.type      = ParticleEnergy")
    a(f"EP.intervals = {red_int}")
    a("PN.type      = ParticleNumber")
    a(f"PN.intervals = {red_int}")
    a("")
    diags = ["diag1"] + (["diag_fields"] if field_int else [])
    nframes = int(cfg["numerics"]["max_step"]) // plot_int if plot_int else 0
    a(f"# diag1: fields + raw particles for (z,uz) phase space (~{nframes} frames).")
    a(f"diagnostics.diags_names = {' '.join(diags)}")
    a(f"diag1.intervals = {plot_int}")
    a("diag1.diag_type = Full")
    if field_int:
        # field-only (write_species = 0) at high cadence: streaks + field structure.
        # Particles are the expensive part; fields are ~1-2 MB/frame, so this is cheap.
        field_vars = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz", "rho"]
        field_vars += [f"rho_{sp}" for sp in species]      # per-species charge density
        nff = int(cfg["numerics"]["max_step"]) // field_int
        a("")
        a(f"# diag_fields: field-only, high cadence (~{nff} frames) for streaks / field")
        a("# structure. No particles (write_species = 0); full grid resolution.")
        a(f"diag_fields.intervals     = {field_int}")
        a("diag_fields.diag_type     = Full")
        a("diag_fields.write_species = 0")
        a(f"diag_fields.fields_to_plot = {' '.join(field_vars)}")
    a("")
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# reverse direction: parse a deck (for post-run verification)
# --------------------------------------------------------------------------- #
def parse_inputs(path: str) -> dict:
    """Parse a WarpX inputs / warpx_used_inputs file into {key: raw_value_string}."""
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            key, val = line.split("=", 1)
            d[key.strip()] = val.strip()
    return d


def _eval(expr: str, ns: dict) -> float:
    """Evaluate a WarpX scalar expression (translate ^ -> **) in a restricted namespace."""
    py = expr.replace("^", "**")
    return float(eval(py, {"__builtins__": {}}, {**ns, **_FUNCS}))


def resolve_constants(d: dict) -> dict:
    """Numerically resolve all my_constants.* expressions (iterative dependency pass)."""
    exprs = {k[len("my_constants."):]: v for k, v in d.items()
             if k.startswith("my_constants.")}
    resolved = dict(CONSTS)
    pending = dict(exprs)
    for _ in range(len(pending) + 2):
        if not pending:
            break
        progressed = False
        for name, expr in list(pending.items()):
            try:
                resolved[name] = _eval(expr, resolved)
                del pending[name]
                progressed = True
            except (NameError, KeyError, ValueError):
                continue
        if not progressed:
            break
    if pending:
        raise ValueError(f"could not resolve my_constants: {list(pending)}")
    return {k: resolved[k] for k in exprs}


def key_params(path: str) -> dict:
    """Resolve the deck at ``path`` to a flat dict of comparable numeric quantities.

    Used to prove two decks are physically equivalent (generated vs authored) and
    to verify ``warpx_used_inputs`` against a config, independent of formatting,
    comments, or whether a value was written as ``20.*de`` or ``2.0*di``.
    """
    d = parse_inputs(path)
    c = resolve_constants(d)
    ns = {**CONSTS, **c}
    out = {f"const:{k}": v for k, v in c.items()}
    out["max_step"] = int(float(d["max_step"]))
    out["n_cell"] = int(float(d["amr.n_cell"].split()[0]))
    out["cfl"] = float(d["warpx.cfl"])
    out["particle_shape"] = int(float(d["algo.particle_shape"]))
    out["dims"] = int(float(d.get("geometry.dims", "1")))
    out["prob_hi"] = _eval(d["geometry.prob_hi"], ns)
    out["prob_lo"] = _eval(d["geometry.prob_lo"], ns)
    for bkey in ("boundary.field_lo", "boundary.field_hi",
                 "boundary.particle_lo", "boundary.particle_hi"):
        if bkey in d:
            out[bkey] = d[bkey].lower()
    out["heater.intervals"] = int(float(d["particle_heater.intervals"]))
    out["injector.intervals"] = int(float(d["target_injector.intervals"]))
    out["injector.tau"] = _eval(d["target_injector.tau"], ns)
    for sp in d.get("particles.species_names", "").split():
        k = f"{sp}.num_particles_per_cell_each_dim"
        if k in d:
            out[f"ppc:{sp}"] = int(float(d[k]))
    cnames = d.get("collisions.collision_names", "").split()
    out["collision_names"] = " ".join(sorted(cnames))
    for cn in cnames:
        out[f"{cn}.species"] = d.get(f"{cn}.species", "").lower()
        out[f"{cn}.type"] = d.get(f"{cn}.type", "pairwisecoulomb").lower()
        out[f"{cn}.CoulombLog"] = _eval(d[f"{cn}.CoulombLog"], ns) \
            if f"{cn}.CoulombLog" in d else -1.0
        out[f"{cn}.ndt_supercycle"] = int(float(d.get(f"{cn}.ndt_supercycle", 1)))
    for diag in ("EP", "PN", "diag1", "diag_fields"):
        k = f"{diag}.intervals"
        if k in d:
            out[k] = int(float(d[k]))
    if "diag_fields.write_species" in d:
        out["diag_fields.write_species"] = int(float(d["diag_fields.write_species"]))
    return out


def verify(cfg: dict, inputs_path: str, rtol: float = 1e-6) -> list[str]:
    """Confirm the deck at ``inputs_path`` matches ``cfg`` (parse -> resolve -> diff).

    Compares against a freshly-rendered deck for the same config, so it catches any
    drift between ``config.yaml`` and what WarpX actually ran. Returns warning
    strings (empty = match).
    """
    import os
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".inputs", delete=False) as fh:
        fh.write(render(cfg))
        gen_path = fh.name
    try:
        want = key_params(gen_path)
        got = key_params(inputs_path)
    finally:
        os.unlink(gen_path)

    warns = []
    for k in sorted(set(want) | set(got)):
        # my_constants are an implementation detail: WarpX prunes unused ones from
        # warpx_used_inputs, and the same value may be written as 20.*de or 2.*di.
        # Value-check only constants present in BOTH decks; the scalar settings
        # (max_step, n_cell, tau, ppc, ...) carry the strict comparison.
        if k not in got:
            if not k.startswith("const:"):
                warns.append(f"{k}: missing from {os.path.basename(inputs_path)}")
            continue
        if k not in want:
            continue  # extra keys in the deck are allowed
        gv, wv = got[k], want[k]
        try:
            if abs(float(gv) - float(wv)) > rtol * max(abs(float(wv)), 1e-30):
                warns.append(f"{k}: deck {gv:.6g} vs config {wv:.6g}")
        except (TypeError, ValueError):
            if gv != wv:
                warns.append(f"{k}: deck {gv!r} vs config {wv!r}")
    return warns
