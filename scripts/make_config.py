#!/usr/bin/env python3
"""Regenerate / verify a run's config.yaml from the WarpX inputs it actually ran.

WarpX writes ``warpx_used_inputs`` (the fully-resolved input deck) into each run
directory. This script parses it, resolves the ``my_constants`` expressions
numerically, maps them to the KinShock2020 config primaries, and either writes
``config.generated.yaml`` or overwrites ``config.yaml`` (``--write``). It always
diffs against the existing ``config.yaml`` so the config provably matches what was
simulated (REPLICATION_PLAN.md §6.0a).

Usage:
    python scripts/make_config.py runs/R1                 # verify + write config.generated.yaml
    python scripts/make_config.py runs/R1 --write         # overwrite config.yaml
    python scripts/make_config.py runs/R1 --inputs deck   # parse an arbitrary deck instead
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import yaml

# --- WarpX physical-constant namespace (names as used in input decks) ---
CONSTS = {
    "m_e": 9.1093837015e-31, "q_e": 1.602176634e-19, "clight": 299792458.0,
    "epsilon0": 8.8541878128e-12, "mu0": 1.25663706212e-6, "kb": 1.380649e-23,
    "pi": math.pi, "m_p": 1.67262192369e-27,
}
_FUNCS = {"sqrt": math.sqrt, "abs": abs, "exp": math.exp, "log": math.log,
          "sin": math.sin, "cos": math.cos, "tan": math.tan, "pow": pow}


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
    # keep only the user constants (drop the base physical constants)
    return {k: resolved[k] for k in exprs}


def _eval(expr: str, ns: dict) -> float:
    """Evaluate a WarpX scalar expression (translate ^ -> **) in a restricted namespace."""
    py = expr.replace("^", "**")
    return float(eval(py, {"__builtins__": {}}, {**ns, **_FUNCS}))


def _num(d: dict, key: str, consts: dict, default=None):
    """Resolve a possibly-symbolic parameter value to a float."""
    if key not in d:
        return default
    try:
        return _eval(d[key], {**CONSTS, **consts})
    except Exception:
        return default


def build_config(d: dict, consts: dict) -> dict:
    """Map parsed inputs + resolved constants to the KinShock2020 config structure."""
    C, clight = consts, CONSTS["clight"]
    me, qe, eps0 = CONSTS["m_e"], CONSTS["q_e"], CONSTS["epsilon0"]
    n0 = C["n0"]
    mass_ratio = C.get("mass_ratio", _num(d, "piston_ions.mass", C, me * 100) / me)
    de = C.get("de", clight / math.sqrt(n0 * qe ** 2 / (eps0 * me)))
    di = C.get("di", de * math.sqrt(mass_ratio))

    prob_lo = _num(d, "geometry.prob_lo", C)
    prob_hi = _num(d, "geometry.prob_hi", C)
    n_cell = int(float(d.get("amr.n_cell", "0").split()[0]))
    dz = (prob_hi - prob_lo) / n_cell if n_cell else float("nan")
    slab = C.get("slab", _num(d, "particle_heater.foil.hi", C))

    theta_e_heat = _num(d, "particle_heater.piston_electrons.theta", C)
    theta0 = C.get("theta0", _num(d, "amb_electrons.ux_std", C, 0) ** 2)

    cfg = {
        "meta": {
            "run_id": os.path.basename(os.path.dirname(os.path.abspath(d.get("__path__", ".")))),
            "generated_from": "scripts/make_config.py <- warpx_used_inputs",
        },
        "reference": {"n0": n0, "mass_ratio": mass_ratio, "charge_state": 1},
        "plasma": {
            "piston": {
                "density_over_n0": C.get("nt", _num(d, "particle_heater.foil.n0", C)) / n0
                if "nt" in C else _num(d, "target_injector.density", C) / n0,
                "theta_e_heat": theta_e_heat,
                "theta_e_init": C.get("theta_e_init",
                                      _num(d, "piston_electrons.ux_std", C, 0) ** 2),
                "theta_i_init": C.get("theta_i_init",
                                      _num(d, "piston_ions.ux_std", C, 0) ** 2),
            },
            "ambient": {
                "density_over_n0": C.get("namb", _num(d, "amb_electrons.density_function", C, n0)) / n0
                if "namb" in C else 0.01,
                "theta_0": theta0,
            },
        },
        "field": {"orientation": "perpendicular", "vA_over_c": C.get("vA", 0) / clight},
        "geometry": {
            "dims": int(float(d.get("geometry.dims", "1"))),
            "normal_axis": "z",
            "slab_halfwidth_di": (slab / di) if (slab and di) else None,
            "domain_halfwidth_de": (prob_hi / de) if de else None,
            "dz_over_de": (dz / de) if de else None,
            "boundary": d.get("boundary.field_lo", "periodic"),
        },
        "numerics": {
            "cfl": _num(d, "warpx.cfl", C, 0.75),
            "particle_shape": int(float(d.get("algo.particle_shape", "2"))),
            "max_step": int(float(d.get("max_step", "0"))),
            "ppc": {
                "piston": int(float(d.get("piston_electrons.num_particles_per_cell_each_dim", "0"))),
                "ambient": int(float(d.get("amb_electrons.num_particles_per_cell_each_dim", "0"))),
            },
        },
        "operators": {
            "heater": {"species": d.get("particle_heater.species", "piston_electrons"),
                       "intervals": int(float(d.get("particle_heater.intervals", "20")))},
            "injector": {"species": d.get("target_injector.species", "piston_electrons"),
                         "neutralizing_species": d.get("target_injector.neutralizing_species",
                                                       "piston_ions"),
                         "intervals": int(float(d.get("target_injector.intervals", "20"))),
                         "tau_over_wpe_inv": _num(d, "target_injector.tau", C, 0)
                         * C.get("wpe", 0) if "wpe" in C else None},
        },
    }
    return cfg


PRIMARY_PATHS = [
    ("reference", "n0"), ("reference", "mass_ratio"),
    ("field", "vA_over_c"),
    ("numerics", "max_step"), ("geometry", "dz_over_de"), ("geometry", "domain_halfwidth_de"),
]


def diff_configs(gen: dict, ref: dict, rtol=1e-3) -> list[str]:
    """Report primary-quantity mismatches between the generated and authored configs."""
    warns = []
    checks = PRIMARY_PATHS + [
        ("plasma", "piston", "theta_e_heat"), ("plasma", "piston", "density_over_n0"),
        ("plasma", "ambient", "theta_0"), ("plasma", "ambient", "density_over_n0"),
    ]
    for path in checks:
        g, r = _dig(gen, path), _dig(ref, path)
        if g is None or r is None:
            continue
        try:
            if abs(float(g) - float(r)) > rtol * max(abs(float(r)), 1e-30):
                warns.append(f"{'/'.join(path)}: generated {g:.6g} vs config {r:.6g}")
        except (TypeError, ValueError):
            if g != r:
                warns.append(f"{'/'.join(path)}: generated {g!r} vs config {r!r}")
    return warns


def _dig(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", help="run directory (containing warpx_used_inputs and config.yaml)")
    ap.add_argument("--inputs", help="parse this inputs file instead of warpx_used_inputs")
    ap.add_argument("--write", action="store_true", help="overwrite config.yaml (else write .generated.yaml)")
    args = ap.parse_args()

    inputs_path = args.inputs or os.path.join(args.run_dir, "warpx_used_inputs")
    if not os.path.isfile(inputs_path):
        sys.exit(f"inputs file not found: {inputs_path}\n"
                 f"(run WarpX first, or pass --inputs <deck>)")

    d = parse_inputs(inputs_path)
    d["__path__"] = os.path.join(args.run_dir, "config.yaml")
    consts = resolve_constants(d)
    gen = build_config(d, consts)

    ref_path = os.path.join(args.run_dir, "config.yaml")
    if os.path.isfile(ref_path):
        with open(ref_path) as fh:
            ref = yaml.safe_load(fh)
        warns = diff_configs(gen, ref)
        print("config verification:",
              "OK (primaries match what was simulated)" if not warns else "MISMATCH")
        for w in warns:
            print("  !", w)

    out = ref_path if args.write else os.path.join(args.run_dir, "config.generated.yaml")
    with open(out, "w") as fh:
        yaml.safe_dump(gen, fh, sort_keys=False, default_flow_style=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
