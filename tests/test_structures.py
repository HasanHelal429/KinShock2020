#!/usr/bin/env python3
"""Structure tests for KinShock2020 — verify the code is correct independently of
any WarpX run. Runnable directly (``python tests/test_structures.py``) or under
pytest. Quenches concerns about: config loading, the units derivation reproducing
Schaeffer 2020 Table I, the config↔deck round-trip (make_config), and the metrics.
"""

from __future__ import annotations

import importlib.util
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import numpy as np
import kinshock
from kinshock import metrics

R1 = os.path.join(ROOT, "runs", "R1")
R0 = os.path.join(ROOT, "runs", "R0")


# --------------------------------------------------------------------------- #
def test_config_loads():
    cfg = kinshock.load(R1)
    for sec in ("reference", "plasma", "field", "geometry", "numerics", "model", "targets"):
        assert sec in cfg, f"missing section {sec}"
    assert cfg["reference"]["mass_ratio"] == 100.0
    assert cfg["plasma"]["piston"]["theta_e_heat"] == 0.092


def test_units_reproduce_table1():
    """units.derive from config primaries must reproduce Table I headline numbers."""
    cfg = kinshock.load(R1)
    sc = kinshock.units.derive(cfg)
    C = kinshock.units.C
    # M_A ~ 14, M_ms ~ 13 (the supercritical-shock headline numbers)
    assert abs(sc.MA - 14) / 14 < 0.15, f"M_A={sc.MA}"
    assert abs(sc.Mms - 13) / 13 < 0.20, f"M_ms={sc.Mms}"
    # sound / shock / piston speeds vs Table I
    assert abs(sc.Cs_ab / C - 0.030) < 0.002, f"Cs_ab/c={sc.Cs_ab/C}"
    assert abs(sc.vsh_model / C - 0.1395) < 0.005, f"vsh/c={sc.vsh_model/C}"
    assert abs(sc.vp_model / C - 0.104) < 1e-6
    # numerics
    assert abs(sc.dt_wpe - 0.225) < 1e-6, f"dt_wpe={sc.dt_wpe}"
    assert sc.n_cell == 50000, f"n_cell={sc.n_cell}"
    assert sc.mass_ratio == 100.0
    # d_i0 = 100 d_e for contrast 0.01 & mass ratio 100; rho_i0 ~ 1000 d_e
    assert abs(sc.di0 / sc.de - 100.0) < 1.0, f"di0/de={sc.di0/sc.de}"
    assert 900 < sc.rho_i0 / sc.de < 1200, f"rho_i0/de={sc.rho_i0/sc.de}"


def test_validation_clean_for_R1():
    cfg = kinshock.load(R1)
    warns = kinshock.config.validate(cfg)
    assert warns == [], f"unexpected validation warnings: {warns}"


def test_R0_matches_R1_physics():
    """R0 must share every physics primary with R1 (only grid/steps differ)."""
    r1 = kinshock.load(R1)
    r0 = kinshock.load(R0)
    for path in (("reference", "n0"), ("reference", "mass_ratio"),
                 ("plasma", "piston", "theta_e_heat"),
                 ("plasma", "piston", "density_over_n0"),
                 ("plasma", "ambient", "density_over_n0"),
                 ("plasma", "ambient", "theta_0"), ("field", "vA_over_c")):
        a, b = r1, r0
        for k in path:
            a, b = a[k], b[k]
        assert a == b, f"R0 differs from R1 at {'/'.join(path)}: {b} != {a}"
    s1, s0 = kinshock.units.derive(r1), kinshock.units.derive(r0)
    assert abs(s0.MA - s1.MA) < 1e-9 and abs(s0.B0 - s1.B0) < 1e-30


def test_make_config_roundtrip():
    """Parsing the deck and rebuilding the config must match the authored config."""
    mc = _load_script("make_config", os.path.join(ROOT, "scripts", "make_config.py"))
    for run, deck in ((R1, "inputs_kinshock_R1"), (R0, "inputs_kinshock_R0")):
        d = mc.parse_inputs(os.path.join(run, deck))
        d["__path__"] = os.path.join(run, "config.yaml")
        consts = mc.resolve_constants(d)
        gen = mc.build_config(d, consts)
        ref = kinshock.load(run)
        warns = mc.diff_configs(gen, ref)
        assert warns == [], f"{run}: round-trip mismatch: {warns}"


def test_metrics_eq1_expansion_speed():
    """Eq. 1: v(n_e,ab) = Cs/2 (1 - ln 1) = Cs/2 at the ablation density."""
    Cs = 0.03 * kinshock.units.C
    v = metrics.expansion_speed(np.array([2.5e16]), 2.5e16, 1e16, Cs)
    assert abs(float(v[0]) - 0.5 * Cs) < 1e-6 * Cs
    # denser than ablation -> slower; log arg <=0 -> NaN
    assert np.isnan(metrics.expansion_speed(np.array([1e16]), 2.5e16, 1e16, Cs))[0]


def test_metrics_rh_ratio():
    """Eq. 2: perpendicular RH speed ratio -> ~4 (strong-shock compression) at high M_A."""
    r = metrics.rh_speed_ratio(14, 0.4)
    assert 3.0 < r < 4.0, f"RH ratio={r}"
    assert metrics.rh_speed_ratio(50, 0.2) > metrics.rh_speed_ratio(5, 0.2)


def test_metrics_reflected_and_front():
    uz = np.array([-1., 0., 0.5, 1.5, 2.0])
    assert abs(metrics.reflected_fraction_G(uz, 1.0) - 2 / 5) < 1e-12
    zc = np.linspace(-100, 100, 201)
    n = np.ones_like(zc) * 1e16
    n[(zc > 40) & (zc < 45)] = 5e16          # a compression at z ~ 42
    zf = metrics.track_front(zc, n, 1e16, threshold=1.5, z_exclude=20.0)
    assert 40 <= zf <= 45, f"front={zf}"


def test_metrics_criteria_shape():
    cfg = kinshock.load(R1)
    sc = kinshock.units.derive(cfg)
    zc = np.linspace(-2000, 2000, 401) * sc.de
    n = np.full_like(zc, sc.namb)
    B = np.full_like(zc, sc.B0)
    j = (zc > 500 * sc.de) & (zc < 560 * sc.de)
    n[j] = 4 * sc.namb            # >2x density compression
    B[j] = 4 * sc.B0              # >2x field compression
    uz = np.array([2.0 * sc.vsh_model, 0.0])   # one reflected ambient ion
    res = metrics.evaluate_criteria(
        zc=zc, n_e=n, Bmag=B, uz_ambient=uz, z_ambient=np.array([540 * sc.de, 0.0]),
        scales=sc, vsh=sc.vsh_model, v_front=sc.vsh_model,
        piston_field_z=0.0, front_z=540 * sc.de)
    assert set(res.flags) == {f"{i}_{n}" for i, n in enumerate(
        ["super_magnetosonic", "collisionless", "density_compression", "field_compression",
         "steep_ramp", "reflected_ions", "piston_separation"], start=1)}
    assert res.flags["3_density_compression"] and res.flags["4_field_compression"]
    assert res.flags["6_reflected_ions"] and res.flags["7_piston_separation"]


# --------------------------------------------------------------------------- #
def _load_script(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_all():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    npass = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            npass += 1
        except Exception as e:  # noqa: BLE001
            print(f"  FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{npass}/{len(tests)} structure tests passed")
    return npass == len(tests)


if __name__ == "__main__":
    print("KinShock2020 structure tests\n" + "=" * 40)
    sys.exit(0 if _run_all() else 1)
