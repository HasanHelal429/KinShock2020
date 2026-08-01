#!/usr/bin/env python3
"""Structure tests for KinShock2020 — verify the code is correct independently of
any WarpX run. Runnable directly (``python tests/test_structures.py``) or under
pytest. Quenches concerns about: config loading, the units derivation reproducing
Schaeffer 2020 Table I, config->deck generation (make_inputs), and the metrics.
"""

from __future__ import annotations

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
R0_HALF = os.path.join(ROOT, "runs", "R0_half")
R1_WARM = os.path.join(ROOT, "runs", "R1_warm")
R1_COLL = os.path.join(ROOT, "runs", "R1_coll")


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
                 ("plasma", "ambient", "theta_0"), ("field", "B0_tesla")):
        a, b = r1, r0
        for k in path:
            a, b = a[k], b[k]
        assert a == b, f"R0 differs from R1 at {'/'.join(path)}: {b} != {a}"
    s1, s0 = kinshock.units.derive(r1), kinshock.units.derive(r0)
    assert abs(s0.MA - s1.MA) < 1e-9 and abs(s0.B0 - s1.B0) < 1e-30


def test_make_inputs_roundtrip():
    """config -> deck generation: the deck rendered from config must (a) resolve
    back to the config primaries and (b) be physically equivalent to the committed
    hand-written deck (same resolved constants + scalar settings)."""
    from kinshock import deck
    for run, name in ((R1, "inputs_kinshock_R1"), (R0, "inputs_kinshock_R0")):
        cfg = kinshock.load(run)
        # render must produce a parseable deck that resolves back to the config
        # (deck.verify renders internally and diffs against the committed deck).
        warns = deck.verify(cfg, os.path.join(run, name))
        assert warns == [], f"{run}: generated deck differs from committed deck: {warns}"
        # sanity: the generated deck contains the expected primaries
        text = deck.render(cfg)
        assert "particle_heater" in text and "target_injector" in text
        assert f"amr.n_cell        = {kinshock.units.derive(cfg).n_cell}" in text


def test_one_sided_half_domain():
    """The one-sided layout (runs/R0_half) must: halve the cell count, put the
    domain at [0, half] with a reflecting z=0 wall + open far boundary, and KEEP
    the foil at [-slab, +slab] so the PSC heating rate (H ~ 1/width) matches the
    full-domain rate — the domain then clips the heated region to [0, slab]. It
    must also round-trip (deck resolves back to config)."""
    from kinshock import deck
    cfg_half = kinshock.load(R0_HALF)
    cfg_full = kinshock.load(R0)
    sc_half = kinshock.units.derive(cfg_half)
    sc_full = kinshock.units.derive(cfg_full)
    # same physics primaries, but half the cells at the same dz / domain half-width
    assert cfg_half["geometry"]["domain_halfwidth_de"] == cfg_full["geometry"]["domain_halfwidth_de"]
    assert sc_half.n_cell == sc_full.n_cell // 2, (sc_half.n_cell, sc_full.n_cell)

    text = deck.render(cfg_half)
    assert "geometry.prob_lo  =  0." in text
    assert "boundary.field_lo    = pec" in text          # z=0 symmetry/foil wall
    assert "boundary.particle_lo = reflecting" in text
    assert "boundary.particle_hi = absorbing" in text     # open far boundary
    assert "boundary.field_hi    = pec" in text           # div-cleaner-safe with B0
    # foil width preserved (NOT clipped to [0, slab]) so the heating rate matches full
    assert "particle_heater.foil.lo          = -slab" in text
    assert "target_injector.lo                   = -slab" in text
    assert deck.verify(cfg_half, os.path.join(R0_HALF, "inputs_kinshock_R0_half")) == []

    # symmetric decks must be unaffected by the layout feature (default = symmetric)
    assert "geometry.prob_lo  = -half" in deck.render(cfg_full)


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
def test_collisional_twin_of_R1_warm():
    """R1_coll must be R1_warm's collisional twin: the ONLY config differences are
    the absolute density scale (n0, so the ambient sits at 1e18 cm^-3) and the
    `collisions` block — every dimensionless quantity is untouched. The paper's Table I
    collisionality (lambda_ab = wce_ab/nu_ei,ab = 20) must come out of the derivation
    exactly, the deck must round-trip, and nu_ei*dt must stay well under 1 for
    Takizuka-Abe.
    """
    from kinshock import deck, units
    warm, coll = kinshock.load(R1_WARM), kinshock.load(R1_COLL)
    sw, sc = units.derive(warm), units.derive(coll)

    # (a) the ambient is at the requested absolute density
    assert abs(sc.namb / 1e6 - 1e18) / 1e18 < 1e-12, f"n_e0 = {sc.namb/1e6:.4g} cm^-3"

    # (b) every dimensionless quantity matches R1_warm exactly
    for k in ("mass_ratio", "MA", "Mms", "beta_ab", "beta_0", "dt_wpe", "n_cell",
              "steps_per_wci0"):
        a, b = getattr(sw, k), getattr(sc, k)
        assert abs(a - b) <= 1e-9 * max(1.0, abs(a)), f"{k}: warm {a} vs coll {b}"
    for k in ("de", "di0", "rho_i0", "dz"):        # lengths scale as 1/sqrt(n0)
        assert abs(getattr(sw, k) / getattr(sc, k) - 1e4) / 1e4 < 1e-9, k

    # (c) the paper's Table I collisionality is hit, and is resolved in time
    assert abs(sc.lambda_ab - 20.0) < 1e-9, f"lambda_ab = {sc.lambda_ab}"
    assert sc.nu_ei_dt < 0.1, f"nu_ei*dt = {sc.nu_ei_dt}"
    assert sw.coulomb_log is None, "R1_warm must stay collisionless"
    # lambda_ab = wce/nu_ei is mfp/rho_e,ab, NOT mfp/d_e: with rho_e,ab = 27.9 d_e it
    # must land at mfp ~ 559 d_e = 5.6 d_i0, i.e. still collisionless at ion scales
    # (the premise of the paper's Fig. 13).
    assert abs(sc.mfp_ei_ab / sc.de - 558.6) < 0.5, f"mfp/de = {sc.mfp_ei_ab/sc.de}"
    assert sc.mfp_ei_ab / sc.di0 > 1.0, f"mfp/di0 = {sc.mfp_ei_ab/sc.di0}"

    # guard the distinction so the two "20"s can never be silently conflated: the
    # mfp = 20 d_e target is a DIFFERENT, ~28x more collisional run.
    alt = kinshock.load(R1_COLL)
    alt["collisions"]["target"] = {"quantity": "mfp_over_de", "value": 20.0}
    sa = units.derive(alt)
    assert abs(sa.mfp_ei_ab / sa.de - 20.0) < 1e-9
    assert sa.nu_ei_ab / sc.nu_ei_ab > 25, "mfp=20 d_e must be far more collisional"

    # (d) deck round-trip, and the collisions actually reach the deck
    warns = deck.verify(coll, os.path.join(R1_COLL, "inputs_kinshock_R1_coll"))
    assert warns == [], f"R1_coll deck differs from config: {warns}"
    text = deck.render(coll)
    names = [ln.split("=", 1)[1].split() for ln in text.splitlines()
             if ln.startswith("collisions.collision_names")][0]
    assert len(names) == 10, f"expected all 10 species pairs, got {names}"
    assert "col_ae_ai.species    = amb_electrons amb_ions" in text
    assert f"my_constants.coulomb_log = {sc.coulomb_log!r}" in text
    # a collisionless config must emit no collisions block at all
    assert "collisions.collision_names" not in deck.render(warm)


# --------------------------------------------------------------------------- #
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


def test_B0_is_primary_and_vA_is_derived():
    """B0 must be independent of the ambient density; v_A must scale as 1/sqrt(namb).

    This is the whole point of the inversion: under the old vA_over_c primary, B0 --
    and therefore wci0, the clock for every t*wci0 plot -- scaled as sqrt(namb), so a
    wrong ambient density silently rescaled time (RESULTS 2026-07-31).
    """
    import copy
    import math

    cfg = kinshock.load(R1)
    base = kinshock.units.derive(cfg)

    hot = copy.deepcopy(cfg)
    hot["plasma"]["ambient"]["density_over_n0"] *= 4.0
    quad = kinshock.units.derive(hot)

    assert quad.B0 == base.B0, "B0 must not move when the ambient density changes"
    assert quad.wci0 == base.wci0, "wci0 must not move when the ambient density changes"
    assert math.isclose(quad.vA, base.vA / 2.0, rel_tol=1e-12), \
        "v_A must be derived as B0/sqrt(mu0*namb*m_i)"


def test_legacy_vA_over_c_is_refused():
    """A stale config must fail loudly with the migration command, not silently."""
    import pytest

    cfg = kinshock.load(R1)
    cfg["field"].pop("B0_tesla")
    cfg["field"]["vA_over_c"] = 0.01
    with pytest.raises(KeyError, match="migrate_field_b0"):
        kinshock.units.derive(cfg)
