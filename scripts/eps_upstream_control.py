#!/usr/bin/env python3
"""Far-upstream electron control for the eps ladder, on a MATCHED ABSOLUTE window.

The committed result (f9f9b9b) records a caveat: the far-upstream control is "not flat
across the ladder", and es_47keV and R1_paper disagree on it at the SAME eps. That caveat
is load-bearing -- if the ambient itself is at a different normalized temperature at each
rung, absolute normalized temperatures are not comparable across rungs.

Part of it looks like a metric bug rather than physics. eps_shape_test.py computes the
shocked-layer speed as p/(m gamma) but the upstream one as p/m:

    ue    = sqrt(p2) / (ME * ge)      # shocked layer, gamma-corrected
    ue_up = sqrt(...) / ME            # upstream, NOT gamma-corrected

At 470 eV that is nothing (v_te/c = 0.0044). At 47 keV the ambient is 1022 eV and
v_te/c = 0.045, and the tail of the distribution carries more, so the two definitions
diverge -- and they diverge BY RUNG, in the same direction as the reported drift. This
recomputes the control both ways on the same frames and says how much of the drift is the
definition and how much is real.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock.units import C, ME  # noqa: E402

T_STAR, TOL = 0.260, 0.030

ROWS = [
    ("470 eV  anchor", "runs/S_phase/ss_dz16_ppc100", 0.0303),
    ("1.5 keV",        "runs/E_phase/es_1p5keV",      0.0539),
    ("4.7 keV",        "runs/E_phase/es_4p7keV",      0.0959),
    ("15 keV",         "runs/E_phase/es_15keV",       0.1705),
    ("47 keV",         "runs/E_phase/es_47keV",       0.3033),
    ("R1_paper 47keV", "runs/R1_phase/R1_paper",      0.3033),
]

rows = []
for label, d, eps in ROWS:
    cfg = kinshock.load(d)
    sc = kinshock.units.derive(cfg)
    pfs = kio.plotfiles(d)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    i = int(np.argmin(np.abs(times - T_STAR)))
    if abs(times[i] - T_STAR) > TOL:
        print(f"  {label}: no frame within {TOL} of t*={T_STAR}")
        continue
    ad = kio.load_frame(pfs[i]).ds.all_data()
    ze = np.asarray(ad[("amb_electrons", "particle_position_x")]) / sc.di0
    px = np.asarray(ad[("amb_electrons", "particle_momentum_x")])
    py = np.asarray(ad[("amb_electrons", "particle_momentum_y")])
    pz = np.asarray(ad[("amb_electrons", "particle_momentum_z")])
    zp = np.percentile(np.asarray(ad[("piston_ions", "particle_position_x")]),
                       99.5) / sc.di0

    m = (ze > 5.0) & (ze < 11.0)   # MATCHED absolute window, fits the 12.02 d_i0 box
    p2 = px[m] ** 2 + py[m] ** 2 + pz[m] ** 2
    ge = np.sqrt(1.0 + p2 / (ME * C) ** 2)
    Te_bug = float(np.var(np.sqrt(p2) / ME) / sc.Cs_ab ** 2)          # as committed
    Te_fix = float(np.var(np.sqrt(p2) / (ME * ge)) / sc.Cs_ab ** 2)   # gamma-corrected
    # the DECK's own ambient temperature in the same normalization, i.e. what an
    # unheated upstream must return. theta_0 is kT/(m_e c^2); C_s,ab^2 = theta_ab c^2/mu.
    th0 = float(cfg["plasma"]["ambient"]["theta_0"])
    Te_deck = th0 * C ** 2 / sc.Cs_ab ** 2
    rows.append((label, eps, times[i], Te_bug, Te_fix, Te_deck, float(ge.mean()),
                 int(m.sum())))
    print(f"  {label:16} t*={times[i]:.3f}  Te_up(bug)={Te_bug:.4g}  "
          f"Te_up(fix)={Te_fix:.4g}  deck={Te_deck:.4g}", flush=True)

print("\n" + "=" * 100)
print(f"{'run':17}{'eps':>8}{'t*':>7}{'Te_up bug':>11}{'Te_up fix':>11}"
      f"{'deck T_0':>10}{'fix/deck':>10}{'<gamma>':>9}{'N':>10}")
print("=" * 100)
for lab, eps, t, tb, tf, td, g, n in rows:
    print(f"{lab:17}{eps:>8.4f}{t:>7.3f}{tb:>11.4g}{tf:>11.4g}{td:>10.4g}"
          f"{tf/td:>10.3f}{g:>9.5f}{n:>10d}")
print("=" * 100)

print("\nWhat this says:")
print("  Te_up(fix)/deck is the FAR-UPSTREAM HEATING FACTOR -- 1.0 means the ambient is")
print("  still at its initial temperature. That, not the raw Te_up, is the control that")
print("  has to be flat across the ladder for cross-rung comparison to be valid.")
if len(rows) >= 2:
    f = np.array([r[4] / r[5] for r in rows])
    b = np.array([r[3] for r in rows])
    fx = np.array([r[4] for r in rows])
    print(f"\n  raw Te_up(bug) spread across rungs : {b.max()/b.min():.2f}x")
    print(f"  raw Te_up(fix) spread across rungs : {fx.max()/fx.min():.2f}x")
    print(f"  heating factor spread              : {f.max()/f.min():.2f}x")
