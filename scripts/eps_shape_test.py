#!/usr/bin/env python3
"""Read the SHAPE of wedge(eps) -- and test whether electron heating is the mediator.

Two problems with the first pass:
  1. the wedge metric stepped in 0.25 d_i0 bins, so a 1.0 d_i0 total range was 4 resolution
     elements -- far too coarse to tell a threshold from a power law. Here the scan step is
     0.02 d_i0 and the crossing is linearly interpolated, giving a continuous depth.
  2. eps is ONE parameter, so BOTH candidate mechanisms are functions of it and the ladder
     alone cannot attribute. What CAN discriminate is the mediator: both candidates act by
     changing how hot the shocked electrons get. So measure T_e,shocked(eps) directly and
     ask whether wedge tracks it -- and whether T_e,shocked itself shows a relativistic
     knee (flat, then a drop once gamma departs from 1) or a smooth power law in eps.

Fits both candidate laws to the continuous wedge and reports which describes it better:
    (a) relativistic capping   wedge ~ a + b/gamma_sh   (knee where gamma departs from 1)
    (b) wave regime            wedge ~ a * (w_pe/w_ce)^b  (pure power law, no knee)
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock.units import C, ME, QE  # noqa: E402

T_STAR = 0.260
TOL = 0.030
STEP = 0.02          # d_i0, wedge scan resolution (was 0.25)
SPREAD = 0.06        # v_z std / v_sh below which the layer counts as undisturbed

ROWS = [
    ("470 eV  anchor", "runs/S_phase/ss_dz16_ppc100", 0.0303, 100.0, 1.0014),
    ("1.5 keV",        "runs/E_phase/es_1p5keV",      0.0539,  56.2, 1.0044),
    ("4.7 keV",        "runs/E_phase/es_4p7keV",      0.0959,  31.6, 1.0138),
    ("15 keV",         "runs/E_phase/es_15keV",       0.1705,  17.8, 1.0436),
    ("47 keV",         "runs/E_phase/es_47keV",       0.3033,  10.0, 1.1380),
    ("R1_paper 47keV", "runs/R1_phase/R1_paper",      0.3033,  10.0, 1.1380),
]

out = []
for label, d, eps, wpwc, gam_th in ROWS:
    cfg = kinshock.load(d)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(d)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    i = int(np.argmin(np.abs(times - T_STAR)))
    if abs(times[i] - T_STAR) > TOL:
        print(f"  {label}: no frame within {TOL} of t*={T_STAR}")
        continue
    fr = kio.load_frame(pfs[i])
    ad = fr.ds.all_data()

    z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
    p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
    g = np.sqrt(1.0 + (p / (mi * C)) ** 2)
    v = p / (mi * g) / sc.vsh_model
    zp = np.percentile(np.asarray(ad[("piston_ions", "particle_position_x")]),
                       99.5) / sc.di0

    # continuous wedge depth: first offset where the local v_z spread drops below SPREAD,
    # linearly interpolated between the two bracketing samples instead of snapped to a bin
    offs = np.arange(0, 6, STEP)
    sig = np.full(offs.size, np.nan)
    for k, lo in enumerate(offs):
        m = (z - zp >= lo) & (z - zp < lo + 0.25)      # window width kept at 0.25
        if m.sum() > 100:
            sig[k] = v[m].std()
    depth = np.nan
    ok = np.isfinite(sig)
    for k in range(1, offs.size):
        if not (ok[k] and ok[k - 1]):
            continue
        if sig[k] < SPREAD <= sig[k - 1]:
            f = (sig[k - 1] - SPREAD) / (sig[k - 1] - sig[k])
            depth = offs[k - 1] + f * STEP
            break

    # shocked-layer electrons: the mediator both candidates act through. Take the layer
    # between the piston front and the wedge tip, in the run's OWN units.
    ze = np.asarray(ad[("amb_electrons", "particle_position_x")]) / sc.di0
    pex = np.asarray(ad[("amb_electrons", "particle_momentum_x")])
    pey = np.asarray(ad[("amb_electrons", "particle_momentum_y")])
    pez = np.asarray(ad[("amb_electrons", "particle_momentum_z")])
    hi = zp + (depth if np.isfinite(depth) else 1.0)
    ms = (ze >= zp) & (ze < hi)
    p2 = (pex[ms] ** 2 + pey[ms] ** 2 + pez[ms] ** 2)
    ge = np.sqrt(1.0 + p2 / (ME * C) ** 2)
    # thermal spread of the ELECTRON population, normalized to the run's own C_s,ab so it
    # is a dimensionless temperature comparable across rungs
    ue = np.sqrt(p2) / (ME * ge)
    Te_norm = float(np.var(ue) / sc.Cs_ab ** 2) if ms.sum() > 100 else np.nan
    gam_sh = float(np.mean(ge)) if ms.sum() > 100 else np.nan
    # Far-upstream control: must stay cold at every rung or the comparison is invalid.
    # The window is ABSOLUTE, not piston-relative. `ze > zp + 4` lands in very different
    # places in a 12.02 d_i0 box than in R1_paper's 80.5 d_i0 one -- the outer half next to
    # the open boundary vs genuinely far material -- and that alone made es_47keV and
    # R1_paper disagree 1.7x at the SAME eps (2026-08-17). On 5-11 d_i0 they agree to 0.5%.
    # gamma-corrected to match the shocked-layer definition above; that correction is only
    # <=4.5% here, so it was NOT the cause of the disagreement, but the two must still use
    # one definition. See scripts/eps_upstream_control.py.
    mu_ = (ze > 5.0) & (ze < 11.0)
    p2u = pex[mu_] ** 2 + pey[mu_] ** 2 + pez[mu_] ** 2
    gu = np.sqrt(1.0 + p2u / (ME * C) ** 2)
    Te_up = float(np.var(np.sqrt(p2u) / (ME * gu)) / sc.Cs_ab ** 2) if mu_.sum() > 100 else np.nan

    out.append((label, eps, wpwc, gam_th, depth, Te_norm, gam_sh, Te_up))
    print(f"  {label:16} depth={depth:.3f}  Te_sh={Te_norm:.3g}  gam_sh={gam_sh:.4f}",
          flush=True)

print("\n" + "=" * 104)
print(f"{'run':17}{'eps':>8}{'wpe/wce':>9}{'gam(Te)':>9}{'wedge':>8}"
      f"{'Te_shock':>11}{'gam_shock':>11}{'Te_upstr':>11}")
print("=" * 104)
for r in out:
    lab, eps, wpwc, gth, dep, Te, gsh, Tup = r
    print(f"{lab:17}{eps:>8.4f}{wpwc:>9.1f}{gth:>9.4f}{dep:>8.3f}"
          f"{Te:>11.4g}{gsh:>11.4f}{Tup:>11.4g}")
print("=" * 104)

lad = [r for r in out if not r[0].startswith("R1_paper")]
e = np.array([r[1] for r in lad])
w = np.array([r[4] for r in lad])
wp = np.array([r[2] for r in lad])
gs = np.array([r[6] for r in lad])
m = np.isfinite(w)

print("\nSHAPE TEST -- eps is one parameter, so both laws are functions of it; what")
print("separates them is whether the curve has a relativistic KNEE or is a pure power law.")


def rms(pred, obs):
    return float(np.sqrt(np.mean((pred - obs) ** 2)))


# (b) pure power law in w_pe/w_ce
A = np.vstack([np.log(wp[m]), np.ones(m.sum())]).T
b_pow, a_pow = np.linalg.lstsq(A, np.log(w[m]), rcond=None)[0]
pred_b = np.exp(a_pow) * wp ** b_pow
print(f"\n  (b) power law   wedge = {np.exp(a_pow):.3f} * (w_pe/w_ce)^{b_pow:+.3f}"
      f"    rms = {rms(pred_b[m], w[m]):.4f}")

# (a) relativistic capping: wedge ~ a + b/gamma_shocked (measured, not assumed)
if np.all(np.isfinite(gs)):
    A2 = np.vstack([1.0 / gs[m], np.ones(m.sum())]).T
    b_rel, a_rel = np.linalg.lstsq(A2, w[m], rcond=None)[0]
    pred_a = a_rel + b_rel / gs
    print(f"  (a) rel. cap    wedge = {a_rel:.3f} + {b_rel:.3f}/gamma_shocked"
          f"    rms = {rms(pred_a[m], w[m]):.4f}")

print("\n  observed vs each fit:")
print(f"    {'run':17}{'obs':>8}{'pow(b)':>9}{'rel(a)':>9}")
for k, r in enumerate(lad):
    pa = pred_a[k] if np.all(np.isfinite(gs)) else np.nan
    print(f"    {r[0]:17}{w[k]:>8.3f}{pred_b[k]:>9.3f}{pa:>9.3f}")

r1 = [r for r in out if r[0].startswith("R1_paper")]
if r1 and np.isfinite(lad[-1][4]) and np.isfinite(r1[0][4]):
    print(f"\n  INDEPENDENT CHECK: es_47keV {lad[-1][4]:.3f} vs R1_paper {r1[0][4]:.3f} d_i0 "
          f"({abs(lad[-1][4]-r1[0][4])/r1[0][4]*100:.1f}% apart) -- same eps, different "
          f"domain / dz/d_e / N_D.")
