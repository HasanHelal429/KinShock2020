#!/usr/bin/env python3
"""Measure the ramp's E_z wavenumber spectrum -- the direct test of the magnetization story.

THE ARGUMENT SO FAR IS INFERENCE. The eps ladder shows the wedge tracks eps; the lever-arm
argument rules out relativistic capping (measured gamma_shocked moves 5.7% against a 44%
wedge change); that leaves the electron-scale wave regime, rho_e/lambda_D = w_pe/w_ce,
which runs 100 -> 10 across the ladder. This measures it instead of arguing it.

THE DISCRIMINATOR, AND WHY THIS k-AXIS. Both candidate branches are driven by the same
cross-field drift v_d (~ v_sh; and v_d/v_te0 is PRESERVED across the ladder, so the drive
is identical at every rung and only its magnetic organisation changes). Doppler matching
k*v_d = omega turns each branch into a prediction on the axis

    k* == k v_d / w_ce

  * ECDI / Bernstein (magnetized): resonances at omega = n*w_ce  ->  k* = 1, 2, 3, ...
    a COMB AT INTEGERS, with ~w_pe/w_ce teeth before it runs out at w_pe.
  * Buneman (unmagnetized):        resonance at omega ~ w_pe     ->  k* = w_pe/w_ce
    a SINGLE BROAD PEAK at 100 / 56 / 32 / 18 / 10 per rung, no harmonic structure.

So the same axis carries both predictions and they do not overlap: teeth at small integers
versus a lone peak at the (rung-dependent) plasma-frequency mark.

RESOLUTION, CHECKED BEFORE RUNNING. Resolving the comb fundamental needs a window longer
than 2*pi*v_d/w_ce, which is 0.14 d_i0 at every rung (that ratio is preserved). A 3 d_i0
window gives dk* = 0.29, about three samples per tooth, and dz gives k*_max = 264 (47 keV)
to 2611 (470 eV) -- both far above the w_pe/w_ce marks. Nothing here is resolution-limited.

Spectra are averaged incoherently over several frames around the matched time, and the
window is Hann-tapered after removing a linear trend, so the DC/ramp profile does not leak
into low k.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io as kio, units as U  # noqa: E402

T_STAR = 0.260
NAVG = 9           # frames averaged, centred on T_STAR
WIN_DI0 = 3.0      # ramp window width [d_i0]

R = "runs"
ROWS = [
    ("470 eV anchor", f"{R}/S_phase/ss_dz16_ppc100"),
    ("1.5 keV",       f"{R}/E_phase/es_1p5keV"),
    ("4.7 keV",       f"{R}/E_phase/es_4p7keV"),
    ("15 keV",        f"{R}/E_phase/es_15keV"),
    ("47 keV",        f"{R}/E_phase/es_47keV"),
]

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
fig, axes = plt.subplots(len(ROWS), 1, figsize=(9.5, 15.0), dpi=190,
                         facecolor=SURFACE, sharex=True)
summary = []

for r, (label, d) in enumerate(ROWS):
    cfg = kinshock.load(d)
    sc = U.derive(cfg)
    qe, me, c = U.QE, U.ME, U.C
    n0 = sc.namb
    wpe0 = math.sqrt(n0 * qe * qe / (U.EPS0 * me))
    wce = qe * sc.B0 / me
    ratio = wpe0 / wce
    v_d = sc.vsh_model
    kstar_unit = v_d / wce                     # k* = k * (v_d/w_ce)

    pfs = kio.field_plotfiles(d) or kio.plotfiles(d)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    i0 = int(np.argmin(np.abs(times - T_STAR)))
    lo_i = max(0, i0 - NAVG // 2)
    sel = list(range(lo_i, min(len(pfs), lo_i + NAVG)))

    P_acc, kk, nused = None, None, 0
    for i in sel:
        fr = kio.load_frame(pfs[i], fields=("Ez",))
        Ez = fr.comps.get("Ez")
        if Ez is None:
            continue
        z = fr.z_centers / sc.di0
        # ramp location: the compressed background field Bx peaks there
        ir = int(np.argmax(np.abs(fr.Bx)))
        # FIXED cell count, not a coordinate mask: the ramp moves between frames and a
        # coordinate window lands on a different number of cells each time, which makes the
        # FFT lengths differ and the incoherent average impossible (and, worse, would give
        # each frame a slightly different k grid if it were silently interpolated).
        nwin = int(round(WIN_DI0 / (z[1] - z[0]))) // 2 * 2
        a = min(max(0, ir - nwin // 2), max(0, Ez.size - nwin))
        seg = Ez[a:a + nwin].astype(float)
        if seg.size < 256:
            continue
        # remove a linear trend, then taper -- keeps the ramp profile out of low k
        x = np.arange(seg.size)
        seg = seg - np.polyval(np.polyfit(x, seg, 1), x)
        seg = seg * np.hanning(seg.size)
        dz_m = (fr.z_centers[1] - fr.z_centers[0])
        F = np.fft.rfft(seg)
        k = 2.0 * np.pi * np.fft.rfftfreq(seg.size, d=dz_m)      # [1/m]
        P = np.abs(F) ** 2
        P_acc = P if P_acc is None else P_acc + P
        kk, nused = k, nused + 1

    if not nused:
        print(f"  {label}: no Ez frames usable", flush=True)
        continue
    P = P_acc / nused
    ks = kk * kstar_unit                       # dimensionless k*
    good = ks > 0
    ks, P = ks[good], P[good]
    P = P / P.max()

    ax = axes[r]
    ax.set_facecolor("#111014")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASE)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_color(INK2)
    ax.loglog(ks, P, color="#ffd166", lw=0.9)
    # ECDI comb: integers
    for n in range(1, 11):
        ax.axvline(n, color="#5ad2ff", lw=0.7, alpha=0.45, ls=(0, (3, 3)))
    # Buneman / unmagnetized mark
    ax.axvline(ratio, color="#ff6b6b", lw=1.6, alpha=0.95)
    ax.text(ratio, 1.4, r"$\omega_{pe}/\omega_{ce}$" + f" = {ratio:.0f}",
            color="#ff6b6b", fontsize=8.5, ha="center")
    ax.set_xlim(0.3, 400)
    ax.set_ylim(1e-6, 3)
    ax.set_ylabel(f"{label}\n" + rf"$\omega_{{pe}}/\omega_{{ce}}$ = {ratio:.0f}"
                  "\n" + r"$|E_z(k)|^2$ (norm.)", color=INK2, fontsize=9)
    if r == len(ROWS) - 1:
        ax.set_xlabel(r"$k^* = k\,v_d/\omega_{ce}$   "
                      "(dashed = ECDI comb at integers, red = Buneman mark)",
                      color=INK2, fontsize=10)

    # where is the spectral peak, and is there power at the integer teeth?
    kpk = float(ks[np.argmax(P)])
    band = (ks > 0.7) & (ks < 10.5)
    lowband = float(P[band].sum()) if band.any() else 0.0
    hib = (ks > 0.5 * ratio) & (ks < 2.0 * ratio)
    hiband = float(P[hib].sum()) if hib.any() else 0.0
    summary.append((label, ratio, kpk, lowband, hiband, nused))
    print(f"  {label:16} frames={nused}  k*_peak={kpk:.2f}  "
          f"w_pe/w_ce={ratio:.0f}", flush=True)

fig.suptitle(r"Ramp $E_z$ wavenumber spectrum vs $\varepsilon$ — ECDI comb or Buneman peak?"
             "\nsame drive at every rung; only the magnetic organisation of the response "
             "changes", color=INK, fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.965])
os.makedirs("media/E_phase", exist_ok=True)
out = "media/E_phase/ramp_Ez_spectrum.png"
fig.savefig(out, facecolor=SURFACE)
print(f"\n  wrote {out}")

print("\n" + "=" * 84)
print(f"{'run':16}{'wpe/wce':>10}{'k*_peak':>10}{'P(k*~1-10)':>13}"
      f"{'P(k*~wpe/wce)':>16}{'ratio hi/lo':>14}")
print("=" * 84)
for lab, ratio, kpk, lo, hi, n in summary:
    print(f"{lab:16}{ratio:>10.0f}{kpk:>10.2f}{lo:>13.4g}{hi:>16.4g}"
          f"{(hi/lo if lo > 0 else float('nan')):>14.3g}")
print("=" * 84)
print("\nPrediction: magnetized (47 keV) -> power concentrated in the integer teeth,")
print("peak at small k*.  Unmagnetized (470 eV) -> peak near w_pe/w_ce, teeth absent.")
