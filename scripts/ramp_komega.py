#!/usr/bin/env python3
"""(k, omega) spectrum of E_z in the shock ramp -- the test the k-only version could not do.

WHY k ALONE FAILED (RESULTS 2026-08-19). The first attempt plotted |E_z(k)|^2 against
k* = k v_d/w_ce and looked for an ECDI comb at integers. It could not work, for three
reasons that are worth keeping written down:

  1. k*_Nyquist/(w_pe/w_ce) = 26.1 IDENTICALLY at every rung -- both scale the same way
     with eps -- so a numerical cutoff at a fixed fraction of Nyquist and a physical one at
     w_pe/w_ce are indistinguishable on that axis.
  2. The summary metrics were artifacts: a "hi/lo band ratio" whose bands scaled with the
     rung (per k-sample it was 1.00 +- 0.05 everywhere), and an argmax on a flat noisy
     plateau.
  3. Worst: the E_phase rungs stop at t* = 0.302 and the shock does not form until
     t* ~ 1.4. There was no ramp in them at all -- argmax|Bx| was finding piston-driven
     compression.

THE omega AXIS IS IMMUNE TO (1). A grid cutoff acts in k, not in omega, so discrete peaks
at omega = n*w_ce cannot be manufactured by the mesh. ECDI/Bernstein predicts power at
integer omega/w_ce; Buneman predicts a single blob at omega/w_ce = w_pe/w_ce. Both lie on
the Doppler ridge omega = k v_d, so the ridge is a consistency check, not the measurement.

TWO THINGS THIS SCRIPT REFUSES TO PAPER OVER:
  * It prints omega_Nyquist/w_ce and, if w_pe/w_ce exceeds it, says outright that the
    Buneman line is ALIASED and that low-order peaks may be aliases rather than harmonics.
    On the existing 470 eV production run that is the case (Nyquist 4.7 vs w_pe/w_ce 100).
  * The ramp window co-moves with the shock by INTEGER cell shifts. A fractional shift
    would need interpolation, which smears exactly the harmonic structure being looked for.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io as kio, units as U  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--t0", type=float, default=1.40, help="window start in t*wci0")
    ap.add_argument("--nframes", type=int, default=64)
    ap.add_argument("--win-di0", type=float, default=3.0, help="ramp window width [d_i0]")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = U.derive(cfg)
    run_id = cfg["meta"]["run_id"]
    wce = U.QE * sc.B0 / U.ME
    n0 = sc.namb
    wpe0 = math.sqrt(n0 * U.QE ** 2 / (U.EPS0 * U.ME))
    ratio = wpe0 / wce
    v_d = sc.vsh_model

    pfs = kio.field_plotfiles(args.run_dir) or kio.plotfiles(args.run_dir)
    times = np.array([kio.load_frame(p).time for p in pfs])
    tstar = times * sc.wci0
    i0 = int(np.argmin(np.abs(tstar - args.t0)))
    sel = list(range(i0, min(len(pfs), i0 + args.nframes)))
    if len(sel) < 8:
        sys.exit(f"only {len(sel)} frames from t*={args.t0}; need more")

    dt_f = float(np.median(np.diff(times[sel])))
    wce_dt = wce * dt_f
    w_nyq = math.pi / wce_dt
    dw = 2.0 * math.pi / (len(sel) * wce_dt)
    print(f"{run_id}:  w_pe/w_ce = {ratio:.1f}")
    print(f"  frames {len(sel)}  t* {tstar[sel[0]]:.3f} -> {tstar[sel[-1]]:.3f}"
          f"  ({len(sel)*wce_dt/(2*math.pi):.2f} electron gyroperiods)")
    print(f"  w_ce*dt_field = {wce_dt:.4f}   w_Nyquist/w_ce = {w_nyq:.2f}"
          f"   dw/w_ce = {dw:.3f}")
    if ratio > w_nyq:
        print(f"  !! w_pe/w_ce = {ratio:.0f} EXCEEDS omega-Nyquist ({w_nyq:.1f}):")
        print( "     the Buneman line is ALIASED. Peaks below Nyquist may be aliases,")
        print( "     NOT cyclotron harmonics. Treat a positive detection as unproven.")
    else:
        print(f"  w_pe/w_ce = {ratio:.1f} is inside Nyquist -- both branches observable.")

    # ---- assemble Ez(t, z_ramp) on a co-moving, integer-shifted window ----------------
    rows, tt = [], []
    nwin = None
    for i in sel:
        fr = kio.load_frame(pfs[i], fields=("Ez",))
        Ez = fr.comps.get("Ez")
        if Ez is None:
            continue
        z = fr.z_centers / sc.di0
        if nwin is None:
            nwin = int(round(args.win_di0 / (z[1] - z[0]))) // 2 * 2
        # RAMP LOCATOR -- a KINEMATIC PRIOR, not a global extremum. Verified 2026-08-19:
        #   argmax|Bx| over the whole domain returns z ~ 0.1 d_i0 every frame -- it locks
        #   onto a 13x B0 spike at the FOIL WALL, not the shock, and the window then sits
        #   still while the shock runs away (median jump 7 cells/frame against an expected
        #   277). That is what made the first (k,omega) attempt come out white.
        #   "outermost |Bx| > 1.5 B0" fails too: past t* ~ 1.36 the precursor has filled
        #   the box and that test returns the far boundary.
        # Restricting argmax to model_front +/- 4 d_i0 tracks at 250-270 cells/frame.
        zf = sc.MA * (fr.time * sc.wci0)
        band = np.nonzero((z > zf - 4.0) & (z < zf + 4.0))[0]
        if band.size < nwin // 4:
            continue
        ir = int(band[int(np.argmax(np.abs(fr.Bx)[band]))])
        a = min(max(0, ir - nwin // 2), max(0, Ez.size - nwin))
        seg = Ez[a:a + nwin].astype(float)
        if seg.size != nwin:
            continue
        x = np.arange(seg.size)
        rows.append(seg - np.polyval(np.polyfit(x, seg, 1), x))
        tt.append(fr.time)
    D = np.asarray(rows)
    if D.shape[0] < 8:
        sys.exit("not enough usable Ez frames")
    nt, nz = D.shape
    print(f"  window {nwin} cells = {args.win_di0} d_i0;  array {nt} x {nz}")

    D = D * np.hanning(nt)[:, None] * np.hanning(nz)[None, :]
    dz_m = float(np.diff(kio.load_frame(pfs[sel[0]]).z_centers)[0])
    F = np.fft.fftshift(np.fft.fft2(D))
    P = np.abs(F) ** 2
    kax = np.fft.fftshift(np.fft.fftfreq(nz, d=dz_m)) * 2 * np.pi
    wax = np.fft.fftshift(np.fft.fftfreq(nt, d=dt_f)) * 2 * np.pi
    ks = kax * (v_d / wce)          # k*
    ws = wax / wce                  # omega/w_ce

    # keep positive omega; power is symmetric
    pos = ws > 0
    ws_p, P_p = ws[pos], P[pos, :]
    wspec = P_p.sum(axis=1)
    wspec = wspec / wspec.max()

    SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 9.5), dpi=190,
                                   facecolor=SURFACE,
                                   gridspec_kw=dict(height_ratios=[2.1, 1.0]))
    for ax in (ax1, ax2):
        ax.set_facecolor("#111014")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(BASE)
        ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            t.set_color(INK2)

    ext = [ks[0], ks[-1], ws_p[0], ws_p[-1]]
    ax1.imshow(np.ma.masked_where(P_p <= 0, P_p), origin="lower", aspect="auto",
               extent=ext, norm=LogNorm(), cmap="magma")
    kk = np.linspace(max(ks[0], 0), ks[-1], 50)
    ax1.plot(kk, kk, color="#5ad2ff", lw=1.0, ls=(0, (5, 4)), alpha=0.8)
    ax1.text(0.98, 0.03, r"dashed: $\omega = k v_d$", transform=ax1.transAxes,
             ha="right", color="#5ad2ff", fontsize=8.5)
    for n in range(1, int(min(w_nyq, 12)) + 1):
        ax1.axhline(n, color="#ffd166", lw=0.6, alpha=0.35)
    if ratio <= w_nyq:
        ax1.axhline(ratio, color="#ff6b6b", lw=1.5)
    ax1.set_xlim(0, min(ks[-1], 40))
    ax1.set_ylim(0, min(ws_p[-1], 32))
    ax1.set_xlabel(r"$k^* = k\,v_d/\omega_{ce}$", color=INK2, fontsize=10)
    ax1.set_ylabel(r"$\omega/\omega_{ce}$", color=INK2, fontsize=10)

    ax2.semilogy(ws_p, wspec, color="#ffd166", lw=1.0)
    for n in range(1, int(min(w_nyq, 32)) + 1):
        ax2.axvline(n, color="#5ad2ff", lw=0.6, alpha=0.4, ls=(0, (3, 3)))
    if ratio <= w_nyq:
        ax2.axvline(ratio, color="#ff6b6b", lw=1.6)
        ax2.text(ratio, 1.3, r"$\omega_{pe}/\omega_{ce}$", color="#ff6b6b",
                 fontsize=8.5, ha="center")
    ax2.set_xlim(0, min(ws_p[-1], 32))
    ax2.set_xlabel(r"$\omega/\omega_{ce}$   (dashed = cyclotron harmonics)",
                   color=INK2, fontsize=10)
    ax2.set_ylabel(r"$\int |E_z(k,\omega)|^2 dk$", color=INK2, fontsize=10)

    warn = ("  ALIASED: w_pe/w_ce > Nyquist" if ratio > w_nyq else "")
    fig.suptitle(f"Ramp $E_z$ $(k,\\omega)$ — {run_id}\n"
                 rf"$\omega_{{pe}}/\omega_{{ce}}$ = {ratio:.0f},  "
                 rf"$\omega_{{Nyq}}/\omega_{{ce}}$ = {w_nyq:.1f},  "
                 rf"$\Delta\omega/\omega_{{ce}}$ = {dw:.2f}{warn}",
                 color=INK, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("media/E_phase", exist_ok=True)
    out = f"media/E_phase/komega_{args.tag or run_id}.png"
    fig.savefig(out, facecolor=SURFACE)
    print(f"  wrote {out}")

    print(f"\n  {'n':>4}{'omega/w_ce':>12}{'P(n)/P_median':>16}")
    med = float(np.median(wspec))
    for n in range(1, int(min(w_nyq, 12)) + 1):
        j = int(np.argmin(np.abs(ws_p - n)))
        print(f"  {n:>4}{ws_p[j]:>12.2f}{wspec[j]/med:>16.2f}")
    print(f"\n  a comb shows as P(n)/median >> 1 at consecutive integers.")


if __name__ == "__main__":
    main()
