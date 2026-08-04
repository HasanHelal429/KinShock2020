#!/usr/bin/env python3
"""Characterize a run's magnetic-field fluctuations — physical vs numerical.

A supercritical perpendicular shock produces *real* electromagnetic turbulence in
the reflected-ion foot/ramp; a PIC run can *also* grow grid-scale numerical field
noise (e.g. from an under-resolved Debye length, dz/lambda_D > ~pi). This driver
separates the two for a completed run using four independent tests and writes a
multi-panel figure to ``media/<run_id>/bfield_diagnostic.png``:

  1. Spatial power spectrum per zone — physical modes peak at a plasma scale and
     decline toward the grid; grid noise rises toward the Nyquist/filter scale.
  2. Polarization (hodogram corr(Bx,By)) — a coherent EM wave is elliptically
     polarized; broadband/random ~ 0 correlation.
  3. Particle response — a real dB/B~1 fluctuation scatters/heats the plasma;
     if ions AND electrons sit at the t=0 thermal floor while dB~B0, the field
     is not exchanging energy with the plasma (numerical field noise).
  4. Debye resolution — reports dz/lambda_D (the enabling condition for grid heating).

Usage:
    python scripts/bfield_diagnostic.py runs/R1_phase/R1_core            # last frame
    python scripts/bfield_diagnostic.py runs/R1_phase/R1_core --twci 1.4 # nearest frame to t*wci
    python scripts/bfield_diagnostic.py runs/R1_phase/R1_core --frame 25 --out media/testing/bf.png
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import io  # noqa: E402
from kinshock.units import C  # noqa: E402


def _smooth(a, k=15):
    return np.convolve(a, np.ones(k) / k, mode="same")


def _front(fr, sc, de):
    """Forward (+z) shock front: furthest z where smoothed ambient density > 1.5 n0."""
    zc = np.asarray(fr.z_centers) / de
    na = io.species_density(fr, "amb_ions") / sc.namb
    m = zc > 50
    above = zc[m][_smooth(na[m]) > 1.5]
    return above.max() if above.size else np.nan


def _thermal_spread(frame, sp, de, me):
    """RMS momentum per component (in c) and count, over all particles of species."""
    ad = frame.ds.all_data()
    try:
        u = np.stack([np.asarray(ad[(sp, f"particle_momentum_{c}")]) / (me * C) for c in "xyz"])
    except Exception:
        return None
    return np.sqrt((u ** 2).sum(0).mean() / 3.0), u.shape[1]


def _spread_in(frame, sp, de, me, zlo, zhi):
    ad = frame.ds.all_data()
    try:
        z = np.asarray(ad[(sp, "particle_position_x")]) / de
        u = np.stack([np.asarray(ad[(sp, f"particle_momentum_{c}")]) / (me * C) for c in "xyz"])
    except Exception:
        return None
    m = (z > zlo) & (z < zhi)
    if m.sum() < 50:
        return None
    return np.sqrt((u[:, m] ** 2).sum(0).mean() / 3.0), int(m.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=os.path.join(ROOT, "runs", "R1_phase", "R1_core"))
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--frame", type=int, help="plotfile index (default: last)")
    g.add_argument("--twci", type=float, help="use frame nearest this t*wci0")
    ap.add_argument("--out", help="figure path (default media/<run_id>/bfield_diagnostic.png)")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    de = sc.de
    rid = cfg["meta"]["run_id"]
    me = sc.mi / sc.mass_ratio
    pf = io.plotfiles(args.run_dir)
    if not pf:
        sys.exit(f"no plotfiles in {args.run_dir}")

    if args.twci is not None:
        idx = int(np.argmin([abs(io.load_frame(p).time * sc.wci0 - args.twci) for p in pf]))
    else:
        idx = args.frame if args.frame is not None else len(pf) - 1
    fr = io.load_frame(pf[idx])
    fr0 = io.load_frame(pf[0])
    zc = np.asarray(fr.z_centers) / de
    dz = float(zc[1] - zc[0])
    bx = fr.Bx / sc.B0
    by = fr.By / sc.B0
    front = _front(fr, sc, de)
    rho_i = sc.rho_i0 / de

    # Debye resolution (enabling condition for grid heating)
    lamD = _thermal_spread(fr0, "amb_electrons", de, me)[0]  # v_the/c == lambda_D/d_e
    print(f"=== B-field diagnostic: {rid}  frame {idx}  t*wci={fr.time*sc.wci0:.2f} ===")
    print(f"front={front:.0f} d_e  rho_i0={rho_i:.0f} d_e  dz={dz:.2f} d_e")
    print(f"lambda_D={lamD:.3f} d_e  ->  dz/lambda_D={dz/lamD:.1f}  "
          f"(finite-grid heating threshold ~pi=3.1)")

    # zones
    zones = {
        "downstream": (front - 800, front - 300, "purple"),
        "foot (reflected ions)": (front + 50, front + 300, "red"),
        "far upstream": (front + 600, min(front + 1400, zc.max() - 20), "green"),
    }

    def psd(seg):
        seg = (seg - seg.mean()) * np.hanning(len(seg))
        k = np.fft.rfftfreq(len(seg), d=dz)
        return k[1:], (np.abs(np.fft.rfft(seg)) ** 2)[1:]

    print("\nzone                      <Bperp>  corr(Bx,By)  peak-lambda[d_e]  %pwr(lambda<2de)")
    spec = {}
    for lbl, (lo, hi, c) in zones.items():
        m = (zc > lo) & (zc < hi)
        if m.sum() < 40:
            print(f"  {lbl}: too few cells"); continue
        k, P = psd(bx[m])
        lam = 1 / k
        cc = np.corrcoef(bx[m], by[m])[0, 1]
        bp = np.hypot(fr.Bx, fr.By)[m].mean() / sc.B0
        spec[lbl] = (k, P, c, bx[m], by[m])
        print(f"  {lbl:24s} {bp:5.2f}   {cc:+.2f}       {lam[np.argmax(P)]:7.1f}         "
              f"{100*P[lam<2].sum()/P.sum():5.1f}%")

    # particle-response test vs distance ahead of front
    print("\nparticle response ahead of front (dB/B~1 should scatter a real wave):")
    print("  dist[rho_i]   dBx_rms   ion u/floor   e- u/floor")
    i0 = _thermal_spread(fr0, "amb_ions", de, me)[0]
    e0 = _thermal_spread(fr0, "amb_electrons", de, me)[0]
    for d0, d1 in [(50, 300), (300, 600), (600, 1000), (1000, 1400)]:
        m = (zc > front + d0) & (zc < front + d1)
        if m.sum() < 20:
            continue
        rms = (bx[m] - bx[m].mean()).std()
        ri = _spread_in(fr, "amb_ions", de, me, front + d0, front + d1)
        re = _spread_in(fr, "amb_electrons", de, me, front + d0, front + d1)
        fi = ri[0] / i0 if ri else float("nan")
        fe = re[0] / e0 if re else float("nan")
        print(f"  {(d0+d1)/2/rho_i:5.2f}        {rms:.3f}     x{fi:.2f}         x{fe:.2f}")

    # ---- figure ----
    out = args.out or os.path.join(ROOT, "media", rid, "bfield_diagnostic.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2)
    # (0,:) profile
    ax = fig.add_subplot(gs[0, :])
    mp = zc > 0
    ax.plot(zc[mp], bx[mp], lw=0.5, color="tab:blue", label="Bx/B0")
    ax.plot(zc[mp], by[mp], lw=0.5, color="tab:orange", alpha=0.7, label="By/B0")
    ax.axvline(front, color="k", ls=":", label="shock front")
    for lbl, (lo, hi, c) in zones.items():
        ax.axvspan(max(lo, 0), hi, color=c, alpha=0.08)
    ax.set_xlim(0, zc.max()); ax.set_xlabel("z [d_e]"); ax.set_ylabel("B/B0")
    ax.legend(fontsize=8, ncol=3)
    ax.set_title(f"{rid}  t*wci={fr.time*sc.wci0:.2f}  (dz/lambda_D={dz/lamD:.1f})")
    # (1,0) spectra
    ax = fig.add_subplot(gs[1, 0])
    for lbl, (k, P, c, *_ ) in spec.items():
        ax.loglog(1 / k, P / P.max(), color=c, lw=1.3, label=lbl)
    ax.axvline(2 * dz, color="k", ls="--", lw=1, label=f"grid Nyquist ({2*dz:.1f} d_e)")
    ax.axvline(1.0, color="gray", ls=":", lw=1, label="d_e")
    ax.set_xlabel("wavelength [d_e]"); ax.set_ylabel("norm power |Bx_k|^2")
    ax.set_xlim(0.4, 300); ax.legend(fontsize=7)
    ax.set_title("spectrum: peak at plasma scale=physical; rising to grid=numerical")
    # (1,1) hodogram of far-upstream + foot
    ax = fig.add_subplot(gs[1, 1])
    for lbl in ("foot (reflected ions)", "far upstream"):
        if lbl in spec:
            _, _, c, bxs, bys = spec[lbl]
            cc = np.corrcoef(bxs, bys)[0, 1]
            ax.plot(bxs - bxs.mean(), bys - bys.mean(), lw=0.35, color=c,
                    label=f"{lbl}: corr={cc:+.2f}")
    ax.set_xlabel("dBx/B0"); ax.set_ylabel("dBy/B0"); ax.set_aspect("equal")
    ax.legend(fontsize=8); ax.set_title("polarization (ellipse=wave, blob=noise)")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
