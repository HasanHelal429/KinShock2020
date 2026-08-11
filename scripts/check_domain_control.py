#!/usr/bin/env python
"""Does a truncated-domain run reproduce the full one over the window they share?

    python scripts/check_domain_control.py runs/S_phase/ss_dz1_ppc100 \
        runs/R1_phase/R1_paper_470eV --tmax 0.30

WHY THIS EXISTS. runs/S_phase shrinks R1_paper_470eV from 80.5 d_i0 to 12.0 d_i0 so a
resolution sweep of the early shock costs minutes instead of days. That is not free:
`field_hi = open` is pec, which REFLECTS fields (only particles are absorbed; Silver-Mueller
is unavailable because the B-field divergence cleaner runs whenever a background B is set,
see kinshock.deck._BC_MAP). Light crosses 12 d_i0 in ~6000 steps, so the piston turn-on
precursor makes ~25 round trips inside the S_phase window against ~4 in the full domain.
Until that is shown not to matter, every S_phase result is uninterpretable.

WHY NOT plot_ez.py's UPSTREAM NUMBER. It defines "far upstream" as the outer 5% of the
domain, which is z = 11.4-12.0 d_i0 in the truncated run and 76.5-80.5 d_i0 in the full
one -- different physical plasma, so the comparison would be meaningless. This script
fixes an ABSOLUTE z-window in d_i0 that is upstream in BOTH runs (default 8-12 d_i0: the
model shock only reaches 4.19 d_i0 by t*wci0 = 0.30) and measures both runs there.

WHAT IS COMPARED, and why these two:
  * rms E_z in that window, in v_A*B0. This is the quantity the whole 470 eV diagnosis
    turns on -- 30.2 against R1_paper's 2.25 at t*wci0 = 0.02 (RESULTS 2026-08-05). If the
    truncated domain inflates it, the sweep would be measuring its own boundary.
  * coherent max |B_perp|/B0, smoothed over 0.2 d_i0. A shock BUILDS its barrier, so this
    is the shape that says whether the same thing is forming at the same time.

BOUNDARY EXCLUSION IS NOT OPTIONAL for the B_perp number. An earlier "11x instant magnetic
barrier" was retracted because every maximum sat at z ~ 0.10 d_i0 -- the foil wall, not the
shock (RESULTS 2026-08-05). The default cuts 0.5 d_i0 from each end, and at 12 d_i0 that
matters far more than it did at 80.5.

A PASS is agreement inside the run-to-run scatter, NOT an exact match: the two runs differ
in domain AND in RNG stream, so some difference is expected. Read the traces, not a single
ratio.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402

import kinshock                                       # noqa: E402
from kinshock import io as kio, plotting as P         # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BASELINE = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834"]


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=9, length=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK2)


def boxcar(a, w):
    return a if w < 2 else np.convolve(a, np.ones(w) / w, mode="same")


def measure(run_dir, tmax, zlo, zhi, edge, smooth_di0):
    """(run_id, t*wci0[], rms Ez in [zlo,zhi] / (vA*B0)[], coherent max|Bperp|/B0[])."""
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    E_scale = sc.vA * sc.B0
    w = max(1, int(round(smooth_di0 * sc.di0 / sc.dz)))
    ne = max(1, int(round(edge * sc.di0 / sc.dz)))
    t, ez, bp = [], [], []
    for p in kio.field_plotfiles(run_dir):
        # Bx/By are first-class Frame attributes; Ez is not -- it only appears in
        # Frame.comps when requested, and comp_opt drops it silently if the plotfile
        # lacks it. Hence the explicit None check rather than a .get() default: a
        # missing Ez must stop the comparison, not be laundered into zeros.
        fr = kio.load_frame(p, fields=("Ez",))
        tt = fr.time * sc.wci0
        if tt > tmax:
            break
        E = fr.comps.get("Ez")
        if E is None:
            raise SystemExit(f"{run_dir}: {os.path.basename(p)} carries no Ez — is this "
                             f"a diag_fields plotfile?")
        E = np.asarray(E)
        zc = (np.arange(len(E)) + 0.5) * sc.dz / sc.di0
        m = (zc >= zlo) & (zc <= zhi)
        if not m.any():                    # window outside this run's domain -> say so
            raise SystemExit(f"{run_dir}: z-window {zlo}-{zhi} d_i0 is outside the "
                             f"domain (0-{zc[-1]:.1f} d_i0); pass a --zwindow inside both")
        t.append(tt)
        ez.append(float(np.sqrt(np.mean((E[m] / E_scale) ** 2))))
        sm = boxcar(np.hypot(np.asarray(fr.Bx), np.asarray(fr.By)), w)[ne:-ne]
        bp.append(float(np.max(sm)) / sc.B0)
    return cfg["meta"]["run_id"], np.array(t), np.array(ez), np.array(bp), sc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("truncated")
    ap.add_argument("reference")
    ap.add_argument("--tmax", type=float, default=0.30, help="t*wci0 limit (default 0.30)")
    ap.add_argument("--zwindow", type=float, nargs=2, default=(8.0, 12.0),
                    metavar=("LO", "HI"),
                    help="absolute upstream window in d_i0, must be inside BOTH domains")
    ap.add_argument("--edge", type=float, default=0.5,
                    help="d_i0 cut from each end before the B_perp max (default 0.5)")
    ap.add_argument("--smooth-di0", type=float, default=0.2)
    ap.add_argument("--at", type=float, nargs="*", default=[0.02, 0.10, 0.20, 0.30])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    runs = [measure(d, args.tmax, *args.zwindow, args.edge, args.smooth_di0)
            for d in (args.truncated, args.reference)]
    if any(len(r[1]) == 0 for r in runs):
        print("one of the runs has no field frames in the window yet")
        return 1

    print(f"\nupstream window {args.zwindow[0]:g}-{args.zwindow[1]:g} d_i0, "
          f"t*wci0 <= {args.tmax:g}, B_perp edge cut {args.edge:g} d_i0\n")
    lbl = [r[0] for r in runs]
    print(f"{'t*wci0':>8} | {'rms Ez / vA B0':^{2*13+3}} | {'coherent |Bperp|/B0':^{2*13+3}}")
    print(f"{'':>8} | {lbl[0][:13]:>13} {lbl[1][:13]:>13} {'ratio':>7} |"
          f" {lbl[0][:13]:>13} {lbl[1][:13]:>13} {'ratio':>7}")
    print("-" * 96)
    for tt in args.at:
        cells = []
        for col in (2, 3):
            vals = []
            for _rid, t, ez, bp, _sc in runs:
                arr = ez if col == 2 else bp
                if tt > t.max() * 1.15:
                    vals.append(float("nan"))
                    continue
                vals.append(float(arr[int(np.argmin(np.abs(t - tt)))]))
            cells.append(vals + [vals[0] / vals[1] if vals[1] else float("nan")])
        print(f"{tt:8.2f} | {cells[0][0]:13.3g} {cells[0][1]:13.3g} {cells[0][2]:7.2f} |"
              f" {cells[1][0]:13.3g} {cells[1][1]:13.3g} {cells[1][2]:7.2f}")
    print()
    for rid, t, ez, bp, sc in runs:
        print(f"  {rid}: {len(t)} field frames, t*wci0 {t.min():.4f}..{t.max():.4f}, "
              f"domain {sc.domain_halfwidth/sc.di0:.1f} d_i0, v_A B0 = {sc.vA*sc.B0:.4g} V/m")

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), dpi=200, facecolor=SURFACE)
    for ax in axes:
        style(ax)
    for i, (rid, t, ez, bp, sc) in enumerate(runs):
        tag = f"{rid}  ({sc.domain_halfwidth/sc.di0:.0f} $d_{{i0}}$)"
        axes[0].plot(t, ez, "-o", color=SERIES[i], lw=1.7, ms=3, label=tag)
        axes[1].plot(t, bp, "-o", color=SERIES[i], lw=1.7, ms=3, label=tag)
    axes[0].set_ylabel(rf"rms $E_z$ / $v_A B_0$  ({args.zwindow[0]:g}–{args.zwindow[1]:g} $d_{{i0}}$)",
                       color=INK2, fontsize=10)
    axes[1].set_ylabel(r"coherent max $|B_\perp| / B_0$", color=INK2, fontsize=10)
    axes[1].axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 4)))
    for ax, ttl in zip(axes, ("upstream longitudinal field", "magnetic pileup")):
        ax.set_xlabel(r"$t\,\omega_{ci0}$", color=INK2, fontsize=10)
        ax.set_title(ttl, color=INK, fontsize=11.5, loc="left", pad=10)
        leg = ax.legend(frameon=False, fontsize=8.5)
        for x in leg.get_texts():
            x.set_color(INK2)
    fig.text(0.01, 0.015,
             "Truncated vs full domain over the window they share. Traces that separate "
             "mean the pec far boundary is shaping the result, and the S_phase domain has "
             "to grow before its sweep can be read.",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.165, wspace=0.24)
    out = args.out or os.path.join(P.ROOT, "media", "S_phase", "domain_control.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
