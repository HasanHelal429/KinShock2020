#!/usr/bin/env python
"""Coherent magnetic pileup |B_perp|/B0 vs time, for one or more runs side by side.

    python scripts/plot_bperp_pileup.py runs/R1_phase/R1_paper runs/implicit_phase/i0_...

WHAT THIS MEASURES AND WHY IT EXISTS. A perpendicular shock forms by sweeping up ambient
plasma into a magnetic barrier: |B_perp| should START near B0 and GROW as the piston plows
through the ambient. R1_paper does exactly that (3.5 -> 4.9 -> 6.1 over t*wci0 = 0.14 ->
0.44). R1_paper_470eV instead sits at ~11x from the FIRST output frame and stays flat --
the barrier is already there before any sweeping could have built it, and the run
accelerates ambient ions immediately instead of first sweeping them up (RESULTS 2026-08-05).
This script is the discriminating measurement between candidate causes.

WHY SMOOTHING IS NOT OPTIONAL. These runs are Debye-under-resolved (dz/lambda_D,amb = 6.07
in the 470 eV run), so the raw B_perp trace carries grid noise comparable to the signal.
The bare max over the domain is therefore a noise statistic, not a pileup statistic -- it
was the first thing I checked and discarded. Smoothing over 0.2 d_i0 keeps the barrier
(which is ~1 d_i0 wide) and averages the noise down.

The NOISE RATIO column is the honesty check: max|raw| / max|smoothed|. If a run's apparent
pileup were merely noise, that ratio would be LARGE where the pileup looks large. On the
two reference runs it goes the other way (1.47-1.67 for the 11x run against 1.66-2.54 for
the physical one), which is what rules out "the 11x is a noise artifact".
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
# categorical, colour-blind safe; first two are the validated diverging pair's hues
SERIES = ["#2a78d6", "#eb6834", "#3fa34d", "#8e5ea2", "#c0392b"]


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=1.0)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK2)


def boxcar(a, w):
    return a if w < 2 else np.convolve(a, np.ones(w) / w, mode="same")


def measure(run_dir, smooth_di0):
    """(label, t*wci0[], coherent max/B0[], noise ratio[]) for every field frame."""
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    w = max(1, int(round(smooth_di0 * sc.di0 / sc.dz)))
    t, coh, ratio = [], [], []
    for p in kio.field_plotfiles(run_dir):
        fr = kio.load_frame(p)
        bp = np.hypot(fr.Bx, fr.By)            # perpendicular to the shock normal z
        sm = boxcar(bp, w)
        mx_s = float(np.max(sm))
        t.append(fr.time * sc.wci0)
        coh.append(mx_s / sc.B0)
        ratio.append(float(np.max(bp)) / mx_s if mx_s > 0 else np.nan)
    return cfg["meta"]["run_id"], np.array(t), np.array(coh), np.array(ratio), w


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--smooth-di0", type=float, default=0.2,
                    help="boxcar width in d_i0 (default 0.2)")
    ap.add_argument("--at", type=float, nargs="*", default=[0.14, 0.29, 0.44],
                    help="t*wci0 values to tabulate (nearest frame)")
    ap.add_argument("--out", default=None, help="output png (default media/bperp_pileup.png)")
    args = ap.parse_args()

    runs = []
    for d in args.run_dirs:
        try:
            runs.append(measure(d, args.smooth_di0))
        except Exception as e:                  # a missing run must not look like a flat line
            print(f"  SKIP {d}: {type(e).__name__}: {e}")
    if not runs:
        return 1

    print(f"\ncoherent |B_perp|/B0, smoothed over {args.smooth_di0:g} d_i0\n")
    hdr = f"{'t*wci0':>8}" + "".join(f"{r[0][:22]:>24}" for r in runs)
    print(hdr)
    print("-" * len(hdr))
    for target in args.at:
        row = f"{target:8.2f}"
        for _, t, coh, ratio, _w in runs:
            if len(t) == 0 or target > t.max() * 1.15:
                row += f"{'--':>24}"           # do not extrapolate past what ran
                continue
            i = int(np.argmin(np.abs(t - target)))
            row += f"{coh[i]:>13.2f} (n {ratio[i]:.2f})"
        print(row)
    print()
    for rid, t, coh, ratio, w in runs:
        print(f"  {rid}: {len(t)} frames, t*wci0 {t.min():.3f}..{t.max():.3f}, "
              f"smoothing {w} cells, coherent {coh.min():.2f}..{coh.max():.2f}")

    fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=200, facecolor=SURFACE)
    style(ax)
    for i, (rid, t, coh, _r, _w) in enumerate(runs):
        ax.plot(t, coh, "-o", color=SERIES[i % len(SERIES)], lw=1.8, ms=3.5, label=rid)
    ax.axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 4)))
    ax.text(ax.get_xlim()[1], 1.0, " B$_0$ ", va="center", fontsize=8, color=MUTED)
    ax.set_xlabel(r"$t\,\omega_{ci0}$", color=INK2, fontsize=10)
    ax.set_ylabel(r"max coherent $|B_\perp| / B_0$", color=INK2, fontsize=10)
    ax.set_title(rf"magnetic pileup, smoothed over {args.smooth_di0:g} $d_{{i0}}$",
                 color=INK, fontsize=12, loc="left", pad=10)
    leg = ax.legend(frameon=False, fontsize=9)
    for x in leg.get_texts():
        x.set_color(INK2)
    fig.text(0.01, 0.015,
             "A shock BUILDS its barrier: |B_perp| should start near B0 and grow. A run that "
             "starts high and stays flat never swept anything up.",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.10, right=0.97, top=0.91, bottom=0.145)
    out = args.out or os.path.join(P.ROOT, "media", "bperp_pileup.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
