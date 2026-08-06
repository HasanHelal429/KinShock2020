#!/usr/bin/env python
"""Side-by-side ambient-ion phase space for two runs at matched t*wci0.

    python scripts/compare_phase.py runs/R1_phase/R1_paper runs/R1_phase/R1_paper_470eV

WHY. R1_paper and R1_paper_470eV are the SAME dimensionless setup -- identical n_amb, dz,
dt, rho_i0/dz, beta_0 ~ 0.2, M_ms ~ 12.75 -- differing only in that every velocity is 10x
smaller and every temperature 100x smaller in the 470 eV run, with dz and dt held FIXED.
That single change takes dz/lambda_D,amb from 0.60 to 6.07. R1_paper passes through the
stages of shock formation; the 470 eV run appears to accelerate ambient ions immediately
instead of sweeping them up first. This plots the two together so the difference is visible
rather than inferred from scalars.

NORMALISATION. Velocities are in v_sh (0.1430c vs 0.0165c, both by-eye fits) and positions
in d_i0 (identical in metres). Absolute units would put the two runs on axes 10x apart and
show nothing. Both rows share limits, so a feature at the same place means the same thing.

READ IT AS: how far ahead of the piston does the accelerated ambient population reach? A
forming shock keeps its reflected ions in a foot ~1-2 d_i0 deep; ions running many d_i0
upstream are not a foot.
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
from matplotlib.colors import LogNorm                 # noqa: E402

import kinshock                                       # noqa: E402
from kinshock import io as kio, metrics, plotting as P  # noqa: E402
from kinshock.units import C, ME                      # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BASELINE = "#c3c2b7"


def style(ax):
    ax.set_facecolor("#111014")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK2)


def load(d, targets, tmax):
    cfg = kinshock.load(d)
    sc = kinshock.units.derive(cfg)
    fit = metrics.load_shock_fit(d, sc)
    vsh = fit.v_sh if fit is not None else sc.vsh_model
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(d)
    times = []
    for p in pfs:
        times.append(kio.load_frame(p).time * sc.wci0)
        if times[-1] > tmax:
            break
    times = np.array(times)
    out = []
    for tt in targets:
        i = int(np.argmin(np.abs(times - tt)))
        fr = kio.load_frame(pfs[i])
        ad = fr.ds.all_data()
        z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
        p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
        w = np.asarray(ad[("amb_ions", "particle_weight")])
        g = np.sqrt(1.0 + (p / (mi * C)) ** 2)
        v = p / (mi * g) / vsh
        zp = np.percentile(
            np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5) / sc.di0
        out.append((times[i], z, v, w, zp))
    return cfg["meta"]["run_id"], sc, vsh, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs=2)
    ap.add_argument("--at", type=float, nargs="+", default=[0.11, 0.45, 1.01])
    ap.add_argument("--zmax", type=float, default=30.0, help="z/d_i0 axis limit")
    ap.add_argument("--vmax", type=float, default=2.5, help="v_z/v_sh axis limit")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    tmax = max(args.at) * 1.2
    A = load(args.run_dirs[0], args.at, tmax)
    B = load(args.run_dirs[1], args.at, tmax)
    n = len(args.at)

    fig, axes = plt.subplots(2, n, figsize=(4.3 * n, 7.2), dpi=200, facecolor=SURFACE,
                            sharex=True, sharey=True)
    zb = np.linspace(0, args.zmax, 320)
    vb = np.linspace(-args.vmax * 0.4, args.vmax, 260)

    for r, (rid, sc, vsh, frames) in enumerate((A, B)):
        for c, (t, z, v, w, zp) in enumerate(frames):
            ax = axes[r, c]
            style(ax)
            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            H = np.ma.masked_where(H <= 0, H)
            ax.pcolormesh(zb, vb, H.T, norm=LogNorm(), cmap="magma", shading="auto")
            ax.axhline(0, color="#ffffff", lw=0.5, alpha=0.25)
            ax.axhline(1.0, color="#5ad2ff", lw=0.9, ls=(0, (5, 4)), alpha=0.9)
            ax.axvline(zp, color="#ffd166", lw=1.2)
            # how far the accelerated population runs ahead of the piston
            acc = v > 0.3
            lead = (np.percentile(z[acc], 99.5) - zp) if acc.sum() else np.nan
            ax.annotate(f"lead {lead:.1f} $d_{{i0}}$", xy=(0.97, 0.95),
                        xycoords="axes fraction", ha="right", va="top",
                        fontsize=9, color="#ffd166")
            if r == 0:
                ax.set_title(rf"$t\,\omega_{{ci0}}$ = {t:.2f}", color=INK, fontsize=11)
            if c == 0:
                ax.set_ylabel(rf"{rid}" "\n" r"$v_z / v_{sh}$", color=INK2, fontsize=9.5)
            if r == 1:
                ax.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=10)
            ax.set_xlim(0, args.zmax)
            ax.set_ylim(vb[0], vb[-1])

    fig.suptitle("ambient ions — yellow = piston front, dashed = $v_{sh}$",
                 color=INK, fontsize=12.5)
    fig.text(0.5, 0.012,
             "Same dimensionless setup: identical n_amb, dz, dt, rho_i0/dz, beta_0, M_ms. "
             "Only the velocity/temperature scale differs (10x/100x), which takes "
             "dz/lambda_D from 0.60 to 6.07.",
             ha="center", fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.925, bottom=0.10,
                        wspace=0.06, hspace=0.06)
    out = args.out or os.path.join(P.ROOT, "media", "compare_phase.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
