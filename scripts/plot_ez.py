#!/usr/bin/env python
"""Plot E_z, the longitudinal electric field, for a KinShock2020 run.

    python scripts/plot_ez.py runs/R1_phase/R1_paper_470eV [--stride 4]

WHAT E_z IS HERE. B0 lies along x and the shock propagates along z, so the MOTIONAL field
is E_y = -v_z B_x and E_z is the longitudinal, essentially electrostatic field -- the one
that carries the cross-shock potential. It is the interesting component for shock structure
and the one most exposed to grid noise.

WHY THIS PLOTS TWO THINGS. In R1_paper_470eV the raw E_z has rms ~3.6e9 V/m against a
motional scale v_sh B0 ~ 3.5e7 V/m -- about 100x GLOBALLY. That global number is dominated
by the shocked plasma: the far upstream sits at only ~4-6x until the shock arrives, and it
is the upstream value that accompanies the grid heating. The figure reports both, because
quoting the global one alone would badly overstate the noise in the region the heating
measurement actually samples. A raw streak is therefore hash, and
showing only that would misrepresent the run as having no field structure; showing only a
smoothed version would hide that the noise dwarfs the signal. Both panels, always:

  * the streak is smoothed over ~1 d_i0 so the coherent ramp field is visible at all;
  * the line-outs draw the RAW trace behind the smoothed one, so the noise it was extracted
    from stays on the page.

That noise level is not incidental -- it is the same Debye under-resolution (dz/lambda_D,amb
= 6.07 at t=0) that drives the grid heating in RESULTS 2026-08-05. Reading the two together
is the point.

Colour: a diverging blue<->red map with a NEUTRAL GREY midpoint (#f0efec), symmetric about
zero so the midpoint is genuinely "no field". Diverging is required because E_z is signed;
a sequential ramp would put zero at one end and make sign changes invisible. Blue<->red
rather than blue<->aqua because two cool hues do not read as opposite, and it stays
distinguishable under red-green CVD.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

import kinshock                                       # noqa: E402
from kinshock import io as kio, metrics, plotting as P  # noqa: E402
from kinshock.units import C                          # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
BLUE, RED = "#2a78d6", "#c0392b"

# blue <-> grey <-> red, equal steps per arm (palette.md "Diverging pair")
DIVERGING = LinearSegmentedColormap.from_list("kin_div", [
    "#0d366b", "#184f95", "#256abf", "#3987e5", "#86b6ef", "#cde2fb",
    "#f0efec",
    "#fbd5d5", "#f2a0a0", "#e56b6b", "#d63f3f", "#ab2c2c", "#7a1d1d",
])


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
    if w < 2:
        return a
    k = np.ones(w) / w
    return np.convolve(a, k, mode="same")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--stride", type=int, default=4,
                    help="field-frame stride for the streak (default 4)")
    ap.add_argument("--smooth-di0", type=float, default=1.0,
                    help="boxcar width for the streak, in d_i0 (default 1.0)")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    fit = metrics.load_shock_fit(args.run_dir, sc)
    v_sh = fit.v_sh if fit is not None else sc.vsh_model
    E_scale = v_sh * sc.B0                        # motional field, the natural unit

    pfs = kio.field_plotfiles(args.run_dir)[::args.stride]
    t, rows, raws = [], [], []
    w = max(1, int(round(args.smooth_di0 * sc.di0 / sc.dz)))
    for p in pfs:
        fr = kio.load_frame(p, fields=("Ez",))
        ez = fr.comps.get("Ez")
        if ez is None:
            continue
        t.append(fr.time * sc.wci0)
        raws.append(ez)
        rows.append(boxcar(ez, w))
    t = np.array(t)
    Es = np.array(rows) / E_scale
    zc = np.asarray(kio.load_frame(pfs[0]).z_centers) / sc.di0
    R = np.array(raws) / E_scale
    rms_raw = float(np.sqrt(np.mean(R ** 2)))
    # The GLOBAL rms is dominated by the shocked plasma and says nothing about the far
    # upstream, where the grid heating was measured. Quote both: on R1_paper_470eV the
    # global value reaches ~157 while the outer 5% is only ~6 until the shock arrives.
    nz = R.shape[1]
    up = R[:, int(0.95 * nz):]
    rms_up_early = float(np.sqrt(np.mean(up[: max(1, len(up) * 3 // 4)] ** 2)))

    fig = plt.figure(figsize=(12.4, 5.2), dpi=200, facecolor=SURFACE)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.26)
    ax, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    style(ax)
    style(ax2)

    lim = float(np.nanpercentile(np.abs(Es), 99.5))    # symmetric -> grey sits at zero
    pc = ax.pcolormesh(t, zc, Es.T, shading="auto", cmap=DIVERGING, vmin=-lim, vmax=lim)
    if fit is not None or True:
        tt = np.linspace(t.min(), t.max(), 50)
        ax.plot(tt, (v_sh * (tt / sc.wci0)) / sc.di0, "-", color=INK, lw=1.3, alpha=0.8,
                label=rf"$v_{{sh}}$ = {v_sh/C:.4f}c" +
                      ("  (by-eye fit)" if fit is not None else "  (model)"))
        leg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
        for x in leg.get_texts():
            x.set_color(INK)
    ax.set_ylim(0, zc.max())
    ax.set_xlabel(r"$t\,\omega_{ci0}$", color=INK2, fontsize=10)
    ax.set_ylabel(r"$z / d_{i0}$", color=INK2, fontsize=10)
    ax.set_title(rf"$E_z$ smoothed over {args.smooth_di0:g} $d_{{i0}}$ "
                 rf"({w} cells) — the coherent field", color=INK, fontsize=11.5,
                 loc="left", pad=10)
    cb = fig.colorbar(pc, ax=ax, label=r"$E_z\,/\,(v_{sh}B_0)$")
    cb.outline.set_edgecolor(BASELINE)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.ax.yaxis.label.set_color(INK2)

    # ---- line-outs -------------------------------------------------------------------
    # The raw trace is ~20x the smoothed one, so overlaying it at full amplitude buries the
    # signal in a hairball and makes every smoothed curve look flat. Draw the noise as a
    # +/-1 sigma BAND around the smoothed line instead: the amplitude stays on the page,
    # the coherent field stays readable, and the offsets can be scaled to the band.
    picks = [int(f * (len(pfs) - 1)) for f in (0.25, 0.55, 0.85)]
    sig = []
    for i in picks:
        resid = raws[i] / E_scale - Es[i]
        sig.append(boxcar(resid ** 2, w) ** 0.5)
    span = 2.6 * float(np.nanpercentile(np.concatenate(sig), 90))
    off = 0.0
    for j, i in enumerate(picks):
        col = BLUE if j % 2 == 0 else RED
        ax2.fill_between(zc, Es[i] - sig[j] + off, Es[i] + sig[j] + off,
                         color=MUTED, alpha=0.28, lw=0, zorder=2,
                         label=r"$\pm1\sigma$ of the raw field" if j == 0 else None)
        ax2.plot(zc, Es[i] + off, "-", color=col, lw=1.7, zorder=3,
                 label=r"smoothed over 1 $d_{i0}$" if j == 0 else None)
        ax2.axhline(off, color=BASELINE, lw=0.8, ls=(0, (4, 4)), zorder=1)
        ax2.text(zc.max() * 0.985, off + span * 0.34,
                 rf"$t\,\omega_{{ci0}}$ = {t[i]:.2f}", ha="right", fontsize=8.5, color=INK)
        off -= span
    ax2.set_xlim(0, zc.max())
    ax2.set_ylim(off + span * 0.45, span * 0.62)
    ax2.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=10)
    ax2.set_ylabel(r"$E_z\,/\,(v_{sh}B_0)$   (offset per time)", color=INK2, fontsize=10)
    ax2.set_yticks([])
    ax2.set_title(r"line-outs: coherent field on the noise band", color=INK,
                  fontsize=11.5, loc="left", pad=10)
    leg2 = ax2.legend(frameon=False, fontsize=8.5, loc="lower left")
    for x in leg2.get_texts():
        x.set_color(INK2)

    fig.text(0.008, 0.105,
             f"{cfg['meta']['run_id']}    {len(pfs)} field frames (stride {args.stride})"
             rf"    $v_{{sh}}B_0$ = {E_scale:.3g} V/m    raw rms $E_z$ = {rms_raw:.0f}"
             rf" $v_{{sh}}B_0$ (ALL z)    far upstream only: {rms_up_early:.1f}"
             rf"    smoothing {w} cells = {args.smooth_di0:g} $d_{{i0}}$",
             fontsize=8, color=MUTED)
    fig.text(0.008, 0.018,
             f"E_z is Debye-noise dominated, but by very different amounts either side of "
             f"the front: ~{rms_raw:.0f}x the motional scale over the whole domain against "
             f"~{rms_up_early:.0f}x in the far upstream.\nThe global figure is the shocked "
             "plasma; it is the upstream one that accompanies the grid heating, and it "
             "roughly doubles over the run as T_0 rises (E ~ T/e.lambda_D ~ sqrt(T)).",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.055, right=0.975, top=0.90, bottom=0.235,
                        wspace=0.26)
    out = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "efield_ez.png")
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  v_sh B0 = {E_scale:.4g} V/m")
    print(f"  raw rms E_z: {rms_raw:.1f} v_sh B0 over all z, "
          f"{rms_up_early:.1f} in the far upstream (pre-arrival)")
    print(f"  smoothing {w} cells = {args.smooth_di0:g} d_i0   colour limit +-{lim:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
