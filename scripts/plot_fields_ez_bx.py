#!/usr/bin/env python3
"""E_z and B_x on SHARED spacetime axes, so the two can be read against each other.

    python scripts/plot_fields_ez_bx.py runs/W_phase/ws_470eV_t14

WHICH COMPONENTS, AND WHY THESE ARE THE "PARALLEL E" AND "PERPENDICULAR B". B0 lies along x
and the shock propagates along z, so components are named against the shock NORMAL:

  * E_z is PARALLEL to the normal -- the longitudinal, essentially electrostatic field that
    carries the cross-shock potential. (E_y = -v_z B_x is the motional field and is not the
    interesting one.) It is what accelerates ambient ions along z, i.e. what builds the C.
  * B_x is PERPENDICULAR to the normal -- the swept-up magnetic barrier, the field that
    piles up as the driver plows through the ambient.

`plot_ez.py` already draws E_z alone. This exists because "is the field the root cause?" is
a question about CO-LOCATION: whether an E_z feature sits at the same z and t as a B_x
feature, or leads it, or has no counterpart. That cannot be read off two figures with
independent axes, so both panels here share one z and one t axis and are drawn from the same
frames.

B_x IS PLOTTED SIGNED, NOT AS |B_x|. Taking the absolute value first rectifies upstream
fluctuation into a DC offset that any subsequent smoothing preserves -- measured 2026-09-04,
upstream <|Bx|>/B0 reads 1.67 and 10.92 on two runs where the truth is 1.00. Signed avoids
the trap entirely and shows field REVERSALS, which |.| would fold away.

Smoothing differs per field on purpose: 0.2 d_i0 on B_x, which RESULTS 2026-08-11 requires
(a raw max runs 2.2-2.8x the smoothed value), and 1.0 d_i0 on E_z, following plot_ez.py --
these runs are Debye-under-resolved and a raw E_z streak is hash, not structure. The
line-out panels draw the RAW trace behind the smoothed one so the noise stays on the page.

Colour is diverging with a neutral grey midpoint, symmetric about 0 for E_z and about 1 for
B_x/B0, so grey means "no longitudinal field" and "still at B0" respectively.
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib                                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
import kinshock                                        # noqa: E402
from kinshock import io as kio                         # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
RAW = "#c9c7bf"
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
        ax.spines[s].set_color(BASELINE); ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=1.0)


def boxcar(a, w):
    return a if w < 2 else np.convolve(a, np.ones(w) / w, mode="same")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--ez-smooth-di0", type=float, default=1.0)
    ap.add_argument("--bx-smooth-di0", type=float, default=0.2)
    ap.add_argument("--zmax", type=float, default=None)
    ap.add_argument("--times", type=float, nargs="*", default=None,
                    help="t*wci0 values for the line-outs (default: 5 spread over the run)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import yt; yt.set_log_level(50)
    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    rid = cfg["meta"]["run_id"]
    pfs = kio.field_plotfiles(args.run_dir)[:: args.stride]
    if not pfs:
        print("no field plotfiles"); return 1

    zc = np.asarray(kio.load_frame(pfs[0]).z_centers) / sc.di0
    dz = float(np.median(np.diff(zc)))
    wE = max(int(round(args.ez_smooth_di0 / dz)), 1)
    wB = max(int(round(args.bx_smooth_di0 / dz)), 1)
    Escale = sc.vsh_model * sc.B0            # motional scale v_sh*B0

    t, Ez, Bx, Ezr, Bxr = [], [], [], [], []
    for p in pfs:
        fr = kio.load_frame(p, fields=("Ez", "Bx"))
        e = fr.comps.get("Ez"); b = fr.comps.get("Bx")
        if e is None or b is None:
            continue
        t.append(fr.time * sc.wci0)
        Ezr.append(e / Escale); Bxr.append(b / sc.B0)
        Ez.append(boxcar(e / Escale, wE)); Bx.append(boxcar(b / sc.B0, wB))
    t = np.asarray(t); Ez = np.asarray(Ez); Bx = np.asarray(Bx)
    Ezr = np.asarray(Ezr); Bxr = np.asarray(Bxr)
    print(f"{rid}: {len(t)} field frames, t* {t[0]:.3f}..{t[-1]:.3f}, "
          f"boxcar E {wE} cells / B {wB} cells", flush=True)

    zmax = args.zmax or zc.max()
    m = zc <= zmax
    times = args.times or list(np.linspace(t[0], t[-1], 5))

    fig, ax = plt.subplots(2, 2, figsize=(13.6, 8.4), dpi=130, facecolor=SURFACE,
                           gridspec_kw=dict(width_ratios=[1.35, 1]))
    for a in ax.ravel():
        style(a)

    for row, (F, Fr, lab, unit, ctr, sm) in enumerate([
            (Ez, Ezr, r"$E_z$", r"$E_z\,/\,(v_{sh}B_0)$", 0.0, args.ez_smooth_di0),
            (Bx, Bxr, r"$B_x$", r"$B_x\,/\,B_0$",         1.0, args.bx_smooth_di0)]):
        lim = float(np.nanpercentile(np.abs(F[:, m] - ctr), 99.0)) or 1.0
        a = ax[row][0]
        pc = a.pcolormesh(t, zc[m], F[:, m].T, shading="auto", cmap=DIVERGING,
                          vmin=ctr - lim, vmax=ctr + lim)
        cb = fig.colorbar(pc, ax=a, pad=0.015, fraction=0.046)
        cb.set_label(unit, color=INK2, fontsize=9)
        cb.ax.tick_params(colors=MUTED, labelsize=8)
        cb.outline.set_edgecolor(BASELINE)
        a.set_xlabel(r"$t\,\omega_{ci0}$", color=INK2, fontsize=10)
        a.set_ylabel(r"$z / d_{i0}$", color=INK2, fontsize=10)
        a.set_title(f"{lab}  —  {sm} $d_{{i0}}$ boxcar", color=INK, fontsize=11, loc="left")
        for T in times:
            a.axvline(T, color=MUTED, lw=0.6, ls=(0, (3, 3)), alpha=.7)

        a2 = ax[row][1]
        # Offsets and y-limits come from the SMOOTHED traces. Scaling to the raw ones
        # instead lets E_z's near-wall spikes (|E_z| ~ 10x the coherent field) set the range
        # and squash all five line-outs into flat lines -- which is what the first version
        # of this figure did. The raw trace is still drawn, and simply clips.
        span = float(np.nanpercentile(np.abs(F[:, m] - ctr), 99.5)) or 1.0
        off = 2.6 * span
        for k, T in enumerate(times):
            i = int(np.argmin(np.abs(t - T)))
            a2.plot(zc[m], np.clip(Fr[i][m], ctr - 1.4 * span, ctr + 1.4 * span) + k * off,
                    color=RAW, lw=0.5, zorder=1)
            a2.plot(zc[m], F[i][m] + k * off, color="#1a1a1a", lw=1.2, zorder=2)
            a2.axhline(ctr + k * off, color=GRID, lw=0.7, zorder=0)
            a2.annotate(rf"$t\omega_{{ci0}}$={t[i]:.2f}", xy=(0.995, ctr + k * off + 0.42 * off),
                        xycoords=("axes fraction", "data"), ha="right",
                        color=MUTED, fontsize=8)
        a2.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=10)
        a2.set_ylabel(unit + "   (offset per time)", color=INK2, fontsize=10)
        a2.set_yticks([])
        a2.set_ylim(ctr - 1.7 * span, ctr + (len(times) - 1) * off + 1.7 * span)
        a2.set_title("line-outs (raw behind smoothed, raw clipped)",
                     color=INK, fontsize=11, loc="left")

    fig.suptitle(f"{rid} — longitudinal $E_z$ and perpendicular $B_x$ "
                 f"(B$_0\\parallel x$, normal $\\parallel z$)", color=INK, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = args.out or f"media/fields/{rid}_ez_bx.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
