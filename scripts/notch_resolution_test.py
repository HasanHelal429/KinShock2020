#!/usr/bin/env python3
"""Is the notch a function of resolution?  A visual ladder plus a number, at 470 eV.

THE NOTCH IS A STEP IN THE MEAN, NOT A SPREAD. It is the bump at the front of the C's
LOWER branch, where the ambipolar-accelerated ambient ions rejoin the undisturbed
upstream. Healthy runs roll over smoothly there; the 470 eV runs build a dense plateau at
v_z ~ 0.9 v_sh and drop off it vertically. Every earlier metric in this repo measured the
wedge -- an offset where the v_z STANDARD DEVIATION falls below a threshold -- which
integrates over the whole C and is structurally blind to a step on one branch. Hence a
separate diagnostic here:

    plateau    max of the weighted-median v_z profile in [piston, piston + 0.4 rho_i0],
               taken over the lower branch only (v_z < VCUT, which excludes the reflected
               population riding above it)
    drop       the largest NEGATIVE gradient of that profile, in v_sh per rho_i0 --
               the cliff steepness, which is the thing that looks like a shock
    z_drop     where that gradient sits, relative to the piston front

WHY THESE RUNS. S_phase is a 2-D grid in (dz/lambda_D, N_D = ppc*lambda_D/dz), built so
lines of constant N_D cross lines of constant dz -- refining dz improves BOTH aliasing and
particle noise, while ppc improves only noise, so a 1-D dz scan cannot say which matters.
Four points reach past the t* ~ 0.48 where the cliff appears:

    ss_dz1_ppc100_int1   dz/lD 6.07  N_D  16.5   baseline
    ss_dz1_ppc1600       dz/lD 6.07  N_D 264     ppc x16 at FIXED dz  -> particle noise
    ss_dz4_ppc400        dz/lD 1.52  N_D 264     dz /4 at FIXED N_D   -> aliasing
    R1_470eV_s3f8_t2     dz/lD 0.76  N_D 132     finest, + shape 3, filter 8

If the plateau and drop survive across an 8x span in dz and a 16x span in ppc, the notch
is not a resolution artifact. If they weaken as either is refined, it is.

    python scripts/notch_resolution_test.py
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                    # noqa: E402
from matplotlib.colors import Normalize            # noqa: E402
import kinshock                                    # noqa: E402
from kinshock import io as kio                     # noqa: E402
from kinshock.units import C, ME                   # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
ACC, LINE, BAD = "#ffd166", "#5ad2ff", "#ff6b57"
WLO, WHI = -0.15, 0.40      # rho_i0 about the piston front
VLO, VHI = -0.30, 1.45
VCUT = 1.30                 # exclude the reflected population above the lower branch
TOL = 0.05

RUNS = [("ss_dz1_ppc100_int1", "runs/S_phase/ss_dz1_ppc100_int1", 6.07, 16.5),
        ("ss_dz1_ppc1600",     "runs/S_phase/ss_dz1_ppc1600",     6.07, 264.0),
        ("ss_dz4_ppc400",      "runs/S_phase/ss_dz4_ppc400",      1.52, 264.0),
        ("R1_470eV_s3f8_t2",   "runs/R1_phase/R1_470eV_s3f8_t2",  0.76, 132.0)]


def profile(z, v, w, zp, nb=160):
    """Weighted-median v_z of the LOWER branch vs (z - z_piston)/rho_i0."""
    edges = np.linspace(WLO, WHI, nb + 1)
    rel = z - zp
    med = np.full(nb, np.nan)
    for k in range(nb):
        m = (rel >= edges[k]) & (rel < edges[k + 1]) & (v < VCUT)
        if m.sum() < 40:
            continue
        vv, ww = v[m], w[m]
        o = np.argsort(vv)
        cw = np.cumsum(ww[o]) / ww.sum()
        med[k] = vv[o][np.searchsorted(cw, 0.5)]
    return 0.5 * (edges[1:] + edges[:-1]), med


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--times", type=float, nargs="+", default=[0.40, 0.52])
    ap.add_argument("--out", default="media/comparison/notch_resolution.png")
    args = ap.parse_args()
    import yt; yt.set_log_level(50)

    fig, ax = plt.subplots(len(RUNS), len(args.times),
                           figsize=(6.2 * len(args.times), 2.9 * len(RUNS)),
                           dpi=160, facecolor=SURFACE, squeeze=False)
    print(f"{'run':22}{'dz/lD':>7}{'N_D':>7}{'t*':>7}{'plateau':>9}"
          f"{'drop':>9}{'z_drop':>8}")
    table = []
    for r, (lab, d, dzl, nd) in enumerate(RUNS):
        cfg = kinshock.load(d); sc = kinshock.units.derive(cfg)
        mi = cfg["reference"]["mass_ratio"] * ME
        pfs = kio.plotfiles(d)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        for c, T in enumerate(args.times):
            a = ax[r][c]; a.set_facecolor("#111014")
            i = int(np.argmin(np.abs(times - T)))
            if abs(times[i] - T) > TOL:
                a.annotate(f"no frame within {TOL} of {T}", xy=(.5, .5),
                           xycoords="axes fraction", ha="center", color=MUTED, fontsize=9)
                continue
            ad = kio.load_frame(pfs[i]).ds.all_data()
            z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.rho_i0
            p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
            g = np.sqrt(1 + (p / (mi * C)) ** 2)
            v = p / (mi * g) / sc.vsh_model
            w = np.asarray(ad[("amb_ions", "particle_weight")])
            zp = np.percentile(
                np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5) / sc.rho_i0
            zb = np.linspace(zp + WLO, zp + WHI, 260)
            vb = np.linspace(VLO, VHI, 200)
            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            nz = H[H > 0]
            a.pcolormesh(np.linspace(WLO, WHI, 259), np.linspace(VLO, VHI, 199),
                         np.ma.masked_where(H <= 0, H).T,
                         norm=Normalize(0, nz.max() if nz.size else 1),
                         cmap="magma", shading="auto", zorder=0)
            x, med = profile(z, v, w, zp)
            a.plot(x, med, color=LINE, lw=1.6, zorder=4)
            ok = np.isfinite(med)
            gx, gm = x[ok], med[ok]
            grad = np.gradient(gm, gx) if gm.size > 3 else np.array([np.nan])
            j = int(np.nanargmin(grad)) if np.isfinite(grad).any() else 0
            drop, zdrop = grad[j], gx[j]
            # Plateau is the TOP OF THE CLIFF -- the branch just upstream of the steepest
            # drop -- NOT "the max ahead of the piston front". The earlier definition
            # searched z > 0 only, and in ss_dz4_ppc400 the whole feature sits at
            # z - z_piston ~ -0.08, so it stepped over a plateau of 0.88 and reported the
            # 0.18 left behind it. The 99.5th-percentile piston front is a rough anchor,
            # not a reliable origin for the notch; locate the notch by its own gradient.
            lo = max(j - 25, 0)
            plateau = np.nanmax(gm[lo:j + 1]) if j > lo else (gm[j] if gm.size else np.nan)
            a.axvline(0, color=ACC, lw=1.2, zorder=3)
            a.axvline(zdrop, color=BAD, lw=1.1, ls=(0, (3, 3)), zorder=3)
            a.set_xlim(WLO, WHI); a.set_ylim(VLO, VHI)
            a.tick_params(colors=MUTED, labelsize=8)
            for s in a.spines.values():
                s.set_color(BASE)
            a.annotate(f"plateau {plateau:.2f}   drop {drop:.1f}", xy=(.975, .06),
                       xycoords="axes fraction", ha="right", color=LINE, fontsize=9)
            if r == 0:
                a.set_title(rf"$t\,\omega_{{ci0}}$ = {times[i]:.2f}", color=INK, fontsize=11)
            if c == 0:
                a.set_ylabel(f"{lab}\n" + rf"dz/$\lambda_D$={dzl}  $N_D$={nd:.0f}"
                             "\n" r"$v_z/v_{sh}$", color=INK2, fontsize=8.5)
            if r == len(RUNS) - 1:
                a.set_xlabel(r"$(z-z_{piston})/\rho_{i0}$", color=INK2, fontsize=10)
            print(f"{lab:22}{dzl:>7.2f}{nd:>7.0f}{times[i]:>7.3f}"
                  f"{plateau:>9.3f}{drop:>9.1f}{zdrop:>8.3f}", flush=True)
            table.append((lab, dzl, nd, times[i], plateau, drop))
    fig.suptitle("Is the notch a resolution artifact?  470 eV, 8x in dz, 16x in ppc\n"
                 "cyan = weighted-median $v_z$ of the lower branch; red dash = steepest drop",
                 color=INK, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"\n  wrote {args.out}")
    late = [t for t in table if t[3] > 0.45]
    if len(late) >= 2:
        pl = [t[4] for t in late]; dr = [t[5] for t in late]
        print(f"\n  across the ladder at t* ~ 0.5: plateau {min(pl):.2f}-{max(pl):.2f} "
              f"(spread {max(pl)/min(pl):.2f}x), steepest drop {min(dr):.0f} to {max(dr):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
