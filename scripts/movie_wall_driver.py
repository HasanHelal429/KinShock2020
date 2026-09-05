#!/usr/bin/env python3
"""Phase-space + profile movie for a REFLECTING-WALL run (W_phase).

Separate from `movie_three_run_compare.py` on purpose: that script keys every position off
`z_piston`, the 99.5th percentile of `piston_ions`, and a wall-driven run HAS NO PISTON
SPECIES. There is nothing to take a percentile of. The wall at z = 0 is the origin instead,
which is a better one anyway -- it is exact, fixed, and free of the estimator ambiguity
that made the notch metric misread `ss_dz4_ppc400`.

WHAT IT SHOWS, top to bottom, sharing one z axis:
  1. ambient-ion phase space v_z/v_sh -- the C and the notch on its lower branch
  2. n/n_0 and |Bx|/B0 -- the compression and field pileup, so a real ramp is visible as
     both at once rather than inferred from the ions alone
The wall is drawn at z = 0 and v_z = 0 and -v_sh are marked: the inflow arrives at -v_sh
and specular reflection puts the reflected population at +v_sh, so a shock forming shows up
as material appearing above the v_z = 0 line and a compression peak walking to the right.

|Bx|/B0 is smoothed with a 0.2 d_i0 boxcar. RESULTS 2026-08-11 requires it -- a raw max runs
2.2-2.8x the smoothed value and that trap has already produced one retracted result -- and
a raw line also paints the panel solid when the cell count is high.

Linear colour scale, fixed across all frames (the run's own max), for the reason the notch
work established: the notch is a step in the DENSE lower branch, and a log stretch lifts the
sparse tails until they compete with it, while a per-frame scale makes something that grows
look constant.

    python scripts/movie_wall_driver.py runs/W_phase/ws_smoke_470eV [--fps 8]
"""
import argparse, os, shutil, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.colors import Normalize, LogNorm     # noqa: E402
import kinshock                                      # noqa: E402
from kinshock import io as kio, plotting as P        # noqa: E402
from kinshock.units import C, ME                     # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
WALL, GOOD, DENS, FIELD = "#ffd166", "#5ad2ff", "#ff6b57", "#e8e6ef"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--norm", choices=("linear", "log"), default="linear")
    ap.add_argument("--zmax", type=float, default=None, help="z limit in d_i0")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rid = os.path.basename(args.run.rstrip("/"))
    out = args.out or f"media/W_phase/{rid}_wall.mp4"

    import yt; yt.set_log_level(50)
    cfg = kinshock.load(args.run)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(args.run)
    zmax = args.zmax or sc.domain_halfwidth / sc.di0
    nz, nv = 320, 230
    zb = np.linspace(0, zmax, nz + 1)
    vb = np.linspace(-1.8, 1.8, nv + 1)

    # one load pass, so the colour scale can be fixed across the whole movie
    frames = []
    for p in pfs:
        fr = kio.load_frame(p, fields=("Bx",))
        ad = fr.ds.all_data()
        z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
        pm = np.asarray(ad[("amb_ions", "particle_momentum_z")])
        g = np.sqrt(1 + (pm / (mi * C)) ** 2)
        v = pm / (mi * g) / sc.vsh_model
        w = np.asarray(ad[("amb_ions", "particle_weight")])
        H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
        n, _ = np.histogram(z, bins=zb, weights=w)
        far = np.median(n[-40:])
        Bx = fr.comps.get("Bx")
        zf = fr.z_centers / sc.di0
        if Bx is not None:
            kw = max(int(round(0.2 / float(np.median(np.diff(zf))))), 1)
            b = np.convolve(np.abs(Bx) / sc.B0, np.ones(kw) / kw, mode="same")
        else:
            zf, b = np.array([]), np.array([])
        frames.append(dict(t=fr.time * sc.wci0, H=H, n=n / far if far > 0 else n,
                           zf=zf, b=b))
    print(f"  loaded {len(frames)} frames", flush=True)
    vmax = max(f["H"].max() for f in frames)
    nmax = max(3.0, min(np.nanmax([f["n"].max() for f in frames]) * 1.1, 12))
    bmax = max(3.0, min(np.nanmax([f["b"].max() for f in frames if f["b"].size]) * 1.1, 12))

    scratch = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))
    fdir = os.path.join(scratch, f"wall_{rid}")
    shutil.rmtree(fdir, ignore_errors=True)
    os.makedirs(fdir, exist_ok=True)
    zc = 0.5 * (zb[1:] + zb[:-1])
    for k, f in enumerate(frames):
        # even pixel dims: libx264 + yuv420p refuses odd ones, and only says so at encode
        # time -- after the whole load pass has run (cost one encode on M_phase).
        fig, ax = plt.subplots(2, 1, figsize=(9.2, 7.04), dpi=125, facecolor=SURFACE,
                               sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
        a = ax[0]; a.set_facecolor("#111014")
        nrm = (Normalize(0, vmax) if args.norm == "linear"
               else LogNorm(vmin=max(vmax * 1e-4, 1e-30), vmax=vmax))
        a.pcolormesh(zc, 0.5 * (vb[1:] + vb[:-1]),
                     np.ma.masked_where(f["H"] <= 0, f["H"]).T, norm=nrm,
                     cmap="magma", shading="auto", zorder=0)
        a.axhline(0, color="#fff", lw=.5, alpha=.25, zorder=2)
        a.axhline(-1, color=GOOD, lw=.9, ls=(0, (5, 4)), alpha=.8, zorder=2)
        a.axhline(1, color=GOOD, lw=.9, ls=(0, (5, 4)), alpha=.8, zorder=2)
        a.annotate("inflow  $-v_{sh}$", xy=(0.985, 0.10), xycoords="axes fraction",
                   ha="right", color=GOOD, fontsize=9)
        a.annotate("reflected  $+v_{sh}$", xy=(0.985, 0.86), xycoords="axes fraction",
                   ha="right", color=GOOD, fontsize=9)
        a.set_ylabel(r"$v_z/v_{sh}$   (ambient ions)", color=INK2, fontsize=10)
        a.annotate(rf"$t\,\omega_{{ci0}}$ = {f['t']:.3f}", xy=(0.015, 0.90),
                   xycoords="axes fraction", color=INK, fontsize=11, family="monospace")
        b = ax[1]; b.set_facecolor("#111014")
        b.plot(zc, f["n"], color=DENS, lw=1.4, label=r"$n/n_0$")
        if f["b"].size:
            b.plot(f["zf"], f["b"], color=FIELD, lw=1.3, label=r"$|B_x|/B_0$")
        b.axhline(1, color=MUTED, lw=.7, ls=(0, (4, 4)))
        b.set_ylim(0, max(nmax, bmax)); b.set_ylabel("compression", color=INK2, fontsize=10)
        b.set_xlabel(r"$z / d_{i0}$        (wall at $z=0$)", color=INK2, fontsize=11)
        lg = b.legend(loc="upper right", fontsize=8.5, framealpha=.25,
                      facecolor="#111014", edgecolor=BASE)
        for tx in lg.get_texts():
            tx.set_color(INK2)
        for axx in ax:
            axx.axvline(0, color=WALL, lw=2.6, zorder=4)
            axx.set_xlim(0, zmax)
            axx.tick_params(colors=MUTED, labelsize=8.5)
            for s in axx.spines.values():
                s.set_color(BASE)
        fig.suptitle(f"{rid} — reflecting-wall driver, no piston species", color=INK,
                     fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(fdir, f"frame_{k:03d}.png"), facecolor=SURFACE)
        plt.close(fig)
    print(f"  rendered {len(frames)} frames", flush=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(P.encode(fdir, out, fps=args.fps))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
