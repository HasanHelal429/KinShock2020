#!/usr/bin/env python3
"""Three-panel phase-space movie: the false shock, isolated.

    ek_47keV_t2        w_pe0/w_ce =  10   c/v_ti =  226   healthy
    mb_r100_t2         w_pe0/w_ce = 100   c/v_ti =  226   healthy
    R1_470eV_s3f8_t2   w_pe0/w_ce = 100   c/v_ti = 2261   grows the false shock

The middle and bottom runs share w_pe0/w_ce = 100 AND B0 = 7.03 T, so anything that
appears in the bottom panel and not the middle one cannot be a w_pe0/w_ce effect. What
differs across that pair is c/v_ti, by 10x.

WHAT TO WATCH, ON THE LOWER BRANCH. The C is ambient ions swept up by the piston's
ambipolar field (Sec. III C; forms with no B field at all). The NOTCH is the step at the
front of the C's LOWER branch, where that population rejoins the undisturbed upstream.
Healthy: the branch rolls over smoothly and a reflected population detaches above it (the
paper's feature II, whose appearance is t*_1). Failing: the branch steepens into a standing
vertical cliff a few tenths of a rho_i0 ahead of the piston, and never detaches.

THREE THINGS THAT MAKE THE PANELS COMPARABLE, none of them optional:

  * The window is piston-RELATIVE and measured in rho_i0. The three boxes are 5.29, 6.00
    and 3.29 rho_i0 and the pistons run at different absolute speeds, so a fixed z window
    would frame a different part of the structure in each row. rho_i0/d_i0 also runs
    10.4 -> 104 across these runs, so d_i0 would not do either.
  * Velocities are each run's OWN v_sh, positions each run's OWN rho_i0.
  * Frames are matched on t*w_ci0, not index: the runs write 52, 61 and 52 frames over the
    same t* range. A frame further than TOL from the requested time is left blank and
    reported rather than silently mislabelled.

The colour scale is LINEAR and fixed across all frames of a run (its own max over the whole
movie). Linear because the notch is a step in the DENSE lower branch and a log stretch
lifts the sparse tails until they compete with it; fixed because a per-frame scale makes a
structure that genuinely grows look constant.

    python scripts/movie_three_run_compare.py [--nframes 60] [--fps 10]
"""
import argparse, os, shutil, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.colors import Normalize, LogNorm      # noqa: E402
import kinshock                                       # noqa: E402
from kinshock import io as kio, plotting as P         # noqa: E402
from kinshock.units import C, ME                      # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
ACC, GOOD, BAD = "#ffd166", "#5ad2ff", "#ff6b57"
WIN_LO, WIN_HI = -0.22, 0.52      # rho_i0 about the piston front
VLO, VHI = -0.35, 1.85            # v_z / v_sh
NZ, NV = 300, 220
TOL = 0.05

RUNS = [
    ("ek_47keV_t2",      "runs/E_phase/ek_47keV_t2",        10, 226,  False),
    ("mb_r100_t2",       "runs/M_phase/mb_r100_t2",        100, 226,  False),
    ("R1_470eV_s3f8_t2", "runs/R1_phase/R1_470eV_s3f8_t2", 100, 2261, True),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nframes", type=int, default=60)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--norm", choices=("linear", "log"), default="linear")
    ap.add_argument("--out", default="media/comparison/three_run_phase.mp4")
    args = ap.parse_args()

    import yt; yt.set_log_level(50)
    grid = np.linspace(0.04, 2.00, args.nframes)

    # ---- one load pass: cache a histogram per (run, frame) ------------------
    cache, meta = {}, {}
    for lab, d, wpwc, cvti, bad in RUNS:
        cfg = kinshock.load(d)
        sc = kinshock.units.derive(cfg)
        mi = cfg["reference"]["mass_ratio"] * ME
        pfs = kio.plotfiles(d)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        meta[lab] = dict(sc=sc, wpwc=wpwc, cvti=cvti, bad=bad,
                         box=sc.domain_halfwidth / sc.rho_i0)
        for k, T in enumerate(grid):
            i = int(np.argmin(np.abs(times - T)))
            if abs(times[i] - T) > TOL:
                cache[(lab, k)] = None
                continue
            ad = kio.load_frame(pfs[i]).ds.all_data()
            z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.rho_i0
            p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
            g = np.sqrt(1 + (p / (mi * C)) ** 2)
            v = p / (mi * g) / sc.vsh_model
            w = np.asarray(ad[("amb_ions", "particle_weight")])
            zp = np.percentile(
                np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5) / sc.rho_i0
            zb = np.linspace(zp + WIN_LO, zp + WIN_HI, NZ + 1)
            vb = np.linspace(VLO, VHI, NV + 1)
            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            cache[(lab, k)] = dict(H=H, zp=zp, t=times[i])
        print(f"  loaded {lab}", flush=True)

    # per-run fixed scale over the WHOLE movie
    vmax = {lab: max((c["H"].max() for (l, _), c in cache.items()
                      if l == lab and c is not None), default=1.0)
            for lab, *_ in RUNS}

    # ---- render pass --------------------------------------------------------
    scratch = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))
    fdir = os.path.join(scratch, "three_run")
    shutil.rmtree(fdir, ignore_errors=True)
    os.makedirs(fdir, exist_ok=True)
    rel = np.linspace(WIN_LO, WIN_HI, NZ)
    vb = np.linspace(VLO, VHI, NV)
    for k in range(args.nframes):
        # libx264 + yuv420p REFUSES an odd pixel dimension, and it refuses it only at
        # encode time -- after the whole (expensive) load pass has already run. 9.2 x 9.4
        # at dpi 125 gave 1150 x 1175 and cost an encode. Keep both even.
        fig, ax = plt.subplots(3, 1, figsize=(9.2, 9.44), dpi=125,
                               facecolor=SURFACE, sharex=True)
        for r, (lab, d, wpwc, cvti, bad) in enumerate(RUNS):
            a = ax[r]
            a.set_facecolor("#111014")
            c = cache[(lab, k)]
            col = BAD if bad else INK2
            if c is None:
                a.annotate(f"no frame within {TOL} of t* = {grid[k]:.2f}",
                           xy=(.5, .5), xycoords="axes fraction", ha="center",
                           color=MUTED, fontsize=10)
            else:
                nrm = (Normalize(vmin=0, vmax=vmax[lab]) if args.norm == "linear"
                       else LogNorm(vmin=max(vmax[lab] * 1e-4, 1e-30), vmax=vmax[lab]))
                a.pcolormesh(rel, vb, np.ma.masked_where(c["H"] <= 0, c["H"]).T,
                             norm=nrm, cmap="magma", shading="auto", zorder=0)
                a.axvline(0, color=ACC, lw=1.4, zorder=3)
                a.axhline(0, color="#fff", lw=.5, alpha=.25, zorder=2)
                a.axhline(1, color=GOOD, lw=.9, ls=(0, (5, 4)), alpha=.85, zorder=2)
                a.annotate(rf"$t\,\omega_{{ci0}}$ = {c['t']:.2f}", xy=(.985, .90),
                           xycoords="axes fraction", ha="right", color=INK,
                           fontsize=10.5, family="monospace")
            a.set_xlim(WIN_LO, WIN_HI); a.set_ylim(VLO, VHI)
            a.tick_params(colors=MUTED, labelsize=8.5)
            for s in a.spines.values():
                s.set_color(BAD if bad else BASE)
            a.set_ylabel(f"{lab}\n" + rf"$\omega_{{pe0}}/\omega_{{ce}}$={wpwc}   "
                         rf"$c/v_{{ti}}$={cvti}" "\n" r"$v_z/v_{sh}$",
                         color=col, fontsize=9)
        ax[2].set_xlabel(r"$(z - z_{piston}) / \rho_{i0}$", color=INK2, fontsize=11)
        fig.suptitle("Ambient-ion phase space, piston-relative — the notch is the step on "
                     "the LOWER branch", color=INK, fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.965])
        # plotting.encode() globs frame_%03d.png -- match its contract exactly or
        # ffmpeg exits 254 after the whole (expensive) load pass has run.
        fig.savefig(os.path.join(fdir, f"frame_{k:03d}.png"), facecolor=SURFACE)
        plt.close(fig)
    print(f"  rendered {args.nframes} frames", flush=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    print(P.encode(fdir, args.out, fps=args.fps))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
