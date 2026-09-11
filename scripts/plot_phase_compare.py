#!/usr/bin/env python3
"""Ambient-ion phase space for two or more runs at matched times, PISTON-RELATIVE.

    python scripts/plot_phase_compare.py runs/R1_phase/R1_470eV_s3f8_t2 \
                                         runs/S_phase/ss_dz4_ppc400

WHY PISTON-RELATIVE. The runs being compared have different boxes (34.2 vs 12.0 d_i0) and
their pistons are at different absolute z, so an absolute axis both wastes most of the panel
on undisturbed upstream and puts the two runs' features at different x. The interesting
structure -- the C, and the notch on its lower branch -- lives within a couple of d_i0 of
the piston front, which at full-box scale is close to sub-pixel. x = 0 is the piston front
(99.5th percentile of piston_ions) in every panel, so features are directly comparable
across rows.

LINEAR COLOUR BY DEFAULT. A log stretch lifts the sparse tails until they compete with the
dense lower branch, which is exactly where the notch lives -- so log actively works against
reading this feature. Per-panel scaling by default (each panel to its own max) because the
runs differ in particle count; --scale row or shared when relative density matters more than
shape.

The row labels carry dz/lambda_D and N_D = ppc*lambda_D/dz, because notch ONSET TIME is
ordered by N_D (RESULTS 2026-09-05) and any two runs being compared should be read with that
in view: these two differ 2x in dz/lambda_D and 2x in N_D in opposite directions.
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.colors import Normalize, LogNorm    # noqa: E402
import yaml                                         # noqa: E402
import kinshock                                     # noqa: E402
from kinshock import io as kio                      # noqa: E402
from kinshock.units import C, ME                    # noqa: E402

SURFACE, INK, INK2, MUTED = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798"
BASE, PISTON, VSH = "#3a3745", "#ffd166", "#5ad2ff"


def _find(o, k):
    if isinstance(o, dict):
        if k in o:
            return o[k]
        for v in o.values():
            r = _find(v, k)
            if r is not None:
                return r
    return None


def row_meta(run):
    raw = yaml.safe_load(open(os.path.join(run, "config.yaml")))
    dz = float(raw["geometry"]["dz_over_de"])
    th0 = float(_find(raw, "theta_0"))
    amb = float(raw["plasma"]["ambient"]["density_over_n0"])
    lD = np.sqrt(th0) / np.sqrt(amb)          # lambda_D,amb in d_e,ab
    ppcd = _find(raw, "ppc")
    ppc = float(ppcd["ambient"] if isinstance(ppcd, dict) else ppcd)
    return dz / lD, ppc, ppc * lD / dz


def movie(args):
    """One panel per run, stacked, swept over a COMMON time grid.

    The colour scale is fixed per run across the whole movie, from a first pass. A
    per-frame scale makes a feature that grows look constant, which is precisely the
    thing these two runs are being compared for.

    Runs have different output cadence (0.039 vs 0.112 t* here), so a common grid plus
    nearest-frame lookup is the only way to march them together; each panel is stamped with
    the frame's TRUE t*, so the two never silently drift apart.
    """
    import yt; yt.set_log_level(50)
    from kinshock import plotting as P
    nz, nv = 380, 220
    zb = np.linspace(args.win[0], args.win[1], nz + 1)
    vb = np.linspace(args.vlim[0], args.vlim[1], nv + 1)

    runs = []
    for run in args.runs:
        cfg = kinshock.load(run); sc = kinshock.units.derive(cfg)
        pfs = kio.plotfiles(run)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        runs.append(dict(run=run, cfg=cfg, sc=sc, pfs=pfs, times=times,
                         mi=cfg["reference"]["mass_ratio"] * ME,
                         rid=cfg["meta"]["run_id"], meta=row_meta(run)))
        print(f"  {runs[-1]['rid']}: {len(pfs)} frames to t*={times[-1]:.3f}", flush=True)
    tmax = args.tmax if args.tmax else min(r["times"][-1] for r in runs)
    grid_t = np.linspace(0.0, tmax, args.nframes)
    print(f"  sweeping t* 0 -> {tmax:.3f} in {args.nframes} frames", flush=True)

    frames = [[] for _ in runs]
    for r, R in enumerate(runs):
        seen = {}
        for T in grid_t:
            i = int(np.argmin(np.abs(R["times"] - T)))
            if i not in seen:
                ad = kio.load_frame(R["pfs"][i]).ds.all_data()
                z = np.asarray(ad[("amb_ions", "particle_position_x")]) / R["sc"].di0
                pm = np.asarray(ad[("amb_ions", "particle_momentum_z")])
                g = np.sqrt(1 + (pm / (R["mi"] * C)) ** 2)
                v = pm / (R["mi"] * g) / R["sc"].vsh_model
                w = np.asarray(ad[("amb_ions", "particle_weight")])
                zp = float(np.percentile(
                    np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5)) / R["sc"].di0
                H, _, _ = np.histogram2d(z - zp, v, bins=[zb, vb], weights=w)
                seen[i] = (R["times"][i], H)
            frames[r].append(seen[i])
        print(f"  {R['rid']}: loaded {len(seen)} distinct frames", flush=True)
    vmax = [max(h.max() for _, h in fr) for fr in frames]

    scratch = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))
    fdir = os.path.join(scratch, "phase_cmp_" + "_".join(R["rid"] for R in runs))
    import shutil; shutil.rmtree(fdir, ignore_errors=True); os.makedirs(fdir, exist_ok=True)
    zc = 0.5 * (zb[1:] + zb[:-1]); vc = 0.5 * (vb[1:] + vb[:-1])
    for k in range(args.nframes):
        # even pixel dims: libx264 + yuv420p rejects odd ones, and only at encode time,
        # after the whole load pass has been paid for.
        fig, ax = plt.subplots(len(runs), 1, figsize=(8.8, 4.4 * len(runs)), dpi=125,
                               facecolor=SURFACE, squeeze=False, sharex=True)
        for r, R in enumerate(runs):
            a = ax[r][0]; a.set_facecolor("#111014")
            t, H = frames[r][k]
            nrm = (Normalize(0, vmax[r]) if args.norm == "linear"
                   else LogNorm(vmin=max(vmax[r] * 1e-4, 1e-30), vmax=vmax[r]))
            a.pcolormesh(zc, vc, np.ma.masked_where(H <= 0, H).T, norm=nrm,
                         cmap="magma", shading="auto", zorder=0)
            a.axhline(1.0, color=VSH, lw=1.0, ls=(0, (6, 4)), alpha=.9, zorder=3)
            a.axhline(0.0, color="#ffffff", lw=.5, alpha=.22, zorder=2)
            a.axvline(0.0, color=PISTON, lw=1.6, zorder=4)
            dzl, ppc, nd = R["meta"]
            a.annotate(f"{R['rid']}   " + rf"$dz/\lambda_D$={dzl:.2f}  ppc={ppc:.0f}  "
                       + rf"$N_D$={nd:.0f}", xy=(0.015, 0.93), xycoords="axes fraction",
                       color=INK, fontsize=10)
            a.annotate(rf"$t\omega_{{ci0}}$ = {t:.3f}", xy=(0.985, 0.93),
                       xycoords="axes fraction", ha="right", color=VSH, fontsize=10,
                       family="monospace")
            a.set_ylabel(r"$v_z/v_{sh}$", color=INK2, fontsize=10)
            a.set_xlim(args.win); a.set_ylim(args.vlim)
            a.tick_params(colors=MUTED, labelsize=8.5)
            for sp in a.spines.values():
                sp.set_color(BASE)
        ax[-1][0].set_xlabel(r"$(z - z_{piston}) / d_{i0}$", color=INK2, fontsize=11)
        fig.suptitle(f"piston-relative ambient-ion phase space  |  {args.norm} colour, "
                     f"fixed per run", color=INK, fontsize=11.5)
        fig.tight_layout(rect=[0, 0, 1, 0.965])
        fig.savefig(os.path.join(fdir, f"frame_{k:03d}.png"), facecolor=SURFACE)
        plt.close(fig)
    print(f"  rendered {args.nframes} frames", flush=True)
    tag = "_".join(R["rid"] for R in runs)
    out = args.out or f"media/comparison/phase_{tag}_{args.norm}.mp4"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(P.encode(fdir, out, fps=args.fps))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", type=float, nargs="*", default=[0.15, 0.30, 0.45, 0.60])
    ap.add_argument("--win", type=float, nargs=2, default=[-2.5, 5.5],
                    help="z - z_piston window in d_i0")
    ap.add_argument("--vlim", type=float, nargs=2, default=[-0.8, 2.2])
    ap.add_argument("--norm", choices=("linear", "log"), default="linear")
    ap.add_argument("--scale", choices=("panel", "row", "shared"), default="panel")
    ap.add_argument("--tol", type=float, default=0.06,
                    help="max |t* - requested| to accept a frame. Runs differ in output "
                         "cadence (0.039 vs 0.112 here), so exact matching is impossible; "
                         "panel titles carry the TRUE t* so any mismatch stays visible.")
    ap.add_argument("--out", default=None)
    ap.add_argument("--movie", action="store_true",
                    help="sweep time and encode an mp4 instead of a static grid")
    ap.add_argument("--nframes", type=int, default=60)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--tmax", type=float, default=None,
                    help="movie: last t* to sweep to (default: min over runs)")
    args = ap.parse_args()
    if args.movie:
        return movie(args)

    import yt; yt.set_log_level(50)
    nz, nv = 380, 220
    zb = np.linspace(args.win[0], args.win[1], nz + 1)
    vb = np.linspace(args.vlim[0], args.vlim[1], nv + 1)

    grid = {}
    print(f"{'run':24}{'t*':>7}{'z_piston':>10}{'Hmax':>12}")
    for r, run in enumerate(args.runs):
        cfg = kinshock.load(run); sc = kinshock.units.derive(cfg)
        mi = cfg["reference"]["mass_ratio"] * ME
        pfs = kio.plotfiles(run)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        for c, T in enumerate(args.times):
            i = int(np.argmin(np.abs(times - T)))
            if abs(times[i] - T) > args.tol:
                grid[(r, c)] = None
                print(f"{os.path.basename(run):24}{T:7.2f}   no frame within {args.tol}")
                continue
            ad = kio.load_frame(pfs[i]).ds.all_data()
            z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
            pm = np.asarray(ad[("amb_ions", "particle_momentum_z")])
            g = np.sqrt(1 + (pm / (mi * C)) ** 2)
            v = pm / (mi * g) / sc.vsh_model
            w = np.asarray(ad[("amb_ions", "particle_weight")])
            zp = float(np.percentile(
                np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5)) / sc.di0
            H, _, _ = np.histogram2d(z - zp, v, bins=[zb, vb], weights=w)
            grid[(r, c)] = (times[i], H)
            print(f"{os.path.basename(run):24}{times[i]:7.3f}{zp:10.2f}{H.max():12.3e}",
                  flush=True)

    vals = [g[1].max() for g in grid.values() if g is not None]
    shared = max(vals) if vals else 1.0
    rowmax = {r: max([grid[(r, c)][1].max() for c in range(len(args.times))
                      if grid.get((r, c)) is not None] or [1.0])
              for r in range(len(args.runs))}

    fig, ax = plt.subplots(len(args.runs), len(args.times),
                           figsize=(3.35 * len(args.times), 3.5 * len(args.runs)),
                           dpi=125, facecolor=SURFACE, squeeze=False,
                           sharex=True, sharey=True)
    for r, run in enumerate(args.runs):
        dzl, ppc, nd = row_meta(run)
        rid = kinshock.load(run)["meta"]["run_id"]
        for c in range(len(args.times)):
            a = ax[r][c]; a.set_facecolor("#111014")
            for s in a.spines.values():
                s.set_color(BASE)
            a.tick_params(colors=MUTED, labelsize=8.5)
            g = grid.get((r, c))
            if g is None:
                a.annotate("no frame", xy=(.5, .5), xycoords="axes fraction",
                           ha="center", color=MUTED, fontsize=9)
                continue
            t, H = g
            vmax = {"panel": H.max(), "row": rowmax[r], "shared": shared}[args.scale]
            nrm = (Normalize(0, vmax) if args.norm == "linear"
                   else LogNorm(vmin=max(vmax * 1e-4, 1e-30), vmax=vmax))
            a.pcolormesh(0.5 * (zb[1:] + zb[:-1]), 0.5 * (vb[1:] + vb[:-1]),
                         np.ma.masked_where(H <= 0, H).T, norm=nrm,
                         cmap="magma", shading="auto", zorder=0)
            a.axhline(1.0, color=VSH, lw=1.0, ls=(0, (6, 4)), alpha=.9, zorder=3)
            a.axhline(0.0, color="#ffffff", lw=.5, alpha=.22, zorder=2)
            a.axvline(0.0, color=PISTON, lw=1.6, alpha=.95, zorder=4)
            if r == 0:
                a.set_title(rf"$t\,\omega_{{ci0}}$ = {t:.2f}", color=INK, fontsize=11.5)
            if c == 0:
                a.set_ylabel(f"{rid}\n" + rf"$dz/\lambda_D$={dzl:.2f}  ppc={ppc:.0f}"
                             + f"\n$N_D$={nd:.0f}\n\n" + r"$v_z/v_{sh}$",
                             color=INK2, fontsize=8.5)
            if r == len(args.runs) - 1:
                a.set_xlabel(r"$(z - z_{piston}) / d_{i0}$", color=INK2, fontsize=10)

    fig.suptitle("Ambient-ion phase space, piston-relative  |  "
                 f"{args.norm} colour ({args.scale})  |  "
                 "yellow = piston front, dashed = $v_{sh}$",
                 color=INK, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    tag = "_".join(kinshock.load(r)["meta"]["run_id"] for r in args.runs)
    out = args.out or f"media/comparison/phase_{tag}_{args.norm}.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
