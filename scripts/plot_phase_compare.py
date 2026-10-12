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
    """Stacked comparison movie in the R1-run style: `plotting.phase_distribution`.

    This is the SAME renderer the shock_phase.mp4 movies in the R1 runs use -- an additive
    two-colour image with the ambient in blue (C_AMBIENT) and the piston in orange
    (C_PISTON), overlapping toward white -- not a single-species density map. Both species
    are drawn, which is the point: the piston is what drives the C, so seeing where it is
    relative to the ambient is most of the reason to watch.

    STATIONARY WINDOW by default here (--frame box): one fixed absolute-z axis shared by
    both panels, so the distribution moves and the axes do not. Both runs are shown over the
    same 0..34 d_i0 whatever their own box (34.2 and 80.5), so features are at the same x.

    PEAKS ARE FIXED PER SPECIES PER RUN over the whole movie. phase_distribution's default
    is each species' own peak IN THAT FRAME, which its docstring flags as wrong for a movie:
    it rescales every frame independently so brightness carries no information across time.
    Per-species normalisation is kept -- the ambient is ~250x rarer than the piston and a
    single shared scale would erase it.
    """
    import yt; yt.set_log_level(50)
    from kinshock import io as _io, plotting as P
    from kinshock.units import C as _C

    z_edges = np.linspace(args.win[0], args.win[1], 421)
    v_edges = np.linspace(args.vlim[0], args.vlim[1], 241)

    runs = []
    for run in args.runs:
        cfg = kinshock.load(run); sc = kinshock.units.derive(cfg)
        pfs = kio.plotfiles(run)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        import yaml as _y
        _raw = _y.safe_load(open(os.path.join(run, "config.yaml")))
        _th = float(_find(_raw, "theta_e_heat"))
        runs.append(dict(run=run, cfg=cfg, sc=sc, pfs=pfs, times=times,
                         rid=cfg["meta"]["run_id"], meta=row_meta(run),
                         Tab=_th * 511000.0, B0=sc.B0,
                         vsh_c=sc.vsh_model / C,
                         cs_c=np.sqrt(_th / cfg["reference"]["mass_ratio"])))
        print(f"  {runs[-1]['rid']}: {len(pfs)} frames to t*={times[-1]:.3f}", flush=True)
    if args.match == "piston":
        # piston front vs frame, per run, then a shared travel grid
        for R in runs:
            tr = []
            for pf in R["pfs"]:
                ad = kio.load_frame(pf).ds.all_data()
                try:
                    tr.append(float(np.percentile(
                        np.asarray(ad[("piston_ions", "particle_position_x")]),
                        99.5)) / R["sc"].di0)
                except Exception:
                    tr.append(np.nan)
            R["travel"] = np.asarray(tr)
            print(f"  {R['rid']}: piston {np.nanmin(R['travel']):.2f} -> "
                  f"{np.nanmax(R['travel']):.2f} d_i0", flush=True)
        dmax = args.tmax if args.tmax else min(np.nanmax(R["travel"]) for R in runs)
        grid_key = np.linspace(float(min(np.nanmin(R["travel"]) for R in runs)),
                               dmax, args.nframes)
        for R in runs:
            R["key"] = R["travel"]
        print(f"  marching on PISTON TRAVEL 0 -> {dmax:.2f} d_i0 in {args.nframes} frames "
              f"({args.frame} frame, {args.norm})", flush=True)
    else:
        dmax = args.tmax if args.tmax else min(r["times"][-1] for r in runs)
        grid_key = np.linspace(0.0, dmax, args.nframes)
        for R in runs:
            R["key"] = R["times"]
        print(f"  sweeping t* 0 -> {dmax:.3f} in {args.nframes} frames "
              f"({args.frame} frame, {args.norm})", flush=True)
    grid_t = grid_key

    frames, peaks = [], []
    for R in runs:
        sc = R["sc"]; seen = {}; pk = {}
        for T in grid_t:
            i = int(np.nanargmin(np.abs(R["key"] - T)))
            if i in seen:
                continue
            fr = kio.load_frame(R["pfs"][i])
            sd = {}
            zp = np.nan
            for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
                z, uz, w = _io.species_phase_weighted(fr, sp, sc, mass=sc.mi)
                if not len(z):
                    continue
                if sp == "piston_ions":
                    zp = float(np.percentile(z, 99.5)) / sc.di0
                x = z / sc.di0
                if args.frame == "piston":
                    x = x - zp
                sd[key] = (x, uz * _C / sc.vsh_model, w)
                H, _, _ = np.histogram2d(sd[key][0], sd[key][1],
                                         bins=[z_edges, v_edges], weights=w)
                nzv = H[H > 0]
                if nzv.size:
                    pk.setdefault(key, []).append(nzv)
            seen[i] = (R["times"][i], sd, zp)
        frames.append([seen[int(np.nanargmin(np.abs(R["key"] - T)))] for T in grid_t])
        # EXPOSE ON A PERCENTILE OF POOLED CELLS, not the global max. The max is set by the
        # single brightest frame -- early, when the distribution is still concentrated -- so
        # normalising to it leaves every later frame nearly black once the population
        # spreads. A high percentile of all nonzero cells pooled over the movie is still
        # ONE fixed number per species (so brightness stays comparable across frames, which
        # is the whole reason not to use per-frame peaks) but exposes typical frames
        # properly and lets the few brightest cells clip.
        pk = {k: float(np.percentile(np.concatenate(v), args.expose))
              for k, v in pk.items()}
        peaks.append(pk)
        print(f"  {R['rid']}: {len(seen)} distinct frames, p{args.expose} peaks "
              + ", ".join(f"{k}={v:.3e}" for k, v in pk.items()), flush=True)

    scratch = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))
    fdir = os.path.join(scratch, "phase_cmp_" + "_".join(R["rid"] for R in runs))
    import shutil; shutil.rmtree(fdir, ignore_errors=True); os.makedirs(fdir, exist_ok=True)
    for k in range(args.nframes):
        # even pixel dims: libx264 + yuv420p rejects odd ones, and only at encode time,
        # after the whole load pass has been paid for.
        fig, ax = plt.subplots(len(runs), 1, figsize=(8.8, 4.4 * len(runs)), dpi=125,
                               squeeze=False, sharex=True)
        for r, R in enumerate(runs):
            a = ax[r][0]
            t, sd, zp = frames[r][k]
            # phase_distribution draws on the v_edges grid, which is in v_sh. Rescaling to
            # c is a pure relabelling of that axis by this run's own v_sh/c -- the data are
            # untouched, so the 43x span between these runs lands in the tick labels.
            sc_y = R["vsh_c"] if args.yunit == "c" else 1.0
            vline = {"vsh": 1.0, "cs_ab": 1.0 / 4.6, "none": None}[args.vline]
            vlab = {"vsh": r"$v_z=v_{sh}$", "cs_ab": r"$v_z=C_{s,ab}$",
                    "none": ""}[args.vline]
            P.phase_distribution(a, sd, z_edges, v_edges, vline=vline, vline_label=vlab,
                                 legend=(r == 0), norm=args.norm, peaks=peaks[r])
            if args.yunit == "c":
                import matplotlib.ticker as _mt
                a.yaxis.set_major_formatter(
                    _mt.FuncFormatter(lambda y, _p, s=sc_y: f"{y*s:.3f}"))
            if args.frame == "box" and np.isfinite(zp):
                a.axvline(zp, color="#ffd166", lw=1.0, alpha=.55, zorder=5)
            if args.label == "minimal":
                _T = f"{R['Tab']:,.1f}".removesuffix(".0")   # 47,012.1 / 470 / 25.6,
                a.set_title(f"{_T} eV      {R['B0']:.2f} T"     # never scientific notation
                            + rf"          $t\omega_{{ci0}}$ = {t:.3f}", fontsize=11)
            else:
                dzl, ppc, nd = R["meta"]
                a.set_title(f"{R['rid']}   " + rf"$dz/\lambda_D$={dzl:.2f}  ppc={ppc:.0f}  "
                            + rf"$N_D$={nd:.0f}" + rf"      $t\omega_{{ci0}}$ = {t:.3f}",
                            fontsize=10)
            a.set_ylabel(r"$v_z / c$" if args.yunit == "c" else r"$v_z / v_{sh}$",
                         fontsize=10)
        ax[-1][0].set_xlabel(r"$z / d_{i0}$" if args.frame == "box"
                             else r"$(z - z_{piston}) / d_{i0}$", fontsize=11)
        fig.suptitle(("stationary window" if args.frame == "box" else "piston-relative")
                     + "   |   " + ("dashed = $C_{s,ab}$" if args.vline == "cs_ab"
                                    else "dashed = $v_{sh}$")
                     + "   |   y axis per panel, "
                     + ("$v_z/c$" if args.yunit == "c" else "$v_z/v_{sh}$"), fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        fig.savefig(os.path.join(fdir, f"frame_{k:03d}.png"))
        plt.close(fig)
    print(f"  rendered {args.nframes} frames", flush=True)
    tag = "_".join(R["rid"] for R in runs)
    out = args.out or f"media/comparison/phase_{tag}_{args.frame}_{args.norm}.mp4"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(P.encode(fdir, out, fps=args.fps))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", type=float, nargs="*", default=[0.15, 0.30, 0.45, 0.60])
    ap.add_argument("--frame", choices=("piston", "box"), default="piston",
                    help="piston: x = z - z_piston, window follows the driver. "
                         "box: STATIONARY window in absolute z, the piston sweeps through "
                         "it -- use when the point is where things are, not where they are "
                         "relative to the driver.")
    ap.add_argument("--win", type=float, nargs=2, default=None,
                    help="window in d_i0. Default -2.5..5.5 about the piston (frame=piston) "
                         "or 0..34 absolute (frame=box).")
    ap.add_argument("--yunit", choices=("vsh", "c"), default="vsh",
                    help="v_z normalisation. 'c' rescales each panel by that run's own "
                         "v_sh/c, so per-panel limits differ by the same factor the runs "
                         "differ in v/c -- the disparity shows in the TICK LABELS, not the "
                         "geometry. --vlim is then read in units of v_sh and converted.")
    ap.add_argument("--vline", choices=("vsh", "cs_ab", "none"), default="vsh",
                    help="what the dashed line marks. cs_ab = the ABLATION sound speed "
                         "sqrt(Z kT_e,ab/m_i) = v_sh/4.6, the partner to a T_ab label; "
                         "the upstream C_s0 is ~10x smaller again and would sit on zero.")
    ap.add_argument("--label", choices=("full", "minimal"), default="full",
                    help="minimal = ablation temperature and B0 only, no run id or grid "
                         "metadata. Run ids then live only in the filename/caption.")
    ap.add_argument("--cmap", default="magma",
                    help="matplotlib colormap; cividis is the repo's CVD-safe "
                         "blue->orange phase-space map (make_figures FIG7_CMAP)")
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
    ap.add_argument("--expose", type=float, default=99.5,
                    help="percentile of pooled nonzero cells used as each species' fixed "
                         "normalisation peak; lower = brighter. 100 = the strict max.")
    ap.add_argument("--match", choices=("tstar", "piston"), default="tstar",
                    help="what to march the runs on. tstar: common t* grid. piston: common "
                         "PISTON-TRAVEL grid -- required whenever v_p differs between runs, "
                         "e.g. V_phase, where v_p ~ sqrt(T_ab) is 4.28x slower than the "
                         "reference while w_ci0 is unchanged, so equal t* is NOT equal stage.")
    ap.add_argument("--tmax", type=float, default=None,
                    help="movie: last t* to sweep to (default: min over runs)")
    args = ap.parse_args()
    if args.win is None:
        args.win = [-2.5, 5.5] if args.frame == "piston" else [0.0, 34.0]
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
            x = (z - zp) if args.frame == "piston" else z
            H, _, _ = np.histogram2d(x, v, bins=[zb, vb], weights=w)
            grid[(r, c)] = (times[i], H, zp)
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
            t, H, zp = g
            vmax = {"panel": H.max(), "row": rowmax[r], "shared": shared}[args.scale]
            nrm = (Normalize(0, vmax) if args.norm == "linear"
                   else LogNorm(vmin=max(vmax * 1e-4, 1e-30), vmax=vmax))
            a.pcolormesh(0.5 * (zb[1:] + zb[:-1]), 0.5 * (vb[1:] + vb[:-1]),
                         np.ma.masked_where(H <= 0, H).T, norm=nrm,
                         cmap=args.cmap, shading="auto", zorder=0)
            a.axhline(1.0, color=VSH, lw=1.0, ls=(0, (6, 4)), alpha=.9, zorder=3)
            a.axhline(0.0, color="#ffffff", lw=.5, alpha=.22, zorder=2)
            a.axvline(0.0 if args.frame == "piston" else zp,
                      color=PISTON, lw=1.6, alpha=.95, zorder=4)
            if r == 0:
                a.set_title(rf"$t\,\omega_{{ci0}}$ = {t:.2f}", color=INK, fontsize=11.5)
            if c == 0:
                a.set_ylabel(f"{rid}\n" + rf"$dz/\lambda_D$={dzl:.2f}  ppc={ppc:.0f}"
                             + f"\n$N_D$={nd:.0f}\n\n" + r"$v_z/v_{sh}$",
                             color=INK2, fontsize=8.5)
            if r == len(args.runs) - 1:
                a.set_xlabel(r"$(z - z_{piston}) / d_{i0}$" if args.frame == "piston"
                             else r"$z / d_{i0}$", color=INK2, fontsize=10)

    fig.suptitle(f"Ambient-ion phase space, "
                 f"{'piston-relative' if args.frame == 'piston' else 'stationary window'}  |  "
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
