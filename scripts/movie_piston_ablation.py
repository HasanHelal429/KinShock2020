#!/usr/bin/env python3
"""Piston-plasma density and velocity, overplotted, for N runs at matched gyrotime.

    python scripts/movie_piston_ablation.py runs/E_phase/ek_47keV_t2 \
        runs/R1_phase/R1_470eV_s3f8_t2 runs/V_phase/vA_lowTab --tmax 2.0

One panel per run, stacked. Each panel carries TWO profiles of `piston_ions` on a shared
z axis: n/n_t on a LOG left axis, and the weighted-mean v_z/c on a linear right axis with a
+/-1 sigma band. The point is the ablation itself -- a rarefaction should show density
falling while velocity rises, and the band says whether that flow is cold and directed or
just hot and smeared.

n IS NORMALISED TO n_t, NOT n_amb. The piston starts at n_t and rarefies, so n/n_t runs
1 -> 0 and reads directly as "how much of the target is left here". n/n_amb would span
0-250 (n_t/n_amb = 250 exactly) and bury the structure. Log, because the rarefaction spans
orders of magnitude and a linear axis shows only the dense body.

v IS IN UNITS OF c, PER PANEL. These runs differ 43x in speed (peak piston v_z = 0.202 c,
0.0219 c, 0.0054 c) so one shared axis would flatten the two colder panels. Limits are fixed
across the WHOLE movie per run -- computed in a first pass -- so brightness and extent carry
information across frames rather than rescaling every step. In units of the ablation sound
speed all three would instead collapse onto 0-7.6 C_s,ab; that is the comparable view, and
this is the absolute one.

Runs march on a common t* grid (gyrotime). Nearest-frame lookup means a run whose data ends
early simply HOLDS its last frame rather than vanishing.
"""
import argparse, os, shutil, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import yaml                                         # noqa: E402
import kinshock                                     # noqa: E402
from kinshock import io as kio, plotting as P       # noqa: E402
from kinshock.units import C, ME                    # noqa: E402

SURFACE, PANEL, INK, INK2, MUTED, BASE = "#0d0c10", "#111014", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
C_N, C_V, C_BAND = "#ff8c42", "#5ad2ff", "#5ad2ff"
NB = 300


def _find(o, k):
    if isinstance(o, dict):
        if k in o:
            return o[k]
        for v in o.values():
            r = _find(v, k)
            if r is not None:
                return r
    return None


def profiles(fr, sc, mi, n_t, zb):
    ad = fr.ds.all_data()
    z = np.asarray(ad[("piston_ions", "particle_position_x")]) / sc.di0
    w = np.asarray(ad[("piston_ions", "particle_weight")])
    pm = np.asarray(ad[("piston_ions", "particle_momentum_z")])
    g = np.sqrt(1 + (pm / (mi * C)) ** 2)
    v = pm / (mi * g) / C
    dz_m = (zb[1] - zb[0]) * sc.di0
    sw, _ = np.histogram(z, bins=zb, weights=w)
    n = sw / dz_m / n_t                      # n/n_t
    swv, _ = np.histogram(z, bins=zb, weights=w * v)
    swv2, _ = np.histogram(z, bins=zb, weights=w * v * v)
    cnt, _ = np.histogram(z, bins=zb)
    ok = cnt >= 20
    mean = np.where(ok, np.divide(swv, sw, out=np.zeros_like(sw), where=sw > 0), np.nan)
    var = np.where(ok, np.divide(swv2, sw, out=np.zeros_like(sw), where=sw > 0) - mean**2,
                   np.nan)
    return np.where(n > 0, n, np.nan), mean, np.sqrt(np.clip(var, 0, None))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--win", type=float, nargs=2, default=[0.0, 24.0])
    ap.add_argument("--tmax", type=float, default=2.0)
    ap.add_argument("--nframes", type=int, default=60)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--nfloor", type=float, default=1e-4, help="log density floor, n/n_t")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import yt; yt.set_log_level(50)
    zb = np.linspace(args.win[0], args.win[1], NB + 1)
    zc = 0.5 * (zb[1:] + zb[:-1])

    R = []
    for run in args.runs:
        cfg = kinshock.load(run); sc = kinshock.units.derive(cfg)
        raw = yaml.safe_load(open(os.path.join(run, "config.yaml")))
        n_t = float(raw["plasma"]["piston"]["density_over_n0"]) * float(raw["reference"]["n0"])
        pfs = kio.plotfiles(run)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        R.append(dict(run=run, sc=sc, pfs=pfs, times=times, n_t=n_t,
                      mi=cfg["reference"]["mass_ratio"] * ME,
                      Tab=float(_find(raw, "theta_e_heat")) * 511000.0, B0=sc.B0))
        print(f"  {os.path.basename(run)}: {len(pfs)} frames to t*={times[-1]:.3f}",
              flush=True)
    grid = np.linspace(0.0, args.tmax, args.nframes)

    for r in R:
        seen, out = {}, []
        for T in grid:
            i = int(np.argmin(np.abs(r["times"] - T)))
            if i not in seen:
                seen[i] = profiles(kio.load_frame(r["pfs"][i]), r["sc"], r["mi"],
                                   r["n_t"], zb) + (r["times"][i],)
            out.append(seen[i])
        r["fr"] = out
        vs = np.concatenate([np.concatenate([f[1] - f[2], f[1] + f[2]]) for f in out])
        vs = vs[np.isfinite(vs)]
        lo, hi = np.nanpercentile(vs, 0.2), np.nanpercentile(vs, 99.8)
        pad = 0.08 * (hi - lo)
        r["vlim"] = (lo - pad, hi + pad)
        print(f"  {os.path.basename(r['run'])}: {len(seen)} distinct frames, "
              f"v_z/c fitted to [{r['vlim'][0]:+.4f}, {r['vlim'][1]:+.4f}]", flush=True)

    scratch = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))
    fdir = os.path.join(scratch, "piston_abl_" + "_".join(
        kinshock.load(x)["meta"]["run_id"] for x in args.runs))
    shutil.rmtree(fdir, ignore_errors=True); os.makedirs(fdir, exist_ok=True)

    for k in range(args.nframes):
        # even pixel dims: libx264+yuv420p rejects odd ones, and only at encode time.
        fig, ax = plt.subplots(len(R), 1, figsize=(9.0, 3.4 * len(R)), dpi=120,
                               facecolor=SURFACE, squeeze=False, sharex=True)
        for r, Rr in enumerate(R):
            a = ax[r][0]; a.set_facecolor(PANEL)
            n, mean, sd, t = Rr["fr"][k]
            a.semilogy(zc, n, color=C_N, lw=1.6, zorder=3)
            a.set_ylim(args.nfloor, 1.6)
            a.set_ylabel(r"$n/n_t$", color=C_N, fontsize=10)
            a.tick_params(axis="y", colors=C_N, labelsize=8.5)
            a.tick_params(axis="x", colors=MUTED, labelsize=8.5)
            b = a.twinx()
            b.fill_between(zc, mean - sd, mean + sd, color=C_BAND, alpha=0.22,
                           lw=0, zorder=1)
            b.plot(zc, mean, color=C_V, lw=1.5, zorder=2)
            b.axhline(0.0, color="#ffffff", lw=0.5, alpha=0.25, zorder=0)
            b.set_ylim(*Rr["vlim"])
            b.set_ylabel(r"$v_z/c$", color=C_V, fontsize=10)
            b.tick_params(axis="y", colors=C_V, labelsize=8.5)
            for s in list(a.spines.values()) + list(b.spines.values()):
                s.set_color(BASE)
            _T = f"{Rr['Tab']:,.1f}".removesuffix(".0")
            a.set_title(f"{_T} eV      {Rr['B0']:.2f} T"
                        + rf"          $t\omega_{{ci0}}$ = {t:.3f}",
                        color=INK, fontsize=11)
            a.set_xlim(*args.win)
        ax[-1][0].set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=11)
        fig.suptitle("piston ions  |  "
                     r"$n/n_t$ (log, left)  and  $v_z/c$ with $\pm1\sigma$ (right)  |  "
                     "velocity axis fitted per panel", color=INK, fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        fig.savefig(os.path.join(fdir, f"frame_{k:03d}.png"), facecolor=SURFACE)
        plt.close(fig)
    print(f"  rendered {args.nframes} frames", flush=True)
    tag = "_".join(kinshock.load(x)["meta"]["run_id"] for x in args.runs)
    out = args.out or f"media/comparison/piston_ablation_{tag}.mp4"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    print(P.encode(fdir, out, fps=args.fps))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
