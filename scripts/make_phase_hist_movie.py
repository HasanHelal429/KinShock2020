#!/usr/bin/env python3
"""Phase-space movie from the BINNED diagnostic (WarpX ParticleHistogram2D, PS_*_zuz).

    python scripts/make_phase_hist_movie.py <run_dir> [--fps 15] [--norm asinh|linear]
                                            [--fixed-scale] [--vsh-c 0.0326] [--out DIR]

Two panels on a shared z axis, one frame per histogram dump (deck.render writes them at
~10x the raw-particle cadence, so this is the smooth-in-time movie; make_movies.py's
shock_phase.mp4 is the particle-dump one):

  * ions      piston + ambient, v_z / v_sh       (dotted: v_z = v_sh, reflected-ion line)
  * electrons piston + ambient, v_z / v_te,ab

The rendering is plotting.phase_image -- the same stretch, per-species normalisation and
two-colour tint as shock_phase.mp4 -- so the two movies read the same way. Default is
asinh normalised per species per frame (make_movies' convention); --fixed-scale uses each
species' peak over ALL frames instead, so brightness is comparable across time.

The histogram ordinate is u = gamma*beta. Every species here is non-relativistic
(gamma - 1 < 1e-2), so u*c is v to that accuracy.

Reads openPMD (ADIOS2 .bp5) with openpmd_api -- use the WarpX venv. The shock speed
follows make_movies: --vsh-c > shock_fit.yaml > the config's model value.
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import metrics  # noqa: E402
from kinshock import plotting as P  # noqa: E402

PANELS = (  # (title, axis label, {role: histogram name}, velocity unit key)
    ("ions", r"$v_z / v_{sh}$", {"piston": "PS_pi_zuz", "ambient": "PS_ai_zuz"}, "vsh"),
    ("electrons", r"$v_z / v_{te,ab}$", {"piston": "PS_pe_zuz", "ambient": "PS_ae_zuz"}, "vte"),
)


def open_series(run_dir, name):
    import openpmd_api as opmd
    path = os.path.join(run_dir, "diags", "reducedfiles", name, "openpmd_%06T.bp5")
    return opmd.Series(path, opmd.Access.read_only)


def load_hist(series, it_num):
    """(H[nv, nz], extent [z0, z1, u0, u1] in m and gamma*beta, time [s]) of one dump."""
    import openpmd_api as opmd
    it = series.iterations[it_num]
    m = it.meshes["data"]
    rc = m[opmd.Mesh_Record_Component.SCALAR]
    H = rc.load_chunk()
    series.flush()
    du, dz = m.grid_spacing
    u0, z0 = m.grid_global_offset
    nv, nz = H.shape
    return np.array(H), [z0, z0 + nz * dz, u0, u0 + nv * du], it.time * it.time_unit_SI


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--norm", choices=("asinh", "linear"), default="asinh")
    ap.add_argument("--fixed-scale", action="store_true",
                    help="normalise each species to its peak over all frames")
    ap.add_argument("--vsh-c", type=float, default=None)
    ap.add_argument("--out", default=None, help="output dir (default media/<phase>/<id>)")
    ap.add_argument("--every", type=int, default=1, help="use every Nth dump")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    C = kinshock.units.C
    run_id = cfg["meta"]["run_id"]
    if args.vsh_c:
        vsh, src = args.vsh_c * C, "--vsh-c"
    else:
        fit = metrics.load_shock_fit(args.run_dir, sc)
        vsh, src = ((fit.v_sh, "by-eye fit (shock_fit.yaml)") if fit is not None
                    else (sc.vsh_model, "config MODEL value -- no shock_fit.yaml"))
    vunit = {"vsh": vsh, "vte": sc.vte_ab}

    series = {h: open_series(args.run_dir, h)
              for _, _, roles, _ in PANELS for h in roles.values()}
    its = sorted(set.intersection(*(set(s.iterations) for s in series.values())))
    its = its[::args.every]
    print(f"{run_id}: {len(its)} histogram dumps -> movie "
          f"(v_sh = {vsh / C:.4f} c from {src}; norm {args.norm}"
          f"{', fixed scale' if args.fixed_scale else ', per-frame scale'})")

    peaks = {}
    if args.fixed_scale:
        for h, s in series.items():
            peaks[h] = max(float(load_hist(s, i)[0].max()) for i in its)

    framedir = os.path.join(os.environ.get("MOVIE_SCRATCH", os.path.join(args.run_dir, ".frames")),
                            run_id, "phase_hist")
    os.makedirs(framedir, exist_ok=True)
    for f in os.listdir(framedir):
        if f.startswith("frame_"):
            os.remove(os.path.join(framedir, f))

    for k, it in enumerate(its):
        fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.6), sharex=True)
        t = None
        for ax, (title, ylab, roles, vk) in zip(axes, PANELS):
            hists, ext = {}, None
            for role, h in roles.items():
                H, ext, t = load_hist(series[h], it)
                hists[role] = H
            extent = [ext[0] / sc.di0, ext[1] / sc.di0,
                      ext[2] * C / vunit[vk], ext[3] * C / vunit[vk]]
            pk = ({role: peaks[h] for role, h in roles.items()} if peaks else None)
            P.phase_image(ax, hists, extent, vline=1.0 if vk == "vsh" else None,
                          legend=True, norm=args.norm, peaks=pk)
            ax.set_ylabel(ylab)
            ax.text(0.985, 0.965, title, transform=ax.transAxes, va="top", ha="right",
                    fontsize=9, color="w")
        axes[-1].set_xlabel(r"$z / d_{i0}$")
        axes[0].set_title(rf"{run_id}   $t\,\omega_{{ci0}} = {t * sc.wci0:.3f}$", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{k:03d}.png"), dpi=110)
        plt.close(fig)

    out_dir = args.out or P.media_dir(run_id=run_id)
    os.makedirs(out_dir, exist_ok=True)
    name = f"phase_hist_{args.norm}{'_fixed' if args.fixed_scale else ''}.mp4"
    return P.encode(framedir, os.path.join(out_dir, name), fps=args.fps)


if __name__ == "__main__":
    main()
