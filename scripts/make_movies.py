#!/usr/bin/env python3
"""Movies for a KinShock2020 run -> media/<run_id>/ (config-driven, yt + ffmpeg).

  * shock_ni.mp4     : ion density line-out n_i/n_e0(z) vs time
  * shock_phase.mp4  : ambient+piston ion (z, v_z/v_sh) phase space vs time
                       (v_z > v_sh dotted line = reflected-ion threshold)

Frames are rendered to $MOVIE_SCRATCH (default the session scratchpad) then
encoded. Uses the same shock speed the figures use (model value, unless a
measured one is passed with --vsh-c).

Usage:  python scripts/make_movies.py [runs/R1] [--fps 8] [--vsh-c 0.14]
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io, plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

SCRATCH = os.environ.get(
    "MOVIE_SCRATCH",
    "/tmp/claude-1024/-home-hhelal/23de6041-80fd-48f5-8fec-f183097726ec/scratchpad/kinshock_frames",
)


def movie_ni(pfs, cfg, sc, fps):
    framedir = os.path.join(SCRATCH, cfg["meta"]["run_id"], "ni")
    os.makedirs(framedir, exist_ok=True)
    for i, pf in enumerate(pfs):
        fr = io.load_frame(pf)
        zc = np.asarray(fr.z_centers) / sc.di0
        n = io.species_density(fr, cfg["ion_species"])
        fig, ax = plt.subplots(figsize=(7.4, 4.0))
        ax.plot(zc, np.where(n > 0, n / sc.namb, np.nan), color=P.C_PISTON, lw=1.1)
        ax.axhline(1.0, color=P.C_REF, ls="--", lw=0.8)
        ax.set_yscale("log")
        ax.set_ylim(0.3, max(10.0, 5.0 * (np.nanmax(n) / sc.namb if n.max() else 1)))
        ax.set_xlabel(r"$z / d_{i0}$")
        ax.set_ylabel(r"$n_i / n_{e0}$")
        ax.set_title(rf"$t\,\omega_{{ci0}}={fr.time*sc.wci0:.2f}$")
        P.style_axes(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{i:03d}.png"), dpi=110)
        plt.close(fig)
    out = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "shock_ni.mp4")
    return P.encode(framedir, out, fps=fps)


def movie_phase(pfs, cfg, sc, vsh, fps):
    framedir = os.path.join(SCRATCH, cfg["meta"]["run_id"], "phase")
    os.makedirs(framedir, exist_ok=True)
    c = kinshock.units.C
    for i, pf in enumerate(pfs):
        fr = io.load_frame(pf)
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        for sp, col, lbl in (("piston_ions", P.C_PISTON, "piston ions"),
                             ("amb_ions", P.C_AMBIENT, "ambient ions")):
            z, uz = io.species_phase(fr, sp, sc, mass=sc.mi)
            if len(z):
                ax.plot(z / sc.di0, uz * c / vsh, ".", ms=0.5, alpha=0.25,
                        color=col, rasterized=True, label=lbl)
        ax.axhline(1.0, color=P.C_REF, ls=":", lw=0.8)
        ax.set_ylim(-2.0, 2.5)
        ax.set_xlabel(r"$z / d_{i0}$")
        ax.set_ylabel(r"$v_z / v_{sh}$")
        ax.set_title(rf"$t\,\omega_{{ci0}}={fr.time*sc.wci0:.2f}$")
        leg = ax.legend(frameon=False, fontsize=8, markerscale=15, loc="upper right")
        for lh in leg.legend_handles:
            lh.set_alpha(1.0)
        P.style_axes(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{i:03d}.png"), dpi=110)
        plt.close(fig)
    out = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "shock_phase.mp4")
    return P.encode(framedir, out, fps=fps)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", nargs="?", default=os.path.join(ROOT, "runs", "R1"))
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--vsh-c", type=float, default=None,
                    help="shock speed in units of c (default: model value from config)")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    vsh = args.vsh_c * kinshock.units.C if args.vsh_c else sc.vsh_model
    pfs = io.plotfiles(args.run_dir)
    print(f"{cfg['meta']['run_id']}: {len(pfs)} plotfiles -> movies (v_sh={vsh/kinshock.units.C:.4f} c)")
    movie_ni(pfs, cfg, sc, args.fps)
    movie_phase(pfs, cfg, sc, vsh, args.fps)


if __name__ == "__main__":
    main()
