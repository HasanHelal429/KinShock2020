#!/usr/bin/env python3
"""Phase-space movie with a selectable colour stretch -- LINEAR by default.

WHY THIS EXISTS ALONGSIDE make_movies.py. `make_movies.py` renders the phase space with an
asinh stretch normalised to each species' peak IN THE CURRENT FRAME. That is the right
default for a single figure -- it makes the ~250x-rarer ambient population and its faint
reflected-ion beam visible at once -- but it does two things that get in the way when the
question is "how much is actually there":

  1. asinh flatters faint tails by construction, so a sparse beam can read as bright as a
     dense core.
  2. per-FRAME normalisation rescales every frame independently, so brightness carries no
     information across time: a structure whose density genuinely doubles looks unchanged,
     and one that fades looks constant.

This script fixes both. One pass caches the weighted histograms and takes each species'
peak over ALL frames; the second renders against that fixed scale, so brightness is
comparable frame to frame. Per-species normalisation is KEPT: the ambient is ~250x rarer
than the piston and a single shared scale would erase it entirely.

Linear is the honest view of where the bulk sits. Pass --norm asinh to get the tails back,
still on the fixed cross-frame scale.

    python scripts/make_phase_movie_norm.py runs/R1_phase/<ID> [--norm linear] [--fps 8]
"""
import argparse
import os
import sys

import numpy as np

ROOT = "/pscratch/sd/h/hhelal/KinShock2020"
sys.path.insert(0, os.path.join(ROOT, "src"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io, plotting as P  # noqa: E402

SCRATCH = os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, ".movie_scratch"))


def resolve_vsh(run_dir, sc, override):
    """Same precedence as make_movies / make_figures: --vsh-c > shock_fit.yaml > model."""
    c = kinshock.units.C
    if override:
        return override * c, "--vsh-c"
    try:
        fit = kinshock.metrics.load_shock_fit(run_dir)
    except Exception:
        fit = None
    if fit and fit.get("vsh_over_c"):
        return float(fit["vsh_over_c"]) * c, "shock_fit.yaml"
    return sc.vsh_model, "config MODEL value (no shock_fit.yaml)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--norm", choices=("linear", "asinh"), default="linear")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--vsh-c", type=float, default=None,
                    help="shock speed in units of c (default: shock_fit.yaml, else model)")
    ap.add_argument("--out", default=None,
                    help="output filename (default shock_phase_<norm>.mp4)")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    c = kinshock.units.C
    run_id = cfg["meta"]["run_id"]

    vsh, src = resolve_vsh(args.run_dir, sc, args.vsh_c)
    print(f"{run_id}: v_sh = {vsh / c:.4f} c from {src}")

    pfs = io.plotfiles(args.run_dir)
    Hd = sc.domain_halfwidth / sc.di0
    z_edges = np.linspace(0, Hd, 421)
    v_edges = np.linspace(-2.5, 3.0, 241)

    # ---- pass 1: cache histograms, find each species' peak over ALL frames ----------
    print(f"  pass 1/2: histogramming {len(pfs)} frames", flush=True)
    frames, peaks = [], {}
    for i, pf in enumerate(pfs):
        fr = io.load_frame(pf)
        sd = {}
        for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
            z, uz, w = io.species_phase_weighted(fr, sp, sc, mass=sc.mi)
            if len(z):
                H, _, _ = np.histogram2d(z / sc.di0, uz * c / vsh,
                                         bins=[z_edges, v_edges], weights=w)
                sd[key] = H
                peaks[key] = max(peaks.get(key, 0.0), float(H.max()))
        frames.append((fr.time * sc.wci0, sd))
        if (i + 1) % 10 == 0:
            print(f"    {i + 1}/{len(pfs)}", flush=True)
    print("  peaks over all frames: " +
          "  ".join(f"{k}={v:.4e}" for k, v in peaks.items()))

    # ---- pass 2: render against the fixed scale -------------------------------------
    framedir = os.path.join(SCRATCH, run_id, f"phase_{args.norm}")
    os.makedirs(framedir, exist_ok=True)
    print(f"  pass 2/2: rendering ({args.norm}, fixed cross-frame scale)", flush=True)
    colors = P.PHASE_COLORS
    for i, (t, sd) in enumerate(frames):
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        rgb = np.zeros((len(v_edges) - 1, len(z_edges) - 1, 3))
        for key, H in sd.items():
            pk = peaks[key]
            if args.norm == "linear":
                inten = (np.clip(H / pk, 0.0, 1.0) if pk > 0 else H).T
            else:
                inten = P._asinh_norm_to(H, pk).T
            rgb += inten[..., None] * np.asarray(colors.get(key, (1, 1, 1)))
        ax.imshow(np.clip(rgb, 0, 1), origin="lower", aspect="auto",
                  interpolation="nearest",
                  extent=[z_edges[0], z_edges[-1], v_edges[0], v_edges[-1]])
        ax.set_facecolor(P.PHASE_BG)
        ax.axhline(1.0, color="w", ls=":", lw=0.9, alpha=0.7)
        ax.text(0.985, 1.0, r"$v_z=v_{sh}$", transform=ax.get_yaxis_transform(),
                va="bottom", ha="right", fontsize=7, color="w", alpha=0.7)
        y = 0.965
        for key in sd:
            ax.text(0.025, y, key, transform=ax.transAxes, va="top", ha="left",
                    fontsize=8, fontweight="bold", color=colors.get(key, (1, 1, 1)))
            y -= 0.085
        ax.set_xlim(z_edges[0], z_edges[-1])
        ax.set_ylim(v_edges[0], v_edges[-1])
        ax.set_xlabel(r"$z / d_{i0}$")
        ax.set_ylabel(r"$v_z / v_{sh}$")
        ax.set_title(rf"$t\,\omega_{{ci0}}={t:.2f}$   ({args.norm}, fixed scale)")
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{i:03d}.png"), dpi=110)
        plt.close(fig)

    out = os.path.join(P.media_dir(run_id=run_id),
                       args.out or f"shock_phase_{args.norm}.mp4")
    print(P.encode(framedir, out, fps=args.fps))


if __name__ == "__main__":
    main()
