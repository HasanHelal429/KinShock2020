#!/usr/bin/env python3
"""Phase-space movie for a HYBRID run (kinetic ions, fluid electrons) made by H-PICShock.

    python scripts/hybrid_phase_movie.py runs/hybrid_phase/H3_470eV_dense [--fps 12]
            [--vsh-kms 4182] [--field-smooth-de 2] [--max-frames N] [--out DIR]

H-PICShock decks use their own config schema, which kinshock.load cannot read, so the scales
come from the run's warpx_used_inputs (my_constants: di, vA, B0, Mi, namb), resolved with
kinshock.deck.parse_inputs/resolve_constants. Three panels on a shared z axis, one frame per
`diag_phase` particle dump (matched to the `diag_fields` dump at the same or latest earlier
step):

  * ions      piston + ambient (z, v_z / v_sh), the same renderer as the kinetic movies
  * fields    B_x/B0 (left) and E_z/(v_A B0) (right), running mean --field-smooth-de
  * electrons the fluid has no phase space: T_e [eV] (left) and n_e/n_e0 per species (right)

v_sh defaults to the kinetic 470 eV family's model value (4182 km/s = 13.95 v_A), so the
v_z = v_sh line sits where it does in the kinetic movies; pass --vsh-kms to change it.
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

from kinshock import deck  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock import plotting as P  # noqa: E402

QE = 1.602176634e-19


def scales_from_used_inputs(run_dir):
    d = deck.parse_inputs(os.path.join(run_dir, "warpx_used_inputs"))
    ns = deck.resolve_constants(d, allow_unresolved=True)
    need = ("di", "vA", "B0", "Mi", "namb")
    missing = [k for k in need if k not in ns]
    if missing:
        raise SystemExit(f"warpx_used_inputs lacks my_constants {missing}")
    s = {k: float(ns[k]) for k in need}
    s["wci0"] = QE * s["B0"] / s["Mi"]
    return s


def dumps(run_dir, prefix):
    dd = os.path.join(run_dir, "diags")
    paths = [os.path.join(dd, p) for p in os.listdir(dd)
             if p.startswith(prefix) and p[len(prefix):].isdigit()]
    return {kio._step_key(p)[0]: p for p in paths}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--vsh-kms", type=float, default=4182.0)
    ap.add_argument("--field-smooth-de", type=float, default=2.0)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--ion-vlim", type=float, nargs=2, default=(-0.5, 2.5),
                    metavar=("LO", "HI"), help="ion panel v_z/v_sh axis (default -0.5 2.5)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    run_id = os.path.basename(run_dir)
    s = scales_from_used_inputs(run_dir)
    vsh = args.vsh_kms * 1e3
    de = s["di"] / np.sqrt(s["Mi"] / 9.1093837015e-31)          # d_e at n_amb
    EB = s["vA"] * s["B0"]

    phase = dumps(run_dir, "diag_phase")
    fields = dumps(run_dir, "diag_fields")
    steps = sorted(phase)
    if args.max_frames:
        steps = steps[:args.max_frames]
    fsteps = np.array(sorted(fields))
    print(f"{run_id}: {len(steps)} particle dumps, d_i0 = {s['di'] * 1e6:.2f} um, "
          f"v_A = {s['vA'] / 1e3:.1f} km/s, v_sh = {args.vsh_kms:.0f} km/s "
          f"({vsh / s['vA']:.2f} v_A), 1/w_ci0 = {1e9 / s['wci0']:.4f} ns")

    # fixed axes over the movie (99.9th percentile), smoothing over field_smooth_de
    loaded = {}
    for st in steps:
        fs = int(fsteps[max(0, np.searchsorted(fsteps, st, side="right") - 1)])
        if fs not in loaded:
            fr = kio.load_frame(fields[fs], add_applied_field=False, fields=("Ez", "Te", "rho_piston_ions",
                                                    "rho_amb_ions"))
            dz = fr.z_edges[1] - fr.z_edges[0]
            nk = max(1, int(round(args.field_smooth_de * de / dz)))
            ker = np.ones(nk) / nk

            def sm(a, ker=ker):
                return np.convolve(a, ker, mode="same") if len(ker) > 1 else a

            zero = np.zeros_like(fr.Bx)
            loaded[fs] = dict(
                z=fr.z_centers / s["di"], t=fr.time, bx=sm(fr.Bx / s["B0"]),
                ez=sm(fr.comps.get("Ez", zero) / EB), te=sm(fr.comps.get("Te", zero)),
                npi=sm(fr.comps.get("rho_piston_ions", zero) / (QE * s["namb"])),
                nai=sm(fr.comps.get("rho_amb_ions", zero) / (QE * s["namb"])))
    vals = loaded.values()
    blim = 1.1 * max(np.percentile(v["bx"], 99.9) for v in vals)
    bfloor = min(-0.05 * blim, 1.1 * min(np.percentile(v["bx"], 0.1) for v in vals))
    elim = 1.1 * max(np.percentile(np.abs(v["ez"]), 99.9) for v in vals)
    tlim = 1.1 * max(np.percentile(v["te"], 99.9) for v in vals)
    nlim = 1.2 * max(np.percentile(np.maximum(v["npi"], v["nai"]), 99.9) for v in vals)

    Lz = max(v["z"][-1] for v in vals) + 0.5 * (loaded[next(iter(loaded))]["z"][1]
                                                 - loaded[next(iter(loaded))]["z"][0])
    z_edges = np.linspace(0.0, Lz, 1601)
    v_edges = np.linspace(args.ion_vlim[0], args.ion_vlim[1], 241)
    framedir = os.path.join(os.environ.get("MOVIE_SCRATCH", os.path.join(run_dir, ".frames")),
                            run_id, "hybrid_phase")
    os.makedirs(framedir, exist_ok=True)
    for f in os.listdir(framedir):
        if f.startswith("frame_"):
            os.remove(os.path.join(framedir, f))

    for k, st in enumerate(steps):
        fr = kio.load_frame(phase[st], add_applied_field=False)  # B0 is on the grid here
        sd = {}
        for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
            z, uz, w = kio.species_phase_weighted(fr, sp, None, mass=s["Mi"])
            if len(z):
                sd[key] = (z / s["di"], uz * kio._C() / vsh, w)
        fs = int(fsteps[max(0, np.searchsorted(fsteps, st, side="right") - 1)])
        F = loaded[fs]
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 0.8, 0.7]})
        P.phase_distribution(axes[0], sd, z_edges, v_edges, vline=1.0, legend=True)
        axes[0].set_ylabel(r"$v_z / v_{sh}$")
        axes[0].text(0.985, 0.965, "ions (kinetic)", transform=axes[0].transAxes,
                     va="top", ha="right", fontsize=9, color="w")
        axf, axe = axes[1], axes[1].twinx()
        axe.plot(F["z"], F["ez"], color=P.C_AMBIENT, lw=0.6, alpha=0.8, label=r"$E_z$")
        axe.set_ylim(-elim, elim)
        axe.set_ylabel(r"$E_z / v_A B_0$")
        axf.plot(F["z"], F["bx"], color="k", lw=1.2, label=r"$B_x/B_0$")
        axf.set_ylim(bfloor, blim)
        axf.axhline(0.0, color="0.6", lw=0.6)
        axf.set_ylabel(r"$B_x / B_0$")
        axf.set_zorder(axe.get_zorder() + 1)
        axf.patch.set_visible(False)
        P.style_axes(axf)
        h1, l1 = axf.get_legend_handles_labels()
        h2, l2 = axe.get_legend_handles_labels()
        axf.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=2, frameon=False)
        note = f"running mean {args.field_smooth_de:g} $d_e$"
        if fs != st:
            note += rf"   fields at $t\,\omega_{{ci0}} = {F['t'] * s['wci0']:.3f}$"
        axf.text(0.01, 0.97, note, transform=axf.transAxes, va="top", ha="left",
                 fontsize=7, color="0.35")
        axt, axn = axes[2], axes[2].twinx()
        axn.plot(F["z"], F["npi"], color=P.C_PISTON, lw=0.8, label=r"$n_{pi}/n_0$")
        axn.plot(F["z"], F["nai"], color=P.C_AMBIENT, lw=0.8, label=r"$n_{ai}/n_0$")
        axn.set_yscale("log")
        axn.set_ylim(1e-2, max(nlim, 10))
        axn.set_ylabel(r"$n_i / n_{e0}$")
        axt.plot(F["z"], F["te"], color="k", lw=1.1, label=r"$T_e$")
        axt.set_ylim(0, tlim)
        axt.set_ylabel(r"$T_e$ [eV] (fluid)")
        axt.set_zorder(axn.get_zorder() + 1)
        axt.patch.set_visible(False)
        P.style_axes(axt)
        h1, l1 = axt.get_legend_handles_labels()
        h2, l2 = axn.get_legend_handles_labels()
        axt.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=3, frameon=False)
        axes[-1].set_xlabel(r"$z / d_{i0}$")
        axes[0].set_title(rf"{run_id} (hybrid: fluid e$^-$)   "
                          rf"$t\,\omega_{{ci0}} = {fr.time * s['wci0']:.3f}$", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{k:03d}.png"), dpi=110)
        plt.close(fig)

    out_dir = args.out or P.media_dir(run_id=run_id)
    os.makedirs(out_dir, exist_ok=True)
    return P.encode(framedir, os.path.join(out_dir, "hybrid_phase_fields.mp4"), fps=args.fps)


if __name__ == "__main__":
    main()
