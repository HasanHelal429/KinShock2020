#!/usr/bin/env python3
"""Side-by-side phase-space movie of TWO hybrid runs on SHARED, PHYSICAL axes.

    python scripts/hybrid_pair_movie.py runs/hybrid_phase/A runs/hybrid_phase/B \
        [--labels "A text" "B text"] [--out media/hybrid_phase/pair_A_vs_B.mp4]

`hybrid_phase_movie.py` normalises each run by its own B0, v_A and w_ci0, so two runs with
different B0 cannot be compared frame by frame from their separate movies (different time
labels, different axis scales). Here both columns share every axis, in physical units:
  * ions      (z, v_z) with v_z in units of ONE fixed speed (--vsh-kms, default 4182 km/s)
  * fields    B_x [T] (total, B0 is on the grid) and E_z [V/m]
  * electrons T_e [eV] and n_i / n_amb (log)
Frames are matched by STEP (the runs must share dt and the diag_phase cadence); the title
gives the physical time in ns and each run's own t*w_ci0. Panel limits are the 99.9th
percentile over BOTH runs and all frames. Frame rate = n_frames / --seconds (default 15 s).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

from kinshock import io as kio  # noqa: E402
from kinshock import plotting as P  # noqa: E402
from hybrid_phase_movie import QE, dumps, scales_from_used_inputs  # noqa: E402

ME = 9.1093837015e-31


def load_fields(run, s, steps, smooth_de):
    fields = dumps(run, "diag_fields")
    fsteps = np.array(sorted(fields))
    de = s["di"] / np.sqrt(s["Mi"] / ME)
    out = {}
    for st in steps:
        fs = int(fsteps[max(0, np.searchsorted(fsteps, st, side="right") - 1)])
        if fs in out:
            continue
        fr = kio.load_frame(fields[fs], add_applied_field=False,
                            fields=("Ez", "Te", "rho_piston_ions", "rho_amb_ions"))
        dz = fr.z_edges[1] - fr.z_edges[0]
        nk = max(1, int(round(smooth_de * de / dz)))
        ker = np.ones(nk) / nk
        sm = (lambda a: np.convolve(a, ker, mode="same")) if nk > 1 else (lambda a: a)
        zero = np.zeros_like(fr.Bx)
        out[fs] = dict(z=fr.z_centers / s["di"], t=fr.time, bx=sm(fr.Bx),
                       ez=sm(fr.comps.get("Ez", zero)), te=sm(fr.comps.get("Te", zero)),
                       npi=sm(fr.comps.get("rho_piston_ions", zero) / (QE * s["namb"])),
                       nai=sm(fr.comps.get("rho_amb_ions", zero) / (QE * s["namb"])))
    return fsteps, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    ap.add_argument("--labels", nargs=2, default=None)
    ap.add_argument("--vsh-kms", type=float, default=4182.0)
    ap.add_argument("--field-smooth-de", type=float, default=2.0)
    ap.add_argument("--ion-vlim", type=float, nargs=2, default=(-0.5, 2.5))
    ap.add_argument("--seconds", type=float, default=15.0, help="movie length; fps = frames/seconds")
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    runs = [os.path.abspath(r) for r in (args.run_a, args.run_b)]
    ids = [os.path.basename(r) for r in runs]
    labels = args.labels or ids
    sc = [scales_from_used_inputs(r) for r in runs]
    vsh = args.vsh_kms * 1e3
    steps = sorted(set(dumps(runs[0], "diag_phase")) & set(dumps(runs[1], "diag_phase")))
    if args.max_frames:
        steps = steps[:args.max_frames]
    if not steps:
        raise SystemExit("the two runs share no diag_phase steps")
    F = [load_fields(r, s, steps, args.field_smooth_de) for r, s in zip(runs, sc)]
    allv = [v for _, d in F for v in d.values()]
    # Axis limits skip the outer 2 d_i0: late in a run the open hi boundary throws a B/E spike
    # ~10x the real ramp (CLAUDE.md, B_compression gotcha), which would flatten every panel.
    inner = lambda v: v["z"] < v["z"][-1] - 2.0
    pct = lambda key, q, f=lambda a: a: [np.percentile(f(v[key][inner(v)]), q) for v in allv]
    blim = 1.1 * max(pct("bx", 99.9))
    bfloor = min(-0.05 * blim, 1.1 * min(pct("bx", 0.1)))
    elim = 1.1 * max(pct("ez", 99.9, np.abs))
    tlim = 1.1 * max(pct("te", 99.9))
    nlim = 1.2 * max(np.percentile(np.maximum(v["npi"], v["nai"])[inner(v)], 99.9) for v in allv)
    Lz = max(v["z"][-1] for v in allv)
    z_edges = np.linspace(0.0, Lz, 1601)
    v_edges = np.linspace(args.ion_vlim[0], args.ion_vlim[1], 241)
    phase = [dumps(r, "diag_phase") for r in runs]
    print(f"{len(steps)} shared frames; B0 = {sc[0]['B0']:.3g} / {sc[1]['B0']:.3g} T; "
          f"limits B {bfloor:.3g}..{blim:.3g} T, |E| {elim:.3g} V/m, T_e {tlim:.0f} eV")

    tag = f"pair_{ids[0]}_vs_{ids[1]}"
    framedir = os.path.join(os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, "media", ".frames")), tag)
    os.makedirs(framedir, exist_ok=True)
    for f in os.listdir(framedir):
        if f.startswith("frame_"):
            os.remove(os.path.join(framedir, f))

    for k, st in enumerate(steps):
        fig, axes = plt.subplots(3, 2, figsize=(15.0, 9.0), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 0.8, 0.7]})
        for c in range(2):
            s, (fsteps, loaded) = sc[c], F[c]
            fr = kio.load_frame(phase[c][st], add_applied_field=False)
            sd = {}
            for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
                z, uz, w = kio.species_phase_weighted(fr, sp, None, mass=s["Mi"])
                if len(z):
                    sd[key] = (z / s["di"], uz * kio._C() / vsh, w)
            fs = int(fsteps[max(0, np.searchsorted(fsteps, st, side="right") - 1)])
            D = loaded[fs]
            ax = axes[0, c]
            P.phase_distribution(ax, sd, z_edges, v_edges, vline=1.0, legend=(c == 0))
            ax.set_ylabel(rf"$v_z$ / {args.vsh_kms:.0f} km/s")
            ax.set_title(f"{labels[c]}\n" rf"$t\,\omega_{{ci0}}$ = {fr.time * s['wci0']:.2f}"
                         f"  (B0 = {s['B0']:.2f} T)", fontsize=10)
            axf, axe = axes[1, c], axes[1, c].twinx()
            axe.plot(D["z"], D["ez"] / 1e6, color=P.C_AMBIENT, lw=0.6, alpha=0.8, label=r"$E_z$")
            axe.set_ylim(-elim / 1e6, elim / 1e6)
            axe.set_ylabel(r"$E_z$ [MV/m]")
            axf.plot(D["z"], D["bx"], color="k", lw=1.2, label=r"$B_x$")
            axf.set_ylim(bfloor, blim)
            axf.axhline(0.0, color="0.6", lw=0.6)
            axf.set_ylabel(r"$B_x$ [T]")
            axf.set_zorder(axe.get_zorder() + 1); axf.patch.set_visible(False)
            P.style_axes(axf)
            h1, l1 = axf.get_legend_handles_labels(); h2, l2 = axe.get_legend_handles_labels()
            axf.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=2, frameon=False)
            axt, axn = axes[2, c], axes[2, c].twinx()
            axn.plot(D["z"], D["npi"], color=P.C_PISTON, lw=0.8, label=r"$n_{pi}/n_0$")
            axn.plot(D["z"], D["nai"], color=P.C_AMBIENT, lw=0.8, label=r"$n_{ai}/n_0$")
            axn.set_yscale("log"); axn.set_ylim(1e-2, max(nlim, 10))
            axn.set_ylabel(r"$n_i / n_{amb}$")
            axt.plot(D["z"], D["te"], color="k", lw=1.1, label=r"$T_e$")
            axt.set_ylim(0, tlim); axt.set_ylabel(r"$T_e$ [eV]")
            axt.set_zorder(axn.get_zorder() + 1); axt.patch.set_visible(False)
            P.style_axes(axt)
            h1, l1 = axt.get_legend_handles_labels(); h2, l2 = axn.get_legend_handles_labels()
            axt.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=3, frameon=False)
            axes[2, c].set_xlabel(r"$z / d_{i0}$")
            t_ns = fr.time * 1e9
        fig.suptitle(f"t = {t_ns:.3f} ns  (step {st})", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{k:03d}.png"), dpi=100)
        plt.close(fig)

    out = args.out or os.path.join(ROOT, "media", "hybrid_phase", f"{tag}.mp4")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fps = max(1, int(round(len(steps) / args.seconds)))
    res = P.encode(framedir, out, fps=fps)
    shutil.rmtree(framedir, ignore_errors=True)   # frames are intermediate (user, 2026-11-18)
    return res


if __name__ == "__main__":
    main()
