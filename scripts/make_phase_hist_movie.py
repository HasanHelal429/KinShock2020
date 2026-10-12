#!/usr/bin/env python3
"""Phase-space movie from the BINNED diagnostic (WarpX ParticleHistogram2D, PS_*_zuz).

    python scripts/make_phase_hist_movie.py <run_dir> [--fps 15] [--norm asinh|linear]
                                            [--fixed-scale] [--vsh-c 0.0326] [--out DIR]
                                            [--fields]

Panels on a shared z axis, one frame per histogram dump (deck.render writes them at
~10x the raw-particle cadence, so this is the smooth-in-time movie; make_movies.py's
shock_phase.mp4 is the particle-dump one):

  * ions      piston + ambient, v_z / v_sh       (dotted: v_z = v_sh, reflected-ion line)
  * electrons piston + ambient, v_z / v_te,ab
  * fields    (--fields) B_x/B0 on the left axis, E_z and E_y / (v_A B0) on the right,
              from the diag_fields plotfile at the same step. B0 is along x, so B_x is
              the compressed component. Smoothed by a running mean (--field-smooth-de,
              default 10 d_e ~ d_i,ab). Axis limits are fixed over the whole movie
              (99.9th percentile of every frame), so growth reads as growth.
              E_y stays noise-dominated upstream: it is the E half of the light-wave
              noise that also wiggles B_x, and in v_A B0 units it carries c/v_A (~4e3
              at 470 eV), so delta B_x/B0 ~ 0.1 reads as E_y ~ 300 v_A B0.

The rendering is plotting.phase_image -- the same stretch, per-species normalisation and
two-colour tint as shock_phase.mp4 -- so the two movies read the same way. Default is
asinh normalised per species per frame (make_movies' convention); --fixed-scale uses each
species' peak over ALL frames instead, so brightness is comparable across time.

The histogram ordinate is u = gamma*beta. Where a panel's window stays below u = 0.1
u*c is v to <0.5% and the axis says v_z; wider windows (relativistic electrons) say gamma*v_z.

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
from kinshock import io as kio  # noqa: E402
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
    ap.add_argument("--fields", action="store_true",
                    help="add a B_x / E_z / E_y panel from the diag_fields plotfiles")
    ap.add_argument("--field-smooth-de", type=float, default=10.0,
                    help="running-mean width for the field panel [d_e] (0 = raw)")
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

    fields = {}
    if args.fields:
        by_step = {kio._step_key(p)[0]: p for p in kio.field_plotfiles(args.run_dir)}
        fsteps = np.array(sorted(by_step))
        # Field dumps can be sparser than histogram dumps (PoC: every 630 vs 126 steps).
        # Every phase-space frame is kept; its field panel shows the latest field dump at
        # or before that step (the first one before any exist), labelled with its own time.
        field_of = {i: int(fsteps[max(0, np.searchsorted(fsteps, i, side="right") - 1)])
                    for i in its}
        n_held = sum(1 for i in its if field_of[i] != i)
        if n_held:
            print(f"  {n_held}/{len(its)} frames show the latest earlier field dump "
                  f"(field cadence sparser than histogram cadence)")
        EB = sc.vA * sc.B0
        # Running mean over field_smooth_de: at Debye-resolved dz (~0.01 d_e) the raw
        # cell-scale E noise is ~1e4 v_A B0 and buries the shock. 2 d_e (the histogram
        # z bin) still leaves E_z noise comparable to the ramp spike; the default 10 d_e
        # (~ d_i,ab, the smallest ion scale) keeps the ramp and quiets the upstream.
        nk = max(1, int(round(args.field_smooth_de * sc.de / sc.dz)))
        ker = np.ones(nk) / nk

        def sm(a):
            return np.convolve(a, ker, mode="same") if nk > 1 else a

        # A particle-applied B0(z) is not on the grid. Branches whose io.load_frame knows
        # about it (io._applied_bx) return the TOTAL B_x; the applied profile is also
        # drawn dashed so background and response can be told apart.
        applied = getattr(kio, "_applied_bx", lambda path: None)
        bg_line = None
        ftime = {}
        for s in sorted(set(field_of.values())):
            fr = kio.load_frame(by_step[s], fields=("Ey", "Ez"))
            zero = np.zeros_like(fr.Bx)
            fields[s] = (fr.z_centers / sc.di0, sm(fr.Bx / sc.B0),
                         sm(fr.comps.get("Ez", zero) / EB),
                         sm(fr.comps.get("Ey", zero) / EB))
            ftime[s] = fr.time
            if bg_line is None and applied(by_step[s]) is not None:
                bg_line = applied(by_step[s])(fr.z_centers) / sc.B0
        blim = 1.1 * max(np.percentile(f[1], 99.9) for f in fields.values())
        # B_x reverses sign behind the piston front, so the floor comes from the data too
        bfloor = min(-0.05 * blim, 1.1 * min(np.percentile(f[1], 0.1) for f in fields.values()))
        elim = 1.1 * max(np.percentile(np.abs(np.r_[f[2], f[3]]), 99.9)
                         for f in fields.values())
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
        if fields:
            fig, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True,
                                     gridspec_kw={"height_ratios": [1, 1, 0.8]})
        else:
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
            # The ordinate is u = gamma*beta. Where the histogram window reaches u > 0.1
            # (gamma - 1 > 0.5%, e.g. PoC electrons to u = 1.54) say gamma*v, not v.
            ax.set_ylabel(ylab.replace("v_z", r"\gamma v_z")
                          if max(abs(ext[2]), abs(ext[3])) > 0.1 else ylab)
            ax.text(0.985, 0.965, title, transform=ax.transAxes, va="top", ha="right",
                    fontsize=9, color="w")
        if fields:
            zc, bx, ez, ey = fields[field_of[it]]
            axf = axes[2]
            axe = axf.twinx()
            axe.plot(zc, ez, color=P.C_AMBIENT, lw=0.6, alpha=0.8, label=r"$E_z$")
            axe.plot(zc, ey, color=P.C_PISTON, lw=0.6, alpha=0.8, label=r"$E_y$")
            axe.set_ylim(-elim, elim)
            axe.set_ylabel(r"$E / v_A B_0$")
            axf.plot(zc, bx, color="k", lw=1.2, label=r"$B_x/B_0$")
            if bg_line is not None:
                axf.plot(zc, bg_line, color="0.45", lw=1.0, ls="--",
                         label=r"applied $B_0(z)/B_0$")
            axf.set_ylim(bfloor, blim)
            axf.axhline(0.0, color="0.6", lw=0.6)
            axf.set_ylabel(r"$B_x / B_0$")
            axf.set_zorder(axe.get_zorder() + 1)       # B line drawn over the E traces
            axf.patch.set_visible(False)
            P.style_axes(axf)
            h1, l1 = axf.get_legend_handles_labels()
            h2, l2 = axe.get_legend_handles_labels()
            axf.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=3,
                       frameon=False)
            note = []
            if args.field_smooth_de > 0:
                note.append(f"running mean {args.field_smooth_de:g} $d_e$")
            if field_of[it] != it:
                note.append(rf"fields at $t\,\omega_{{ci0}} = {ftime[field_of[it]] * sc.wci0:.3f}$")
            if note:
                axf.text(0.01, 0.97, "   ".join(note), transform=axf.transAxes, va="top",
                         ha="left", fontsize=7, color="0.35")
        axes[-1].set_xlabel(r"$z / d_{i0}$")
        axes[0].set_title(rf"{run_id}   $t\,\omega_{{ci0}} = {t * sc.wci0:.3f}$", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{k:03d}.png"), dpi=110)
        plt.close(fig)

    out_dir = args.out or P.media_dir(run_id=run_id)
    os.makedirs(out_dir, exist_ok=True)
    name = (f"phase_hist_{args.norm}{'_fixed' if args.fixed_scale else ''}"
            f"{'_fields' if fields else ''}.mp4")
    return P.encode(framedir, os.path.join(out_dir, name), fps=args.fps)


if __name__ == "__main__":
    main()
