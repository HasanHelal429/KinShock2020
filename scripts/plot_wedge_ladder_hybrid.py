#!/usr/bin/env python3
"""The wedge ladder + hybrid comparison — five drivers/resolutions, two times.

    python scripts/plot_wedge_ladder_hybrid.py                 # linear colour (default)
    python scripts/plot_wedge_ladder_hybrid.py --norm log      # the original scaling

Recreates `media/comparison/wedge_ladder_hybrid.png`, whose generator was never committed
(the figure dates from 2026-08-14 and only RESULTS references it). Rebuilt here as a script
so it can be regenerated and so the wedge numbers cannot drift.

THE POINT OF THE FIGURE (RESULTS 2026-08-14). Two unrelated brakes on electron heating --
relativity in the 47 keV reference, and a polytropic fluid closure in the hybrid -- give the
SAME precursor to ~6% at every time, while unconstrained kinetic 470 eV electrons give a
different one, and give it ~1.9x faster. The resolution ladder in between (dz1 -> dz8 ->
dz16, an 16x span) shows the wedge is not a resolution artifact.

LINEAR COLOUR IS THE DEFAULT HERE. A log stretch lifts the sparse tails until they compete
visually with the dense core, which for this figure overstates the thin high-|v_z| arcs and
makes rows with very different particle counts look more alike than they are. Linear shows
where the ions actually are. `--norm log` restores the original if the tails are the point.

`wedge_depth` is IMPORTED from plot_mag_ladder rather than reimplemented, so this figure and
the M_phase ones can never disagree about what a wedge is: the first offset ahead of the
piston where the local v_z spread drops below 0.06 v_sh, linearly interpolated.

v_sh IS THE MODEL VALUE FOR EVERY ROW. Only R1_paper carries a shock_fit.yaml; using it for
that row alone would silently put one row on a fitted v_sh and four on the model one, which
is exactly the figure-vs-movie inconsistency CLAUDE.md warns about. Uniform beats accurate
here because the figure is a comparison.
"""
import argparse, importlib.util, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib                                   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.colors import Normalize, LogNorm    # noqa: E402
import kinshock                                     # noqa: E402
from kinshock import io as kio                      # noqa: E402
from kinshock.units import C, ME                    # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "_mag", os.path.join(ROOT, "scripts", "plot_mag_ladder.py"))
_mag = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mag)
wedge_depth = _mag.wedge_depth

SURFACE, INK, INK2, MUTED = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798"
BASE, PISTON, VSH, BAND = "#3a3745", "#ffd166", "#5ad2ff", "#5ad2ff"

# (run, row label, dz/lambda_D, scales borrowed from)
# The hybrid row borrows its scales. Its config is a DIFFERENT SCHEMA -- hybrid-solver keys
# (solver/operators/gates/acceptance), geometry in d_i not d_e, and `field.vA_over_c`, the
# legacy key `units.derive` refuses outright -- so derive() cannot read it at all. Borrowing
# is not a fudge here: its reference block gives n_amb = 4.8e24 and m_i/m_e = 100, IDENTICAL
# to the 470 eV runs, so d_i0 is the same number; and its own targets (M_A 14, M_ms 13,
# beta_0 0.2, vp_over_c 0.0103986833) are the same Table I point, with
# M_A x v_A/c = 0.0139982 against the sibling's v_sh = 0.0139507 -- 0.34% apart. Putting the
# row on its own 0.34%-different v_sh would make it less comparable, not more.
ROWS = [
    ("runs/R1_phase/R1_paper",             "R1_paper  (47 keV ref)",  "0.60", None),
    ("runs/S_phase/ss_dz1_ppc100",         "dz1  production",         "6.07", None),
    ("runs/S_phase/ss_dz8_ppc100_s3f8",    "dz8 + shape3 + filter8",  "0.76", None),
    ("runs/S_phase/ss_dz16_ppc100",        "dz16  finest",            "0.38", None),
    ("runs/hybrid_phase/H3_470eV_dense",   "HYBRID  (fluid e-)",      "n/a",
     "runs/S_phase/ss_dz1_ppc100"),
]
VLO, VHI = -1.0, 2.6


def panel_data(run, t_target, zmax, scales_from=None):
    if scales_from:
        sc = kinshock.units.derive(kinshock.load(scales_from))
        mi = 100.0 * ME                       # hybrid reference.ion_mass_over_me
    else:
        cfg = kinshock.load(run)
        sc = kinshock.units.derive(cfg)
        mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(run)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    i = int(np.argmin(np.abs(times - t_target)))
    ad = kio.load_frame(pfs[i]).ds.all_data()
    z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
    pm = np.asarray(ad[("amb_ions", "particle_momentum_z")])
    g = np.sqrt(1 + (pm / (mi * C)) ** 2)
    v = pm / (mi * g) / sc.vsh_model
    w = np.asarray(ad[("amb_ions", "particle_weight")])
    zp = float(np.percentile(
        np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5)) / sc.di0
    rho = sc.rho_i0 / sc.di0
    box = sc.domain_halfwidth / sc.di0
    # wedge_depth works in rho_i0 about the piston, and scans to the end of the box
    d_rho = wedge_depth(z / rho, v, zp / rho, scan=max((box - zp) / rho, 0.1))
    return times[i], z, v, w, zp, (d_rho * rho if np.isfinite(d_rho) else np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--norm", choices=("linear", "log"), default="linear")
    ap.add_argument("--times", type=float, nargs="*", default=[0.13, 0.26])
    ap.add_argument("--zmax", type=float, default=13.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import yt; yt.set_log_level(50)
    nb_z, nb_v = 460, 210
    zb = np.linspace(0, args.zmax, nb_z + 1)
    vb = np.linspace(VLO, VHI, nb_v + 1)

    fig, ax = plt.subplots(len(ROWS), len(args.times),
                           figsize=(11.5, 3.05 * len(ROWS)), dpi=125,
                           facecolor=SURFACE, sharex=True, squeeze=False)
    print(f"{'run':26}{'t*':>7}{'z_piston':>10}{'wedge d_i0':>12}")
    for r, (run, label, dzl, borrow) in enumerate(ROWS):
        for c, T in enumerate(args.times):
            a = ax[r][c]; a.set_facecolor("#111014")
            try:
                t, z, v, w, zp, wd = panel_data(run, T, args.zmax, borrow)
            except Exception as e:
                a.annotate(f"{type(e).__name__}", xy=(.5, .5), xycoords="axes fraction",
                           ha="center", color=MUTED, fontsize=9)
                print(f"{os.path.basename(run):26}{T:7.2f}   ERROR {e}")
                continue
            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            nz = H[H > 0]
            vmax = nz.max() if nz.size else 1.0
            nrm = (Normalize(0, vmax) if args.norm == "linear"
                   else LogNorm(vmin=max(vmax * 1e-4, 1e-30), vmax=vmax))
            a.pcolormesh(0.5 * (zb[1:] + zb[:-1]), 0.5 * (vb[1:] + vb[:-1]),
                         np.ma.masked_where(H <= 0, H).T, norm=nrm,
                         cmap="magma", shading="auto", zorder=0)
            a.axhline(1.0, color=VSH, lw=1.1, ls=(0, (6, 4)), alpha=.9, zorder=3)
            a.axhline(0.0, color="#ffffff", lw=.5, alpha=.25, zorder=2)
            a.axvline(zp, color=PISTON, lw=1.8, zorder=4)
            if np.isfinite(wd):
                a.axvspan(zp, zp + wd, color=BAND, alpha=0.13, zorder=1)
                a.annotate(rf"wedge {wd:.2f} $d_{{i0}}$", xy=(0.975, 0.06),
                           xycoords="axes fraction", ha="right", color=VSH, fontsize=9.5)
            else:
                a.annotate("wedge: no crossing", xy=(0.975, 0.06),
                           xycoords="axes fraction", ha="right", color=MUTED, fontsize=9.5)
            a.set_xlim(0, args.zmax); a.set_ylim(VLO, VHI)
            a.tick_params(colors=MUTED, labelsize=8.5)
            for s in a.spines.values():
                s.set_color(BASE)
            if r == 0:
                a.set_title(rf"$t\,\omega_{{ci0}}$ = {T:.2f}", color=INK, fontsize=12,
                            pad=10)
            if c == 0:
                a.set_ylabel(f"{label}\n" + r"$dz/\lambda_D$ = " + dzl
                             + "\n\n" + r"$v_z/v_{sh}$", color=INK2, fontsize=9)
            if r == len(ROWS) - 1:
                a.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=11)
            print(f"{os.path.basename(run):26}{t:7.3f}{zp:10.2f}"
                  f"{wd if np.isfinite(wd) else float('nan'):12.2f}", flush=True)

    fig.suptitle("Ambient-ion phase space: the wedge survives the full resolution ladder\n"
                 f"each row in its own $d_{{i0}}$ and MODEL $v_{{sh}}$  |  "
                 f"yellow = piston front, dashed = $v_{{sh}}$, shaded = wedge  |  "
                 f"{args.norm} colour, per panel",
                 color=INK, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = args.out or f"media/comparison/wedge_ladder_hybrid_{args.norm}.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
