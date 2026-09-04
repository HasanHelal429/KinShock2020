#!/usr/bin/env python3
"""Schaeffer 2020 Fig. 7 layout, for any run in this repo.

Three rows -- ambient-ion, piston-ion and electron v_z phase space, with |Bx|/B0 and
n_e/n_e0 overplotted on the electron row -- across five timesteps, in the paper's own
units (z/d_i0, v_z/v_sh). Built so a run can be put beside the published figure and read
the same way, feature for feature.

WHAT TO LOOK FOR, AND WHERE. Sec. III C is explicit that two ambient-ion populations must
not be confused, and that they are told apart by WHERE they sit:

  the C      swept-up ambient ions accelerated by the piston's AMBIPOLAR field, arcing off
             the undisturbed band and back. Forms with no magnetic field at all. Labelled
             I in the paper's Fig. 7.
  the notch  the bump at the front of the C's LOWER branch, where the accelerated
             population rejoins the undisturbed upstream. In a healthy run the branch rolls
             over smoothly through it; the reflected population that grows out of this
             region is the paper's II, and its appearance is t*_1.

In the 470 eV runs the notch instead steepens into a standing vertical cliff just ahead of
the piston -- a false shock. `ek_47keV_t2` rolls over smoothly at the same t*;
`R1_470eV_s3f8_t2` drops discontinuously from v_z ~ 0.9 to ~0.2 at z ~ 7 d_i0.

THE WINDOW IS SET IN rho_i0, NOT d_i0. The paper's panels run from about the piston front
minus 0.19 rho_i0 to plus 0.48 rho_i0 -- it just happens to look like d_i0 there because
E_phase pins rho_i0/d_i0 = 10.40. That ratio is NOT fixed across M_phase (10.4 -> 104), so
a window in d_i0 would frame a different fraction of the structure at every run. The axis
is still LABELLED in d_i0, as the paper's is.

    python scripts/plot_fig7_style.py runs/M_phase/mb_r100_t2
    python scripts/plot_fig7_style.py runs/E_phase/ek_47keV_t2 --times 0.15 0.49 0.73 0.98 1.25
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.colors import LogNorm, Normalize     # noqa: E402
import kinshock                                      # noqa: E402
from kinshock import io as kio                       # noqa: E402
from kinshock.units import C, ME                     # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
ACC, FIELD, DENS = "#ffd166", "#e8e6ef", "#ff6b57"
WIN_LO, WIN_HI = -0.22, 0.50        # window in rho_i0 about the piston front


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--times", type=float, nargs="+",
                    default=[0.15, 0.49, 0.73, 0.98, 1.25])
    ap.add_argument("--norm", choices=("log", "linear"), default="log",
                    help="colour stretch. log (default) matches the paper's own rendering "
                         "and shows the faint reflected population; linear shows where the "
                         "BULK of each population actually sits, at the cost of the tails. "
                         "The notch is a step in the DENSE lower branch, so linear reads it "
                         "more honestly than log does.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rid = os.path.basename(args.run.rstrip("/"))
    out = args.out or (f"media/comparison/fig7_{rid}"
                       f"{'_linear' if args.norm == 'linear' else ''}.png")

    import yt; yt.set_log_level(50)
    cfg = kinshock.load(args.run)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(args.run)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])

    ROWS = [("Ambient Ions", ["amb_ions"], mi, 2.0, -0.5),
            ("Piston Ions", ["piston_ions"], mi, 1.6, -0.3),
            ("Electrons", ["amb_electrons", "piston_electrons"], ME, 6.0, -6.0)]
    fig, ax = plt.subplots(3, len(args.times),
                           figsize=(3.5 * len(args.times), 8.4), dpi=170,
                           facecolor=SURFACE, squeeze=False)

    for c, T in enumerate(args.times):
        i = int(np.argmin(np.abs(times - T)))
        fr = kio.load_frame(pfs[i], fields=("Bx",))
        ad = fr.ds.all_data()
        zp = np.percentile(np.asarray(ad[("piston_ions", "particle_position_x")]),
                           99.5) / sc.di0
        rho = sc.rho_i0 / sc.di0
        lo, hi = zp + WIN_LO * rho, zp + WIN_HI * rho
        for r, (name, species, m, vhi, vlo) in enumerate(ROWS):
            a = ax[r][c]
            a.set_facecolor("#111014")
            zs, vs, ws = [], [], []
            for sp in species:
                try:
                    zz = np.asarray(ad[(sp, "particle_position_x")]) / sc.di0
                except Exception:
                    continue
                pp = np.asarray(ad[(sp, "particle_momentum_z")])
                g = np.sqrt(1 + (pp / (m * C)) ** 2)
                zs.append(zz); vs.append(pp / (m * g) / sc.vsh_model)
                ws.append(np.asarray(ad[(sp, "particle_weight")]))
            z = np.concatenate(zs); v = np.concatenate(vs); w = np.concatenate(ws)
            zb = np.linspace(lo, hi, 300); vb = np.linspace(vlo, vhi, 230)
            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            # Clip the low end at a percentile of the OCCUPIED bins. The electron
            # distribution spans several decades in one panel, and a bare LogNorm puts
            # almost every bin at the top of the colormap -- the row renders as a flat
            # white fill that hides the structure and both overlays.
            nz = H[H > 0]
            if args.norm == "linear":
                # Linear puts the weight where the PARTICLES are. The notch is a step in
                # the dense lower branch, so it reads more honestly here than under log,
                # which lifts the sparse tails until they compete with the branch itself.
                nrm = Normalize(vmin=0, vmax=(nz.max() if nz.size else 1))
            else:
                vmin = np.quantile(nz, 0.02) if nz.size else None
                nrm = LogNorm(vmin=vmin, vmax=(nz.max() if nz.size else None))
            a.pcolormesh(zb, vb, np.ma.masked_where(H <= 0, H).T, norm=nrm,
                         cmap="magma", shading="auto", zorder=0)
            a.axvline(zp, color=ACC, lw=1.1)
            a.set_xlim(lo, hi); a.set_ylim(vlo, vhi)
            a.tick_params(colors=MUTED, labelsize=8)
            for s in a.spines.values():
                s.set_color(BASE)
            if r == 2:                                   # Bx/B0 and n_e/n_e0 overlay
                zc = fr.z_centers / sc.di0
                Bx = fr.comps.get("Bx")
                m2 = (zc >= lo) & (zc <= hi)
                tw = a.twinx()
                tw.set_ylim(0, 12); tw.set_xlim(lo, hi)
                tw.tick_params(colors=MUTED, labelsize=7.5)
                for s in tw.spines.values():
                    s.set_color(BASE)
                if Bx is not None:
                    # SMOOTH before plotting, for two reasons. Physically, RESULTS
                    # 2026-08-11 requires a 0.2 d_i0 boxcar on this quantity -- a raw max
                    # runs 2.2-2.8x the smoothed value and that trap already produced one
                    # retracted result. Practically, a run with a small B0 and a wide
                    # window (mb_r100_t2: B0 = 7 T, 186k cells) has such large relative
                    # upstream fluctuations that the raw line paints the whole panel solid
                    # and hides the electron phase space underneath it.
                    zz, bb = zc[m2], np.abs(Bx[m2]) / sc.B0
                    if zz.size > 4:
                        dzc = float(np.median(np.diff(zz)))
                        kw = max(int(round(0.2 / dzc)), 1)
                        if kw > 1:
                            bb = np.convolve(bb, np.ones(kw) / kw, mode="same")
                    tw.plot(zz, bb, color=FIELD, lw=1.3, zorder=5,
                            label=r"$|B_x|/B_0$ (0.2 $d_{i0}$ boxcar)")
                # n_e/n_e0 straight from the electron weights, normalised to the far
                # upstream of THIS window -- the plotfiles carry no rho field
                ne, edges = np.histogram(z, bins=np.linspace(lo, hi, 240), weights=w)
                far = np.median(ne[-40:])
                if far > 0:
                    tw.plot(0.5 * (edges[1:] + edges[:-1]), ne / far,
                            color=DENS, lw=1.3, ls=(0, (4, 3)), zorder=5,
                            label=r"$n_e/n_{e0}$")
                if c == 0:
                    lg = tw.legend(loc="upper right", fontsize=7.5, framealpha=.25,
                                   facecolor="#111014", edgecolor=BASE)
                    for t_ in lg.get_texts():
                        t_.set_color(INK2)
                if c == len(args.times) - 1:
                    tw.set_ylabel(r"$|B_x|/B_0$ ,  $n_e/n_{e0}$", color=INK2, fontsize=9)
            if r == 0:
                a.set_title(rf"$t\,\omega_{{ci0}}$ = {times[i]:.2f}", color=INK, fontsize=11)
            if c == 0:
                a.set_ylabel(f"{name}\n" + r"$v_z/v_{sh}$", color=INK2, fontsize=9.5)
            if r == 2:
                a.set_xlabel(r"$z/d_{i0}$", color=INK2, fontsize=10)
        print(f"  column done: t* = {times[i]:.2f}", flush=True)

    wpe0 = np.sqrt(sc.namb * 1.602176634e-19 ** 2 / (8.8541878128e-12 * ME))
    fig.suptitle(f"{rid} — Fig. 7 layout ({args.norm})    "
                 rf"$\omega_{{pe0}}/\omega_{{ce}}$ = {wpe0/(1.602176634e-19*sc.B0/ME):.0f}    "
                 rf"$c/v_{{ti}}$ = {np.sqrt(sc.mass_ratio/(sc.Cs0/C)**2/sc.mass_ratio):.0f}    "
                 rf"$M_A$ = {sc.MA:.0f}",
                 color=INK, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
