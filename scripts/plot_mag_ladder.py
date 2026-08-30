#!/usr/bin/env python3
"""M_phase summary figure: does the ambient wedge depend on w_pe0/w_ce at FIXED eps?

The E_phase analogue (`plot_eps_ladder.py`) with the scan axis swapped. Rows are the B0
scan ordered by `w_pe0/w_ce`; rung 1 IS `runs/E_phase/es_47keV`, which is why it is read
from E_phase rather than duplicated. Columns are matched `t*w_ci0`.

WHY THIS FIGURE CARRIES TWO EXTRA DIAGNOSTICS. As submitted, M_phase rungs 2-3 produced an
INVALID measurement and this plot has to show that rather than hide it (RESULTS /
runs/M_phase/README.md, 2026-08-24). `w_ci0 ~ B0`, so holding `t*w_ci0` fixed stretches
ABSOLUTE time as `1/B0` while `v_p` and `d_i0` stay frozen -- the piston sweeps `1/B0`
further through a box left at a fixed 12.02 d_i0. So each panel also draws

  * the DOMAIN EDGE, and the far-upstream control window (5-11 d_i0) the wedge comparison
    depends on staying cold, and
  * `T_e,upstream` in that window, in the run's own `C_s,ab^2`.

A panel whose control window is hot is stamped INVALID and its wedge number is drawn
struck-through, because with a contaminated upstream the `v_z` spread threshold is crossed
far downstream of the real wedge tip -- the number is an artifact, not a weak signal.

Frame matching keeps E_phase's tolerance: a panel whose nearest frame is further than TOL
from the request is left blank and reported rather than quietly mislabelled.

TWO WAYS TO MATCH THE PANELS, AND THEY DISAGREE.  `--match clock` (default) matches
`t*w_ci0`, the gyro clock shock formation runs on; that is the designed comparison and it
is the one the box failure breaks.  `--match piston` instead picks, per run, the frame
whose piston front sits nearest the anchor's -- same geometry, same box occupancy, both
upstreams cold -- at the cost of a `t*` that differs by the same `1/B0`, so the higher
rungs have had proportionally fewer gyroperiods to form a shock.

Neither is clean, and the wedge is reported in BOTH `d_i0` and `rho_i0` because at fixed
eps they decouple (`rho_i0/d_i0 = v_p/v_A ~ 1/B0`) and **they give opposite answers**.
E_phase never had to choose: it pinned `rho_i0/d_i0` = 10.40.

    python scripts/plot_mag_ladder.py [--match clock|piston] [--times ...] [--out PATH]
"""
import argparse
import os
import sys

import numpy as np

ROOT = "/pscratch/sd/h/hhelal/KinShock2020"
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock.units import C, ME  # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
GOOD, WARN, ACC = "#5ad2ff", "#ff6b57", "#ffd166"
ZMAX, VMAX = 13.5, 2.6
TOL = 0.030
UPSTREAM = (5.0, 11.0)     # the absolute control window, as eps_shape_test.py uses
T_UP_COLD = 5.0            # T_e,upstream above this => the comparison is void
SPREAD, STEP = 0.06, 0.02

# (label, run dir, B0 [T], w_pe0/w_ce, M_A, beta_0, rho_i0/d_i0)
ROWS = [
    ("rung 1  es_47keV", "runs/E_phase/es_47keV",  70.273,  10.0,  13.95, 0.196, 10.40),
    ("rung 2  mb_r20",   "runs/M_phase/mb_r20",    35.137,  20.0,  27.9,  0.783, 20.80),
    ("rung 3  mb_r32",   "runs/M_phase/mb_r32",    22.222,  31.6,  44.1,  1.957, 32.89),
    ("rung 4  mb_r56",   "runs/M_phase/mb_r56",    12.497,  56.2,  78.5,  6.190, 58.48),
    ("rung 5  mb_r100",  "runs/M_phase/mb_r100",    7.027, 100.0, 139.5, 19.57, 104.0),
]


def wedge_depth(z, v, zp):
    """E_phase's continuous depth: first 0.02 d_i0 offset where the local v_z spread
    drops below SPREAD, linearly interpolated between bracketing samples."""
    offs = np.arange(0, 6, STEP)
    sig = np.full(offs.size, np.nan)
    for k, lo in enumerate(offs):
        m = (z - zp >= lo) & (z - zp < lo + 0.25)
        if m.sum() > 100:
            sig[k] = v[m].std()
    ok = np.isfinite(sig)
    for k in range(1, offs.size):
        if ok[k] and ok[k - 1] and sig[k] < SPREAD <= sig[k - 1]:
            f = (sig[k - 1] - SPREAD) / (sig[k - 1] - sig[k])
            return offs[k - 1] + f * STEP
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--times", type=float, nargs="+", default=[0.130, 0.260])
    ap.add_argument("--match", choices=("clock", "piston"), default="clock",
                    help="clock: matched t*w_ci0. piston: matched piston front position "
                         "(taken from the anchor at the LAST requested time).")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.out is None:
        args.out = (f"media/M_phase/wedge_vs_B0"
                    f"{'_matchpiston' if args.match == 'piston' else ''}.png")

    rows = [r for r in ROWS if os.path.isdir(os.path.join(r[1], "diags"))]
    missing = [r[0] for r in ROWS if r not in rows]
    if missing:
        print(f"  not run yet, omitted: {', '.join(missing)}")
    if not rows:
        print("  nothing to plot"); return 1

    z_anchor = None
    if args.match == "piston":
        d0 = rows[0][1]
        sc0 = kinshock.units.derive(kinshock.load(d0))
        pfs0 = kio.plotfiles(d0)
        t0 = np.array([kio.load_frame(p).time * sc0.wci0 for p in pfs0])
        j = int(np.argmin(np.abs(t0 - args.times[-1])))
        ad0 = kio.load_frame(pfs0[j]).ds.all_data()
        z_anchor = np.percentile(
            np.asarray(ad0[("piston_ions", "particle_position_x")]), 99.5) / sc0.di0
        print(f"  matching on piston front z = {z_anchor:.3f} d_i0 "
              f"({rows[0][0]} at t*={t0[j]:.3f})")

    zb = np.linspace(0, ZMAX, 340)
    vb = np.linspace(-VMAX * 0.4, VMAX, 260)
    cols = [args.times[-1]] if args.match == "piston" else args.times
    fig, axes = plt.subplots(len(rows), len(cols),
                             figsize=(6.8 if len(cols) == 1 else 11.5,
                                      3.6 * len(rows) + 1.2), dpi=200,
                             facecolor=SURFACE, sharex=True, sharey=True, squeeze=False)
    table = []
    for r, (label, d, B0, wpwc, MA, b0, rho) in enumerate(rows):
        cfg = kinshock.load(d)
        sc = kinshock.units.derive(cfg)
        mi = cfg["reference"]["mass_ratio"] * ME
        box = sc.domain_halfwidth / sc.di0
        pfs = kio.plotfiles(d)
        times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
        fronts = None
        if args.match == "piston":
            fronts = np.array([
                np.percentile(np.asarray(kio.load_frame(p).ds.all_data()[
                    ("piston_ions", "particle_position_x")]), 99.5) / sc.di0
                for p in pfs])
        for c, T in enumerate(cols):
            ax = axes[r][c]
            ax.set_facecolor("#111014")
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_color(BASE)
            ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
            for t in ax.get_xticklabels() + ax.get_yticklabels():
                t.set_color(INK2)
            ax.set_xlim(0, ZMAX); ax.set_ylim(vb[0], vb[-1])
            if r == len(rows) - 1 or len(cols) == 1 and r == len(rows) - 1:
                ax.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=11)

            if args.match == "piston":
                i = int(np.argmin(np.abs(fronts - z_anchor)))
            else:
                i = int(np.argmin(np.abs(times - T)))
            if args.match == "clock" and abs(times[i] - T) > TOL:
                ax.annotate(f"no frame within {TOL:.3f}\nof $t^*$={T:.3f}\n"
                            f"(nearest {times[i]:.3f})", xy=(0.5, 0.5),
                            xycoords="axes fraction", ha="center", va="center",
                            color=MUTED, fontsize=9)
                print(f"  {label} @ t*={T}: REFUSED, nearest {times[i]:.4f}")
                continue

            ad = kio.load_frame(pfs[i]).ds.all_data()
            z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
            p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
            w = np.asarray(ad[("amb_ions", "particle_weight")])
            g = np.sqrt(1.0 + (p / (mi * C)) ** 2)
            v = p / (mi * g) / sc.vsh_model
            zp = np.percentile(
                np.asarray(ad[("piston_ions", "particle_position_x")]), 99.5) / sc.di0
            depth = wedge_depth(z, v, zp)

            # far-upstream control, in the run's own C_s,ab^2 -- the validity gate
            ze = np.asarray(ad[("amb_electrons", "particle_position_x")]) / sc.di0
            pe = np.stack([np.asarray(ad[("amb_electrons", f"particle_momentum_{a}")])
                           for a in "xyz"])
            mu_ = (ze > UPSTREAM[0]) & (ze < UPSTREAM[1])
            if mu_.sum() > 100:
                p2 = (pe[:, mu_] ** 2).sum(0)
                gu = np.sqrt(1.0 + p2 / (ME * C) ** 2)
                T_up = float(np.var(np.sqrt(p2) / (ME * gu)) / sc.Cs_ab ** 2)
            else:
                T_up = float("nan")
            valid = np.isfinite(T_up) and T_up < T_UP_COLD

            H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
            ax.pcolormesh(zb, vb, np.ma.masked_where(H <= 0, H).T, norm=LogNorm(),
                          cmap="magma", shading="auto")
            ax.axhline(0, color="#fff", lw=.5, alpha=.25)
            ax.axhline(1, color=GOOD, lw=.9, ls=(0, (5, 4)), alpha=.9)
            ax.axvline(zp, color=ACC, lw=1.2)
            # the two diagnostics this figure exists to show
            ax.axvspan(*UPSTREAM, color=GOOD if valid else WARN, alpha=0.05)
            ax.axvline(box, color=WARN, lw=1.4, ls=(0, (2, 2)))
            ax.annotate("domain edge", xy=(box - 0.12, vb[-1] * 0.93), ha="right",
                        color=WARN, fontsize=7.5, rotation=90, va="top")
            if np.isfinite(depth):
                ax.axvspan(zp, zp + depth, color=GOOD, alpha=0.07)
            txt = (f"wedge {depth:.2f} $d_{{i0}}$ = {depth / rho:.4f} " r"$\rho_{i0}$"
                   if np.isfinite(depth) else "wedge: no crossing")
            ax.annotate(txt, xy=(0.975, 0.06), xycoords="axes fraction", ha="right",
                        fontsize=9, color=GOOD if valid else MUTED,
                        bbox=None if valid else dict(boxstyle="round,pad=0.15", fc="none",
                                                     ec=WARN, lw=0.8))
            ax.annotate(rf"$z_{{piston}}$={zp:.2f}   $T_{{e,up}}$={T_up:.1f}",
                        xy=(0.02, 0.06), xycoords="axes fraction", fontsize=8.5,
                        color=INK2 if valid else WARN,
                        bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none",
                                  alpha=0.75))
            if not valid:
                ax.annotate("INVALID — upstream contaminated", xy=(0.5, 0.955),
                            xycoords="axes fraction", ha="center", fontsize=9.5,
                            color=WARN, weight="bold")
            if r == 0:
                ax.set_title(rf"$t\,\omega_{{ci0}}$ = {times[i]:.3f}"
                             + ("   (matched piston front)" if args.match == "piston"
                                else ""), color=INK, fontsize=12)
            elif args.match == "piston":
                ax.set_title(rf"$t\,\omega_{{ci0}}$ = {times[i]:.3f}", color=MUTED,
                             fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{label}\n" + rf"$\omega_{{pe0}}/\omega_{{ce}}$ = {wpwc:.0f}"
                              "\n" rf"$B_0$={B0:.1f} T   $M_A$={MA:.0f}" "\n"
                              rf"$\beta_0$={b0:.2f}   $\rho_{{i0}}/d_{{i0}}$={rho:.1f}"
                              "\n" r"$v_z/v_{sh}$", color=INK2, fontsize=9)
            table.append((label, wpwc, MA, b0, T, times[i], zp, depth,
                          depth / rho, T_up, valid))
        print(f"  row done: {label}", flush=True)

    held = (r"HELD: $\varepsilon$, $\gamma$=1.138, $T_{e,ab}$, $T_0$, $n_0$," "\n"
            r"$\lambda_{D,0}/d_{i0}$=4.42e-3, $dz/\lambda_D$=0.379, $N_D$=264"
            if len(cols) == 1 else
            r"every row HELD: $\varepsilon$, $\gamma$=1.138, $T_{e,ab}$, $T_0$, $n_0$, "
            r"$\lambda_{D,0}/d_{i0}$=4.42e-3, $dz/\lambda_D$=0.379, $N_D$=264")
    fig.suptitle(r"Ambient-ion phase space vs $\omega_{pe0}/\omega_{ce}$" "\n"
                 r"at FIXED $\varepsilon$ = 0.3033" "\n" + held,
                 color=INK, fontsize=11 if len(cols) == 1 else 12.5)
    fig.tight_layout(rect=[0, 0, 1, 1 - (0.075 if len(cols) == 1 else 0.045)
                           * (5.0 / len(rows))])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, facecolor=SURFACE)
    print(f"\n  wrote {args.out}")

    print("\n" + "=" * 96)
    print(f"{'run':20}{'wpe/wce':>9}{'M_A':>7}{'beta_0':>8}{'t* got':>8}"
          f"{'z_pist':>8}{'wedge d_i0':>12}{'wedge rho_i0':>14}{'T_e,up':>9}{'valid':>7}")
    print("=" * 96)
    for t in table:
        print(f"{t[0]:20}{t[1]:>9.1f}{t[2]:>7.0f}{t[3]:>8.2f}{t[5]:>8.3f}"
              f"{t[6]:>8.2f}{t[7]:>12.3f}{t[8]:>14.4f}{t[9]:>9.1f}"
              f"{'yes' if t[10] else 'NO':>7}")
    print("=" * 96)
    ok = [t for t in table if t[10]]
    if len(ok) >= 2:
        a, b = ok[0], ok[-1]
        print(f"\n  valid points only, {a[0].strip()} -> {b[0].strip()} "
              f"(w_pe0/w_ce {a[1]:.0f} -> {b[1]:.0f}):")
        print(f"    in d_i0   : {a[7]:.3f} -> {b[7]:.3f}   ({b[7]/a[7]-1:+.0%} change)")
        print(f"    in rho_i0 : {a[8]:.4f} -> {b[8]:.4f}   ({b[8]/a[8]-1:+.0%} change)")
        if (b[7] - a[7]) * (b[8] - a[8]) < 0:
            print("    -> THE TWO NORMALIZATIONS DISAGREE IN SIGN. rho_i0/d_i0 moves with")
            print("       B0, so the unit choice decides the answer. Settle it first.")
    bad = [t for t in table if not t[10]]
    if bad:
        print(f"\n  {len(bad)} of {len(table)} panels INVALID: the far-upstream control "
              f"({UPSTREAM[0]:g}-{UPSTREAM[1]:g} d_i0) is above {T_UP_COLD:g}, so the wedge")
        print("  numbers there are artifacts of a piston that has swept the box, not "
              "measurements.\n  See runs/M_phase/README.md for the cause and the fix.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
