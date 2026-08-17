#!/usr/bin/env python3
"""E_phase summary figure: does the ambient wedge depend on eps = v_te,ab/c?

Rows  = the five ladder rungs, ordered by eps (470 eV -> 47 keV), plus R1_paper as the
        independent 47 keV reference (different domain and sampling, same physics point).
Cols  = matched t*wci0 inside the S_phase window.

Each panel is the ambient-ion phase space in that run's OWN v_sh and d_i0, so the axes are
the same dimensionless picture at every rung and only eps moves. Wedge depth is MEASURED
per panel by the same rule the resolution ladder used (first 0.25 d_i0 bin ahead of the
piston front whose v_z spread falls below 0.06 v_sh), and printed as a table at the end --
the curve's SHAPE is the result:

    flat-then-drop  -> relativistic capping of electron heating (~eps^2)
    straight decade -> electron-scale wave regime, rho_e/lambda_D = w_pe/w_ce (~eps^-1)
    no trend        -> not an eps effect; the difference is elsewhere

FRAME MATCHING HAS A TOLERANCE. The earlier comparison scripts used a bare argmin and
silently rendered the SAME frame for two different requested times when a run's sampling
was sparse. A panel whose nearest frame is further than TOL from the request is left blank
and reported, rather than quietly mislabelled.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock.units import C, ME  # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
ZMAX, VMAX = 13.5, 2.6
TIMES = [0.130, 0.260]
TOL = 0.030          # max |t_frame - t_request| in 1/w_ci0 before a panel is refused

# (label, run dir, T_e,ab [eV], eps, w_pe0/w_ce, gamma(T_e,ab))
ROWS = [
    ("470 eV   (anchor)", "runs/S_phase/ss_dz16_ppc100", 470, 0.0303, 100.0, 1.0014),
    ("1.5 keV",           "runs/E_phase/es_1p5keV",     1486, 0.0539,  56.2, 1.0044),
    ("4.7 keV",           "runs/E_phase/es_4p7keV",     4700, 0.0959,  31.6, 1.0138),
    ("15 keV",            "runs/E_phase/es_15keV",     14860, 0.1705,  17.8, 1.0436),
    ("47 keV",            "runs/E_phase/es_47keV",     47012, 0.3033,  10.0, 1.1380),
    ("R1_paper (47 keV)", "runs/R1_phase/R1_paper",    47012, 0.3033,  10.0, 1.1380),
]

zb = np.linspace(0, ZMAX, 340)
vb = np.linspace(-VMAX * 0.4, VMAX, 260)
fig, axes = plt.subplots(len(ROWS), len(TIMES), figsize=(11.5, 18.5), dpi=200,
                         facecolor=SURFACE, sharex=True, sharey=True)
table = []
for r, (label, d, TeV, eps, wpwc, gam) in enumerate(ROWS):
    if not os.path.isdir(os.path.join(d, "diags")):
        print(f"  SKIP {label}: no diags/ yet", flush=True)
        for c in range(len(TIMES)):
            axes[r, c].set_facecolor("#111014")
        continue
    cfg = kinshock.load(d)
    sc = kinshock.units.derive(cfg)
    vsh = sc.vsh_model
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(d)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    for c, T in enumerate(TIMES):
        ax = axes[r, c]
        ax.set_facecolor("#111014")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(BASE)
        ax.tick_params(colors=MUTED, labelsize=8.5, length=3)
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            t.set_color(INK2)
        ax.set_xlim(0, ZMAX)
        ax.set_ylim(vb[0], vb[-1])
        if r == len(ROWS) - 1:
            ax.set_xlabel(r"$z / d_{i0}$", color=INK2, fontsize=11)

        i = int(np.argmin(np.abs(times - T)))
        if abs(times[i] - T) > TOL:
            ax.annotate(f"no frame within {TOL:.3f}\nof $t^*$={T:.3f}\n"
                        f"(nearest {times[i]:.3f})", xy=(0.5, 0.5),
                        xycoords="axes fraction", ha="center", va="center",
                        color=MUTED, fontsize=9)
            print(f"  {label} @ t*={T}: REFUSED, nearest {times[i]:.4f}", flush=True)
            continue

        fr = kio.load_frame(pfs[i])
        ad = fr.ds.all_data()
        z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
        p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
        w = np.asarray(ad[("amb_ions", "particle_weight")])
        g = np.sqrt(1.0 + (p / (mi * C)) ** 2)
        v = p / (mi * g) / vsh
        zp = np.percentile(np.asarray(ad[("piston_ions", "particle_position_x")]),
                           99.5) / sc.di0
        depth = float("nan")
        for lo in np.arange(0, 6, 0.25):
            m = (z - zp >= lo) & (z - zp < lo + 0.25)
            if m.sum() > 100 and v[m].std() < 0.06:
                depth = lo
                break

        H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
        ax.pcolormesh(zb, vb, np.ma.masked_where(H <= 0, H).T, norm=LogNorm(),
                      cmap="magma", shading="auto")
        ax.axhline(0, color="#fff", lw=.5, alpha=.25)
        ax.axhline(1, color="#5ad2ff", lw=.9, ls=(0, (5, 4)), alpha=.9)
        ax.axvline(zp, color="#ffd166", lw=1.2)
        if np.isfinite(depth):
            ax.axvspan(zp, zp + depth, color="#5ad2ff", alpha=0.07)
            ax.annotate(f"wedge {depth:.2f} $d_{{i0}}$", xy=(0.975, 0.06),
                        xycoords="axes fraction", ha="right", fontsize=9, color="#5ad2ff")
        if r == 0:
            ax.set_title(rf"$t\,\omega_{{ci0}}$ = {times[i]:.2f}", color=INK, fontsize=12)
        if c == 0:
            ax.set_ylabel(f"{label}\n" + rf"$\varepsilon$ = {eps:.4f}   "
                          rf"$\omega_{{pe}}/\omega_{{ce}}$ = {wpwc:.0f}" "\n"
                          rf"$\gamma$ = {gam:.4f}" "\n" r"$v_z/v_{sh}$",
                          color=INK2, fontsize=9)
        table.append((label, eps, wpwc, gam, T, times[i], depth))
    print(f"  row done: {label}", flush=True)

fig.suptitle(r"Ambient-ion phase space vs $\varepsilon = v_{te,ab}/c$ — the ONLY physics "
             "parameter that differs\n"
             "every row: $M_A$=13.95, $M_{ms}$=12.76, "
             r"$\beta_{ab}$=1150, $\rho_{i0}/d_{i0}$=10.4, $dz/\lambda_D$=0.379, "
             "$N_D$=264 (R1_paper: 0.60 / 167)",
             color=INK, fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.955])
os.makedirs("media/E_phase", exist_ok=True)
out = "media/E_phase/wedge_vs_eps.png"
fig.savefig(out, facecolor=SURFACE)
print(f"\n  wrote {out}")

print("\n" + "=" * 80)
print(f"{'run':20}{'eps':>9}{'wpe/wce':>9}{'gamma':>9}{'t* req':>9}{'t* got':>9}"
      f"{'wedge d_i0':>12}")
print("=" * 80)
for row in table:
    lab, eps, wpwc, gam, Treq, Tgot, depth = row
    dv = f"{depth:.2f}" if np.isfinite(depth) else "n/a"
    print(f"{lab:20}{eps:>9.4f}{wpwc:>9.1f}{gam:>9.4f}{Treq:>9.3f}{Tgot:>9.3f}{dv:>12}")
print("=" * 80)
