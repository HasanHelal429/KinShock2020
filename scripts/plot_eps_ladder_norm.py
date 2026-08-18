#!/usr/bin/env python3
"""The eps ladder replotted under three colour normalizations, to test whether the visual
difference in the wedge is a PLOTTING artifact rather than a physical one.

THE CONCERN (raised 2026-08-17). Each panel in plot_eps_ladder.py gets its OWN auto-scaled
LogNorm. The weighted 2-D histogram is a physical phase-space density and is itself
ppc-invariant -- the macroparticle weight scales as n*dz/ppc while the macroparticle count
per unit length scales inversely, so sum-of-weights is the same physical number either way.
But the FLOOR is not invariant: the faintest thing a panel can show is one macroparticle,
i.e. its own w_macro, and w_macro differs by ~10x across this ladder because dz does. A run
with a 10x lower floor renders 10x fainter structure, and on a per-panel log scale that
reads as "more filled in" even when the underlying density profile is identical.

So this script renders the same frames three ways:

  1. SHARED LOG    -- one vmin/vmax for every panel. Removes per-panel autoscaling.
  2. MATCHED FLOOR -- shared log, but vmin pinned to the COARSEST run's single-macroparticle
                      weight, so every panel is clipped at the noise floor of the worst-
                      sampled run. This is the honest like-for-like comparison: no panel is
                      allowed to show structure that the coarsest run could not have
                      resolved.
  3. LINEAR        -- shared linear scale. Log scales exaggerate faint tails by
                      construction; linear shows where the bulk of the phase-space density
                      actually is.

If the wedge difference is a rendering artifact, 2 and 3 collapse the rows together. If it
survives all three, the difference is in the phase-space density itself.

NOTE the scalar wedge metric is computed from v_z spread per z-window and never touches the
colour map, so it is unaffected by any of this. Whether it is affected by SAMPLING is a
separate question -- see scripts/eps_wedge_subsample.py.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
os.chdir("/pscratch/sd/h/hhelal/KinShock2020")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock.units import C, ME  # noqa: E402

SURFACE, INK, INK2, MUTED, BASE = "#0d0c10", "#f2f0f6", "#c9c5d4", "#8b8798", "#3a3745"
ZMAX, VMAX = 13.5, 2.6
TIMES = [0.130, 0.260]
TOL = 0.030

ROWS = [
    ("470 eV   (anchor)", "runs/S_phase/ss_dz16_ppc100", 0.0303, 100.0),
    ("1.5 keV",           "runs/E_phase/es_1p5keV",      0.0539,  56.2),
    ("4.7 keV",           "runs/E_phase/es_4p7keV",      0.0959,  31.6),
    ("15 keV",            "runs/E_phase/es_15keV",       0.1705,  17.8),
    ("47 keV",            "runs/E_phase/es_47keV",       0.3033,  10.0),
    ("R1_paper (47 keV)", "runs/R1_phase/R1_paper",      0.3033,  10.0),
]

zb = np.linspace(0, ZMAX, 340)
vb = np.linspace(-VMAX * 0.4, VMAX, 260)

# ---- load once, cache the histograms -------------------------------------------------
panels = {}     # (row, col) -> dict(H, zp, t)
wmacro = {}     # row -> median macroparticle weight
for r, (label, d, eps, wpwc) in enumerate(ROWS):
    cfg = kinshock.load(d)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(d)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])
    for c, T in enumerate(TIMES):
        i = int(np.argmin(np.abs(times - T)))
        if abs(times[i] - T) > TOL:
            print(f"  {label} @ t*={T}: REFUSED, nearest {times[i]:.4f}", flush=True)
            continue
        ad = kio.load_frame(pfs[i]).ds.all_data()
        z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
        p = np.asarray(ad[("amb_ions", "particle_momentum_z")])
        w = np.asarray(ad[("amb_ions", "particle_weight")])
        g = np.sqrt(1.0 + (p / (mi * C)) ** 2)
        v = p / (mi * g) / sc.vsh_model
        zp = np.percentile(np.asarray(ad[("piston_ions", "particle_position_x")]),
                           99.5) / sc.di0
        H, _, _ = np.histogram2d(z, v, bins=[zb, vb], weights=w)
        panels[(r, c)] = dict(H=H, zp=zp, t=times[i])
        wmacro[r] = float(np.median(w))
    print(f"  loaded: {label}   w_macro = {wmacro.get(r, float('nan')):.4e}", flush=True)

pos = np.concatenate([p["H"][p["H"] > 0].ravel() for p in panels.values()])
VMIN_SHARED = float(np.percentile(pos, 1.0))
VMAX_SHARED = float(max(p["H"].max() for p in panels.values()))
FLOOR = float(max(wmacro.values()))       # the COARSEST run's single-macroparticle weight

print(f"\n  macroparticle weights: " +
      "  ".join(f"{ROWS[r][0].split()[0]}={w:.3e}" for r, w in sorted(wmacro.items())))
print(f"  floor ratio (coarsest/finest) = {max(wmacro.values())/min(wmacro.values()):.2f}x")
print(f"  shared log vmin={VMIN_SHARED:.3e}  vmax={VMAX_SHARED:.3e}  "
      f"matched floor={FLOOR:.3e}")


def render(kind, norm_for, out, subtitle):
    fig, axes = plt.subplots(len(ROWS), len(TIMES), figsize=(11.5, 18.5), dpi=200,
                             facecolor=SURFACE, sharex=True, sharey=True)
    for r, (label, d, eps, wpwc) in enumerate(ROWS):
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
            pan = panels.get((r, c))
            if pan is None:
                continue
            H = pan["H"]
            ax.pcolormesh(zb, vb, np.ma.masked_where(H <= 0, H).T, norm=norm_for(),
                          cmap="magma", shading="auto")
            ax.axhline(0, color="#fff", lw=.5, alpha=.25)
            ax.axhline(1, color="#5ad2ff", lw=.9, ls=(0, (5, 4)), alpha=.9)
            ax.axvline(pan["zp"], color="#ffd166", lw=1.2)
            if r == 0:
                ax.set_title(rf"$t\,\omega_{{ci0}}$ = {pan['t']:.2f}", color=INK,
                             fontsize=12)
            if c == 0:
                ax.set_ylabel(f"{label}\n" + rf"$\varepsilon$ = {eps:.4f}   "
                              rf"$\omega_{{pe}}/\omega_{{ce}}$ = {wpwc:.0f}" "\n"
                              rf"$w_{{macro}}$ = {wmacro.get(r, float('nan')):.2e}" "\n"
                              r"$v_z/v_{sh}$", color=INK2, fontsize=9)
    fig.suptitle(f"Ambient-ion phase space vs $\\varepsilon$ — {kind}\n{subtitle}",
                 color=INK, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    os.makedirs("media/E_phase", exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"  wrote {out}")


render("SHARED LOG SCALE",
       lambda: LogNorm(vmin=VMIN_SHARED, vmax=VMAX_SHARED),
       "media/E_phase/wedge_vs_eps_sharedlog.png",
       "one vmin/vmax for every panel — removes per-panel autoscaling")

render("MATCHED NOISE FLOOR",
       lambda: LogNorm(vmin=FLOOR, vmax=VMAX_SHARED),
       "media/E_phase/wedge_vs_eps_matchedfloor.png",
       f"log, vmin pinned to the COARSEST run's single-macroparticle weight "
       f"({FLOOR:.2e}) — no panel may show what the worst-sampled run could not resolve")

render("SHARED LINEAR SCALE",
       lambda: Normalize(vmin=0.0, vmax=VMAX_SHARED),
       "media/E_phase/wedge_vs_eps_linear.png",
       "linear — shows where the bulk of the phase-space density actually is")
