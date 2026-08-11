#!/usr/bin/env python
"""Grid-heating rate of a uniform ambient box, for one or more H_phase runs.

    python scripts/heating_rate.py runs/H_phase/*
    python scripts/heating_rate.py runs/H_phase/* --ref hs_dz1_ppc100 --window 5 30

WHY THIS IS A DIFFERENT SCRIPT FROM grid_heating.py. That one exists because in the
full piston-driven deck the ParticleEnergy reduced diagnostic is useless: it reports one
mean energy per species over the WHOLE domain, and that number is dominated by real
physical heating of the plasma the piston has already swept up (2.2x within 3.4 t_ab in
R1_paper). It therefore has to bin particles by z out of the plotfiles and keep only the
far upstream. In an H_phase box there IS no piston, no injector and no shock, so nothing
physical can change the temperature and EP is exactly the right instrument -- domain-wide,
unwindowed, ~1000 rows per run, and free. Do not point this script at a production run.

WHAT IT REPORTS, AND WHAT TO QUOTE. The headline is dT over the run's own window, which
is a measurement and needs no extrapolation -- 30 t_ab by construction, the same window
h0_baseline measured (+29.7 eV) in the full deck. The projection to the 190 t_ab
production window is printed too, but it is a LINEAR extrapolation 6.3x past its data and
is there only to say whether a point is obviously survivable; quote rate RATIOS against
the reference point instead. That is the trap the h0_baseline fitted asymptote fell into
(RESULTS 2026-08-05: 132 -> 87 eV, extrapolated 2-3x past its data, not quotable).

    T = (2/3) <KE> / k_B      three degrees of freedom, non-relativistic here
                              (v_te/c = 4.4e-3 at 10 eV, so gamma - 1 -> u^2/2 exactly)

The check that makes the number trustworthy is `N drift`: with periodic boundaries and no
injector the macroparticle count must be constant to the digit. Anything else means
particles are being created or lost and the mean energy is not a temperature.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import kinshock                                       # noqa: E402
from kinshock import io as kio                        # noqa: E402
from kinshock.units import C, EPS0, ME, ME_C2_EV, QE  # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BASELINE = "#c3c2b7"
# categorical, colour-blind safe; the first two are the validated diverging pair's hues
SERIES = ["#2a78d6", "#eb6834", "#3fa34d", "#8e5ea2", "#c0392b", "#0f8b8d",
          "#d4a017", "#6b6b6b"]


def resolution(cfg, sc):
    """(dz/lambda_D, N_D, lambda_D [m]) for the AMBIENT plasma of this run.

    lambda_D is built from the ambient temperature and density, not the ablation ones:
    it is the upstream that grid-heats, and lambda_D ~ sqrt(T) makes the two differ by
    10x in this parameter set (RESULTS 2026-08-04).
    """
    T0_eV = float(cfg["plasma"]["ambient"]["theta_0"]) * ME_C2_EV
    lam = math.sqrt(EPS0 * T0_eV * QE / (sc.namb * QE * QE))
    ppc = int(cfg["numerics"]["ppc"]["ambient"])
    return sc.dz / lam, ppc * lam / sc.dz, lam


def temperatures(run_dir, species=None):
    """(t [s], {species: T [eV]}, N_drift) from the EP/PN reduced diagnostics.

    ``<species>_mean(J)`` is WarpX's mean kinetic energy per physical particle, so the
    conversion is exact and needs no particle loop.
    """
    hdr, dat = kio.reduced_diag(run_dir, "EP")
    if dat.ndim == 1:
        dat = dat[None, :]
    cols = {}
    for i, h in enumerate(hdr):
        name = h.split("]", 1)[-1]
        if name.endswith("_mean(J)"):
            cols[name[: -len("_mean(J)")]] = i
    cols.pop("total", None)
    t = dat[:, 1]
    T = {s: (2.0 / 3.0) * dat[:, i] / QE for s, i in cols.items()
         if species is None or s in species}
    try:
        _, pn = kio.reduced_diag(run_dir, "PN")
        n = pn[:, 2] if pn.ndim > 1 else np.array([pn[2]])
        drift = float(n.max() / n.min() - 1.0) if n.min() > 0 else math.nan
    except FileNotFoundError:
        drift = math.nan
    return t, T, drift


def fit_rate(t_tab, T, lo, hi):
    """Least-squares dT/dt [eV per t_ab] over the window, plus the endpoint delta.

    The delta is the measurement; the slope is what projects. They disagree when the
    heating is not linear, which is itself worth seeing -- so both are returned.
    """
    m = (t_tab >= lo) & (t_tab <= hi)
    if m.sum() < 3:
        return math.nan, math.nan, math.nan
    slope = float(np.polyfit(t_tab[m], T[m], 1)[0])
    return slope, float(T[m][0]), float(T[m][-1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--species", default="amb_ions,amb_electrons",
                    help="comma-separated; default the two ambient species")
    ap.add_argument("--window", type=float, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="fit window in t_ab (default: 10%%..100%% of what ran)")
    ap.add_argument("--project-to", type=float, default=190.0,
                    help="t_ab of the production window for the projection (default 190)")
    ap.add_argument("--ref", default=None,
                    help="run_id whose rate the ratio column is against "
                         "(default: the first run given)")
    ap.add_argument("--out", default=None, help="figure path (default media/H_phase/heating_rate.png)")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    want = [s.strip() for s in args.species.split(",") if s.strip()]
    rows = []
    for d in args.run_dirs:
        if not os.path.isfile(os.path.join(d, "config.yaml")):
            continue
        try:
            cfg = kinshock.load(d)
            sc = kinshock.units.derive(cfg)
            t, T, drift = temperatures(d, want)
        except Exception as e:                # a missing run must not look like a flat line
            print(f"  SKIP {d}: {type(e).__name__}: {e}")
            continue
        if not T:
            print(f"  SKIP {d}: EP has none of {want}")
            continue
        dzlam, ND, lam = resolution(cfg, sc)
        t_tab = t / sc.t_ab
        lo, hi = args.window or (0.10 * t_tab.max(), t_tab.max())
        rows.append(dict(rid=cfg["meta"]["run_id"], dzlam=dzlam, ND=ND, t_tab=t_tab,
                         T=T, drift=drift, lo=lo, hi=hi, ran=t_tab.max(),
                         fits={s: fit_rate(t_tab, T[s], lo, hi) for s in T}))

    if not rows:
        print("no runs with EP output yet")
        return 1

    ref = next((r for r in rows if r["rid"] == args.ref), rows[0])
    key = "amb_ions" if "amb_ions" in ref["T"] else sorted(ref["T"])[0]

    print(f"\ngrid heating, ambient box — fit window {rows[0]['lo']:.1f}..{rows[0]['hi']:.1f} t_ab"
          f", projection to {args.project_to:g} t_ab, ratios vs {ref['rid']}\n")
    hdr = (f"{'run':>18} {'dz/lD':>6} {'N_D':>6} {'ran':>7} {'species':>15} "
           f"{'T0':>7} {'T_end':>8} {'dT':>8} {'eV/t_ab':>9} {'proj':>9} {'ratio':>6} {'Ndrift':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        for s in sorted(r["T"]):
            slope, T0, Tend = r["fits"][s]
            rslope = ref["fits"].get(s, (math.nan,))[0]
            ratio = slope / rslope if rslope not in (0.0, math.nan) else math.nan
            proj = T0 + slope * args.project_to
            print(f"{r['rid']:>18} {r['dzlam']:6.2f} {r['ND']:6.1f} {r['ran']:7.1f} {s:>15} "
                  f"{T0:7.2f} {Tend:8.2f} {Tend-T0:8.2f} {slope:9.4f} {proj:9.1f} "
                  f"{ratio:6.2f} {r['drift']:8.1e}")
    print("\n  dT is the MEASUREMENT over the window; proj is a linear extrapolation "
          f"{args.project_to/max(1e-9, rows[0]['hi']):.1f}x past it — quote the ratio, not proj.")
    bad = [r["rid"] for r in rows if not (r["drift"] < 1e-9 or math.isnan(r["drift"]))]
    if bad:
        print(f"  ! macroparticle count is NOT constant in {bad} — the mean energy is "
              f"not a temperature there; check the boundaries.")

    if args.no_plot:
        return 0

    import matplotlib                                  # noqa: E402
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt                    # noqa: E402
    from kinshock import plotting as P                 # noqa: E402

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), dpi=200, facecolor=SURFACE,
                             gridspec_kw=dict(width_ratios=[1.35, 1.0]))
    for ax in axes:
        ax.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(BASELINE)
        ax.tick_params(colors=MUTED, labelsize=9, length=3)
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_color(INK2)

    for i, r in enumerate(rows):
        c = SERIES[i % len(SERIES)]
        axes[0].plot(r["t_tab"], r["T"][key], "-", color=c, lw=1.6,
                     label=f"{r['rid']}  (dz/$\\lambda_D$ {r['dzlam']:.2f}, $N_D$ {r['ND']:.0f})")
    axes[0].axhline(10.0, color=BASELINE, lw=1.0, ls=(0, (4, 4)))
    axes[0].text(axes[0].get_xlim()[1], 10.0, " $T_0$ ", va="center", fontsize=8, color=MUTED)
    axes[0].set_xlabel(r"$t / t_{ab}$", color=INK2, fontsize=10)
    axes[0].set_ylabel(f"{key}  T [eV]", color=INK2, fontsize=10)
    axes[0].set_title("uniform ambient box: nothing but the grid can heat this",
                      color=INK, fontsize=11.5, loc="left", pad=10)
    leg = axes[0].legend(frameon=False, fontsize=8)
    for x in leg.get_texts():
        x.set_color(INK2)

    # the decision panel: dT against aliasing, one line per noise level
    byND = {}
    for r in rows:
        byND.setdefault(round(r["ND"]), []).append(r)
    for j, (nd, group) in enumerate(sorted(byND.items())):
        group.sort(key=lambda r: r["dzlam"])
        x = [g["dzlam"] for g in group]
        y = [g["fits"][key][2] - g["fits"][key][1] for g in group]
        axes[1].plot(x, y, "-o", color=SERIES[j % len(SERIES)], lw=1.8, ms=5,
                     label=f"$N_D$ = {nd}")
    axes[1].set_xscale("log")
    axes[1].axvline(6.07, color=BASELINE, lw=1.0, ls=(0, (4, 4)))
    axes[1].text(6.07, axes[1].get_ylim()[1], " production ", ha="right", va="top",
                 fontsize=8, color=MUTED, rotation=90)
    axes[1].set_xlabel(r"$dz / \lambda_{D,amb}$", color=INK2, fontsize=10)
    axes[1].set_ylabel(r"$\Delta T$ over the window [eV]", color=INK2, fontsize=10)
    axes[1].set_title("aliasing (x) vs noise (colour)", color=INK, fontsize=11.5,
                      loc="left", pad=10)
    leg = axes[1].legend(frameon=False, fontsize=9)
    for x in leg.get_texts():
        x.set_color(INK2)

    fig.text(0.01, 0.015,
             "Falling along a coloured line = aliasing, fixed by dz. Separation BETWEEN "
             "lines at fixed x = particle noise, fixed 4x more cheaply by ppc.",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.90, bottom=0.135, wspace=0.22)
    out = args.out or os.path.join(P.ROOT, "media", "H_phase", "heating_rate.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
