#!/usr/bin/env python
"""Measure numerical (grid) heating of the FAR UPSTREAM, spatially resolved.

WHY NOT JUST USE EP.txt. The reduced ParticleEnergy diagnostic reports one mean
energy per species over the WHOLE domain, and that number is dominated by real
physical heating of the ambient plasma the piston has already swept up -- in
runs/R1_paper it rises 2.2x within the first 3.4 t_ab, long before the shock forms.
It is therefore useless as a grid-heating diagnostic. Numerical heating has to be
read where nothing physical has happened yet: far ahead of the piston.

WHAT THIS DOES. For each plotfile, bins the chosen species by z, keeps only the far
upstream (default: outer 25% of the domain), and computes the weight-averaged kinetic
energy per particle there, exactly as EP does but restricted spatially:

    KE = (gamma - 1) m c^2,   gamma = sqrt(1 + u^2),   u = p/(m c)
    T   = (2/3) <KE> / k      (3 degrees of freedom)

It then fits dT/dt in eV per t_ab, projects to the full 220 t_ab, and -- if a second
run is given -- reports the RATIO of heating rates, which is the robust number: it
cancels most of the extrapolation risk and both runs' absolute normalizations.

It also reports how far the piston has actually reached (from the piston-ion density),
so you can confirm the "far upstream" window really is pristine. If `piston reach`
ever approaches `window lo`, the measurement is contaminated and the window must move.

    python scripts/grid_heating.py runs/R1_paper_470eV_pilot runs/R1_paper
    python scripts/grid_heating.py runs/R1_paper --species amb_ions --upstream-frac 0.85
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import kinshock                                    # noqa: E402
from kinshock import io as kio                      # noqa: E402
from kinshock.units import C, ME, QE                # noqa: E402


def species_mean_energy(frame, species, mass, z_lo, z_hi):
    """Weight-averaged kinetic energy [J] per particle of ``species`` in [z_lo, z_hi].

    Uses all three momentum components so the result is directly comparable to the
    ParticleEnergy reduced diagnostic. Returns (mean_KE, total_weight, n_macro).
    """
    ad = frame.ds.all_data()
    try:
        z = np.asarray(ad[(species, "particle_position_x")])
        w = np.asarray(ad[(species, "particle_weight")])
        u2 = np.zeros_like(z)
        for comp in ("x", "y", "z"):
            p = np.asarray(ad[(species, f"particle_momentum_{comp}")]) / (mass * C)
            u2 += p * p
    except Exception:
        return math.nan, 0.0, 0
    m = (z >= z_lo) & (z <= z_hi)
    if not m.any() or w[m].sum() <= 0:
        return math.nan, 0.0, 0
    ke = (np.sqrt(1.0 + u2[m]) - 1.0) * mass * C * C
    return float(np.average(ke, weights=w[m])), float(w[m].sum()), int(m.sum())


def piston_reach(frame, scales, frac=1e-3):
    """Rightmost z [m] where the piston-ion density exceeds ``frac`` of the ambient.

    Cheap contamination check for the upstream window.
    """
    try:
        n = kio.species_density(frame, "piston_ions")
    except Exception:
        return math.nan
    zc = np.linspace(0.0, scales.domain_halfwidth, len(n))
    hit = np.nonzero(n > frac * scales.namb)[0]
    return float(zc[hit[-1]]) if len(hit) else 0.0


def analyse(run_dir, species, frac, tab_max=None):
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    mass = sc.mi if "ion" in species else ME
    T0_eV = float(cfg["plasma"]["ambient"]["theta_0"]) * 511000.0
    z_lo, z_hi = frac * sc.domain_halfwidth, sc.domain_halfwidth
    rows = []
    for p in kio.plotfiles(run_dir):
        fr = kio.load_frame(p)
        tab = fr.time / sc.t_ab
        if tab_max is not None and tab > tab_max:
            continue
        ke, _, nmac = species_mean_energy(fr, species, mass, z_lo, z_hi)
        if math.isnan(ke):
            continue
        rows.append((tab, (2.0 / 3.0) * ke / QE, piston_reach(fr, sc) / sc.domain_halfwidth,
                     nmac))
    return cfg, sc, T0_eV, np.array([r[:3] for r in rows]), [r[3] for r in rows], (z_lo, z_hi)


def report(run_dir, species, frac, tab_max, tau_sim_tab=220.0):
    cfg, sc, T0, a, nmac, (z_lo, z_hi) = analyse(run_dir, species, frac, tab_max)
    rid = cfg["meta"]["run_id"]
    print(f"\n=== {rid} ===  species={species}  t_ab={sc.t_ab*1e12:.4f} ps  "
          f"T_0(init)={T0:.3g} eV")
    print(f"    upstream window: z in [{z_lo*1e3:.3f}, {z_hi*1e3:.3f}] mm "
          f"(outer {100*(1-frac):.0f}% of {sc.domain_halfwidth*1e3:.3f} mm), "
          f"dz/lambda_D,ab = {sc.dz/(math.sqrt(float(cfg['plasma']['piston']['theta_e_heat']))*sc.de):.3g}")
    if len(a) < 2:
        print("    not enough frames yet")
        return None
    print(f"    {'t/t_ab':>9}{'t*wci0':>9}{'T [eV]':>12}{'T/T_0':>9}"
          f"{'piston reach':>14}{'n_macro':>10}")
    for (tab, T, reach), nm in zip(a, nmac):
        print(f"    {tab:9.3f}{tab*sc.t_ab/sc.wci0_inv:9.4f}{T:12.4g}{T/T0:9.3f}"
              f"{reach:13.3%}{nm:10d}")
    slope = np.polyfit(a[:, 0], a[:, 1], 1)[0]          # eV per t_ab
    print(f"    dT/dt = {slope:.4g} eV/t_ab = {slope/T0*100:.4g} %T_0/t_ab"
          f"  ({slope*sc.wci0_inv/sc.t_ab:.4g} eV per 1/wci0)")
    proj = a[0, 1] + slope * tau_sim_tab
    print(f"    linear projection to {tau_sim_tab:g} t_ab: {proj:.4g} eV "
          f"= {proj/T0:.3g} x T_0")
    return dict(run=rid, T0=T0, slope=slope, frac_slope=slope / T0, proj=proj,
                proj_over_T0=proj / T0, sc=sc, arr=a)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run dirs; the first is the run under test")
    ap.add_argument("--species", default="amb_electrons")
    ap.add_argument("--upstream-frac", type=float, default=0.75,
                    help="keep z > frac*domain as 'far upstream' (default 0.75)")
    ap.add_argument("--tab-max", type=float, default=None,
                    help="only use frames up to this t/t_ab (match windows across runs)")
    ap.add_argument("--tau-sim-tab", type=float, default=220.0)
    args = ap.parse_args(argv)

    out = [report(r, args.species, args.upstream_frac, args.tab_max, args.tau_sim_tab)
           for r in args.runs]
    out = [o for o in out if o]
    if len(out) >= 2:
        a, b = out[0], out[1]
        print(f"\n=== {a['run']} vs {b['run']} ===")
        print(f"  heating rate  [eV/t_ab]      {a['slope']:.4g}   vs {b['slope']:.4g}")
        print(f"  as %T_0 per t_ab             {a['frac_slope']*100:.4g}   vs "
              f"{b['frac_slope']*100:.4g}    RATIO = {a['frac_slope']/b['frac_slope']:.3g}x")
        print(f"  projected T/T_0 at 220 t_ab  {a['proj_over_T0']:.3g}   vs "
              f"{b['proj_over_T0']:.3g}")
        print("\n  The RATIO of fractional rates is the robust number: it cancels both\n"
              "  runs' absolute T_0 and most of the linear-extrapolation risk.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
