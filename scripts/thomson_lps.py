"""Synthetic Thomson (EPW, IAW, ion streak) for LaserProdShock runs, at a fixed probe distance
from the target. These runs are in PHYSICAL units already, so no family mapping is applied, with
one exception: P7's heavy electrons (m_e x me_scale at real c) are put back on real-electron
velocities at the same temperature (v * sqrt(me_scale)).

Reuses scripts/thomson_physical.py (forward model, hybrid fluid electrons, parametric twin,
plotting). The Schaeffer PlasmaPy fork must be on PYTHONPATH.

Ion labels must be real PlasmaPy particles. The effective CH ion of the *_CH hybrid runs
(Z = 3, A = 6.502) has no exact counterpart, so --ions takes an explicit stand-in label (e.g.
"Li 3+", 6.94 u, 7% heavy); prefer the *_CpH twin (real C6+ / H+).

Example:
  PYTHONPATH=~/Schaeffer_PlasmaPy/src python scripts/thomson_lps.py \\
      /pscratch/sd/h/hhelal/LaserProdShock/runs/P8/P8_prod_dz0p5_CpH --hybrid \\
      --ions targ_ions_C:"C 6+":1.9921e-26 targ_ions_H:p+:1.6726e-27 \\
             amb_ions_C:"C 6+":1.9921e-26 amb_ions_H:p+:1.6726e-27 --position-mm 3
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thomson_physical as tp  # noqa: E402

ME, QE = tp.ME, tp.QE


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--hybrid", action="store_true",
                    help="fluid electrons from diag_fields (rho, Te); ions from --ion-prefix")
    ap.add_argument("--electrons", nargs="*", default=["targ_electrons", "amb_electrons"])
    ap.add_argument("--me-scale", type=float, default=1.0,
                    help="simulated / real electron mass (P7: 18.3615)")
    ap.add_argument("--ions", nargs="+", required=True, metavar="NAME:LABEL:MASS_KG")
    ap.add_argument("--particle-prefix", default=None,
                    help="particle plotfile prefix (default diag1 kinetic, diag_phase hybrid)")
    ap.add_argument("--position-mm", type=float, default=3.0)
    ap.add_argument("--probe-wavelength", type=float, default=532.0)
    ap.add_argument("--angle", type=float, default=90.0)
    ap.add_argument("--te-ref", type=float, default=470.0,
                    help="eV, sizes the IAW doublet part of the window")
    ap.add_argument("--iaw-halfwidths", type=float, default=2.5)
    ap.add_argument("--iaw-half-min", type=float, default=5.0)
    ap.add_argument("--epw-half", type=float, default=200.0)
    ap.add_argument("--wavelength-bins", type=int, default=600)
    ap.add_argument("--velocity-bins", type=int, default=2048)
    ap.add_argument("--smoothing-window", type=int, default=15)
    ap.add_argument("--v-ion-max-kms", type=float, default=4000.0)
    ap.add_argument("--v-e-max-c", type=float, default=0.3)
    ap.add_argument("--parametric", action="store_true")
    ap.add_argument("--vacuum-density", type=float, default=9e21,
                    help="m^-3; frames below it are vacuum (NaN). Default = 1%% of the OMEGA 2019 "
                         "ambient (n_e0 = 9e23). The pipeline's own default is 1%% of the run's "
                         "PEAK n_e at the probe, which blanks the upstream behind a dense front")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--note", default="")
    args = ap.parse_args()

    from plasmapy.diagnostics import pic_thomson as pt
    import astropy.units as u

    run = os.path.abspath(args.run_dir)
    rid = os.path.basename(run)
    diags = os.path.join(run, "diags")
    lps = run.split("/runs/")[0]
    phase = os.path.basename(os.path.dirname(run))
    cache = os.path.join(os.environ.get("THOMSON_CACHE", "/global/homes/h/hhelal/.claude/jobs/"
                                        "32b8a20b/tmp/thomson_lps_cache"), rid)
    os.makedirs(cache, exist_ok=True)          # never inside the run dir (owner's request)
    prefix = args.particle_prefix or ("diag_phase" if args.hybrid else "diag1")
    z_p = args.position_mm * 1e-3
    vi = args.v_ion_max_kms * 1e3
    ve = min(args.v_e_max_c / np.sqrt(args.me_scale), 0.999) * tp.C   # raw (heavy-e) units

    def read(name, mass, label, is_e, rng):
        nb = args.velocity_bins
        return pt.read_warpx_phase_space(
            diags, name, mass=mass, label=label, is_electron=is_e, prefix=prefix,
            n_velocity_bins=nb, n_position_bins=512, velocity_range=rng,
            cache=os.path.join(cache, f"{name}_{prefix}_v{nb}.npz"))

    ions = []
    for spec in args.ions:
        name, label, mass = spec.split(":")
        ions.append(read(name, float(mass), label, False, (-0.25 * vi, vi)))
    if args.hybrid:
        electrons = [tp.hybrid_electrons(pt, diags, ions, ions[0].x)]
    else:
        electrons = [read(n, ME * args.me_scale, None, True, (-ve, ve)) for n in args.electrons]
        # heavy electron -> real electron at the same T: v * sqrt(me_scale)
        electrons = [tp.scaled(e, 1.0 / np.sqrt(args.me_scale)) for e in electrons]

    ang = np.deg2rad(args.angle)
    lam0 = args.probe_wavelength
    cs = np.sqrt(args.te_ref * QE / tp.MP)
    iaw_half = max(args.iaw_halfwidths * tp.doppler_nm(cs, lam0, ang), args.iaw_half_min)
    win = lambda h: np.linspace(lam0 - min(h, lam0 - 40), lam0 + min(h, lam0 - 40),
                                args.wavelength_bins)
    cond = {"smoothing_window": args.smoothing_window, "smoothing_iterations": 2,
            "max_taper_bins": 20}

    def spectra(els, ios):
        return pt.spectra_from_phase_spaces(
            els, ios, position=z_p, reference_density=1 * u.m**-3,
            vacuum_density=args.vacuum_density * u.m**-3,
            probe_wavelength=lam0 * u.nm,
            epw_wavelengths=win(args.epw_half) * u.nm, iaw_wavelengths=win(iaw_half) * u.nm,
            epw_notches=[[lam0 - iaw_half, lam0 + iaw_half]] * u.nm, iaw_notches=None,
            probe_vec=(1.0, 0.0, 0.0), scatter_vec=(0.0, 1.0, 0.0),
            electron_conditioning=cond, ion_conditioning=cond,
            velocity_scale_factor=None, scattered_power=False, progress=True)

    vgrid = np.linspace(-0.25 * vi, vi, 600)

    def streak_of(ios):
        return sum(np.stack([np.interp(vgrid, i.v, row) for row in tp.at_position(i, z_p)])
                   for i in ios)

    t_ns = ions[0].t * 1e9
    rows = [("kinetic VDFs", spectra(electrons, ions), streak_of(ions))]
    if args.parametric:
        p_el = electrons if args.hybrid else [tp.maxwellianized(pt, e, ME) for e in electrons]
        p_io = [tp.maxwellianized(pt, i, float(s.split(":")[2])) for i, s in zip(ions, args.ions)]
        rows.append(("parametric: moment-matched Maxwellian per species",
                     spectra(p_el, p_io), streak_of(p_io)))
    args.xlabel = "t (ns)"
    args.streak_title = f"ion phase space at z = {args.position_mm:g} mm (log f)"
    args.position_di0 = 0
    args.streak_ylabel = "ion v_z (km/s)"
    args.title = (f"{rid} — Thomson at z = {args.position_mm:g} mm from the target "
                  f"({'hybrid: fluid e⁻' if args.hybrid else 'fully kinetic'}"
                  f"{'; e⁻ v×' + format(np.sqrt(args.me_scale), '.3g') if args.me_scale != 1 else ''}"
                  f"); {lam0:g} nm, {args.angle:g}°{'; ' + args.note if args.note else ''}")
    args.out_dir = os.path.join(lps, "media", phase, rid)
    args.tag = args.tag or f"thomson_{args.position_mm:g}mm".replace("thomson_", "")
    tp.plot(rid, run, rows, t_ns, vgrid, lam0, args, 1.0, 1.0)


if __name__ == "__main__":
    main()
