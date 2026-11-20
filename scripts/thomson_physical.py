"""Synthetic Thomson (EPW + IAW) of a run, mapped onto a PHYSICAL 470 eV hydrogen plasma, plus
the ion phase-space streak at the probe.

The runs sit at different points of the same dimensionless family, so each is put on
physical 470 eV hydrogen scales by rescaling velocities per species before the forward model:
  electrons: v / s_e,   s_e = sqrt(T_e,ab / 470 eV)        (10 for the 47 keV runs, 1 at 470 eV)
  ions:      v / s_i,   s_i = s_e * sqrt((m_p/m_e) / mu)    (light ion mu = 100 -> proton)
f is multiplied by the same factor so its zeroth moment, the density, is unchanged; n_e stays
the run's own (n_e,ab = 6e26 m^-3 for every run here). Ions are modelled as protons ("p+").
Time is shown as t*w_ci0, the invariant of the family.

Kinetic runs read every species from the particle dumps (prefix diag1). Hybrid runs
(H-PICShock: fluid electrons) read ions from `diag_phase`. The electron VDF is a Maxwellian
built from the fluid moments in `diag_fields`: n_e = rho/e, T_e = Te [eV], drift = ion bulk
v_z at the probe. That drift is exact in 1D, where J_z = 0. Their scales (d_i0, w_ci0,
T_e,ab, mu) come from the kinetic twin's config.

Usage (WarpX venv; fork on PYTHONPATH):
  PYTHONPATH=~/Schaeffer_PlasmaPy/src python scripts/thomson_physical.py runs/R1_phase/R1_470eV_s3f8_t2
  ... scripts/thomson_physical.py runs/hybrid_phase/hc_470eV_t2 --twin runs/R1_phase/R1_470eV_s3f8_t2
Output: media/<phase>/<id>/thomson_phys470.png (+ .npz)
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from dataclasses import replace

import numpy as np
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from kinshock import units  # noqa: E402
from kinshock.units import ME, MP  # noqa: E402

QE, C = 1.602176634e-19, 299792458.0
T_PHYS = 470.0  # eV: the physical regime every run is mapped onto


def doppler_nm(v, probe_nm, angle):
    return 2 * probe_nm * (v / C) * np.sin(angle / 2)


def scaled(ps, s):
    """Velocities divided by s, f multiplied by s: density unchanged."""
    return replace(ps, v=ps.v / s, f=ps.f * s)


def at_position(ps, z):
    """f(t, v) at the bin nearest z."""
    j = int(np.argmin(np.abs(ps.x - z)))
    return ps.f[:, :, j]


def maxwellianized(pt, ps, mass):
    """Drifting Maxwellian with the same n, u, T as `ps` at every (t, x), on its own grid."""
    import astropy.units as u

    v = ps.v
    n = np.trapezoid(ps.f, v, axis=1)
    safe = np.where(n > 0, n, 1.0)
    um = np.trapezoid(ps.f * v[None, :, None], v, axis=1) / safe
    var = np.trapezoid(ps.f * (v[None, :, None] - um[:, None, :]) ** 2, v, axis=1) / safe
    T = np.maximum(mass * var / QE, 1e-3)          # eV
    return pt.from_moments(np.maximum(n, 1.0), T * u.eV, np.where(n > 0, um, 0.0),
                           t=ps.t, x=ps.x, v=v, label=ps.label, is_electron=ps.is_electron,
                           mass=mass, meta={"source": "parametric (moment-matched Maxwellian)"})


def hybrid_electrons(pt, diags, ions, z_axis):
    """Maxwellian electrons from the hybrid fluid moments, on the ions' (t, x) grid."""
    import yt

    yt.set_log_level(40)
    t_ion = ions[0].t
    frames = sorted(glob.glob(os.path.join(diags, "diag_fields*")),
                    key=lambda p: int(os.path.basename(p)[len("diag_fields"):]))
    tf, ne, te = [], [], []
    for p in frames:
        ds = yt.load(p)
        cg = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
        nz = int(ds.domain_dimensions[0])
        z = np.linspace(float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0]), nz + 1)
        zc = 0.5 * (z[1:] + z[:-1])
        rho = np.asarray(cg["boxlib", "rho"]).reshape(nz)
        Te = np.asarray(cg["boxlib", "Te"]).reshape(nz)
        tf.append(float(ds.current_time))
        ne.append(np.interp(z_axis, zc, rho / QE))
        te.append(np.interp(z_axis, zc, Te))
    tf, ne, te = np.array(tf), np.array(ne), np.array(te)
    idx = [int(np.argmin(np.abs(tf - t))) for t in t_ion]   # field frame matching each ion dump
    ne, te = ne[idx], te[idx]
    # electron drift = charge-weighted ion bulk velocity, sum(Z n u) / sum(Z n): J_z = 0 in 1D
    from plasmapy.particles import Particle
    Zs = [abs(Particle(s.label).charge_number) for s in ions]
    num = sum(Z * np.trapezoid(s.f * s.v[None, :, None], s.v, axis=1) for Z, s in zip(Zs, ions))
    den = sum(Z * np.trapezoid(s.f, s.v, axis=1) for Z, s in zip(Zs, ions))
    drift = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    import astropy.units as u
    return pt.from_moments(np.maximum(ne, 1.0), te * u.eV, drift, t=t_ion, x=z_axis,
                           label="e-", is_electron=True, mass=ME,
                           meta={"source": "hybrid fluid moments"})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--twin", default=None,
                    help="kinetic run supplying scales (required for hybrid runs)")
    ap.add_argument("--position-di0", type=float, default=10.0)
    ap.add_argument("--probe-wavelength", type=float, default=532.0)
    ap.add_argument("--angle", type=float, default=90.0)
    ap.add_argument("--iaw-halfwidths", type=float, default=2.5)
    ap.add_argument("--wavelength-bins", type=int, default=600)
    ap.add_argument("--smoothing-window", type=int, default=15,
                    help="boxcar width in velocity bins (15 on the 2048-bin grids adds <~7%% "
                         "to the upstream thermal width)")
    ap.add_argument("--epw-half", type=float, default=200.0,
                    help="EPW window half-width in nm; the collective satellites downstream "
                         "(n_e ~ 1e26, alpha ~ 6) sit near +/-120 nm")
    ap.add_argument("--velocity-bins", type=int, default=2048)
    ap.add_argument("--iaw-half-min", type=float, default=5.0,
                    help="minimum IAW window half-width in nm (covers the bulk-flow shift)")
    ap.add_argument("--vacuum-fraction", type=float, default=1e-2,
                    help="frames with n_e below this fraction of the AMBIENT density count as "
                         "vacuum (NaN). The pipeline default is a fraction of the run's peak "
                         "n_e at the probe, which blanks the whole upstream once the "
                         "downstream exceeds 100 n_amb (H3_470eV_dense: 116 n_amb)")
    ap.add_argument("--tag", default="phys470")
    ap.add_argument("--parametric", action="store_true",
                    help="also model moment-matched Maxwellians of every species and write a "
                         "two-row figure (kinetic over parametric), tag <tag>_param")
    args = ap.parse_args()

    from plasmapy.diagnostics import pic_thomson as pt
    import astropy.units as u

    run = os.path.abspath(args.run_dir)
    rid = os.path.basename(run)
    hybrid = args.twin is not None
    cfg = yaml.safe_load(open(os.path.join(os.path.abspath(args.twin) if hybrid else run,
                                           "config.yaml")))
    sc = units.derive(cfg)
    mu = float(cfg["reference"]["mass_ratio"])
    s_e = np.sqrt(sc.Te_ab_eV / T_PHYS)
    s_i = s_e * np.sqrt((MP / ME) / mu)
    z_p = args.position_di0 * sc.di0
    print(f"{rid}: T_e,ab {sc.Te_ab_eV:.4g} eV, mu {mu:g} -> s_e {s_e:.4g}, s_i {s_i:.4g}; "
          f"probe at {args.position_di0:g} d_i0 = {z_p:.4g} m")

    diags = os.path.join(run, "diags")
    cache = os.path.join(run, "thomson_cache")
    os.makedirs(cache, exist_ok=True)
    prefix = "diag_phase" if hybrid else "diag1"

    # Bounded velocity ranges in RAW (simulation) units: a sparse hot tail otherwise sets
    # the grid (+-0.76 c for the 470 eV ambient electrons), leaving the 10 eV upstream with
    # ~2 bins per thermal width, which the smoothing then destroys.
    v_e = min(0.3 * s_e, 0.999) * C
    v_rng = {True: (-v_e, v_e), False: (-1.0 * sc.vsh_model, 3.0 * sc.vsh_model)}

    def read(name, mass, is_e):
        nb = args.velocity_bins
        return pt.read_warpx_phase_space(
            diags, name, mass=mass, label=None if is_e else "p+", is_electron=is_e,
            prefix=prefix, n_velocity_bins=nb, n_position_bins=256,
            velocity_range=v_rng[is_e],
            cache=os.path.join(cache, f"{name}_{prefix}_v{nb}_r{v_rng[is_e][1]:.4g}.npz"))

    ions = [read(n, mu * ME, False) for n in ("piston_ions", "amb_ions")]
    if hybrid:
        electrons = [hybrid_electrons(pt, diags, ions, ions[0].x)]
    else:
        electrons = [read(n, ME, True) for n in ("piston_electrons", "amb_electrons")]
    electrons = [scaled(e, s_e) for e in electrons]
    ions = [scaled(i, s_i) for i in ions]

    # Windows sized for the PHYSICAL 470 eV hydrogen plasma.
    ang = np.deg2rad(args.angle)
    cs = np.sqrt(T_PHYS * QE / MP)
    epw_half = args.epw_half
    # The downstream flow (~600 km/s physical) Doppler-shifts the whole ion feature by
    # ~1.5 nm, so the window must cover the drift, not just the +/- k C_s doublet.
    iaw_half = max(args.iaw_halfwidths * doppler_nm(cs, args.probe_wavelength, ang),
                   args.iaw_half_min)
    lam0 = args.probe_wavelength
    win = lambda h: np.linspace(lam0 - min(h, lam0 - 40), lam0 + min(h, lam0 - 40),
                                args.wavelength_bins)
    print(f"  EPW window +/-{epw_half:g} nm, IAW window +/-{iaw_half:.3g} nm (= EPW notch)")
    cond = {"smoothing_window": args.smoothing_window, "smoothing_iterations": 2,
            "max_taper_bins": 20}
    def spectra(els, ios):
        return pt.spectra_from_phase_spaces(
            els, ios, position=z_p, reference_density=1 * u.m**-3,
            vacuum_density=args.vacuum_fraction * sc.namb * u.m**-3,
            probe_wavelength=lam0 * u.nm,
            epw_wavelengths=win(epw_half) * u.nm, iaw_wavelengths=win(iaw_half) * u.nm,
            epw_notches=[[lam0 - iaw_half, lam0 + iaw_half]] * u.nm, iaw_notches=None,
            probe_vec=(1.0, 0.0, 0.0), scatter_vec=(0.0, 1.0, 0.0),
            electron_conditioning=cond, ion_conditioning=cond,
            velocity_scale_factor=None, scattered_power=False, progress=True)

    # Ion streak f_i(v_z, t) at the probe, all ion species on one velocity grid.
    vgrid = np.linspace(min(i.v.min() for i in ions), max(i.v.max() for i in ions), 600)

    def streak_of(ios):
        return sum(np.stack([np.interp(vgrid, i.v, row) for row in at_position(i, z_p)])
                   for i in ios)

    tw = ions[0].t * sc.wci0
    rows = [("kinetic VDFs", spectra(electrons, ions), streak_of(ions))]
    if args.parametric:
        # Hybrid electrons are already a Maxwellian; re-fitting them is the identity.
        p_el = electrons if hybrid else [maxwellianized(pt, e, ME) for e in electrons]
        p_io = [maxwellianized(pt, i, MP) for i in ions]
        rows.append(("parametric: moment-matched Maxwellian per species",
                     spectra(p_el, p_io), streak_of(p_io)))
    plot(rid, run, rows, tw, vgrid, lam0, args, s_e, s_i)


def _get(spec, *names):
    for n in names:
        if hasattr(spec, n):
            return getattr(spec, n)
    raise AttributeError(f"spectrogram has none of {names}: {[a for a in dir(spec) if not a.startswith('_')]}")


def plot(rid, run, rows, tw, vgrid, lam0, args, s_e, s_i):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(len(rows), 3, figsize=(16, 4.8 * len(rows)), constrained_layout=True,
                            squeeze=False)
    saved = {}
    for r, (rlab, spec, streak) in enumerate(rows):
        epw, iaw, lw_e, lw_i, tws = _row(fig, axs[r], spec, tw, vgrid, streak, lam0, args,
                                         rlab if len(rows) > 1 else "")
        key = "" if r == 0 else "_param"
        saved.update({f"epw{key}": epw, f"iaw{key}": iaw, f"ion_f{key}": streak})
    title = getattr(args, "title", None) or (
        f"{rid} — Thomson mapped to physical 470 eV hydrogen "
        f"(electrons v/{s_e:.3g}, ions v/{s_i:.3g}); {lam0:g} nm probe, {args.angle:g}°")
    fig.suptitle(title, fontsize=11)
    phase = os.path.basename(os.path.dirname(run))
    out = getattr(args, "out_dir", None) or os.path.join(ROOT, "media", phase, rid)
    out = out.replace(os.path.join(".claude", "worktrees", "psc-schaeffer2020") + os.sep, "")
    os.makedirs(out, exist_ok=True)
    tag = args.tag + ("_param" if len(rows) > 1 else "")
    png = os.path.join(out, f"thomson_{tag}.png")
    fig.savefig(png, dpi=130)
    np.savez_compressed(png[:-4] + ".npz", t_wci0=tws, epw_nm=lw_e, iaw_nm=lw_i, ion_v=vgrid,
                        ion_t_wci0=tw, s_e=s_e, s_i=s_i, **saved)
    print("wrote", png)


def _row(fig, ax, spec, tw, vgrid, streak, lam0, args, rlab):
    epw = np.asarray(_get(spec, "epw", "epw_spectrogram", "S_epw"))
    iaw = np.asarray(_get(spec, "iaw", "iaw_spectrogram", "S_iaw"))
    lw_e = np.asarray(_get(spec, "epw_wavelengths", "wavelengths_epw")) * 1e9
    lw_i = np.asarray(_get(spec, "iaw_wavelengths", "wavelengths_iaw")) * 1e9
    ts = np.asarray(_get(spec, "t", "times"))
    tws = tw if len(ts) == len(tw) else np.interp(ts, np.linspace(ts.min(), ts.max(), len(tw)), tw)
    pre = f"{rlab}\n" if rlab else ""
    for a, S, lw, ttl in ((ax[0], epw, lw_e, pre + "EPW (electron feature)"),
                          (ax[1], iaw, lw_i, pre + "IAW (ion feature)")):
        S2 = S if S.shape == (len(lw), len(tws)) else S.T
        # Each time column normalised to its own 99.5th percentile: the collective EPW
        # resonance throws single-pixel spikes (extrapolated tail), and the downstream is
        # ~50x denser than the upstream, so a global max hides everything else.
        ref = np.nanpercentile(np.where(S2 > 0, S2, np.nan), 99.5, axis=0)
        S2 = np.where(np.isfinite(ref) & (ref > 0), S2 / np.where(ref > 0, ref, 1), np.nan)
        S2 = np.log10(np.clip(S2, 1e-3, 1.0))
        m = a.pcolormesh(tws, lw, S2, shading="nearest", cmap="inferno", vmin=-3, vmax=0)
        fig.colorbar(m, ax=a, label="log10 S (per-time normalised)")
        a.axhline(lam0, color="c", lw=0.6, ls="--")
        a.set_xlabel(getattr(args, "xlabel", "t·ω_ci0"))
        a.set_ylabel("wavelength (nm)")
        a.set_title(ttl)
    # Same 4-decade floor for kinetic and parametric: a Maxwellian is nonzero everywhere, so
    # without a floor its 1e-300 tails colour the whole panel.
    fmax = np.nanmax(streak) if np.any(streak > 0) else 1.0
    flog = np.log10(np.clip(streak / fmax, 1e-4, None))
    flog = np.where(flog <= -4, np.nan, flog)
    ax[2].pcolormesh(tw, vgrid / 1e3, flog.T, shading="nearest", cmap="viridis", vmin=-4, vmax=0)
    ax[2].set_xlabel(getattr(args, "xlabel", "t·ω_ci0"))
    ax[2].set_ylabel(getattr(args, "streak_ylabel", "ion v_z (km/s, physical H⁺)"))
    ax[2].set_title(pre + (getattr(args, "streak_title", None)
                           or f"ion phase space at z = {args.position_di0:g} d_i0 (log f)"))
    return epw, iaw, lw_e, lw_i, tws


if __name__ == "__main__":
    main()
