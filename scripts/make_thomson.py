"""Synthetic Thomson scattering spectra (EPW + IAW) from a run's particle plotfiles.

Forward-models the Thomson-scattered spectrum a real diagnostic would record, at one
point in the domain, for every particle frame of a run. Writes two spectrograms and
the underlying arrays into ``media/<run_id>/``:

    thomson_epw.png      electron feature (Doppler-broadened by v_te)
    thomson_iaw.png      ion-acoustic feature (the +/- k*C_s doublet)
    thomson_spectra.npz  t, both spectrograms, both wavelength axes, alpha, n_e

Config-driven like every other script here: species come from the ``species:`` block,
the ion mass from ``reference.mass_ratio``, and C_s,ab (which sizes the IAW window)
from ``kinshock.units.derive``. Nothing about the run is hard-coded.

WHETHER THERE IS AN ION FEATURE AT ALL DEPENDS ON THE DENSITY. The collectivity
parameter alpha = 1/(k lambda_D) scales as sqrt(n)/T, so the same deck gives
alpha ~ 1e-5 at n0 = 1e18 m^-3 and ~0.5 at Table I's 6e26 (RESULTS 2026-08-03). Below
alpha ~ 1 the electron susceptibility vanishes, the ion feature disappears, and the
spectrum reduces to the Doppler-mapped electron VDF. The script reports alpha per
frame rather than assuming a regime -- read it before interpreting the IAW panel.

TWO ENVIRONMENTS. The forward model needs PyTorch and the plotfile reader needs yt,
and no single conda env here has both:

    /opt/anaconda3/envs/physics   yt, no torch     -> can read, cannot model
    /opt/anaconda3/envs/tsnn      torch, no yt     -> can model, cannot read

So the binned phase spaces are cached to ``runs/<ID>/thomson_cache/`` and the two
halves run separately. ``--stage`` defaults to ``auto``, which does as much as the
current interpreter allows and prints the exact command for the other half:

    python scripts/make_thomson.py runs/R1_paper                    # in physics: reads
    /opt/anaconda3/envs/tsnn/bin/python scripts/make_thomson.py runs/R1_paper   # models

The cache key includes the particle mass, so it is deliberately taken from
``kinshock.units`` rather than from astropy: astropy 8.0.0 (physics) and 7.0.1 (tsnn)
ship different CODATA m_e, and a cache written by one env would be rejected by the
other. Delete ``thomson_cache/`` to force a re-read.

The PlasmaPy fork carrying the pipeline is found at ``$KINSHOCK_PLASMAPY`` or
``~/Schaeffer_PlasmaPy/src`` (branch ``feature/pic-thomson-pipeline``).
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from kinshock import config as kconfig  # noqa: E402
from kinshock import units  # noqa: E402
from kinshock.units import ME  # noqa: E402

PLASMAPY_SRC = os.environ.get(
    "KINSHOCK_PLASMAPY", os.path.expanduser("~/Schaeffer_PlasmaPy/src")
)

# The pipeline calls np.trapezoid (numpy >= 2). tsnn still has numpy 1.26, where the
# same function is np.trapz -- identical signature and semantics, renamed in 2.0.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz


def _pt():
    """Import the fork's pic_thomson, with a pointed message if it is not there."""
    if PLASMAPY_SRC not in sys.path:
        sys.path.insert(0, PLASMAPY_SRC)
    try:
        from plasmapy.diagnostics import pic_thomson
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise SystemExit(
            f"cannot import plasmapy.diagnostics.pic_thomson from {PLASMAPY_SRC}.\n"
            "Point KINSHOCK_PLASMAPY at a checkout of the Schaeffer PlasmaPy fork "
            "(branch feature/pic-thomson-pipeline)."
        ) from exc
    return pic_thomson


def _have(module: str) -> bool:
    import importlib.util

    return importlib.util.find_spec(module) is not None


def species_of(cfg, kind):
    """Species names of one ``kind`` ('electron' / 'ion'), in config order."""
    return [n for n, s in cfg["species"].items() if s.get("kind") == kind]


def read_phase_spaces(pt, run_dir, cfg, args, cache_dir):
    """Bin every species' phase space, memoising to *cache_dir*."""
    os.makedirs(cache_dir, exist_ok=True)
    mass_ratio = float(cfg["reference"]["mass_ratio"])
    diags = os.path.join(run_dir, "diags")

    def read(name, mass, *, is_electron):
        return pt.read_warpx_phase_space(
            diags,
            name,
            mass=mass,
            label=None if is_electron else args.ion_label,
            is_electron=is_electron,
            cache=os.path.join(cache_dir, f"{name}.npz"),
            n_velocity_bins=args.velocity_bins,
            n_position_bins=args.position_bins,
            timesteps=None,
        )

    electrons = [read(n, ME, is_electron=True) for n in species_of(cfg, "electron")]
    ions = [read(n, mass_ratio * ME, is_electron=False) for n in species_of(cfg, "ion")]
    return electrons, ions


def peak_sigma(pt, phase_spaces):
    """Widest velocity sigma across species, each at its own densest cell."""
    widest = 0.0
    for ps in phase_spaces:
        density = pt.number_density(ps.f, ps.v)
        step, cell = np.unravel_index(np.argmax(density), density.shape)
        f = ps.f[step, :, cell]
        total = np.trapezoid(f, ps.v)
        if not np.isfinite(total) or total <= 0:
            continue
        mean = np.trapezoid(f * ps.v, ps.v) / total
        var = np.trapezoid(f * (ps.v - mean) ** 2, ps.v) / total
        widest = max(widest, float(np.sqrt(max(var, 0.0))))
    return widest


def doppler_nm(speed, probe_nm, angle_rad, c):
    """Wavelength shift at *probe_nm* for a scatterer moving at *speed*."""
    return probe_nm * 2.0 * np.sin(angle_rad / 2.0) * speed / c


def window(half_nm, probe_nm, bins):
    half = min(half_nm, probe_nm - 40.0)
    return np.linspace(probe_nm - half, probe_nm + half, bins)


def model(pt, electrons, ions, cfg, scales, args):
    """Forward-model both features. Returns the pipeline's spectrogram object."""
    import astropy.units as u

    c = 299792458.0
    angle = np.deg2rad(args.angle)

    e_sigma = peak_sigma(pt, electrons)
    i_sigma = peak_sigma(pt, ions)
    print(f"  electron sigma {e_sigma:.4g} m/s ({e_sigma / c:.3f} c)")
    print(f"  ion      sigma {i_sigma:.4g} m/s ({i_sigma / c:.4f} c)")

    epw_half = 2.0 * doppler_nm(e_sigma, args.probe_wavelength, angle, c)
    # Size the IAW window from the ion-acoustic speed, NOT the ion thermal sigma:
    # the piston drift dominates that sigma, and using it widens the window ~6x
    # until the doublet spans 2-3 pixels.
    cs_shift = doppler_nm(scales.Cs_ab, args.probe_wavelength, angle, c)
    iaw_half = args.iaw_halfwidths * cs_shift
    print(f"  C_s,ab {scales.Cs_ab:.4g} m/s -> IAW doublet at +/-{cs_shift:.1f} nm")
    print(f"  EPW window +/-{min(epw_half, args.probe_wavelength - 40):.1f} nm")
    print(f"  IAW window +/-{min(iaw_half, args.probe_wavelength - 40):.1f} nm")

    position = args.position
    if position is None:
        position = float(np.mean(electrons[0].x))
        print(f"  no --position; sampling the domain centre, {position:.5g} m")

    conditioning = {
        "smoothing_window": args.smoothing_window,
        "smoothing_iterations": args.smoothing_iterations,
        "max_taper_bins": 20,
    }

    return pt.spectra_from_phase_spaces(
        electrons,
        ions,
        position=position,
        # The reader already scales f so its zeroth moment is a density in m^-3.
        # This is a unit scale, not the run's density -- passing the real density
        # double-counts it and drives n_e to ~1e51 and alpha to NaN.
        reference_density=1 * u.m**-3,
        probe_wavelength=args.probe_wavelength * u.nm,
        epw_wavelengths=window(epw_half, args.probe_wavelength, args.wavelength_bins) * u.nm,
        iaw_wavelengths=window(iaw_half, args.probe_wavelength, args.wavelength_bins) * u.nm,
        probe_vec=(1.0, 0.0, 0.0),
        scatter_vec=(0.0, 1.0, 0.0),
        electron_conditioning=conditioning,
        ion_conditioning=conditioning,
        # S(k, omega) rather than power per unit wavelength.
        scattered_power=False,
        progress=True,
    )


def panel(out_dir, kind, data, wavelengths_m, alpha, t, probe_nm, title, fname):
    """Spectrogram: absolute, per-timestep normalised, and alpha."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # spectra_from_phase_spaces returns the wavelength axis in SI metres, not the
    # nm the window was given in.
    nm = np.asarray(wavelengths_m) * 1e9
    data = np.asarray(data)
    t_ps = np.asarray(t) * 1e12

    fig, axes = plt.subplots(
        1, 3, figsize=(15.5, 4.4), gridspec_kw={"width_ratios": [1.25, 1.25, 0.8]}
    )

    finite = data[np.isfinite(data).all(axis=1)]
    im = axes[0].imshow(
        data.T, origin="lower", aspect="auto",
        extent=[t_ps[0], t_ps[-1], nm[0], nm[-1]], cmap="inferno",
        vmax=np.percentile(finite, 99) if finite.size else None,
    )
    fig.colorbar(im, ax=axes[0], label="S(k, w)")
    axes[0].set_title("absolute")

    # The scattered power climbs ~2 decades as the piston arrives, which would
    # otherwise saturate the late frames and black out the early ones.
    with np.errstate(invalid="ignore"):
        norm = data / np.nanmax(data, axis=1, keepdims=True)
    im = axes[1].imshow(
        norm.T, origin="lower", aspect="auto",
        extent=[t_ps[0], t_ps[-1], nm[0], nm[-1]], cmap="inferno", vmin=0.0, vmax=1.0,
    )
    fig.colorbar(im, ax=axes[1], label="S(k, w) / max per timestep")
    axes[1].set_title("normalised per timestep")

    for ax in axes[:2]:
        ax.axhline(probe_nm, color="cyan", ls="--", lw=0.8, alpha=0.8)
        ax.set_xlabel("time (ps)")
        ax.set_ylabel("wavelength (nm)")

    alpha = np.asarray(alpha)
    axes[2].semilogy(t_ps, alpha, lw=1.4)
    axes[2].axhline(1.0, color="0.4", ls=":", lw=1.0)
    axes[2].set_xlabel("time (ps)")
    axes[2].set_ylabel(r"$\alpha$")
    axes[2].set_title(f"scattering parameter (median {np.nanmedian(alpha):.2f})")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = os.path.join(out_dir, fname)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=os.path.join(ROOT, "runs", "R1"))
    ap.add_argument("--stage", choices=["auto", "read", "model"], default="auto",
                    help="'read' bins the plotfiles (needs yt), 'model' forward-models "
                         "them (needs torch). Default 'auto' does what this interpreter can.")
    ap.add_argument("--probe-wavelength", type=float, default=532.0, metavar="NM")
    ap.add_argument("--angle", type=float, default=90.0, metavar="DEG",
                    help="scattering angle (default 90)")
    ap.add_argument("--position", type=float, default=None, metavar="M",
                    help="sampling position in m (default: domain centre)")
    ap.add_argument("--ion-label", default="p+",
                    help="physical species the simulation ions stand for")
    ap.add_argument("--iaw-halfwidths", type=float, default=2.5, metavar="N",
                    help="IAW window half-width in units of the C_s Doppler shift")
    ap.add_argument("--wavelength-bins", type=int, default=400)
    ap.add_argument("--velocity-bins", type=int, default=512)
    ap.add_argument("--position-bins", type=int, default=256)
    ap.add_argument("--smoothing-window", type=int, default=31,
                    help="boxcar width in velocity bins. Above alpha ~ 1 the spectrum "
                         "carries |1 - chi_e/eps|^2, which diverges at the EPW resonance "
                         "and turns VDF shot noise into speckle, so collective runs need "
                         "markedly more smoothing than the pipeline's default of 9.")
    ap.add_argument("--smoothing-iterations", type=int, default=2)
    ap.add_argument("--cache-dir", default=None,
                    help="default runs/<ID>/thomson_cache")
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    cfg = kconfig.load(run_dir)
    run_id = cfg["meta"]["run_id"]
    scales = units.derive(cfg)
    cache_dir = args.cache_dir or os.path.join(run_dir, "thomson_cache")

    pt = _pt()
    has_yt, has_torch = _have("yt"), _have("torch")
    cached = all(
        os.path.exists(os.path.join(cache_dir, f"{n}.npz"))
        for n in species_of(cfg, "electron") + species_of(cfg, "ion")
    )

    stage = args.stage
    if stage == "auto":
        stage = "model" if (has_torch and (cached or has_yt)) else "read"
    if stage == "read" and not has_yt:
        raise SystemExit("stage 'read' needs yt; run this under the physics env.")
    if stage == "model" and not has_torch:
        raise SystemExit("stage 'model' needs torch; run this under the tsnn env.")

    print(f"{run_id}: {stage} stage  (yt={has_yt} torch={has_torch} cached={cached})")
    electrons, ions = read_phase_spaces(pt, run_dir, cfg, args, cache_dir)

    if stage == "read":
        print(f"  cached {len(electrons) + len(ions)} phase spaces in {cache_dir}")
        print("\nNow forward-model them under an env that has torch:\n"
              f"  /opt/anaconda3/envs/tsnn/bin/python scripts/make_thomson.py {args.run_dir}")
        return

    spec = model(pt, electrons, ions, cfg, scales, args)

    print(f"\n  alpha_epw: {np.nanmin(spec.alpha_epw):.3e} to {np.nanmax(spec.alpha_epw):.3e}")
    if spec.alpha_iaw is not None:
        print(f"  alpha_iaw: {np.nanmin(spec.alpha_iaw):.3e} to {np.nanmax(spec.alpha_iaw):.3e}")
    print(f"  n_e:       {np.nanmin(spec.electron_density):.3e} to "
          f"{np.nanmax(spec.electron_density):.3e} m^-3")
    if np.nanmedian(spec.alpha_epw) < 1.0:
        print("  NOTE alpha < 1: sub-collective, so the ion feature is weak and the IAW "
              "panel is dominated by the electron feature and any bulk drift.")

    out_dir = os.path.join(ROOT, "media", run_id)
    os.makedirs(out_dir, exist_ok=True)
    probe = args.probe_wavelength
    stamp = f"{probe:.0f} nm probe, {args.angle:.0f} deg"

    panel(out_dir, "epw", spec.epw, spec.epw_wavelengths, spec.alpha_epw, spec.t,
          probe, f"{run_id} — Thomson EPW (electron feature), {stamp}", "thomson_epw.png")
    if spec.iaw is not None:
        panel(out_dir, "iaw", spec.iaw, spec.iaw_wavelengths, spec.alpha_iaw, spec.t,
              probe, f"{run_id} — Thomson IAW (ion-acoustic feature), {stamp}",
              "thomson_iaw.png")

    npz = os.path.join(out_dir, "thomson_spectra.npz")
    np.savez_compressed(
        npz,
        t=np.asarray(spec.t),
        epw=np.asarray(spec.epw),
        epw_wavelengths=np.asarray(spec.epw_wavelengths),
        alpha_epw=np.asarray(spec.alpha_epw),
        iaw=np.asarray(spec.iaw) if spec.iaw is not None else np.zeros(0),
        iaw_wavelengths=(np.asarray(spec.iaw_wavelengths)
                         if spec.iaw_wavelengths is not None else np.zeros(0)),
        alpha_iaw=(np.asarray(spec.alpha_iaw)
                   if spec.alpha_iaw is not None else np.zeros(0)),
        electron_density=np.asarray(spec.electron_density),
        position=float(spec.position),
        probe_wavelength_nm=probe,
        scattering_angle_deg=args.angle,
    )
    print(f"  wrote {npz}")


if __name__ == "__main__":
    main()
