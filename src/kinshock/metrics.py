"""Shock-formation metrics for KinShock2020 (Schaeffer 2020, Sec. III).

Implements the paper's quantitative diagnostics, operating on already-loaded
profiles / phase-space arrays (see :mod:`kinshock.io`) plus a
:class:`kinshock.units.Scales`. Nothing here reads files or knows about yt.

Covers: the plasma-speed model (Eqs. 1-2), the seven shock-formation criteria
(Sec. III B), the reflected-ambient-ion functions F/G and the t*/z* timescales
(Figs. 4, 6), and compression ratios (Sec. III C).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .units import Scales


# --------------------------------------------------------------------------- #
# Plasma-speed model (Sec. III A)
# --------------------------------------------------------------------------- #
def expansion_speed(n_e, n_e_ab, n_e0, Cs_ab):
    """Eq. 1: ablation expansion speed at density ``n_e`` (scale-free profile).

    v = (C_s,ab / 2) [ 1 - ln( (n_e - n_e0) / (n_e,ab - n_e0) ) ].
    Vectorized; returns NaN where the log argument is non-positive.
    """
    n_e = np.asarray(n_e, dtype=float)
    num = n_e - n_e0
    den = n_e_ab - n_e0
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = np.where((num > 0) & (den > 0), num / den, np.nan)
        return 0.5 * Cs_ab * (1.0 - np.log(arg))


def rh_speed_ratio(MA, beta0):
    """Eq. 2: perpendicular RH upstream/downstream speed ratio v1/v2.

    v1/v2 = -[M_A^2 + b0] + sqrt( (M_A^2 + b0)^2 + 8 M_A^2 ),  b0 = (5/2)(1 + beta0).
    """
    b0 = 2.5 * (1.0 + beta0)
    x = MA * MA + b0
    return -x + math.sqrt(x * x + 8.0 * MA * MA)


# --------------------------------------------------------------------------- #
# Compression ratios (Sec. III C)
# --------------------------------------------------------------------------- #
def compression(profile, upstream_value):
    """Peak compression ratio max(profile)/upstream_value (density or |B|)."""
    profile = np.asarray(profile, dtype=float)
    if upstream_value == 0:
        return np.nan
    return float(np.nanmax(profile) / upstream_value)


# --------------------------------------------------------------------------- #
# Front tracking / speeds
# --------------------------------------------------------------------------- #
def track_front(zc, n, upstream_density, threshold=1.5, z_exclude=None, side=+1):
    """Outermost position where n exceeds ``threshold * upstream_density``.

    ``zc`` and returned position share units. ``z_exclude`` masks out the piston
    slab (|z| < z_exclude); ``side`` selects the +z (default) or -z front.
    Returns np.nan if no cell qualifies.
    """
    zc = np.asarray(zc, dtype=float)
    n = np.asarray(n, dtype=float)
    mask = n > threshold * upstream_density
    if z_exclude is not None:
        mask &= np.abs(zc) > z_exclude
    mask &= (zc > 0) if side > 0 else (zc < 0)
    if not mask.any():
        return np.nan
    return float(zc[mask].max() if side > 0 else zc[mask].min())


def speed_from_trajectory(t, z, use_second_half=True):
    """Linear-fit speed dz/dt (units follow inputs). Fits the second half by default
    (after the formation transient)."""
    t = np.asarray(t, dtype=float)
    z = np.asarray(z, dtype=float)
    good = np.isfinite(t) & np.isfinite(z)
    t, z = t[good], z[good]
    if t.size < 2:
        return np.nan
    if use_second_half:
        h = t.size // 2
        t, z = t[h:], z[h:]
    return float(np.polyfit(t, z, 1)[0])


# --------------------------------------------------------------------------- #
# Reflected-ambient-ion functions F, G and timescales (Figs. 4, 6)
# --------------------------------------------------------------------------- #
def reflected_fraction_G(uz_ambient, vsh):
    """G(t) = N_reflected / N_total: fraction of ambient ions with v_z > vsh
    (faster than the shock, i.e. magnetically reflected in the upstream frame)."""
    uz = np.asarray(uz_ambient, dtype=float)
    if uz.size == 0:
        return 0.0
    return float(np.count_nonzero(uz > vsh) / uz.size)


def reflected_profile_F(z_ambient, uz_ambient, vsh, edges):
    """F(z,t) = spatial distribution of reflected ambient ions, normalized by the
    total ambient-ion count (histogram over ``edges``). Returns (F, bin_centers)."""
    z = np.asarray(z_ambient, dtype=float)
    uz = np.asarray(uz_ambient, dtype=float)
    edges = np.asarray(edges, dtype=float)
    ntot = max(z.size, 1)
    sel = uz > vsh
    counts, _ = np.histogram(z[sel], bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return counts.astype(float) / ntot, centers


def onset_time_from_G(t, G):
    """t*_1: time of maximum dG/dt (onset of shock formation). Returns (t_star, index)."""
    t = np.asarray(t, dtype=float)
    G = np.asarray(G, dtype=float)
    if t.size < 3:
        return np.nan, -1
    dG = np.gradient(G, t)
    i = int(np.nanargmax(dG))
    return float(t[i]), i


def onset_location_from_F(z, F):
    """z*_1: location of maximum dF/dz at the onset time. Returns (z_star, index)."""
    z = np.asarray(z, dtype=float)
    F = np.asarray(F, dtype=float)
    if z.size < 3:
        return np.nan, -1
    dF = np.gradient(F, z)
    i = int(np.nanargmax(dF))
    return float(z[i]), i


# --------------------------------------------------------------------------- #
# The seven shock-formation criteria (Sec. III B)
# --------------------------------------------------------------------------- #
@dataclass
class CriteriaResult:
    values: dict          # measured numbers behind each criterion
    flags: dict           # bool per criterion 1..7
    is_precursor: bool    # criteria 1-6
    is_shock: bool        # criteria 1-7


def evaluate_criteria(zc, n_e, Bmag, uz_ambient, z_ambient, scales: Scales,
                      vsh, v_front, piston_field_z, front_z,
                      lambda_ii_over_di0=None):
    """Evaluate the seven criteria at one time (all inputs in SI; zc/front_z/piston_field_z
    in the same length unit).

    Parameters mirror the paper:
      * ``vsh``/``v_front`` : shock speed and local front speed [m/s]
      * ``piston_field_z``  : position of the piston's peak field
      * ``front_z``         : position of the shock/precursor field peak
      * ``lambda_ii_over_di0``: ion-ion mfp / d_i0 (defaults to config-scale ~350 if None)
    """
    n_e = np.asarray(n_e, dtype=float)
    Bmag = np.asarray(Bmag, dtype=float)

    n_comp = compression(n_e, scales.namb)
    b_comp = compression(Bmag, scales.B0)
    G = reflected_fraction_G(uz_ambient, vsh)

    # steep ramp: characteristic gradient scale of |B| near the front, in d_i0
    ramp_di0 = _ramp_scale(zc, Bmag) / scales.di0 if scales.di0 else np.inf
    Mms_front = v_front / math.sqrt(scales.vA ** 2 + scales.Cs0 ** 2)
    lam = lambda_ii_over_di0 if lambda_ii_over_di0 is not None else 350.0
    sep = (abs(front_z - piston_field_z) >= 0.25 * scales.rho_i0)

    values = {
        "M_ms_front": Mms_front,
        "lambda_ii_over_di0": lam,
        "n_compression": n_comp,
        "B_compression": b_comp,
        "ramp_scale_over_di0": ramp_di0,
        "reflected_fraction_G": G,
        "piston_separation_over_rho_i0": abs(front_z - piston_field_z) / scales.rho_i0,
    }
    flags = {
        "1_super_magnetosonic": Mms_front > 1.0,
        "2_collisionless": lam > 1.0,
        "3_density_compression": n_comp > 2.0,
        "4_field_compression": b_comp > 2.0,
        "5_steep_ramp": ramp_di0 <= 1.0,
        "6_reflected_ions": G > 0.0,
        "7_piston_separation": sep,
    }
    is_precursor = all(flags[k] for k in list(flags)[:6])
    is_shock = is_precursor and flags["7_piston_separation"]
    return CriteriaResult(values=values, flags=flags,
                          is_precursor=is_precursor, is_shock=is_shock)


def _ramp_scale(zc, field):
    """Gradient length L = field_peak / max|dfield/dz| near the peak (same units as zc)."""
    zc = np.asarray(zc, dtype=float)
    f = np.asarray(field, dtype=float)
    if zc.size < 3:
        return np.inf
    g = np.gradient(f, zc)
    gmax = np.nanmax(np.abs(g))
    fmax = np.nanmax(np.abs(f))
    if gmax == 0:
        return np.inf
    return float(fmax / gmax)
