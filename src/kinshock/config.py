"""Load and validate a per-run KinShock2020 config (see REPLICATION_PLAN.md §6.0a).

The config holds only PRIMARY quantities; :mod:`kinshock.units` derives the rest.
:func:`validate` checks the derived dimensionless numbers against the run's Table I
``targets`` and returns human-readable warnings (it never raises on a physics
mismatch — the deck may legitimately explore off-Table-I points).
"""

from __future__ import annotations

import os
from typing import Any

import yaml

from . import units


def load(path: str) -> dict:
    """Load a config from a YAML file or a run directory (which must contain config.yaml)."""
    if os.path.isdir(path):
        path = os.path.join(path, "config.yaml")
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_path"] = os.path.abspath(path)
    cfg["_run_dir"] = os.path.dirname(cfg["_path"])
    _require(cfg, ["reference", "plasma", "field", "geometry", "numerics", "model"])
    return cfg


def derive(cfg: dict) -> units.Scales:
    """Convenience: config dict -> derived :class:`kinshock.units.Scales`."""
    return units.derive(cfg)


def validate(cfg: dict, scales: units.Scales | None = None,
             rel_tol: float = 0.15) -> list[str]:
    """Compare derived dimensionless numbers to the config's Table I ``targets``.

    Returns a list of warning strings (empty = all within tolerance). M_A and M_ms
    use ``rel_tol``; the betas are checked only to within a factor of ~2.5, since
    the paper's normalization (n_e,ab = 1.25, reduced c) shifts them by ~2x
    relative to our SI formula while keeping the same regime (see config comments).
    """
    if scales is None:
        scales = units.derive(cfg)
    warns: list[str] = []
    tgt = cfg.get("targets", {})

    def check_rel(name, got, want, tol):
        if want in (None, 0):
            return
        rel = abs(got - want) / abs(want)
        if rel > tol:
            warns.append(f"{name}: derived {got:.4g} vs target {want:.4g} "
                         f"(rel diff {rel*100:.0f}% > {tol*100:.0f}%)")

    def check_factor(name, got, want, factor):
        if want in (None, 0):
            return
        r = got / want
        if r < 1.0 / factor or r > factor:
            warns.append(f"{name}: derived {got:.4g} vs target {want:.4g} "
                         f"(factor {r:.2g}, outside {factor:g}x)")

    check_rel("M_A", scales.MA, tgt.get("M_A"), rel_tol)
    check_rel("M_ms", scales.Mms, tgt.get("M_ms"), max(rel_tol, 0.2))
    check_factor("beta_ab", scales.beta_ab, tgt.get("beta_ab"), 2.5)
    check_factor("beta_0", scales.beta_0, tgt.get("beta_0"), 2.5)

    # numerical sanity
    if scales.dt_wpe > 0.5:
        warns.append(f"dt*wpe = {scales.dt_wpe:.3g} > 0.5 (resolve omega_pe: lower CFL or dz)")
    front_full = 6.0 * scales.rho_i0
    if scales.domain_halfwidth < front_full:
        warns.append(f"domain half-width {scales.domain_halfwidth/scales.de:.0f} d_e may be too "
                     f"small: shock reaches ~6 rho_i0 = {front_full/scales.de:.0f} d_e by t*_3")

    # geometry <-> boundary consistency
    geo = cfg["geometry"]
    b = geo.get("boundary", "periodic")
    faces = (b.get("lo", "periodic"), b.get("hi", "periodic")) if isinstance(b, dict) else (b, b)
    n_periodic = sum(f == "periodic" for f in faces)
    if n_periodic == 1:
        warns.append("boundary: exactly one face is 'periodic' — WarpX requires both faces "
                     "periodic together or neither; set a non-periodic pair (e.g. "
                     "lo: reflecting, hi: absorbing)")
    if geo.get("layout", "symmetric") == "one_sided" and faces[0] == "periodic":
        warns.append("one-sided layout with a periodic z=0 face: the center should be a "
                     "reflecting (symmetry/foil) wall — set boundary.lo: reflecting")
    # collisions: Takizuka-Abe/Perez scattering assumes the collision time is
    # resolved (nu*dt << 1) and needs enough pairs per cell to be a fair sample.
    if scales.nu_ei_dt is not None:
        if scales.nu_ei_dt > 0.1:
            warns.append(f"collisions: nu_ei*dt = {scales.nu_ei_dt:.3g} > 0.1 — the collision "
                         f"time is under-resolved; raise ndt_subcycle or lower the target "
                         f"collisionality")
        if scales.mfp_ei_ab < scales.dz:
            warns.append(f"collisions: mfp_ei,ab = {scales.mfp_ei_ab/scales.dz:.3g} dz — the "
                         f"mean free path is below the cell size (unresolved transport)")
        units.collision_pairs(cfg)   # raises on an unknown species name

    has_bfield = (cfg["field"].get("orientation", "perpendicular") != "none"
                  and float(cfg["field"].get("B0_tesla", 0) or 0) != 0)
    if has_bfield and "absorbing" in faces:
        warns.append("boundary 'absorbing' (Silver-Mueller) is incompatible with the B-field "
                     "divergence cleaner that runs when a background B is set; use 'open' "
                     "(pec fields + absorbing particles) instead")
    return warns


def _require(cfg: dict, keys: list[str]) -> None:
    missing = [k for k in keys if k not in cfg]
    if missing:
        raise KeyError(f"config missing required section(s): {missing}")
