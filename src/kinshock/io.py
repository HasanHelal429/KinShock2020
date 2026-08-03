"""WarpX diagnostic I/O for KinShock2020 (yt-backed).

Designed after ``warpx-cda/heating_operator/scripts/make_shock_figures.py`` (same
yt loader, per-species macroparticle histogramming, field read) but standalone and
config-driven — it takes the run's :class:`kinshock.units.Scales` instead of module-
level constants. ``yt`` is imported lazily so the rest of the package imports without
it (e.g. for the units↔Table-I test).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from .units import QE, Scales


@dataclass
class Frame:
    """One plotfile: grid + fields + time (SI)."""
    path: str
    time: float                # [s]
    z_edges: np.ndarray        # cell edges along the propagation axis [m]
    z_centers: np.ndarray      # cell centers [m]
    Bx: np.ndarray
    By: np.ndarray
    Bz: np.ndarray
    ds: object                 # underlying yt dataset (for particle access)
    comps: dict = field(default_factory=dict)   # extra grid components (see load_frame)

    @property
    def Bperp(self):
        """Perpendicular field magnitude (the compressed component; deck uses Bx)."""
        return np.sqrt(self.Bx ** 2 + self.By ** 2)


def _yt():
    import yt
    yt.set_log_level(40)
    return yt


def plotfiles(run_dir, subdirs=("diags_movies", "diags")):
    """Sorted plotfile paths under a run dir. Tries movie diags first, then diags."""
    for sub, prefix in ((subdirs[0], "plt"), (subdirs[1], "diag1")):
        d = os.path.join(run_dir, sub)
        if os.path.isdir(d):
            pfs = sorted(os.path.join(d, p) for p in os.listdir(d)
                         if p.startswith(prefix))
            if pfs:
                return pfs
    raise FileNotFoundError(f"no plotfiles under {run_dir}/{{{','.join(subdirs)}}}")


def field_plotfiles(run_dir):
    """High-cadence field-only plotfiles (``diag_fields*``) if the run wrote them,
    else fall back to the Full ``diag1`` plotfiles. Lets field diagnostics (e.g. the
    B_perp streak) use the dense field series when available without breaking runs
    that predate the field-only diagnostic."""
    d = os.path.join(run_dir, "diags")
    if os.path.isdir(d):
        pfs = sorted(os.path.join(d, p) for p in os.listdir(d)
                     if p.startswith("diag_fields"))
        if pfs:
            return pfs
    return plotfiles(run_dir)


def load_frame(path, fields=()) -> Frame:
    """Load one plotfile into a :class:`Frame` (1D covering grid along z).

    ``fields`` names extra grid components to read into ``Frame.comps`` (e.g.
    ``("rho_amb_electrons",)``). Components the plotfile does not carry are simply
    absent from the dict, so callers can probe availability and fall back.
    """
    yt = _yt()
    ds = yt.load(path)
    dims = ds.domain_dimensions
    nz = int(dims[0])
    lo = float(ds.domain_left_edge[0])
    hi = float(ds.domain_right_edge[0])
    edges = np.linspace(lo, hi, nz + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # A covering_grid spanning the full domain can round up to one ULP past
    # domain_right_edge, and yt refuses to read outside a NON-PERIODIC domain
    # (RuntimeError, "attempted to read outside the boundaries"). That silently
    # zeroed Bx/By/Bz for every run whose right edge happened to round that way
    # -- and a zero B field looks like a physics result, not a read failure
    # (RESULTS 2026-08-03). force_periodicity only affects ghost-cell lookups at
    # the edge, never the cell values themselves.
    ds.force_periodicity()
    cg = ds.covering_grid(0, ds.domain_left_edge, dims)

    def comp_opt(name):
        """The component, or None if the plotfile genuinely does not carry it.

        Only a missing-field error may fall through to None. Anything else --
        notably the boundary RuntimeError above -- must raise, or a read failure
        gets laundered into a plausible-looking array of zeros.
        """
        if not any(f[1] == name for f in ds.field_list):
            return None
        return np.asarray(cg["boxlib", name]).reshape(nz)

    def comp(name):
        a = comp_opt(name)
        return np.zeros(nz) if a is None else a

    comps = {}
    for name in fields:
        a = comp_opt(name)
        if a is not None:
            comps[name] = a

    return Frame(path=path, time=float(ds.current_time),
                 z_edges=edges, z_centers=centers,
                 Bx=comp("Bx"), By=comp("By"), Bz=comp("Bz"), ds=ds, comps=comps)


def species_density(frame: Frame, species) -> np.ndarray:
    """Number density [m^-3] histogrammed from macroparticles for one or more species
    (1D cell volume = dz * 1 m^2)."""
    if isinstance(species, str):
        species = [species]
    edges = frame.z_edges
    nz = len(edges) - 1
    dz = (edges[-1] - edges[0]) / nz
    ad = frame.ds.all_data()
    n = np.zeros(nz)
    for sp in species:
        try:
            z = np.asarray(ad[(sp, "particle_position_x")])
            w = np.asarray(ad[(sp, "particle_weight")])
        except Exception:
            continue
        h, _ = np.histogram(z, bins=edges, weights=w)
        n += h / dz
    return n


def rho_field_names(species):
    """Per-species charge-density component names (``rho_<species>``) for load_frame."""
    if isinstance(species, str):
        species = [species]
    return tuple(f"rho_{sp}" for sp in species)


def field_species_density(frame: Frame, species, charge_states=None) -> np.ndarray:
    """Number density [m^-3] from the *deposited* ``rho_<species>`` grid components
    instead of the macroparticles (requires ``load_frame(..., fields=...)``).

    This is what lets a density streak run at the field-only cadence: ``diag_fields``
    sets ``write_species = 0``, so it carries no particles but does carry per-species
    rho. Values are shape-factor (and ``filter_npass``) smoothed, so they differ
    slightly from :func:`species_density` at the sharpest gradients.

    Raises KeyError if the frame carries none of the requested rho components, so
    callers can fall back to :func:`species_density` on the particle plotfiles.
    """
    if isinstance(species, str):
        species = [species]
    charge_states = charge_states or {}
    n = np.zeros(len(frame.z_centers))
    found = False
    for sp in species:
        rho = frame.comps.get(f"rho_{sp}")
        if rho is None:
            continue
        n += np.abs(rho) / (charge_states.get(sp, 1) * QE)
        found = True
    if not found:
        raise KeyError(f"no rho_<species> components for {species} in {frame.path}")
    return n


def species_phase(frame: Frame, species, scales: Scales, mass=None):
    """(z [m], u_z = p_z / (m c)) for one species. ``mass`` defaults to the ion mass
    for ion species and m_e otherwise (pass explicitly to be safe)."""
    from .units import ME
    m = mass if mass is not None else (scales.mi if "ion" in species else ME)
    ad = frame.ds.all_data()
    try:
        z = np.asarray(ad[(species, "particle_position_x")])
        pz = np.asarray(ad[(species, "particle_momentum_z")]) / (m * _C())
    except Exception:
        z, pz = np.array([]), np.array([])
    return z, pz


def species_phase_weighted(frame: Frame, species, scales: Scales, mass=None):
    """Like :func:`species_phase` but also returns the macroparticle weights, so a
    2D histogram over (z, u_z) is a true phase-space *distribution* f(z, u_z) rather
    than a raw macroparticle count. Returns (z [m], u_z = p_z/(m c), weight)."""
    from .units import ME
    m = mass if mass is not None else (scales.mi if "ion" in species else ME)
    ad = frame.ds.all_data()
    try:
        z = np.asarray(ad[(species, "particle_position_x")])
        pz = np.asarray(ad[(species, "particle_momentum_z")]) / (m * _C())
        w = np.asarray(ad[(species, "particle_weight")])
    except Exception:
        z, pz, w = np.array([]), np.array([]), np.array([])
    return z, pz, w


def _C():
    from .units import C
    return C


def reduced_diag(run_dir, name):
    """Load a WarpX reduced-diagnostic table (whitespace columns, '#'-commented header)
    e.g. ParticleEnergy 'EP' or ParticleNumber 'PN'. Returns (header_list, data ndarray)."""
    path = None
    for cand in (os.path.join(run_dir, "diags", "reducedfiles", name + ".txt"),
                 os.path.join(run_dir, "diags", name + ".txt"),
                 os.path.join(run_dir, name + ".txt")):
        if os.path.isfile(cand):
            path = cand
            break
    if path is None:
        raise FileNotFoundError(f"reduced diag '{name}' not found under {run_dir}")
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    return header, data
