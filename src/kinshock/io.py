"""WarpX diagnostic I/O for KinShock2020 (yt-backed).

Designed after ``warpx-cda/heating_operator/scripts/make_shock_figures.py`` (same
yt loader, per-species macroparticle histogramming, field read) but standalone and
config-driven — it takes the run's :class:`kinshock.units.Scales` instead of module-
level constants. ``yt`` is imported lazily so the rest of the package imports without
it (e.g. for the units↔Table-I test).
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

from .units import Scales


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


def load_frame(path) -> Frame:
    """Load one plotfile into a :class:`Frame` (1D covering grid along z)."""
    yt = _yt()
    ds = yt.load(path)
    dims = ds.domain_dimensions
    nz = int(dims[0])
    lo = float(ds.domain_left_edge[0])
    hi = float(ds.domain_right_edge[0])
    edges = np.linspace(lo, hi, nz + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    cg = ds.covering_grid(0, ds.domain_left_edge, dims)

    def comp(name):
        try:
            return np.asarray(cg["boxlib", name]).reshape(nz)
        except Exception:
            return np.zeros(nz)

    return Frame(path=path, time=float(ds.current_time),
                 z_edges=edges, z_centers=centers,
                 Bx=comp("Bx"), By=comp("By"), Bz=comp("Bz"), ds=ds)


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
