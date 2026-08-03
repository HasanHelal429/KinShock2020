#!/usr/bin/env python3
"""Reproduce the Schaeffer 2020 shock diagnostics for a KinShock2020 run.

Config-driven (all scales from runs/RXX/config.yaml via kinshock.units); writes
to media/<run_id>/. Implements analyses A–E of REPLICATION_PLAN.md §6:

  A. B_perp(z,t) streak plot with piston/shock speed lines  -> shock_streak.png
     (uses the high-cadence field-only diag_fields series when present, else diag1)
  A. shock-front trajectory + measured v_sh, M_A, M_ms       -> shock_trajectory.png
  B. ion (z, u_z) phase-space distribution at several times  -> shock_phase.png (Fig. 5)
     + density / B_perp line-outs                            -> shock_lineouts.png
     + species-resolved 3-row phase space (amb/piston/e-)    -> shock_fig7.png (Fig. 7)
  C. seven shock-formation criteria per frame                -> criteria.json (+ stdout)
  D. reflected-ambient-ion fraction G(t), F(z), t*_1, z*_1   -> shock_reflected.png

Usage:  python scripts/make_figures.py [runs/R1] [--nframes 5]
"""

from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import kinshock  # noqa: E402
from kinshock import io, metrics, plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def _json_native(o):
    """Cast numpy scalars (bool_/float64/int64) to native Python types for json."""
    if hasattr(o, "item"):
        return o.item()
    raise TypeError(f"not JSON serializable: {type(o).__name__}")


def _select(pfs, n):
    if len(pfs) <= n:
        return list(range(len(pfs)))
    return list(np.linspace(0, len(pfs) - 1, n).astype(int))


def _ion_density(frame, cfg):
    return io.species_density(frame, cfg["ion_species"])


def load_series(pfs):
    """Load a list of plotfile paths into Frames (fields + time; particles lazy)."""
    return [io.load_frame(pf) for pf in pfs]


# --- the single source of truth for the shock speed + front trajectory ---
class _Shock:
    """v_sh and z_front(frame) for the whole figure suite. Wraps the by-eye fit
    (runs/<ID>/shock_fit.yaml via scripts/tune_shock.py) when present, else the
    automatic track_front + speed_from_trajectory fallback. Every diagnostic reads
    v_sh and the front through this one object, so they can no longer drift apart."""

    def __init__(self, v_sh, sc, cfg, fit=None):
        self.v_sh, self.sc, self.cfg, self.fit = v_sh, sc, cfg, fit
        self.from_fit = fit is not None
        self.MA = v_sh / sc.vA
        self.Mms = v_sh / np.sqrt(sc.vA ** 2 + sc.Cs0 ** 2)

    def front_m(self, fr):
        """Forward (+z) shock-front position [m] for frame ``fr`` (np.nan if none)."""
        if self.fit is not None:
            return self.fit.z_front(fr.time, fr.time * self.sc.wci0)
        zexcl = self.cfg["geometry"]["slab_halfwidth_di"] * self.sc.di
        return metrics.track_front(fr.z_centers, _ion_density(fr, self.cfg),
                                   self.sc.namb, threshold=1.5, z_exclude=zexcl, side=+1)


def resolve_shock(run_dir, frames, cfg, sc):
    """Build the run's :class:`_Shock`: prefer the by-eye fit, else auto-fallback."""
    fit = metrics.load_shock_fit(run_dir, sc)
    if fit is not None:
        print(f"  shock: by-eye fit (shock_fit.yaml)  v_sh={fit.v_sh/kinshock.units.C:.4f}c  "
              f"M_A={fit.MA:.2f}  M_ms={fit.Mms:.2f}")
        return _Shock(fit.v_sh, sc, cfg, fit=fit)
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di
    ts, zs = [], []
    for fr in frames:
        zf = metrics.track_front(fr.z_centers, _ion_density(fr, cfg), sc.namb,
                                 threshold=1.5, z_exclude=zexcl)
        if np.isfinite(zf):
            ts.append(fr.time); zs.append(zf)
    if len(ts) >= 2:
        z_edge = cfg["geometry"]["domain_halfwidth_de"] * sc.de
        v = metrics.speed_from_trajectory(np.array(ts), np.array(zs),
                                          z_edge=z_edge, use_second_half=False)
    else:
        v = sc.vsh_model
    v = v if np.isfinite(v) else sc.vsh_model
    print(f"  shock: NO shock_fit.yaml -> auto track_front fallback (v_sh={v/kinshock.units.C:.4f}c). "
          f"Run scripts/tune_shock.py for a by-eye fit.")
    return _Shock(v, sc, cfg, fit=None)


# --- A: streak + trajectory ---
def fig_streak(frames, cfg, sc):
    t = np.array([fr.time for fr in frames]) * sc.wci0
    zc = np.asarray(frames[0].z_centers) / sc.di0
    B = np.array([fr.Bx / sc.B0 for fr in frames])  # (nt, nz)
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    pc = ax.pcolormesh(t, zc, B.T, shading="auto", cmap="viridis",
                       vmin=0, vmax=max(2.0, np.nanpercentile(B, 99)))
    # overlay model piston/shock speed lines (z = v t), converted to (t*wci0, z/d_i0)
    tt = np.linspace(t.min(), t.max(), 50)
    for v, lbl, ls in ((sc.vp_model, "v_p", "--"), (sc.vsh_model, "v_sh", "-")):
        z_di0 = (v * (tt / sc.wci0)) / sc.di0
        ax.plot(tt, z_di0, ls, color="w", lw=1.2, label=lbl)
    ax.set_ylim(0, zc.max())
    ax.set_xlabel(r"$t\,\omega_{ci0}$")
    ax.set_ylabel(r"$z / d_{i0}$")
    ax.legend(frameon=False, labelcolor="w", loc="upper left")
    fig.colorbar(pc, ax=ax, label=r"$B_\perp / B_0$")
    P.stamp(ax, cfg, sc)
    fig.suptitle(f"{cfg['meta']['run_id']}: magnetic streak (piston + shock)")
    fig.tight_layout()
    return P.savefig(fig, "shock_streak.png", run_id=cfg["meta"]["run_id"])


def fig_trajectory(frames, cfg, sc, shock):
    """Empirical front points (track_front) with the shock model overlaid. v_sh / M_A
    come from ``shock`` -- the by-eye fit (shock_fit.yaml) or the auto fallback -- so
    this figure reports the SAME speed every other diagnostic uses."""
    ts, zs = [], []
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di   # exclude piston slab
    for fr in frames:
        n = _ion_density(fr, cfg)
        zf = metrics.track_front(fr.z_centers, n, sc.namb, threshold=1.5, z_exclude=zexcl)
        if np.isfinite(zf):
            ts.append(fr.time)
            zs.append(zf)
    ts, zs = np.array(ts), np.array(zs)
    C = kinshock.units.C
    vsh, MA, Mms = shock.v_sh, shock.MA, shock.Mms
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    if ts.size:
        ax.plot(ts * sc.wci0, zs / sc.di0, "o", color=P.C_PISTON, ms=3,
                label="tracked front (n>1.5 n_e0)")
    # overlay the shock MODEL used everywhere else: fit -> z0 + v_sh*t (+ overrides);
    # fallback -> a line of the fitted slope through the tracked points.
    tw_all = np.array([fr.time * sc.wci0 for fr in frames])
    tgrid = np.linspace(tw_all.min(), tw_all.max(), 100)
    if shock.from_fit:
        zmodel = (shock.fit.z0 + vsh * (tgrid / sc.wci0)) / sc.di0
        if shock.fit.fronts:
            fk = np.array(sorted(shock.fit.fronts))
            ax.plot(fk, np.array([shock.fit.fronts[k] for k in fk]) / sc.di0,
                    "x", color="cyan", ms=9, mew=2.0, label="per-time override")
        lbl = "by-eye fit"
    else:
        z0 = float(np.median(zs - vsh * ts)) if ts.size else 0.0
        zmodel = (z0 + vsh * (tgrid / sc.wci0)) / sc.di0
        lbl = "auto fit (no shock_fit.yaml)"
    ax.plot(tgrid, zmodel, "-", color=P.C_REF, lw=1.6, label=lbl)
    ax.set_xlabel(r"$t\,\omega_{ci0}$")
    ax.set_ylabel(r"shock front $z / d_{i0}$")
    P.style_axes(ax)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    txt = (rf"$v_{{sh}}={vsh/C:.4f}\,c$" + "\n"
           rf"$M_A={MA:.1f}$ (target {cfg['targets']['M_A']})" + "\n"
           rf"$M_{{ms}}={Mms:.1f}$ (target {cfg['targets']['M_ms']})" + "\n"
           rf"$v_{{sh}}/v_p={vsh/sc.vp_model:.2f}$ (RH $\to$ 4/3)")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec=P.C_REF, alpha=0.9))
    src = "by-eye fit" if shock.from_fit else "AUTO fallback"
    fig.suptitle(f"{cfg['meta']['run_id']}: shock-front trajectory ({src})")
    fig.tight_layout()
    P.savefig(fig, "shock_trajectory.png", run_id=cfg["meta"]["run_id"])
    print(f"  v_sh = {vsh/C:.4f} c ; M_A = {MA:.2f} ; M_ms = {Mms:.2f}  [{src}]")


# --- B: line-outs + phase space ---
def fig_lineouts(frames, cfg, sc, nframes):
    idx = _select(frames, nframes)
    cmap = plt.get_cmap("viridis")
    fig, (a_n, a_b) = plt.subplots(2, 1, figsize=(7.6, 6.6), sharex=True)
    for k, i in enumerate(idx):
        fr = frames[i]
        zc = np.asarray(fr.z_centers) / sc.di0
        n = _ion_density(fr, cfg)
        col = cmap(k / max(len(idx) - 1, 1))
        lbl = rf"$t\,\omega_{{ci0}}={fr.time*sc.wci0:.1f}$"
        a_n.plot(zc, np.where(n > 0, n / sc.namb, np.nan), color=col, lw=1.1, label=lbl)
        a_b.plot(zc, fr.Bperp / sc.B0, color=col, lw=1.1)
    a_n.axhline(1.0, color=P.C_REF, ls="--", lw=0.8)
    a_n.set_yscale("log")
    a_n.set_ylabel(r"$n_i / n_{e0}$")
    a_b.set_ylabel(r"$B_\perp / B_0$")
    a_b.set_xlabel(r"$z / d_{i0}$")
    for ax in (a_n, a_b):
        P.style_axes(ax)
    a_n.legend(frameon=False, fontsize=8, ncol=2)
    P.stamp(a_n, cfg, sc)
    fig.suptitle(f"{cfg['meta']['run_id']}: density compression + magnetic ramp")
    fig.tight_layout()
    return P.savefig(fig, "shock_lineouts.png", run_id=cfg["meta"]["run_id"])


# phase panels follow the FORWARD (+z) shock: the window is anchored on the tracked
# front and sized to the feature -- extents scale with the front position so earlier
# (more compact) times are zoomed tighter -- and clipped to z >= 0. rho_i0 floors keep
# the earliest panels from collapsing to a degenerate width.
PHASE_WIN_BACK_FRAC = 0.5       # downstream extent behind the front ~ frac * z_front
PHASE_WIN_AHEAD_FRAC = 0.20     # upstream (foot) extent ahead of the front
PHASE_WIN_BACK_MIN_RHO = 0.6    # floors, in rho_i0
PHASE_WIN_AHEAD_MIN_RHO = 0.05
# left (v_z/v_sh) axis: NOT centred on 0 -- the forward shock is dominated by positive
# v_z, so the range is weighted upward to best capture the piston/reflected dynamics.
PHASE_V_LO, PHASE_V_HI = -1.0, 3.0


def _electron_species(cfg):
    return [n for n, s in cfg["species"].items() if s.get("kind") == "electron"]


def _smooth(a, k=7):
    """Light boxcar smoothing for the overlaid line-outs (beats per-cell shot noise)."""
    a = np.asarray(a, dtype=float)
    if a.size < k or k < 2:
        return a
    return np.convolve(a, np.ones(k) / k, mode="same")


def _select_nonzero(frames, n, sc):
    """Indices of ``n`` frames evenly spanning the run but EXCLUDING t=0 -- nothing has
    happened at t=0, so Fig-5-style panels start once the piston is moving."""
    cand = [i for i, fr in enumerate(frames) if fr.time * sc.wci0 > 1e-6]
    if not cand:
        cand = list(range(len(frames)))
    if len(cand) <= n:
        return cand
    return [cand[j] for j in np.linspace(0, len(cand) - 1, n).astype(int)]


def _select_times(frames, times, sc):
    """Indices of the frames NEAREST each requested time (in t*wci0), preserving the
    requested order. Used by --phase-times to hand-pick the phase-space panels."""
    tw = np.array([fr.time * sc.wci0 for fr in frames])
    idx = []
    for t in times:
        j = int(np.argmin(np.abs(tw - t)))
        if j not in idx:            # skip duplicates when two requests snap to one frame
            idx.append(j)
        print(f"  phase panel: requested t*wci0={t:.2f} -> nearest frame t*wci0={tw[j]:.2f}")
    return idx


def _shock_window(fr, cfg, sc, shock):
    """(z_lo, z_hi, z_front) in d_i0: a window that FOLLOWS the forward (+z) shock front,
    sized to the feature (extents grow with the front position; rho_i0 floors) and clipped
    to z >= 0. The front position comes from ``shock`` (the by-eye fit or auto fallback),
    so the phase panels track the identical front the trajectory/reflected diagnostics use."""
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di
    zf = shock.front_m(fr)
    rho = sc.rho_i0 / sc.di0
    zhi = float(np.asarray(fr.z_centers).max()) / sc.di0
    if np.isfinite(zf):
        c = zf / sc.di0
        back = max(PHASE_WIN_BACK_MIN_RHO * rho, PHASE_WIN_BACK_FRAC * c)
        ahead = max(PHASE_WIN_AHEAD_MIN_RHO * rho, PHASE_WIN_AHEAD_FRAC * c)
        lo, hi, zfront = c - back, c + ahead, c
    else:
        lo, hi, zfront = 0.0, min(zexcl / sc.di0 + rho, zhi), None
    return max(lo, 0.0), min(hi, zhi), zfront


def fig_phase(frames, cfg, sc, shock, nframes, times=None):
    idx = _select_times(frames, times, sc) if times else _select_nonzero(frames, nframes, sc)
    fig, axes = plt.subplots(1, len(idx), figsize=(3.4 * len(idx), 4.0), sharey=True)
    axes = np.atleast_1d(axes)
    vnorm = shock.v_sh
    v_edges = np.linspace(PHASE_V_LO, PHASE_V_HI, 201)
    e_species = _electron_species(cfg)
    ks = max(5, int(round(0.30 * sc.di0 / sc.dz)))   # smooth overlays over ~0.30 d_i0

    wins = {i: _shock_window(frames[i], cfg, sc, shock) for i in idx}

    # Right-axis scale for the overlaid B_x/B0 + n_e/n_e0 profiles: anchor on the
    # (smoothed) magnetic ramp B_x, O(4-6) at the shock, using panels whose window has
    # cleared the dense piston. The total electron density (piston ejecta filling the
    # downstream) then clips at the top, so its jump across the shock stays legible.
    slab_di0 = (cfg["geometry"]["slab_halfwidth_di"] * sc.di) / sc.di0
    piston_margin = max(3.0, 4.0 * slab_di0)

    def _win_bmax(i):
        fr = frames[i]
        lo, hi, _ = wins[i]
        zc = np.asarray(fr.z_centers) / sc.di0
        w = (zc >= lo) & (zc <= hi)
        return float(np.nanmax(_smooth((fr.Bx / sc.B0)[w], ks))) if w.any() else 0.0

    clean = [i for i in idx if wins[i][0] > piston_margin]
    pool = clean if clean else list(idx)
    prof_max = float(np.clip(1.3 * max(_win_bmax(i) for i in pool), 4.0, 14.0))
    # align the profile baseline (value 0) with the phase-space v_z = 0 line
    f0 = (0.0 - PHASE_V_LO) / (PHASE_V_HI - PHASE_V_LO)
    r_lo = -f0 / (1.0 - f0) * prof_max

    for k, (ax, i) in enumerate(zip(axes, idx)):
        fr = frames[i]
        lo, hi, zfront = wins[i]
        z_edges = np.linspace(lo, hi, 241)
        sd = {}
        for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
            z, uz, w = io.species_phase_weighted(fr, sp, sc, mass=sc.mi)
            if len(z):
                sd[key] = (z / sc.di0, uz * kinshock.units.C / vnorm, w)
        P.phase_distribution(ax, sd, z_edges, v_edges, vline=1.0, legend=(k == 0))

        # overlay B_x/B0 and n_e/n_e0 line-outs on a twin axis (Schaeffer 2020 Fig. 5),
        # semi-transparent + baseline aligned to v_z=0 so they guide without dominating
        zc = np.asarray(fr.z_centers) / sc.di0
        m = (zc >= lo) & (zc <= hi)
        ne = io.species_density(fr, e_species)
        ax2 = ax.twinx()
        ax2.plot(zc[m], _smooth((fr.Bx / sc.B0)[m], ks), color="#eaeaff", lw=1.1,
                 alpha=0.6, label=r"$B_x/B_0$")
        ax2.plot(zc[m], _smooth(np.where(ne > 0, ne / sc.namb, 0.0)[m], ks),
                 color="#7CFC5A", lw=1.1, alpha=0.6, label=r"$n_e/n_{e0}$")
        ax2.set_ylim(r_lo, prof_max)
        ax2.set_yticks(np.arange(0, prof_max + 0.1, 2.0))
        if zfront is not None:
            ax.axvline(zfront, color="w", ls="-", lw=0.6, alpha=0.35)
        if k == len(idx) - 1:
            ax2.set_ylabel(r"$B_x/B_0,\ \ n_e/n_{e0}$")
            ax2.legend(loc="upper right", fontsize=7, framealpha=0.35,
                       facecolor=P.PHASE_BG, labelcolor="w", edgecolor="none")
        else:
            ax2.set_yticklabels([])
        ax.set_title(rf"$t\,\omega_{{ci0}}={fr.time*sc.wci0:.1f}$", fontsize=9)
        ax.set_xlabel(r"$z / d_{i0}$")
    axes[0].set_ylabel(r"$v_z / v_{sh}$")
    fig.suptitle(f"{cfg['meta']['run_id']}: ion phase space following the shock "
                 r"($v_z>v_{sh}$ = reflected; $B_x$, $n_e$ overlaid)")
    fig.tight_layout()
    return P.savefig(fig, "shock_phase.png", run_id=cfg["meta"]["run_id"])


# --- B (Fig. 7): species-resolved phase-space rows ---
# Schaeffer 2020 Fig. 7 is a 3-row x N-timestep grid: row 1 ambient-ion (z, v_z),
# row 2 piston-ion (z, v_z) -- both with v_z relative to v_sh -- and row 3 electron
# (z, v_z) with B_x (black) and total n_e (red), relative to their upstream values,
# overlaid. Columns share a per-time window that follows the forward shock (identical
# to fig_phase, so the two figures register). Electrons are hot & bidirectional, so
# their row gets its own symmetric, data-sized v-axis.
#
# Unlike fig_phase (additive tints on a black panel), Fig. 7 shows ONE species per
# panel, so it is rendered on a WHITE background through a perceptually-uniform
# colormap (empty bins stay white); this also lets B_x be black, as in the paper.
FIG7_E_VPCTL = 99.5     # electron v-axis half-range from this |v_z|/v_sh percentile
FIG7_E_VMIN = 3.0       # ... clipped to a sane band (units of v_sh)
FIG7_E_VMAX = 10.0
FIG7_CMAP = "cividis"   # perceptually-uniform, CVD-safe (dark blue -> yellow)

# Horizontal-axis normalizations Fig. 7 can be drawn in. The windows (and every other
# figure) are computed in d_i0, so each entry only needs the length in metres, the axis
# label, and a filename suffix; the panels are rescaled by d_i0/length at plot time.
#   d_i0   -- upstream ion inertial length  c/omega_pi0      (default, as in fig_phase)
#   rho_i0 -- upstream ion gyroradius       v_p/omega_ci0    (the paper's z* normalization)
# ``note`` (when set) adds a normalization line to the suptitle -- used for the
# non-default axis, so the default d_i0 figure stays exactly as before.
FIG7_XUNITS = {
    "d_i0":   dict(attr="di0",    label=r"$z / d_{i0}$",     suffix="", note=None),
    "rho_i0": dict(attr="rho_i0", label=r"$z / \rho_{i0}$",  suffix="_rho_i0",
                   note=r"$\rho_{i0}=v_p/\omega_{ci0}$"),
}


def _phase_hist_cmap(ax, z, v, w, z_edges, v_edges, cmap=FIG7_CMAP):
    """Render a single-species phase-space density f(z, v) through a sequential colormap
    on a WHITE background: empty bins stay white, an asinh stretch keeps both the dense
    core and the sparse reflected beam visible. Returns nothing (draws into ``ax``)."""
    ax.set_facecolor("white")
    ax.set_xlim(z_edges[0], z_edges[-1])
    ax.set_ylim(v_edges[0], v_edges[-1])
    if len(z) == 0:
        return
    H, _, _ = np.histogram2d(z, v, bins=[z_edges, v_edges], weights=w)
    inten = np.ma.masked_less_equal(P._asinh_norm(H).T, 0.0)   # (nv, nz), empty -> masked
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("white")
    ax.imshow(inten, origin="lower", aspect="auto", interpolation="nearest",
              extent=[z_edges[0], z_edges[-1], v_edges[0], v_edges[-1]],
              cmap=cm, vmin=0.0, vmax=1.0)


def _electron_phase_weighted(fr, cfg, sc):
    """Combined (z, u_z, w) phase-space arrays for ALL electron species (row 3 shows the
    total electron distribution, matching the total n_e overlay)."""
    from kinshock.units import ME
    zs, us, ws = [], [], []
    for sp in _electron_species(cfg):
        z, uz, w = io.species_phase_weighted(fr, sp, sc, mass=ME)
        if len(z):
            zs.append(z); us.append(uz); ws.append(w)
    if not zs:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(zs), np.concatenate(us), np.concatenate(ws)


def fig_fig7(frames, cfg, sc, shock, nframes, times=None, xunits="d_i0", vranges=None):
    """Fig. 7 (3 rows x N times). ``xunits`` selects the horizontal normalization from
    :data:`FIG7_XUNITS` ("d_i0" -> shock_fig7.png, "rho_i0" -> shock_fig7_rho_i0.png);
    only the x scaling/label/filename change, so the panels stay directly comparable.

    ``vranges`` optionally pins the v_z/v_sh axis of each row, as a dict with any of
    ``ambient`` / ``piston`` / ``electron`` -> ``(lo, hi)``. Rows left out keep their
    defaults: the two ion rows share fig_phase's band, and the electron row is sized
    symmetrically from a percentile of |v_z|/v_sh. Pinning them makes panels comparable
    across runs, which the auto-sized electron axis cannot be."""
    xu = FIG7_XUNITS[xunits]
    xf = sc.di0 / getattr(sc, xu["attr"])      # d_i0-valued windows -> requested units
    idx = _select_times(frames, times, sc) if times else _select_nonzero(frames, nframes, sc)
    vnorm = shock.v_sh
    ks = max(5, int(round(0.30 * sc.di0 / sc.dz)))    # smooth overlays over ~0.30 d_i0
    C = kinshock.units.C
    ncol = len(idx)
    wins = {i: _shock_window(frames[i], cfg, sc, shock) for i in idx}

    vr = vranges or {}

    # Ion rows default to fig_phase's shared band; each may be pinned separately.
    v_amb = np.linspace(*vr.get("ambient", (PHASE_V_LO, PHASE_V_HI)), 201)
    v_pis = np.linspace(*vr.get("piston", (PHASE_V_LO, PHASE_V_HI)), 201)

    if "electron" in vr:
        v_e = np.linspace(*vr["electron"], 201)
    else:
        # electron v-axis: symmetric range from a high percentile of |v_z|/v_sh over the
        # selected frames (electrons are hot and not a directed beam), clipped to a band.
        evabs = []
        for i in idx:
            _, uz, _ = _electron_phase_weighted(frames[i], cfg, sc)
            if len(uz):
                evabs.append(np.nanpercentile(np.abs(uz * C / vnorm), FIG7_E_VPCTL))
        ve = float(np.clip(max(evabs) if evabs else FIG7_E_VMIN, FIG7_E_VMIN, FIG7_E_VMAX))
        v_e = np.linspace(-ve, ve, 201)

    # electron-row profile (B_x/B0, n_e/n_e0) scale: anchor on the smoothed magnetic ramp
    # in windows that have cleared the dense piston (same recipe as fig_phase).
    slab_di0 = (cfg["geometry"]["slab_halfwidth_di"] * sc.di) / sc.di0
    piston_margin = max(3.0, 4.0 * slab_di0)

    def _win_bmax(i):
        fr = frames[i]; lo, hi, _ = wins[i]
        zc = np.asarray(fr.z_centers) / sc.di0
        w = (zc >= lo) & (zc <= hi)
        return float(np.nanmax(_smooth((fr.Bx / sc.B0)[w], ks))) if w.any() else 0.0

    clean = [i for i in idx if wins[i][0] > piston_margin]
    pool = clean if clean else list(idx)
    prof_max = float(np.clip(1.3 * max(_win_bmax(i) for i in pool), 4.0, 14.0))

    # Align the overplotted B_x / n_e baseline (value 0) with the phase-space v_z = 0
    # line, matching fig_phase's convention: the twin axis' lower limit is set so that
    # its 0 falls at the same panel fraction as v_z = 0 on the electron velocity axis.
    f0 = (0.0 - v_e[0]) / (v_e[-1] - v_e[0])
    rlo = -f0 / (1.0 - f0) * prof_max

    fig, axes = plt.subplots(3, ncol, figsize=(3.25 * ncol, 8.6), squeeze=False,
                             sharey="row")
    e_species = _electron_species(cfg)

    for c, i in enumerate(idx):
        fr = frames[i]
        lo, hi, zfront = wins[i]
        z_edges = np.linspace(lo * xf, hi * xf, 241)

        # Row 1: ambient ions ; Row 2: piston ions  (v_z / v_sh)
        for r, (sp, v_row) in enumerate((("amb_ions", v_amb), ("piston_ions", v_pis))):
            z, uz, w = io.species_phase_weighted(fr, sp, sc, mass=sc.mi)
            axc = axes[r][c]
            _phase_hist_cmap(axc, z / sc.di0 * xf, uz * C / vnorm, w, z_edges, v_row)
            axc.axhline(1.0, color="0.25", ls=":", lw=0.9)          # v_z = v_sh

        # Row 3: total electron phase space + B_x (black) & n_e (red) overlays
        ax = axes[2][c]
        ze, ue, we = _electron_phase_weighted(fr, cfg, sc)
        _phase_hist_cmap(ax, ze / sc.di0 * xf, ue * C / vnorm, we, z_edges, v_e)
        ax.axhline(0.0, color="0.25", ls=":", lw=0.8)
        zc = np.asarray(fr.z_centers) / sc.di0
        m = (zc >= lo) & (zc <= hi)
        ne = io.species_density(fr, e_species)
        ax2 = ax.twinx()
        # white background -> B_x black and n_e red, exactly as the paper's Fig. 7.
        ax2.plot(zc[m] * xf, _smooth((fr.Bx / sc.B0)[m], ks), color="k", lw=1.3,
                 label=r"$B_x/B_0$")
        ax2.plot(zc[m] * xf, _smooth(np.where(ne > 0, ne / sc.namb, 0.0)[m], ks),
                 color="#d62728", lw=1.3, label=r"$n_e/n_{e0}$", ls='--')
        ax2.set_ylim(rlo, prof_max)
        ax2.set_yticks(np.arange(0, prof_max + 0.1, 2.0))
        if c == ncol - 1:
            ax2.set_ylabel(r"$B_x/B_0,\ \ n_e/n_{e0}$")
            ax2.legend(loc="upper right", fontsize=7, framealpha=0.85,
                       facecolor="white", labelcolor="k", edgecolor="0.7")
        else:
            ax2.set_yticklabels([])

        # shock-front guide on all three rows of the column (region marker)
        for r in range(3):
            if zfront is not None:
                axes[r][c].axvline(zfront * xf, color="0.4", ls="-", lw=0.7, alpha=0.6)
        axes[0][c].set_title(rf"$t\,\omega_{{ci0}}={fr.time*sc.wci0:.2f}$", fontsize=9)
        axes[2][c].set_xlabel(xu["label"])

    axes[0][0].set_ylabel(r"ambient ion  $v_z / v_{sh}$")
    axes[1][0].set_ylabel(r"piston ion  $v_z / v_{sh}$")
    axes[2][0].set_ylabel(r"electron  $v_z / v_{sh}$")
    title = (f"{cfg['meta']['run_id']}: initial shock formation "
             r"(Fig. 7 — ambient-ion / piston-ion / electron phase space)")
    if xu["note"]:
        title += ("\n" + rf"$z$ in {xu['note']}"
                  + rf"$\ = {getattr(sc, xu['attr'])/sc.de:.0f}\,d_e"
                  + rf" = {getattr(sc, xu['attr'])/sc.di0:.1f}\,d_{{i0}}$")
    fig.suptitle(title)
    fig.tight_layout()
    return P.savefig(fig, "shock_fig7" + xu["suffix"] + ".png",
                     run_id=cfg["meta"]["run_id"])


# --- D: reflected-ion fraction + timescales ---
def fig_reflected(frames, cfg, sc, shock):
    vnorm = shock.v_sh
    t, G = [], []
    for fr in frames:
        z, uz = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        G.append(metrics.reflected_fraction_G(uz * kinshock.units.C, vnorm))
        t.append(fr.time)
    t, G = np.array(t), np.array(G)
    tstar, i1 = metrics.onset_time_from_G(t, G)

    dG = np.gradient(G)

    fig, ((a_g, a_f), (a_dg, a_df)) = plt.subplots(2, 2, figsize=(12.4, 6.2))
    a_g.plot(t * sc.wci0, G, color=P.C_AMBIENT, ms=3, lw=1.0)
    if i1 >= 0:
        a_g.axvline(tstar * sc.wci0, color=P.C_PISTON, ls="--",
                    label=rf"$t^*_1={tstar*sc.wci0:.2f}\,\omega_{{ci0}}^{{-1}}$ (target ~1)")
        a_g.legend(frameon=False, fontsize=9)
    a_g.set_xlabel(r"$t\,\omega_{ci0}$")
    a_g.set_ylabel(r"$G=N_{a,refl}/N_{a,tot}$")
    P.style_axes(a_g)

    a_dg.plot(t * sc.wci0, dG, color=P.C_AMBIENT, ms=3, lw=1.0)
    if i1 >= 0:
        a_dg.axvline(tstar * sc.wci0, color=P.C_PISTON, ls="--",
                    label=rf"$t^*_1={tstar*sc.wci0:.2f}\,\omega_{{ci0}}^{{-1}}$ (target ~1)")
        a_dg.legend(frameon=False, fontsize=9)
    a_dg.set_xlabel(r"$t\,\omega_{ci0}$")
    a_dg.set_ylabel(r"$dG/dt$")
    P.style_axes(a_dg)

    if i1 >= 0:
        fr = frames[i1]
        z, uz = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        edges = np.asarray(fr.z_edges)
        F, centers = metrics.reflected_profile_F(z * 1.0, uz * kinshock.units.C, vnorm, edges)
        for _ in range(6):
            F = _smooth(F, k=50)
        zc = centers / sc.rho_i0
        zstar, _ = metrics.onset_location_from_F(zc, F)
        dF = np.gradient(F)
        a_f.plot(zc, F, color=P.C_AMBIENT, lw=1.0)
        a_f.axvline(zstar, color=P.C_PISTON, ls="--",
                    label=rf"$z^*_1={zstar:.1f}\,\rho_{{i0}}$ (target ~1)")
        a_f.legend(frameon=False, fontsize=9)
        a_f.set_xlim(0, max(2.0, np.nanmax(zc[F > 0]) if (F > 0).any() else 2.0))
        a_df.plot(zc, dF, color=P.C_AMBIENT, lw=1.0)
        a_df.axvline(zstar, color=P.C_PISTON, ls="--",
                    label=rf"$z^*_1={zstar:.1f}\,\rho_{{i0}}$ (target ~1)")
        a_df.legend(frameon=False, fontsize=9)
        a_df.set_xlim(0, max(2.0, np.nanmax(zc[F > 0]) if (F > 0).any() else 2.0))
    a_f.set_xlabel(r"$z / \rho_{i0}$")
    a_f.set_ylabel(r"$F(z,t^*_1)$")
    a_df.set_xlabel(r"$z / \rho_{i0}$")
    a_df.set_ylabel(r"$dF/dt$")
    P.style_axes(a_f)
    P.stamp(a_g, cfg, sc)
    P.style_axes(a_df)
    P.stamp(a_dg, cfg, sc)
    fig.suptitle(f"{cfg['meta']['run_id']}: reflected ambient ions & onset")
    fig.tight_layout()
    P.savefig(fig, "shock_reflected.png", run_id=cfg["meta"]["run_id"])
    return tstar


# --- C: seven criteria ---
def criteria_table(frames, cfg, sc, shock):
    vnorm = shock.v_sh
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di
    out = []
    for fr in frames:
        n = _ion_density(fr, cfg)
        z_a, uz_a = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        front_z = shock.front_m(fr)
        # piston peak-field position (within the slab)
        zc = np.asarray(fr.z_centers)
        pmask = np.abs(zc) <= zexcl
        piston_z = zc[pmask][np.argmax(fr.Bperp[pmask])] if pmask.any() else 0.0
        if not np.isfinite(front_z):
            continue
        res = metrics.evaluate_criteria(
            zc=fr.z_centers, n_e=n, Bmag=fr.Bperp,
            uz_ambient=uz_a * kinshock.units.C, z_ambient=z_a,
            scales=sc, vsh=vnorm, v_front=vnorm,
            piston_field_z=piston_z, front_z=front_z)
        out.append({"t_wci0": fr.time * sc.wci0, "flags": res.flags,
                    "values": res.values, "is_precursor": res.is_precursor,
                    "is_shock": res.is_shock})
    path = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "criteria.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=_json_native)
    print("wrote", path)
    prec = next((r["t_wci0"] for r in out if r["is_precursor"]), None)
    shock = next((r["t_wci0"] for r in out if r["is_shock"]), None)
    print(f"  first precursor (crit 1-6) at t*wci0 = {prec}; "
          f"first shock (crit 1-7) at t*wci0 = {shock} (paper t*_2 ~ 2.5)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir", nargs="?", default=os.path.join(ROOT, "runs", "R1"))
    ap.add_argument("--nframes", type=int, default=5)
    ap.add_argument("--phase-times", type=float, nargs="+", default=None, metavar="T",
                    help="times (in t*wci0) for the ion phase-space panels; snaps to the "
                         "nearest available frame. Default: --nframes evenly-spaced frames.")
    ap.add_argument("--fig7-xunits", nargs="+", choices=sorted(FIG7_XUNITS), metavar="U",
                    default=["d_i0"],
                    help="horizontal-axis normalization(s) for Fig. 7: d_i0 (default, "
                         "-> shock_fig7.png) and/or rho_i0, the upstream ion gyroradius "
                         "v_p/omega_ci0 (-> shock_fig7_rho_i0.png).")
    ap.add_argument("--v-ambient", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                    help="fig7 ambient-ion row v_z/v_sh limits (default: the fig_phase "
                         "band %g %g)" % (PHASE_V_LO, PHASE_V_HI))
    ap.add_argument("--v-piston", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                    help="fig7 piston-ion row v_z/v_sh limits (default: same band)")
    ap.add_argument("--v-electron", nargs=2, type=float, default=None, metavar=("LO", "HI"),
                    help="fig7 electron row v_z/v_sh limits. Default is auto-sized and "
                         "SYMMETRIC, from a percentile of |v_z|/v_sh, so it varies between "
                         "runs; pin it to compare runs panel-for-panel.")
    ap.add_argument("--only", nargs="+", default=None, metavar="FIG",
                    choices=["streak", "trajectory", "lineouts", "phase", "fig7",
                             "reflected", "criteria"],
                    help="build only these figures (default: all). e.g. --only fig7.")
    args = ap.parse_args()
    want = (lambda name: args.only is None or name in args.only)

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    print(f"{cfg['meta']['run_id']}: d_e={sc.de:.3e} m  d_i0={sc.di0:.3e} m  "
          f"B0={sc.B0:.3e} T  rho_i0={sc.rho_i0/sc.de:.0f} d_e")

    # particle frames (diag1): phase space, density, criteria, trajectory
    particle_pfs = io.plotfiles(args.run_dir)
    frames = load_series(particle_pfs)
    print(f"{len(frames)} particle frames, t*wci0 in [{frames[0].time*sc.wci0:.2f}, "
          f"{frames[-1].time*sc.wci0:.2f}]")

    # field frames (diag_fields if present, else the same diag1 frames): streak
    if want("streak"):
        field_pfs = io.field_plotfiles(args.run_dir)
        if field_pfs == particle_pfs:
            field_frames = frames
        else:
            field_frames = load_series(field_pfs)
            print(f"{len(field_frames)} field-only frames for the streak "
                  f"({len(field_frames)/max(len(frames),1):.0f}x the particle cadence)")

    # SINGLE SOURCE OF TRUTH for v_sh + z_front(t): the by-eye fit (shock_fit.yaml via
    # scripts/tune_shock.py) if present, else the automatic track_front fallback.
    shock = resolve_shock(args.run_dir, frames, cfg, sc)

    if want("streak"):
        fig_streak(field_frames, cfg, sc)
    if want("trajectory"):
        fig_trajectory(frames, cfg, sc, shock)
    if want("lineouts"):
        fig_lineouts(frames, cfg, sc, args.nframes)
    if want("phase"):
        fig_phase(frames, cfg, sc, shock, args.nframes, times=args.phase_times)
    if want("fig7"):
        vranges = {k: tuple(v) for k, v in (("ambient", args.v_ambient),
                                            ("piston", args.v_piston),
                                            ("electron", args.v_electron)) if v}
        for xu in args.fig7_xunits:
            fig_fig7(frames, cfg, sc, shock, args.nframes, times=args.phase_times,
                     xunits=xu, vranges=vranges)
    if want("reflected"):
        fig_reflected(frames, cfg, sc, shock)
    if want("criteria"):
        criteria_table(frames, cfg, sc, shock)


if __name__ == "__main__":
    main()
