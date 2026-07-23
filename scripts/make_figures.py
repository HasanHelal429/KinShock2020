#!/usr/bin/env python3
"""Reproduce the Schaeffer 2020 shock diagnostics for a KinShock2020 run.

Config-driven (all scales from runs/RXX/config.yaml via kinshock.units); writes
to media/<run_id>/. Implements analyses A–E of REPLICATION_PLAN.md §6:

  A. B_perp(z,t) streak plot with piston/shock speed lines  -> shock_streak.png
     (uses the high-cadence field-only diag_fields series when present, else diag1)
  A. shock-front trajectory + measured v_sh, M_A, M_ms       -> shock_trajectory.png
  B. ion (z, u_z) phase-space distribution at several times  -> shock_phase.png
     + density / B_perp line-outs                            -> shock_lineouts.png
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


# --- A: streak + trajectory ---
def fig_streak(frames, cfg, sc):
    t = np.array([fr.time for fr in frames]) * sc.wci0
    zc = np.asarray(frames[0].z_centers) / sc.di0
    B = np.array([fr.Bperp / sc.B0 for fr in frames])  # (nt, nz)
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


def fig_trajectory(frames, cfg, sc):
    ts, zs = [], []
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di   # exclude piston slab
    for fr in frames:
        n = _ion_density(fr, cfg)
        zf = metrics.track_front(fr.z_centers, n, sc.namb, threshold=1.5, z_exclude=zexcl)
        if np.isfinite(zf):
            ts.append(fr.time)
            zs.append(zf)
    if len(ts) < 2:
        print("  shock front not detected; skipping trajectory")
        return None, np.nan
    ts, zs = np.array(ts), np.array(zs)
    vsh = metrics.speed_from_trajectory(ts, zs)
    MA, Mms = vsh / sc.vA, vsh / np.sqrt(sc.vA ** 2 + sc.Cs0 ** 2)
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    ax.plot(ts * sc.wci0, zs / sc.di0, "o-", color=P.C_PISTON, ms=3, lw=1.0)
    ax.set_xlabel(r"$t\,\omega_{ci0}$")
    ax.set_ylabel(r"shock front $z / d_{i0}$")
    P.style_axes(ax)
    txt = (rf"$v_{{sh}}={vsh/kinshock.units.C:.4f}\,c$" + "\n"
           rf"$M_A={MA:.1f}$ (target {cfg['targets']['M_A']})" + "\n"
           rf"$M_{{ms}}={Mms:.1f}$ (target {cfg['targets']['M_ms']})" + "\n"
           rf"$v_{{sh}}/v_p={vsh/sc.vp_model:.2f}$ (RH $\to$ 4/3)")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec=P.C_REF, alpha=0.9))
    fig.suptitle(f"{cfg['meta']['run_id']}: shock-front trajectory")
    fig.tight_layout()
    P.savefig(fig, "shock_trajectory.png", run_id=cfg["meta"]["run_id"])
    print(f"  measured v_sh = {vsh/kinshock.units.C:.4f} c ; M_A = {MA:.2f} ; M_ms = {Mms:.2f}")
    return vsh, MA


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


def _shock_window(fr, cfg, sc):
    """(z_lo, z_hi, z_front) in d_i0: a window that FOLLOWS the forward (+z) shock front,
    sized to the feature (extents grow with the front position; rho_i0 floors) and clipped
    to z >= 0. Falls back to the near-piston region before a front is detectable."""
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di
    zf = metrics.track_front(fr.z_centers, _ion_density(fr, cfg), sc.namb,
                             threshold=1.5, z_exclude=zexcl, side=+1)
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


def fig_phase(frames, cfg, sc, vsh, nframes):
    idx = _select_nonzero(frames, nframes, sc)
    fig, axes = plt.subplots(1, len(idx), figsize=(3.4 * len(idx), 4.0), sharey=True)
    axes = np.atleast_1d(axes)
    vnorm = vsh if (vsh and np.isfinite(vsh)) else sc.vsh_model
    v_edges = np.linspace(PHASE_V_LO, PHASE_V_HI, 201)
    e_species = _electron_species(cfg)
    ks = max(5, int(round(0.30 * sc.di0 / sc.dz)))   # smooth overlays over ~0.30 d_i0

    wins = {i: _shock_window(frames[i], cfg, sc) for i in idx}

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


# --- D: reflected-ion fraction + timescales ---
def fig_reflected(frames, cfg, sc, vsh):
    vnorm = vsh if (vsh and np.isfinite(vsh)) else sc.vsh_model
    t, G = [], []
    for fr in frames:
        z, uz = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        G.append(metrics.reflected_fraction_G(uz * kinshock.units.C, vnorm))
        t.append(fr.time)
    t, G = np.array(t), np.array(G)
    tstar, i1 = metrics.onset_time_from_G(t, G)

    fig, (a_g, a_f) = plt.subplots(2, 1, figsize=(7.2, 6.2))
    a_g.plot(t * sc.wci0, G, "o-", color=P.C_AMBIENT, ms=3, lw=1.0)
    if i1 >= 0:
        a_g.axvline(tstar * sc.wci0, color=P.C_PISTON, ls="--",
                    label=rf"$t^*_1={tstar*sc.wci0:.2f}\,\omega_{{ci0}}^{{-1}}$ (target ~1)")
        a_g.legend(frameon=False, fontsize=9)
    a_g.set_xlabel(r"$t\,\omega_{ci0}$")
    a_g.set_ylabel(r"$G=N_{a,refl}/N_{a,tot}$")
    P.style_axes(a_g)

    if i1 >= 0:
        fr = frames[i1]
        z, uz = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        edges = np.asarray(fr.z_edges)
        F, centers = metrics.reflected_profile_F(z * 1.0, uz * kinshock.units.C, vnorm, edges)
        zc = centers / sc.rho_i0
        zstar, _ = metrics.onset_location_from_F(zc, F)
        a_f.plot(zc, F, color=P.C_AMBIENT, lw=1.0)
        a_f.axvline(zstar, color=P.C_PISTON, ls="--",
                    label=rf"$z^*_1={zstar:.1f}\,\rho_{{i0}}$ (target ~1)")
        a_f.legend(frameon=False, fontsize=9)
        a_f.set_xlim(0, max(6.0, np.nanmax(zc[F > 0]) if (F > 0).any() else 6.0))
    a_f.set_xlabel(r"$z / \rho_{i0}$")
    a_f.set_ylabel(r"$F(z,t^*_1)$")
    P.style_axes(a_f)
    P.stamp(a_g, cfg, sc)
    fig.suptitle(f"{cfg['meta']['run_id']}: reflected ambient ions & onset")
    fig.tight_layout()
    P.savefig(fig, "shock_reflected.png", run_id=cfg["meta"]["run_id"])
    return tstar


# --- C: seven criteria ---
def criteria_table(frames, cfg, sc, vsh):
    vnorm = vsh if (vsh and np.isfinite(vsh)) else sc.vsh_model
    zexcl = cfg["geometry"]["slab_halfwidth_di"] * sc.di
    out = []
    for fr in frames:
        n = _ion_density(fr, cfg)
        z_a, uz_a = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        front_z = metrics.track_front(fr.z_centers, n, sc.namb, 1.5, z_exclude=zexcl)
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
    args = ap.parse_args()

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
    field_pfs = io.field_plotfiles(args.run_dir)
    if field_pfs == particle_pfs:
        field_frames = frames
    else:
        field_frames = load_series(field_pfs)
        print(f"{len(field_frames)} field-only frames for the streak "
              f"({len(field_frames)/max(len(frames),1):.0f}x the particle cadence)")

    fig_streak(field_frames, cfg, sc)
    vsh, _ = fig_trajectory(frames, cfg, sc)
    fig_lineouts(frames, cfg, sc, args.nframes)
    fig_phase(frames, cfg, sc, vsh, args.nframes)
    fig_reflected(frames, cfg, sc, vsh)
    criteria_table(frames, cfg, sc, vsh)


if __name__ == "__main__":
    main()
