#!/usr/bin/env python3
"""scripts/tune_shock.py -- fit a run's shock front BY EYE (KinShock2020).

The analysis suite (make_figures, crosscheck) needs one consistent shock speed and
front trajectory. Deriving them automatically (track_front + speed_from_trajectory)
drifted out of sync between scripts. This tool instead lets you fit the front by eye
against the B_perp / n_e streaks and writes the result to

    runs/<ID>/shock_fit.yaml          # SINGLE SOURCE OF TRUTH (kinshock.metrics.ShockFit)
        shock.v_sh_over_c, shock.z0_de          -> z_front(t) = z0 + v_sh*t
        shock.fronts_de: {t*wci0: z_front}      -> optional per-time overrides (regions mode)

which every downstream diagnostic then reads. Display is a PNG you refresh in your
editor (robust over SSH, no X11): each command re-renders media/<ID>/tune_*.png.

Two modes
---------
trajectory (default) -- fit shock.v_sh_over_c / shock.z0_de against the |B_perp| and
    n_e streaks. Both streaks run at the high-cadence field-only diagnostic's rate
    (diags/diag_fields*, ~20x diag1), with n_e taken from the deposited rho_<species>
    fields since that diagnostic writes no particles. The trial line z = z0 + v_sh*t is
    overlaid; adjust until it rides the forward front. Commands:
        v <val>     set trial v_sh [c]            (e.g. v 0.14)
        z <val>     set trial z0   [d_e]          (front intercept at t=0)
        ma <val>    set v_sh from a target M_A    (v_sh = val * v_A)
        save        write shock.{v_sh_over_c,z0_de} to shock_fit.yaml (asks y/N)
        q           quit

regions -- refine the front position at ONE time against its ion phase space + line-outs.
    Pick the time with --time (t*wci0). Commands:
        shock <z>   set trial front position [d_i0]
        save        write shock.fronts_de[t*wci0] to shock_fit.yaml (asks y/N)
        q           quit

Examples:
    python scripts/tune_shock.py runs/R1_warm
    python scripts/tune_shock.py runs/R1_warm --mode regions --time 2.4
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import io, metrics, plotting as P  # noqa: E402
from kinshock.units import C  # noqa: E402


def _electron_species(cfg):
    return [n for n, s in cfg["species"].items() if s.get("kind") == "electron"]


def _block_mean(a, k):
    """Block-average a 1D profile over k cells (display decimation along z)."""
    a = np.asarray(a, dtype=float)
    if k <= 1:
        return a
    n = (len(a) // k) * k                      # drops a <k-cell remainder at the top
    return a[:n].reshape(-1, k).mean(axis=1)


def _confirm_write(run_dir, sc, v_sh, z0, fronts, no_write):
    """Preview the shock_fit.yaml edits and (unless --no-write) ask y/N before saving."""
    print(f"  -> shock_fit.yaml: v_sh_over_c={v_sh/C:.6f}  z0_de={z0/sc.de:.3f}  "
          f"fronts_de={{{', '.join(f'{k:.2f}:{v/sc.de:.1f}' for k, v in sorted(fronts.items()))}}}")
    if no_write:
        print("  (--no-write) not saved.")
        return
    try:
        ans = input("  write shock_fit.yaml? [y/N] ").strip().lower()
    except EOFError:
        ans = "n"
    if ans == "y":
        path = metrics.save_shock_fit(run_dir, sc, v_sh, z0, fronts)
        print(f"  wrote {path}")
    else:
        print("  not saved.")


# --------------------------------------------------------------------------- #
# Trajectory mode
# --------------------------------------------------------------------------- #
class TrajectoryTuner:
    """|B_perp| + n_e streaks with a movable trial front line z = z0 + v_sh*t."""

    Z_DISPLAY = 1600      # target z bins for the streak (see the decimation note below)

    def __init__(self, run_dir, cfg, sc, args):
        self.run_dir, self.cfg, self.sc = run_dir, cfg, sc
        self.png = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "tune_trajectory.png")
        esp = _electron_species(cfg)
        rho_fields = io.rho_field_names(esp)
        charges = {n: s.get("charge_state", 1) for n, s in cfg["species"].items()}

        # Prefer the dense field-only series (diag_fields, ~20x the particle cadence).
        # It has write_species = 0, so n_e must come from the deposited rho_<species>
        # fields rather than macroparticles -- that is exactly what makes the full field
        # cadence usable for BOTH panels. Runs whose plotfiles carry no per-species rho
        # (e.g. diag1 alone) fall back to the particle diags + macroparticle histograms.
        particle_pfs = io.plotfiles(run_dir)
        pfs = io.field_plotfiles(run_dir)
        use_rho = bool(io.load_frame(pfs[0], fields=rho_fields).comps)
        if not use_rho and pfs != particle_pfs:
            print("  field plotfiles carry no per-species rho -> using the particle diags")
            pfs = particle_pfs
        if len(pfs) < 2:
            raise RuntimeError(f"need >=2 plotfiles for a streak; got {len(pfs)}")
        stride = max(1, args.stride)
        pfs = pfs[::stride]
        src = "deposited rho fields" if use_rho else "macroparticle histograms"
        print(f"loading {len(pfs)} frames for the streak (stride {stride}; n_e from {src}; "
              f"{len(particle_pfs)} particle frames exist)...")

        # Stream the series: hold one yt dataset at a time and keep only the decimated
        # profiles, so a 1000-frame x 25000-cell streak stays in a few hundred MB.
        zbin, t_s, B, ne = None, [], [], []
        for i, p in enumerate(pfs):
            fr = io.load_frame(p, fields=rho_fields if use_rho else ())
            if zbin is None:
                nz = len(fr.z_centers)
                zbin = args.zbin if args.zbin > 0 else max(1, nz // self.Z_DISPLAY)
                dz = float(fr.z_centers[1] - fr.z_centers[0])
                self.zc_di0 = _block_mean(fr.z_centers, zbin) / sc.di0
                if zbin > 1:
                    # The streak is drawn into ~1200 px, so plotting all nz cells costs
                    # pcolormesh dearly and resolves nothing extra; block-averaging also
                    # suppresses the far-upstream numerical hash (lambda < 2-3 d_e).
                    print(f"  z: {nz} cells -> {len(self.zc_di0)} bins "
                          f"(block-mean x{zbin} = {zbin*dz/sc.de:.1f} d_e)")
            t_s.append(fr.time)
            B.append(_block_mean(fr.Bperp, zbin) / sc.B0)
            d = (io.field_species_density(fr, esp, charges) if use_rho
                 else io.species_density(fr, esp))
            ne.append(_block_mean(d, zbin) / sc.namb)
            del fr                                  # release the yt dataset
            if (i + 1) % 50 == 0 or i + 1 == len(pfs):
                print(f"  {i+1}/{len(pfs)} frames", end="\r", flush=True)
        print()
        self.t_s = np.array(t_s)
        self.t_wci0 = self.t_s * sc.wci0
        self.B = np.array(B, dtype=np.float32)                                 # (nt, nz)
        self.ne = np.array(ne, dtype=np.float32)
        # seed from an existing fit, else the model
        fit = metrics.load_shock_fit(run_dir, sc)
        if fit is not None:
            self.v_sh, self.z0 = fit.v_sh, fit.z0
            self.fronts = fit.fronts
            print(f"seeded from shock_fit.yaml: v_sh={self.v_sh/C:.4f}c z0={self.z0/sc.de:.1f} d_e")
        else:
            self.v_sh, self.z0, self.fronts = sc.vsh_model, 0.0, {}
            print(f"no shock_fit.yaml; seeded from model v_sh={self.v_sh/C:.4f}c")

    def render(self):
        MA, Mms = self.v_sh / self.sc.vA, self.v_sh / np.hypot(self.sc.vA, self.sc.Cs0)
        z_trial = (self.z0 + self.v_sh * self.t_s) / self.sc.di0
        fig, (aB, aN) = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
        for ax, S, lbl, cmap, vmax in (
            (aB, self.B, r"$B_\perp/B_0$", "viridis", max(2.0, np.nanpercentile(self.B, 99))),
            # cap n_e well below the dense piston (~250 n_e0) so the ~2-5x shock
            # compression -- the feature to fit the front against -- stays visible.
            (aN, np.where(self.ne > 0, self.ne, np.nan), r"$n_e/n_{e0}$", "magma", 8.0)):
            pc = ax.pcolormesh(self.t_wci0, self.zc_di0, S.T, shading="auto",
                               cmap=cmap, vmin=0, vmax=vmax)
            ax.plot(self.t_wci0, z_trial, "-", color="w", lw=2.0,
                    label=rf"trial: $v_{{sh}}$={self.v_sh/C:.4f}c, $z_0$={self.z0/self.sc.de:.0f}$d_e$, $M_A$={MA:.2f}")
            for k, zf in self.fronts.items():           # per-time overrides (regions mode)
                ax.plot(k, zf / self.sc.di0, "x", color="cyan", ms=9, mew=2.0)
            ax.set_ylabel(r"$z / d_{i0}$")
            ax.set_ylim(0, self.zc_di0.max())
            ax.legend(frameon=False, labelcolor="w", loc="upper left", fontsize=9)
            fig.colorbar(pc, ax=ax, label=lbl)
        aN.set_xlabel(r"$t\,\omega_{ci0}$")
        fig.suptitle(f"{self.cfg['meta']['run_id']}: trajectory fit  "
                     f"(v_sh={self.v_sh/C:.4f}c, M_A={MA:.2f}, M_ms={Mms:.2f})")
        fig.tight_layout()
        fig.savefig(self.png, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  v_sh={self.v_sh/C:.4f}c  M_A={MA:.2f}  M_ms={Mms:.2f}  "
              f"v_sh/v_p={self.v_sh/self.sc.vp_model:.2f}")
        print(f"  ↻ wrote {self.png} -- refresh in your IDE")

    def loop(self, no_write):
        print("\ntrajectory mode -- commands: v <c> | z <d_e> | ma <M_A> | save | q")
        self.render()
        while True:
            try:
                raw = input("tune> ").strip()
            except EOFError:
                break
            if not raw:
                continue
            cmd, *rest = raw.split()
            cmd = cmd.lower()
            if cmd in ("q", "quit", "exit"):
                break
            elif cmd == "v" and rest:
                self.v_sh = float(rest[0]) * C; self.render()
            elif cmd == "z" and rest:
                self.z0 = float(rest[0]) * self.sc.de; self.render()
            elif cmd == "ma" and rest:
                self.v_sh = float(rest[0]) * self.sc.vA; self.render()
            elif cmd == "save":
                _confirm_write(self.run_dir, self.sc, self.v_sh, self.z0, self.fronts, no_write)
            else:
                print("  ? commands: v <c> | z <d_e> | ma <M_A> | save | q")


# --------------------------------------------------------------------------- #
# Regions mode
# --------------------------------------------------------------------------- #
class RegionsTuner:
    """One frame's ion phase space + n_e/B line-outs with a movable front marker."""

    def __init__(self, run_dir, cfg, sc, args):
        self.run_dir, self.cfg, self.sc = run_dir, cfg, sc
        pfs = io.plotfiles(run_dir)
        frames_t = [(p, io.load_frame(p)) for p in pfs]
        tw = np.array([fr.time * sc.wci0 for _, fr in frames_t])
        i = int(np.argmin(np.abs(tw - args.time)))
        self.p, self.fr = frames_t[i]
        self.t_wci0 = float(tw[i])
        print(f"regions: requested t*wci0={args.time:.2f} -> frame t*wci0={self.t_wci0:.2f}")
        self.png = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]),
                                f"tune_regions_t{self.t_wci0:.2f}.png")
        # seed / require an existing trajectory fit (v_sh, z0) so save is well-defined
        self.fit = metrics.load_shock_fit(run_dir, sc)
        if self.fit is None:
            raise RuntimeError("run trajectory mode first (need v_sh/z0 in shock_fit.yaml)")
        self.fronts = dict(self.fit.fronts)
        self.z_front = self.fronts.get(self.t_wci0,
                                       self.fit.z_front(self.fr.time, self.t_wci0))

    def render(self):
        sc, fr = self.sc, self.fr
        z, uz, w = io.species_phase_weighted(fr, "amb_ions", sc, mass=sc.mi)
        esp = _electron_species(self.cfg)
        ne = io.species_density(fr, esp)
        zc = np.asarray(fr.z_centers) / sc.di0
        fig, (ap, al) = plt.subplots(2, 1, figsize=(11, 8), sharex=True,
                                     gridspec_kw=dict(height_ratios=[2, 1]))
        if len(z):
            ap.hist2d(z / sc.di0, uz * C / self.fit.v_sh, bins=[220, 160],
                      range=[[0, zc.max()], [-1, 3]], cmap="cividis", cmin=1)
        ap.axhline(1.0, color="0.4", ls=":", lw=1.0)
        ap.set_ylabel(r"ambient ion $v_z/v_{sh}$")
        al.plot(zc, np.where(ne > 0, ne / sc.namb, np.nan), color=P.C_AMBIENT, lw=1.1, label=r"$n_e/n_{e0}$")
        al.plot(zc, fr.Bperp / sc.B0, color=P.C_PISTON, lw=1.1, label=r"$B_\perp/B_0$")
        al.axhline(2.0, color="0.6", ls="--", lw=0.7)     # criterion-3/4 (2x) guide
        al.set_ylim(0, 10)                                 # cap: hide the dense piston, show the shock jump
        al.set_ylabel(r"$n_e/n_{e0},\ B_\perp/B_0$"); al.set_xlabel(r"$z / d_{i0}$")
        al.legend(frameon=False, fontsize=9)
        for ax in (ap, al):
            ax.axvline(self.z_front / sc.di0, color="cyan", ls="-", lw=1.8,
                       label=f"front {self.z_front/sc.di0:.2f} d_i0")
        ap.legend(frameon=False, fontsize=9, loc="upper right")
        fig.suptitle(f"{self.cfg['meta']['run_id']}: front @ t*wci0={self.t_wci0:.2f}  "
                     f"(z_front={self.z_front/sc.di0:.2f} d_i0 = {self.z_front/sc.rho_i0:.2f} rho_i0)")
        fig.tight_layout()
        fig.savefig(self.png, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"  z_front={self.z_front/sc.di0:.2f} d_i0 = {self.z_front/sc.rho_i0:.2f} rho_i0")
        print(f"  ↻ wrote {self.png} -- refresh in your IDE")

    def loop(self, no_write):
        print(f"\nregions mode t*wci0={self.t_wci0:.2f} -- commands: shock <d_i0> | save | q")
        self.render()
        while True:
            try:
                raw = input("tune> ").strip()
            except EOFError:
                break
            if not raw:
                continue
            cmd, *rest = raw.split()
            cmd = cmd.lower()
            if cmd in ("q", "quit", "exit"):
                break
            elif cmd == "shock" and rest:
                self.z_front = float(rest[0]) * self.sc.di0; self.render()
            elif cmd == "save":
                self.fronts[self.t_wci0] = self.z_front
                _confirm_write(self.run_dir, self.sc, self.fit.v_sh, self.fit.z0,
                               self.fronts, no_write)
            else:
                print("  ? commands: shock <d_i0> | save | q")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--mode", choices=("trajectory", "regions"), default="trajectory")
    ap.add_argument("--time", type=float, default=None,
                    help="(regions) frame time t*wci0 to refine; default: last frame.")
    ap.add_argument("--stride", type=int, default=1,
                    help="(trajectory) plotfile stride for the streak (default 1 = every "
                         "field frame; raise it to cut load time on long runs).")
    ap.add_argument("--zbin", type=int, default=0,
                    help="(trajectory) block-average the streak over this many z cells "
                         "for display (default 0 = auto, ~1600 bins).")
    ap.add_argument("--no-write", action="store_true",
                    help="preview edits but never modify shock_fit.yaml.")
    args = ap.parse_args()

    cfg = kinshock.load(args.run_dir)
    sc = kinshock.units.derive(cfg)
    print(f"{cfg['meta']['run_id']}: v_A={sc.vA/C:.4f}c  Cs_ab={sc.Cs_ab/C:.4f}c  "
          f"d_e={sc.de:.3e}m  d_i0={sc.di0:.3e}m")

    if args.mode == "trajectory":
        TrajectoryTuner(args.run_dir, cfg, sc, args).loop(args.no_write)
    else:
        if args.time is None:
            args.time = float(io.load_frame(io.plotfiles(args.run_dir)[-1]).time * sc.wci0)
            print(f"(no --time; using last frame t*wci0={args.time:.2f})")
        RegionsTuner(args.run_dir, cfg, sc, args).loop(args.no_write)


if __name__ == "__main__":
    main()
