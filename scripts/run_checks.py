#!/usr/bin/env python3
"""Bring-up / progress checks for a KinShock2020 run (writes to media/testing/).

Always emits a config-summary figure (derived scales vs Table I targets) so
progress is visible even before any WarpX output exists. If the run has produced
plotfiles / reduced diagnostics, also emits:
  * loaded-state sanity: initial per-species density and B_perp profiles, with the
    measured C_s,ab / v_A / M_A / M_ms / beta_0 stamped vs their Table I targets;
  * operator balance: energy-conservation and piston-inventory histories (heater
    <-> injector balance), inherited from the flatfoil validation.

Usage:  python scripts/run_checks.py [runs/R1_phase/R1]
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import kinshock  # noqa: E402
from kinshock import plotting as P  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def config_summary_figure(cfg, sc, warns):
    """Text figure of derived scales vs targets -> media/testing (no run data needed)."""
    rid = cfg.get("meta", {}).get("run_id", "?")
    tgt = cfg.get("targets", {})
    rows = [
        ("M_A", sc.MA, tgt.get("M_A")),
        ("M_ms", sc.Mms, tgt.get("M_ms")),
        ("beta_ab", sc.beta_ab, tgt.get("beta_ab")),
        ("beta_0", sc.beta_0, tgt.get("beta_0")),
        ("C_s,ab / c", sc.Cs_ab / kinshock.units.C, None),
        ("v_A / c", sc.vA / kinshock.units.C, None),
        ("v_sh(model) / c", sc.vsh_model / kinshock.units.C, None),
        ("B0 [T]", sc.B0, None),
        ("d_e [m]", sc.de, None),
        ("d_i0 [m]", sc.di0, None),
        ("rho_i0 / d_e", sc.rho_i0 / sc.de, None),
        ("1/wci0 [wpe^-1]", sc.wci0_inv * sc.wpe, None),
        ("dt*wpe", sc.dt_wpe, None),
        ("n_cell", sc.n_cell, None),
        ("domain half [d_e]", sc.domain_halfwidth / sc.de, None),
        ("steps / wci0^-1", sc.steps_per_wci0, None),
    ]
    lines = [f"{'quantity':22s} {'derived':>14s}  {'Table I':>10s}"]
    lines.append("-" * 50)
    for name, got, want in rows:
        w = "" if want is None else f"{want:>10.4g}"
        lines.append(f"{name:22s} {got:>14.5g}  {w}")
    lines.append("")
    lines.append("validation: " + ("OK (within tolerance)" if not warns else "WARNINGS"))
    for wmsg in warns:
        lines.append("  ! " + wmsg)

    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.axis("off")
    ax.text(0.0, 1.0, f"{rid} — derived scales vs Schaeffer 2020 Table I\n\n" + "\n".join(lines),
            transform=ax.transAxes, va="top", ha="left", family="monospace", fontsize=9)
    return P.savefig(fig, f"{cfg['meta']['run_id']}_config_summary.png", testing=True)


def loaded_state_figure(cfg, sc, frame):
    """Initial per-species density + B_perp from the first plotfile."""
    from kinshock import io
    zc_di0 = np.asarray(frame.z_centers) / sc.di0
    fig, (a_n, a_b) = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True)
    for sp, col in (("piston_ions", P.C_PISTON), ("amb_ions", P.C_AMBIENT)):
        n = io.species_density(frame, sp)
        a_n.plot(zc_di0, np.where(n > 0, n / sc.namb, np.nan), color=col, lw=1.1, label=sp)
    a_n.set_yscale("log")
    a_n.set_ylabel(r"$n_i / n_{e0}$")
    a_n.legend(frameon=False, fontsize=8)
    a_b.plot(zc_di0, frame.Bperp / sc.B0 if sc.B0 else frame.Bperp, color=P.C_REF, lw=1.1)
    a_b.set_ylabel(r"$B_\perp / B_0$")
    a_b.set_xlabel(r"$z / d_{i0}$")
    for ax in (a_n, a_b):
        P.style_axes(ax)
    P.stamp(a_n, cfg, sc, extra=f"t·wpe={frame.time*sc.wpe:.1f} (loaded state)")
    fig.suptitle(f"{cfg['meta']['run_id']}: loaded-state sanity")
    fig.tight_layout()
    return P.savefig(fig, f"{cfg['meta']['run_id']}_loaded_state.png", testing=True)


def operator_balance_figure(cfg, sc, run_dir):
    """Energy-conservation + piston-inventory histories from reduced diags EP/PN."""
    from kinshock import io
    try:
        heP, EP = io.reduced_diag(run_dir, "EP")
        hpn, PN = io.reduced_diag(run_dir, "PN")
    except FileNotFoundError as e:
        print("  (skipping operator-balance figure:", e, ")")
        return None
    # EP/PN columns are [0]step [1]time(s) [2]total ... -- reduced_diag keeps the step
    # column, so time is 1 and the total is 2. Indexing these as 0/1 plotted step-vs-time
    # in BOTH panels and made the conservation check vacuous (found 2026-08-01).
    t = EP[:, 1] * sc.wpe
    fig, (a_e, a_n) = plt.subplots(2, 1, figsize=(7.6, 6.0), sharex=True)
    a_e.plot(t, EP[:, 2] / EP[0, 2] if EP[0, 2] else EP[:, 2], color=P.C_PISTON, lw=1.1)
    a_e.set_ylabel("total particle energy / E(0)")
    a_n.plot(PN[:, 1] * sc.wpe, PN[:, 2], color=P.C_AMBIENT, lw=1.1)
    a_n.set_ylabel("total macroparticles")
    a_n.set_xlabel(r"$t\,\omega_{pe}$")
    for ax in (a_e, a_n):
        P.style_axes(ax)
    fig.suptitle(f"{cfg['meta']['run_id']}: operator balance (heater ↔ injector)")
    fig.tight_layout()
    return P.savefig(fig, f"{cfg['meta']['run_id']}_operator_balance.png", testing=True)


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "runs", "R1_phase", "R1")
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    warns = kinshock.config.validate(cfg, sc)

    print(sc.pretty())
    print("\nvalidation:", "OK" if not warns else "WARNINGS")
    for w in warns:
        print("  !", w)

    config_summary_figure(cfg, sc, warns)

    # data-dependent checks (only if the run has produced output)
    from kinshock import io
    try:
        pfs = io.plotfiles(run_dir)
    except FileNotFoundError:
        print("\nno plotfiles yet — wrote config-summary progress figure only.")
        return
    print(f"\n{len(pfs)} plotfiles found; writing loaded-state + operator-balance figures.")
    loaded_state_figure(cfg, sc, io.load_frame(pfs[0]))
    operator_balance_figure(cfg, sc, run_dir)


if __name__ == "__main__":
    main()
