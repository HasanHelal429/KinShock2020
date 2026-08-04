#!/usr/bin/env python3
"""Visualize the 3-way half-domain wall crosscheck (RESULTS.md 2026-07-23).

full = R1_core (z>=0, reference) | spec = R1_core_half (specular, artifact)
sym  = R1_core_half_sym (pi-rotation symmetry wall, the fix)

Four panels:
  A  R1 shock-front trajectory z/rho_i0 vs t*wci0 (does sym track full better than spec?)
  B  R1 reflected-ambient-ion fraction G(t) + onset t*_1 markers
  C  R1 key metrics as deviation-from-full (%) -- shorter sym bar = closer to full = fix works
  D  R0 near-wall B-field components (full/spec/sym) -- the completed smoke-tier confirmation

Writes media/testing/crosscheck_3way.png. Run after tmp/crosscheck_3way.py has data.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import kinshock
from kinshock import io
from crosscheck_3way import analyze, RUNS   # reuse the metric engine

# Fixed-order, CVD-safe categorical colors (never cycled): full=neutral reference,
# spec=vermillion (artifact), sym=blue (fix). Identity also carried by legend + labels.
COL = {"full": "#8a8a8a", "spec": "#d55e00", "sym": "#0072b2"}
NAME = {"full": "full R1_core (ref)", "spec": "half: specular", "sym": "half: symmetry (fix)"}
ORDER = ["full", "spec", "sym"]

R0_RUNS = {"full": "runs/R0_phase/R0", "spec": "runs/R0_phase/R0_half", "sym": "runs/R0_phase/R0_half_sym"}


def _style(ax):
    ax.grid(True, alpha=0.25, linewidth=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def r0_near_wall():
    """{tag: {'Bx':.., 'By':.., 'total':..}} near-wall peak /B0 for the R0 trio."""
    out = {}
    for tag, rd in R0_RUNS.items():
        rd = os.path.join(ROOT, rd)
        if not os.path.isdir(os.path.join(rd, "diags")):
            continue
        cfg = kinshock.load(rd); sc = kinshock.units.derive(cfg)
        nw = 3.0 * float(cfg["geometry"]["slab_halfwidth_di"]) * sc.di / sc.de
        fr = [io.load_frame(p) for p in io.plotfiles(rd)][-1]
        zc = np.asarray(fr.z_centers) / sc.de
        m = (zc > 0) & (zc < nw)
        out[tag] = {"Bx": np.nanmax(np.abs(fr.Bx[m])) / sc.B0,
                    "By": np.nanmax(np.abs(fr.By[m])) / sc.B0,
                    "total": np.nanmax(fr.Bperp[m]) / sc.B0}
    return out


def main():
    tags = [t for t in ORDER if os.path.isdir(os.path.join(ROOT, RUNS[t], "diags"))]
    print("[plot] R1 runs available:", tags, flush=True)
    res = {t: analyze(t, RUNS[t]) for t in tags}
    r0 = r0_near_wall()

    fig = plt.figure(figsize=(13.5, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.24,
                          left=0.07, right=0.985, top=0.90, bottom=0.08)
    axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

    # --- A: front trajectory ---
    edge = max((res[t]["edge_rho"] for t in tags), default=None)
    tcm = res[tags[0]]["t_clean_max"] if tags else 2.25
    for t in tags:
        r = res[t]
        tw = np.array(r["traj_tw"]); fr = np.array(r["traj_front_rho"])
        g = np.isfinite(fr)
        axA.plot(tw[g], fr[g], "-o", color=COL[t], ms=3.2, lw=2.0,
                 label=NAME[t], zorder=3 if t != "full" else 2)
    if edge is not None:
        axA.axhline(edge, color="#555", ls=":", lw=1.2)
        axA.text(0.98, edge, "domain edge", va="bottom", ha="right",
                 fontsize=8, color="#555", transform=axA.get_yaxis_transform())
    axA.axvspan(tcm, 3.0, color="#000", alpha=0.05, zorder=0)
    axA.text(tcm + 0.02, 0.04, "boundary-\ncontaminated", fontsize=7.5, color="#777",
             transform=axA.get_xaxis_transform(), va="bottom")
    axA.set_xlabel(r"$t\,\omega_{ci0}$"); axA.set_ylabel(r"shock front  $z/\rho_{i0}$")
    axA.set_title("A. Shock-front trajectory", fontweight="bold", loc="left")
    axA.set_xlim(0, 2.9); axA.legend(frameon=False, fontsize=9, loc="upper left")
    _style(axA)

    # --- B: reflected fraction G(t) + onset markers ---
    label_y = {"full": 0.62, "spec": 0.50, "sym": 0.38}  # stagger to avoid overlap
    for t in tags:
        r = res[t]
        tw = np.array(r["G_tw"]); G = np.array(r["G"])
        axB.plot(tw, G, "-", color=COL[t], lw=2.0, label=NAME[t])
        ts = r["tstar_wci0"]
        if np.isfinite(ts):
            axB.axvline(ts, color=COL[t], ls="--", lw=1.4, alpha=0.9)
            axB.text(0.03, label_y[t], f"$t^*_1$={ts:.2f}", color=COL[t], fontsize=8.5,
                     va="center", ha="left", fontweight="bold",
                     transform=axB.transAxes)
    axB.set_xlabel(r"$t\,\omega_{ci0}$")
    axB.set_ylabel(r"$G=N_{a,\mathrm{refl}}/N_{a,\mathrm{tot}}$")
    axB.set_title(r"B. Reflected ambient ions & onset $t^*_1$", fontweight="bold", loc="left")
    axB.set_xlim(0, 2.9); axB.legend(frameon=False, fontsize=9, loc="upper left")
    _style(axB)

    # --- C: R1 metrics as deviation from full (%) ---
    METRICS = [("v_sh_Cs", r"$v_{sh}$"), ("n_comp", r"$n$-comp"),
               ("b_comp", r"$B$-comp"), ("near_wall_bperp", "near-wall\n$B_\\perp$"),
               ("tstar_wci0", r"$t^*_1$"), ("zstar_rhoi0", r"$z^*_1$")]
    have_full = "full" in res
    others = [t for t in ("spec", "sym") if t in res]
    x = np.arange(len(METRICS)); w = 0.38
    if have_full:
        for j, t in enumerate(others):
            devs = [100.0 * (res[t][k] / res["full"][k] - 1.0)
                    if res["full"][k] else np.nan for k, _ in METRICS]
            off = (j - (len(others) - 1) / 2) * w
            bars = axC.bar(x + off, devs, w, color=COL[t], label=NAME[t], zorder=3)
            for b, d in zip(bars, devs):
                if np.isfinite(d):
                    axC.annotate(f"{d:+.0f}%", (b.get_x() + b.get_width() / 2,
                                 b.get_height()), fontsize=7.5, ha="center",
                                 va="bottom" if d >= 0 else "top",
                                 xytext=(0, 2 if d >= 0 else -2),
                                 textcoords="offset points", color="#333")
        axC.axhline(0, color=COL["full"], lw=1.6)
        axC.text(len(METRICS) - 0.5, 0, "  full = 0", color=COL["full"], fontsize=8,
                 va="center", ha="left")
        axC.set_xticks(x); axC.set_xticklabels([lbl for _, lbl in METRICS], fontsize=8.5)
        axC.set_ylabel("deviation from full-domain  [%]")
        axC.set_title("C. R1 metrics vs full  (shorter = closer = better)",
                      fontweight="bold", loc="left")
        axC.legend(frameon=False, fontsize=9)
    else:
        axC.text(0.5, 0.5, "full R1_core not available", ha="center", transform=axC.transAxes)
    _style(axC)

    # --- D: R0 near-wall B components (completed smoke-tier confirmation) ---
    comps = [("Bx", r"$|B_x|$"), ("By", r"$|B_y|$"), ("total", r"$B_\perp$")]
    xr = np.arange(len(comps)); wr = 0.26
    r0tags = [t for t in ORDER if t in r0]
    for j, t in enumerate(r0tags):
        vals = [r0[t][c] for c, _ in comps]
        off = (j - (len(r0tags) - 1) / 2) * wr
        bars = axD.bar(xr + off, vals, wr, color=COL[t], label=NAME[t], zorder=3)
        for b, v in zip(bars, vals):
            axD.annotate(f"{v:.1f}", (b.get_x() + b.get_width() / 2, b.get_height()),
                         fontsize=7.5, ha="center", va="bottom",
                         xytext=(0, 2), textcoords="offset points", color="#333")
    axD.set_xticks(xr); axD.set_xticklabels([lbl for _, lbl in comps])
    axD.set_ylabel(r"near-wall peak  $/B_0$   ($z<3\,$slab)")
    axD.set_title("D. R0 near-wall field: specular over-shoots $B_y$, symmetry fixes it",
                  fontweight="bold", loc="left", fontsize=10.5)
    axD.legend(frameon=False, fontsize=8.5)
    _style(axD)

    sym_note = "" if "sym" in res else "  [sym pending]"
    nf = res.get("sym", {}).get("n_frames", "?")
    fig.suptitle("Half-domain z=0 wall: specular vs π-rotation symmetry vs full domain"
                 f"   (3-way crosscheck{sym_note}; sym n_frames={nf})",
                 fontsize=13.5, fontweight="bold")
    out = os.path.join(ROOT, "media", "testing", "crosscheck_3way.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=135, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
