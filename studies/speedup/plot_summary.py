#!/usr/bin/env python
"""Plot the optimization sweep into ``media/opt_phase/``.

Reads the same ``out/*/result.txt`` files as make_summary.py, so the figures cannot drift
from the tables.

DESIGN NOTES (so the choices are auditable rather than taste):
  * Palette is two categorical slots, blue #2a78d6 and orange #eb6834, validated in light
    mode: worst all-pairs CVD separation dE 24.7 (protan) against a >=8 target, normal-vision
    floor 33.6 against >=15, contrast >=3:1 on the #fcfcfb surface. A third slot (aqua
    #1baf7a) was dropped because it lands below 3:1 and would have obliged extra relief.
  * LIGHT MODE ONLY, deliberately. A PNG cannot respond to the viewer's theme, and the
    rest of media/ is light; a half-adapted figure is worse than one that commits.
  * Figure 2 is a LOLLIPOP on a log axis, not bars. The configurations span three orders of
    magnitude (0.68 d to 71 d), and bar length must stay proportional to value -- a log bar
    chart lies about ratios. Dots encode by position, so a log axis is legitimate for them.
  * Grid and axes are recessive (hairline #e1e0d9, muted ink #898781); values are direct-
    labelled in text ink, never in the series colour.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "studies" / "speedup" / "out"
MEDIA = ROOT / "media" / "opt_phase"

BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"

STEPS_FULL, CFL_DECK, DRIFT = 3_224_046, 0.75, 1.279
RESULT = re.compile(r"^(\S+)\s+thr=\s*(\d+)\s+n=\s*(\d+)\s+s/step=\s*([0-9.]+)")
CFL_OVERRIDE = re.compile(r"warpx\.cfl=([0-9.]+)")


def points():
    """label -> (mean s/step, cfl). The cfl is needed because dt is proportional to it, so a
    large-dt point reaches 220 t_ab in STEPS_FULL*0.75/cfl steps, not STEPS_FULL."""
    out = {}
    for d in sorted(OUT.iterdir()) if OUT.is_dir() else []:
        f = d / "result.txt"
        if f.is_file():
            m = RESULT.match(f.read_text().strip())
            if m:
                meta = (d / "meta.txt").read_text() if (d / "meta.txt").is_file() else ""
                cm = CFL_OVERRIDE.search(meta)
                out[m.group(1)] = (float(m.group(4)),
                                   float(cm.group(1)) if cm else CFL_DECK)
    return out


def days(sc):
    s, cfl = sc
    return STEPS_FULL * CFL_DECK / cfl * s * DRIFT / 86400.0


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(BASELINE)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=1.0)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK2)


def fig_threads(p, path):
    thr = [4, 8, 12, 16, 20, 24]
    keys = {4: "l1_thr4", 8: "l1_thr8", 12: "l1_thr12", 16: "l1_thr16",
            20: "l1_thr20_clean", 24: "l1_thr24_clean"}
    have = [t for t in thr if keys[t] in p]
    if len(have) < 3:
        print("  threads: too few points, skipped")
        return
    base = p["l1_thr8"][0]
    spd = [base / p[keys[t]][0] for t in have]

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=200, facecolor=SURFACE)
    style(ax)
    ax.yaxis.grid(True, color=GRID, linewidth=1.0)
    ax.set_axisbelow(True)

    # Ideal linear scaling from the 8-thread baseline: chrome, not a data series.
    ax.plot(have, [t / 8 for t in have], linestyle=(0, (5, 4)), linewidth=2,
            color=MUTED, label="ideal linear (from 8 thr)", zorder=2)
    ax.plot(have, spd, "-o", linewidth=2, markersize=8, color=BLUE,
            markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=2,
            label="measured", zorder=3)

    for t, s in zip(have, spd):
        ax.annotate(f"{s:.2f}x", (t, s), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=9, color=INK)
    if 20 in have:
        i = have.index(20)
        ax.annotate("sweet spot — 24 thr buys only 5% more\nfor 75% of a shared 32-core box",
                    (20, spd[i]), textcoords="offset points", xytext=(-18, -42),
                    ha="center", fontsize=8.5, color=INK2,
                    arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=1))

    ax.set_xlabel("OMP threads", color=INK2, fontsize=10)
    ax.set_ylabel("speedup vs 8 threads", color=INK2, fontsize=10)
    ax.set_title("CPU thread scaling — R1_paper_470eV, idle chablis",
                 color=INK, fontsize=12, pad=14, loc="left")
    ax.set_xticks(have)
    ax.set_ylim(0, max(max(spd), 3.0) * 1.18)
    leg = ax.legend(frameon=False, fontsize=9, loc="upper left")
    for t in leg.get_texts():
        t.set_color(INK2)
    fig.text(0.01, 0.02, "mean s/step, 1500 steps/point. No cliff to 24 threads — the "
             "earlier ~20x collapse needed a busy box and a 5x smaller deck.",
             fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    print(f"  wrote {path}")


def fig_cost(p, path):
    rows = [("explicit GPU, 1 box", "l3B_gpu_1box"),
            ("explicit GPU, 235 boxes (default)", "l3A_gpu_default"),
            ("explicit CPU 20 thr", "l1_thr20_clean"),
            ("explicit CPU 8 thr (baseline)", "l1_thr8"),
            ("implicit GPU, cfl 3.0", "l4_gpu_impl_cfl3.0"),
            ("implicit GPU, cfl 7.5", "l4_gpu_impl_cfl7.5"),
            ("implicit GPU, cfl 0.75", "l4_gpu_impl_cfl0.75"),
            ("implicit CPU, cfl 7.5", "l2D_cfl7.5"),
            ("implicit CPU, Picard cfl 0.75", "l2B_time_picard")]
    rows = [(lab, days(p[k])) for lab, k in rows if k in p]
    if len(rows) < 3:
        print("  cost: too few points, skipped")
        return
    rows.sort(key=lambda r: r[1])
    labels = [r[0] for r in rows]
    dvals = [r[1] for r in rows]
    y = list(range(len(rows)))
    # Orange marks the two RECOMMENDATIONS, which are not the same thing: cheapest overall,
    # and cheapest that also conserves energy exactly (no grid heating by construction).
    best = min(range(len(dvals)), key=lambda i: dvals[i])
    best_impl = min((i for i, l in enumerate(labels) if l.startswith("implicit")),
                    key=lambda i: dvals[i], default=None)
    cols = [ORANGE if i in (best, best_impl) else BLUE for i in y]
    days_list = dvals

    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=200, facecolor=SURFACE)
    style(ax)
    ax.set_xscale("log")
    ax.xaxis.grid(True, color=GRID, linewidth=1.0, which="major")
    ax.set_axisbelow(True)

    for yi, d, c in zip(y, days_list, cols):
        ax.plot([min(days_list) * 0.55, d], [yi, yi], color=GRID, linewidth=2, zorder=2)
        ax.plot([d], [yi], "o", markersize=9, color=c, markeredgecolor=SURFACE,
                markeredgewidth=2, zorder=3)
        txt = f"{d:.2f} d" if d < 100 else f"{d/365:.1f} yr"
        ax.annotate(txt, (d, yi), textcoords="offset points", xytext=(12, 0),
                    va="center", fontsize=9, color=INK)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(y) - 0.3)
    ax.set_xlim(min(days_list) * 0.45, max(days_list) * 4.0)
    ax.set_xlabel("projected full-run wall clock (days, log scale)", color=INK2, fontsize=10)
    ax.set_title("Cost of the full 220 t_ab run — orange = the two recommendations",
                 color=INK, fontsize=12, pad=14, loc="left")
    fig.text(0.01, 0.02, "steps(cfl) x mean s/step x 1.279 drift, steps = 3,224,046 x 0.75/cfl (large dt\n"
             "reaches 220 t_ab in fewer steps). Dots on a log axis, not bars: bar length\n"
             "must stay proportional to value.",
             fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.105, 1, 1))
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    p = points()
    if not p:
        print("no benchmark points found", file=sys.stderr)
        return 1
    MEDIA.mkdir(parents=True, exist_ok=True)
    print(f"{len(p)} points ->  {MEDIA}")
    fig_threads(p, MEDIA / "threads_scaling.png")
    fig_cost(p, MEDIA / "lever_cost.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
