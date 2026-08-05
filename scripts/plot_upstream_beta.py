#!/usr/bin/env python
"""Plot the upstream beta_0 trajectory under grid heating, and where it saturates.

    python scripts/plot_upstream_beta.py runs/R1_phase/R1_paper_470eV

WHY THIS FIGURE. beta_0 is a Table I primary (0.2) but in a run that under-resolves the
ambient Debye length it is an OUTPUT, not an input: grid heating raises T_0, and
beta_0 = mu0 n k T_0 / B^2 rises with it at fixed n and B. For R1_paper_470eV that means
beta_0 leaves 0.2 within the first tens of t_ab, so any conclusion resting on the upstream
beta has to be read against this curve rather than the config value.

TWO SATURATION NUMBERS, and they differ by 1.7x -- which is the point of the figure:

  * CEILING. The finite-grid instability shuts off once the Debye length reaches the cell
    size. lambda_D ~ sqrt(T), so T_ceiling = T_0 (dz/lambda_D,0)^2 and beta scales with it.
  * ASYMPTOTE. The measured rate falls off linearly in sqrt(T) -- exactly what a drive that
    weakens as lambda_D -> dz should do -- and extrapolates to zero BELOW the ceiling,
    because the instability runs out of drive before lambda_D/dz quite reaches 1.

The fit is the honest weak link: a two-parameter extrapolation over a factor-~11 range in T,
so the SHAPE (rate falling as sqrt(T)) is well constrained while the intercept, and hence the
asymptote, is much less so. The figure marks that rather than hiding it.

Measurement window: the far upstream must stay AHEAD of the piston, so this uses the same
outer-fraction window as grid_heating.py and stops where the piston reaches it. Past that
point the "upstream" is shocked plasma and the temperature is not grid heating at all.

Design notes: two categorical colours only (blue #2a78d6 measured, orange #eb6834 fit),
validated in light mode -- worst all-pairs CVD separation dE 24.7, normal-vision floor 33.6,
contrast >= 3:1 on the #fcfcfb surface. Reference levels are chrome (muted grey), not series
colours, because they are annotation rather than data. Light mode only: a PNG cannot follow
the viewer's theme.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import matplotlib                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402

import kinshock                                      # noqa: E402
from kinshock import io as kio, plotting as P         # noqa: E402
from kinshock.units import QE, MU0, EPS0, ME         # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"


def measure(run_dir, frac, species="amb_electrons"):
    """(t/t_ab, T_eV, piston_reach_fraction) per plotfile, in the outer ``frac`` window."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import grid_heating as gh

    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    z_lo, z_hi = frac * sc.domain_halfwidth, sc.domain_halfwidth
    rows = []
    for p in kio.plotfiles(run_dir):
        fr = kio.load_frame(p)
        ke, _, _ = gh.species_mean_energy(fr, species, ME, z_lo, z_hi)
        if math.isnan(ke):
            continue
        rows.append((fr.time / sc.t_ab, (2.0 / 3.0) * ke / QE,
                     gh.piston_reach(fr, sc) / sc.domain_halfwidth))
    a = np.array(sorted(rows))
    return cfg, sc, a[:, 0], a[:, 1], a[:, 2]


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
    ax.grid(True, color=GRID, linewidth=1.0)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--upstream-frac", type=float, default=0.95,
                    help="keep z > frac*domain as far upstream (default 0.95, tighter than "
                         "grid_heating's 0.75 so the window outruns the piston for longer)")
    ap.add_argument("--reach-max", type=float, default=0.95,
                    help="hard backstop: drop frames once the piston reaches this fraction "
                         "of the domain. The primary cut is the rate-break detector below.")
    ap.add_argument("--break-factor", type=float, default=3.0,
                    help="cut the series where the local dT/dt first exceeds this multiple "
                         "of the trailing median rate (precursor arrival; default 3)")
    ap.add_argument("--t-max", type=float, default=900.0,
                    help="how far to extrapolate the fit, in t_ab")
    args = ap.parse_args()

    cfg, sc, t, T, reach = measure(args.run_dir, args.upstream_frac)
    n, dz, B = sc.namb, sc.dz, sc.B0
    T0 = T[0]

    def beta(TT):
        return MU0 * n * TT * QE / B ** 2

    def lamD(TT):
        return math.sqrt(EPS0 * TT * QE / (n * QE ** 2))

    T_ceil = T0 * (dz / lamD(T0)) ** 2                 # lambda_D = dz
    # Trim to the window where the far upstream is still PRISTINE. A reach threshold alone
    # is too blunt: on R1_paper_470eV, reach <= 0.80 admitted the frame at 144.4 t_ab whose
    # local rate is 5.5 eV/t_ab against the trend's 0.55 -- the shock precursor arriving --
    # and that single point flipped the fitted slope positive, producing an "asymptote"
    # BELOW the current temperature (RESULTS 2026-08-05). Grid heating is smooth and
    # decelerating; the precursor is a step. Cut on that contrast instead.
    keep_reach = reach <= args.reach_max
    t, T, reach = t[keep_reach], T[keep_reach], reach[keep_reach]
    rate = np.diff(T) / np.diff(t)
    cut = len(t)
    for i in range(4, len(rate)):
        trail = np.median(rate[max(0, i - 5):i])
        if trail > 0 and rate[i] > args.break_factor * trail:
            cut = i + 1
            break
    tc, Tc = t[:cut], T[:cut]
    if cut < len(t):
        print(f"  precursor break at t/t_ab = {t[cut]:.1f} "
              f"(rate {rate[cut-1]:.2f} vs trailing median "
              f"{np.median(rate[max(0,cut-6):cut-1]):.2f} eV/t_ab); "
              f"using {cut} of {len(t)} frames")

    # Rate model: dT/dt = a + b*sqrt(T), b < 0 -- the drive weakening as lambda_D -> dz.
    dT = np.diff(Tc) / np.diff(tc)
    Tm = 0.5 * (Tc[1:] + Tc[:-1])
    keep = Tm > 2.5 * T0                               # drop the initial transient
    b, a = np.polyfit(np.sqrt(Tm[keep]), dT[keep], 1)
    corr = float(np.corrcoef(np.sqrt(Tm[keep]), dT[keep])[0, 1])
    T_asym = (a / -b) ** 2

    # integrate the fit forward from the last clean point
    tf = [tc[-1]]
    Tf = [Tc[-1]]
    while tf[-1] < args.t_max:
        step = 0.5
        nxt = Tf[-1] + (a + b * math.sqrt(max(Tf[-1], 1e-6))) * step
        tf.append(tf[-1] + step)
        Tf.append(min(nxt, T_asym))
    tf, Tf = np.array(tf), np.array(Tf)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.9), dpi=200,
                                  facecolor=SURFACE, gridspec_kw={"width_ratios": [1.55, 1]})
    style(ax)
    style(ax2)

    # ---- panel A: beta_0(t) -----------------------------------------------------------
    for lvl, lbl in ((beta(T0), f"Table I target  {beta(T0):.2f}"),
                     (beta(T_asym), f"fitted asymptote  {beta(T_asym):.2f}"),
                     (beta(T_ceil), rf"$\lambda_D\!=\!\Delta z$ ceiling  {beta(T_ceil):.2f}")):
        ax.axhline(lvl, color=MUTED, ls=(0, (5, 4)), lw=1.2, zorder=2)
        ax.text(args.t_max * 0.985, lvl, lbl, va="bottom", ha="right",
                fontsize=8.5, color=INK2)

    ax.axvspan(tc[-1], t[-1], color=GRID, alpha=0.45, zorder=1)
    ax.axvspan(t[-1], args.t_max, color=GRID, alpha=0.9, zorder=1)
    ax.text(0.5 * (tc[-1] + t[-1]), beta(T_ceil) * 0.62, "precursor in window",
            ha="center", va="center", fontsize=7.5, color=MUTED, rotation=90)
    ax.plot(tf, beta(Tf), "-", color=ORANGE, lw=2, zorder=4,
            label=r"fit  $dT/dt = a + b\sqrt{T}$, extrapolated")
    ax.plot(tc, beta(Tc), "o", ms=7, color=BLUE, mec=SURFACE, mew=1.6, zorder=5,
            label="measured (far upstream, ahead of piston)")

    ax.annotate(f"end of clean window\n{tc[-1]:.0f} " + r"$t_{ab}$" +
                f"\n" + rf"$\beta_0={beta(Tc[-1]):.2f}$",
                (tc[-1], beta(Tc[-1])), textcoords="offset points", xytext=(14, -30),
                fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.text(0.5 * (t[-1] + args.t_max), beta(T_ceil) * 0.30,
            "extrapolation only\n(run ends at " + f"{t[-1]:.0f}" + r" $t_{ab}$)",
            ha="center", fontsize=8.5, color=INK2)

    ax.set_xlim(0, args.t_max)
    ax.set_ylim(0, beta(T_ceil) * 1.16)
    ax.set_xlabel(r"$t / t_{ab}$", color=INK2, fontsize=10)
    ax.set_ylabel(r"upstream $\beta_0 = \mu_0 n k T_0 / B_0^2$", color=INK2, fontsize=10)
    ax.set_title(r"$\beta_0$ is an OUTPUT here, not the Table I input", color=INK,
                 fontsize=11.5, loc="left", pad=10)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    for txt in leg.get_texts():
        txt.set_color(INK2)

    # ---- panel B: the evidence for the shape ------------------------------------------
    xb = np.sqrt(beta(Tm[keep]))
    yb = np.diff(beta(Tc))[keep] / np.diff(tc)[keep]
    xx = np.linspace(0, math.sqrt(beta(T_asym)) * 1.05, 50)
    slope = b / math.sqrt(MU0 * n * QE / B ** 2) ** 0   # rate in beta units
    kb = MU0 * n * QE / B ** 2
    ax2.plot(xx, (a + b * xx / math.sqrt(kb)) * kb, "-", color=ORANGE, lw=2, zorder=3,
             label="linear fit")
    ax2.plot(xb, yb, "o", ms=7, color=BLUE, mec=SURFACE, mew=1.6, zorder=4,
             label="measured intervals")
    ax2.axhline(0, color=BASELINE, lw=1.2, zorder=2)
    ax2.axvline(math.sqrt(beta(T_asym)), color=MUTED, ls=(0, (5, 4)), lw=1.2, zorder=2)
    ax2.text(math.sqrt(beta(T_asym)) - 0.03, ax2.get_ylim()[1] * 0.55,
             rf"rate$\to$0 at $\beta_0={beta(T_asym):.2f}$ ", fontsize=8.5, color=INK,
             ha="right", va="center")
    ax2.set_xlabel(r"$\sqrt{\beta_0}\ \ (\propto \lambda_D)$", color=INK2, fontsize=10)
    ax2.set_ylabel(r"$d\beta_0 / d(t/t_{ab})$", color=INK2, fontsize=10)
    ax2.set_title(f"the drive falls off linearly in " + r"$\lambda_D$" +
                  f"   (r = {corr:.3f})", color=INK, fontsize=11.5, loc="left", pad=10)
    leg2 = ax2.legend(frameon=False, fontsize=8.5, loc="upper right")
    for txt in leg2.get_texts():
        txt.set_color(INK2)

    foot = (
        rf"{cfg['meta']['run_id']}    $T_0$ {T0:.1f}$\to${Tc[-1]:.0f} eV over {tc[-1]:.0f} "
        rf"$t_{{ab}}$    $\Delta z/\lambda_D$ {dz/lamD(T0):.2f}$\to${dz/lamD(Tc[-1]):.2f}"
        rf"    fit $dT/dt={a:.3f}{b:+.4f}\sqrt{{T}}$  (r={corr:.3f}, "
        rf"{int(keep.sum())} intervals)" "\n"
        f"The asymptote is a 2-parameter extrapolation over a factor-{Tc[-1]/T0:.0f} range "
        "in T. The SHAPE (rate falling as sqrt(T)) is well constrained; the intercept,\n"
        "and hence the asymptote, much less so. The shaded band is where the shock precursor "
        "reaches the measurement window and the 'upstream' stops being pristine."
    )
    fig.text(0.008, 0.012, foot, fontsize=8, color=MUTED, va="bottom", linespacing=1.5)
    fig.tight_layout(rect=(0, 0.155, 1, 1))
    out = os.path.join(P.media_dir(run_id=cfg["meta"]["run_id"]), "upstream_beta.png")
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  T_0 {T0:.2f} -> {Tc[-1]:.1f} eV   beta_0 {beta(T0):.3f} -> {beta(Tc[-1]):.2f}")
    print(f"  fit dT/dt = {a:.4f} {b:+.5f} sqrt(T)   corr {corr:.3f}")
    print(f"  asymptote  T={T_asym:.1f} eV  beta_0={beta(T_asym):.2f}")
    print(f"  ceiling    T={T_ceil:.1f} eV  beta_0={beta(T_ceil):.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
