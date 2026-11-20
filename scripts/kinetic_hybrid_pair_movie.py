#!/usr/bin/env python3
"""Side-by-side movie: a KINETIC run (left) vs a HYBRID run (right), on DIMENSIONLESS shared axes.

    python scripts/kinetic_hybrid_pair_movie.py runs/E_phase/ek_47keV_t2 runs/hybrid_phase/hc_470eV_dense_th63 \
        [--labels "..." "..."] [--teab 47010 470] [--out media/.../x.mp4]

The two runs need not share physical scales (e.g. the 47 keV kinetic twin vs a 470 eV hybrid): every axis
is normalised by the run's own upstream quantities, and frames are paired by each run's own t*w_ci0.
  * ions      (z/d_i0, v_z/v_A), piston + ambient
  * fields    B_x/B0 and E_z/(v_A B0)
  * electrons T_e/T_e,ab: hybrid = fluid T_e; kinetic = T_e,z = m_e var(v_z) of all electrons per bin
              (one velocity component), plus n_i/n_amb per species (log)
Scales come from my_constants (B0, Mi, namb); d_i0 = c/w_pi(n_amb), v_A = B0/sqrt(mu0 n_amb M_i),
w_ci0 = e B0/M_i. Kinetic particle dumps are `diag<N>`, hybrid ones `diag_phase<N>`; fields from
`diag_fields<N>`. Frame rate = n_frames / --seconds (default 15 s).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)
from kinshock import deck  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock import plotting as P  # noqa: E402

QE, ME, MU0, EPS0, C = 1.602176634e-19, 9.1093837015e-31, 4e-7 * np.pi, 8.8541878128e-12, 2.99792458e8


def scales(run):
    d = deck.parse_inputs(os.path.join(run, "warpx_used_inputs"))
    ns = deck.resolve_constants(d, allow_unresolved=True)
    B0, Mi, namb = float(ns["B0"]), float(ns["Mi"]), float(ns["namb"])
    if "warpx.const_dt" in d:
        dt = float(d["warpx.const_dt"].split()[0])
    else:   # CFL-set dt (kinetic decks): take it from the run log's step lines
        import re
        with open(os.path.join(run, "run.log")) as f:
            m = re.search(r"DT = ([0-9.eE+-]+)", f.read())
        dt = float(m.group(1))
    return dict(B0=B0, Mi=Mi, namb=namb, dt=dt, wci=QE * B0 / Mi,
                vA=B0 / np.sqrt(MU0 * namb * Mi), di=C / np.sqrt(namb * QE**2 / (EPS0 * Mi)))


def dumps(run, prefix):
    dd = os.path.join(run, "diags")
    return {int(p[len(prefix):]): os.path.join(dd, p) for p in os.listdir(dd)
            if p.startswith(prefix) and p[len(prefix):].isdigit()}


def fields(run, s, step, kinetic, teab, smooth_de):
    fl = dumps(run, "diag_fields")
    fs = min(fl, key=lambda k: abs(k - step))
    want = ("Ez", "rho_piston_ions", "rho_amb_ions") + (() if kinetic else ("Te",))
    fr = kio.load_frame(fl[fs], add_applied_field=False, fields=want)
    z = fr.z_centers / s["di"]
    dz = fr.z_edges[1] - fr.z_edges[0]
    nk = max(1, int(round(smooth_de * s["di"] / np.sqrt(s["Mi"] / ME) / dz)))
    sm = (lambda a: np.convolve(a, np.ones(nk) / nk, mode="same")) if nk > 1 else (lambda a: a)
    out = dict(z=z, bx=sm(fr.Bx / s["B0"]), ez=sm(fr.comps["Ez"] / (s["vA"] * s["B0"])),
               npi=sm(fr.comps["rho_piston_ions"] / (QE * s["namb"])),
               nai=sm(fr.comps["rho_amb_ions"] / (QE * s["namb"])))
    if not kinetic:
        out["te"] = sm(fr.comps["Te"]) / teab
    return out


def kinetic_te(fr, s, teab, nbins=400):
    zs, ws, us = [], [], []
    for sp in ("piston_electrons", "amb_electrons"):
        z, uz, w = kio.species_phase_weighted(fr, sp, None, mass=ME)
        if len(z):
            zs.append(z); ws.append(w); us.append(uz * C)
    z, w, v = (np.concatenate(a) for a in (zs, ws, us))
    edges = np.linspace(fr.z_edges[0], fr.z_edges[-1], nbins + 1)
    W = np.histogram(z, edges, weights=w)[0]
    m1 = np.histogram(z, edges, weights=w * v)[0] / np.maximum(W, 1e-300)
    m2 = np.histogram(z, edges, weights=w * v * v)[0] / np.maximum(W, 1e-300)
    T = ME * np.maximum(m2 - m1**2, 0) / QE / teab
    T[W <= 0] = np.nan
    return 0.5 * (edges[1:] + edges[:-1]) / s["di"], T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kinetic_run")
    ap.add_argument("hybrid_run")
    ap.add_argument("--labels", nargs=2, default=None)
    ap.add_argument("--teab", nargs=2, type=float, default=(47010.0, 470.0),
                    help="T_e,ab [eV] used to normalise T_e: kinetic, hybrid")
    ap.add_argument("--vlim", nargs=2, type=float, default=(-5.0, 40.0), help="v_z/v_A axis")
    ap.add_argument("--field-smooth-de", type=float, default=2.0)
    ap.add_argument("--seconds", type=float, default=15.0)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    runs = [os.path.abspath(args.kinetic_run), os.path.abspath(args.hybrid_run)]
    kin = (True, False)
    labels = args.labels or [os.path.basename(r) for r in runs]
    sc = [scales(r) for r in runs]
    ph = [dumps(runs[0], "diag1"), dumps(runs[1], "diag_phase")]   # kinetic particle dumps: diag1<step>
    if not ph[0]:
        raise SystemExit("no diag1<step> particle dumps in the kinetic run")
    ts = lambda c, st: st * sc[c]["dt"] * sc[c]["wci"]
    # walk the run that is coarser in t*, pair with the other's nearest dump, over the common t* range
    s0, s1 = sorted(ph[0]), sorted(ph[1])
    d0, d1 = ts(0, s0[1] - s0[0]), ts(1, s1[1] - s1[0])
    lead = 0 if d0 >= d1 else 1
    tmax = min(ts(0, s0[-1]), ts(1, s1[-1]))
    so = s1 if lead == 0 else s0
    pairs = []
    for st in (s0 if lead == 0 else s1):
        t = ts(lead, st)
        if t > tmax + 1e-9:
            break
        o = min(so, key=lambda q: abs(ts(1 - lead, q) - t))
        pairs.append((st, o) if lead == 0 else (o, st))
    if args.max_frames:
        pairs = pairs[:args.max_frames]
    worst = max(abs(ts(0, a) - ts(1, b)) for a, b in pairs)
    print(f"{len(pairs)} frames to t*w_ci0 {tmax:.2f}; worst t* mismatch {worst:.3f}; "
          f"v_A = {sc[0]['vA']/1e3:.0f} / {sc[1]['vA']/1e3:.0f} km/s; d_i0 = {sc[0]['di']*1e6:.2f} / {sc[1]['di']*1e6:.2f} um")

    F = [[fields(runs[c], sc[c], p[c], kin[c], args.teab[c], args.field_smooth_de) for p in pairs] for c in range(2)]
    zmax = min(F[0][0]["z"][-1], F[1][0]["z"][-1])
    inner = lambda v: v["z"] < zmax - 2.0
    allv = F[0] + F[1]
    blim = 1.1 * max(np.percentile(v["bx"][inner(v)], 99.9) for v in allv)
    bfloor = min(-0.05 * blim, 1.1 * min(np.percentile(v["bx"][inner(v)], 0.1) for v in allv))
    elim = 1.1 * max(np.percentile(np.abs(v["ez"][inner(v)]), 99.9) for v in allv)
    nlim = 1.2 * max(np.percentile(np.maximum(v["npi"], v["nai"])[inner(v)], 99.9) for v in allv)
    z_edges = np.linspace(0.0, zmax, 1201)
    v_edges = np.linspace(args.vlim[0], args.vlim[1], 241)
    framedir = os.path.join(os.environ.get("MOVIE_SCRATCH", os.path.join(ROOT, "media", ".frames")),
                            "kin_hyb_" + os.path.basename(runs[1]))
    os.makedirs(framedir, exist_ok=True)
    for f in os.listdir(framedir):
        os.remove(os.path.join(framedir, f))
    tlim = None
    for k, pair in enumerate(pairs):
        fig, axes = plt.subplots(3, 2, figsize=(15.0, 9.0), sharex=True,
                                 gridspec_kw={"height_ratios": [1.1, 0.8, 0.7]})
        for c in range(2):
            s = sc[c]; D = F[c][k]
            fr = kio.load_frame(ph[c][pair[c]], add_applied_field=False)
            sd = {}
            for sp, key in (("piston_ions", "piston"), ("amb_ions", "ambient")):
                z, uz, w = kio.species_phase_weighted(fr, sp, None, mass=s["Mi"])
                if len(z):
                    sd[key] = (z / s["di"], uz * C / s["vA"], w)
            ax = axes[0, c]
            # v_sh,model = 13.95 v_A in both families (4182 km/s at 470 eV, 41820 at 47 keV)
            P.phase_distribution(ax, sd, z_edges, v_edges, vline=13.95, legend=(c == 0))
            ax.set_ylabel(r"$v_z / v_A$")
            ax.set_title(f"{labels[c]}\n" rf"$t\,\omega_{{ci0}}$ = {fr.time * s['wci']:.3f}"
                         f"   (B0 {s['B0']:.1f} T, v_A {s['vA']/1e3:.0f} km/s)", fontsize=10)
            axf, axe = axes[1, c], axes[1, c].twinx()
            axe.plot(D["z"], D["ez"], color=P.C_AMBIENT, lw=0.6, alpha=0.8, label=r"$E_z/v_AB_0$")
            axe.set_ylim(-elim, elim); axe.set_ylabel(r"$E_z / v_A B_0$")
            axf.plot(D["z"], D["bx"], color="k", lw=1.2, label=r"$B_x/B_0$")
            axf.set_ylim(bfloor, blim); axf.axhline(0, color="0.6", lw=0.6); axf.set_ylabel(r"$B_x/B_0$")
            axf.set_zorder(axe.get_zorder() + 1); axf.patch.set_visible(False); P.style_axes(axf)
            h1, l1 = axf.get_legend_handles_labels(); h2, l2 = axe.get_legend_handles_labels()
            axf.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8, ncol=2, frameon=False)
            axt, axn = axes[2, c], axes[2, c].twinx()
            axn.plot(D["z"], D["npi"], color=P.C_PISTON, lw=0.8, label=r"$n_{pi}/n_{amb}$")
            axn.plot(D["z"], D["nai"], color=P.C_AMBIENT, lw=0.8, label=r"$n_{ai}/n_{amb}$")
            axn.set_yscale("log"); axn.set_ylim(1e-2, max(nlim, 10)); axn.set_ylabel(r"$n_i/n_{amb}$")
            if kin[c]:
                zt, te = kinetic_te(fr, s, args.teab[c])
                lab = r"$T_{e,z}/T_{e,ab}$ (kinetic, 1 comp.)"
            else:
                zt, te = D["z"], D["te"]
                lab = r"$T_e/T_{e,ab}$ (fluid)"
            axt.plot(zt, te, color="k", lw=1.1, label=lab)
            axt.set_ylim(0, 3.0); axt.set_ylabel(r"$T_e / T_{e,ab}$")
            axt.set_zorder(axn.get_zorder() + 1); axt.patch.set_visible(False); P.style_axes(axt)
            h1, l1 = axt.get_legend_handles_labels(); h2, l2 = axn.get_legend_handles_labels()
            axt.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=7, ncol=3, frameon=False)
            axes[2, c].set_xlabel(r"$z / d_{i0}$")
            axes[2, c].set_xlim(0, zmax)
            if c == 1:
                tw = fr.time * s["wci"]
        fig.suptitle(rf"$t\,\omega_{{ci0}}$ = {tw:.2f}  (each run's own gyrofrequency; dimensionless axes)", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(framedir, f"frame_{k:03d}.png"), dpi=100)
        plt.close(fig)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fps = max(1, int(round(len(pairs) / args.seconds)))
    res = P.encode(framedir, args.out, fps=fps)
    shutil.rmtree(framedir, ignore_errors=True)
    return res


if __name__ == "__main__":
    main()
