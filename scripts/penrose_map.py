#!/usr/bin/env python3
"""False-shock test A, step 2: is the plasma ahead of the piston ELECTROSTATICALLY UNSTABLE in 1D,
judged from the measured particle distributions alone? (Penrose criterion.)

    python scripts/penrose_map.py runs/F_phase/fs_470eV_dz0379_t1 runs/F_phase/ft_15keV_t2 \
        --times 0.24 0.32 0.40 0.48 [--zwin -0.5 4] [--dzw 0.1] [--out media/F_phase/testA]

For each particle dump and each window [z, z + dzw] d_i0 (z from the piston front), build the
charge-weighted 1D distribution  F(v) = sum_s (q_s^2 / m_s) sum_i w_i delta(v - v_z,i)  over ALL
four species (ambient + piston ions and electrons), on a uniform v grid, Gaussian-smoothed by
--sigma v_sh. Penrose (1960): a 1D Vlasov-Poisson plasma is unstable iff F has a minimum v0 with
    P(v0) = PV int (F(v) - F(v0)) / (v - v0)^2 dv  > 0,
and the unstable waves then have phase velocity near v0. Reported: max over minima of
P(v0) / [int F dv / v_sh^2] (dimensionless) and its v0 / v_sh. No fit, no model.

Output per run: one figure, ambient-ion phase space (from the particles) per time with the
windows where P > 0 marked along the top and coloured by v0, plus F(v) in three windows.
Caveats: 1D electrostatic only (no B, no oblique modes); the smoothing sets the smallest
velocity structure resolved -- rerun with --sigma halved to check a verdict is not set by it.
Shock features (false shock) are for the user to judge by eye.
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock import units as U  # noqa: E402

SPECIES = (("amb_ions", "i"), ("piston_ions", "i"), ("amb_electrons", "e"), ("piston_electrons", "e"))


def load_particles(path, sc):
    ad = kio.load_frame(path, add_applied_field=False).ds.all_data()
    out = {}
    for sp, kind in SPECIES:
        m = sc.mi if kind == "i" else sc.me
        try:
            z = np.asarray(ad[(sp, "particle_position_x")])
            pz = np.asarray(ad[(sp, "particle_momentum_z")])
            w = np.asarray(ad[(sp, "particle_weight")])
        except Exception:
            continue
        u = pz / (m * U.C)
        v = u / np.sqrt(1 + u * u) * U.C
        out[sp] = (z, v, w, m)
    return out


def penrose(vgrid, F, vlo, vhi, floor_frac):
    """max over interior minima of PV int (F - F0)/(v - v0)^2 dv, and the v0 achieving it.

    Only minima with v0 in [vlo, vhi] and F(v0) >= floor_frac * max F over that range count:
    a minimum where F ~ 0 (a sampling gap in a sparse tail) makes P > 0 trivially.
    """
    dv = vgrid[1] - vgrid[0]
    band = (vgrid >= vlo) & (vgrid <= vhi)
    floor = floor_frac * F[band].max()
    i_min = np.where((F[1:-1] < F[:-2]) & (F[1:-1] <= F[2:]))[0] + 1
    i_min = i_min[band[i_min] & (F[i_min] >= floor)]
    best, v0b = -np.inf, np.nan
    for i in i_min:
        d = vgrid - vgrid[i]
        d[i] = np.inf
        p = np.sum((F - F[i]) / d ** 2) * dv
        if p > best:
            best, v0b = p, vgrid[i]
    return best, v0b


def analyse(run_dir, times, zwin, dzw, sigma, floor_frac, esmooth):
    cfg = kinshock.load(run_dir)
    sc = U.derive(cfg)
    vsh = sc.vsh_model
    vmax_e = 12 * math.sqrt(float(cfg["plasma"]["piston"]["theta_e_heat"])) * U.C
    dv = 0.005 * vsh
    vgrid = np.arange(-vmax_e, vmax_e + dv, dv)
    paths = kio.plotfiles(run_dir)
    ptimes = np.array([kio.load_frame(p, add_applied_field=False).time for p in paths]) * sc.wci0
    res = []
    for t in times:
        p = paths[int(np.argmin(np.abs(ptimes - t)))]
        tt = ptimes[int(np.argmin(np.abs(ptimes - t)))]
        P = load_particles(p, sc)
        za, wa = P["piston_ions"][0], P["piston_ions"][2]
        o = np.argsort(za)
        zp = za[o][np.searchsorted(np.cumsum(wa[o]) / wa.sum(), 0.99)]
        edges = zp + np.arange(zwin[0], zwin[1] + 1e-9, dzw) * sc.di0
        rows = []
        rng = np.random.default_rng(1)
        half = {sp: rng.random(len(q[0])) < 0.5 for sp, q in P.items()}

        def build(a, b, part):
            # Each species smoothed on ITS OWN scale: ions by sigma*v_sh, electrons by
            # esmooth * their local rms velocity (finer electron structure is sampling noise
            # at these counts, and it dominates F in the ion band through the 1/m_e weight).
            F = np.zeros_like(vgrid)
            for sp, (z, v, w, m) in P.items():
                sel = (z >= a) & (z < b)
                if part is not None:
                    sel &= half[sp] if part == 0 else ~half[sp]
                if sel.sum() < 2:
                    continue
                h, _ = np.histogram(v[sel], bins=len(vgrid), range=(vgrid[0] - dv / 2, vgrid[-1] + dv / 2),
                                    weights=w[sel])
                if m > 10 * sc.me:
                    s = sigma * vsh
                else:
                    s = max(esmooth * math.sqrt(np.average(v[sel] ** 2, weights=w[sel])), sigma * vsh)
                F += gaussian_filter1d(h / (m * dv), s / dv)
            return F

        for a, b in zip(edges[:-1], edges[1:]):
            F = build(a, b, None)
            norm = F.sum() * dv / vsh ** 2
            pv, v0 = penrose(vgrid, F, -0.5 * vsh, 3.0 * vsh, floor_frac)
            # half-split robustness: both independent halves must agree
            ok = False
            if norm > 0 and pv > 0:
                p0, a0 = penrose(vgrid, build(a, b, 0), -0.5 * vsh, 3.0 * vsh, floor_frac)
                p1, a1 = penrose(vgrid, build(a, b, 1), -0.5 * vsh, 3.0 * vsh, floor_frac)
                ok = p0 > 0 and p1 > 0 and abs(a0 - v0) < 0.1 * vsh and abs(a1 - v0) < 0.1 * vsh
            rows.append(dict(zc=(0.5 * (a + b) - zp) / sc.di0,
                             P=(pv / norm if norm > 0 else np.nan) if ok else -abs(pv / norm) if norm > 0 else np.nan,
                             P_raw=pv / norm if norm > 0 else np.nan, v0=v0 / vsh, F=F))
        ai = P["amb_ions"]
        res.append(dict(t=tt, zp=zp, rows=rows, z_ai=(ai[0] - zp) / sc.di0, v_ai=ai[1] / vsh, w_ai=ai[2]))
    meta = dict(name=os.path.basename(run_dir.rstrip("/")), vgrid=vgrid / vsh,
                eps_ab=math.sqrt(float(cfg["plasma"]["piston"]["theta_e_heat"])))
    return meta, res


def draw(meta, res, zwin, path):
    n = len(res)
    fig, ax = plt.subplots(2, n, figsize=(3.8 * n, 7.5), squeeze=False,
                           gridspec_kw=dict(height_ratios=[1.3, 1]))
    for j, r in enumerate(res):
        a = ax[0, j]
        a.hist2d(r["z_ai"], r["v_ai"], bins=[400, 300], range=[zwin, [-0.5, 2.5]],
                 weights=r["w_ai"], cmap="magma", norm=matplotlib.colors.PowerNorm(0.3))
        zc = np.array([q["zc"] for q in r["rows"]])
        Pv = np.array([q["P"] for q in r["rows"]])
        v0 = np.array([q["v0"] for q in r["rows"]])
        uns = Pv > 0
        sc_ = a.scatter(zc[uns], np.full(uns.sum(), 2.35), c=v0[uns], cmap="cool", vmin=-0.5, vmax=2.5,
                        marker="s", s=14)
        a.scatter(zc[uns], v0[uns], c="cyan", s=6, marker="x", lw=0.6)
        a.set_title(f"t*w_ci0 = {r['t']:.3f}\nrobustly unstable: {uns.sum()}/{len(uns)}", fontsize=9)
        a.set_xlabel("(z - z_piston) / d_i0")
        if j == 0:
            a.set_ylabel("ambient ions v_z / v_sh\n(top strip: Penrose-unstable windows;\n x = phase velocity v0)")
        b = ax[1, j]
        good = [q for q in r["rows"] if np.isfinite(q["P"])]
        picks = sorted(good, key=lambda q: -q["P"])[:1] + [good[len(good) // 4], good[len(good) // 2]]
        for q in picks:
            F = q["F"] / q["F"].max()
            b.plot(meta["vgrid"], F, lw=0.8,
                   label=f"z={q['zc']:.2f}: P={q['P']:+.2g}" + (f", v0={q['v0']:.2f}" if np.isfinite(q['v0']) else ""))
        b.set_xlim(-0.5, 2.5)
        b.set_yscale("log")
        b.set_ylim(1e-4, 1.5)
        b.set_xlabel("v / v_sh")
        b.legend(fontsize=6)
        if j == 0:
            b.set_ylabel("F(v) = sum q^2/m f_s  (norm.)")
    fig.colorbar(sc_, ax=ax[0, :].tolist(), shrink=0.6, label="v0 / v_sh")
    fig.suptitle(f"{meta['name']} (eps_ab = {meta['eps_ab']:.3f}): 1D electrostatic Penrose test from the "
                 "particles, all four species", y=0.995)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", type=float, nargs="+", required=True)
    ap.add_argument("--zwin", type=float, nargs=2, default=(-0.5, 4.0))
    ap.add_argument("--dzw", type=float, default=0.1)
    ap.add_argument("--sigma", type=float, default=0.02, help="velocity smoothing [v_sh]")
    ap.add_argument("--esmooth", type=float, default=0.2,
                    help="electron smoothing as a fraction of the window's electron rms velocity")
    ap.add_argument("--floor", type=float, default=0.01,
                    help="minima count only where F(v0) >= floor * max F in the ion band")
    ap.add_argument("--out", default="media/F_phase/testA")
    args = ap.parse_args()
    for r in args.runs:
        meta, res = analyse(os.path.abspath(r), args.times, args.zwin, args.dzw, args.sigma, args.floor, args.esmooth)
        draw(meta, res, args.zwin, os.path.join(args.out, f"penrose_{meta['name']}_s{args.sigma:g}_f{args.floor:g}_e{args.esmooth:g}.png"))
        for q in res:
            P = np.array([w["P"] for w in q["rows"]])
            zc = np.array([w["zc"] for w in q["rows"]])
            v0 = np.array([w["v0"] for w in q["rows"]])
            u = P > 0
            print(f"  {meta['name']} t*={q['t']:.3f}: unstable {u.sum()}/{len(u)} windows"
                  + (f", z in [{zc[u].min():.2f}, {zc[u].max():.2f}] d_i0, v0/v_sh in "
                     f"[{np.nanmin(v0[u]):.2f}, {np.nanmax(v0[u]):.2f}]" if u.any() else ""))


if __name__ == "__main__":
    main()
