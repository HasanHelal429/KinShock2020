#!/usr/bin/env python3
"""False-shock test A, step 3: LINEAR GROWTH RATE of the 1D electrostatic modes ahead of the
piston, from Maxwellian fits to the measured particle populations.

    python scripts/es_growth_rate.py runs/F_phase/fs_470eV_dz0379_t1 runs/F_phase/ft_15keV_t2 \
        --times 0.32 0.40 0.48 [--zwin 0 4] [--dzw 0.25] [--ncomp 3]

Per window [z, z+dzw] d_i0 from the piston front:
  ions (ambient + piston, mass m_i) -> weighted 1D Gaussian mixture in v_z, --ncomp components (EM);
  electrons -> one drifting Maxwellian per species (ambient, piston) from weighted moments;
  densities from sum(w)/window length (checked against n_amb far upstream);
then the unmagnetized ES kinetic dispersion
    D(k, w) = 1 + sum_s (w_ps^2 / (k^2 vt_s^2)) [1 + zeta_s Z(zeta_s)],  zeta_s = (w - k u_s)/(sqrt2 k vt_s)
is solved for complex w on a k grid (Newton from several phase-velocity seeds), keeping the most
unstable root. Reported: gamma_max / w_pi0 and / w_ci0, k at max (k d_i0, k lambda_D0), and v_ph/v_sh.

M1 (ES ion-ion / ion-beam) predicts gamma/w_pi0 similar in a forming and a clean twin (their
normalised distributions are near self-similar) but gamma/w_ci0 larger by w_pi0/w_ci0 = c/v_A.
Caveats: unmagnetized, 1D ES, Maxwellian fits (a fit can miss a sharp beam edge that Penrose sees).
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
from scipy.special import wofz

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock import units as U  # noqa: E402
from penrose_map import load_particles  # noqa: E402


def gmm1d(v, w, k, iters=200, seed=0):
    rng = np.random.default_rng(seed)
    q = np.quantile(v, np.linspace(0.15, 0.95, k))
    mu = q + 1e-6 * rng.standard_normal(k) * (v.std() + 1e-30)
    sd = np.full(k, v.std() / k + 1e-30)
    pi = np.full(k, 1.0 / k)
    W = w / w.sum()
    for _ in range(iters):
        g = pi / (np.sqrt(2 * np.pi) * sd) * np.exp(-0.5 * ((v[:, None] - mu) / sd) ** 2)
        r = g / np.maximum(g.sum(1, keepdims=True), 1e-300)
        rw = r * W[:, None]
        nk = rw.sum(0) + 1e-300
        mu = (rw * v[:, None]).sum(0) / nk
        sd = np.sqrt((rw * (v[:, None] - mu) ** 2).sum(0) / nk) + 1e-12 * abs(v).max()
        pi = nk / nk.sum()
    return pi, mu, sd


def Zf(z):
    return 1j * math.sqrt(math.pi) * wofz(z)


def D(w, k, comps):
    s = 1.0 + 0j
    for wp2, u, vt in comps:
        zeta = (w - k * u) / (math.sqrt(2) * k * vt)
        s += wp2 / (k * k * vt * vt) * (1 + zeta * Zf(zeta))
    return s


def solve(k, comps, seeds_v, gamma0):
    best = None
    for vs in seeds_v:
        w = k * vs + 1j * gamma0
        for _ in range(60):
            f = D(w, k, comps)
            h = 1e-7 * (abs(w) + 1.0)
            dfd = (D(w + h, k, comps) - f) / h
            if dfd == 0:
                break
            step = f / dfd
            w = w - step
            if abs(step) < 1e-10 * (abs(w) + 1):
                break
        if np.isfinite(w) and abs(D(w, k, comps)) < 1e-6 and w.imag > 0:
            if best is None or w.imag > best.imag:
                best = w
    return best


def analyse(run_dir, times, zwin, dzw, ncomp):
    cfg = kinshock.load(run_dir)
    sc = U.derive(cfg)
    vsh = sc.vsh_model
    wpi0 = math.sqrt(sc.namb * U.QE ** 2 / (U.EPS0 * sc.mi))
    th0 = float(cfg["plasma"]["ambient"]["theta_0"])
    lD0 = math.sqrt(th0 * sc.me * U.C ** 2 * U.EPS0 / (sc.namb * U.QE ** 2))
    paths = kio.plotfiles(run_dir)
    ptimes = np.array([kio.load_frame(p, add_applied_field=False).time for p in paths]) * sc.wci0
    name = os.path.basename(run_dir.rstrip("/"))
    print(f"\n== {name}: w_pi0/w_ci0 = {wpi0 / sc.wci0:.0f}, d_i0/lambda_D0 = {sc.di0 / lD0:.0f}")
    out = []
    for t in times:
        i = int(np.argmin(np.abs(ptimes - t)))
        P = load_particles(paths[i], sc)
        za, wa = P["piston_ions"][0], P["piston_ions"][2]
        o = np.argsort(za)
        zp = za[o][np.searchsorted(np.cumsum(wa[o]) / wa.sum(), 0.99)]
        # density check far upstream
        zf, wf = P["amb_ions"][0], P["amb_ions"][2]
        L = zf.max()
        far = zf > L - 1.0 * sc.di0
        nchk = wf[far].sum() / (1.0 * sc.di0) / sc.namb
        print(f"  t*={ptimes[i]:.3f}: far-upstream n_amb check = {nchk:.3f}")
        for a in np.arange(zwin[0], zwin[1] - 1e-9, dzw):
            lo, hi = zp + a * sc.di0, zp + (a + dzw) * sc.di0
            comps = []
            zi = []
            vi = []
            wi = []
            for sp in ("amb_ions", "piston_ions"):
                z, v, w, m = P[sp]
                s = (z >= lo) & (z < hi)
                zi.append(z[s]); vi.append(v[s]); wi.append(w[s])
            vi = np.concatenate(vi); wi = np.concatenate(wi)
            if wi.sum() <= 0 or len(vi) < 100:
                continue
            ni = wi.sum() / (hi - lo)
            pi_, mu, sd = gmm1d(vi, wi, ncomp)
            for p_, m_, s_ in zip(pi_, mu, sd):
                if p_ < 1e-3:
                    continue
                comps.append((p_ * ni * U.QE ** 2 / (U.EPS0 * sc.mi), m_, max(s_, 1e-3 * vsh)))
            for sp in ("amb_electrons", "piston_electrons"):
                z, v, w, m = P[sp]
                s = (z >= lo) & (z < hi)
                if w[s].sum() <= 0:
                    continue
                ne = w[s].sum() / (hi - lo)
                u = np.average(v[s], weights=w[s])
                vt = math.sqrt(np.average((v[s] - u) ** 2, weights=w[s]))
                comps.append((ne * U.QE ** 2 / (U.EPS0 * m), u, vt))
            ks = np.geomspace(5, 3 * sc.di0 / lD0, 60) / sc.di0
            seeds = sorted({*(m_ for _, m_, _ in comps[:ncomp])} |
                           {0.5 * (x + y) for x in mu for y in mu if x < y})
            best = (0.0, np.nan, np.nan)
            for k in ks:
                w = solve(k, comps, seeds, 0.05 * wpi0)
                if w is not None and w.imag > best[0]:
                    best = (w.imag, k, w.real / k)
            g, kb, vph = best
            out.append(dict(name=name, t=ptimes[i], z=a + dzw / 2, g_pi=g / wpi0, g_ci=g / sc.wci0,
                            comps_norm=[(wp2 / wpi0 ** 2, u / vsh, vt / vsh) for wp2, u, vt in comps],
                            wpi0=wpi0, vsh=vsh, lD0=lD0, di0=sc.di0, wci0=sc.wci0,
                            kdi=kb * sc.di0, klD=kb * lD0, vph=vph / vsh,
                            pops=[(round(p_, 3), round(m_ / vsh, 2), round(s_ / vsh, 3))
                                  for p_, m_, s_ in zip(pi_, mu, sd)]))
            r = out[-1]
            print(f"    z={r['z']:5.2f}  gamma/w_pi0={r['g_pi']:.3f}  gamma/w_ci0={r['g_ci']:7.1f}  "
                  f"k d_i0={r['kdi']:7.0f}  k lD0={r['klD']:.2f}  v_ph/v_sh={r['vph']:.2f}  "
                  f"ions(frac,u,vt)/v_sh={r['pops']}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", type=float, nargs="+", required=True)
    ap.add_argument("--zwin", type=float, nargs=2, default=(0.0, 4.0))
    ap.add_argument("--dzw", type=float, default=0.25)
    ap.add_argument("--ncomp", type=int, default=3)
    ap.add_argument("--dump", default=None, help="pickle every window's fitted populations here")
    args = ap.parse_args()
    allr = []
    for r in args.runs:
        allr += analyse(os.path.abspath(r), args.times, args.zwin, args.dzw, args.ncomp)
    if args.dump:
        import pickle
        with open(args.dump, "wb") as fh:
            pickle.dump(allr, fh)
    print("\nSUMMARY (max over windows per run and time)")
    for name in dict.fromkeys(r["name"] for r in allr):
        for t in dict.fromkeys(r["t"] for r in allr if r["name"] == name):
            rr = [r for r in allr if r["name"] == name and r["t"] == t]
            b = max(rr, key=lambda r: r["g_ci"])
            print(f"  {name:22s} t*={t:.3f}  max gamma/w_pi0={b['g_pi']:.3f}  gamma/w_ci0={b['g_ci']:.1f}"
                  f"  at z={b['z']:.2f}, k lD0={b['klD']:.2f}, v_ph/v_sh={b['vph']:.2f};"
                  f"  e-folds per 0.01 t*: {0.01 * b['g_ci']:.2f}")


if __name__ == "__main__":
    main()
