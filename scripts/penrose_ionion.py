"""Audit test A: is the ion velocity distribution at the front Penrose-unstable to the
electrostatic ion-ion mode, and at what wavelength? (runs/F_phase/FALSE_SHOCK_AUDIT.md)

Per z column of the binned phase space (PS_ai + PS_pi for ions, PS_ae + PS_pe for electrons):
F(v) = normalized ion v_z distribution; T_e,loc from the electron u_z variance (1D); then

    P(v0) = int (F - F(v0)) / (v - v0)^2 dv   at each interior local minimum v0 of F
    S     = max P * C_s^2,   C_s^2 = T_e,loc / m_i     (unstable for some k iff S > 1)
    k_c^2 = omega_pi^2 P - 1/lambda_De^2             (unstable band k < k_c)

i.e. the unmagnetized 1D electrostatic dielectric 1 + 1/(k lambda_De)^2 - (omega_pi/k)^2 int F'/(v-u)
with Boltzmann-like electron shielding. Columns need n_i > --nmin * n_amb. Prints, per dump,
the max S over z, where it sits (z/d_i0 and relative to the ion front) and 2 pi/k_c in d_i0.
Run with the WarpX venv (openpmd_api).

    python scripts/penrose_ionion.py runs/F_phase/fs_470eV_dz0379_t1 [--every 5]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kinshock import units  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from make_phase_hist_movie import load_hist, open_series  # noqa: E402

EPS0, QE, C = 8.8541878128e-12, 1.602176634e-19, 299792458.0


def smooth(f, w=1.5):
    x = np.arange(-4, 5)
    g = np.exp(-0.5 * (x / w) ** 2)
    return np.convolve(f, g / g.sum(), mode="same")


def penrose(F, v):
    """max over interior local minima of int (F - F0)/(v - v0)^2 dv, with F normalized."""
    best = 0.0
    dv = v[1] - v[0]
    for i in range(1, len(F) - 1):
        if F[i] <= F[i - 1] and F[i] <= F[i + 1] and F[:i].max() > F[i] and F[i + 1:].max() > F[i]:
            d = np.delete(v - v[i], i)
            num = np.delete(F - F[i], i)
            best = max(best, float(np.sum(num / d**2) * dv))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--every", type=int, default=5)
    ap.add_argument("--nmin", type=float, default=0.3)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(os.path.join(a.run_dir, "config.yaml")))
    s = units.derive(cfg)
    me, mi = s.me, s.mi
    S = {k: open_series(a.run_dir, f"PS_{k}_zuz") for k in ("ai", "pi", "ae", "pe")}
    its = sorted(S["ai"].iterations)[:: a.every]
    print(f"# {os.path.basename(a.run_dir.rstrip('/'))}: eps_ab={np.sqrt(s.Te_ab_eV / 510998.95 * ME_FAC(cfg)):.3g} "
          f"d_i0={s.di0:.3g} m  n_amb={s.namb:.3g}")
    print("t*wci0 | S_max | z/d_i0 | z - v_sh t [d_i0] | 2pi/k_c [d_i0] | T_e,loc [eV] | n_i/n_amb")
    for it in its:
        Hi = None
        for k in ("ai", "pi"):
            H, ext, t = load_hist(S[k], it)
            Hi = H if Hi is None else Hi + H
        He = None
        for k in ("ae", "pe"):
            H, exte, _ = load_hist(S[k], it)
            He = H if He is None else He + H
        nv, nz = Hi.shape
        dz = (ext[1] - ext[0]) / nz
        z = ext[0] + (np.arange(nz) + 0.5) * dz
        vi = (ext[2] + (np.arange(nv) + 0.5) * (ext[3] - ext[2]) / nv) * C   # u ~ v for ions
        nve = He.shape[0]
        ue = exte[2] + (np.arange(nve) + 0.5) * (exte[3] - exte[2]) / nve
        ni = Hi.sum(0) / dz
        ne = He.sum(0) / dz
        front = s.vsh_model * t                                   # model shock front v_sh t
        best = (0.0, np.nan, np.nan, np.nan, np.nan)
        for j in np.nonzero(ni > a.nmin * s.namb)[0]:
            w = He[:, j]
            if w.sum() <= 0:
                continue
            um = np.sum(w * ue) / w.sum()
            Te = me * C**2 * np.sum(w * (ue - um) ** 2) / w.sum() / QE         # eV
            F = smooth(Hi[:, j].astype(float))
            if F.sum() <= 0:
                continue
            F = F / (F.sum() * (vi[1] - vi[0]))
            P = penrose(F, vi)
            Sv = P * Te * QE / mi
            if Sv > best[0]:
                wpi2 = ni[j] * QE**2 / (EPS0 * mi)
                lde2 = EPS0 * Te * QE / (max(ne[j], 1.0) * QE**2)
                kc2 = wpi2 * P - 1.0 / lde2
                lam = 2 * np.pi / np.sqrt(kc2) / s.di0 if kc2 > 0 else np.nan
                best = (Sv, z[j], lam, Te, ni[j])
        tw = t * s.wci0
        print(f"{tw:.3f} | {best[0]:.3g} | {best[1] / s.di0:.3f} | {(best[1] - front) / s.di0:.3f} | "
              f"{best[2]:.3g} | {best[3]:.4g} | {best[4] / s.namb:.3g}")


def ME_FAC(cfg):
    return float(cfg["reference"].get("electron_mass_factor", 1.0))


if __name__ == "__main__":
    main()
