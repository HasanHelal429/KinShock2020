#!/usr/bin/env python3
"""False-shock test A, step 1: WHERE and AT WHAT WAVENUMBER short-scale E_z lives, relative to
the two upstream ion populations, around the false-shock onset.

    python scripts/notch_ez_spectrum.py runs/F_phase/fs_470eV_dz0379_t1 --times 0.25 0.32 0.40 0.48 \
        [--control runs/F_phase/ft_15keV_t2] [--zwin -1 5] [--out media/F_phase/testA]

Per time (columns), three rows:
  1. ambient-ion phase space from PS_ai_zuz (v_z / v_sh vs z / d_i0), asinh-stretched;
  2. raw E_z / (v_sh B0) from the nearest diag_fields dump, NOT smoothed;
  3. sliding-window E_z power spectrum, log10 |E_k|^2 vs (z, k d_i0): Hann windows of
     --win d_i0, step --step d_i0. Reference lines: k d_i0 = c/v_sh and c/(0.5 v_sh) (the ES
     ion-ion scale omega_pi0/v_rel for v_rel = v_sh, v_sh/2), and k lambda_D0 = 1.
z is measured from the piston front (99th percentile of the piston-ion histogram density).

With --control, a second figure is made for the control run at the same t*w_ci0, same layout
and same normalisation, so the two compare panel for panel. Shock features are for the user to
judge by eye; this script measures nothing about them.
Needs the WarpX venv (openpmd_api + yt): source perlmutter/_common.sh && pm_profile.
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import kinshock  # noqa: E402
from kinshock import io as kio  # noqa: E402
from kinshock import units as U  # noqa: E402
from make_phase_hist_movie import load_hist, open_series  # noqa: E402


def piston_front(run_dir, it):
    s = open_series(run_dir, "PS_pi_zuz")
    H, ext, _ = load_hist(s, it)
    prof = H.sum(axis=0)
    z = np.linspace(ext[0], ext[1], H.shape[1] + 1)
    c = np.cumsum(prof) / max(prof.sum(), 1e-300)
    return z[np.searchsorted(c, 0.99) + 1]


def spectrogram(z, ez, win, step):
    dz = z[1] - z[0]
    n = int(round(win / dz)) // 2 * 2
    hann = np.hanning(n)
    centers, rows = [], []
    for i0 in range(0, len(z) - n, max(1, int(round(step / dz)))):
        seg = ez[i0:i0 + n]
        seg = (seg - seg.mean()) * hann
        rows.append(np.abs(np.fft.rfft(seg)) ** 2)
        centers.append(z[i0 + n // 2])
    k = 2 * np.pi * np.fft.rfftfreq(n, dz)
    return np.array(centers), k, np.array(rows)


def frame_set(run_dir, times, zwin, win, step):
    cfg = kinshock.load(run_dir)
    sc = U.derive(cfg)
    th0 = float(cfg["plasma"]["ambient"]["theta_0"])
    lD0 = math.sqrt(th0 * sc.me * U.C ** 2 * U.EPS0 / (sc.namb * U.QE ** 2))
    mi = sc.mi
    ps = open_series(run_dir, "PS_ai_zuz")
    its = np.array(sorted(ps.iterations))
    ps_t = np.array([ps.iterations[i].time * ps.iterations[i].time_unit_SI for i in its]) * sc.wci0
    fdir = os.path.join(run_dir, "diags")
    fpaths = {kio._step_key(os.path.join(fdir, p))[0]: os.path.join(fdir, p)
              for p in os.listdir(fdir) if p.startswith("diag_fields") and p[11:].isdigit()}
    fsteps = np.array(sorted(fpaths))
    f_t = fsteps * sc.dt * sc.wci0
    out = []
    for t in times:
        ip = int(its[np.argmin(np.abs(ps_t - t))])
        H, ext, tt = load_hist(ps, ip)
        zp = piston_front(run_dir, ip)
        fs = int(fsteps[np.argmin(np.abs(f_t - t))])
        fr = kio.load_frame(fpaths[fs], fields=("Ez",), add_applied_field=False)
        z = fr.z_centers
        ez = fr.comps["Ez"] / (sc.vsh_model * sc.B0)
        m = (z >= zp + zwin[0] * sc.di0) & (z <= zp + zwin[1] * sc.di0)
        zc, k, P = spectrogram(z[m], ez[m], win * sc.di0, step * sc.di0)
        out.append(dict(t_ps=tt * sc.wci0, t_f=fs * sc.dt * sc.wci0, H=H, ext=ext, zp=zp,
                        z=(z[m] - zp) / sc.di0, ez=ez[m], zc=(zc - zp) / sc.di0,
                        kdi=k * sc.di0, P=P))
    meta = dict(name=os.path.basename(run_dir.rstrip("/")), di0=sc.di0, vsh=sc.vsh_model,
                mi=mi, c_over_vsh=U.C / sc.vsh_model, kD=sc.di0 / lD0,
                eps_ab=math.sqrt(float(cfg["plasma"]["piston"]["theta_e_heat"])))
    return meta, out


def draw(meta, frames, zwin, path, ezlim=None, plim=None):
    n = len(frames)
    fig, ax = plt.subplots(3, n, figsize=(3.6 * n, 9.5), sharex=True, squeeze=False,
                           gridspec_kw=dict(height_ratios=[1.2, 0.8, 1.2]))
    for j, f in enumerate(frames):
        z0, z1, u0, u1 = f["ext"]
        vs = U.C / meta["vsh"]
        a = ax[0, j]
        a.imshow(np.arcsinh(f["H"] / max(np.percentile(f["H"][f["H"] > 0], 50), 1e-300)),
                 origin="lower", aspect="auto", cmap="magma",
                 extent=[(z0 - f["zp"]) / meta["di0"], (z1 - f["zp"]) / meta["di0"], u0 * vs, u1 * vs])
        a.set_ylim(-0.5, 2.5)
        a.set_title(f"t*w_ci0 = {f['t_ps']:.3f}" + (f"  (E: {f['t_f']:.3f})" if abs(f['t_f'] - f['t_ps']) > 2e-3 else ""))
        a.axhline(1, color="c", lw=0.6, ls="--")
        a = ax[1, j]
        a.plot(f["z"], f["ez"], lw=0.4, color="k")
        if ezlim:
            a.set_ylim(-ezlim, ezlim)
        a = ax[2, j]
        P = np.log10(np.maximum(f["P"], 1e-300)).T
        vmax = plim if plim else np.percentile(P, 99.9)
        a.pcolormesh(f["zc"], f["kdi"][1:], P[1:], shading="nearest", cmap="viridis",
                     vmin=vmax - 6, vmax=vmax)
        a.set_yscale("log")
        a.set_ylim(f["kdi"][1], f["kdi"][-1])
        for kk, ls, lab in ((meta["c_over_vsh"], "--", "c/v_sh"), (2 * meta["c_over_vsh"], ":", "2c/v_sh"),
                            (meta["kD"], "-", "k lD0 = 1")):
            a.axhline(kk, color="w", ls=ls, lw=0.8)
            if j == 0:
                a.text(zwin[0] + 0.05, kk * 1.08, lab, color="w", fontsize=7)
        a.set_xlabel("(z - z_piston) / d_i0")
    ax[0, 0].set_ylabel("ambient ions  v_z / v_sh")
    ax[1, 0].set_ylabel("E_z / (v_sh B0)  [raw]")
    ax[2, 0].set_ylabel("k d_i0   (log10 |E_k|^2, 6 decades)")
    ax[0, 0].set_xlim(*zwin)
    fig.suptitle(f"{meta['name']}  (eps_ab = {meta['eps_ab']:.3f}):  ion phase space, raw E_z, "
                 "E_z k-spectrum along z", y=0.995)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--times", type=float, nargs="+", required=True)
    ap.add_argument("--control", default=None)
    ap.add_argument("--zwin", type=float, nargs=2, default=(-1.0, 5.0))
    ap.add_argument("--win", type=float, default=0.5, help="window length [d_i0]")
    ap.add_argument("--step", type=float, default=0.05, help="window step [d_i0]")
    ap.add_argument("--out", default="media/F_phase/testA")
    args = ap.parse_args()
    runs = [args.run_dir] + ([args.control] if args.control else [])
    data = [frame_set(os.path.abspath(r), args.times, args.zwin, args.win, args.step) for r in runs]
    ezlim = 1.1 * max(np.percentile(np.abs(f["ez"]), 99.9) for _, fr in data for f in fr)
    plim = max(np.percentile(np.log10(np.maximum(f["P"], 1e-300)), 99.9) for _, fr in data for f in fr)
    for meta, fr in data:
        draw(meta, fr, args.zwin, os.path.join(args.out, f"ez_spectrum_{meta['name']}.png"),
             ezlim=ezlim, plim=plim)
        for f in fr:
            print(f"  {meta['name']} t*={f['t_ps']:.3f}  z_piston={f['zp']/meta['di0']:.2f} d_i0  "
                  f"rms E_z/(v_sh B0)={np.sqrt(np.mean(f['ez']**2)):.3e}")


if __name__ == "__main__":
    main()
