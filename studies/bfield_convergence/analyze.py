#!/usr/bin/env python3
"""Compare far-upstream B-fluctuation across convergence variants.

Reads the variant run dirs produced by run_variants.sh, measures RMS(dBx) in the
cold far-upstream zone at matched physical time, and writes a bar+spectrum figure.
Collapse under filter/shape/finer-dz => numerical; invariance => physical.

Usage: python studies/bfield_convergence/analyze.py <out_dir> [--ref-run runs/R1_phase/R1_core]
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "src"))
import kinshock  # noqa: E402
from kinshock import io  # noqa: E402

def smooth(a, k=15): return np.convolve(a, np.ones(k)/k, mode="same")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir", help="scratch dir with baseline/filt8/shape3/finer_dz run subdirs")
    ap.add_argument("--config", default=os.path.join(ROOT, "runs", "R1_phase", "R1_core"),
                    help="run dir to derive scales from (default runs/R1_phase/R1_core)")
    ap.add_argument("--out", default=os.path.join(ROOT, "media", "testing", "bfield_convergence.png"))
    args = ap.parse_args()
    sc = kinshock.units.derive(kinshock.load(args.config)); de = sc.de
    variants = [("baseline", "tab:red"), ("filt8", "tab:blue"),
                ("shape3", "tab:green"), ("finer_dz", "tab:purple")]
    def front(fr):
        zc = np.asarray(fr.z_centers)/de; na = io.species_density(fr, "amb_ions")/sc.namb
        m = zc > 50; a = zc[m][smooth(na[m]) > 1.5]; return a.max() if a.size else np.nan
    # common latest t*wci across variants
    frames = {}
    for name, _ in variants:
        pf = io.plotfiles(os.path.join(args.out_dir, name))
        if pf: frames[name] = pf
    if not frames: sys.exit(f"no variant runs found under {args.out_dir}")
    tmax = min(io.load_frame(pf[-1]).time*sc.wci0 for pf in frames.values())
    print(f"comparing at t*wci~{tmax:.2f}   (zone: cold far-upstream, front+600..+1400 d_e)")
    print(f"{'variant':10s} {'dBx_rms':>8s}  vs baseline")
    res = {}
    for name, c in variants:
        if name not in frames: continue
        pf = frames[name]
        fr = min((io.load_frame(p) for p in pf), key=lambda f: abs(f.time*sc.wci0 - tmax))
        zc = np.asarray(fr.z_centers)/de; bx = fr.Bx/sc.B0; f0 = front(fr)
        m = (zc > f0+600) & (zc < min(f0+1400, zc.max()-20))
        rms = (bx[m]-bx[m].mean()).std() if m.sum() > 20 else float("nan")
        res[name] = (rms, c, fr, f0)
    base = res.get("baseline", (float("nan"),))[0]
    for name, c in variants:
        if name in res:
            r = res[name][0]; print(f"{name:10s} {r:8.3f}  {100*(r-base)/base:+.0f}%")
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    names = [n for n, _ in variants if n in res]
    vals = [res[n][0] for n in names]; cols = [res[n][1] for n in names]
    b = ax[0].bar(range(len(names)), vals, color=cols, alpha=0.85)
    for bar, v in zip(b, vals):
        ax[0].annotate(f"{v:.2f}\n{100*(v-base)/base:+.0f}%", (bar.get_x()+bar.get_width()/2, v),
                       ha="center", va="bottom", fontsize=9)
    ax[0].set_xticks(range(len(names))); ax[0].set_xticklabels(names)
    ax[0].set_ylabel("RMS dBx/B0 (cold far-upstream)")
    ax[0].set_title(f"Far-upstream fluctuation at t*wci~{tmax:.2f}\ncollapse=numerical, invariant=physical")
    for name in names:
        _, c, fr, f0 = res[name]; zc = np.asarray(fr.z_centers)/de; bx = fr.Bx/sc.B0
        m = (zc > f0+600) & (zc < min(f0+1400, zc.max()-20)); seg = (bx[m]-bx[m].mean())*np.hanning(m.sum())
        dz = float(zc[1]-zc[0]); k = np.fft.rfftfreq(len(seg), d=dz)[1:]
        ax[1].loglog(1/k, (np.abs(np.fft.rfft(seg))**2)[1:], color=c, lw=1.2, label=name)
    ax[1].set_xlabel("wavelength [d_e]"); ax[1].set_ylabel("power |Bx_k|^2"); ax[1].legend(fontsize=8)
    ax[1].set_title("far-upstream spectrum by variant")
    fig.tight_layout(); os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=115); print("wrote", args.out)

if __name__ == "__main__":
    main()
