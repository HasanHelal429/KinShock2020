"""Reader for psc_schaeffer2020's binned phase-space diagnostic (ps.<step>.bin).

    from read_ps import load_header, load_frame
    hdr = load_header(run_dir)
    fr = load_frame(run_dir, step, hdr)
    fr["hist"][kind][axis]  # (nz_bins, nu_bins), axis 0 = (z, u_z), 1 = (z, u_x)
    fr["n"][kind], fr["T"][kind]   # density and temperature profiles along z

Density per z bin = sum w / (nicell * ny * cells_per_zbin), the particle-dump
normalisation. T = <u^2> - <u>^2 averaged over the three axes (m = 1 for e; times
the mass ratio for ions gives T in m_e c^2), non-relativistic moment.
"""
import os
import numpy as np


def load_header(run):
    hdr = {"kinds": []}
    with open(os.path.join(run, "ps_header.txt"), errors="replace") as fh:
        for line in fh:
            w = line.split()
            if not w or w[0] in ("file", "mom"):
                continue
            if w[0] == "kind":
                i = w.index("ulo")
                hdr["kinds"].append(dict(name=" ".join(w[2:i]), ulo=float(w[i + 1]),
                                         uhi=float(w[w.index("uhi") + 1])))
            else:
                hdr[w[0]] = float(w[1]) if "." in w[1] or "e" in w[1] else int(w[1])
    return hdr


def load_frame(run, step, hdr, masses=None):
    nk, nz, nu = len(hdr["kinds"]), hdr["nz_bins"], hdr["nu_bins"]
    raw = np.fromfile(os.path.join(run, f"ps.{step:09d}.bin"), dtype=np.float32)
    nh = nk * 2 * nz * nu
    hist = raw[:nh].reshape(nk, 2, nz, nu)
    mom = raw[nh:].reshape(nk, 7, nz).astype(float)
    norm = hdr["nicell"] * hdr["ny"] * hdr["cells_per_zbin"]
    w = mom[:, 0]
    ws = np.where(w > 0, w, np.nan)
    mean = mom[:, 1:4] / ws[:, None]
    var = mom[:, 4:7] / ws[:, None] - mean ** 2
    masses = masses or [1.0] * nk
    T = var.mean(axis=1) * np.array(masses)[:, None]
    z = hdr["zlo"] + (np.arange(nz) + 0.5) * (hdr["zhi"] - hdr["zlo"]) / nz
    return dict(z=z, hist=hist / norm, n=w / norm, u_mean=mean, T=T)
