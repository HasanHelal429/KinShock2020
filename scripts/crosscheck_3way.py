#!/usr/bin/env python3
"""3-way half-domain wall validation (RESULTS.md 2026-07-23).

Compares the shock diagnostics of three runs on z >= 0:
  full  = runs/R1_core           (symmetric domain, the reference)
  spec  = runs/R1_core_half      (one-sided, SPECULAR wall -- old, artifact)
  sym   = runs/R1_core_half_sym  (one-sided, pi-rotation SYMMETRY wall -- the fix)

Question: does the symmetry wall move the half-domain shock kinematics
(front speed, compression, onset) closer to the full-domain reference than the
specular wall did? Prints a table; the fix "passes" if sym is closer to full
than spec on front speed / field compression / onset, and near-wall B_perp drops.

Usage: python tmp/crosscheck_3way.py            # all three
       python tmp/crosscheck_3way.py full spec  # subset (e.g. before sym finishes)
"""
import os, sys, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
import kinshock
from kinshock import io, metrics
from kinshock.units import C

RUNS = {
    "full": "runs/R1_core",
    "spec": "runs/R1_core_half",
    "sym":  "runs/R1_core_half_sym",
}
T_CLEAN_MAX = 2.25          # exclude boundary-contaminated frames (shock hits edge ~2.5)
AMB_WIN_DE = (-400.0, -50.0)  # ambient shock window rel. to front [d_e]: [front-400, front-50]


def analyze(tag, run_dir):
    run_dir = os.path.join(ROOT, run_dir)
    cfg = kinshock.load(run_dir)
    sc = kinshock.units.derive(cfg)
    slab_de = float(cfg["geometry"]["slab_halfwidth_di"]) * sc.di / sc.de  # piston half-width [d_e]
    near_wall_de = 3.0 * slab_de
    edge_de = float(cfg["geometry"]["domain_halfwidth_de"])
    # piston exclusion for front tracking: MATCH make_figures.fig_trajectory exactly
    # (slab_halfwidth_di * di) so both diagnostics track the identical front -> identical v_sh
    zexcl = float(cfg["geometry"]["slab_halfwidth_di"]) * sc.di

    pfs = io.plotfiles(run_dir)
    frames = [io.load_frame(p) for p in pfs]

    ts, zf_de, ncomp, bcomp, nearwall = [], [], [], [], []
    for fr in frames:
        t = fr.time
        zc_de = np.asarray(fr.z_centers) / sc.de
        n_all = io.species_density(fr, cfg["ion_species"])          # all ions -> front tracking
        n_amb = io.species_density(fr, cfg["ambient_ion_species"])  # ambient only -> compression
        bperp = fr.Bperp
        # +z shock front (total-ion density threshold), piston excluded
        zf = metrics.track_front(fr.z_centers, n_all, sc.namb, threshold=1.5,
                                 z_exclude=zexcl, side=+1)
        ts.append(t); zf_de.append(zf / sc.de if np.isfinite(zf) else np.nan)
        # ambient compression in a window behind the front, EXCLUDING the piston/near-wall
        # zone (z > 3*slab) so piston ions / field pileup do not contaminate the ambient value
        if np.isfinite(zf):
            f_de = zf / sc.de
            w = ((zc_de >= f_de + AMB_WIN_DE[0]) & (zc_de <= f_de + AMB_WIN_DE[1])
                 & (zc_de > near_wall_de))
            ncomp.append(np.nanmax(n_amb[w]) / sc.namb if w.any() else np.nan)
            bcomp.append(np.nanmax(bperp[w]) / sc.B0 if w.any() else np.nan)
        else:
            ncomp.append(np.nan); bcomp.append(np.nan)
        # near-wall B_perp (z in (0, 3*slab))
        nw = (zc_de > 0) & (zc_de < near_wall_de)
        nearwall.append(np.nanmax(bperp[nw]) / sc.B0 if nw.any() else np.nan)

    ts = np.array(ts); zf_de = np.array(zf_de)
    tw = ts * sc.wci0
    ncomp = np.array(ncomp); bcomp = np.array(bcomp); nearwall = np.array(nearwall)

    # front speed: SHARED metrics.speed_from_trajectory with the identical
    # domain-aware clean window make_figures uses (front < 0.94*edge, no second-half)
    # -> the reflected-ion threshold is now consistent across both diagnostics.
    vfit = metrics.speed_from_trajectory(ts, zf_de * sc.de,
                                         z_edge=edge_de * sc.de, use_second_half=False)  # m/s
    # clean mask retained only for the compression-averaging window below
    clean = np.isfinite(zf_de) & (tw <= T_CLEAN_MAX) & (zf_de < 0.94 * edge_de)

    # onset t*_1 (max dG/dt) and z*_1, using the +z ambient ions only (the full
    # symmetric run has ambient ions on both sides -> restrict to z>=0 so all three
    # runs compare the same +z shock)
    vref = vfit if np.isfinite(vfit) else sc.vsh_model

    def amb_phase_zpos(fr):
        z, uz = io.species_phase(fr, "amb_ions", sc, mass=sc.mi)
        m = np.asarray(z) >= 0.0
        return np.asarray(z)[m], np.asarray(uz)[m]

    G = []
    for fr in frames:
        z, uz = amb_phase_zpos(fr)
        G.append(metrics.reflected_fraction_G(uz * C, vref))
    G = np.array(G)
    # Onset t*_1: SHARED metrics.onset_time_from_G (first prominent dG/dt peak) so
    # this matches make_figures and is stable against the bimodal-G argmax flip.
    tstar, i1 = metrics.onset_time_from_G(ts, G)
    zstar = np.nan
    if i1 >= 0:
        fr = frames[i1]
        z, uz = amb_phase_zpos(fr)
        edges = np.asarray(fr.z_edges)
        edges = edges[edges >= 0.0]
        F, centers = metrics.reflected_profile_F(z, uz * C, vref, edges)
        ker = np.ones(50) / 50
        for _ in range(6):
            F = np.convolve(F, ker, mode="same")
        zc = centers / sc.rho_i0
        zstar, _ = metrics.onset_location_from_F(zc, F)

    # clean-window averaged compressions
    cmask = clean & (tw >= 0.5)
    return {
        "tag": tag, "n_frames": len(frames), "n_clean": int(clean.sum()),
        "v_sh_c": vfit / C, "v_sh_Cs": vfit / sc.Cs_ab, "v_sh_vp": vfit / sc.vp_model,
        "n_comp": float(np.nanmean(ncomp[cmask])) if cmask.any() else np.nan,
        "b_comp": float(np.nanmean(bcomp[cmask])) if cmask.any() else np.nan,
        "near_wall_bperp": float(np.nanmax(nearwall[clean])) if clean.any() else np.nan,
        "tstar_wci0": tstar * sc.wci0, "zstar_rhoi0": zstar,
        # arrays for plotting (front trajectory in rho_i0 vs t*wci0, and G(t))
        "traj_tw": (tw).tolist(),
        "traj_front_rho": (zf_de * sc.de / sc.rho_i0).tolist(),
        "G_tw": (tw).tolist(), "G": G.tolist(),
        "edge_rho": edge_de * sc.de / sc.rho_i0, "t_clean_max": T_CLEAN_MAX,
    }


def main():
    which = sys.argv[1:] or list(RUNS)
    res = {}
    for tag in which:
        if not os.path.isdir(os.path.join(ROOT, RUNS[tag], "diags")):
            print(f"[skip] {tag}: {RUNS[tag]}/diags not found (run not done?)")
            continue
        print(f"[load] {tag} <- {RUNS[tag]} ...", flush=True)
        res[tag] = analyze(tag, RUNS[tag])

    cols = [("v_sh (c)", "v_sh_c", "{:.4f}"), ("v_sh (C_s,ab)", "v_sh_Cs", "{:.2f}"),
            ("v_sh / v_p", "v_sh_vp", "{:.2f}"), ("n_comp (amb)", "n_comp", "{:.2f}"),
            ("B_comp (amb)", "b_comp", "{:.2f}"), ("near-wall B/B0", "near_wall_bperp", "{:.1f}"),
            ("t*_1 (wci0^-1)", "tstar_wci0", "{:.2f}"), ("z*_1 (rho_i0)", "zstar_rhoi0", "{:.2f}"),
            ("n_clean", "n_clean", "{:d}")]
    print("\n=== 3-way half-domain wall validation (z >= 0) ===")
    hdr = f"{'quantity':<18}" + "".join(f"{t:>12}" for t in which)
    print(hdr); print("-" * len(hdr))
    for label, key, fmt in cols:
        row = f"{label:<18}"
        for tag in which:
            v = res.get(tag, {}).get(key, float('nan'))
            row += f"{(fmt.format(v) if v == v else 'n/a'):>12}"
        print(row)
    print("\nPaper: v_sh=4.6 C_s,ab=0.138c; v_sh/v_p~1.33; t*_1~1 rho_i0; z*_1~1.")
    if {"full", "spec", "sym"} <= set(res):
        print("\n--- did the symmetry wall move toward full? (|sym-full| < |spec-full| = better) ---")
        for label, key, _ in cols[:-1]:
            f, s, y = res["full"][key], res["spec"][key], res["sym"][key]
            if not (f == f and s == s and y == y):
                continue
            d_spec, d_sym = abs(s - f), abs(y - f)
            verdict = "BETTER" if d_sym < d_spec else ("worse" if d_sym > d_spec else "same")
            print(f"  {label:<18} full={f:8.3f}  spec Δ={s-f:+7.3f}  sym Δ={y-f:+7.3f}  -> {verdict}")
    out = os.path.join(ROOT, "media", "testing", "crosscheck_3way.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, "w"), indent=2)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
