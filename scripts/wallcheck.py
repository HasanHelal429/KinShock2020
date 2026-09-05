#!/usr/bin/env python3
"""Go/no-go check for a reflecting-wall (W_phase) run.

Three questions, in the order they can fail:

  1. DID THE DRIFT LAND?  The ambient must arrive at the wall at -v_sh. The first
     attempt at this smoke test looked perfect -- COMPLETED, zero warnings -- while
     the plasma sat still, because `uz_m` is silently ignored under
     momentum_distribution_type = maxwellian (the key is `uz_mean`, alongside
     maxwellian_u_mean_distribution_type). WarpX's only tell was a line in the
     "Unused ParmParse Variables" block. So this is checked FIRST and against the
     data, not the deck.
  2. DOES IT REFLECT?  Fraction of ambient ions with v_z > 0. Starts at ~0 and must
     grow. If it stays at 0 the boundary is transmitting, not reflecting.
  3. DOES A RAMP FORM AND LEAVE?  Peak n/n_0 and peak |Bx|/B0 must both exceed 1 and
     their z must WALK AWAY from the wall. A pile that sits at z=0 is a boundary
     artefact; a shock propagates.

|Bx|/B0 is smoothed with the 0.2 d_i0 boxcar RESULTS 2026-08-11 requires -- an
unsmoothed peak runs 2.2-2.8x high and has already cost one retracted result.

    python scripts/wallcheck.py runs/W_phase/ws_smoke_470eV
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)
import kinshock                                   # noqa: E402
from kinshock import io as kio                    # noqa: E402
from kinshock.units import C, ME                  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    args = ap.parse_args()

    import yt; yt.set_log_level(50)
    cfg = kinshock.load(args.run)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    want = cfg["species"]["amb_ions"].get("drift_uz_over_c")
    pfs = kio.plotfiles(args.run)
    if not pfs:
        print("NO PLOTFILES -- run produced no output"); return 2

    # the deck's own confession, if it made one
    log = os.path.join(args.run, "run.log")
    if os.path.isfile(log):
        txt = open(log, errors="ignore").read()
        blk = txt.split("Unused ParmParse Variables")
        if len(blk) > 1:
            bad = [l for l in blk[1].splitlines()[:60] if "uz_m" in l or "u_mean" in l]
            if bad:
                print("!! WarpX IGNORED drift keys:")
                for l in bad:
                    print("    " + l.strip())

    print(f"{'t*':>7} {'<vz>/vsh':>9} {'refl':>7} {'n/n0':>7} {'z_n':>7} "
          f"{'|Bx|/B0':>8} {'z_B':>7}")
    rows = []
    # n_0 is fixed at t=0 and reused. A wall run is a FINITE SLAB drifting into the wall
    # with no injector behind it, so the far field EMPTIES: measured here, the trailing
    # edge walks 12.02 -> 10.71 d_i0 by t* = 0.12 and the last-40-bin median hits exactly
    # zero, which sent n/n_0 to 1e19. Total weight is conserved (1.4e21 throughout), so
    # that was the normalisation, not lost particles. Normalising against the instantaneous
    # far field is only valid where an injector holds it up.
    n0_ref = None
    for p in pfs:
        fr = kio.load_frame(p, fields=("Bx",))
        ad = fr.ds.all_data()
        z = np.asarray(ad[("amb_ions", "particle_position_x")]) / sc.di0
        pm = np.asarray(ad[("amb_ions", "particle_momentum_z")])
        g = np.sqrt(1 + (pm / (mi * C)) ** 2)
        v = pm / (mi * g)
        w = np.asarray(ad[("amb_ions", "particle_weight")])
        vbar = np.average(v, weights=w) / sc.vsh_model
        refl = float(np.sum(w[v > 0]) / np.sum(w))
        zb = np.linspace(0, sc.domain_halfwidth / sc.di0, 321)
        n, _ = np.histogram(z, bins=zb, weights=w)
        if n0_ref is None:
            n0_ref = float(np.median(n[-40:]))
        nn = n / n0_ref if n0_ref > 0 else n
        zc = 0.5 * (zb[1:] + zb[:-1])
        # search only where the slab still is -- argmax over emptied bins picks noise
        occ = n > 0
        i = int(np.argmax(np.where(occ, nn, -np.inf)))
        Bx = fr.comps.get("Bx")
        if Bx is not None:
            zf = fr.z_centers / sc.di0
            kw = max(int(round(0.2 / float(np.median(np.diff(zf))))), 1)
            # smooth SIGNED then |.|; abs-first rectifies noise into a DC offset big
            # enough (10.9 on mb_r100) to sail past the ramp test as a false PASS
            b = np.abs(np.convolve(Bx / sc.B0, np.ones(kw) / kw, mode="same"))
            j = int(np.argmax(b)); bmax, zbm = b[j], zf[j]
        else:
            bmax, zbm = float("nan"), float("nan")
        t = fr.time * sc.wci0
        print(f"{t:7.4f} {vbar:9.4f} {refl:7.4f} {nn[i]:7.3f} {zc[i]:7.3f} "
              f"{bmax:8.3f} {zbm:7.3f}")
        rows.append((t, vbar, refl, nn[i], zc[i], bmax, zbm))

    a = np.array(rows)
    print()
    ok = True
    if want is not None:
        got = a[0, 1] * sc.vsh_model / C
        err = abs(got - want) / abs(want)
        good = err < 0.05
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] drift: want {want:+.6f} c, "
              f"got {got:+.6f} c ({err*100:.2f}% off)")
    grew = a[-1, 2] > max(3 * a[0, 2], 0.01)
    ok &= grew
    print(f"[{'PASS' if grew else 'FAIL'}] reflection: {a[0,2]:.4f} -> {a[-1,2]:.4f}")
    ramp = a[:, 3].max() > 1.3 and (np.nanmax(a[:, 5]) > 1.3)
    ok &= ramp
    print(f"[{'PASS' if ramp else 'FAIL'}] ramp: peak n/n0 {a[:,3].max():.2f}, "
          f"peak |Bx|/B0 {np.nanmax(a[:,5]):.2f}")
    # Only frames that HAVE a ramp can say whether the ramp moves. Before one forms, the
    # density argmax is noise (z_n jumps 0.02 -> 0.66 -> 0.43 in consecutive early frames),
    # and averaging that into a first-half mean buried a clean monotonic 1.45 -> 1.82 walk.
    # Require a positive least-squares slope of ramp position vs t* over formed frames.
    formed = a[:, 3] > 1.3
    if formed.sum() >= 4:
        sl_n = float(np.polyfit(a[formed, 0], a[formed, 4], 1)[0])
        sl_B = float(np.polyfit(a[formed, 0], a[formed, 6], 1)[0])
        moved = sl_n > 0 and sl_B > 0
        print(f"[{'PASS' if moved else 'FAIL'}] ramp moves off wall over {int(formed.sum())} "
              f"formed frames: dz_n/dt* = {sl_n:+.2f}, dz_B/dt* = {sl_B:+.2f} d_i0 per t*")
    else:
        moved = False
        print(f"[FAIL] ramp moves off wall: only {int(formed.sum())} frames with n/n0 > 1.3")
    ok &= moved
    print(f"\n{'ALL CHECKS PASS' if ok else 'CHECKS FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
