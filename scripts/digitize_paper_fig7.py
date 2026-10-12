#!/usr/bin/env python3
"""Digitize Schaeffer 2020 Fig. 7 and compare it against `ek_47keV_t2`.

WHY THIS IS THE PSC CROSS-CHECK, WITHOUT PSC. Everything in E_phase and M_phase is WarpX
validating WarpX. The one question none of it answers is whether the eps dependence is
physics or a WarpX artifact -- that needs a second code. PSC is not available here (no
tree on Perlmutter; the sibling laser-deposition project records PSC access as blocked),
but the paper's Fig. 7 IS a PSC run at the point we care about: its caption says "for the
run in Table I", i.e. eps = 0.3033, w_pe0/w_ce = 10, M_A = 14, beta_0 = 0.2, and it plots
ambient-ion phase space in z/d_i0 vs v_z/v_sh at t*w_ci0 = 0.15 - 1.25. `ek_47keV_t2` is
the same dimensionless point run to t* = 2.0, so every panel has a partner.

METHOD. Panel frames and tick positions are found from the image (dark rows/columns and
inward tick stubs); tick VALUES are transcribed in PANELS below from the printed labels.
Data pixels are the CHROMATIC ones -- the colormap is white -> blue -> yellow, so
max(RGB)-min(RGB) > 25 selects data and rejects the white ground, the grey shaded bands
and the black frame and text in one test. Our particles are then rendered onto the SAME
pixel grid with the SAME extent rule, so the two sides are measured identically.

WHAT IS AND IS NOT COMPARABLE. Two extents survive the round trip:
  z_piston   rightmost column carrying piston-ion signal
  z_disturb  rightmost column with ambient signal ABOVE v_z/v_sh = 0.2, i.e. clear of the
             undisturbed band at v_z ~ 0
The repo's actual wedge metric -- first offset ahead of the piston where the local v_z
STANDARD DEVIATION falls below 0.06 v_sh -- does NOT survive it. A printed contour cannot
show a low-amplitude spread enhancement, so the wedge cannot be read off the figure at all
and is not attempted here.

THE THRESHOLD IS THE WHOLE BALL GAME, AND IT IS NOT OPTIONAL. Our particles have no
rendering floor; the printed figure has one. Two earlier attempts at this got opposite answers from the same data purely through the
threshold, so it is matched by inked AREA FRACTION against each printed panel and
`--ink-scale` exists to report how far the answer moves when that match is varied.

MEASURED SENSITIVITY (--ink-scale 0.6 / 1.0 / 1.6):
  z_piston  is INSENSITIVE -- identical to 2 dp at every scale.
  z_disturb is insensitive at t* >= 0.98 (W/P = 1.05 throughout) but swings 1.02-1.27 at
            t* = 0.49, so nothing may be concluded from the early panels.
Quote nothing from this script without the sensitivity beside it.

RESULT. At t* = 0.98 and 1.25 WarpX reproduces PSC to 1% on the piston front and 5% on the
ambient disturbance extent. Early panels differ by 13-38%, where the structures are a few
d_i0 across and the threshold dominates. So WarpX is not grossly wrong at eps = 0.3033 --
which is what the eps ladder needed, since every other check in E_phase and M_phase is
WarpX against itself. It does NOT validate eps = 0.0303, and it does not validate the
wedge metric, which is unreadable from a printed contour.

    python scripts/digitize_paper_fig7.py [--run runs/E_phase/ek_47keV_t2]

Needs `schaeffer2020.pdf` in the repo root (gitignored; OSTI 1643946) or the pre-extracted
`media/p7_0_Im1.png`.
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
os.chdir(ROOT)

# (x0, x1, major tick x-positions, their printed values, t*w_ci0)
PANELS = [
    (185,  758,  [229, 326, 423, 520, 617, 714], [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], 0.15),
    (766,  1340, [794, 924, 1053, 1182, 1311],   [2, 4, 6, 8, 10],               0.49),
    (1348, 1922, [1370, 1476, 1581, 1687, 1793, 1899], [2, 4, 6, 8, 10, 12],     0.73),
    (1930, 2503, [1950, 2047, 2143, 2240, 2337, 2434], [4, 6, 8, 10, 12, 14],    0.98),
    (2511, 3085, [2527, 2607, 2687, 2768, 2848, 2928, 3009],
     [6, 8, 10, 12, 14, 16, 18], 1.25),
]
R1, R2 = (70, 449), (457, 837)      # ambient-ion and piston-ion row frames [px]
Y1_0, Y1_S = 388.0, 160.5           # row 1: v_z = 0 at y = 388; 160.5 px per v_z/v_sh
Y2_0, Y2_S = 761.0, 170.0           # row 2
VCUT = 0.2                          # v_z/v_sh above which ambient counts as disturbed


def _vrange(row, y0, ys):
    return ((row[1] - y0) / -ys, (row[0] - y0) / -ys)


def _rightmost(mask, ze):
    xs = np.where(mask.sum(0) >= 2)[0]
    return ze[xs.max()] if len(xs) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", default="runs/E_phase/ek_47keV_t2")
    ap.add_argument("--img", default="media/p7_0_Im1.png")
    ap.add_argument("--ink-scale", type=float, default=1.0,
                    help="multiply the matched ink fraction, to test sensitivity")
    args = ap.parse_args()

    from PIL import Image
    import kinshock
    from kinshock import io as kio
    from kinshock.units import C, ME
    import yt; yt.set_log_level(50)

    a = np.asarray(Image.open(args.img).convert("RGB")).astype(int)
    SIG = (a.max(2) - a.min(2)) > 25

    cfg = kinshock.load(args.run)
    sc = kinshock.units.derive(cfg)
    mi = cfg["reference"]["mass_ratio"] * ME
    pfs = kio.plotfiles(args.run)
    times = np.array([kio.load_frame(p).time * sc.wci0 for p in pfs])

    print(f"\nSchaeffer 2020 Fig. 7 (PSC, Table I run)  vs  {os.path.basename(args.run)} "
          "(WarpX, same dimensionless point)")
    print("matched grid, identical extent rule, all in d_i0\n")
    print(f"{'t*':>6} | {'z_piston':^22} | {'z_ambient_disturb':^22}")
    print(f"{'':>6} | {'PSC':>7}{'WarpX':>8}{'W/P':>7} | {'PSC':>7}{'WarpX':>8}{'W/P':>7}")
    for x0, x1, tx, tv, ts in PANELS:
        A = np.polyfit(tx, tv, 1)
        zlo, zhi = A[0] * x0 + A[1], A[0] * x1 + A[1]
        W = x1 - x0
        ze = np.linspace(zlo, zhi, W)
        ycut = int(Y1_0 - VCUT * Y1_S) - R1[0]
        mA, mP = SIG[R1[0]:R1[1], x0:x1], SIG[R2[0]:R2[1], x0:x1]
        zpP, zdP = _rightmost(mP, ze), _rightmost(mA[:ycut, :], ze)

        ad = kio.load_frame(pfs[int(np.argmin(np.abs(times - ts)))]).ds.all_data()

        def phase(sp):
            z = np.asarray(ad[(sp, "particle_position_x")]) / sc.di0
            p = np.asarray(ad[(sp, "particle_momentum_z")])
            g = np.sqrt(1 + (p / (mi * C)) ** 2)
            return z, p / (mi * g) / sc.vsh_model

        def rendered(z, v, row, y0, ys, ink):
            """Render to the panel's grid and threshold so the INKED AREA FRACTION equals
            the printed panel's. Without this the comparison is not like-for-like: our
            particles have no rendering floor, so `H > 0` counts single particles the
            figure's lowest contour would never show, and the extents come out long."""
            vlo, vhi = _vrange(row, y0, ys)
            H, _, _ = np.histogram2d(
                z, v, bins=[np.linspace(zlo, zhi, W + 1),
                            np.linspace(vlo, vhi, row[1] - row[0] + 1)])
            H = H.T[::-1]
            return H > np.quantile(H, 1.0 - ink)      # area match, over ALL bins

        za, va = phase("amb_ions")
        zp_, vp_ = phase("piston_ions")
        zpW = _rightmost(rendered(zp_, vp_, R2, Y2_0, Y2_S, min(mP.mean()*args.ink_scale, 1.0)), ze)
        zdW = _rightmost(rendered(za, va, R1, Y1_0, Y1_S, min(mA.mean()*args.ink_scale, 1.0))[:ycut, :], ze)
        print(f"{ts:>6.2f} | {zpP:>7.2f}{zpW:>8.2f}{zpW / zpP:>7.2f} |"
              f" {zdP:>7.2f}{zdW:>8.2f}{zdW / zdP:>7.2f}")
    print("\nW/P < 1 means WarpX is SHORT of PSC. Read the header on why a short WarpX")
    print("extent is the stronger result: the rule biases our side long.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
