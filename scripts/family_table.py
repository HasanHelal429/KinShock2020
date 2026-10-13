#!/usr/bin/env python3
"""Summary table of the configuration families (one representative run each), B0 and T_e first.

    python scripts/family_table.py            # markdown to stdout
    python scripts/family_table.py --csv F    # also write a CSV (one column per family)

Every number is derived from the run's config.yaml through kinshock.units.derive, so the table
cannot drift from what the decks ran. Columns follow the shock workbook's families:
47 keV = "Inflated Te" (reduce m_i), 470 eV = "Reduced ME", Inflated m_e (PSC's convention at
real c), RMR = real mass ratio, and the V_phase speed-corner runs (deflated T_e, deflated B).
The false-shock row is the user's by-eye verdict (RESULTS 2026-10-12).
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import kinshock  # noqa: E402
from kinshock import units  # noqa: E402

# (column label, run dir under runs/, how mu is reached, false-shock verdict)
FAMILIES = [
    ("47 keV", "E_phase/ek_47keV_t2", "light ion", "clean"),
    ("470 eV", "R1_phase/R1_470eV_s3f8_t2", "light ion", "notch (after 0.3)"),
    ("Inflated m_e 470 eV", "IM_phase/im_470eV_t03", "heavy electron", "notch (by 0.3)"),
    ("Inflated m_e 47 keV", "IM_phase/im_47keV_t03", "heavy electron", "clean"),
    ("RMR 470 eV", "U_phase/um1836_470eV", "real (1836)", "notch (by 0.3)"),
    ("RMR 47 keV", "U_phase/um1836_47keV", "real (1836)", "clean"),
    ("Deflated T_e", "V_phase/vA_lowTab", "light ion", "notch (by 0.3)"),
    ("Deflated B", "V_phase/vB_lowB0", "light ion", "no movie"),
    ("Scaled (vC)", "V_phase/vC_scaled", "light ion", "not run"),
]


def _num(x):
    """Plain number: thousands separators at >= 1000, 3 significant figures below."""
    return f"{x:,.0f}" if abs(x) >= 1000 else f"{x:.3g}"


def quantities(run_rel):
    cfg = kinshock.load(os.path.join(ROOT, "runs", run_rel))
    sc = units.derive(cfg)
    C, QE, EPS0 = units.C, units.QE, units.EPS0
    T0 = cfg["plasma"]["ambient"]["theta_0"] * sc.me * C ** 2 / QE
    vte0 = math.sqrt(T0 * QE / sc.me)
    lamD0 = math.sqrt(EPS0 * T0 / (sc.namb * QE))
    wpe0 = math.sqrt(sc.namb * QE ** 2 / (EPS0 * sc.me))
    wce0 = QE * sc.B0 / sc.me
    ppc = cfg["numerics"]["ppc"]["ambient"]
    t_end = int(cfg["numerics"]["max_step"]) / sc.steps_per_wci0
    return [
        ("B0 [T]", f"{sc.B0:.3g}"),
        ("T_e,ab [eV]", _num(sc.Te_ab_eV)),
        ("T_0 upstream [eV]", _num(T0)),
        ("m_i/m_e (mu)", f"{sc.mi / sc.me:.0f}"),
        ("m_e / m_e,real", f"{sc.me / units.ME:.4g}"),
        ("n_e,ab [m^-3]", f"{sc.n0:.2e}"),
        ("n_e0 upstream [m^-3]", f"{sc.namb:.2e}"),
        ("eps_ab = v_te,ab/c", f"{sc.vte_ab / C:.4f}"),
        ("eps_0 = v_te,0/c", f"{vte0 / C:.5f}"),
        ("v_sh (model) [km/s]", _num(sc.vsh_model / 1e3)),
        ("v_A [km/s]", _num(sc.vA / 1e3)),
        ("M_A", f"{sc.MA:.2f}"),
        ("M_ms", f"{sc.Mms:.2f}"),
        ("beta_ab", _num(sc.beta_ab)),
        ("beta_0", f"{sc.beta_0:.3g}"),
        ("w_pe0 / w_ce0", f"{wpe0 / wce0:.3g}"),
        ("d_i0 [um]", _num(sc.di0 * 1e6)),
        ("1/w_ci0 [ns]", f"{1e9 / sc.wci0:.4g}"),
        ("dz / lambda_D,0", f"{sc.dz / lamD0:.3f}"),
        ("N_D = ppc lambda_D/dz", f"{ppc * lamD0 / sc.dz:.0f}"),
        ("runs to t*w_ci0", f"{t_end:.2f}"),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=None)
    args = ap.parse_args()

    cols = [(label, run, how, verdict, quantities(run)) for label, run, how, verdict in FAMILIES]
    names = [k for k, _ in cols[0][4]]
    rows = [["run"] + [f"`{run.split('/')[1]}`" for _, run, _, _, _ in cols],
            ["mu reached by"] + [how for _, _, how, _, _ in cols]]
    rows += [[n] + [dict(q)[n] for *_, q in cols] for n in names]
    rows += [["false shock (by eye)"] + [v for _, _, _, v, _ in cols]]
    header = ["quantity"] + [label for label, *_ in cols]
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for r in rows:
        print("| " + " | ".join(r) + " |")
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            for r in rows:
                w.writerow([c.strip("`") for c in r])
        print(f"\nwrote {args.csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
