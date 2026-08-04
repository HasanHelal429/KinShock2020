#!/usr/bin/env python
"""Schaeffer 2020 Table I in three unit systems: PSC code units, physical, WarpX.

WHY THIS SCRIPT EXISTS
----------------------
A PIC run of this problem is a set of DIMENSIONLESS numbers. It corresponds to a
whole family of real plasmas, and you pick one member of that family by choosing
exactly three free parameters:

    1. the ablation density        n_e,ab   (real 6e20 cm^-3; code value ~1)
    2. the ablation temperature    T_e,ab   (real 470 eV; code value ~0.1, which
                                             IS the reduced speed of light --
                                             Sec. II: c = sqrt(mu_p/T_e,ab) C_s,ab)
    3. the collisionality          lambda_ab (= mfp_e,ab/d_e,ab = 20)

Everything else in Table I follows. The order of operations (and this script's
structure) is:

    step 1  pick the three REAL values
    step 2  pick the corresponding CODE values (density O(1), temperature O(0.1)
            so the run stays non-relativistic)
    step 3  derive the rest of the CODE column -- all of it follows from the three
    step 4  derive the rest of the PHYSICAL column; beta_ab then sets B0

THE ONE TRAP: THE MASS RATIO
----------------------------
The code runs at mu_p = m_i/m_e = 100; the real hydrogen plasma is at 1836. A
dimensionless number transfers between the columns ONLY if it is insensitive to
that. Sec. II says exactly which ones are not (p. 3):

    "This scaling ensures that dimensionless quantities such as the magnetic
     Reynolds number are correct, but electron collisionality relative to global
     scales (e.g. nu_ei,ab t_ab) is only quantitatively matched at physical mass
     ratios."

So `lambda_ab` transfers (it is a purely electron-scale ratio, mfp/d_e) but
`nu_ei,ab * t_ab` does NOT: it equals mu/lambda_ab, hence 5.0 in the code column
and 91.8 in the physical one -- a factor mu_phys/mu_p = 18.4. Deriving the
physical Coulomb logarithm therefore has to go through lambda_ab, never through
nu_ei t_ab, or the answer comes out 18.4x low (0.5 instead of ~9).

Run with --show-work to print the algebra behind the Coulomb logarithm.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kinshock import config as kcfg                     # noqa: E402
from kinshock.units import (C, EPS0, ME, ME_C2_EV, MP, MU0, NU_EI_NRL,  # noqa: E402
                            QE)

MU_PHYS = MP / ME          # 1836.15, the real proton-to-electron mass ratio

# Table I's own numbers, transcribed from schaeffer2020.pdf p.4, for the side-by-side
# check. Two significant figures throughout, so ~10% agreement is a match.
TABLE1 = {
    "Lx":            ("0.5 d_i,ab",     "4.7 um"),
    "Lz":            ("900 d_i,ab",     "8.4 mm"),
    "tau_sim":       ("220 t_ab",       "10.9 ns"),
    "Z_ab":          ("1",              "1"),
    "n_e,ab":        ("1.25",           "6e20 cm^-3"),
    "T_e,ab":        ("0.092 m_e c^2",  "470 eV"),
    "tau_ei,ab":     ("0.009 t_ab",     "0.43 ps"),
    "C_s,ab":        ("0.030 c",        "210 km/s"),
    "v_p":           ("0.104 c",        "730 km/s"),
    "v_sh":          ("4.6 C_s,ab",     "980 km/s"),
    "B_0":           ("0.01 sqrt(m_e c^2)", "7 T"),
    "Z_0":           ("1",              "1"),
    "n_e0":          ("0.01 (code)",    "4.8e18 cm^-3"),
    "T_0":           ("0.002 m_e c^2",  "10 eV"),
    "d_i0":          ("11.2 d_i,ab",    "104 um"),
    "1/w_ci0":       ("33.9 t_ab",      "1.5 ns"),
    "m_i/m_e":       ("100",            ""),
    "c_sim/c_phys":  ("0.02",           ""),
    "beta_ab":       ("1150",           ""),
    "beta_0":        ("0.2",            ""),
    "lambda_ab":     ("20",             ""),
    "mfp/d_i0":      ("350",            ""),
    "M_A":           ("14",             ""),
    "M_ms":          ("13",             ""),
}


def wpe(n):
    """Electron plasma frequency [rad/s] at density n [m^-3]. Real e, m_e, eps0."""
    return math.sqrt(n * QE * QE / (EPS0 * ME))


def nu_ei_formulary(n_cm3, Te_eV, lnL, coeff):
    """nu_ei = coeff * n[cm^-3] * lnLambda * T[eV]^{-3/2}   [s^-1]."""
    return coeff * n_cm3 * lnL * Te_eV ** -1.5


def lnL_nrl(n_cm3, Te_eV):
    """NRL electron-ion Coulomb logarithm, T_e > 10 eV: 24 - ln(sqrt(n)/T)."""
    return 24.0 - math.log(math.sqrt(n_cm3) / Te_eV)


def mfp_directed(v, n, mi, lnL, Z=1.0):
    """Rutherford momentum-transfer mean free path [m] of an ion of speed ``v``
    streaming through density ``n`` [m^-3]:

        lambda = 4 pi eps0^2 mi^2 v^4 / (n Z^4 e^4 lnLambda)

    This is the DIRECTED-ion mfp, ~v^4, and it is what Table I's
    ``lambda_mfp/d_i0 = 350`` row reports -- see the note in main(). The thermal
    ion-ion mfp at T_0 = 10 eV is smaller by (v_ti/v_sh)^4 ~ 1e-6 and is NOT
    what that row means.
    """
    return (4.0 * math.pi * EPS0 ** 2 * mi ** 2 * v ** 4
            / (n * Z ** 4 * QE ** 4 * lnL))


# --------------------------------------------------------------------------- #
# step 2 + 3 : the CODE column
# --------------------------------------------------------------------------- #
def code_column(cfg, lambda_ab):
    """Everything PSC actually integrates, in PSC's normalization.

    Lengths in d_e,ab = c_sim/w_pe(n_ref), times in 1/w_pe(n_ref), velocities in
    c_sim, temperatures in m_e c_sim^2, densities in n_ref, B in m_e c_sim w_pe/e
    (so B_code = w_ce/w_pe).
    """
    ref, pis, amb = cfg["reference"], cfg["plasma"]["piston"], cfg["plasma"]["ambient"]
    mu_p = float(ref["mass_ratio"])
    th_e = float(pis["theta_e_heat"])
    th_0 = float(amb["theta_0"])

    # PSC's reference density is n_ref = n_e,ab/1.25, so n_e,ab = 1.25 code units and
    # the ambient is 0.01 code units. That 0.01 is what Table I's "0.01 n_e,ab" row
    # really means -- as a FRACTION of n_e,ab it is 0.008, which is what reproduces
    # both d_i0/d_i,ab = 11.2 and beta_0 = 0.2 exactly (see the checks below).
    n_ab_code = 1.25
    n_0_code = n_ab_code * float(amb["density_over_n0"])

    c = 1.0                                   # c_sim is the velocity unit
    Cs_ab = math.sqrt(th_e / mu_p) * c        # sqrt(Z T_e/m_i)
    vte_ab = math.sqrt(th_e) * c              # sqrt(T_e/m_e) -- the paper's mfp speed
    Cs_0 = math.sqrt(th_0 / mu_p) * c

    de_ab = 1.0                               # the length unit
    di_ab = de_ab * math.sqrt(mu_p)           # = 10 d_e,ab at mu_p = 100
    di_0 = di_ab * math.sqrt(n_ab_code / n_0_code)
    t_ab = di_ab / Cs_ab                      # ablation time, in 1/w_pe

    # Sec. II: lambda_ab = w_ce,ab/nu_ei,ab = mfp_ab/d_e,ab, with w_ce,ab at the
    # FUNDAMENTAL field B_ab = sqrt(mu0 n_e,ab T_e,ab) and mfp = sqrt(T_ab/m_e)/nu_ei.
    # The two forms are identical because rho_e(B_ab) = c/w_pe = d_e,ab exactly.
    mfp_ab = lambda_ab * de_ab
    nu_ei_ab = vte_ab / mfp_ab

    # B_code = w_ce/w_pe(n_ref). Derived from the config's own SI primary so this
    # column is provably the same run as the deck, not a re-typed Table I.
    n0_si = float(ref["n0"])
    B_code = (QE * float(cfg["field"]["B0_tesla"]) / ME) / wpe(n0_si / n_ab_code)

    # Identity (see derivation in the module docstring of scripts/README): for PSC's
    # B normalization, mu0 n kT/B^2 = theta * n_code / B_code^2. This is Table I's
    # beta convention -- NOT 2*mu0 n T/B^2. It reproduces 1150 and 0.2 exactly, and
    # is independently confirmed by Table I's own 1/w_ci0 = sqrt(beta_ab) t_ab row.
    beta_ab = th_e * n_ab_code / B_code ** 2
    beta_0 = th_0 * n_0_code / B_code ** 2

    # p.3: 1/w_ci0 = (Z_ab/Z_0) sqrt(beta_ab) t_ab
    wci0_inv = math.sqrt(beta_ab) * t_ab
    vA = B_code / math.sqrt(n_0_code * mu_p)      # code-unit Alfven speed

    vp = float(cfg["model"]["vp_over_c"]) * c
    vsh = float(cfg["model"]["vsh_over_Csab"]) * Cs_ab

    return dict(
        mu=mu_p, th_e=th_e, th_0=th_0, n_ab=n_ab_code, n_0=n_0_code,
        Cs_ab=Cs_ab, vte_ab=vte_ab, Cs_0=Cs_0, de_ab=de_ab, di_ab=di_ab, di_0=di_0,
        t_ab=t_ab, mfp_ab=mfp_ab, nu_ei_ab=nu_ei_ab, nu_t_ab=nu_ei_ab * t_ab,
        B=B_code, beta_ab=beta_ab, beta_0=beta_0, wci0_inv=wci0_inv, vA=vA,
        vp=vp, vsh=vsh, MA=vsh / vA, Mms=vsh / math.hypot(vA, Cs_0),
        lambda_ab=lambda_ab, c_over_Csab=c / Cs_ab,
    )


# --------------------------------------------------------------------------- #
# step 4 : the PHYSICAL column
# --------------------------------------------------------------------------- #
def phys_column(code, n_ab_cm3, Te_ab_eV, coeff):
    """The real HED plasma the run represents: real c, real proton mass.

    Only DIMENSIONLESS numbers come across from `code`, and only those that are
    mass-ratio-safe (Sec. II). Ratios used: n_0/n_e,ab, T_0/T_e,ab, beta_ab,
    beta_0, v_p/c_sim, v_sh/C_s,ab, lambda_ab, tau_sim/t_ab, L/d_i,ab.
    """
    mu = MU_PHYS
    n_ab = n_ab_cm3 * 1e6                       # -> m^-3
    n_0 = n_ab * (code["n_0"] / code["n_ab"])
    Te_0 = Te_ab_eV * (code["th_0"] / code["th_e"])

    w_pe = wpe(n_ab)
    de_ab = C / w_pe                            # REAL c here, not c_sim
    di_ab = de_ab * math.sqrt(mu)
    di_0 = di_ab * math.sqrt(n_ab / n_0)

    Cs_ab = math.sqrt(Te_ab_eV / (mu * ME_C2_EV)) * C
    vte_ab = math.sqrt(Te_ab_eV / ME_C2_EV) * C
    Cs_0 = math.sqrt(Te_0 / (mu * ME_C2_EV)) * C
    t_ab = di_ab / Cs_ab

    # The reduced speed of light, read off the velocity anchor (Sec. II:
    # c_sim = sqrt(mu_p/T_e,ab) C_s,ab). This is Table I's c_sim/c_phys row.
    c_sim = code["c_over_Csab"] * Cs_ab

    # beta_ab is dimensionless and mass-ratio-safe -> it SETS B0. (User's step 4.)
    p_ab = n_ab * Te_ab_eV * QE                 # electron pressure [Pa]
    B0 = math.sqrt(MU0 * p_ab / code["beta_ab"])
    vA = B0 / math.sqrt(MU0 * n_0 * mu * ME)
    wci0_inv = mu * ME / (QE * B0)

    vp = code["vp"] / 1.0 * c_sim                # v_p is quoted in units of c_sim
    vsh = code["vsh"] / code["Cs_ab"] * Cs_ab    # v_sh is quoted in units of C_s,ab

    # ---- collisionality: lambda_ab is the invariant that DOES transfer ----
    mfp_ab = code["lambda_ab"] * de_ab
    nu_ei_ab = vte_ab / mfp_ab
    lnL = nu_ei_ab / nu_ei_formulary(n_ab_cm3, Te_ab_eV, 1.0, coeff)

    return dict(
        mu=mu, n_ab=n_ab, n_0=n_0, Te_ab=Te_ab_eV, Te_0=Te_0,
        de_ab=de_ab, di_ab=di_ab, di_0=di_0, Cs_ab=Cs_ab, vte_ab=vte_ab, Cs_0=Cs_0,
        t_ab=t_ab, c_sim=c_sim, B0=B0, vA=vA, wci0_inv=wci0_inv, vp=vp, vsh=vsh,
        MA=vsh / vA, Mms=vsh / math.hypot(vA, Cs_0),
        beta_ab=code["beta_ab"], beta_0=code["beta_0"],
        mfp_ab=mfp_ab, nu_ei_ab=nu_ei_ab, nu_t_ab=nu_ei_ab * t_ab,
        tau_ei=1.0 / nu_ei_ab, lnL=lnL, lnL_nrl=lnL_nrl(n_ab_cm3, Te_ab_eV),
        w_pe=w_pe,
    )


# --------------------------------------------------------------------------- #
# the WarpX column
# --------------------------------------------------------------------------- #
def warpx_column(cfg, code, coeff):
    """What our SI deck actually contains: REAL c, but the code column's mu_p and theta.

    WarpX has no reduced-c option, so the only way to run PSC's dimensionless
    problem is to keep theta_e,ab = 0.092 at the real c. Every dimensionless row
    is then exact, at the price of absolute values that look unphysical
    (T_e,ab = 47 keV, B0 = 70 T) and a Coulomb logarithm that must be a dial.
    """
    ref = cfg["reference"]
    mu_p, n0 = float(ref["mass_ratio"]), float(ref["n0"])
    Te_ab = code["th_e"] * ME_C2_EV
    Te_0 = code["th_0"] * ME_C2_EV

    w_pe = wpe(n0)
    de_ab = C / w_pe
    di_ab = de_ab * math.sqrt(mu_p)
    n_0 = n0 * (code["n_0"] / code["n_ab"])
    di_0 = di_ab * math.sqrt(n0 / n_0)

    Cs_ab = code["Cs_ab"] * C
    vte_ab = code["vte_ab"] * C
    t_ab = di_ab / Cs_ab
    B0 = float(cfg["field"]["B0_tesla"])
    vA = B0 / math.sqrt(MU0 * n_0 * mu_p * ME)

    mfp_ab = code["lambda_ab"] * de_ab
    nu_ei_ab = vte_ab / mfp_ab
    lnL = nu_ei_ab / nu_ei_formulary(n0 / 1e6, Te_ab, 1.0, coeff)

    return dict(
        mu=mu_p, n_ab=n0, n_0=n_0, Te_ab=Te_ab, Te_0=Te_0,
        de_ab=de_ab, di_ab=di_ab, di_0=di_0, Cs_ab=Cs_ab, vte_ab=vte_ab,
        Cs_0=code["Cs_0"] * C, t_ab=t_ab, B0=B0, vA=vA,
        wci0_inv=mu_p * ME / (QE * B0),
        vp=code["vp"] * C, vsh=code["vsh"] * C,
        mfp_ab=mfp_ab, nu_ei_ab=nu_ei_ab, nu_t_ab=nu_ei_ab * t_ab,
        lnL=lnL, lnL_nrl=lnL_nrl(n0 / 1e6, Te_ab), w_pe=w_pe,
    )


# --------------------------------------------------------------------------- #
# rendering
# --------------------------------------------------------------------------- #
def sci(x, unit="1/s", digits=4):
    """Plain scientific notation. Used for rates, where an SI prefix on "1/s"
    would read as a nonsense unit ("2.1 T1/s")."""
    return f"{x:.{digits}g} {unit}"


def eng(x, unit, digits=3):
    """Format a number with an SI prefix chosen so the mantissa is 1-999."""
    if x == 0:
        return f"0 {unit}"
    prefixes = [(1e-15, "f"), (1e-12, "p"), (1e-9, "n"), (1e-6, "u"), (1e-3, "m"),
                (1.0, ""), (1e3, "k"), (1e6, "M"), (1e9, "G"), (1e12, "T"), (1e15, "P")]
    for scale, p in reversed(prefixes):
        if abs(x) >= scale:
            return f"{x / scale:.{digits}g} {p}{unit}"
    return f"{x:.{digits}g} {unit}"


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Schaeffer 2020 Table I in PSC code units, physical units and WarpX SI.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", nargs="?", default="runs/R1_paper",
                    help="run directory holding config.yaml (default: runs/R1_paper)")
    ap.add_argument("--n-ab-cm3", type=float, default=6.0e20,
                    help="FREE PARAMETER 1: real ablation density [cm^-3] (Table I: 6e20)")
    ap.add_argument("--Te-ab-eV", type=float, default=470.0,
                    help="FREE PARAMETER 2: real ablation temperature [eV]; this is the "
                         "reduced speed of light (Table I: 470)")
    ap.add_argument("--lambda-ab", type=float, default=None,
                    help="FREE PARAMETER 3: collisionality mfp_e,ab/d_e,ab "
                         "(default: from the config's collisions.target)")
    ap.add_argument("--nu-coeff", type=float, default=3.95e-6,
                    help="coefficient in nu_ei = C*n*lnL*T^-1.5. Default 3.95e-6 is what "
                         "reproduces Table I's tau_ei = 0.43 ps at lnLambda = 10; NRL's "
                         "value is 2.91e-6 (pass it to see the 1.36x sensitivity)")
    ap.add_argument("--show-work", action="store_true",
                    help="print the Coulomb-logarithm algebra step by step")
    args = ap.parse_args(argv)

    cfg = kcfg.load(Path(args.run) / "config.yaml")
    lam = args.lambda_ab
    if lam is None:
        tgt = (cfg.get("collisions") or {}).get("target") or {}
        if tgt.get("quantity") != "lambda_ab":
            ap.error("config has no lambda_ab target; pass --lambda-ab explicitly")
        lam = float(tgt["value"])

    co = code_column(cfg, lam)
    ph = phys_column(co, args.n_ab_cm3, args.Te_ab_eV, args.nu_coeff)
    wx = warpx_column(cfg, co, args.nu_coeff)

    w = (34, 22, 18, 20, 14)
    print(f"\nSchaeffer 2020 Table I  --  run {cfg['meta']['run_id']}")
    print(f"free parameters:  n_e,ab = {args.n_ab_cm3:.3g} cm^-3   "
          f"T_e,ab = {args.Te_ab_eV:g} eV   lambda_ab = {lam:g}\n")
    print(f"{'Parameter':<{w[0]}}{'PSC (code)':<{w[1]}}{'Physical':<{w[2]}}"
          f"{'WarpX (SI)':<{w[3]}}{'Table I':<{w[4]}}")
    print("-" * sum(w))

    def row(name, cv, pv, wv, key=None, phys_ref=True):
        t1 = TABLE1.get(key or name, ("", ""))
        ref = t1[1] if phys_ref and t1[1] else t1[0]
        print(f"{name:<{w[0]}}{cv:<{w[1]}}{pv:<{w[2]}}{wv:<{w[3]}}{ref:<{w[4]}}")

    def sec(t):
        print(f"\n{t}")

    sec("Ablation")
    row("n_e,ab", f"{co['n_ab']:g}", f"{ph['n_ab']/1e6:.3g} cm^-3",
        f"{wx['n_ab']:.3g} m^-3", "n_e,ab")
    row("T_e,ab", f"{co['th_e']:g} m_e c^2", eng(ph["Te_ab"], "eV"),
        eng(wx["Te_ab"], "eV"), "T_e,ab")
    row("C_s,ab", f"{co['Cs_ab']:.4f} c", eng(ph["Cs_ab"], "m/s"),
        eng(wx["Cs_ab"], "m/s"), "C_s,ab")
    row("v_te,ab", f"{co['vte_ab']:.4f} c", eng(ph["vte_ab"], "m/s"),
        eng(wx["vte_ab"], "m/s"))
    row("d_e,ab", f"{co['de_ab']:g}", eng(ph["de_ab"], "m"), eng(wx["de_ab"], "m"))
    row("d_i,ab", f"{co['di_ab']:g} d_e,ab", eng(ph["di_ab"], "m"),
        eng(wx["di_ab"], "m"))
    row("t_ab = d_i,ab/C_s,ab", f"{co['t_ab']:.4g} /w_pe", eng(ph["t_ab"], "s"),
        eng(wx["t_ab"], "s"))
    row("1/w_pe", "1", eng(1 / ph["w_pe"], "s"), eng(1 / wx["w_pe"], "s"))
    row("v_p", f"{co['vp']:g} c", eng(ph["vp"], "m/s"), eng(wx["vp"], "m/s"), "v_p")
    row("v_sh", f"{co['vsh']/co['Cs_ab']:g} C_s,ab", eng(ph["vsh"], "m/s"),
        eng(wx["vsh"], "m/s"), "v_sh")

    sec("Upstream")
    row("B_0", f"{co['B']:g} sqrt(m_e c^2)", eng(ph["B0"], "T"),
        eng(wx["B0"], "T"), "B_0")
    row("n_e0", f"{co['n_0']:g}", f"{ph['n_0']/1e6:.3g} cm^-3",
        f"{wx['n_0']:.3g} m^-3", "n_e0")
    row("T_0", f"{co['th_0']:g} m_e c^2", eng(ph["Te_0"], "eV"),
        eng(wx["Te_0"], "eV"), "T_0")
    row("d_i0", f"{co['di_0']/co['di_ab']:.4g} d_i,ab", eng(ph["di_0"], "m"),
        eng(wx["di_0"], "m"), "d_i0")
    row("v_A", f"{co['vA']:.4g} c", eng(ph["vA"], "m/s"), eng(wx["vA"], "m/s"))
    row("1/w_ci0", f"{co['wci0_inv']/co['t_ab']:.4g} t_ab",
        eng(ph["wci0_inv"], "s"), eng(wx["wci0_inv"], "s"), "1/w_ci0")

    sec("Dimensionless (transfer between columns)")
    row("m_i/m_e", f"{co['mu']:g}", f"{ph['mu']:.4g}", f"{wx['mu']:g}",
        "m_i/m_e", phys_ref=False)
    row("c_sim/c_phys", f"{co['c_over_Csab']:.4g} C_s,ab",
        f"{ph['c_sim']/C:.4f}", "1 (real c)", "c_sim/c_phys", phys_ref=False)
    row("beta_ab = mu0 n T/B^2", f"{co['beta_ab']:.4g}", f"{ph['beta_ab']:.4g}",
        f"{co['beta_ab']:.4g}", "beta_ab", phys_ref=False)
    row("beta_0", f"{co['beta_0']:.4g}", f"{ph['beta_0']:.4g}",
        f"{co['beta_0']:.4g}", "beta_0", phys_ref=False)
    row("M_A", f"{co['MA']:.4g}", f"{ph['MA']:.4g}", f"{wx['vsh']/wx['vA']:.4g}",
        "M_A", phys_ref=False)
    row("M_ms", f"{co['Mms']:.4g}", f"{ph['Mms']:.4g}", "", "M_ms", phys_ref=False)

    sec("Collisions  (lambda_ab transfers; nu_ei t_ab does NOT -- Sec. II)")
    row("lambda_ab = mfp/d_e,ab", f"{co['lambda_ab']:g}", f"{co['lambda_ab']:g}",
        f"{co['lambda_ab']:g}", "lambda_ab", phys_ref=False)
    row("mfp_e,ab", f"{co['mfp_ab']:g} d_e,ab", eng(ph["mfp_ab"], "m"),
        eng(wx["mfp_ab"], "m"))
    row("nu_ei,ab", f"{co['nu_ei_ab']:.4g} w_pe", sci(ph["nu_ei_ab"]),
        sci(wx["nu_ei_ab"]))
    row("tau_ei,ab", f"{1/co['nu_ei_ab']/co['t_ab']:.4g} t_ab",
        eng(1 / ph["nu_ei_ab"], "s"), eng(1 / wx["nu_ei_ab"], "s"), "tau_ei,ab")
    row("nu_ei,ab * t_ab  (mu/lambda)", f"{co['nu_t_ab']:.4g}",
        f"{ph['nu_t_ab']:.4g}", f"{wx['nu_t_ab']:.4g}", "tau_ei,ab", phys_ref=False)
    print()
    row(f"ln Lambda (coeff {args.nu_coeff:.3g})", "n/a (Takizuka-Abe)",
        f"{ph['lnL']:.3g}", f"{wx['lnL']:.4g}  (DIAL)")
    # units.py -- and therefore every deck -- uses the NRL coefficient, so print the
    # dial under it too or the deck value looks like it disagrees with this table.
    scale = args.nu_coeff / NU_EI_NRL
    row("ln Lambda (NRL coeff 2.91e-6)", "n/a", f"{ph['lnL']*scale:.3g}",
        f"{wx['lnL']*scale:.4g}  <- the deck")
    row("ln Lambda (NRL 24-ln(sqrt n/T))", "", f"{ph['lnL_nrl']:.3g}",
        f"{wx['lnL_nrl']:.3g}")

    # Table I's last collision row. It is ~1e6 x larger than the THERMAL upstream
    # ion-ion mfp, which is what the repo had been comparing it against; the row is
    # the DIRECTED-ion mfp (~v^4) of the piston/shocked ions streaming into the
    # ambient, i.e. the statement "the experiment is globally collisionless".
    sec("Table I's lambda_mfp/d_i0 = 350  (directed ions, ~v^4)")
    mi_p = ph["mu"] * ME
    v_ti0 = math.sqrt(ph["Te_0"] / (ph["mu"] * ME_C2_EV)) * C   # upstream ion thermal
    for label, v in (("at v_p", ph["vp"]), ("at v_sh", ph["vsh"]),
                     ("at v_ti(T_0)  [thermal]", v_ti0)):
        lam_m = mfp_directed(v, ph["n_0"], mi_p, ph["lnL"])
        over_di0 = f'{lam_m / ph["di_0"]:.4g} d_i0'
        print(f"{'  mfp ' + label:<{w[0]}}{'':<{w[1]}}{eng(lam_m, 'm'):<{w[2]}}"
              f"{over_di0:<{w[3]}}{'350' if 'v_sh' in label else '':<{w[4]}}")

    if args.show_work:
        show_work(co, ph, args)

    print()
    return 0


def show_work(co, ph, args):
    lam, mu, muc = co["lambda_ab"], ph["mu"], co["mu"]
    rate1 = nu_ei_formulary(args.n_ab_cm3, args.Te_ab_eV, 1.0, args.nu_coeff)
    print(f"""
--- the Coulomb logarithm, step by step ---------------------------------------

Sec. II defines the collisionality as a purely ELECTRON-scale ratio,

    lambda_ab = w_ce,ab/nu_ei,ab = mfp_ab/d_e,ab = {lam:g},
    mfp_ab = sqrt(T_ab/m_e)/nu_ei,ab,     d_e,ab = c/w_pe,

so it is independent of the ion mass and transfers between columns unchanged.
That is the ONLY collision quantity that does. Note

    nu_ei,ab t_ab = (v_te/mfp)(d_i/C_s) = sqrt(mu)*sqrt(mu)/lambda_ab = mu/lambda_ab,

which carries the ion mass explicitly: {muc:g}/{lam:g} = {muc/lam:.3g} in code units but
{mu:.4g}/{lam:g} = {mu/lam:.3g} physically. Sec. II flags exactly this quantity as
"only quantitatively matched at physical mass ratios", so the physical lnLambda
must be derived from lambda_ab, not from nu_ei t_ab.

  1. real electron thermal speed at T_e,ab = {args.Te_ab_eV:g} eV
        v_te,ab = sqrt(T_e,ab/m_e)              = {eng(ph['vte_ab'], 'm/s')}
  2. real electron skin depth at n_e,ab = {args.n_ab_cm3:.3g} cm^-3 (REAL c)
        d_e,ab  = c/w_pe                        = {eng(ph['de_ab'], 'm')}
  3. mean free path from the collisionality
        mfp_ab  = lambda_ab * d_e,ab            = {eng(ph['mfp_ab'], 'm')}
  4. collision rate
        nu_ei,ab = v_te,ab/mfp_ab               = {sci(ph['nu_ei_ab'])}
        tau_ei   = 1/nu_ei,ab                   = {eng(ph['tau_ei'], 's')}
                                                = {ph['tau_ei']/ph['t_ab']:.4g} t_ab
     Table I quotes tau_ei,ab = 0.009 t_ab = 0.43 ps.
  5. invert the formulary at the REAL (n, T)
        nu_ei = {args.nu_coeff:.3g} * n[cm^-3] * lnL * T[eV]^-1.5 = {rate1:.4g} * lnL
        lnLambda = {ph['nu_ei_ab']:.4g} / {rate1:.4g}      = {ph['lnL']:.3g}

  cross-check: the NRL expression 24 - ln(sqrt(n)/T) gives {ph['lnL_nrl']:.3g} at the same
  (n, T), so lambda_ab = {lam:g} corresponds to a genuinely PHYSICAL Coulomb logarithm --
  no dial is needed in the physical column.

  the 18.4x trap: routing through nu_ei t_ab = {muc/lam:.3g} (the CODE value) instead of
  lambda_ab gives nu_ei = {muc/lam:.3g}/t_ab = {sci(muc/lam/ph['t_ab'])} and hence
  lnLambda = {(muc/lam/ph['t_ab'])/rate1:.3g}, low by mu_phys/mu_p = {mu/muc:.4g}. A lnLambda below 1
  is also outside the validity of the Coulomb logarithm itself.
""")


if __name__ == "__main__":
    raise SystemExit(main())
