#!/usr/bin/env python
"""Schaeffer 2020 Table I from a self-consistent set of PHYSICAL scales.

THE INPUT IS A CHOICE OF REAL PLASMA
------------------------------------
A PIC run of this problem is a set of dimensionless numbers and corresponds to a
whole family of real plasmas. You pick one member by choosing, in physical units:

    n_e,ab      ablation electron density        6e26 m^-3
    T_e,ab      ablation electron temperature    470 eV
    lambda_ab   collisionality, mfp_e,ab/d_e,ab  20
    mu          m_i/m_e                          100
    beta_ab     ablation beta, mu0 n kT/B^2      1150
    T_0         upstream temperature             10 eV
    n_e0/n_e,ab upstream density fraction        0.008

Everything else follows, and B0 follows from beta_ab. Note mu = 100 is a *physical*
choice here, not just a code convenience: the represented plasma really is a
light-ion plasma. That is what makes the set self-consistent, and it is why the ion
rows below do NOT reproduce Table I's own SI column, which is a real hydrogen
(mu = 1836) plasma -- the paper's caption calls that "one possible set of
experimentally-relevant physical values", i.e. an illustration rather than a unit map.

beta_0 IS NOT FREE. Given n_e0, T_0 and B0 (itself set by beta_ab) it is determined:
mu0 n_e0 kT_0/B0^2 = 0.196 here. Table I quotes 0.2, which agrees.

THE ONE DEGREE OF FREEDOM LEFT: THE SPEED OF LIGHT
--------------------------------------------------
Sec. II sets c = sqrt(mu_p/T_e,ab) C_s,ab, so the code's theta_e,ab IS its speed of
light. PSC picks theta_e,ab = 0.092 (O(0.1), to stay non-relativistic), which for
T_e,ab = 470 eV means c_sim/c_phys = sqrt(theta_e_phys/0.092) = 0.100. WarpX has no
reduced-c option, so it must use the *physical* theta_e = T_e,ab/(m_e c^2) = 9.2e-4.

Both codes then represent the SAME physical plasma and reproduce every dimensionless
row -- beta_ab, beta_0, M_A, M_ms, lambda_ab, n and T ratios -- except the one that
IS the speed of light, c/C_s,ab: 33.0 for PSC against 329.7 for WarpX. Two
consequences, both reported below and neither optional:

  * WarpX needs 10x the timesteps for the same 220 t_ab, because dt is CFL-locked to
    dz/c while t_ab ~ c.
  * dz/lambda_D,ab goes from ~1 (resolved) to ~10 (under-resolved), so the paper's
    dz = 0.3 d_e,ab grid will grid-heat at real c. Fixing that costs another 10x in
    both dz and dt, i.e. ~100x overall.

Run with --show-work for the Coulomb-logarithm algebra, --deck for the config diff.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from kinshock import config as kcfg                                    # noqa: E402
from kinshock.units import (C, EPS0, ME, ME_C2_EV, ME_C2_J, MP, MU0,   # noqa: E402
                            NU_EI_NRL, QE)

MU_PHYS = MP / ME          # 1836.15, for reference only -- NOT used in the columns

# Table I as printed (schaeffer2020.pdf p.4), for the check column.
TABLE1 = {
    "n_e,ab": ("1.25", "6e20 cm^-3"), "T_e,ab": ("0.092 m_e c^2", "470 eV"),
    "tau_ei,ab": ("0.009 t_ab", "0.43 ps"), "C_s,ab": ("0.030 c", "210 km/s"),
    "v_p": ("0.104 c", "730 km/s"), "v_sh": ("4.6 C_s,ab", "980 km/s"),
    "B_0": ("0.01", "7 T"), "n_e0": ("0.01 (code)", "4.8e18 cm^-3"),
    "T_0": ("0.002 m_e c^2", "10 eV"), "d_i0": ("11.2 d_i,ab", "104 um"),
    "1/w_ci0": ("33.9 t_ab", "1.5 ns"), "m_i/m_e": ("100", ""),
    "c_sim/c_phys": ("0.02", ""), "beta_ab": ("1150", ""), "beta_0": ("0.2", ""),
    "lambda_ab": ("20", ""), "M_A": ("14", ""), "M_ms": ("13", ""),
    "Lz": ("900 d_i,ab", "8.4 mm"), "tau_sim": ("220 t_ab", "10.9 ns"),
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
    """Rutherford momentum-transfer mfp [m] of an ion of speed v in density n.

    lambda = 4 pi eps0^2 mi^2 v^4 / (n Z^4 e^4 lnLambda). This is Table I's
    lambda_mfp/d_i0 = 350 row; the thermal ion-ion mfp is ~1e6 smaller.
    """
    return 4.0 * math.pi * EPS0 ** 2 * mi ** 2 * v ** 4 / (n * Z ** 4 * QE ** 4 * lnL)


# --------------------------------------------------------------------------- #
# the one derivation, run once per choice of speed of light
# --------------------------------------------------------------------------- #
def derive(P, c_used, geo):
    """All scales for the plasma ``P`` as simulated with speed of light ``c_used``.

    ``P`` holds the physical choices; ``geo`` the grid/time ratios from the config.
    Returns SI values plus the code-unit normalizations Table I tabulates.

    The only c-dependent quantities are the inertial lengths (d_e, d_i ~ c), the
    times built on them (t_ab ~ c), and theta_e = kT/(m_e c^2). Densities,
    temperatures, all speeds (C_s, v_te, v_A, v_p, v_sh), all frequencies
    (w_pe, w_ce, w_ci) and beta are c-INDEPENDENT, so they are identical in every
    column by construction.
    """
    n_ab, mu, lam = P["n_ab"], P["mu"], P["lambda_ab"]
    mi = mu * ME
    Te_ab_J = P["Te_ab_eV"] * QE

    theta_e = Te_ab_J / (ME * c_used ** 2)        # THE reduced speed of light
    theta_0 = P["T_0_eV"] * QE / (ME * c_used ** 2)

    w_pe = wpe(n_ab)                              # c-independent
    de_ab = c_used / w_pe
    di_ab = de_ab * math.sqrt(mu)
    n_0 = P["n0_frac"] * n_ab
    de_0, di_0 = c_used / wpe(n_0), c_used / wpe(n_0) * math.sqrt(mu)

    Cs_ab = math.sqrt(Te_ab_J / mi)               # c-independent
    vte_ab = math.sqrt(Te_ab_J / ME)
    Cs_0 = math.sqrt(P["T_0_eV"] * QE / mi)
    lam_D_ab = vte_ab / w_pe                      # c-independent
    t_ab = di_ab / Cs_ab                          # ~ c

    # beta_ab SETS B0 (Table I's convention, no factor of 2 -- see RESULTS 2026-08-03)
    B0 = math.sqrt(MU0 * n_ab * Te_ab_J / P["beta_ab"])
    beta_0 = MU0 * n_0 * P["T_0_eV"] * QE / B0 ** 2       # DERIVED, not free
    vA = B0 / math.sqrt(MU0 * n_0 * mi)                   # c-independent
    # PSC's field normalization, B/sqrt(mu0 n_ref m_e c^2): this is the row Table I
    # prints as "0.01 sqrt(m_e c^2)", and it is c-dependent, hence 0.01 for PSC and
    # 1e-3 for WarpX even though B0 is the same 7.03 T.
    n_ref = n_ab / P["n_ab_code"]
    B_code = B0 / math.sqrt(MU0 * n_ref * ME * c_used ** 2)
    # w_ci0 must be taken in the SAME normalization as the lengths and times of this
    # column, i.e. from B_code, NOT as q B0/m_i. Mixing them (physical w_ci0 against
    # a reduced-c t_ab) breaks Sec. II's exact identity 1/w_ci0 = sqrt(beta_ab) t_ab
    # and reports 339 t_ab instead of Table I's 33.9. For the real-c column the two
    # agree identically, which the assertion below pins.
    w_ci0 = B_code * wpe(n_ref) / mu
    if abs(c_used - C) < 1.0:
        assert abs(w_ci0 / (QE * B0 / mi) - 1.0) < 1e-9, "w_ci0 normalization mismatch"

    vsh = P["vsh_over_Csab"] * Cs_ab
    vp = P["vp_over_Csab"] * Cs_ab
    rho_i0 = vp / w_ci0

    # grid + time: dz is a fixed number of d_e,ab, dt is CFL-locked to dz/c_used
    dz = geo["dz_over_de"] * de_ab
    dt = geo["cfl"] * dz / c_used
    n_cell = int(round(geo["halfwidth_de"] / geo["dz_over_de"]))
    tau_sim = geo["tau_sim_over_tab"] * t_ab
    max_step = int(round(tau_sim / dt))

    # collisionality: lambda_ab = mfp/d_e,ab fixes nu_ei, which fixes lnLambda
    mfp_ab = lam * de_ab
    nu_ei_ab = vte_ab / mfp_ab
    return dict(
        c=c_used, mu=mu, mi=mi, theta_e=theta_e, theta_0=theta_0,
        n_ab=n_ab, n_0=n_0, Te_ab=P["Te_ab_eV"], T_0=P["T_0_eV"],
        w_pe=w_pe, de_ab=de_ab, di_ab=di_ab, de_0=de_0, di_0=di_0,
        Cs_ab=Cs_ab, vte_ab=vte_ab, Cs_0=Cs_0, lam_D_ab=lam_D_ab, t_ab=t_ab,
        B0=B0, B_code=B_code, beta_ab=P["beta_ab"], beta_0=beta_0,
        w_ci0=w_ci0, wci0_inv=1.0 / w_ci0, vA=vA, vp=vp, vsh=vsh, rho_i0=rho_i0,
        MA=vsh / vA, Mms=vsh / math.hypot(vA, Cs_0),
        dz=dz, dt=dt, dt_wpe=dt * w_pe, n_cell=n_cell, tau_sim=tau_sim,
        max_step=max_step, dz_over_lamD=dz / lam_D_ab,
        mfp_ab=mfp_ab, nu_ei_ab=nu_ei_ab, nu_t_ab=nu_ei_ab * t_ab,
        lambda_ab=lam, c_over_Csab=c_used / Cs_ab,
    )


def coulomb_log(d, coeff):
    """lnLambda that makes WarpX's dimensional operator reproduce d's nu_ei,ab."""
    return d["nu_ei_ab"] / nu_ei_formulary(d["n_ab"] / 1e6, d["Te_ab"], 1.0, coeff)


# --------------------------------------------------------------------------- #
# rendering
# --------------------------------------------------------------------------- #
def eng(x, unit, digits=4):
    if x == 0:
        return f"0 {unit}"
    for scale, p in reversed([(1e-15, "f"), (1e-12, "p"), (1e-9, "n"), (1e-6, "u"),
                              (1e-3, "m"), (1.0, ""), (1e3, "k"), (1e6, "M"),
                              (1e9, "G"), (1e12, "T"), (1e15, "P")]):
        if abs(x) >= scale:
            return f"{x / scale:.{digits}g} {p}{unit}"
    return f"{x:.{digits}g} {unit}"


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Table I from a self-consistent set of physical scales.")
    ap.add_argument("run", nargs="?", default="runs/R1_phase/R1_paper")
    ap.add_argument("--n-ab", type=float, default=6.0e26, help="[m^-3] (6e26)")
    ap.add_argument("--Te-ab-eV", type=float, default=470.0, help="[eV] (470)")
    ap.add_argument("--T0-eV", type=float, default=10.0, help="[eV] (10)")
    ap.add_argument("--mu", type=float, default=100.0, help="m_i/m_e (100)")
    ap.add_argument("--beta-ab", type=float, default=1150.0, help="(1150)")
    ap.add_argument("--n0-frac", type=float, default=0.008, help="n_e0/n_e,ab (0.008)")
    ap.add_argument("--beta-0", type=float, default=None,
                    help="OVER-DETERMINED: beta_0 follows from n_e0, T_0 and B0. Pass a "
                         "value only to check it against the derived one.")
    ap.add_argument("--lambda-ab", type=float, default=20.0, help="mfp/d_e,ab (20)")
    ap.add_argument("--theta-e-psc", type=float, default=0.092,
                    help="PSC's code temperature, i.e. its reduced c (0.092)")
    ap.add_argument("--nu-coeff", type=float, default=NU_EI_NRL,
                    help=f"nu_ei coefficient (default NRL {NU_EI_NRL:g}; the paper's "
                         "3.95e-6 reproduces Table I's 0.43 ps at lnLambda = 10)")
    ap.add_argument("--show-work", action="store_true")
    ap.add_argument("--deck", action="store_true",
                    help="print the config.yaml values this implies")
    args = ap.parse_args(argv)

    cfg = kcfg.load(Path(args.run) / "config.yaml")
    P = dict(n_ab=args.n_ab, Te_ab_eV=args.Te_ab_eV, T_0_eV=args.T0_eV, mu=args.mu,
             beta_ab=args.beta_ab, n0_frac=args.n0_frac, lambda_ab=args.lambda_ab,
             n_ab_code=1.25, vsh_over_Csab=float(cfg["model"]["vsh_over_Csab"]),
             # Table I's v_p = 0.104 c_sim, expressed as the c-free ratio v_p/C_s,ab
             vp_over_Csab=float(cfg["model"]["vp_over_c"])
             / math.sqrt(args.theta_e_psc / args.mu))
    geo = dict(dz_over_de=float(cfg["geometry"]["dz_over_de"]),
               cfl=float(cfg["numerics"]["cfl"]),
               halfwidth_de=float(cfg["geometry"]["domain_halfwidth_de"]),
               tau_sim_over_tab=220.0)

    # the physical theta_e is what WarpX must use; PSC's 0.092 defines its reduced c
    theta_e_phys = args.Te_ab_eV * QE / ME_C2_J
    c_sim = C * math.sqrt(theta_e_phys / args.theta_e_psc)

    psc = derive(P, c_sim, geo)
    wx = derive(P, C, geo)

    w = (32, 21, 19, 21, 13)
    print(f"\nTable I from self-consistent physical scales  --  {cfg['meta']['run_id']}")
    print(f"  CHOSEN: n_e,ab = {args.n_ab:.3g} m^-3   T_e,ab = {args.Te_ab_eV:g} eV   "
          f"lambda_ab = {args.lambda_ab:g} d_e   m_i/m_e = {args.mu:g}")
    print(f"          beta_ab = {args.beta_ab:g}   T_0 = {args.T0_eV:g} eV   "
          f"n_e0 = {args.n0_frac:g} n_e,ab")
    print(f"  DERIVED: B0 = {wx['B0']:.4f} T   beta_0 = {wx['beta_0']:.4f}   "
          f"c_sim/c_phys = {c_sim / C:.4f}\n")
    if args.beta_0 is not None and abs(args.beta_0 / wx["beta_0"] - 1.0) > 0.05:
        # beta_0 = mu0 n_e0 kT_0/B0^2 with B0 already fixed by beta_ab, so it cannot be
        # chosen independently. Report what the requested value would actually require.
        n_need = args.beta_0 * wx["B0"] ** 2 / (MU0 * args.T0_eV * QE)
        T_need = args.beta_0 * wx["B0"] ** 2 / (MU0 * wx["n_0"] * QE)
        print(f"  !! beta_0 = {args.beta_0:g} was requested but beta_0 is DERIVED, and "
              f"these inputs give {wx['beta_0']:.4f} (Table I: 0.2).")
        print(f"     {args.beta_0:g} would need n_e0 = {n_need / wx['n_ab']:.5f} n_e,ab "
              f"(not {args.n0_frac:g}) at T_0 = {args.T0_eV:g} eV,")
        print(f"     or T_0 = {T_need:.4g} eV (not {args.T0_eV:g}) at "
              f"n_e0 = {args.n0_frac:g} n_e,ab. Proceeding with {wx['beta_0']:.4f}.\n")
    print(f"{'Parameter':<{w[0]}}{'PSC (code)':<{w[1]}}{'Physical':<{w[2]}}"
          f"{'WarpX deck (SI)':<{w[3]}}{'Table I':<{w[4]}}")
    print("-" * sum(w))

    def row(name, cv, pv, wv, key=None, ref="phys"):
        t1 = TABLE1.get(key or name, ("", ""))
        r = t1[1] if ref == "phys" and t1[1] else (t1[0] if ref != "skip" else "")
        print(f"{name:<{w[0]}}{cv:<{w[1]}}{pv:<{w[2]}}{wv:<{w[3]}}{r:<{w[4]}}")

    print("\nAblation  (the physical column is c-independent unless marked ~c)")
    row("n_e,ab", f"{P['n_ab_code']:g}", f"{wx['n_ab']:.3g} m^-3",
        f"{wx['n_ab']:.3g} m^-3", "n_e,ab")
    row("T_e,ab", f"{psc['theta_e']:.4g} m_e c^2", f"{wx['Te_ab']:g} eV",
        f"{wx['Te_ab']:g} eV  (th={wx['theta_e']:.3g})", "T_e,ab")
    row("C_s,ab", f"{psc['Cs_ab'] / psc['c']:.4f} c", eng(wx["Cs_ab"], "m/s"),
        eng(wx["Cs_ab"], "m/s"), "C_s,ab")
    row("v_te,ab", f"{psc['vte_ab'] / psc['c']:.4f} c", eng(wx["vte_ab"], "m/s"),
        eng(wx["vte_ab"], "m/s"))
    row("d_e,ab", "1 d_e,ab", eng(wx["de_ab"], "m"), eng(wx["de_ab"], "m"))
    row("d_i,ab", f"{math.sqrt(P['mu']):g} d_e,ab", eng(wx["di_ab"], "m"),
        eng(wx["di_ab"], "m"))
    row("lambda_D,ab", f"{psc['lam_D_ab'] / psc['de_ab']:.4f} d_e,ab",
        eng(wx["lam_D_ab"], "m"), eng(wx["lam_D_ab"], "m"))
    row("t_ab", f"{psc['t_ab'] * psc['w_pe']:.4g} /w_pe (PSC clock)",
        eng(wx["t_ab"], "s"), eng(wx["t_ab"], "s"))
    row("v_p", f"{psc['vp'] / psc['c']:.4f} c", eng(wx["vp"], "m/s"),
        eng(wx["vp"], "m/s"), "v_p")
    row("v_sh", f"{P['vsh_over_Csab']:g} C_s,ab", eng(wx["vsh"], "m/s"),
        eng(wx["vsh"], "m/s"), "v_sh")

    print("\nUpstream")
    row("B_0", f"{psc['B_code']:.4g} sqrt(m_e c^2)", eng(wx["B0"], "T"),
        eng(wx["B0"], "T"), "B_0")
    row("n_e0", f"{P['n_ab_code'] * P['n0_frac']:g}", f"{wx['n_0']:.3g} m^-3",
        f"{wx['n_0']:.3g} m^-3", "n_e0")
    row("T_0", f"{psc['theta_0']:.4g} m_e c^2", f"{wx['T_0']:g} eV",
        f"{wx['T_0']:g} eV  (th={wx['theta_0']:.3g})", "T_0")
    row("d_i0", f"{wx['di_0'] / wx['di_ab']:.4g} d_i,ab",
        eng(wx["di_0"], "m"), eng(wx["di_0"], "m"), "d_i0")
    row("v_A", f"{psc['vA'] / psc['c']:.4g} c", eng(wx["vA"], "m/s"),
        eng(wx["vA"], "m/s"))
    row("1/w_ci0", f"{psc['wci0_inv'] / psc['t_ab']:.4g} t_ab",
        eng(wx["wci0_inv"], "s"), eng(wx["wci0_inv"], "s"), "1/w_ci0")
    row("rho_i0", f"{wx['rho_i0'] / wx['de_ab']:.4g} d_e,ab",
        eng(wx["rho_i0"], "m"), eng(wx["rho_i0"], "m"))

    print("\nDimensionless")
    row("m_i/m_e", f"{P['mu']:g}", f"{P['mu']:g}", f"{P['mu']:g}", "m_i/m_e", "code")
    row("c/C_s,ab   <- THE reduced c", f"{psc['c_over_Csab']:.4g}",
        f"{psc['c_over_Csab']:.4g}", f"{wx['c_over_Csab']:.4g}", ref="skip")
    row("c_sim/c_phys", f"{c_sim / C:.4f}", f"{c_sim / C:.4f}", "1 (real c)",
        "c_sim/c_phys", "code")
    row("beta_ab = mu0 n T/B^2", f"{psc['beta_ab']:.5g}", f"{wx['beta_ab']:.5g}",
        f"{wx['beta_ab']:.5g}", "beta_ab", "code")
    row("beta_0  (DERIVED)", f"{psc['beta_0']:.4f}", f"{wx['beta_0']:.4f}",
        f"{wx['beta_0']:.4f}", "beta_0", "code")
    row("M_A", f"{psc['MA']:.4g}", f"{wx['MA']:.4g}", f"{wx['MA']:.4g}", "M_A", "code")
    row("M_ms", f"{psc['Mms']:.4g}", f"{wx['Mms']:.4g}", f"{wx['Mms']:.4g}",
        "M_ms", "code")

    print("\nCollisions")
    row("lambda_ab = mfp/d_e,ab", f"{P['lambda_ab']:g}", f"{P['lambda_ab']:g}",
        f"{P['lambda_ab']:g}", "lambda_ab", "code")
    row("mfp_e,ab", f"{P['lambda_ab']:g} d_e,ab", eng(wx["mfp_ab"], "m"),
        eng(wx["mfp_ab"], "m"))
    row("nu_ei,ab", f"{psc['nu_ei_ab'] / psc['w_pe']:.4g} w_pe",
        f"{wx['nu_ei_ab']:.4g} 1/s", f"{wx['nu_ei_ab']:.4g} 1/s")
    row("tau_ei,ab", f"{1 / psc['nu_ei_ab'] / psc['t_ab']:.4g} t_ab",
        eng(1 / wx["nu_ei_ab"], "s"), eng(1 / wx["nu_ei_ab"], "s"), "tau_ei,ab")
    row("nu_ei,ab t_ab = mu/lambda", f"{psc['nu_t_ab']:.4g}", f"{psc['nu_t_ab']:.4g}",
        f"{wx['nu_t_ab']:.4g}", ref="skip")
    lnL_w, lnL_p = coulomb_log(wx, args.nu_coeff), coulomb_log(psc, args.nu_coeff)
    row(f"ln Lambda (coeff {args.nu_coeff:.3g})", "n/a (Takizuka-Abe)",
        f"{lnL_w:.4g}", f"{lnL_w:.4g}   <- PHYSICAL", ref="skip")
    row("ln Lambda  24-ln(sqrt n/T)", "", f"{lnL_nrl(wx['n_ab'] / 1e6, wx['Te_ab']):.4g}",
        f"{lnL_nrl(wx['n_ab'] / 1e6, wx['Te_ab']):.4g}", ref="skip")
    lam350 = mfp_directed(wx["vsh"], wx["n_0"], wx["mi"], lnL_w)
    row("mfp(v_sh)/d_i0  [Table I 350]", "", f"{lam350 / wx['di_0']:.4g} d_i0",
        f"{lam350 / wx['di_0']:.4g} d_i0", ref="skip")

    print("\nGrid  (identical cell COUNT; the cell SIZE differs because d_e,ab ~ c)")
    row("dz", f"{geo['dz_over_de']:g} d_e,ab", eng(wx["dz"], "m"), eng(wx["dz"], "m"))
    row("n_cell", f"{psc['n_cell']}", f"{psc['n_cell']}", f"{wx['n_cell']}")
    row("domain", f"{geo['halfwidth_de']:g} d_e,ab",
        eng(wx["dz"] * wx["n_cell"], "m"), eng(wx["dz"] * wx["n_cell"], "m"), "Lz")
    row("tau_sim", "220 t_ab", eng(wx["tau_sim"], "s"), eng(wx["tau_sim"], "s"),
        "tau_sim")
    row("dt", f"{psc['dt_wpe']:.4g} /w_pe", eng(wx["dt"], "s"), eng(wx["dt"], "s"))

    # The reduced c is the ONLY difference between the two codes, so it is worth
    # isolating what it buys and what it costs rather than burying it in the rows
    # above. PSC and WarpX represent the same physical plasma; PSC reaches it with a
    # 10x-slow light speed, which shrinks d_e,ab and hence both dz and t_ab by 10x.
    print("\nReduced-c ledger  (PSC vs WarpX -- the one place they differ)")
    led = (("theta_e,ab = kT_e/(m_e c^2)", f"{psc['theta_e']:.5g}",
            f"{wx['theta_e']:.5g}", "sets the code's c"),
           ("c/C_s,ab", f"{psc['c_over_Csab']:.4g}", f"{wx['c_over_Csab']:.4g}",
            "Sec. II's reduced c"),
           ("d_e,ab [m]", eng(psc["de_ab"], "m"), eng(wx["de_ab"], "m"), "~ c"),
           ("dz [m]", eng(psc["dz"], "m"), eng(wx["dz"], "m"), "= 0.3 d_e,ab"),
           ("dz/lambda_D,ab", f"{psc['dz_over_lamD']:.3g}",
            f"{wx['dz_over_lamD']:.3g}", "<= 1 or it grid-heats"),
           ("t_ab [s]", eng(psc["t_ab"], "s"), eng(wx["t_ab"], "s"), "~ c"),
           ("dt [s]", eng(psc["dt"], "s"), eng(wx["dt"], "s"), "CFL: 0.75 dz/c"),
           ("max_step for 220 t_ab", f"{psc['max_step']:,}", f"{wx['max_step']:,}",
            f"{wx['max_step'] / psc['max_step']:.0f}x more"))
    print(f"  {'':<30}{'PSC':<18}{'WarpX':<18}note")
    for name, a, b, note in led:
        print(f"  {name:<30}{a:<18}{b:<18}{note}")
    print(f"\n  PSC resolves lambda_D at 0.3 d_e,ab because its d_e,ab is only "
          f"{psc['de_ab'] / psc['lam_D_ab']:.1f} lambda_D;")
    print(f"  at real c d_e,ab is {wx['de_ab'] / wx['lam_D_ab']:.1f} lambda_D, so the SAME 0.3 d_e,ab "
          f"cell is {wx['dz_over_lamD']:.1f} lambda_D wide.")
    print(f"  Resolving it needs dz ~{wx['dz_over_lamD']:.0f}x smaller and dt with it: "
          f"~{wx['dz_over_lamD'] ** 2:.0f}x cost on top of the {wx['max_step'] / psc['max_step']:.0f}x.")

    if args.show_work:
        show_work(wx, psc, args, lnL_w)
    if args.deck:
        print_deck(wx, P, geo, lnL_w, args)
    print()
    return 0


def show_work(wx, psc, args, lnL_w):
    rate1 = nu_ei_formulary(wx["n_ab"] / 1e6, wx["Te_ab"], 1.0, args.nu_coeff)
    print(f"""
--- the Coulomb logarithm ------------------------------------------------------
lambda_ab = mfp_e,ab/d_e,ab is a pure ELECTRON ratio, so it is independent of mu and
of the ion physics; what it DOES depend on is the speed of light, through
d_e,ab = c/w_pe. WarpX runs at real c, so its d_e,ab is 1/{psc['c'] / C:.4f} = {C / psc['c']:.1f}x PSC's, and
the mean free path it must produce is that much longer:

  1. v_te,ab = sqrt(T_e,ab/m_e), T_e,ab = {args.Te_ab_eV:g} eV = {eng(wx['vte_ab'], 'm/s')}
  2. d_e,ab  = c/w_pe at n_e,ab = {wx['n_ab']:.3g} m^-3     = {eng(wx['de_ab'], 'm')}
  3. mfp_ab  = lambda_ab * d_e,ab                  = {eng(wx['mfp_ab'], 'm')}
  4. nu_ei,ab = v_te,ab/mfp_ab                     = {wx['nu_ei_ab']:.4g} 1/s
     tau_ei   = {eng(1 / wx['nu_ei_ab'], 's')} = {1 / wx['nu_ei_ab'] / wx['t_ab']:.4g} t_ab   (Table I: 0.009 t_ab)
  5. invert nu_ei = {args.nu_coeff:.3g} n[cm^-3] lnL T[eV]^-1.5 = {rate1:.4g} lnL
     lnLambda = {wx['nu_ei_ab']:.4g}/{rate1:.4g} = {lnL_w:.4g}

  Cross-check: 24 - ln(sqrt(n)/T) = {lnL_nrl(wx['n_ab'] / 1e6, wx['Te_ab']):.4g} at the same (n, T). So at
  T_e,ab = {args.Te_ab_eV:g} eV the deck needs NO dial -- lambda_ab = {args.lambda_ab:g} is what a real plasma
  at these conditions actually does. That is the payoff of choosing the physical
  temperature: nu_ei ~ T^-3/2, so dropping T_e,ab from 47 keV (theta_e = 0.092 at real
  c) to {args.Te_ab_eV:g} eV raises nu_ei by ~1000x and turns lnLambda = 1.2e5 into {lnL_w:.1f}.
""")


def print_deck(wx, P, geo, lnL_w, args):
    print(f"""
--- config.yaml values this implies -------------------------------------------
reference:
  n0: {P['n_ab']:.6g}                 # [m^-3] unchanged
  mass_ratio: {P['mu']:g}
plasma:
  piston:
    theta_e_heat: {wx['theta_e']:.6g}      # [kT_e/m_e c^2] = {args.Te_ab_eV:g} eV at REAL c (was 0.092 = 47 keV)
  ambient:
    density_over_n0: {P['n0_frac']:g}
    theta_0: {wx['theta_0']:.6g}          # [kT/m_e c^2] = {P['T_0_eV']:g} eV     (was 0.002 = 1.02 keV)
field:
  B0_tesla: {wx['B0']:.10f}       # [T] from beta_ab = {P['beta_ab']:g}   (was 70.273)
numerics:
  max_step: {wx['max_step']}              # [steps] = 220 t_ab      (was 322400, i.e. {wx['max_step'] / 322400:.0f}x)
collisions:
  target:
    quantity: lambda_ab          # UNCHANGED -- and that is the point: at 470 eV this
    value: 20.0                  # now RESOLVES to lnLambda = {lnL_w:.4g} on its own,
                                 # a physical value, where at 47 keV it needed 1.22e5.
                                 # (Do NOT switch to `value: physical`: that uses
                                 # 24-ln(sqrt n/T) = {lnL_nrl(wx['n_ab'] / 1e6, wx['Te_ab']):.4g}, which is a
                                 # different quantity and would not give lambda_ab = 20.)

geometry/ppc/operators unchanged: d_e,ab = c/w_pe depends only on n0 and real c, so
dz = {eng(wx['dz'], 'm')}, {wx['n_cell']} cells and the {eng(wx['dz'] * wx['n_cell'], 'm')} domain are all identical.

  ** two costs, both consequences of real c, neither avoidable in WarpX **
  1. {wx['max_step'] / 322400:.0f}x the timesteps for the same 220 t_ab (dt is CFL-locked to dz/c, t_ab ~ c).
  2. dz/lambda_D,ab = {wx['dz_over_lamD']:.3g} vs ~1 for PSC, so the paper's 0.3 d_e,ab grid
     UNDER-RESOLVES the Debye length and will grid-heat. Resolving it needs dz ~10x
     smaller and dt with it, i.e. ~100x cost on top.
""")


if __name__ == "__main__":
    raise SystemExit(main())
