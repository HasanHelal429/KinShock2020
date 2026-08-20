#!/usr/bin/env python3
"""Magnetization, Mach numbers, upstream beta and magnetic Reynolds number per eps rung.

Every quantity is computed from the config primaries, with the definition written out, so
none of it depends on which script quoted what. Two conventions matter here and are the
project's, not this script's:

  * beta = mu0 n kT / B^2, NO factor of 2 (pinned 2026-08-03; Table I's own code primaries
    reproduce 1150 and 0.2 exactly under this convention and would give 2300 / 0.4 under
    the other). Electron and ion betas are reported separately since T_i = T_e here.

  * RESISTIVITY IS SET BY ELECTRON-ION COLLISIONS, so R_m uses nu_ei -- not the ion-ion
    rate that criterion 2 uses for "is it collisionless". CLAUDE.md warns those two are
    different quantities that disagree by ~10^3; conflating them is the documented trap.

R_m has a neat exact form here. Spitzer resistivity eta = m_e nu_ei/(n e^2), so the
magnetic diffusivity is

    eta_m = eta/mu0 = nu_ei * m_e/(mu0 n e^2) = nu_ei * d_e^2

because d_e^2 = m_e/(mu0 n e^2) exactly. So R_m = L V / (nu_ei d_e^2), and at the ablation
scale (L = d_e,ab, V = C_s,ab) that collapses to

    R_m,ab = C_s,ab/(nu_ei d_e,ab) = lambda_ab / sqrt(mu)

i.e. Table I's lambda_ab = 20 at mu = 100 IS a magnetic Reynolds number of 2, by
construction. That is the sense in which Sec. II's "the lambda_ab scaling ensures that
dimensionless quantities such as the magnetic Reynolds number are correct" is true.
"""
import math
import sys

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
import kinshock  # noqa: E402
from kinshock import units as U  # noqa: E402

R = "/pscratch/sd/h/hhelal/KinShock2020/runs"
ROWS = [
    ("470 eV anchor", f"{R}/S_phase/ss_dz16_ppc100"),
    ("1.5 keV",       f"{R}/E_phase/es_1p5keV"),
    ("4.7 keV",       f"{R}/E_phase/es_4p7keV"),
    ("15 keV",        f"{R}/E_phase/es_15keV"),
    ("47 keV",        f"{R}/E_phase/es_47keV"),
    ("R1_paper",      f"{R}/R1_phase/R1_paper"),
]


def numbers(path):
    cfg = kinshock.load(path)
    sc = U.derive(cfg)
    mu = sc.mass_ratio
    B0, c, qe, me, mi = sc.B0, U.C, U.QE, U.ME, sc.mi
    th_ab = float(cfg["plasma"]["piston"]["theta_e_heat"])
    th_0 = float(cfg["plasma"]["ambient"]["theta_0"])
    n0 = sc.namb
    T0_eV = th_0 * U.ME_C2_EV
    lnL = sc.coulomb_log

    d = {}
    d["eps"] = math.sqrt(th_ab)

    # ---- MAGNETIZATION ------------------------------------------------------------
    # sigma = B^2/(mu0 n m c^2) = (v_A/c)^2 with the matching mass; the relativistic
    # "cold magnetization" of the upstream.
    vA_c = sc.vA / c
    d["sigma_i"] = vA_c ** 2                       # ion (v_A/c)^2
    d["sigma_e"] = vA_c ** 2 * mu                  # electron
    # the NON-relativistic magnetization: electron gyro vs plasma response
    wpe0 = math.sqrt(n0 * qe * qe / (U.EPS0 * me))
    wce = qe * B0 / me
    d["wce/wpe0"] = wce / wpe0                     # = 1/(w_pe/w_ce)
    d["wpe0/wce"] = wpe0 / wce
    d["rho_e/lD"] = wpe0 / wce                     # identical group, stated plainly

    # ---- MACH NUMBERS --------------------------------------------------------------
    vA, Cs0 = sc.vA, sc.Cs0
    # fast magnetosonic at 90 deg (perpendicular shock): v_f^2 = v_A^2 + C_s^2
    d["M_A"] = sc.vsh_model / vA
    d["M_s"] = sc.vsh_model / Cs0
    d["M_ms"] = sc.vsh_model / math.sqrt(vA * vA + Cs0 * Cs0)
    d["M_A,piston"] = sc.vp_model / vA

    # ---- UPSTREAM BETA -------------------------------------------------------------
    # project convention: beta = mu0 n kT/B^2 (no 2). T_i = T_e = T_0 upstream.
    beta_e = U.MU0 * n0 * (T0_eV * qe) / (B0 * B0)
    d["beta_0,e"] = beta_e
    d["beta_0,i"] = beta_e                          # T_i = T_e here
    d["beta_0,tot"] = 2.0 * beta_e
    d["beta_ab"] = sc.beta_ab

    # ---- MAGNETIC REYNOLDS ---------------------------------------------------------
    # eta_m = nu_ei * d_e^2 exactly (see module docstring).
    # (a) ablation scale, the one lambda_ab controls
    d["R_m,ab"] = sc.Cs_ab / (sc.nu_ei_ab * sc.de)
    d["lambda_ab/sqrt(mu)"] = sc.lambda_ab / math.sqrt(mu)   # must equal R_m,ab
    # (b) shock scale, upstream conditions, SAME deck lnLambda
    n0_cgs = n0 * 1e-6
    nu_ei_0 = U.NU_EI_NRL * n0_cgs * lnL * T0_eV ** -1.5
    d["nu_ei,0 [1/s]"] = nu_ei_0
    d["R_m,shock"] = sc.vsh_model * sc.di0 / (nu_ei_0 * sc.de0 ** 2)
    # (c) global: domain scale
    d["R_m,domain"] = sc.vsh_model * sc.domain_halfwidth / (nu_ei_0 * sc.de0 ** 2)
    d["lnLambda"] = lnL
    d["T_0 [eV]"] = T0_eV
    d["B0 [T]"] = B0
    return d


data = [(lab, numbers(p)) for lab, p in ROWS]
keys = [
    ("--- MAGNETIZATION ---", None),
    ("sigma_i = B^2/(mu0 n_i m_i c^2)", "sigma_i"),
    ("sigma_e = B^2/(mu0 n_e m_e c^2)", "sigma_e"),
    ("w_ce/w_pe0", "wce/wpe0"),
    ("w_pe0/w_ce = rho_e/lambda_D", "wpe0/wce"),
    ("--- MACH NUMBERS ---", None),
    ("M_A    = v_sh/v_A", "M_A"),
    ("M_ms   = v_sh/sqrt(vA^2+Cs0^2)", "M_ms"),
    ("M_s    = v_sh/C_s0", "M_s"),
    ("M_A,piston = v_p/v_A", "M_A,piston"),
    ("--- UPSTREAM BETA (mu0 n kT/B^2) ---", None),
    ("beta_0,e", "beta_0,e"),
    ("beta_0,i", "beta_0,i"),
    ("beta_0,total", "beta_0,tot"),
    ("beta_ab (reference)", "beta_ab"),
    ("--- MAGNETIC REYNOLDS (eta_m = nu_ei d_e^2) ---", None),
    ("R_m,ab    = C_s,ab d_e,ab/eta_m", "R_m,ab"),
    ("  check: lambda_ab/sqrt(mu)", "lambda_ab/sqrt(mu)"),
    ("R_m,shock = v_sh d_i0/eta_m,0", "R_m,shock"),
    ("R_m,domain= v_sh L/eta_m,0", "R_m,domain"),
    ("--- inputs, for reference ---", None),
    ("eps = v_te,ab/c", "eps"),
    ("lnLambda (deck)", "lnLambda"),
    ("T_0 [eV]", "T_0 [eV]"),
    ("B0 [T]", "B0 [T]"),
    ("nu_ei,0 [1/s]", "nu_ei,0 [1/s]"),
]

W = 13
print()
print(f"{'quantity':36}" + "".join(f"{lab:>{W}}" for lab, _ in data))
print("=" * (36 + W * len(data)))
for label, k in keys:
    if k is None:
        print(f"\n{label}")
        continue
    vals = [d[k] for _, d in data]
    row = "".join(f"{v:>{W}.4g}" for v in vals)
    lad = vals[:5]
    spread = max(lad) / min(lad) if min(lad) > 0 else float("inf")
    tag = "  [=] held" if spread < 1.01 else f"  [x] {spread:.3g}x"
    print(f"{label:36}{row}{tag}")
print()
print("[=]/[x] compares the FIVE LADDER RUNGS only (R1_paper is an independent 47 keV run")
print("on a different domain, so it is shown for reference, not included in the spread).")
