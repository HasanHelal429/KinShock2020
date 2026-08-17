#!/usr/bin/env python3
"""Full dimensionless audit: R1_paper (47 keV) vs R1_paper_470eV_ppc400 (470 eV).

The question is not "are these the same problem" -- they were BUILT to be the same
dimensionless problem (Table I). It is: which groups actually survive that construction,
and which ones silently do not, so we know how many independent knobs the difference in
shock structure could be hiding behind.

Every group is computed from the config primaries only. Groups are tagged
  [=] preserved to <1%     [x] differs
so the count of [x] rows IS the dimensionality of the difference.
"""
import math
import os
import sys

sys.path.insert(0, "/pscratch/sd/h/hhelal/KinShock2020/src")
import yaml  # noqa: E402
from kinshock import units as U  # noqa: E402

ROOT = "/pscratch/sd/h/hhelal/KinShock2020"
# Two run dirs, A then B. Override from the command line:
#     python scripts/dimensionless_audit.py runs/R1_phase/R1_paper runs/E_phase/es_47keV
DEFAULT = [f"{ROOT}/runs/R1_phase/R1_paper",
           f"{ROOT}/runs/R1_phase/R1_paper_470eV_ppc400"]
_args = [a for a in sys.argv[1:] if not a.startswith("-")]
_paths = _args if len(_args) == 2 else DEFAULT
RUNS = [(os.path.basename(p.rstrip("/")), os.path.abspath(p)) for p in _paths]


def groups(path):
    cfg = yaml.safe_load(open(path + "/config.yaml"))
    sc = U.derive(cfg)
    g = {}

    mu = sc.mass_ratio
    me, mi, c, qe = U.ME, sc.mi, U.C, U.QE

    # --- temperatures as theta = kT/(m_e c^2) ---
    th_ab = float(cfg["plasma"]["piston"]["theta_e_heat"])
    th_0 = float(cfg["plasma"]["ambient"]["theta_0"])
    n_ab, n_0 = sc.n0, sc.namb
    B0 = sc.B0

    # thermal speeds (non-rel definition sqrt(kT/m))
    vte_ab = math.sqrt(th_ab) * c
    vte_0 = math.sqrt(th_0) * c
    vti_0 = math.sqrt(th_0 / mu) * c

    # Debye lengths
    lD_ab = vte_ab / sc.wpe
    wpe0 = math.sqrt(n_0 * qe * qe / (U.EPS0 * me))
    lD_0 = vte_0 / wpe0

    # gyro / plasma frequencies
    wce = qe * B0 / me
    wci = qe * B0 / mi
    wpi0 = math.sqrt(n_0 * qe * qe / (U.EPS0 * mi))

    # ---------------- the groups that DEFINE the problem (Table I) -------------
    g["mass ratio  m_i/m_e"] = mu
    g["density     n_e0/n_e,ab"] = n_0 / n_ab
    g["density     n_t/n_e,ab"] = sc.nt / n_ab
    g["temperature T_0/T_e,ab"] = th_0 / th_ab
    g["beta_ab   (2 mu0 n kT_e,ab/B^2)"] = sc.beta_ab
    g["beta_0    (upstream e-beta)"] = sc.beta_0
    g["M_A       v_sh/v_A"] = sc.MA
    g["M_ms      v_sh/sqrt(vA^2+Cs0^2)"] = sc.Mms
    g["M_s       v_sh/C_s0"] = sc.vsh_model / sc.Cs0
    g["v_sh/C_s,ab"] = sc.vsh_model / sc.Cs_ab
    g["v_p /C_s,ab"] = sc.vp_model / sc.Cs_ab
    g["v_p /v_A   (piston Alfven Mach)"] = sc.vp_model / sc.vA

    # ---------------- scale separations ---------------------------------------
    g["d_i,ab/d_e,ab  = sqrt(mu)"] = sc.di / sc.de
    g["d_i0 /d_i,ab"] = sc.di0 / sc.di
    g["d_i0 /d_e0"] = sc.di0 / sc.de0
    g["rho_i0/d_i0   (v_p/v_A)"] = sc.rho_i0 / sc.di0
    g["rho_sh/d_i0   (v_sh/v_A)"] = (sc.vsh_model / wci) / sc.di0
    g["L_domain/d_i0"] = sc.domain_halfwidth / sc.di0
    g["L_target/d_i0"] = 2.0 * sc.di / sc.di0

    # ---------------- time separations ----------------------------------------
    g["1/w_ci0 / t_ab  = sqrt(beta_ab)"] = sc.wci0_inv / sc.t_ab
    g["w_pe,ab / w_ci0"] = sc.wpe / sc.wci0
    g["t_run * w_ci0"] = sc.dt * float(cfg["numerics"]["max_step"]) * sc.wci0

    # ---------------- collisionality ------------------------------------------
    g["lambda_ab = mfp_e/d_e,ab"] = sc.lambda_ab
    g["lnLambda (used in deck)"] = sc.coulomb_log
    g["lnLambda (NRL physical)"] = sc.coulomb_log_nrl
    g["mfp_ii,amb / d_i0"] = sc.mfp_ii_amb / sc.di0
    g["nu_ei,ab / w_ce"] = (sc.nu_ei_ab / wce) if sc.nu_ei_ab else float("nan")
    g["nu_ei,ab / w_pe"] = (sc.nu_ei_ab / sc.wpe) if sc.nu_ei_ab else float("nan")

    # ---------------- THE c-GROUPS (everything below scales with 1/c) ----------
    g["*theta_e,ab = kT_e,ab/m_e c^2"] = th_ab
    g["*v_te,ab / c"] = vte_ab / c
    g["*v_te,0  / c"] = vte_0 / c
    g["*C_s,ab  / c"] = sc.Cs_ab / c
    g["*v_A     / c"] = sc.vA / c
    g["*v_sh    / c"] = sc.vsh_model / c
    g["*v_p     / c"] = sc.vp_model / c
    g["*lambda_D,ab / d_e,ab"] = lD_ab / sc.de
    g["*lambda_D,0  / d_e0"] = lD_0 / sc.de0
    g["*lambda_D,0  / d_i0"] = lD_0 / sc.di0
    g["*w_pe,ab / w_ce  (ablation)"] = sc.wpe / wce
    g["*w_pe,0  / w_ce  (upstream)"] = wpe0 / wce
    g["*w_pi,0  / w_ci  (upstream)"] = wpi0 / wci
    g["*sigma_i = B^2/(mu0 n_i0 m_i c^2)"] = (sc.vA / c) ** 2
    g["*sigma_e = B^2/(mu0 n_e0 m_e c^2)"] = (sc.vA / c) ** 2 * mu
    g["*rho_e,ab / d_e,ab (v_te/c)"] = (vte_ab / wce) / sc.de
    g["*gamma of a T_e,ab electron"] = 1.0 + 1.5 * th_ab
    g["*gamma at v_sh (ion frame)"] = 1.0 / math.sqrt(1 - (sc.vsh_model / c) ** 2)
    g["*B0 [T]  (dimensional, for ref)"] = B0
    g["*T_e,ab [eV]"] = th_ab * U.ME_C2_EV
    g["*T_0 [eV]"] = th_0 * U.ME_C2_EV

    # ---------------- numerics (not physics, listed to separate the confound) --
    ppc = float(cfg["numerics"]["ppc"]["ambient"])
    dz = sc.dz
    g["#dz / d_e,ab"] = dz / sc.de
    g["#dz / lambda_D,ab"] = dz / lD_ab
    g["#dz / lambda_D,0"] = dz / lD_0
    g["#ppc (ambient)"] = ppc
    g["#N_D,PIC = ppc*lambda_D0/dz"] = ppc * lD_0 / dz
    g["#N_D,phys = n_0 lambda_D0^3"] = n_0 * lD_0 ** 3
    g["#steps / (1/w_ci0)"] = sc.steps_per_wci0
    return g, sc


tables = []
for label, path in RUNS:
    g, sc = groups(path)
    tables.append((label, g, sc))

keys = list(tables[0][1].keys())
(la, ga, _), (lb, gb, _) = tables

print()
print("=" * 96)
print(f"{'dimensionless group':38}{la:>20}{lb:>20}{'B/A':>10}   tag")
print("=" * 96)
section = None
for k in keys:
    a, b = ga[k], gb[k]
    r = (b / a) if (a not in (0, None) and a == a) else float("nan")
    same = abs(r - 1.0) < 0.01
    tag = "[=]" if same else "[x]"
    kind = k[0]
    lab = k.lstrip("*#")
    newsec = {"*": "C-DEPENDENT", "#": "NUMERICS"}.get(kind, "INVARIANT")
    if newsec != section:
        section = newsec
        print("-" * 96)
        print(f"  == {section} ==")
        print("-" * 96)
    print(f"{lab:38}{a:>20.5g}{b:>20.5g}{r:>10.4g}   {tag}")
print("=" * 96)

TH = "*theta_e,ab = kT_e,ab/m_e c^2"
nd = sum(1 for k in keys if not (abs(gb[k] / ga[k] - 1) < 0.01) and not k.startswith("#"))
print(f"\nphysics groups that DIFFER: {nd}")
eps = math.sqrt(gb[TH] / ga[TH])
# The power column only means anything when the two runs actually differ in eps. When
# they do not (e.g. R1_paper vs the top ladder rung, both 47 keV) log(eps) -> 0 and every
# exponent blows up to a meaningless six-figure number, so suppress it rather than print it.
show_pow = abs(eps - 1.0) > 0.01
print("all differing physics groups" + (", as a power of the single ratio "
      f"v_te,ab/c  (B/A = {eps:.4g}):" if show_pow
      else " (eps is the SAME in both runs -- power column suppressed):"))
for k in keys:
    if k.startswith("#"):
        continue
    r = gb[k] / ga[k]
    if abs(r - 1) < 0.01:
        continue
    if not show_pow:
        print(f"   {k.lstrip('*'):38} B/A = {r:>12.5g}")
        continue
    p = math.log(abs(r)) / math.log(eps) if r > 0 else float("nan")
    print(f"   {k.lstrip('*'):38} B/A = {r:>12.5g}   = (v_te/c ratio)^{p:+.3f}")
print()
