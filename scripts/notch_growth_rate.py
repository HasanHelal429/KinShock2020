#!/usr/bin/env python3
"""Notch GROWTH RATE across M_phase -- the discriminator c/v_ti cannot explain.

THE DEFECT IS A RATE, NOT A DEPTH. The 470 eV runs do not fail to form a shock; they form
a notch/precursor too early and grow it too fast. RESULTS quantifies it: dz8+s3f8 (kinetic,
470 eV) reaches wedge 2.25 d_i0 at t* = 0.26, hybrid and R1_paper at t* = 0.50. As a rate
over t* = 0.13 -> 0.26 that is 9.62 d_i0 per t* against R1_paper's 5.77 -- a factor 1.67.

TWO CANDIDATES E_phase CANNOT SEPARATE, both going as 1/eps_0:
    c/v_ti     = sqrt(mu)/eps_0      2261 at 470 eV vs 226 at 47 keV
    w_pe0/w_ce = sqrt(beta_0)/eps_0   100 at 470 eV vs  10 at 47 keV
Their ratio sqrt(mu/beta_0) = 22.6 is FIXED in E_phase, which is why RESULTS 2026-08-11
could pin both defects on c/v_ti while E_phase pinned the wedge on eps, with no conflict.

WHY M_phase SETTLES IT. c/v_ti depends on T_0 and mu ONLY, both frozen across M_phase, so
it sits at 226.0 at every rung while w_pe0/w_ce moves 10 -> 100. The grid is frozen too --
dz, dt, n_cell, ppc, N_D and dz/lambda_D are identical. A rate that still moves cannot be
c/v_ti and cannot be numerics.

AND THE PAIR OF SCANS EXONERATES M_A. It is pinned at 14 in E_phase but runs 14 -> 140 in
M_phase over the same decade in w_pe0/w_ce. If M_A drove the rate the scans would disagree.

UNITS: rates are rho_i0 per t*. Comparing the RATIO to E_phase's d_i0 rates is legitimate
because rho_i0/d_i0 = 10.40 is fixed in E_phase; it is not fixed across M_phase, so only
the ratio is compared.
"""
import sys, numpy as np
sys.path.insert(0,'src')
import kinshock
from kinshock import io as kio
from kinshock.units import C, ME
import yt; yt.set_log_level(50)
SPREAD,STEP,WIN=0.06,0.002,0.024
def wedge(d,T):
    cfg=kinshock.load(d); sc=kinshock.units.derive(cfg); mi=cfg["reference"]["mass_ratio"]*ME
    pfs=kio.plotfiles(d); times=np.array([kio.load_frame(p).time*sc.wci0 for p in pfs])
    i=int(np.argmin(np.abs(times-T))); ad=kio.load_frame(pfs[i]).ds.all_data()
    z=np.asarray(ad[("amb_ions","particle_position_x")])/sc.rho_i0
    p=np.asarray(ad[("amb_ions","particle_momentum_z")]); g=np.sqrt(1+(p/(mi*C))**2)
    v=p/(mi*g)/sc.vsh_model
    zp=np.percentile(np.asarray(ad[("piston_ions","particle_position_x")]),99.5)/sc.rho_i0
    box=sc.domain_halfwidth/sc.rho_i0
    offs=np.arange(0,max(box-zp-WIN,4*WIN),STEP); sig=np.full(offs.size,np.nan)
    for k,lo in enumerate(offs):
        m=(z-zp>=lo)&(z-zp<lo+WIN)
        if m.sum()>100: sig[k]=v[m].std()
    ok=np.isfinite(sig)
    for k in range(1,offs.size):
        if ok[k] and ok[k-1] and sig[k]<SPREAD<=sig[k-1]:
            f=(sig[k-1]-SPREAD)/(sig[k-1]-sig[k]); return offs[k-1]+f*STEP, times[i]
    return float('nan'), times[i]
RUNS=[("es_47keV","runs/E_phase/es_47keV",10.0),("mb_r20","runs/M_phase/mb_r20",20.0),
      ("mb_r32","runs/M_phase/mb_r32",31.6),("mb_r56","runs/M_phase/mb_r56",56.2),
      ("mb_r100","runs/M_phase/mb_r100",100.0)]
print("NOTCH GROWTH RATE across M_phase -- c/v_ti FROZEN at 226.0 at every rung")
print(f"{'run':10}{'wpe0/wce':>10}{'w(.13)':>9}{'w(.26)':>9}{'dw/dt*':>10}{'vs rung1':>10}")
base=None
for n,d,w in RUNS:
    w1,t1=wedge(d,0.131); w2,t2=wedge(d,0.262)
    r=(w2-w1)/(t2-t1)
    if base is None: base=r
    print(f"{n:10}{w:>10.1f}{w1:>9.4f}{w2:>9.4f}{r:>10.3f}{r/base:>10.2f}", flush=True)
