#!/usr/bin/env python3
"""Test A step 4: swap beam speed and electron temperature between the 470 eV and 15 keV
leading-edge fits (pickle from es_growth_rate.py --dump) and re-solve the ES dispersion.
FIT path is hard-coded to the 2026-11-18 job tmp; pass your own by editing R=."""
import pickle, sys, warnings, numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0,'scripts')
from es_growth_rate import solve
R=pickle.load(open('/global/homes/h/hhelal/.claude/jobs/2222e93a/tmp/fits_t032.pkl','rb'))
def gmax(comps):
    best=(0,None,None)
    seeds=sorted({u for _,u,_ in comps}|{0.5*(x[1]+y[1]) for x in comps for y in comps})
    for k in np.geomspace(1,2000,120):
        w=solve(k,comps,seeds,0.05)
        if w is not None and w.imag>best[0]: best=(w.imag,k,w.real/k)
    return best
def mod(cn, ub=None, vte=None):
    ib=max(range(3),key=lambda i:cn[i][1]); out=[]
    for i,(a,u,v) in enumerate(cn):
        if i==ib and ub is not None: u=ub
        if a>=5 and vte is not None: v=vte
        out.append((a,u,v))
    return out
for name,z in (("fs_470eV_dz0379_t1",2.62),("ft_15keV_t2",1.88)):
    cn=[r for r in R if r["name"]==name and abs(r["z"]-z)<0.01][0]["comps_norm"]
    print(f"{name} z={z}:")
    for lab,ub,vte in (("as measured",None,None),("beam 1.55",1.55,None),("beam 1.78",1.78,None),
                       ("v_te 0.36",None,0.36),("v_te 0.51",None,0.51),
                       ("beam 1.78 + v_te 0.36",1.78,0.36),("beam 1.55 + v_te 0.51",1.55,0.51)):
        g=gmax(mod(cn,ub,vte))
        print(f"   {lab:24s} gamma/w_pi0 = {g[0]:.3f}" + (f"   v_ph/v_sh {g[2]:.2f}" if g[1] else ""), flush=True)
