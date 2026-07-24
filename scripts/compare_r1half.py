import os, sys, numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT,"src"))
import kinshock
from kinshock import io, metrics
from kinshock.units import C

run=os.path.join(ROOT,"runs","R1_core_half")
cfg=kinshock.load(run); sc=kinshock.units.derive(cfg)

# --- scales / paper references ---
Cs = getattr(sc,"C_s_ab",None) or getattr(sc,"Cs_ab",None) or getattr(sc,"cs_ab",None)
print("=== scales ===")
print(f"C_s,ab        = {getattr(sc,'C_s_ab','?')}")
for a in ("C_s_ab","cs_ab","Cs_ab","vp_model","vsh_model","rho_i0","di0","de","wci0","namb","B0","mi"):
    if hasattr(sc,a): print(f"  {a:12s}= {getattr(sc,a)}")

pfs=io.plotfiles(run)
frames=[io.load_frame(p) for p in pfs]
print(f"\n{len(frames)} frames, t*wci0 in [{frames[0].time*sc.wci0:.3f},{frames[-1].time*sc.wci0:.3f}]")

# --- shock front trajectory & speed (same as fig_trajectory) ---
zexcl = 2.0*sc.di0  # exclude piston slab; matches figure default region
ts,zs=[],[]
for fr in frames:
    n=io.species_density(fr, cfg["ion_species"])
    zf=metrics.track_front(fr.z_centers, n, sc.namb, threshold=1.5, z_exclude=zexcl)
    ts.append(fr.time); zs.append(zf)
ts=np.array(ts); zs=np.array(zs)
vsh=metrics.speed_from_trajectory(ts, zs)   # m/s
print("\n=== shock front / speed (second-half linear fit of density front) ===")
print(f"  v_sh = {vsh:.4e} m/s = {vsh/C:.4f} c = {vsh/1e3:.0f} km/s")
Cs_c = sc.C_s_ab/C if hasattr(sc,"C_s_ab") else 0.030
print(f"  v_sh / C_s,ab = {vsh/ (sc.C_s_ab if hasattr(sc,'C_s_ab') else 0.030*C):.2f}   (paper: 4.6)")
print(f"  v_sh / v_p(model) = {vsh/sc.vp_model:.2f}   (paper: ~4/3 at low B)")
MA  = vsh/sc.vA if hasattr(sc,"vA") else np.nan
print(f"  vsh_model = {sc.vsh_model/C:.4f} c ; vp_model = {sc.vp_model/C:.4f} c")

# --- onset time t*_1 from G(t) ---
tt,G=[],[]
for fr in frames:
    z,uz=io.species_phase(fr,"amb_ions",sc,mass=sc.mi)
    G.append(metrics.reflected_fraction_G(uz*C, vsh if np.isfinite(vsh) else sc.vsh_model))
    tt.append(fr.time)
tt=np.array(tt); G=np.array(G)
tstar,i1=metrics.onset_time_from_G(tt,G)
print("\n=== onset time (max dG/dt) ===")
print(f"  t*_1 = {tstar*sc.wci0:.3f} wci0^-1   (paper t*_1 ~ 1)")

# --- onset position z*_1 from F(z) at t*_1 ---
if i1>=0:
    fr=frames[i1]
    z,uz=io.species_phase(fr,"amb_ions",sc,mass=sc.mi)
    edges=np.asarray(fr.z_edges)
    F,centers=metrics.reflected_profile_F(z*1.0, uz*C, vsh if np.isfinite(vsh) else sc.vsh_model, edges)
    def smooth(a,k=50):
        ker=np.ones(k)/k; return np.convolve(a,ker,mode="same")
    for _ in range(6): F=smooth(F,50)
    zc=centers/sc.rho_i0
    zstar,_=metrics.onset_location_from_F(zc,F)
    print("\n=== onset position (max dF/dz at t*_1) ===")
    print(f"  z*_1 = {zstar:.2f} rho_i0   (paper z*_1 ~ 1);  front z(t*_1) = {zs[i1]/sc.rho_i0:.2f} rho_i0")

# --- first precursor / first shock from criteria flags for reference ---
print(f"\n  rho_i0 = {sc.rho_i0/sc.de:.0f} d_e = {sc.rho_i0*1e3:.3f} mm ; wci0^-1 = {1/sc.wci0*1e9:.2f} ns")

# --- clean-window speed (exclude boundary-contaminated frames, t*wci0 <= 2.25) ---
tw = ts*sc.wci0
half_len = 3600*sc.de   # domain edge
clean = np.isfinite(zs) & (tw <= 2.25) & (zs < 0.94*half_len)
tc, zc2 = ts[clean], zs[clean]
if tc.size>=2:
    vclean=float(np.polyfit(tc, zc2, 1)[0])
    print("\n=== clean-window shock speed (t*wci0<=2.25, front<0.94*edge) ===")
    print(f"  n_clean={tc.size}  t*wci0 in [{tc[0]*sc.wci0:.2f},{tc[-1]*sc.wci0:.2f}]")
    print(f"  v_sh(clean) = {vclean/C:.4f} c = {vclean/sc.Cs_ab:.2f} C_s,ab = {vclean/sc.vp_model:.2f} v_p")
