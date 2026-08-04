# Run `R1_warm`

> Full R1 (Schaeffer 2020 Table I, M_A~14 perpendicular) with warm ablative piston ions (T_i = T_e,ab -> v_th,i = C_s,ab) and theta_e recalibrated to 0.062 to hold M_A ~ 14. Reproduces the paper's broad piston-ion phase space (Fig. 7) that cold-ion R1_half missed.

**tier** Full · **deck** `inputs_kinshock_R1_warm` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I

## Notes

<!-- prose:begin -->
**The current full reference run.** Warm ablative piston ions (`theta_i_init =
theta_e_heat/mass_ratio`, so T_i = T_e,ab and v_th,i = C_s,ab) reproduce the paper's broad
Fig. 7 row-2 ion ridge that cold-ion `R1_half` missed. `theta_e_heat` was then recalibrated
0.092 -> 0.062 -> 0.078 by eye so the *settled* M_A lands near 14 (RESULTS 2026-07-24); the
`M_A` row below is the *model* value from `model.vsh_over_Csab`, not the by-eye fit, which is
why it reads 12.8.

**Known deviations from Table I** (see the comparison table, and CLAUDE.md):
- `density_over_n0: 0.01` is **25% high** — Table I's ambient is 0.008 n_e,ab (its 0.01 is in
  code units where n_e,ab = 1.25). Because B0 is now primary this no longer moves the clock,
  but it still makes d_i0 = 10 d_i,ab instead of 11.2 and v_A/c = 0.01 instead of 0.0089.
- The gyroperiod is 27.93 t_ab vs the paper's 33.9 — 1.218x short, from that density (1.118x)
  and the theta_e recal (1.086x). Every `t*wci0` number here is biased late by that factor;
  it is the leading candidate for the onset t*_1 = 1.35-1.41 vs the paper's ~1.
- Domain is 7500 d_e vs the paper's 9000 d_e,ab; 250k steps vs 400k.
- beta rows are 2x the paper's because `units.py` uses beta = 2*mu0*n*T/B^2 while Table I
  tabulates mu0*n*T/B^2. Convention, not physics.

Shock kinematics come from `shock_fit.yaml` (by eye), never auto-detection.
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/R1_warm`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 1.0e18 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2.5 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.078 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.001 | `plasma.piston.theta_e_init` |
| piston init theta_i | 0.00078 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.01 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.002 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 0.00320753 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 7500 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.3 | `geometry.dz_over_de` |
| boundary lo / hi | lo symmetry, hi open | `geometry.boundary` |
| CFL | 0.75 | `numerics.cfl` |
| particle shape | 2 | `numerics.particle_shape` |
| max_step | 250000 | `numerics.max_step` |
| ppc | piston 100, ambient 100 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 40 | `operators.injector.tau_over_wpe_inv` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 5.64146e+10 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 0.00531409 | c/omega_pe |
| d_i,ab [m] | 0.0531409 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.0279285 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 6.34689e-09 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 1e+16 | density_over_n0 * n0 |
| d_i0 [m] | 0.531409 | c/omega_pi(n_amb) |
| d_i0 / d_e | 100 | derived |
| **omega_ci0 [rad/s]** | 5.64146e+06 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 1.77259e-07 | 1/omega_ci0 |
| **v_A / c** | 0.01 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 1040 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.104 | config model.vp_over_c |
| v_sh / c (model) | 0.128471 | model.vsh_over_Csab * C_s,ab |
| M_A | 12.8471 | v_sh/v_A |
| M_ms | 11.7277 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 780 | mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.2 | mu0*n_amb*T_0/B0^2 |
| dz [m] | 0.00159423 | dz_over_de * d_e |
| dt [s] | 3.98833e-12 | CFL-limited |
| dt*omega_pe | 0.225 | dt * omega_pe |
| n_cell | 25000 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 44444.4 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 5.625 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 39857.9 | theta_e,ab * m_e c^2 |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.0279285 | 0.03 | 0.931x (-6.9%) ~ | theta_e_heat recalibrated off the paper's 0.092 |
| piston speed v_p / c | 0.104 | 0.104 | 1.000x (+0.0%) ok |  |
| Alfven Mach v_sh/v_A | 12.8471 | 14 | 0.918x (-8.2%) ~ | model M_A from model.vsh_over_Csab; the by-eye settled value is in shock_fit.yaml |
| magnetosonic Mach | 11.7277 | 13 | 0.902x (-9.8%) ~ |  |
| ablation beta, mu0*n*T/B^2 (Table I 1150) | 780 | 1150 | 0.678x (-32.2%) **OFF** | n_amb 0.01 vs Table I's 0.008, and/or theta_e off 0.092 |
| upstream beta, mu0*n*T/B^2 (Table I 0.2) | 0.2 | 0.2 | 1.000x (-0.0%) ok |  |
| d_i0 / d_i,ab | 10 | 11.18 | 0.894x (-10.6%) ~ | n_amb is 0.01 n0; Table I is 0.008 n0 |
| gyroperiod in ablation times | 27.9285 | 33.9 | 0.824x (-17.6%) ~ | n_amb 25% high (1.118x) and/or theta_e recal (1.086x) |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_R1_warm` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
