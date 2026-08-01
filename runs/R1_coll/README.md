# Run `R1_coll`

> Collisional twin of R1_warm: identical dimensionless setup (M_A~14 perpendicular, warm ablative piston ions) pinned to an absolute ambient density n_e0 = 1e18 cm^-3, with pairwise Coulomb collisions on all species pairs tuned to the paper's Table I collisionality lambda_ab = omega_ce,ab/nu_ei,ab = 20 (mfp = 559 d_e,ab = 5.6 d_i0).

**tier** Full · **deck** `inputs_kinshock_R1_coll` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I + Fig. 13

## Notes

<!-- prose:begin -->
**R1_warm's collisional twin.** Same dimensionless setup, but n0 pinned to 1e18 cm^-3 (1e20
m^-3 -> B0 = 32.1 T, 10x R1_warm's, because B0 in tesla scales with the density it was chosen
at) and pairwise Coulomb collisions on. Built to test the paper's Fig. 13 claim that
collisional and collisionless perpendicular shock formation look the same. It does: R1_coll
reproduces R1_warm's compression, front speed, ramp steepness, piston separation and reflected
-ion fraction (peak B_perp 12% lower).

**Its collisionality target is wrong, twice over** (RESULTS 2026-07-31):
- `lambda_ab` was set from a mis-definition. The paper says outright that
  lambda_ab = omega_ce,ab/nu_ei,ab = **mfp/d_e,ab** with omega_ce,ab at B_ab = sqrt(mu0 n T),
  not at B0. So Table I's 20 means mfp = 20 d_e,ab; this run targets 559 d_e,ab, i.e. it is
  ~28x *less* collisional than the paper.
- Even so it sits at upstream lambda_ii/d_i0 = 0.52 vs Table I's 350, so **criterion 2
  ("collisionless") fails in all 51 frames**. lambda_ab and lambda_mfp/d_i0 cannot both be
  matched at mu_p = 100 — the paper concedes this in §II and names the physical mass ratio,
  not a reduced c, as the only fix.

So this is *not* "the paper's Table I collisionality"; it is a much more collisional run that
behaves the same anyway — which tests Fig. 13 harder than the paper did.
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/R1_coll`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 1.0e26 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2.5 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.078 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.001 | `plasma.piston.theta_e_init` |
| piston init theta_i | 0.00078 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.01 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.002 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 32.0753 | `field.B0_tesla` |
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
| collisions target | quantity coulomb_log, value 7713.3 | `collisions.target` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 5.64146e+14 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 5.31409e-07 | c/omega_pe |
| d_i,ab [m] | 5.31409e-06 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.0279285 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 6.34689e-13 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 1e+24 | density_over_n0 * n0 |
| d_i0 [m] | 5.31409e-05 | c/omega_pi(n_amb) |
| d_i0 / d_e | 100 | derived |
| **omega_ci0 [rad/s]** | 5.64146e+10 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 1.77259e-11 | 1/omega_ci0 |
| **v_A / c** | 0.01 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 1040 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.104 | config model.vp_over_c |
| v_sh / c (model) | 0.128471 | model.vsh_over_Csab * C_s,ab |
| M_A | 12.8471 | v_sh/v_A |
| M_ms | 11.7277 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 1560 | 2*mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.4 | 2*mu0*n_amb*T_0/B0^2 |
| dz [m] | 1.59423e-07 | dz_over_de * d_e |
| dt [s] | 3.98833e-16 | CFL-limited |
| dt*omega_pe | 0.225 | dt * omega_pe |
| n_cell | 25000 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 44444.4 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 5.625 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 39857.9 | theta_e,ab * m_e c^2 |
| lnLambda (used) | 7713.3 | units.coulomb_log_for(collisions.target) |
| lnLambda (physical) | 11.5672 | NRL at (n0, T_e,ab) |
| nu_ei,ab [1/s] | 2.82073e+11 | NRL electron-ion |
| nu_ei*dt | 0.0001125 | must be << 1 |
| mfp_ei,ab / d_e | 558.57 | v_te/nu_ei / d_e |
| lambda_ab | 558.57 | omega_ce,ab/nu_ei,ab = mfp/d_e,ab |
| mfp_ii,amb / d_i0 | 0.519572 | upstream ion-ion — what criterion 2 tests |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.0279285 | 0.03 | 0.931x (-6.9%) ~ | theta_e_heat recalibrated off the paper's 0.092 |
| piston speed v_p / c | 0.104 | 0.104 | 1.000x (+0.0%) ok |  |
| Alfven Mach v_sh/v_A | 12.8471 | 14 | 0.918x (-8.2%) ~ | model M_A from model.vsh_over_Csab; the by-eye settled value is in shock_fit.yaml |
| magnetosonic Mach | 11.7277 | 13 | 0.902x (-9.8%) ~ |  |
| ablation beta (Table I 1150, x2 convention) | 1560 | 2300 | 0.678x (-32.2%) **OFF** | n_amb 0.01 vs Table I's 0.008, and/or theta_e off 0.092 |
| upstream beta (Table I 0.2, x2 convention) | 0.4 | 0.4 | 1.000x (-0.0%) ok |  |
| d_i0 / d_i,ab | 10 | 11.18 | 0.894x (-10.6%) ~ | n_amb is 0.01 n0; Table I is 0.008 n0 |
| gyroperiod in ablation times | 27.9285 | 33.9 | 0.824x (-17.6%) ~ | n_amb 25% high (1.118x) and/or theta_e recal (1.086x) |
| collisionality mfp/d_e,ab | 558.57 | 20 | 27.928x (+2692.8%) **OFF** |  |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_R1_coll` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
