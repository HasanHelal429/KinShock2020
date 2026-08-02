# Run `R1_paper_phys`

> R1_paper with the physical (NRL) lnLambda = 10.836 in place of the 1.223e5 dial that Table I's lambda_ab = 20 requires, truncated to t*wci0 = 0.389 as a quick test. Identical n0, B0, geometry, dz, dt, ppc and operators, so collisionality is the only variable; plotfiles land on R1_paper's own frame steps, making it a matched control. Expected: effectively collisionless (mfp_ei,ab = 2.258e5 d_e,ab, 25x the domain; nu_ei*dt = 3.0e-7), with upstream mfp_ii/d_i0 = 168.8 vs Table I's 350 -- so criterion 2 should flip FAIL -> PASS while lambda_ab misses 20 by 1.1e4x.

**tier** Full · **deck** `inputs_kinshock_R1_paper_phys` · **paper** Schaeffer et al., Phys. Plasmas 27, 042901 (2020), Table I + Sec. II

## Notes

<!-- prose:begin -->
_No hand-written notes yet. Anything written between the
`prose:begin` / `prose:end` markers survives regeneration — put the run's story here
(why it exists, what it showed, what to distrust)._
<!-- prose:end -->

## Primaries — `config.yaml` is the only source of truth

Edit these; never the deck. Regenerate with `python scripts/make_inputs.py runs/R1_paper_phys`.

| Quantity | Value | config key |
|---|---|---|
| reference density n0 [m^-3] | 6.0e26 | `reference.n0` |
| mass ratio m_i/m_e | 100 | `reference.mass_ratio` |
| charge state Z | 1 | `reference.charge_state` |
| target density / n0 | 2 | `plasma.piston.density_over_n0` |
| heater theta_e,ab | 0.092 | `plasma.piston.theta_e_heat` |
| piston init theta_e | 0.002 | `plasma.piston.theta_e_init` |
| piston init theta_i | 2e-05 | `plasma.piston.theta_i_init` |
| ambient density / n0 | 0.008 | `plasma.ambient.density_over_n0` |
| ambient theta_0 | 0.002 | `plasma.ambient.theta_0` |
| field orientation | perpendicular | `field.orientation` |
| **B0 [T]** | 70.2734 | `field.B0_tesla` |
| dims | 1 | `geometry.dims` |
| layout | one_sided | `geometry.layout` |
| slab halfwidth [d_i] | 2 | `geometry.slab_halfwidth_di` |
| domain halfwidth [d_e] | 9000 | `geometry.domain_halfwidth_de` |
| dz [d_e] | 0.3 | `geometry.dz_over_de` |
| boundary lo / hi | lo symmetry, hi open | `geometry.boundary` |
| CFL | 0.75 | `numerics.cfl` |
| particle shape | 2 | `numerics.particle_shape` |
| max_step | 19344 | `numerics.max_step` |
| ppc | piston 25, ambient 25 | `numerics.ppc` |
| heater intervals | 20 | `operators.heater.intervals` |
| injector intervals | 20 | `operators.injector.intervals` |
| injector tau [1/wpe] | 40 | `operators.injector.tau_over_wpe_inv` |
| collisions target | quantity coulomb_log, value physical | `collisions.target` |

## Derived — computed by `units.derive`, not stored

| Quantity | Value | From |
|---|---|---|
| omega_pe [rad/s] | 1.38187e+15 | sqrt(n0 q^2/(eps0 m_e)) |
| d_e [m] | 2.16947e-07 | c/omega_pe |
| d_i,ab [m] | 2.16947e-06 | d_e*sqrt(m_i/m_e) |
| C_s,ab / c | 0.0303315 | sqrt(theta_e,ab/mass_ratio) |
| t_ab [s] | 2.38583e-13 | d_i,ab / C_s,ab |
| n_amb [m^-3] | 4.8e+24 | density_over_n0 * n0 |
| d_i0 [m] | 2.42554e-05 | c/omega_pi(n_amb) |
| d_i0 / d_e | 111.803 | derived |
| **omega_ci0 [rad/s]** | 1.23598e+11 | **q_e*B0/m_i — B0 only, no n_amb** |
| 1/omega_ci0 [s] | 8.09073e-12 | 1/omega_ci0 |
| **v_A / c** | 0.01 | **B0/sqrt(mu0*n_amb*m_i) — DERIVED** |
| rho_i0 / d_e | 1162.76 | v_p/omega_ci0 / d_e |
| v_p / c (model) | 0.104 | config model.vp_over_c |
| v_sh / c (model) | 0.139525 | model.vsh_over_Csab * C_s,ab |
| M_A | 13.9525 | v_sh/v_A |
| M_ms | 12.7368 | v_sh/sqrt(v_A^2+C_s0^2) |
| beta_ab | 2300 | 2*mu0*n0*T_e,ab/B0^2 |
| beta_0 | 0.4 | 2*mu0*n_amb*T_0/B0^2 |
| dz [m] | 6.50841e-08 | dz_over_de * d_e |
| dt [s] | 1.62823e-16 | CFL-limited |
| dt*omega_pe | 0.225 | dt * omega_pe |
| n_cell | 30000 | domain / dz (halved when one_sided) |
| steps per 1/omega_ci0 | 49690.4 | 1/(omega_ci0*dt) |
| run length [1/omega_ci0] | 0.38929 | max_step * dt * omega_ci0 |
| T_e,ab [eV] | 47011.9 | theta_e,ab * m_e c^2 |
| lnLambda (used) | 10.8364 | units.coulomb_log_for(collisions.target) |
| lnLambda (physical) | 10.8364 | NRL at (n0, T_e,ab) |
| nu_ei,ab [1/s] | 1.85617e+09 | NRL electron-ion |
| nu_ei*dt | 3.02227e-07 | must be << 1 |
| mfp_ei,ab / d_e | 225810 | v_te/nu_ei / d_e |
| lambda_ab | 225810 | omega_ce,ab/nu_ei,ab = mfp/d_e,ab |
| mfp_ii,amb / d_i0 | 168.803 | upstream ion-ion — what criterion 2 tests |

## vs Schaeffer 2020 Table I

| Quantity | This run | Table I | ratio | known cause if off |
|---|---|---|---|---|
| mass ratio m_i/m_e | 100 | 100 | 1.000x (+0.0%) ok |  |
| C_s,ab / c | 0.0303315 | 0.03 | 1.011x (+1.1%) ok |  |
| piston speed v_p / c | 0.104 | 0.104 | 1.000x (+0.0%) ok |  |
| Alfven Mach v_sh/v_A | 13.9525 | 14 | 0.997x (-0.3%) ok |  |
| magnetosonic Mach | 12.7368 | 13 | 0.980x (-2.0%) ok |  |
| ablation beta (Table I 1150, x2 convention) | 2300 | 2300 | 1.000x (-0.0%) ok |  |
| upstream beta (Table I 0.2, x2 convention) | 0.4 | 0.4 | 1.000x (-0.0%) ok |  |
| d_i0 / d_i,ab | 11.1803 | 11.18 | 1.000x (+0.0%) ok |  |
| gyroperiod in ablation times | 33.9116 | 33.9 | 1.000x (+0.0%) ok |  |
| collisionality mfp/d_e,ab | 225810 | 20 | 11290.483x (+1128948.3%) **OFF** |  |

## Files

| Path | What |
|---|---|
| `config.yaml` | the primaries (tracked) |
| `inputs_kinshock_R1_paper_phys` | generated deck (tracked, **never hand-edit**) |
| `warpx_used_inputs` | what WarpX actually ran (tracked) |
| `shock_fit.yaml` | by-eye v_sh + front, from `tune_shock.py` |
| `diags/`, `*.log` | output (gitignored, regenerable) |

---

_Tables generated by `scripts/make_run_readme.py` — edit `config.yaml`, not the tables. Prose between the `prose:` markers is preserved._
