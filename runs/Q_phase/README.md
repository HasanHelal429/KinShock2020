# Q_phase — quasi-1D (2D x–z) 470 eV runs

**Question (user, 2026-11-18).** Tests A and B (`runs/F_phase/FALSE_SHOCK_AUDIT.md`) put the
false shock on an ion-acoustic instability at the leading edge of the fast-ion beam. That mode
is unstable at 470 eV, where front v_beam/v_te ≈ 4.9, and Landau-damped at 15 keV (≈ 3.0).
1D forces its wave vector along z. **Does the false shock survive when the mode can go
oblique / k‖B?** If it does, it is robust physics of the real-c plasma. If it does not, it
is largely a 1D artifact, and the honest production run is quasi-1D at 470 eV.

**Design.** These are twins of `R1_470eV_s3f8_t2`:
- physics and numerics are the production run's: dz = dx = 0.0375 d_e, dz/λ_D,amb 0.758,
  shape 3, filter 8, collisions, heater, injector;
- cut to 12 d_i0 and t·ω_ci0 = 0.6 (the 1D notch appears after ~0.36);
- 2D x–z, periodic in x, which is along B₀, so k‖B is allowed;
- ppc 4 × 4 = 16 per cell;
- the 2D CFL gives √2 more steps (3,373,524);
- binary: tag `A2D` (`warpx-cda/build_pm_a2d/bin/warpx.2d`, the same acc2d6621 source and
  flags as build A, with `WarpX_DIMS=2`);
- the deck generator's 2D path is `geometry.dims: 2` plus `geometry.transverse_cells`,
  `numerics.ppc_transverse`, a list `max_grid_size`, `diagnostics.particle_fraction` and
  `diagnostics.field_coarsening`. 1D decks are byte-identical (checked on fs_470eV_dz0202_t1).

| run | width | macroparticles | GPUs | est. cost | role |
|---|---|---|---|---|---|
| `q2_470eV_nx8` | 8 cells = 65 nm, narrower than the ~340 nm mode | 9.2×10⁶ | 1 (shared) | ~10–20 GPU-h | option 0: must reproduce 1D; validates the 2D deck (watch PN for the injector's slab density) |
| `q2_470eV_nx128` | 128 cells = 1.04 µm (~3 mode λ, 0.43 d_e0) | 1.47×10⁸ | 4 (regular, one node, 4 z boxes) | ~35–85 h wall, ~35–85 node-h | option A |

Launch (from the main checkout, after it has these files):

    perlmutter/submit_chain.sh runs/Q_phase/q2_470eV_nx8   --binary A2D --segments 1
    perlmutter/submit_chain.sh runs/Q_phase/q2_470eV_nx128 --binary A2D --gpus 4 --qos regular --segments 3

## Verdicts

| run | status | false shock (by eye) |
|---|---|---|
| `q2_470eV_nx8` | staged | — |
| `q2_470eV_nx128` | staged | — |
