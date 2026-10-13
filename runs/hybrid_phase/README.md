# hybrid_phase — kinetic ions, fluid electrons (H-PICShock decks)

These runs use **H-PICShock's config schema** (github.com/HasanHelal429/H-PICShock), not
kinshock's: `kinshock.load` cannot read them, and `make_inputs.py`, `submit_chain.sh` and
`find_runs`-based tools skip or fail on them. Decks are generated with H-PICShock's
`scripts/make_inputs.py <run_dir>` (commit `59139fe` reproduces `H3_470eV_dense`'s deck
key-for-key apart from diagnostics) and run on the commit-A WarpX binary, which carries the
hybrid solver and the electron-fluid heater. Launch with `perlmutter/job_hybrid.sbatch`.
Render with `scripts/hybrid_phase_movie.py`, which reads scales from `warpx_used_inputs`.

| run | what | false shock (by eye) |
|---|---|---|
| `H3_470eV_eheat_cal` | heater-calibration point, θ = 1.993e-4, 80 d_i0, dz 0.1 d_i0 | — |
| `H3_470eV_dense` | same at dz 0.025 d_i0 (dz/4, dt/4); piston front 6.2–6.4 d_i0/t* (kinetic: 12.1–13.7) | **no** (2026-10-12) |
| `hy_470eV_t2` | hybrid twin of `R1_470eV_s3f8_t2` (34.2 d_i0, to t·ω_ci0 2.0), θ = 6.3e-4 (estimate for ~12.5 d_i0/t*) | — |
| `hy_470eV_t2_lo` / `_hi` | same, θ = 5.0e-4 / 7.5e-4 (brackets) | — |

The `hy_*` twins change only box, clock, diagnostic cadence and heater amplitude relative to
`H3_470eV_dense`. Heater sizing: T_ss ∝ θ^1.21 and v_front ∝ √T give v ∝ θ^0.6. Matching the
kinetic piston speed puts the fluid slab near 1.5–2.5 keV rather than 470 eV; the closure
cannot match both, and the drive speed sets the shock's M_A.
