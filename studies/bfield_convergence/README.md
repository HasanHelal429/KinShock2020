# B-field fluctuation: physical vs numerical

**Question.** The runs show large-amplitude B oscillations. Which are real
collisionless-shock turbulence, and which are a PIC artifact?

**Background (see `RESULTS.md`).** The reflected-ion foot of a supercritical
perpendicular shock is *genuinely* turbulent — so B fluctuations where reflected
ions live (within ~0.3 rho_i of the front) are expected to be physical. Separately,
the Debye length is under-resolved (dz/lambda_D ~ 6.7 at dz=0.3 d_e; the finite-grid
heating threshold is ~pi), which can grow *grid-scale numerical* field noise.

**Method.** `run_variants.sh` runs the same physics with one numerical knob changed
per variant, long enough (t*wci ~ 0.56) for the far-upstream fluctuation to develop:

- `baseline`  — dz=0.3 d_e, filter_npass=1, particle_shape=2
- `filt8`     — heavier current filter (`warpx.filter_npass_each_dir = 8`)
- `shape3`    — cubic particle shape (`algo.particle_shape = 3`)
- `finer_dz`  — resolve the Debye length (`amr.n_cell` x2 -> dz/2; dt halves, so 2x steps)

Then `analyze.py` (and `scripts/bfield_diagnostic.py`) compare RMS(dBx) in the cold
far-upstream zone across variants. **Filter/shape/finer-dz suppress grid-scale numerical
noise but not physical waves at resolved scales**, so:

- fluctuation **collapses** under the knobs  => numerical (grid instability)
- fluctuation **invariant**                  => physical

**Run it:**

```bash
WARPX=/path/to/warpx.1d studies/bfield_convergence/run_variants.sh \
    runs/R1_phase/R1_core_half/inputs_kinshock_R1_core_half     # -> ./scratch/{baseline,filt8,shape3,finer_dz}
python studies/bfield_convergence/analyze.py studies/bfield_convergence/scratch
python scripts/bfield_diagnostic.py studies/bfield_convergence/scratch/finer_dz --twci 0.5
```

**Status / findings.** A short (t*wci~0.11) pass already confirmed the *near-shock/foot*
turbulence is **physical** (spectrum invariant to filter/shape/dz). The long-run test of
the *far-upstream* component is the decisive one; results are logged in `RESULTS.md`
under the B-field entry. The particle-response test (`bfield_diagnostic.py`) is
independent evidence: a dB/B~1 fluctuation that scatters neither ions nor electrons
(both at the t=0 thermal floor) is not exchanging energy with the plasma.
