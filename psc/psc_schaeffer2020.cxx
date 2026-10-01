// ======================================================================
// psc_schaeffer2020 -- the reference run of Schaeffer et al., Phys. Plasmas 27,
// 042901 (2020), written FROM THE PAPER, not adapted from a PSC test case.
//
// Every physics number below cites its source in the paper. PSC's own
// psc_flatfoil_yz.cxx supplied only the API (how a heating operator, foil
// injector and collision object are wired together); where flatfoil's VALUES
// differ from the paper, the paper wins -- most importantly the target
// half-width (flatfoil 1 d_i, paper 2 d_i,ab), which matters because the
// heating rate is proportional to 1/width.
//
// UNITS: PSC code units. Density unit n = 1 (the paper's "n_e,ab = 1", Sec. II:
// "1000 [macroparticles per cell] at n_e,ab = 1"); length unit d_e at n = 1;
// time unit 1/omega_pe at n = 1; temperature m_e c^2; B in sqrt(mu0 n m_e c^2).
// The Table I ablation density n_e,ab = 1.25 is EMERGENT -- the plasma that a
// 2.5 target rarefies to -- and is NOT an input. The inputs are the target
// (2.5) and the ambient (0.01) densities, both stated in Sec. II / Table I.
//
// Usage:  psc_schaeffer2020 params.txt   (every value has a paper default)
// ======================================================================

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "output_fields.hxx"
#include "psc_config.hxx"
#include "input_params.hxx"

#include "../libpsc/psc_heating/psc_heating_impl.hxx"
#include "heating_spot_foil.hxx"

#include <fstream>

// Two species, as in the paper's single-species runs: H ions and ONE electron
// population. (flatfoil adds a hot-electron population; the paper has none.)
// The last kind is the neutralizing population.
enum
{
  MY_ELECTRON,
  MY_ION,
  N_MY_KINDS,
};

namespace
{
struct Sch2020
{
  // -- plasma (Table I, Sec. II.2)
  double mass_ratio;  // mu_p = 100                       Table I, Sec. II.2
  double Zi;          // Z = 1, H                          Table I
  double Te_ab;       // 0.092 m_e c^2                     Table I (heating target)
  double target_n;    // 2.5                               Sec. II.2
  double target_T0;   // initial target T (inferred)       see params.txt note
  double amb_n;       // n_e0 = 0.01                       Table I
  double amb_T;       // T_0 = T_e0 = T_i0 = 0.002 m_e c^2 Table I, Sec. II.1
  double B0;          // 0.01 sqrt(m_e c^2), along y here  Table I
  double lambda_ab;   // 20                                Table I, Sec. II.2
  double target_hw_di;// 2 d_i,ab (0 < z <= 2 d_i,ab)       Sec. II.2
  // -- grid / numerics
  int ny, nz;         // 12 x (see params.txt)             Sec. II.2
  double Ly, Lz;      // 5 x (see params.txt) d_e          Sec. II.2
  int npz;            // patches along z (decomposition)
  int nicell;         // 1000 at n = 1                     Sec. II.2
  double cfl;         // 0.75 (reproduces 400,000 steps)   derived, see params.txt
  int nmax;           // 400,000                           Sec. II.2
  // -- operator cadences (not given in the paper; PSC flatfoil values)
  int inject_interval, heating_interval, collision_interval;
  double inject_tau;
  // -- output
  int field_every, prt_every, energy_every;
  // -- derived
  double d_i;
} g;

std::string params_path;
PscParams psc_params;
} // namespace

using Dim = dim_yz; // heating_spot_foil supports dim_yz/dim_xyz only; B along y,
                    // ablation along z -- the paper's x-z plane with B along x,
                    // relabelled x->y. Physically identical.
#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecSingle<Dim>;
#endif
using Writer = WriterDefault;
using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
// Every 20th particle: ~1e8 particles late in the run would otherwise make each
// phase-space dump tens of GB. Sparse ambient is thinned in proportion.
using OutputParticles =
  OutputParticlesHdf5<Mparticles, ParticleSelectorEveryNth<20>>;
// Host density moment, also on GPU. PSC calls the injection hook after the push
// but before particles are re-sorted into blocks, and both CUDA moment kernels
// read the stale block offsets: they return ~0 (domain sum 680 vs 9550 on CPU)
// and the target is refilled from empty at every injection (flatfoil's CUDA
// selector has this problem). The host moment is correct but copies every
// particle off the GPU, so injection runs every 200 steps (0.11 t_ab) rather
// than flatfoil's 20, where it was 78% of the step time.
using Moment_n = Moment_n_1st<MfieldsSingle::Storage, Dim>;
using Heating = typename HeatingSelector<Mparticles>::Heating;

// ----------------------------------------------------------------------
// the target: 0 < |z| <= 2 d_i,ab, uniform in y (quasi-1D, Sec. II.1)

struct Target
{
  double zh, n, T;
  bool is_inside(const double crd[3]) const { return std::abs(crd[2]) <= zh; }
};

void setupParameters(int argc, char** argv)
{
  if (argc != 2) {
    LOG_ERROR("Usage: %s path/to/params.txt\n", argv[0]);
  }
  params_path = argv[1];
  InputParams p(params_path);

  g.mass_ratio   = p.getOrDefault<double>("mass_ratio", 100.);
  g.Zi           = p.getOrDefault<double>("Zi", 1.);
  g.Te_ab        = p.getOrDefault<double>("Te_ab", 0.092);
  g.target_n     = p.getOrDefault<double>("target_n", 2.5);
  g.target_T0    = p.getOrDefault<double>("target_T0", 0.002);
  g.amb_n        = p.getOrDefault<double>("amb_n", 0.01);
  g.amb_T        = p.getOrDefault<double>("amb_T", 0.002);
  g.B0           = p.getOrDefault<double>("B0", 0.01);
  g.lambda_ab    = p.getOrDefault<double>("lambda_ab", 20.);
  g.target_hw_di = p.getOrDefault<double>("target_hw_di", 2.);

  g.ny     = p.getOrDefault<int>("ny", 12);
  g.nz     = p.getOrDefault<int>("nz", 30000);
  g.Ly     = p.getOrDefault<double>("Ly", 5.);
  g.Lz     = p.getOrDefault<double>("Lz", 9000.);
  g.npz    = p.getOrDefault<int>("npz", 150); // nz/npz must be a multiple of 4 (CUDA BS144)
  g.nicell = p.getOrDefault<int>("nicell", 1000);
  g.cfl    = p.getOrDefault<double>("cfl", 0.75);
  g.nmax   = p.getOrDefault<int>("nmax", 400000);

  g.inject_interval    = p.getOrDefault<int>("inject_interval", 200);
  g.heating_interval   = p.getOrDefault<int>("heating_interval", 20);
  g.collision_interval = p.getOrDefault<int>("collision_interval", 10);
  g.inject_tau         = p.getOrDefault<double>("inject_tau", 40.);

  g.field_every  = p.getOrDefault<int>("field_every", 800);
  g.prt_every    = p.getOrDefault<int>("prt_every", 8000);
  // 0 = off. PSC main's DiagEnergiesField builds an empty view for an
  // invariant x (ibn[0] = 0) and throws; flatfoil disables it the same way.
  g.energy_every = p.getOrDefault<int>("energy_every", 0);

  psc_params.nmax = g.nmax;
  psc_params.cfl = g.cfl;
  psc_params.stats_every = p.getOrDefault<int>("stats_every", 100);
  psc_params.write_checkpoint_every_step =
    p.getOrDefault<int>("checkpoint_every", 0) /* needs ADIOS2 */;
}

Grid_t* setupGrid()
{
  Grid_t::Real3 LL = {1., g.Ly, g.Lz};
  Int3 gdims = {1, g.ny, g.nz};
  Int3 np = {1, 1, g.npz};
  // PSC centres the domain on the origin: z in [-Lz/2, Lz/2]. The target sits
  // at z = 0 and ablates in BOTH directions, as the paper states.
  Grid_t::Domain domain{gdims, LL, -.5 * LL, np};

  // periodic in both x and z (Sec. II.2)
  psc::grid::BC bc{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                   {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  Grid_t::Kinds kinds(N_MY_KINDS);
  kinds[MY_ION] = {g.Zi, g.mass_ratio * g.Zi, "i"};
  kinds[MY_ELECTRON] = {-1., 1., "e"};
  g.d_i = std::sqrt(kinds[MY_ION].m / kinds[MY_ION].q); // = 10 at mu = 100

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = g.nicell;
  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  // --- print the derived scales so they can be checked against Table I
  double Cs_ab = std::sqrt(g.Zi * g.Te_ab / g.mass_ratio);
  double t_ab = g.d_i / Cs_ab;
  double wci0 = g.Zi * g.B0 / (g.mass_ratio * g.Zi);
  double d_i0 = std::sqrt(g.mass_ratio / (g.Zi * g.amb_n));
  mpi_printf(MPI_COMM_WORLD, "=== schaeffer2020 derived scales (code units) ===\n");
  mpi_printf(MPI_COMM_WORLD, "d_e = 1  d_i,ab = %g  d_i0 = %g\n", g.d_i, d_i0);
  mpi_printf(MPI_COMM_WORLD, "dz = %g  dy = %g  dt = %g\n", g.Lz / g.nz,
             g.Ly / g.ny, dt);
  mpi_printf(MPI_COMM_WORLD, "C_s,ab/c = %g  c/C_s,ab = %g  t_ab = %g\n", Cs_ab,
             1. / Cs_ab, t_ab);
  mpi_printf(MPI_COMM_WORLD, "1/w_ci0 = %g  = %g t_ab\n", 1. / wci0,
             1. / wci0 / t_ab);
  mpi_printf(MPI_COMM_WORLD, "lambda_D,amb = %g  dz/lambda_D,amb = %g\n",
             std::sqrt(g.amb_T / g.amb_n), (g.Lz / g.nz) /
             std::sqrt(g.amb_T / g.amb_n));
  mpi_printf(MPI_COMM_WORLD, "run: nmax = %d -> t = %g = %g t_ab = %g /w_ci0\n",
             g.nmax, g.nmax * dt, g.nmax * dt / t_ab, g.nmax * dt * wci0);

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();
  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

void run(int argc, char** argv)
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up schaeffer2020...\n");
  setupParameters(argc, argv);
  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  // keep the exact params used next to the output (provenance)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      std::ifstream src(params_path, std::ios::binary);
      std::ofstream dst("params_record.txt", std::ios::binary);
      dst << src.rdbuf();
    }
  }

  Mparticles mprts(grid);
  MfieldsState mflds(grid);

  psc_params.balance_interval = 500;
  Balance balance{3};
  psc_params.sort_interval = 10;

  // -- Collisions: Takizuka-Abe (Sec. II.1), lambda_ab = w_ce,ab/nu_ei,ab = 20.
  // PSC's operator scales nu internally by n/v^3 ~ n/T^1.5, and the paper
  // evaluates lambda_ab "at T_ab and n_e,ab = 1", where w_ce,ab = B_ab = sqrt(T)
  // so nu_ei = sqrt(T)/lambda_ab. Hence nu ~ T^2/lambda_ab -- the same form PSC's
  // flatfoil uses, with the paper's T_ab and lambda_ab.
  double collision_nu = 3.76 * std::pow(g.Te_ab, 2.) / g.Zi / g.lambda_ab;
  Collision collision{grid, g.collision_interval, collision_nu};
  mpi_printf(grid.comm(), "collision_nu = %g (lambda_ab = %g)\n", collision_nu,
             g.lambda_ab);

  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 0;
  checks_params.gauss.check_interval = 1000;
  checks_params.gauss.err_threshold = 1e-4;
  checks_params.gauss.print_max_err_always = true;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  psc_params.marder_interval = 100;
  Marder marder(grid, 0.9, 3, false);

  // -- output
  OutputFields<MfieldsState, Mparticles, Writer> out_fields;
  out_fields.pfield.out_interval = g.field_every;
  OutputMoments<MfieldsState, Mparticles, Dim, Writer> out_moments{out_fields};
  OutputParticlesParams outp_params{};
  outp_params.every_step = g.prt_every;
  OutputParticles outp{grid, outp_params};
  DiagEnergies<Mparticles, MfieldsState> oute{grid.comm(), g.energy_every};

  // -- the target and its ablation
  Target target{g.target_hw_di * g.d_i, g.target_n, g.target_T0};

  // Heating operator over the WHOLE target, uniform in y (rH = 0 -> "uniform
  // heating, not a spot"): the paper's "uniform driving conditions in the
  // transverse direction". Electrons only, heated to T_e,ab; ions are not
  // heated, they are accelerated by the ambipolar field (Sec. II.1).
  HeatingSpotFoilParams hp{};
  hp.zl = -target.zh;
  hp.zh = target.zh;
  hp.xc = 0.;
  hp.yc = 0.;
  hp.rH = 0.;
  hp.T[MY_ELECTRON] = g.Te_ab;
  hp.T[MY_ION] = 0.;
  hp.Mi = grid.kinds[MY_ION].m;
  hp.n_kinds = N_MY_KINDS;
  HeatingSpotFoil<Dim> heating_spot{grid, hp};
  auto& heating = *new Heating{grid, g.heating_interval, heating_spot};

  SetupParticles<Mparticles> setup_particles(grid);
  setup_particles.fractional_n_particles_per_cell = true;
  setup_particles.neutralizing_population = MY_ION;

  // "the heating operator continuously added new particles to the target to
  // maintain this target density" (Sec. II.2): relax toward target_n.
  double inject_fac = (g.inject_interval * grid.dt / g.inject_tau) /
                      (1. + g.inject_interval * grid.dt / g.inject_tau);

  auto lf_inject_heat = [&](Mparticles& mprts, MfieldsState& mflds) {
    const Grid_t& grid = mprts.grid();
    auto timestep = grid.timestep();

    if (g.inject_interval > 0 && timestep % g.inject_interval == 0) {
      Moment_n moment_n{grid};
      auto d_n = psc::mflds::interior(grid, moment_n(mprts));
      auto&& h_n = gt::host_mirror(d_n);
      gt::copy(d_n, h_n);
      setup_particles.setupParticles(
        mprts,
        [&](int kind, Double3 pos, int p, Int3 idx, psc_particle_npt& npt) {
          double crd[3] = {pos[0], pos[1], pos[2]};
          if (!target.is_inside(crd)) {
            npt.n = 0.;
            return;
          }
          npt.n = target.n - h_n(idx[0], idx[1], idx[2], kind, p);
          npt.T[0] = npt.T[1] = npt.T[2] =
            (kind == MY_ION) ? g.amb_T : target.T;
          if (npt.n < 0) {
            npt.n = 0;
          }
          npt.n *= inject_fac;
        });
    }
    if (g.heating_interval > 0 && timestep % g.heating_interval == 0) {
      heating(mprts);
    }
  };

  // -- initial condition: ambient everywhere, replaced by the target inside it.
  // Ambient plasma AND field occupy |z| > 2 d_i,ab only (Sec. II.2).
  partitionAndSetupParticles(
    setup_particles, balance, grid_ptr, mprts,
    // 5-arg form on purpose: PSC's 3-arg InitNptFunc adapter captures its
    // argument by reference and dangles (setup_particles.hxx), segfaulting here.
    [&](int kind, Double3 crd, int /*p*/, Int3 /*idx*/, psc_particle_npt& npt) {
      double c[3] = {crd[0], crd[1], crd[2]};
      if (target.is_inside(c)) {
        npt.n = target.n;
        npt.T[0] = npt.T[1] = npt.T[2] =
          (kind == MY_ION) ? g.amb_T : target.T;
      } else {
        npt.n = g.amb_n;
        npt.T[0] = npt.T[1] = npt.T[2] = g.amb_T;
      }
    });

  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HY: return std::abs(crd[2]) > target.zh ? g.B0 : 0.;
      default: return 0.;
    }
  });

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks);
  psc.add_gauss_corrector(&marder);
  psc.add_diagnostic(&out_fields);
  psc.add_diagnostic(&out_moments);
  psc.add_diagnostic(&outp);
  psc.add_diagnostic(&oute);
  psc.add_injector(
    new InjectFromLambda<Mparticles, MfieldsState>(lf_inject_heat));

  psc.integrate();
}

int main(int argc, char** argv)
{
  // psc_init parses argv as "--option value" pairs and aborts on a bare path,
  // so hand it argc = 1 and read the params file ourselves -- the same thing
  // PSC's own psc_shock.cxx does.
  int argc_psc = 1;
  psc_init(argc_psc, argv);
  run(argc, argv);
  psc_finalize();
  return 0;
}
