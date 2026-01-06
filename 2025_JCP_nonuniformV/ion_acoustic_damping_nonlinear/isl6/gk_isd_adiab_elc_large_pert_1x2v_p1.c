#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>

#include <rt_arg_parse.h>

struct gk_app_ctx {
  int cdim, vdim; // Dimensionality.

  double charge_elc; // electron charge
  double mass_elc; // electron mass
  double charge_ion; // ion charge
  double mass_ion; // ion mass
  double n0; // reference density
  double Te; // electron temperature
  double Ti; // ion temperature
  double B0; // reference magnetic field
  double nu_ion; // ion collision frequency
  double Lz; // Box size in z.
  double wavek; // Wave number of ion acoustic wave.

  double vpar_min_ion; // Minimum vparallel for ions.
  double vpar_max_ion; // Maximum vparallel for ions.
  double mu_max_ion; // Velocity space extents in mu for ions
  // Computational velocity space limits.
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  // Grid.
  int Nz; // Number of cells in z.
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

// Initial conditions.
void eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  double wavek = app->wavek;

  double alpha = 0.5;
  double perturb = alpha*cos(wavek*z);

  fout[0] = n0*(1.+perturb);
}

void eval_upar_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti = app->Ti;
  fout[0] = Ti;
}

// Collision frequency.
void eval_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void bfield_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = app->B0;
}

struct gk_app_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  double mi = 1.0; // ion mass
  double me = mi/1836.16; // electron mass
  double qi = 1.0; // ion charge
  double qe = -1.0; // electron charge

  // Reference density and temperature.
  double n0 = 1.0;
  double Te = 1.0;
  double Ti = 1.0;
  double B0 = 1.0;
  double nu_ion = 0.0;
  double wavek = 0.125;

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);

  // Simulation box size (m).
  double Lz = 2.*M_PI/wavek;

  // Physical velocity space limits.
  double vpar_max_ion =  6.0*vtIon;
  double vpar_min_ion = -vpar_max_ion;
  double mu_max_ion = mi*pow(6.0*vtIon,2)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = vpar_min_ion;
  double vpar_max_ion_c = vpar_max_ion;
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = mu_max_ion;

  int Nz = 64; // Number of cells in z, originally 8.
  int Nvpar = 64; // Number of cells in vpar.
  int Nmu = 4; // Number of cells in mu.

  double t_end = 1000.0;
  double num_frames = 1000;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_app_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .charge_elc = qe,  .mass_elc = me,
    .charge_ion = qi,  .mass_ion = mi,
    .n0 = n0,
    .Te = Te,  .Ti = Ti,
    .B0 = B0,
    .nu_ion = nu_ion,
    .Lz = Lz,
    .wavek = wavek,
    // Physical velocity space limits.
    .vpar_min_ion = vpar_min_ion,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    // Computational velocity space limits.
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    .Nz = Nz,  .Nvpar = Nvpar,  .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = eval_density_ion,
      .upar = eval_upar_ion,
      .temp = eval_temp_ion,      
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
    },

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN, },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN, },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ADIABATIC,

    .electron_mass = ctx.mass_elc,
    .electron_charge = ctx.charge_elc,
    .electron_density = ctx.n0,
    .electron_temp = ctx.Te,
  };

  // GK app
  struct gkyl_gk app_inp = {
    .name = "gk_isd_adiab_elc_large_pert_1x2v_p1",

    .cdim = ctx.cdim,
    .lower = { -ctx.Lz/2.0 },
    .upper = {  ctx.Lz/2.0 },
    .cells = { ctx.cells[0] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 1,
    .species = { ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
    .print_verbosity = {
      .enabled = true,
      .frequency = 0.01,
    },
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
