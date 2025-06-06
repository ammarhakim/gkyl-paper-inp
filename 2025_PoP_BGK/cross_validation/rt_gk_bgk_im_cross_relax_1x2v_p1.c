#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>

struct bgk_cross_relax_ctx
{ 
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double charge_ion; // Proton charge.

  double T_par_elc; // Parallel electron temperature.
  double T_par_ion; // Parallel ion temperature.
  double alpha; // Ratio of perpendicular to parallel temperatures.

  double n0; // Reference number density (1 / m^3).
  double B0; // Reference magnetic field strength (Tesla).

  // Derived physical quantities (using non-normalized physical units).
  double T_perp_elc; // Perpendicular electron temperature.
  double T_perp_ion; // Perpendicular ion temperature.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double vte_par; // Parallel electron thermal velocity.
  double vti_par; // Parallel ion thermal velocity.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.

  double upar_elc; // Parallel electron velocity.
  double upar_ion; // Parallel ion velocity.

  // Simulation parameters.
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct bgk_cross_relax_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double T_par_elc = 300.0 * GKYL_ELEMENTARY_CHARGE; // Parallel electron temperature.
  double T_par_ion = 200.0 * GKYL_ELEMENTARY_CHARGE; // Parallel ion temperature.
  double alpha = 1.3; // Ration of perpendicular to parallel temperatures.

  double n0 = 7.0e19; //  Reference number density (1 / m^3).
  double B0 = 1.0; // Reference magnetic field strength (Tesla).

  // Derived physical quantities (using non-normalized physical units).
  double T_perp_elc = alpha * T_par_elc; // Perpendicular electron temperature.
  double T_perp_ion = alpha * T_par_ion; // Perpendicular ion temperature.

  double Te = (T_par_elc + (2.0 * T_perp_elc)) / 3.0; // Electron temperature.
  double Ti = (T_par_ion + (2.0 * T_perp_ion)) / 3.0; // Ion temperature.

  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion); // Electron Coulomb logarithm.
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion); // Ion Coulomb logarithm.
  double nu_elc = log_lambda_elc * pow(charge_ion, 4.0) * n0 /
    (6.0 * sqrt(2.0) * pow(pi, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_elc) * pow(Te, 3.0 / 2.0)); // Electron collision frequency.
  double nu_ion = log_lambda_ion * pow(charge_ion, 4.0) * n0 /
    (12.0 * pow(pi, 3.0 / 2.0) * pow(epsilon0, 2.0) * sqrt(mass_ion) * pow(Ti, 3.0 / 2.0)); // Ion collision frequency.
  
  double vte_par = sqrt(T_par_elc / mass_elc); // Parallel electron thermal velocity.
  double vti_par = sqrt(T_par_ion / mass_ion); // Parallel ion thermal velocity.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.

  double upar_elc = 0.5 * sqrt(mass_elc / mass_ion) * vte; // Parallel electron velocity.
  double upar_ion = 50.0 * (mass_elc / mass_ion) * vti; // Parallel ion velocity.

  // Simulation parameters.
  int Nz = 1; // Cell count (configuration space: z-direction).
  int Nvpar = 16; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 16; // Cell count (velocity space: magnetic moment direction).
  double Lz = 4.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 5.0 * vte_par; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc * pow(5.0 * vte_par, 2.0) / (2.0 * B0); // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 5.0 * vti_par; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = mass_ion * pow(5.0 * vti_par, 2.0) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 500.0 / nu_ion; // Final simulation time.
  int num_frames = 50; // Number of output frames.
  int field_energy_calcs = INT_MAX; // Number of times to calculate field energy.
  int integrated_mom_calcs = INT_MAX; // Number of times to calculate integrated moments.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct bgk_cross_relax_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .T_par_elc = T_par_elc,
    .T_par_ion = T_par_ion,
    .n0 = n0,
    .B0 = B0,
    .T_perp_elc = T_perp_elc,
    .T_perp_ion = T_perp_ion,
    .Te = Te,
    .Ti = Ti,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .vte_par = vte_par,
    .vti_par = vti_par,
    .vte = vte,
    .vti = vti,
    .upar_elc = upar_elc,
    .upar_ion = upar_ion,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .Lz = Lz,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double n0 = app->n0;

  // Set electron total number density.
  fout[0] = n0;
}

void
evalElcTparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_par_elc = app->T_par_elc;

  // Set electron parallel temperature.
  fout[0] = T_par_elc;
}

void
evalElcTperpInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_perp_elc = app->T_perp_elc;

  // Set electron perpendicular temperature.
  fout[0] = T_perp_elc;
}

void
evalElcUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double upar_elc = app->upar_elc;

  // Set electron parallel velocity.
  fout[0] = upar_elc;
}

void
evalIonDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double n0 = app->n0;

  // Set ion total number density.
  fout[0] = n0;
}

void
evalIonTparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_par_ion = app->T_par_ion;

  // Set ion parallel temperature.
  fout[0] = T_par_ion;
}

void
evalIonTperpInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double T_perp_ion = app->T_perp_ion;

  // Set ion perpendicular temperature.
  fout[0] = T_perp_ion;
}

void
evalIonUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;
  
  double upar_ion = app->upar_ion;

  // Set ion parallel velocity.
  fout[0] = upar_ion;
}

void
evalElcNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalIonNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bmag_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct bgk_cross_relax_ctx *app = ctx;

  double B0 = app->B0;

  // Set magnetic field strength.
  fout[0] = B0;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_gyrokinetic_app_write(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct bgk_cross_relax_ctx ctx = create_ctx(); // Context for initialization functions.

  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nz);
  int NVPAR = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nmu);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif

  int ccells[] = { NZ };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif

  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vpar_max_elc, 0.0 },
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc },
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,
    .no_by = true,

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = evalElcDensityInit,
      .ctx_density = &ctx,
      .temppar = evalElcTparInit,
      .ctx_temppar = &ctx,
      .tempperp = evalElcTperpInit,
      .ctx_tempperp = &ctx,
      .upar = evalElcUparInit,
      .ctx_upar = &ctx,
      .correct_all_moms = true,
    },
    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    },   
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Te,
      .self_nu = evalElcNu,
      .ctx = &ctx,
      .has_implicit_coll_scheme = true,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .num_diag_moments = 6,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "BiMaxwellianMoments" },
  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vpar_max_ion, 0.0 },
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion },
    .cells = { NVPAR, NMU },
    .polarization_density = ctx.n0,
    .no_by = true,

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = evalIonDensityInit,
      .ctx_density = &ctx,
      .temppar = evalIonTparInit,
      .ctx_temppar = &ctx,
      .tempperp = evalIonTperpInit,
      .ctx_tempperp = &ctx,
      .upar = evalIonUparInit,
      .ctx_upar = &ctx,
      .correct_all_moms = true,
    },
    .correct = {
      .correct_all_moms = true, 
      .use_last_converged = true, 
      .iter_eps = 1e-12,
      .max_iter = 10,
    },      
    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Ti,
      .self_nu = evalIonNu,
      .ctx = &ctx,
      .has_implicit_coll_scheme = true,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    
    .num_diag_moments = 6,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "BiMaxwellianMoments" },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {
    .name = "gk_bgk_im_cross_relax_1x2v_p1",

    .cdim = 1, .vdim = 2,
    .lower = { -0.5 * ctx.Lz },
    .upper = { 0.5 * ctx.Lz },
    .cells = { NZ },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = { 0.0, 0.0 },

      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr, false);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr, false);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr, false);
    calc_integrated_mom(&im_trig, app, t_curr, false);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

        calc_field_energy(&fe_trig, app, t_curr, true);
        calc_integrated_mom(&im_trig, app, t_curr, true);
        write_data(&io_trig, app, t_curr, true);

        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }
  
  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);
  gkyl_gyrokinetic_app_stat_write(app);
  
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
