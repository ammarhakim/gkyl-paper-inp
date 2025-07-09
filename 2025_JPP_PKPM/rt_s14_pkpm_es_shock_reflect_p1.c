#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>
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

struct es_shock_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  // Physical constants (using normalized code units).
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Ion mass.
  double charge_ion; // Ion charge.

  double n0; // Reference number density.

  double Te_over_Ti; // Electron temperature / ion temperature.
  double B0; // Reference magnetic field strength.

  // Derived physical quantities (using normalized code units).
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double cs; // Sound speed.

  double Vx_drift; // Drift velocity (x-direction).

  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  // Simulation parameters.
  int Nx; // Cell count (configuration space: x-direction).
  int Nvx; // Cell count (velocity space: vx-direction).
  double Lx; // Domain size (configuration space: x-direction).
  double vx_max_elc; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion; // Domain boundary (ion velocity space: vx-direction).
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct es_shock_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using normalized code units).
  double epsilon0 = 1.0; // Permittivity of free space.
  double mu0 = 1.0; // Permeability of free space.
  double mass_elc = 1.0; // Electron mass.
  double charge_elc = -1.0; // Electron charge.
  double mass_ion = 1836.153; // Ion mass.
  double charge_ion = 1.0; // Ion charge.

  double n0 = 1.0; // Reference number density.

  double Te_over_Ti = 4.0; // Electron temperature / ion temperature.
  double B0 = 1.0; // Reference magnetic field strength.

  // Derived physical quantities (using normalized code units).
  double vte = 1.0; // Electron thermal velocity.
  double vti = vte / sqrt(Te_over_Ti * mass_ion); // Ion thermal velocity.
  double cs = vte / sqrt(mass_ion); // Sound speed.

  double Vx_drift = 3.0 * cs; // Drift velocity (x-direction).
  
  double nu_elc = 1.0e-6; // Electron collision frequency.
  double nu_ion = 1.0e-6 / sqrt(mass_ion) * (Te_over_Ti * sqrt(Te_over_Ti)); // Ion collision frequency.

  // Simulation parameters.
  int Nx = 1792; // Cell count (configuration space: x-direction).
  int Nvx = 256; // Cell count (velocity space: vx-direction).
  double Lx = 3584.0; // Domain size (configuration space: x-direction).
  double vx_max_elc = 8.0 * vte; // Domain boundary (electron velocity space: vx-direction).
  double vx_max_ion = 64.0 * vti; // Domain boundary (ion velocity space: vx-direction).
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 1500.0; // Final simulation time.
  int num_frames = 15; // Number of output frames.
  int field_energy_calcs = 15000; // Number of times to calculate field energy.
  int integrated_mom_calcs = 15000; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = 15000; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct es_shock_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .charge_ion = charge_ion,
    .n0 = n0,
    .Te_over_Ti = Te_over_Ti,
    .B0 = B0,
    .vte = vte,
    .vti = vti,
    .cs = cs,
    .Vx_drift = Vx_drift,
    .nu_elc = nu_elc,
    .nu_ion = nu_ion,
    .Nx = Nx,
    .Nvx = Nvx,
    .Lx = Lx,
    .vx_max_elc = vx_max_elc,
    .vx_max_ion = vx_max_ion,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
    .t_end = t_end,
    .num_frames = num_frames,
    .field_energy_calcs = field_energy_calcs,
    .integrated_mom_calcs = integrated_mom_calcs,
    .integrated_L2_f_calcs = integrated_L2_f_calcs,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };

  return ctx;
}

void
evalElcDistInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;
  double vx = xn[1];

  double pi = app->pi;

  double n0 = app->n0;
  double vte = app->vte;

  double F0 = (n0 / sqrt(2.0 * pi * vte * vte)) * (exp(-(vx * vx) / (2.0 * vte * vte))); // Electron distribution function (F0).
  double G = (vte * vte) * F0; // Electron distribution function (G).

  // Set electron distribution function.
  fout[0] = F0; fout[1] = G;
}

void
evalElcFluidInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;
  double x = xn[0];

  double mass_elc = app->mass_elc;
  
  double n0 = app->n0;
  double Vx_drift = app->Vx_drift;

  double mom_x = -n0 * mass_elc * Vx_drift; // Electron total momentum density (x-direction).
  double mom_y = 0.0; // Electron total momentum density (y-direction).
  double mom_z = 0.0; // Electron total momentum density (z-direction).

  // Set electron total momentum density.
  fout[0] = mom_x; fout[1] = mom_y; fout[2] = mom_z;
}

void
evalIonDistInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;
  double vx = xn[1];

  double pi = app->pi;

  double n0 = app->n0;
  double vti = app->vti;

  double F0 = (n0 / sqrt(2.0 * pi * vti * vti)) * (exp(-(vx * vx) / (2.0 * vti * vti))); // Ion distribution function (F0).
  double G = (vti * vti) * F0; // Ion distribution function (G).

  // Set electron distribution function.
  fout[0] = F0; fout[1] = G;
}

void
evalIonFluidInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;
  double x = xn[0];

  double mass_ion = app->mass_ion;

  double n0 = app->n0;
  double Vx_drift = app->Vx_drift;

  double mom_x = -n0 * mass_ion * Vx_drift; // Ion total momentum density (x-direction).
  double mom_y = 0.0; // Ion total momentum density (y-direction).
  double mom_z = 0.0; // Ion total momentum density (z-direction).

  // Set ion total momentum density.
  fout[0] = mom_x; fout[1] = mom_y; fout[2] = mom_z;
}

void
evalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double Ex = 0.0; // Total electric field (x-direction).
  double Ey = 0.0; // Total electric field (y-direction).
  double Ez = 0.0; // Total electric field (z-direction).

  double Bx = 0.0; // Total magnetic field (x-direction).
  double By = 0.0; // Total magnetic field (y-direction).
  double Bz = 0.0; // Total magnetic field (z-direction).
  
  // Set electric field.
  fout[0] = Ex; fout[1] = Ey, fout[2] = Ez;
  // Set magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
  // Set correction potentials.
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalExternalFieldInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;

  double B0 = app->B0;

  double Ex = 0.0; // External electric field (x-direction).
  double Ey = 0.0; // External electric field (y-direction).
  double Ez = 0.0; // External electric field (z-direction).

  double Bx = B0; // External magnetic field (x-direction).
  double By = 0.0; // External magnetic field (y-direction).
  double Bz = 0.0; // External magnetic field (z-direction).

  // Set external electric field.
  fout[0] = Ex; fout[1] = Ey; fout[2] = Ez;
  // Set external magnetic field.
  fout[3] = Bx; fout[4] = By; fout[5] = Bz;
}

void
evalElcNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalIonNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct es_shock_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr)) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr)) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr)) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  struct es_shock_ctx ctx = create_ctx(); // Context for initialization functions.

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
  int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NX };
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
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vx_max_elc },
    .upper = { ctx.vx_max_elc },
    .cells = { NVX },

    .init_dist = evalElcDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalElcFluidInit,
    .ctx_fluid = &ctx,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalElcNu,
      .ctx = &ctx,
    },

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  // Ions.
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -ctx.vx_max_ion },
    .upper = { ctx.vx_max_ion },
    .cells = { 8*NVX },

    .init_dist = evalIonDistInit,
    .ctx_dist = &ctx,
    .init_fluid = evalIonFluidInit,
    .ctx_fluid = &ctx,

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = evalIonNu,
      .ctx = &ctx,
    },

    .bcx = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_COPY },
  };

  // Field.
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .init = evalFieldInit,
    .ctx = &ctx,

    .ext_em = evalExternalFieldInit,
    .ext_em_ctx = &ctx,

    .bcx = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_COPY },
  };

  // PKPM app.
  struct gkyl_pkpm app_inp = {
    .name = "s14_pkpm_es_shock_reflect_p1",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0},
    .upper = { ctx.Lx },
    .cells = { NX },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .use_explicit_source = true,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

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
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&app_inp);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr);

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
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr);
    calc_integrated_mom(&im_trig, app, t_curr);
    calc_integrated_L2_f(&l2f_trig, app, t_curr);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr);
  calc_integrated_mom(&im_trig, app, t_curr);
  calc_integrated_L2_f(&l2f_trig, app, t_curr);
  write_data(&io_trig, app, t_curr, false);
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid species RHD calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
