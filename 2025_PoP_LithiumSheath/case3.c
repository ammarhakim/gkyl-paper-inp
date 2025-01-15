#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_bc_emission.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

struct sheath_ctx {
  int cdim;
  int vdim;
  double epsilon0;
  double mu0;
  double q0;
  double charge_elc; // electron charge
  double mass_elc; // electron mass
  double charge_ion; // ion charge
  double mass_ion; // ion mass
  double n0;
  double Te; // electron to ion temperature ratio
  double Ti;
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double lambda_D;
  double Lx; // size of the box
  double Ls;
  double omega_pe;
  double phi;
  double deltahat_ts;
  double Ehat_ts;
  double t1;
  double t2;
  double t3;
  double t4;
  double s;
  double E_f;
  double phi_r;
  double nu_ee;
  double nu_ii;
  int Nx;
  int Nv;
  int cells[GKYL_MAX_DIM];
  int num_emission_species;
  double t_end;
  int num_frames;
  double dt_failure_tol;
  int num_failures_max;
};

static inline double sq(double x) { return x*x; }

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte;
  double n = app->n0;
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncElcSource(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte;
  double Ls = app->Ls;
  double fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  if(fabs(x) < Ls) {
    fout[0] = 2*(Ls - fabs(x))/Ls*fv;
  } else {
    fout[0] = 0.0;
  }
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti;
  double n = app->n0;
  double fv = n/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncIonSource(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti;
  double Ls = app->Ls;
  double fv = 1.0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v)/(2*sq(vt))));
  if(fabs(x) < Ls) {
    fout[0] = 2*(Ls - fabs(x))/Ls*fv;
  } else {
    fout[0] = 0.0;
  }
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalElcNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double nu = app->nu_ee;
  double lambda_D = app->lambda_D;

  // Set collision frequency.
  fout[0] = nu/(1 + exp(fabs(x)/(6.0*lambda_D) - 8.0/1.5));
}

void
evalIonNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double nu = app->nu_ii;
  double lambda_D = app->lambda_D;

  // Set collision frequency.   
  fout[0] = nu/(1 + exp(fabs(x)/(6.0*lambda_D) - 8.0/1.5));
}

struct sheath_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 1; // Dimensionality
  double mass_elc = 9.109e-31;
  double mass_ion = 1836.153*mass_elc;
  double q0 = 1.602e-19;
  double epsilon0 = 8.854e-12;

  double n0 = 1.0e16;
  double Te = 150.0*q0;
  double Ti = 150.0*q0;
  double vte = sqrt(Te/mass_elc);
  double vti = sqrt(Ti/mass_ion);
  double c = 6.0*vte;
  double mu0 = 1.0/(epsilon0*c*c);

  double lambda_D = sqrt(epsilon0*Te/(n0*q0*q0));
  double Lx = 128.0*lambda_D;
  double Ls = 40.0*lambda_D;
  double omega_pe = sqrt(n0*q0*q0/(epsilon0*mass_elc));

  double nu_ee = vte/(50.0*lambda_D);
  double nu_ii = vti/(50.0*lambda_D);

  // SEE parameters
  double phi = 3.0;
  double deltahat_ts = 4.208;
  double Ehat_ts = 354.52;
  double t1 = 0.66;
  double t2 = 0.8;
  double t3 = 0.7;
  double t4 = 1.0;
  double s = 1.79;
  double E_f = 290.31;
  double phi_r = 144.49;

  int Nx = 960;
  int Nv = 512;
  struct sheath_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .q0 = q0,
    .charge_elc = -q0,
    .mass_elc = mass_elc,
    .charge_ion = q0,
    .mass_ion = mass_ion,
    .n0 = n0,
    .Te = Te,
    .Ti = Ti,
    .vte = vte,
    .vti = vti,
    .lambda_D = lambda_D,
    .Lx = Lx,
    .Ls = Ls,
    .omega_pe = omega_pe,
    .phi = phi,
    .deltahat_ts = deltahat_ts,
    .Ehat_ts = Ehat_ts,
    .t1 = t1,
    .t2 = t2,
    .t3 = t3,
    .t4 = t4,
    .s = s,
    .E_f = E_f,
    .phi_r = phi_r,
    .nu_ee = nu_ee,
    .nu_ii = nu_ii,
    .Nx = Nx,
    .Nv = Nv,
    .cells = {Nx, Nv},
    .num_emission_species = 1,
    .t_end = 10000.0/omega_pe,
    .num_frames = 100,
    .dt_failure_tol = 1.0e-4,
    .num_failures_max = 20,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app,
  double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}


void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr,
  bool force_write)
{
  gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, iot->curr - 1);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, iot->curr - 1);
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

  struct sheath_ctx ctx = create_ctx(); // Context for initialization functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_vlasov_comms_new(app_args.use_mpi,
    app_args.use_gpu, stderr);

  char in_species[1][128] = { "elc" };
  /* struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_secondary_electron_copper_new(ctx.num_emission_species, 0.0, in_species, app_args.use_gpu); */
  struct gkyl_emission_spectrum_model *spectrum_model[1];
  spectrum_model[0] = gkyl_emission_spectrum_chung_everhart_new(ctx.q0, ctx.phi,
    app_args.use_gpu);
  struct gkyl_emission_yield_model *yield_model[1];
  yield_model[0] = gkyl_emission_yield_furman_pivi_new(ctx.q0, ctx.deltahat_ts,
    ctx.Ehat_ts, ctx.t1, ctx.t2, ctx.t3, ctx.t4, ctx.s, app_args.use_gpu);
  struct gkyl_emission_elastic_model *elastic_model = gkyl_emission_elastic_cazaux_new(ctx.q0, ctx.E_f, ctx.phi_r, app_args.use_gpu);
  struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_new(ctx.num_emission_species,
    1000.0/ctx.omega_pe, true, spectrum_model, yield_model, elastic_model, in_species);

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -4.0*ctx.vte},
    .upper = { 4.0*ctx.vte}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncElc,
      .ctx_func = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.Ls,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFuncElcSource,
        .ctx_func = &ctx,
      },
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalElcNu,
      .ctx = &ctx,
      .fixed_temp_relax = true,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_EMISSION,
                 .aux_ctx = bc_ctx, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -6.0*ctx.vti},
    .upper = { 6.0*ctx.vti}, 
    .cells = { cells_v[0] },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncIon,
      .ctx_func = &ctx,
    },

    .source = {
      .source_id = GKYL_BFLUX_SOURCE,
      .source_length = ctx.Ls,
      .source_species = "ion",
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_FUNC,
        .func = evalDistFuncIonSource,
        .ctx_func = &ctx,
      },
    },

    .collisions =  {
      .collision_id = GKYL_BGK_COLLISIONS,
      .self_nu = evalIonNu,
      .ctx = &ctx,
      .fixed_temp_relax = true,
    },

    .bcx = {
      .lower = { .type = GKYL_SPECIES_REFLECT, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_NEUMANN },
      .up_type = { GKYL_POISSON_DIRICHLET },
      .lo_value = { 0.0 }, .up_value = { 0.0 }
    },
  };

  // VM app
  struct gkyl_vm app_inp = {
    .name = "sheath",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0 },
    .upper = { ctx.Lx },
    .cells = { cells_x[0] },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,
    .is_electrostatic = true,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

  // Initial & final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_vlasov_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  write_data(&trig_write, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&trig_write, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  write_data(&trig_write, app, t_curr, false);
  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_stat_write(app);

  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  gkyl_vlasov_app_cout(app, stdout, "\n");
  gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.nio);
  gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_vlasov_app_release(app);
  for (int i=0; i<ctx.num_emission_species; ++i) {
    gkyl_emission_spectrum_model_release(spectrum_model[i]);
    gkyl_emission_yield_model_release(yield_model[i]);
  }
  gkyl_emission_elastic_model_release(elastic_model);
  gkyl_bc_emission_release(bc_ctx);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
