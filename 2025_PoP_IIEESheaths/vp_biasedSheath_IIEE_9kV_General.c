#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_vlasov.h>
#include <gkyl_bc_emission.h>
#include <rt_arg_parse.h>

struct sheath_ctx {
  int cdim;
  int vdim;
  double epsilon0;
  double q0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double n0;
  double Te; // electron to ion temperature ratio
  double Ti;
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double lambda_D;
  double mfp;
  double nu_elc;
  double nu_ion;
  double Lx; // size of the box
  double Ls;
  double omega_pe;
  double gauss_E0;
  double gauss_tau;
  double intWall;  
  double A2;
  double A3; 
  double A4;
  double A5; 
  double nw;   // Wall number density
  double phi_bias;  // bias potential
  int Nx;
  int Nvi;
  int Nve;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
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
    fout[0] = (Ls - fabs(x))/Ls*fv;
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
    fout[0] = (Ls - fabs(x))/Ls*fv;
  } else {
    fout[0] = 0.0;
  }
}

// Functions setting the collision frequency profile.
double
nu_fk(double x, double lambda_D)
{
  return 1.0/(1.0 + exp(x/(12.0*lambda_D) - 16.0/3.0));
}

double
nu_profile(double x, double lambda_D, double Lx)
{
  return 0.5*(nu_fk(x-0.5*Lx+128*lambda_D, lambda_D)+nu_fk(-x-0.5*Lx+128*lambda_D, lambda_D));
}

// Electron and ion collision frequency profiles.
void
eval_nu_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];

  struct sheath_ctx *app = ctx;

  double nu_elc = app->nu_elc;
  double lambda_D = app->lambda_D;
  double Lx = app->Lx;

  fout[0] = nu_elc * nu_profile(x, lambda_D, Lx);
}

void
eval_nu_ion(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  double x = xn[0];

  struct sheath_ctx *app = ctx;

  double nu_ion = app->nu_ion;
  double lambda_D = app->lambda_D;
  double Lx = app->Lx;

  fout[0] = nu_ion * nu_profile(x, lambda_D, Lx);
}

struct sheath_ctx
create_ctx(void)
{
  struct sheath_ctx ctx = {
    .cdim = 1,
    .vdim = 1,
    .epsilon0 = GKYL_EPSILON0,
    .chargeElc = -GKYL_ELEMENTARY_CHARGE,
    .massElc = GKYL_ELECTRON_MASS,
    .chargeIon = GKYL_ELEMENTARY_CHARGE,
    .massIon = GKYL_PROTON_MASS,
    .n0 = 1.1e23,
    .Te = 2000.0*GKYL_ELEMENTARY_CHARGE,
    .Ti = 2000.0*GKYL_ELEMENTARY_CHARGE,
    .phi_bias = 9.0e3,   // V
    .vte = sqrt(ctx.Te/ctx.massElc),
    .vti = sqrt(ctx.Ti/ctx.massIon),
    .lambda_D = sqrt(ctx.epsilon0*ctx.Te/(ctx.n0*GKYL_ELEMENTARY_CHARGE*GKYL_ELEMENTARY_CHARGE)),
    .mfp = 50.0*ctx.lambda_D,
    .nu_elc = ctx.vte/ctx.mfp,
    .nu_ion = ctx.vti/ctx.mfp,
    .Lx = 256.0*ctx.lambda_D,
    .Ls = 100.0*ctx.lambda_D,
    .omega_pe = sqrt(ctx.n0*GKYL_ELEMENTARY_CHARGE*GKYL_ELEMENTARY_CHARGE/(ctx.epsilon0*ctx.massElc)),
    .gauss_E0 = 2.807,      // eV
    .gauss_tau = 1.210,
    .intWall = 63366509.0954913, // m/J
    .A2 = 4.194,
    .A3 = 4.649e3,
    .A4 = 8.113e1,
    .A5 = 2.242e-2,
    .nw = 8.491e+28,
    .Nx = 1024,
    .Nvi = 64,
    .Nve = 512,
    .num_emission_species = 1,
    .t_end = 10000.0/ctx.omega_pe,
    .num_frames = 100,
    .dt_failure_tol = 1.0e-6,
    .num_failures_max = 20,
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger* iot, gkyl_vlasov_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_vlasov_app_write(app, t_curr, iot->curr - 1);

    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, t_curr, iot->curr - 1);

    gkyl_vlasov_app_write_integrated_mom(app);
    gkyl_vlasov_app_write_field_energy(app);
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
  struct gkyl_comm *comm = gkyl_vlasov_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  char in_species[1][128] = { "ion" };
  struct gkyl_emission_spectrum_model *spectrum_model[1];
  spectrum_model[0] = gkyl_emission_spectrum_gaussian_new(ctx.chargeIon, ctx.gauss_E0, ctx.gauss_tau, app_args.use_gpu);
  struct gkyl_emission_yield_model *yield_model[1];
  yield_model[0] = gkyl_emission_yield_schou_new(ctx.chargeIon, ctx.intWall, ctx.A2, ctx.A3, ctx.A4, ctx.A5, ctx.nw, app_args.use_gpu);
  struct gkyl_bc_emission_ctx *bc_ctx = gkyl_bc_emission_new(ctx.num_emission_species, 1000.0/ctx.omega_pe, false, spectrum_model, yield_model, NULL, in_species);

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -4.0*ctx.vte},
    .upper = { 4.0*ctx.vte}, 
    .cells = { ctx.Nve },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncElc,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = eval_nu_elc,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
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

    .bcx = {
      .lower = { .type = GKYL_SPECIES_EMISSION,
                 .aux_ctx = bc_ctx, },
      .upper = { .type = GKYL_SPECIES_EMISSION,
                 .aux_ctx = bc_ctx, },
    },
    
    .num_diag_moments = 4,
    .diag_moments = { "M0", "M1i", "M2", "M3i", "intM0" },
  };

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0*ctx.vti},
    .upper = { 6.0*ctx.vti}, 
    .cells = { ctx.Nvi },

    .num_init = 1,
    .projection[0] = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalDistFuncIon,
      .ctx_func = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .self_nu = eval_nu_ion,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
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

    .bcx = {
      .lower = { .type = GKYL_SPECIES_ABSORB, },
      .upper = { .type = GKYL_SPECIES_ABSORB, },
    },
    
    .num_diag_moments = 4,
    .diag_moments = { "M0", "M1i", "M2", "M3i", "intM0" },
  };

  // Field.
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0,
    .poisson_bcs = {
      .lo_type = { GKYL_POISSON_DIRICHLET },
      .up_type = { GKYL_POISSON_DIRICHLET },
      .lo_value = {ctx.phi_bias }, .up_value = { 0.0 }
    },
  };

  // VM app
  struct gkyl_vm app_inp = {
    .name = "vp_biasedSheath_IIEE_9kV_1x1v_p2",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { ctx.Nx },
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
    }

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
  gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  gkyl_vlasov_app_calc_field_energy(app, t_curr);
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

    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

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
	write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

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
  gkyl_vlasov_app_release(app);
  gkyl_vlasov_comms_release(comm);
  for (int i=0; i<ctx.num_emission_species; ++i) {
    gkyl_emission_spectrum_model_release(spectrum_model[i]);
    gkyl_emission_yield_model_release(yield_model[i]);
  }
  gkyl_bc_emission_release(bc_ctx);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
