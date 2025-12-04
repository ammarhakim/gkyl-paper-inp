#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>

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
  double vpar_lin_fac_inv; // Inverse factor of where linear mapping ends.
  double vpar_pow;
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

  double alpha = 0.01;
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

void bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->B0;
}

// Velocity space mappings.
void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double vpar_min = app->vpar_min_ion;
  double vpar_max = app->vpar_max_ion;
  double mu_max = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
//  // Uniform vpar.
//  vp[0] = cvpar;
//  // Quadratic vpar.
//  if (cvpar < 0.)
//    vp[0] = vpar_min*pow(cvpar,2);
//  else
//    vp[0] = vpar_max*pow(cvpar,2);
//  // Linear map up to vpar_max/lin_frac_inv, then a power grid.
  double vpar_lin_fac_inv = app->vpar_lin_fac_inv;
  double vpar_pow = app->vpar_pow;
  if (fabs(cvpar) <= 1.0/vpar_lin_fac_inv)
    vp[0] = vpar_max*cvpar;
  else if (cvpar < -1.0/vpar_lin_fac_inv)
    vp[0] = -vpar_max*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);
  else
    vp[0] =  vpar_max*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);
//  // Quadratic vpar centered at +/- 1.
//  double inflex_pt = (2.0-sqrt(2.0))/2.0;
//  if (cvpar < -inflex_pt)
//    vp[0] = -1.0-8.0*pow(cvpar+inflex_pt,2);
//  else if (cvpar < 0.0)
//    vp[0] = -1.0+pow(1.0/inflex_pt,2)*pow(cvpar+inflex_pt,2);
//  else if (cvpar < inflex_pt)
//    vp[0] =  1.0-pow(1.0/inflex_pt,2)*pow(cvpar-inflex_pt,2);
//  else
//    vp[0] = 1.0+8.0*pow(cvpar-inflex_pt,2);
//  if (cvpar < -inflex_pt)
//    vp[0] = -2.0-6.0*pow(cvpar+inflex_pt,2);
//  else if (cvpar < 0.0)
//    vp[0] = -2.0+2.0*pow(1.0/inflex_pt,2)*pow(cvpar+inflex_pt,2);
//  else if (cvpar < inflex_pt)
//    vp[0] =  2.0-2.0*pow(1.0/inflex_pt,2)*pow(cvpar-inflex_pt,2);
//  else
//    vp[0] =  2.0+6.0*pow(cvpar-inflex_pt,2);

  vp[1] = mu_max*cmu;
//  vp[1] = mu_max*pow(cmu,2);
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
  double vpar_min_ion = -5.0*vtIon;
  double vpar_max_ion =  5.0*vtIon;
  double mu_max_ion = mi*pow(5.0*vtIon,2)/(2.0*B0);

  // Computational velocity space limits.
//  double vpar_min_ion_c = -1.0;
//  double vpar_max_ion_c =  1.0;
  double vpar_lin_fac_inv = 2.0;
  double vpar_pow = 2;
  double vpar_min_ion_c = -1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double vpar_max_ion_c =  1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;

  int Nz = 32; // Number of cells in z, originally 8.
  int Nvpar = 16; // Number of cells in vpar.
  int Nmu = 4; // Number of cells in mu.

  double t_end = 50.0;
  double num_frames = 100;
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
    .vpar_lin_fac_inv = vpar_lin_fac_inv,
    .vpar_pow = vpar_pow,
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

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
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
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = eval_density_ion,
      .upar = eval_upar_ion,
      .temp = eval_temp_ion,      
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
    },

//    .collisions =  {
//      .collision_id = GKYL_LBO_COLLISIONS,
//      .ctx = &ctx,
//      .self_nu = eval_nu_ion,
//    },
    
    .num_diag_moments = 4,
    .diag_moments = { "M1", "M2par", "M2perp", "BiMaxwellianMoments" },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "HamiltonianMoments" },
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
  struct gkyl_gk gk = {
    .name = "gk_ion_sound_adiabatic_elc_1x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { -ctx.Lz/2.0 },
    .upper = {  ctx.Lz/2.0 },
    .cells = { ctx.cells[0] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .bmag_func = bmag_func, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_ctx = &ctx
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

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
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

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false);

  // start, end and initial time-step
  double dt = t_end-t_curr;
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, t_curr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);
  
  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
