#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct lbo_relax_ctx
{ 
  int cdim, vdim; // Dimensionality.

  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.

  double n0; // Reference number density (1 / m^3).
  double Te; // Electron temperature.
  double B0; // Reference magnetic field strength (Tesla).

  double upar_elc; // Parallel electron drift speed.
  double T_par_elc; // Parallel electron temperature.
  double T_perp_elc; // Perpendicular electron temperature.

  double log_lambda_elc; // Electron Coulomb logarithm.
  double nu_elc; // Electron collision frequency.

  double vte; // Electron thermal velocity.

  int Nz; // Number of cells along the magnetic field.
  int Nvpar; // Number of cells in v_parallel.
  int Nmu; // Number of cells in mu.
  int poly_order; // Polynomial order.
  int cells[GKYL_MAX_DIM]; // Number of cells in each direction.

  double Lz; // Domain size along the field.

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;

  double cfl_frac; // CFL coefficient.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

double coulombLog_sr(double qs, double qr, double ms, double mr, double ns, double nr, double vts, double vtr, double B0)
{
  // Coulomb logarithm in Gkeyll (see online documentation):
  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double hbar = GKYL_PLANCKS_CONSTANT_H/(2.*M_PI);
  double m_sr = ms*mr/(ms+mr);
  double u_sr = sqrt(3.*pow(vts,2)+3.*pow(vtr,2));
  double omega_ps = sqrt(ns*pow(qs,2)/(ms*eps0));
  double omega_pr = sqrt(nr*pow(qr,2)/(mr*eps0));
  double omega_cs = fabs(qs*B0/ms);
  double omega_cr = fabs(qr*B0/mr);
  double rMax = 1./sqrt( (pow(omega_ps,2)+pow(omega_cs,2))/(pow(vts,2)+3.*pow(vts,2))
                        +(pow(omega_pr,2)+pow(omega_cr,2))/(pow(vtr,2)+3.*pow(vts,2)) );
  double rMin = fmax(fabs(qs*qr)/(4.*M_PI*eps0*m_sr*pow(u_sr,2)), hbar/(2.*exp(0.5)*u_sr*m_sr));
  return 0.5*log(1. + pow(rMax/rMin,2));
}

struct lbo_relax_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double eV = GKYL_ELEMENTARY_CHARGE; // Proton charge.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double mass_ion = 2.014*GKYL_PROTON_MASS; // Ion mass.
  double charge_elc = -eV; // Electron charge.

  double n0 = 7.0e19; //  Reference number density (1 / m^3).
  double Te = 300.0*eV; // Electron temperature.
  double B0 = 1.0; // Magnetic field axis (simple toroidal coordinates).

  double alpha = 1.3;
  double T_par_elc = Te; // Parallel electron temperature.
  double T_perp_elc = alpha*Te; // Perpendicular electron temperature.

  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.

  // Coulomb logarithm and collision frequency.
  double log_lambda_elc = coulombLog_sr(charge_elc, charge_elc, mass_elc, mass_elc, n0, n0, vte, vte, B0);
  double nu_elc = (1/sqrt(2.0)) * log_lambda_elc * pow(fabs(charge_elc), 4.0) * n0 /
    (3.0 * pow(2.0*M_PI, 3.0 / 2.0) * pow(epsilon0, 2.0) * pow(mass_elc,2) * pow(pow(vte,2), 3.0 / 2.0));
  
  // Bulk flow speed along B-field in terms of reference temperatures.
  double upar_elc = 0.5*sqrt(mass_elc/mass_ion)*vte;

  int Nz = 1; // Number of cells in z, originally 8.
  int Nvpar = 12; // Number of cells in vpar.
  int Nmu = 16; // Number of cells in mu.

  // Physical velocity space limits.
  double vpar_max_elc =  5.0*vte;
  double mu_max_elc = mass_elc*pow(5.0*vte,2)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  double Lz = 4.0; // Domain size along the field.
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 7.0/nu_elc;
  double num_frames = 30;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct lbo_relax_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .n0 = n0,
    .Te = Te,
    .B0 = B0,
    .upar_elc = upar_elc,
    .T_par_elc = T_par_elc,
    .T_perp_elc = T_perp_elc,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .vte = vte,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .Lz = Lz,
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
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
init_density_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double n0 = app->n0;

  // Set electron total number density.
  fout[0] = n0;
}

void
init_upar_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double upar_elc = app->upar_elc;
  // Set electron parallel velocity.
  fout[0] = upar_elc;
}

void
init_Tpar_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double T_par_elc = app->T_par_elc;

  // Set electron parallel temperature.
  fout[0] = T_par_elc;
}

void
init_Tperp_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double T_perp_elc = app->T_perp_elc;

  // Set electron perpendicular temperature.
  fout[0] = T_perp_elc;
}

void
init_nu_elc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
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
  struct lbo_relax_ctx *app = ctx;

  double B0 = app->B0;

  // Set magnetic field strength.
  fout[0] = B0;
}

static inline void
mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct lbo_relax_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
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
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct lbo_relax_ctx ctx = create_ctx(); // Context for initialization functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,
    .no_collisionless_terms = true,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN,
      .density = init_density_elc,
      .ctx_density = &ctx,
      .temppar = init_Tpar_elc,
      .ctx_temppar = &ctx,
      .tempperp = init_Tperp_elc,
      .ctx_tempperp = &ctx,
      .upar = init_upar_elc,
      .ctx_upar = &ctx,
      .correct_all_moms = true,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .self_nu = init_nu_elc,
      .ctx = &ctx,
    },

    .num_diag_moments = 4,
    .diag_moments = { "M1", "M2", "M2par", "BiMaxwellianMoments" },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {
    .name = "gk_lbo_relax_nonuniv_1x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { -0.5 * ctx.Lz },
    .upper = {  0.5 * ctx.Lz },
    .cells = { cells_x[0] },

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

    .num_species = 1,
    .species = { elc },

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
