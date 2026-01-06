#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_gyrokinetic_run.h>
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

double nu_ss(double qs, double ms, double ns, double us, double vts, double B0)
{
  // Like-species collision frequency.
  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double logLambda = coulombLog_sr(qs, qs, ms, ms, ns, ns, vts, vts, B0);
  return (1./sqrt(2.0))*pow(qs,4)*ns*logLambda/(3.0*(pow(2.0*M_PI,3.0/2.0))*pow(eps0,2)*pow(ms,2)*(pow(vts,3.0)));
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

  // Bulk flow speed along B-field in terms of reference temperatures.
  double upar_elc = 0.5*sqrt(mass_elc/mass_ion)*vte;

  // Coulomb logarithm and collision frequency.
  double nu_elc = nu_ss(charge_elc, mass_elc, n0, upar_elc, vte, B0);
  
  int Nz = 1; // Number of cells in z, originally 8.
  int Nvpar = 4; // Number of cells in vpar.
  int Nmu = 12; // Number of cells in mu.

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

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  // Set physical coordinates (X, Y, Z) from computational coordinates (x, y, z).
  xp[0] = zc[0]; xp[1] = zc[1]; xp[2] = zc[2];
}

void
bfield_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct lbo_relax_ctx *app = ctx;

  double B0 = app->B0;

  // Set magnetic field strength.
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = B0;
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
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

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
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
    },

    .num_diag_moments = 3,
    .diag_moments = { GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_BIMAXWELLIAN, },
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .zero_init_field = true, // Don't compute the field at t = 0.
    .is_static = true, // Don't evolve the field in time.
  };

  // Gyrokinetic app.
  struct gkyl_gk app_inp = {
    .name = "gk_lbo_relax_nonuniv_1x2v_p1",

    .cdim = ctx.cdim,
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
      .bfield_func = bfield_func,
      .bfield_ctx = &ctx
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
  };

  gkyl_gyrokinetic_run_simulation(&run_inp);

  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
