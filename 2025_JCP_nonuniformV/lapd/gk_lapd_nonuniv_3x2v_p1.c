#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_util.h>

#include <rt_arg_parse.h>

struct gk_lapd_ctx {
  int cdim, vdim; // Dimensionality.

  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double n_max; // Reference density.
  double Te_max; // Reference electron temperature.
  double Ti_max; // Reference ion temperature.
  double n_ref; // Reference density.
  double Te_ref; // Reference electron temperature.
  double Ti_ref; // Reference ion temperature.
  double B0; // Reference magnetic field.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.
  double nu_frac; // Fraction to multiply nu_sr by.

  // Source parameters
  double S0_src; // Source amplitude.
  double Te_src; // Electron source temperature.
  double Ti_src; // Ion source temperature.
  double r_src; // Radial extent of source.
  double L_src; // Length of source.
  double floor_src; // Minimum value of the source.

  double Lperp; // Perpendicular length.
  double Rmin; // Minimum radius.
  double Rmax; // Maximum radius.
  double Lz; // Box size in z.

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.
  double vpar_max_ion; // Parallel velocity extents for ions.
  double mu_max_ion; // Maximum magnetic moment for ions.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  int Nr; // Number of cells in radial direction.
  int Ntheta; // Number of cells in azimuthal direction.
  int Nz; // Number of cells in paralle direction.
  int Nvpar; // Number of cells in parallel velocity direction.
  int Nmu; // Number of cells in magnetic moment direction.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

static inline double
profileA(double r, double c_edge, double Lperp)
{
  if (r < Lperp/2.0)
    return (1.0-c_edge)*pow(1.0-pow(r/(Lperp/2.),2.0), 3.0) + c_edge;
  else
    return c_edge;
}

void
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct gk_lapd_ctx *app = ctx;
  double n_max = app->n_max;
  double Lperp = app->Lperp;

  pcg32_random_t rng = gkyl_pcg32_init(0);; // RNG for use in IC
  double perturb = 2e-3*(1.0 - 0.5*gkyl_pcg32_rand_double(&rng));

  fout[0] = n_max*profileA(r, 1.0/20.0, Lperp)*(1.0 + perturb);
}

void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct gk_lapd_ctx *app = ctx;
  double Lperp = app->Lperp;
  double Te_max = app->Te_max;

  fout[0] = Te_max*profileA(r, 1.0/5.0, Lperp);
}

void
init_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct gk_lapd_ctx *app = ctx;

  fout[0] = app->Ti_max;
}

void
init_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct gk_lapd_ctx *app = ctx;
  double Lz = app->Lz;
  double r_src = app->r_src;
  double L_src = app->L_src;

  double floor = app->floor_src;
  double S0 = app->S0_src;

  fout[0] = S0*(floor + (1.0-floor)*0.5*(1.0 - tanh((r-r_src)/L_src)));
}

void
init_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  struct gk_lapd_ctx *app = ctx;
  double Te_src = app->Te_src;
  double Lperp = app->Lperp;
  double floor = 0.01;
  double r_src = app->r_src;
  double L_src = app->L_src;

  fout[0] = Te_src*profileA(r, 1.0/2.5, Lperp);
//  fout[0] = Te_src*(floor + (1.0-floor)*0.5*(1.0 - tanh((r-r_src)/L_src)));
}

void
init_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double Ti_src = app->Ti_src;

  fout[0] = Ti_src;
}

void
init_nu_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  fout[0] = app->nu_elc;
}

void
init_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

// Computational to physical map (r,theta) -> (X,Y).
void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], theta = xc[1], z = xc[2];
  xp[0] = r*cos(theta); xp[1] = r*sin(theta); xp[2] = z;
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_lapd_ctx *gc = ctx;
  fout[0] = gc->B0;
}

void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
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

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_lapd_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

struct gk_lapd_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 3.973*GKYL_PROTON_MASS; // Ion mass.
  double mRat = 400.0;
  double me = mi/mRat;
  double qi = eV; // Ion charge.
  double qe = -eV; // Electron charge.

  double n_max = 2.0e18; // Maximum particle density.
  double Te_max = 4.0*eV; // Maximum electron temperature.
  double Ti_max = 1.0*eV; // Maximum ion temperature.

  double n_ref = n_max/2.0; // Reference density.
  double Te_ref = Te_max/2.0; // Reference electron temperature.
  double Ti_ref = Ti_max/2.0; // Reference ion temperature.

  double B0 = 0.0398; // Magnetic field magnitude.

  // Bias parameters.
  double V_bias = 0;  // Biasing voltage (V).
  double r_bias = 0.25; // Radius of the biased plate (m).

  // Simulation box size (m).
  double diam = 1.2; // Diameter. Schaffner says vessel diameter is 1 m w/ a 0.6 diameter plasma, Fisher, Shi, Frei all used 1.4 m.
  double Lz = 18.0;
  double Rmax = diam/2.0;
  double Rmin = 0.05*Rmax;
  double Lperp = diam*0.8; // Scale length of ICs and sources.

  // Physical velocity space limits
  double vt_ion_ref = sqrt(Ti_ref/mi);
  double vt_elc_ref = sqrt(Te_ref/me);

  double vpar_max_elc = 6.0*vt_elc_ref;
  double mu_max_elc = 0.5*me*pow(6.0*vt_elc_ref,2)/B0;

  double vpar_max_ion = 6.0*vt_ion_ref;
  double mu_max_ion = 0.5*mi*pow(6.0*vt_ion_ref,2)/B0;

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  // Number of cells.
  int Nr = 36, Ntheta = 36, Nz = 8;
  int Nvpar = 8, Nmu = 4;
                 
  // Source parameters.
  double c_s_max = sqrt(Te_max/mi);
  double S0_src = 1.08*n_max*c_s_max/Lz;
  double Te_src = 12.0*eV;
  double Ti_src = 1.0*eV;
  double L_src = 2*0.625/100.0;
  double r_src = 0.25;
  double floor_src = 0.01;

  // Coulomb logarithms.
  double logLambdaElc = 6.6 - 0.5*log(n_ref/1e20) + 1.5*log(Te_ref/eV);
  double logLambdaIon = 6.6 - 0.5*log(n_ref/1e20) + 1.5*log(Ti_ref/eV);

  // Collision frequencies.
  double nu_frac = 1.0/10.0;
  double nu_elc = logLambdaElc*pow(eV, 4.0)*n_ref
    /(6.0*sqrt(2.0)*pow(M_PI,1.5)*pow(eps0,2)*sqrt(me)*pow(Te_ref,1.5));

  double nu_ion = logLambdaIon*pow(eV, 4.0)*n_ref
    /(12.0*pow(M_PI,1.5)*pow(eps0,2)*sqrt(mi)*pow(Ti_ref,1.5));

  double t_end = 6.0e-3;
  double num_frames = 600;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Compute the omega_H mode frequency. The time step has to stay below it.
  int pOrder      = 1;
  double dzEff    = Lz/(Nz*(pOrder+1));
  double kParMax  = 2.0*M_PI/(2.0*dzEff);
  double Lr       = Rmax-Rmin;
  double Ly_max   = 2.0*M_PI*Rmax;
  double kperpMin = 2.0*M_PI/(2.0*GKYL_MAX2(Lr,Ly_max));
  double omega_H  = kParMax*qi*B0/(kperpMin*sqrt(me*mi));
  printf("omega_H = %g\n", omega_H);

  struct gk_lapd_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .charge_elc = qe, 
    .mass_elc = me, 
    .charge_ion = qi, 
    .mass_ion = mi,
    .n_max = n_max, 
    .Te_max = Te_max, 
    .Ti_max = Ti_max, 
    .n_ref = n_ref, 
    .Te_ref = Te_ref, 
    .Ti_ref = Ti_ref, 
    .B0 = B0, 
    .nu_elc = nu_elc, 
    .nu_ion = nu_ion, 
    .nu_frac = nu_frac,
    .S0_src = S0_src,
    .Te_src = Te_src, 
    .Ti_src = Ti_src, 
    .r_src = r_src, 
    .L_src = L_src, 
    .floor_src = floor_src,
    .Lperp = Lperp, 
    .Rmin = Rmin, 
    .Rmax = Rmax, 
    .Lz = Lz, 
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    // Number of cells.
    .Nr = Nr,  .Ntheta = Ntheta,  .Nz = Nz,
    .Nvpar = Nvpar,  .Nmu = Nmu,
    .cells = {Nr, Ntheta, Nz, Nvpar, Nmu},
    // Time stepping and I/O parameters.
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
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_app_write_field_energy(app);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_write_dt(app);
    }
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
  if (app_args.use_mpi) MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_lapd_ctx ctx = create_ctx(); // Context for init functions.

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
    .polarization_density = ctx.n_ref,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = init_density,
      .upar = init_upar,
      .temp = init_temp_elc,      
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n_ref, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te_ref, // Temperature used to calculate coulomb logarithm
      .self_nu = init_nu_elc,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
      .nuFrac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = init_density_source,
        .upar= init_upar_source,
        .temp = init_temp_elc_source,      
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
      }, 
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { "HamiltonianMoments" },
      }
    },
    
    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 4,
    .diag_moments = { "M1", "M2par", "M2perp", "BiMaxwellianMoments" },

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "HamiltonianMoments" },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { "HamiltonianMoments" },
    },

  };

  // Ions.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n_ref,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = init_density,
      .upar = init_upar,
      .temp = init_temp_ion,      
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temp = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n_ref, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti_ref, // Temperature used to calculate coulomb logarithm
      .self_nu = init_nu_ion,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
      .nuFrac = ctx.nu_frac,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = init_density_source,
        .upar= init_upar_source,
        .temp = init_temp_ion_source,      
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
      }, 
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { "HamiltonianMoments" },
      }
    },
    
    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 4,
    .diag_moments = { "M1", "M2par", "M2perp", "BiMaxwellianMoments" },

    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "HamiltonianMoments" },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { "HamiltonianMoments" },
    },

  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC}, 
                    .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC}, 
                    .lo_value = {0.0}, .up_value = {0.0}}, 
    .time_rate_diagnostics = true,
  };

  // GK app.
  struct gkyl_gk gk = {
    .name = "gk_lapd_nonuniv_3x2v_p1",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { ctx.Rmin, -M_PI, -ctx.Lz/2.0 },
    .upper = { ctx.Rmax,  M_PI,  ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
//    .cfl_frac_omegaH = 1e12,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // mapping of computational to physical space
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
      .comm = comm,
    },
  };

  // create app object
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
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, app_args.is_restart, false, -1.0);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, app_args.is_restart, false);

  // start, end and initial time-step
  double dt = t_end-t_curr;
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (step == 1 || step % 10 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);

    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step == 1 || step % 10 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

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
	calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
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
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
