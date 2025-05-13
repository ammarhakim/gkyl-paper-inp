#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.

  // Plasma parameters
  double mi; // Ion mass.
  double me; // Electron mass.
  double qi; // Ion charge.
  double qe; // Electron charge.
  double Te0; // Electron temperature.
  double Ti0; // Ion temperature.
  double n0; // Density.
  double B_p; // Plasma magnetic field (mirror center).
  double beta; // Plasma beta in the center.
  double tau; // Temperature ratio.

  // Parameters controlling initial conditions.
  double alim;
  double alphaIC0;
  double alphaIC1;
  double nuFrac;

  double logLambdaIon; // Ion Coulomb logarithm.
  double nuIon; // Ion-ion collision freq.
  
  double vti; // Ion thermal speed.
  double vte; // Electron thermal speed.
  double c_s; // Ion sound speed.
  double omega_ci; // Ion gyrofrequency.
  double rho_s; // Ion sound gyroradius.

  double RatZeq0; // Radius of the field line at Z=0.
  double Z_min; // Minimum axial coordinate Z.
  double Z_max; // Maximum axial coordinate Z.
  double z_min; // Minimum value of the position along the field line.
  double z_max; // Maximum value of the position along the field line.
  double psi_eval; // Psi (poloidal flux) of the field line.
  double psi_in, z_in; // Auxiliary psi and z.

  // Magnetic equilibrium model parameters.
  double mcB;
  double gamma;
  double Z_m;

  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double B_bt; // Magnetic field at the banana tip.
  double R_bt; // Radius at the banana tip.
  double Z_bt; // Axial coordinate at the banana tip.
  double z_bt; // Length along the field line at the banana tip.

  // Physics parameters at mirror throat
  double R_m; // Mirror ratio.
  double B_m; // Magnetic field at the throat.
  double z_m; // Length along the field line at the banana tip.
  double n_m; // Density at the throat.
  double Ti_m; // Ion temperature
  double Ti_perp0; // Reference ion perp temperature.
  double Ti_par0; // Reference ion par temperature.
  double Ti_perp_m; // Ion perp temperature at the throat.
  double Ti_par_m; // Ion par temperature at the throat.
  double cs_m; // Ion sound speed at the throat.

  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;

  // Physical velocity space limits.
  double vpar_min_ion, vpar_max_ion;
  double mu_max_ion;
  // Computational velocity space limits.
  double vpar_lin_fac_inv, mu_lin_fac_inv; // Inverse factor of where linear mapping ends.
  double vpar_pow, mu_pow; // Power of the velocity grid.
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

  // Grid DOF.
  int num_cell_vpar;
  int num_cell_mu;
  int num_cell_z;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

double
psi_RZ(double RIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double psi = 0.5 * pow(RIn, 2.) * mcB *
               (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
                1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))));
  return psi;
}

double
R_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rout = sqrt(2. * psiIn / (mcB * 
    (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
     1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))))));
  return Rout;
}

void
Bfield_psiZ(double psiIn, double ZIn, void *ctx, double *BRad, double *BZ, double *Bmag)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rcoord = R_psiZ(psiIn, ZIn, ctx);

  BRad[0] = -(1. / 2.) * Rcoord * mcB *
          (-2. * (ZIn - Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn - Z_m) / gamma, 2.), 2.)))
           -2. * (ZIn + Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn + Z_m) / gamma, 2.), 2.))));

  BZ[0] = mcB *
        ( 1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.)))
         +1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))) );

  Bmag[0] = sqrt(pow(BRad[0], 2) + pow(BZ[0], 2));
}

double
integrand_z_psiZ(double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_in;
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, ZIn, ctx, &BRad, &BZ, &Bmag);
  return Bmag / BZ;
}

double
z_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double eps = 0.0;
  struct gkyl_qr_res integral;
  if (eps <= ZIn)
  {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, eps, ZIn, 7, 1e-14);
  }
  else
  {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, ZIn, eps, 7, 1e-14); 
    integral.res = -integral.res;
  }
  return integral.res;
}

// Invert z(Z) via root-finding.
double
root_Z_psiz(double Z, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  return app->z_in - z_psiZ(app->psi_in, Z, ctx);
}

double
Z_psiz(double psiIn, double zIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double maxL = app->Z_max - app->Z_min;
  double eps = maxL / app->num_cell_z;   // Interestingly using a smaller eps yields larger errors in some geo quantities.
  app->psi_in = psiIn;
  app->z_in = zIn;
  struct gkyl_qr_res Zout;
  if (0.0 <= zIn)
  {
    double fl = root_Z_psiz(-eps, ctx);
    double fr = root_Z_psiz(app->Z_max + eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, -eps, app->Z_max + eps, fl, fr, 1000, 1e-14);
  }
  else
  {
    double fl = root_Z_psiz(app->Z_min - eps, ctx);
    double fr = root_Z_psiz(eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, app->Z_min - eps, eps, fl, fr, 1000, 1e-14);
  }
  return Zout.res;
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;

  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx); // Cylindrical axial coordinate.

  if (fabs(Z) <= app->Z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                              exp(-pow(z - zSrc, 2) / (2.0 * pow(sigSrc, 2))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double sigSrc = app->sigSrcIon;
  double TSrc0 = app->TSrc0Ion;
  double Tfloor = app->TSrcFloorIon;

  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = (fabs(z) / z) * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m));
  }
}

void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double Ti_par0 = app->Ti_par0;
  double Ti_par_m = app->Ti_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_par_m+(Ti_par0-Ti_par_m)*tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_par_m;
  }
}

void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_mirror_ctx *app = ctx;
  double z_m = app->z_m;
  double Ti_perp0 = app->Ti_perp0;
  double Ti_perp_m = app->Ti_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_perp_m - Ti_perp0*tanh(3.*fabs(z_m-fabs(z)));
  }
  else
  {
    fout[0] = Ti_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  eval_temp_par_ion(t, xn, fout, ctx);
  double Tpar = fout[0];
  eval_temp_perp_ion(t, xn, fout, ctx);
  double Tperp = fout[0];
  fout[0] = (Tpar + 2. * Tperp) / 3.;
}

void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

// Geometry evaluation functions for the gk app
void
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  double psi = xc[0], theta = xc[1], z = xc[2];

  double Z = Z_psiz(psi, z, ctx);
  double R = R_psiZ(psi, Z, ctx);

  // Cartesian coordinates on plane perpendicular to Z axis.
  double x = R * cos(theta);
  double y = R * sin(theta);

  xp[0] = x;  xp[1] = y;  xp[2] = Z;
}

void
bmag_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  double psi = xc[0], theta = xc[1], z = xc[2];

  struct gk_mirror_ctx *app = ctx;
  double Z = Z_psiz(psi, z, ctx);
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);
  fout[0] = Bmag;
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/lin_frac_inv, then a power grid.
  double vpar_lin_fac_inv = app->vpar_lin_fac_inv;
  double vpar_pow = app->vpar_pow;
  if (fabs(cvpar) <= 1.0/vpar_lin_fac_inv)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -1.0/vpar_lin_fac_inv)
    vp[0] = -vpar_max_ion*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);
  else
    vp[0] =  vpar_max_ion*pow(vpar_lin_fac_inv,vpar_pow-1)*pow(fabs(cvpar),vpar_pow);

//  // Quadratic mu.
//  vp[1] = mu_max_ion*pow(cmu,2.0);
  // Linear map up to mu_max/lin_frac_inv, then a power grid.
  double mu_lin_fac_inv = app->mu_lin_fac_inv;
  double mu_pow = app->mu_pow;
//  if (cmu <= 1.0/mu_lin_fac_inv)
//    vp[0] = mu_max_ion*cmu;
//  else
//    vp[0] = mu_max_ion*pow(mu_lin_fac_inv,mu_pow-1)*pow(cmu,mu_pow);
  double w = 0.3;
  double f = 0.012;
  double a = mu_max_ion*(f-1.0)/(w*w-1.0);
  double b = mu_max_ion*(w*w-f)/(w*w-1.0);
  if (cmu <= w)
    vp[1] = (f*mu_max_ion/w)*cmu;
  else
    vp[1] = a*pow(cmu,2)+b;

}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0,  mu0 = GKYL_MU0;
  double mp = GKYL_PROTON_MASS,  me = GKYL_ELECTRON_MASS;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double qi =  eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;

  // Parameters controlling initial conditions.
  double alim = 0.125;
  double alphaIC0 = 2;
  double alphaIC1 = 10;

  double nuFrac = 1.0;
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
    (12. * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Perpendicular wavenumber in SI units:
  // Geometry parameters.
  double RatZeq0 = 0.10; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  // the boundary of a cell (due to AD errors).
  double Z_min = -2.5;
  double Z_max =  2.5;

  // Parameters controlling the magnetic equilibrium model.
  double mcB = 6.51292;
  double gamma = 0.124904;
  double Z_m = 0.98;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = Z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;

  // Physical velocity space limits.
  double vpar_max_ion =  16.0*vti;
  double vpar_min_ion = -vpar_max_ion;
  double mu_max_ion = mi*pow(3.0*vti,2)/(2.0*B_p);

  // Computational velocity space limits.
  double vpar_lin_fac_inv = 4;
  double vpar_pow = 3;
  double vpar_min_ion_c = -1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double vpar_max_ion_c =  1.0/pow(vpar_lin_fac_inv,(vpar_pow-1)/vpar_pow);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  double mu_lin_fac_inv = 1.0/0.012;
  double mu_pow = 2;
//  double mu_min_ion_c = 0.0;
//  double mu_max_ion_c = 1.0/pow(mu_lin_fac_inv,(mu_pow-1)/mu_pow);

  // Grid DOF:
  int num_cell_z = 192;
  int num_cell_vpar = 48; // Number of cells in the paralell velocity direction 96
  int num_cell_mu = 16;  // Number of cells in the mu direction 192
  int poly_order = 1;

  // Banana tip and throat info.
  double B_bt = 1.058278;
  double R_bt = 0.071022;
  double Z_bt = 0.467101;
  double z_bt = 0.468243;
  double R_m = 0.017845;
  double B_m = 16.662396;
  double z_m = 0.982544;

  // Physics parameters at mirror throat.
  double n_m = 1.105617e19;
  double Ti_m = 3081.437703 * eV;
  double cs_m = 4.037740e5;

  // Initial conditions parameter.s
  double Ti_perp0 = 10000 * eV;
  double Ti_par0 = 7500 * eV;
  double Ti_perp_m = 15000 * eV;
  double Ti_par_m = 1000 * eV;

  // Temporal inputs.
  double t_end = 100.0e-6;
  double num_frames = 400;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,  .vdim = vdim,
    .mi = mi,  .qi = qi,
    .me = me,  .qe = qe,
    .Te0 = Te0,  .Ti0 = Ti0,  .n0 = n0,
    .B_p = B_p,  .beta = beta,  .tau = tau,
    .alim = alim,
    .alphaIC0 = alphaIC0,
    .alphaIC1 = alphaIC1,
    .nuFrac = nuFrac,  .logLambdaIon = logLambdaIon,  .nuIon = nuIon,
    .vti = vti,  .vte = vte,  .c_s = c_s,
    .omega_ci = omega_ci,  .rho_s = rho_s,
    .RatZeq0 = RatZeq0,
    .Z_min = Z_min,  .Z_max = Z_max,
    // Parameters controlling the magnetic equilibrium model.
    .mcB = mcB,  .gamma = gamma,
    // Banana tip and throat info.
    .B_bt = B_bt,  .R_bt = R_bt,
    .Z_bt = Z_bt,  .z_bt = z_bt,
    .Z_m = Z_m,  .R_m = R_m,  .B_m = B_m,
    .z_m = z_m,  .n_m = n_m,  .Ti_m = Ti_m,
    // Initial condition parameters.
    .Ti_perp0 = Ti_perp0,  .Ti_par0 = Ti_par0,
    .Ti_perp_m = Ti_perp_m,  .Ti_par_m = Ti_par_m,  .cs_m = cs_m,
    // Source parameters
    .NSrcIon = NSrcIon,  .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,  .TSrcFloorIon = TSrcFloorIon,
    .lineLengthSrcIon = lineLengthSrcIon,  .sigSrcIon = sigSrcIon,
    // Physical velocity space limits.
    .vpar_min_ion = vpar_min_ion,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    // Computational velocity space limits.
    .vpar_lin_fac_inv = vpar_lin_fac_inv,
    .vpar_pow = vpar_pow,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_lin_fac_inv = mu_lin_fac_inv,
    .mu_pow = mu_pow,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    // Grid DOF.
    .num_cell_z = num_cell_z,
    .num_cell_vpar = num_cell_vpar,
    .num_cell_mu = num_cell_mu,
    .cells = {num_cell_z, num_cell_vpar, num_cell_mu},
    .poly_order = poly_order,
    // Temporal inputs.
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

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_mirror_ctx ctx = create_ctx(); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Populate a couple more values in the context.
  ctx.psi_eval = psi_RZ(ctx.RatZeq0, 0., &ctx);
  ctx.z_min    = z_psiZ(ctx.psi_eval, ctx.Z_min, &ctx);
  ctx.z_max    = z_psiZ(ctx.psi_eval, ctx.Z_max, &ctx);
  printf("\n");
  printf(" psi_eval = %e\n",ctx.psi_eval);
  printf(" z_min    = %e\n",ctx.z_min   );
  printf(" z_max    = %e\n",ctx.z_max   );
  printf("\n");

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,  .mass = ctx.mi,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN, 
      .density = eval_density_ion,
      .upar = eval_upar_ion,
      .temppar = eval_temp_par_ion,      
      .tempperp = eval_temp_perp_ion,   
      .ctx_density = &ctx,
      .ctx_upar = &ctx,
      .ctx_temppar = &ctx,
      .ctx_tempperp = &ctx,
    },

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = eval_density_ion_source,
        .upar = eval_upar_ion_source,
        .temp = eval_temp_ion_source,      
        .ctx_density = &ctx,
        .ctx_upar = &ctx,
        .ctx_temp = &ctx,
      }, 
    },

    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    .num_diag_moments = 4,
    .diag_moments = {"M1", "M2par", "M2perp", "BiMaxwellianMoments"},
  };

  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
  };

  struct gkyl_gk gk = {  // GK app
    .name = "gz57_1x2v_p1",
    .cdim = ctx.cdim,  .vdim = ctx.vdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {ctx.psi_eval, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .bmag_func = bmag_func, // magnetic field magnitude
      .c2p_ctx = &ctx,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {},

    .num_species = 1,
    .species = {ion},

    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
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
    if (step % 100 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);

    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step % 100 == 0)
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

  // fetch simulation statistics
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
