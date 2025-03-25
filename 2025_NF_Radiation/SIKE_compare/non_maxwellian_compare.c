#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
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
#define BUFFER_LEN 100

struct sheath_ctx
{
  // Mathematical constants (dimensionless).
  double pi;

  double *f_from_file;
  double *v_from_file;
  double *te_from_file;
  size_t num_te;
  size_t num_v;
  int te_ind;
  // Physical constants (using non-normalized physical units).
  double epsilon0; // Permittivity of free space.
  double mass_elc; // Electron mass.
  double charge_elc; // Electron charge.
  double mass_ion; // Proton mass.
  double massAr; // Argon mass
  double charge_ion; // Proton charge.

  double Te; // Electron temperature.
  double Ti; // Ion temperature.
  double TAr; // Argon temperature.
  double n0; // Reference number density (1 / m^3).
  double n0high; // Reference number density (1 / m^3).
  double n0low; // Reference number density (1 / m^3).

  double nu_frac; // Collision frequency fraction.
  double k_perp_rho_s; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Derived physical quantities (using non-normalized physical units).
  double R_outer; // Radial coordinate (simple toroidal coordinates).
  double B0; // Reference magnetic field strength (Tesla).

  double log_lambda_elc; // Electron Coulomb logarithm.
  double log_lambda_ion; // Ion Coulomb logarithm.
  double nu_elc; // Electron collision frequency.
  double nu_ion; // Ion collision frequency.

  double c_s; // Sound speed.
  double vte; // Electron thermal velocity.
  double vti; // Ion thermal velocity.
  double vtAr; // Argon thermal velocity.
  double omega_ci; // Ion cyclotron frequency.
  double rho_s; // Ion-sound gyroradius.
  double k_perp; // Perpendicular wavenumber (for Poisson solver).

  double n_src; // Source number density.
  double T_src; // Source temperature.
  double cx; // Source mean position (x-direction).
  double cz; // Source standard deviation (x-direction).
  double x_center; // Source center (x-direction).

  // Velocity mapping parameters
  double vscale; // How much smaller cells below vcut and mucut are compared to above
  double vcut; // magnitude of parralel velocity to switch from small to large cells
  double mucut; // mu velocity to switch from small to large cells
  
  // Simulation parameters.
  int cdim;
  int vdim;
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM];
  double Lx; // Domain size (configuration space: x-direction).
  double Lz; // Domain size (configuration space: z-direction).
  double vpar_max_elc; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion; // Domain boundary (ion velocity space: magnetic moment direction).
  double vpar_max_Ar; // Velocity space extents in vparallel for Ar
  double mu_max_Ar; // Velocity space extents in mu for Ar

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct sheath_ctx
create_ctx(void)
{
  // Mathematical constants (dimensionless).
  double pi = M_PI;

  // Physical constants (using non-normalized physical units).
  double epsilon0 = GKYL_EPSILON0; // Permittivity of free space.
  double mass_elc = GKYL_ELECTRON_MASS; // Electron mass.
  double charge_elc = -GKYL_ELEMENTARY_CHARGE; // Electron charge.
  double mass_ion = 2.014 * GKYL_PROTON_MASS; // Proton mass.
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double charge_ion = GKYL_ELEMENTARY_CHARGE; // Proton charge.

  double Te = 70.0 * GKYL_ELEMENTARY_CHARGE; // Electron temperature.
  double Ti = 120.0 * GKYL_ELEMENTARY_CHARGE; // Ion temperature.
  double n0 = 3.0e19; //  Reference number density (1 / m^3).
  double TAr = 10.0*GKYL_ELEMENTARY_CHARGE;
  double n0high = n0*250.0; // Particle density in 1/m^3
  double n0low = n0 * 0.000001/5; // Particle density in 1/m^3

  double B0 = 2.51; // Magnetic field axis (simple toroidal coordinates).
  double R_outer = 5.6;
  double cx = 0.00159/2.16;
  double cz = 7.22285;
  double x_center = 0.03;

  double nu_frac = 1.0; // Collision frequency fraction.
  double k_perp_rho_s = 0.3; // Product of perpendicular wavenumber and ion-sound gyroradius.

  // Coulomb logarithms.
  double log_lambda_elc = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Te / charge_ion);
  double log_lambda_ion = 6.6 - 0.5 * log(n0 / 1.0e20) + 1.5 * log(Ti / charge_ion);

  // Collision frequencies.
  double nu_elc = nu_frac * log_lambda_elc * pow(charge_ion,4) * n0 /
    (6.0 * sqrt(2.0) * pow(pi,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_elc) * pow(Te,3.0/2.0));
  double nu_ion = nu_frac * log_lambda_ion * pow(charge_ion,4) * n0 /
    (12.0 * pow(pi,3.0/2.0) * pow(epsilon0,2) * sqrt(mass_ion) * pow(Ti,3.0/2.0));
  
  double c_s = sqrt(Te / mass_ion); // Sound speed.
  double vte = sqrt(Te / mass_elc); // Electron thermal velocity.
  double vti = sqrt(Ti / mass_ion); // Ion thermal velocity.
  double vtAr = sqrt(TAr/mAr); // Argon thermal velocity
  double omega_ci = fabs(charge_ion * B0 / mass_ion); // Ion cyclotron frequency.
  double rho_s = c_s / omega_ci; // Ion-sound gyroradius.
  double k_perp = k_perp_rho_s / rho_s; // Perpendicular wavenumber (for Poisson solver).

  double n_src = 2.1e23; // Source number density.
  //double T_src = 200*5.0/2.0*0.3*1.6251586572438161*1.17*5.0/8.0 * GKYL_ELEMENTARY_CHARGE; // Source Temperature
  double T_src = 1805.0 * GKYL_ELEMENTARY_CHARGE; // Source Temperature

  // Velocity mapping parameters
  double vscale = 0.2;
  double vcut = 5e6;
  double mucut = 5e-17;
  
  // Simulation parameters.
  int cdim = 1;
  int vdim = 2;
  int Nz = 2; // Cell count (configuration space: z-direction).
  int Nvpar = 256; // Cell count (velocity space: parallel velocity direction).
  int Nmu = 128; // Cell count (velocity space: magnetic moment direction).
  double Lx = 0.06; // Domain size (configuration space: x-direction).
  double Lz = 2.0; // Domain size (configuration space: z-direction).
  double vpar_max_elc = 2.29628446e7;//6.0 * vte + vcut/vscale; // Domain boundary (electron velocity space: parallel velocity direction).
  double mu_max_elc = mass_elc*vpar_max_elc*vpar_max_elc/(2*B0);//3*(3.0 / 2.0) * 0.5 * mass_elc * pow(4.0 * vte,2) / (2.0 * B0);// + mucut/vscale; // Domain boundary (electron velocity space: magnetic moment direction).
  double vpar_max_ion = 6.0 * vti; // Domain boundary (ion velocity space: parallel velocity direction).
  double mu_max_ion = 3*(3.0 / 2.0) * 0.5 * mass_ion * pow(4.0 * vti,2) / (2.0 * B0); // Domain boundary (ion velocity space: magnetic moment direction).
  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 18.*mAr*vtAr*vtAr/(2.0*B0);


  double t_end = 2.0e-11; // Final simulation time.
  int num_frames = 1; // Number of output frames.
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct sheath_ctx ctx = {
    .pi = pi,
    .epsilon0 = epsilon0,
    .mass_elc = mass_elc,
    .charge_elc = charge_elc,
    .mass_ion = mass_ion,
    .massAr = mAr,
    .charge_ion = charge_ion,
    .Te = Te,
    .Ti = Ti,
    .TAr = TAr, 
    .n0 = n0,
    .n0high = n0high,
    .n0low = n0low,
    .nu_frac = nu_frac,
    .k_perp_rho_s = k_perp_rho_s,
    .B0 = B0,
    .R_outer = R_outer,
    .x_center = x_center,
    .cx = cx,
    .cz = cz,
    .log_lambda_elc = log_lambda_elc,
    .nu_elc = nu_elc,
    .log_lambda_ion = log_lambda_ion,
    .nu_ion = nu_ion,
    .c_s = c_s,
    .vte = vte,
    .vti = vti,
    .vtAr = vtAr,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .k_perp = k_perp,
    .n_src = n_src,
    .T_src = T_src,
    .cdim = cdim,
    .vdim = vdim,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .Lx = Lx,
    .Lz = Lz,
    .vscale = vscale,
    .vcut = vcut,
    .mucut = mucut,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_elc = mu_max_elc,
    .vpar_max_ion = vpar_max_ion,
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  const char f_dist_filename[20] = "example_dist/fe.npy"; //Has to be binary
  const char v_grid_filename[23] = "example_dist/vgrid.npy"; //Has to be binary
  const char te_filename[20] = "example_dist/Te.npy"; //Has to be binary
  size_t num_te = 199;
  size_t num_v = 80;
  FILE *feptr = fopen(f_dist_filename, "rb");
  FILE *vptr = fopen(v_grid_filename, "rb");
  FILE *teptr = fopen(te_filename, "rb");
  if (feptr == NULL) {
    printf("Feptr = NULL\n");
    exit(EXIT_FAILURE);
  }
  if (teptr == NULL) {
    printf("teptr = NULL\n");
    exit(EXIT_FAILURE);
  }
  if (vptr == NULL) {
    printf("vptr = NULL\n");
    exit(EXIT_FAILURE);
  }
  ctx.num_te = num_te;
  ctx.num_v = num_v;
  ctx.te_ind = 0;
  ctx.f_from_file = (double *)malloc(num_te*num_v*sizeof(double));
  ctx.v_from_file = (double *)malloc(num_v*sizeof(double));
  ctx.te_from_file = (double *)malloc(num_te*sizeof(double));
  long res_sz1 = fread(ctx.f_from_file, 1, num_te*num_v*sizeof(double), feptr);
  long res_sz2 = fread(ctx.v_from_file, 1, num_v*sizeof(double), vptr);
  long res_sz3 = fread(ctx.te_from_file, 1, num_te*sizeof(double), teptr);
  /*for (int i=0; i<num_te; i++)
    printf("i=%d, fe[v=0, Te=%e]=%e\n",i, ctx.te_from_file[i],ctx.f_from_file[i]);

  for (int i=0; i<num_v; i++)
    printf("i=%d, fe[v=%e, Te=6eV]=%e\n",i, ctx.v_from_file[i],ctx.f_from_file[i*num_te]);
  */
  fclose(feptr);  // Close file
  fclose(teptr);  // Close file
  fclose(vptr);  // Close file
  
  //free(f_dist_filename);
  return ctx;
}

void mapc2p_vel_elc(double t, const double* vc, double* GKYL_RESTRICT vp, void* ctx)
{
    struct sheath_ctx* app = ctx;

    double scale = app->vscale; //0.2;
    double vcut = app->vcut/scale; //5e6/scale;
    double mucut = app->mucut/scale; //5e-17/scale;
    double cvpar = vc[0], cmu = vc[1];
    //  vp[0] = cvpar;
    if (app->vpar_max_elc * cvpar < -vcut) {
      vp[0] = -fabs(app->vpar_max_elc * pow(cvpar, 1)) + (1-scale)*vcut;
    } else if (cvpar < 0.0) {
      vp[0] = -fabs(app->vpar_max_elc * scale * pow(cvpar, 1));
    }else if (app->vpar_max_elc * cvpar < vcut) {
      vp[0] = app->vpar_max_elc * scale * pow(cvpar, 1);
    } else {
      vp[0] = fabs(app->vpar_max_elc * cvpar) - (1-scale)*vcut;
    }
    vp[1] = app->mu_max_elc * pow(cmu, 2);// - mucut*(1-scale);
}

void mapc2p_vel_ion(double t, const double* vc, double* GKYL_RESTRICT vp, void* ctx)
{
    struct sheath_ctx* app = ctx;

    double cvpar = vc[0], cmu = vc[1];
    //  vp[0] = cvpar;
    if (cvpar < 0.)
      vp[0] = -fabs(app->vpar_max_ion * pow(cvpar, 1));
    else
      vp[0] = app->vpar_max_ion * pow(cvpar, 1);

    //  vp[1] = cmu;
    vp[1] = app->mu_max_ion * pow(cmu, 1);
}

void
evalSourceDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n_src = app->n_src;
  double cx= app->cx;
  double cz= app->cz;
  double x_center= app->x_center;
  double Lz = app->Lz;

  double n = 0.0;

  n = exp( -z*z/2/cz/cz );
  if (n < 1e-5)
    n = 1e-5;
  n = n*n_src;
  fout[0] = n;
}

void
evalSourceUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set source parallel velocity.
  fout[0] = 0.0;
}

void
evalSourceTempInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double T_src = app->T_src;
  fout[0] = T_src;
}

void
evalDensityInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[0];

  double n0= app->n0;
  double cx= app->cx;
  double cz= app->cz;
  double x_center= app->x_center;
  double Lz = app->Lz;

  double n = 0.0;

  n = exp( -z*z/2/cz/cz );
  if (n < 1e-1)
    n = 1e-1;
  n = n*n0;
  fout[0] = n;
}

void
eval_density_low(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double z = xn[0];
  double n0 = app->n0low;
  double cz = app->cz/1.4/2.0;
  double zcenter = app->Lz/2;
  double n = 0.0;
  if (z>0)
    n = n0 * exp(-(z-zcenter)*(z-zcenter)/(2.0*cz*cz));
  else
    n = n0 * exp(-(z+zcenter)*(z+zcenter)/(2.0*cz*cz));
  if (n < 1.0e8)
    n = 1.0e8;

  fout[0] = n;
}

void
eval_density_high(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
    struct sheath_ctx* app = ctx;
    double z = xn[0];
    double n0 = app->n0high;
    double cz = app->cz / 1.4 / 2.0;
    double zcenter = app->Lz/2;
    double n = 0.0;
    if (z > 0)
        n = n0 * exp(-(z - zcenter) * (z - zcenter) / (2.0 * cz * cz));
    else
        n = n0 * exp(-(z + zcenter) * (z + zcenter) / (2.0 * cz * cz));
    if (n < 1.0e8)
        n = 1.0e8;

    fout[0] = n;
}

void
eval_density_arion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0e5;
}

void
evalUparInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Set parallel velocity.
  fout[0] = 0.0;
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

double interpolate1D(double x, double x1, double x2, double f1, double f2){
  return (x2-x)/(x2-x1)*f1 + (x-x1)/(x2-x1)*f2;
}

double interpolate2D(double vmag, double Te, struct sheath_ctx *app) {
  int num_te = app->num_te;
  int num_v = app->num_v;
  int te_ind = app->te_ind;
  double *f = app->f_from_file;
  double *vFile = app->v_from_file;
  double *teFile = app->te_from_file;

  if (vmag > vFile[num_v-1])     
    return 0.0;

  int x1=0, x2=0, y1=0, y2=0;
  for (int i=num_te-1; i>=0; i--) {    
    if (Te > teFile[i]) {
      x1 = i-1;
      x2 = i;
    }
    //printf("i=%d, Te=%e, x1=%d, x2=%d, num_te=%d, teFile[i]=%e\n ", i, Te, x1,x2,num_te,teFile[i]);
  }
  if (vmag < vFile[0]) {
    printf("vmag=%e, vfile=%e, f=%e\n",vmag, vFile[0],f[num_te*x1]);
    return f[num_te*te_ind];
  }
  
  for (int i=1; i<num_v; i++) {
    if (vmag > vFile[i-1]) {
      y1 = i-1;
      y2 = i;
    }
    //    printf("i=%d, vmag=%e, y1=%d, y2=%d, num_v=%d, vFile[i]=%e, f[%d]=%e\n ", i, vmag, y1,y2,num_v,vFile[i],num_te*te_ind+i,f[num_te*i]);
  }
  /*  if(x1==x2) 
    return interpolate1D(vmag, vFile[y1], vFile[y2], f[num_te*y1+x1], f[num_te*y1+x2]);
    
  double z1=interpolate1D(Te, teFile[x1], teFile[x2], f[num_te*y1+x1], f[num_te*y2+x1]);
  double z2=interpolate1D(Te, teFile[x1], teFile[x2], f[num_te*y1+x2], f[num_te*y2+x2]);
  double ffinal = interpolate1D(vmag, vFile[y1], vFile[y2], z1, z2);*/
  double ffinal = interpolate1D(vmag, vFile[y1], vFile[y2], f[num_te*y1+te_ind], f[num_te*y2+te_ind]);
  return fmax(ffinal,0.0);
}


void
evalFuncElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x=xn[0], v=xn[1], mu=xn[2];
  double mass = app->mass_elc;
  double vmag = pow(v*v+2*mu*app->B0/mass,1.0/2.0);
  double f;
  //
  f = interpolate2D(vmag, x, app); // Te is set to x
  if (f!=f)
    printf("x=%f, vmag=%e, v=%e, mu=%e, f=%e\n",x, vmag, v, mu, f);
  fout[0]=f;
}

void
evalTempIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double Ti = app->Ti;
  fout[0] = Ti;
}

void
eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct sheath_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void
evalNuElcInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_elc = app->nu_elc;

  // Set electron collision frequency.
  fout[0] = nu_elc;
}

void
evalNuIonInit(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;

  double nu_ion = app->nu_ion;

  // Set ion collision frequency.
  fout[0] = nu_ion;
}

double plasma_frequency(double n, double m)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  return sqrt(n*eV*eV/m/eps0);
}
double coulomb_log(double ns, double nr, double ms, double mr, double Ts, double Tr, double qs, double qr)
{

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double hbar = GKYL_PLANCKS_CONSTANT_H/2/M_PI;
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  double wps = plasma_frequency(ns,ms);
  double wpr = plasma_frequency(nr,mr);
  double inner1 = wps*wps/(Ts/ms + 3*Ts/ms) + wpr*wpr/(Tr/mr + 3*Ts/ms);
  double u = 3*(vts*vts + vtr*vtr);
  double msr = ms*mr/(ms+mr);
  double inner2 = fmax(fabs(qs*qr)/(4*M_PI*eps0*msr*u*u), hbar/(2*sqrt(eV)*msr*u));
  double inner = (1/inner1)*(1/inner2/inner2) + 1;
  return 0.5*log(inner);
}

double norm_nu_func(double nu_frac, double ns, double nr, double ms, double mr, double qs, double qr, double Ts, double Tr)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double clog = coulomb_log(ns,nr,ms,mr,Ts, Tr, qs, qr);
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  return nu_frac/ms*(1/mr+1/ms)*qs*qs*qr*qr*clog/(6*pow(M_PI,1.5)*eps0*eps0);
}

static inline void
mapc2p(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT xp, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double x = zc[0], y = zc[1], z = zc[2];
  xp[0] = x; xp[1] = y; xp[2] = z;
}

void bmag_func(double t, const double* GKYL_RESTRICT zc, double* GKYL_RESTRICT fout, void* ctx)
{
  struct sheath_ctx *app = ctx;
  double B0 = app->B0;
  fout[0] = B0;
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
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);

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

  ctx.te_ind = APP_ARGS_CHOOSE(app_args.vcells[2], ctx.te_ind);
  printf("Using the %dth temperature in example_dist/te.cxv\n",ctx.te_ind);
  
  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);
  
  // Electron species.
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .no_by = true,

    /*    .mapc2p = {
       .mapping = mapc2p_vel_elc,
       .ctx = &ctx,
       },*/

    .projection = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = evalFuncElcInit,
      .ctx_func = &ctx,
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nu_frac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuElcInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    /* .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar = evalSourceUparInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempInit,
        .ctx_temp = &ctx,
      }, 
      },*/

    .bcx = {
      .lower = { .type = GKYL_SPECIES_GK_SHEATH, },
      .upper = { .type = GKYL_SPECIES_GK_SHEATH, },
    },

    .radiation = {
      .radiation_id = GKYL_GK_RADIATION, 
      .num_cross_collisions = 6, 
      .collide_with = { "Li0", "Li1", "Li2", "Ar8", "Ar9", "Ar11"},
      .z = {3, 3, 3, 18, 18, 18},
      .charge_state = {0, 1, 2, 8, 9, 11},
      .num_of_densities = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // Must be 1 for now
      },

    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ion species.
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .lower = { -1.0, 0.0},
    .upper = {  1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0, 
    .no_by = true, 

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    
    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .density = evalDensityInit,
      .ctx_density = &ctx,
      .upar = evalUparInit,
      .ctx_upar = &ctx,
      .temp = evalTempIonInit,
      .ctx_temp = &ctx,
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nu_frac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .self_nu = evalNuIonInit,
      .ctx = &ctx,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },
    /*    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = evalSourceDensityInit,
        .ctx_density = &ctx,
        .upar = evalSourceUparInit,
        .ctx_upar = &ctx,
        .temp = evalSourceTempInit,
        .ctx_temp = &ctx,
      }, 
      },*/

    .bcx = {
      .lower = { .type = GKYL_SPECIES_GK_SHEATH, },
      .upper = { .type = GKYL_SPECIES_GK_SHEATH, },
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ar+8
  struct gkyl_gyrokinetic_neut_species Ar8 = {
    .name = "Ar8", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_low,
      .ctx_upar = &ctx,
      .udrift = eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };

  // Ar+9
  struct gkyl_gyrokinetic_neut_species Ar9 = {
    .name = "Ar9", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_low,
      .ctx_upar = &ctx,
      .udrift = eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };

  // Ar +11
  struct gkyl_gyrokinetic_neut_species Ar11 = {
    .name = "Ar11", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_high,
      .ctx_upar = &ctx,
      .udrift= eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };
  
  // neutral Li
  struct gkyl_gyrokinetic_neut_species Li0 = {
    .name = "Li0", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_high,
      .ctx_upar = &ctx,
      .udrift = eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };

  // Li+1
  struct gkyl_gyrokinetic_neut_species Li1 = {
    .name = "Li1", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_low,
      .ctx_upar = &ctx,
      .udrift = eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };

  // Li+2
  struct gkyl_gyrokinetic_neut_species Li2 = {
    .name = "Li2", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { 16, 16, 16},
    .is_static = true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = eval_density_low,
      .ctx_upar = &ctx,
      .udrift = eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,
    },

    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},
  };

  // Field.
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .kperpSq = ctx.k_perp * ctx.k_perp,
  };
  
  // GK app.
  struct gkyl_gk app_inp = {
    .name = "non_maxwellian_elc",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .lower = { -ctx.Lz/2 },
    .upper = { ctx.Lz/2 },
    .cells = { cells_x[0] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {ctx.Lx/2.0, 0.0},
      .mapc2p = mapc2p,
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func,
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = { 1 },

    .num_species = 2,
    .species = { elc, ion},
    .num_neut_species = 6,
    .neut_species = {Li0, Li1, Li2, Ar8, Ar9, Ar11},
    .field = field,


    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1], app_args.cuts[2] },
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
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

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
        write_data(&trig_write, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

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

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_comm_release(comm);
  free(ctx.f_from_file);
  free(ctx.v_from_file);

  mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
