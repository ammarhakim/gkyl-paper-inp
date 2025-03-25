#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>

/*
This is a test to see if our ionization/recombination operators produce coronal equilibrium
Need at least 1e21 argon density otr distribution of charge states may be slightly different
Need ne*time > 10^20

*/

struct gk_app_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massAr; // Li mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double c_s; // sound speed
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac;
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // argon reference density
  double n0Ar1; // argon1 reference density
  double n0Ar2; // argon2 reference density
  double n0Ar3; // argon3 reference density
  double n0Ar4; // argon4 reference density
  double n0Ar5; // argon5 reference density
  double n0Ar6; // argon6 reference density
  double n0Ar7; // argon7 reference density
  double n0Ar8; // argon8 reference density
  double n0Ar9; // argon9 reference density
  double n0Ar10; // argon10 reference density
  double n0Ar11; // argon11 reference density
  double n0Ar12; // argon12 reference density
  double n0Ar13; // argon13 reference density
  double n0Ar14; // argon14 reference density
  double n0Ar15; // argon15 reference density
  double n0Ar16; // argon16 reference density
  double n0Ar17; // argon17 reference density
  double n0Ar18; // argon18 reference density

  // Simulation parameters
  double Lz; // Box size in z.
  double kperp; // perpendicular wave number used in Poisson solve
  double n_src; // Source density.
  double T_src; // Temp density.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Li ions
  double mu_max_Ar; // Velocity space extents in mu for Li ions  
  double finalTime; // end time
  int numFrames; // number of output frames
};

// Source profiles.
void eval_source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_src;
  double Ls = app->Lz/4.;

  if (fabs(z) < Ls) {
    fout[0] = 1.;
  } else {
    fout[0] = 1.e-40;
  }
  fout[0] = n_src*fout[0];
}
void eval_source_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n_src = app->n_src;
  double Ls = app->Lz/4.;

  if (fabs(z) < Ls) {
    fout[0] = 1.;
  } else {
    fout[0] = 1.e-40;
  }
  fout[0] = n_src*fout[0];
}

void eval_source_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double T_src = app->T_src;
  fout[0] = T_src;
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

// Initial conditions.
void eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];

  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0 = app->n0;
  double n0Ararr[9] = {0};
  n0Ararr[0] = app->n0Ar;
  n0Ararr[1] = app->n0Ar1;
  n0Ararr[2] = app->n0Ar2;
  n0Ararr[3] = app->n0Ar3;
  n0Ararr[4] = app->n0Ar4;
  n0Ararr[5] = app->n0Ar5;
  n0Ararr[6] = app->n0Ar6;
  n0Ararr[7] = app->n0Ar7;
  n0Ararr[8] = app->n0Ar8;
  n0Ararr[9] = app->n0Ar9;
  n0Ararr[10] = app->n0Ar10;
  n0Ararr[11] = app->n0Ar11;
  n0Ararr[12] = app->n0Ar12;
  n0Ararr[13] = app->n0Ar13;
  n0Ararr[14] = app->n0Ar14;
  n0Ararr[15] = app->n0Ar15;
  n0Ararr[16] = app->n0Ar16;
  n0Ararr[17] = app->n0Ar17;
  n0Ararr[18] = app->n0Ar18;
      
  // fout[0] = n0 - 10*n0Ar;
  for (int i=0; i<=8; i++)
    n0 -= n0Ararr[i] * i;
  fout[0] = n0;
}

void eval_density_ar0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar;
  fout[0] = n0Ar;
}

void eval_density_ar1(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar1;
  fout[0] = n0Ar;
}

void eval_density_ar2(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar2;
  fout[0] = n0Ar;
}

void eval_density_ar3(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar3;
  fout[0] = n0Ar;
}

void eval_density_ar4(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar4;
  fout[0] = n0Ar;
}

void eval_density_ar5(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar5;
  fout[0] = n0Ar;
}

void eval_density_ar6(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar6;
  fout[0] = n0Ar;
}

void eval_density_ar7(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar7;
  fout[0] = n0Ar;
}

void eval_density_ar8(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar8;
  fout[0] = n0Ar;
}

void eval_density_ar9(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar9;
  fout[0] = n0Ar;
}

void eval_density_ar10(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar10;
  fout[0] = n0Ar;
}

void eval_density_ar11(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar11;
  fout[0] = n0Ar;
}

void eval_density_ar12(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar12;
  fout[0] = n0Ar;
}

void eval_density_ar13(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar13;
  fout[0] = n0Ar;
}

void eval_density_ar14(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar14;
  fout[0] = n0Ar;
}

void eval_density_ar15(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar15;
  fout[0] = n0Ar;
}

void eval_density_ar16(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar17;
  fout[0] = n0Ar;
}

void eval_density_ar17(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar17;
  fout[0] = n0Ar;
}

void eval_density_ar18(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  struct gk_app_ctx *app = ctx;
  double n0Ar = app->n0Ar18;
  fout[0] = n0Ar;
}

void eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Te = app->Te;
  fout[0] = Te;
}

void eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double Ti = app->Ti;
  fout[0] = Ti;
}

void eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_app_ctx *app = ctx;
  fout[0] = app->nuIon;
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

double norm_nu_func(double nuFrac, double ns, double nr, double ms, double mr, double qs, double qr, double Ts, double Tr)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double clog = coulomb_log(ns,nr,ms,mr,Ts, Tr, qs, qr);
  double vts = sqrt(Ts/ms);
  double vtr = sqrt(Tr/mr);
  return nuFrac/ms*(1/mr+1/ms)*qs*qs*qr*qr*clog/(6*pow(M_PI,1.5)*eps0*eps0);
}

struct gk_app_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // Proton mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.

  double mi = 2.014*mp; // D ion mass
  double mAr = 14*GKYL_PROTON_MASS; // Ar ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  // Reference density and temperature.
  double Te = 2*eV;
  double Ti = 2*eV;
  double TAr = 2*eV;
  double n0 = 1.0e27;
  //double n0Ar = n0 * 1e-18;
  double n0Ar = n0 * 1e-24;
  double n0Ar1 = n0Ar * 1e3;
  double n0Ar2 = n0Ar * 2e6;
  double n0Ar3 = n0Ar * 1e8;
  double n0Ar4 = n0Ar * 8e9;
  double n0Ar5 = n0Ar * 3e9;
  double n0Ar6 = n0Ar * 2e9;
  double n0Ar7 = n0Ar * 3e8;
  double n0Ar8 = n0Ar * 2e8;
  double n0Ar9 = n0Ar * 2e5;
  double n0Ar10 = n0Ar * 2e3;
  double n0Ar11 = n0Ar * 2e1;
  double n0Ar12 = n0Ar * 2e0;
  double n0Ar13 = n0Ar;
  double n0Ar14 = n0Ar;
  double n0Ar15 = n0Ar;
  double n0Ar16 = n0Ar;
  double n0Ar17 = n0Ar;
  double n0Ar18 = n0Ar;

  // Geometry and magnetic field.
  double B_axis = 0.5;
  double R0     = 0.85;
  double a0     = 0.15;
  double R      = R0 + a0;
  double B0     = B_axis*(R0/R);

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(Ti/mAr);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq
  
  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lz = 4.;

  // Perpendicular wavenumber in SI units:
  double kperpRhos = 0.3;
  double kperp = kperpRhos / rho_s;

  // Source parameters.
  //double n_src = 2.870523e+21;
  //double n_src = 2.870523e2*n0;
  double n_src = 0.0;
  double T_src = 2.*Te;

  double vpar_max_elc = 6.0*vtElc*1.4;
  double mu_max_elc = (3./2.)*0.5*me*pow(4.0*vtElc,2)/(2.0*B0)*2;

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = (3./2.)*0.5*mi*pow(4.0*vtIon,2)/(2.0*B0);

  double vpar_max_Ar = 6.0*vtAr;
  double mu_max_Ar = (3./2.)*0.5*mAr*pow(4.0*vtAr,2)/(2.0*B0);
  
  double finalTime = 2e-7; 
  double numFrames = 100;

  struct gk_app_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti,
    .TAr = TAr,
    .c_s = c_s, 
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac, 
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar,
    .n0Ar1 = n0Ar1,
    .n0Ar2 = n0Ar2,
    .n0Ar3 = n0Ar3,
    .n0Ar4 = n0Ar4,
    .n0Ar5 = n0Ar5,
    .n0Ar6 = n0Ar6,
    .n0Ar7 = n0Ar7,
    .n0Ar8 = n0Ar8,
    .n0Ar9 = n0Ar9,
    .n0Ar10 = n0Ar10,
    .n0Ar11 = n0Ar11,
    .n0Ar12 = n0Ar12,
    .n0Ar13 = n0Ar13,
    .n0Ar14 = n0Ar14,
    .n0Ar15 = n0Ar15,
    .n0Ar16 = n0Ar16,
    .n0Ar17 = n0Ar17,
    .n0Ar18 = n0Ar18,
    .Lz = Lz, 
    .kperp = kperp, 
    .n_src = n_src,
    .T_src = T_src,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar, 
    .finalTime = finalTime, 
    .numFrames = numFrames,
  };
  return ctx;
}

struct all_of_1species_reactions {
  struct gkyl_gyrokinetic_react_type recombination[GKYL_MAX_SPECIES];
  struct gkyl_gyrokinetic_react_type ionization[GKYL_MAX_SPECIES];
};

struct all_of_1species_reactions*
create_isonuclear_react(double ion_mass, double electron_mass, int max_z, enum gkyl_ion_type ion_id, char* electron_name, char* ion_name)
{
  static struct all_of_1species_reactions reactions[4];  // index = gkyl_react_self_type
  // 0 = elc, 1=ion, 2 = donor, 3 = recvr
  char ion_nm[128];
  char donor_nm[128];  // Same as recvr_nm
  for (int j=0; j<4; j++) {
    for (int i=0; i<max_z; i++) {
      snprintf(ion_nm, sizeof(ion_nm), "%s%d", ion_name, i+1);
      snprintf(donor_nm, sizeof(donor_nm), "%s%d", ion_name, i);
      reactions[j].recombination[i].react_id = GKYL_REACT_RECOMB;
      reactions[j].recombination[i].ion_id = ion_id;
      strcpy(reactions[j].recombination[i].elc_nm, electron_name);
      reactions[j].recombination[i].charge_state = i;
      reactions[j].recombination[i].ion_mass = ion_mass;
      reactions[j].recombination[i].elc_mass = electron_mass;
      strcpy(reactions[j].recombination[i].ion_nm, ion_nm);
      strcpy(reactions[j].recombination[i].recvr_nm, donor_nm);
      reactions[j].recombination[i].type_self = j;
      
      reactions[j].ionization[i].react_id = GKYL_REACT_IZ;
      reactions[j].ionization[i].ion_id = ion_id;
      strcpy(reactions[j].ionization[i].elc_nm, electron_name);
      reactions[j].ionization[i].charge_state = i;
      reactions[j].ionization[i].ion_mass = ion_mass;
      reactions[j].ionization[i].elc_mass = electron_mass;
      strcpy(reactions[j].ionization[i].ion_nm, ion_nm);
      strcpy(reactions[j].ionization[i].donor_nm, donor_nm);
      reactions[j].ionization[i].type_self = j;
    }
  }
  return reactions;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr-1);
    gkyl_gyrokinetic_app_calc_mom(app); gkyl_gyrokinetic_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int  
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_app_ctx ctx = create_ctx(); // context for init functions

  // Reactions for argon isonuclear sequence
  struct all_of_1species_reactions *reactions = create_isonuclear_react(ctx.massAr, ctx.massElc, 6, GKYL_ION_C, "elc", "Ar");

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 2);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 6);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 4);

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { 4*NV, 4*NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massElc, ctx.chargeElc, ctx.chargeElc, ctx.Te, ctx.Te),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massIon, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.Ti), 
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar1, ctx.massElc, ctx.massAr, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar2, ctx.massElc, ctx.massAr, ctx.chargeElc, 2*ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar3, ctx.massElc, ctx.massAr, ctx.chargeElc, 3*ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar4, ctx.massElc, ctx.massAr, ctx.chargeElc, 4*ctx.chargeIon, ctx.Te, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar5, ctx.massElc, ctx.massAr, ctx.chargeElc, 5*ctx.chargeIon, ctx.Te, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar6, ctx.massElc, ctx.massAr, ctx.chargeElc, 6*ctx.chargeIon, ctx.Te, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 7,
      .collide_with = { "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6"},
      //.num_cross_collisions = 3,
      //.collide_with = { "ion", "Ar1", "Ar2" },
    },

    /*    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
      },*/

    .react_neut = {
      .num_react = 2,
      .react_type = {reactions[GKYL_SELF_ELC].ionization[0], 
		     reactions[GKYL_SELF_ELC].recombination[0]},
    }, 

    .react= { 
      .num_react = 12, 
      .react_type = {reactions[GKYL_SELF_ELC].ionization[1], 
		     reactions[GKYL_SELF_ELC].recombination[1],
		     reactions[GKYL_SELF_ELC].ionization[2], 
		     reactions[GKYL_SELF_ELC].recombination[2],
		     reactions[GKYL_SELF_ELC].ionization[3], 
		     reactions[GKYL_SELF_ELC].recombination[3],
		     reactions[GKYL_SELF_ELC].ionization[4], 
		     reactions[GKYL_SELF_ELC].recombination[4],
		     reactions[GKYL_SELF_ELC].ionization[5], 
		     reactions[GKYL_SELF_ELC].recombination[5],
		     reactions[GKYL_SELF_ELC].ionization[6], 
		     reactions[GKYL_SELF_ELC].recombination[6],
		     },      
    },  
    /*.radiation = {
       .radiation_id = GKYL_GK_RADIATION,
       .num_cross_collisions = 18,
       .collide_with = { "Ar0", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17"},
       .z = { 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 },
       .charge_state = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 },
       .num_of_densities = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
       },*/

    .radiation = {
       .radiation_id = GKYL_GK_RADIATION,
       .num_cross_collisions = 6,
       .collide_with = { "Ar0", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5"},
       .z = { 6, 6, 6, 6, 6, 6},
       .charge_state = { 0, 1, 2, 3, 4, 5},
       .num_of_densities = { 1, 1, 1, 1, 1, 1},
       },
    
    /*.bcx = {
          .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
          .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
	  },*/
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massIon, ctx.chargeIon, ctx.chargeIon, ctx.Ti, ctx.Ti),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.Ti, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar1, ctx.massIon, ctx.massAr, ctx.chargeIon, ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar2, ctx.massIon, ctx.massAr, ctx.chargeIon, 2*ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar3, ctx.massIon, ctx.massAr, ctx.chargeIon, 3*ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar4, ctx.massIon, ctx.massAr, ctx.chargeIon, 4*ctx.chargeIon, ctx.Ti, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar5, ctx.massIon, ctx.massAr, ctx.chargeIon, 5*ctx.chargeIon, ctx.Ti, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar6, ctx.massIon, ctx.massAr, ctx.chargeIon, 6*ctx.chargeIon, ctx.Ti, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6"},
      //.num_cross_collisions = 3,
      //.collide_with = { "elc", "Ar1", "Ar2" },
    },

    /*    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_source_density_ion,
        .ctx_upar = &ctx,
        .upar= eval_source_upar,
        .ctx_temp = &ctx,
        .temp = eval_source_temp,      
      }, 
      },*/
    /*    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar1+ ions
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar1,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar1,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar1, ctx.massAr, ctx.massAr, ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0, ctx.massAr, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0, ctx.massAr, ctx.massIon, ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar2, ctx.massAr, ctx.massAr, ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar3, ctx.massAr, ctx.massAr, ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar4, ctx.massAr, ctx.massAr, ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar5, ctx.massAr, ctx.massAr, ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar1, ctx.n0Ar6, ctx.massAr, ctx.massAr, ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6"},
      //.num_cross_collisions = 3,
      //.collide_with = { "elc", "ion", "Ar2" },
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {reactions[GKYL_SELF_ION].ionization[0], 
		     reactions[GKYL_SELF_ION].recombination[0]},
    },
    .react = {
      .num_react = 2,
      .react_type = {reactions[GKYL_SELF_DONOR].ionization[1], 
		     reactions[GKYL_SELF_RECVR].recombination[1]},
    },
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },


    /*    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar2+ ions
  struct gkyl_gyrokinetic_species Ar2 = {
    .name = "Ar2",
    .charge = 2*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar2,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar2,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar2, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0, ctx.massAr, ctx.massElc, 2*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0, ctx.massAr, ctx.massIon, 2*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar1, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar3, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar4, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar5, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar2, ctx.n0Ar6, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar1", "Ar3", "Ar4", "Ar5", "Ar6"},
      //.num_cross_collisions = 3,
      //.collide_with = { "elc", "ion", "Ar1" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[1], 
		     reactions[GKYL_SELF_ION].recombination[1],
                     reactions[GKYL_SELF_DONOR].ionization[2],
		     reactions[GKYL_SELF_RECVR].recombination[2]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
    },
    */
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar3+ ions
  struct gkyl_gyrokinetic_species Ar3 = {
    .name = "Ar3",
    .charge = 3*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar3,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar3,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar3, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0, ctx.massAr, ctx.massElc, 3*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0, ctx.massAr, ctx.massIon, 3*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar1, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar2, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar4, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar5, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar3, ctx.n0Ar6, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar4", "Ar5", "Ar6"}
    },
    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[2], 
		     reactions[GKYL_SELF_ION].recombination[2],
                     reactions[GKYL_SELF_DONOR].ionization[3],
		     reactions[GKYL_SELF_RECVR].recombination[3]},
    },

	 
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar4+ ions
  struct gkyl_gyrokinetic_species Ar4 = {
    .name = "Ar4",
    .charge = 4*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar4,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar4,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar4, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0, ctx.massAr, ctx.massElc, 4*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0, ctx.massAr, ctx.massIon, 4*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar1, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar2, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar3, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar5, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar4, ctx.n0Ar6, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar5", "Ar6"},
    },
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[3], 
		     reactions[GKYL_SELF_ION].recombination[3] ,
                     reactions[GKYL_SELF_DONOR].ionization[4],
		     reactions[GKYL_SELF_RECVR].recombination[4]},
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar5+ ions
  struct gkyl_gyrokinetic_species Ar5 = {
    .name = "Ar5",
    .charge = 5*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar5,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar5,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar5, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0, ctx.massAr, ctx.massElc, 5*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0, ctx.massAr, ctx.massIon, 5*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar1, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar2, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar3, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar4, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar5, ctx.n0Ar6, ctx.massAr, ctx.massAr, 5*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar6"},
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[4], 
		     reactions[GKYL_SELF_ION].recombination[4] ,
                     reactions[GKYL_SELF_DONOR].ionization[5],
		     reactions[GKYL_SELF_RECVR].recombination[5]},
    },
     
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar6+ ions
  struct gkyl_gyrokinetic_species Ar6 = {
    .name = "Ar6",
    .charge = 6*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar6,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar6,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar6, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0, ctx.massAr, ctx.massElc, 6*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0, ctx.massAr, ctx.massIon, 6*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar1, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar2, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar3, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar4, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar6, ctx.n0Ar5, ctx.massAr, ctx.massAr, 6*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 7,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5"},
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[5], 
	             reactions[GKYL_SELF_ION].recombination[5]} 
    },
	
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar7+ ions
  struct gkyl_gyrokinetic_species Ar7 = {
    .name = "Ar7",
    .charge = 7*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar7,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar7,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar7, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0, ctx.massAr, ctx.massElc, 7*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0, ctx.massAr, ctx.massIon, 7*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar1, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar2, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar3, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar4, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar5, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar6, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar8, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar9, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar10, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar11, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar12, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar13, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar14, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar15, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar16, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar17, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar7, ctx.n0Ar18, ctx.massAr, ctx.massAr, 7*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[6], 
		     reactions[GKYL_SELF_ION].recombination[6] ,
                     reactions[GKYL_SELF_DONOR].ionization[7],
		     reactions[GKYL_SELF_RECVR].recombination[7]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },

    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar8+ ions
  struct gkyl_gyrokinetic_species Ar8 = {
    .name = "Ar8",
    .charge = 8*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar8,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar8,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar8, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0, ctx.massAr, ctx.massElc, 8*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0, ctx.massAr, ctx.massIon, 8*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar1, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar2, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar3, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar4, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar5, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar6, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar7, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar9, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar10, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar11, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar12, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar13, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar14, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar15, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar16, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar17, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar8, ctx.n0Ar18, ctx.massAr, ctx.massAr, 8*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[7], 
		     reactions[GKYL_SELF_ION].recombination[7] ,
                     reactions[GKYL_SELF_DONOR].ionization[8],
		     reactions[GKYL_SELF_RECVR].recombination[8]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar9+ ions
  struct gkyl_gyrokinetic_species Ar9 = {
    .name = "Ar9",
    .charge = 9*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar9,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar9,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar9, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0, ctx.massAr, ctx.massElc, 9*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0, ctx.massAr, ctx.massIon, 9*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar1, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar2, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar3, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar4, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar5, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar6, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar7, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar8, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar10, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar11, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar12, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar13, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar14, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar15, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar16, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar17, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar9, ctx.n0Ar18, ctx.massAr, ctx.massAr, 9*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[8], 
		     reactions[GKYL_SELF_ION].recombination[8] ,
                     reactions[GKYL_SELF_DONOR].ionization[9],
		     reactions[GKYL_SELF_RECVR].recombination[9]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar10+ ions
  struct gkyl_gyrokinetic_species Ar10 = {
    .name = "Ar10",
    .charge = 10*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar10,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar10,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar10, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0, ctx.massAr, ctx.massElc, 10*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0, ctx.massAr, ctx.massIon, 10*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar1, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar2, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar3, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar4, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar5, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar6, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar7, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar8, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar9, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar11, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar12, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar13, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar14, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar15, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar16, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar17, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar10, ctx.n0Ar18, ctx.massAr, ctx.massAr, 10*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[9], 
		     reactions[GKYL_SELF_ION].recombination[9] ,
                     reactions[GKYL_SELF_DONOR].ionization[10],
		     reactions[GKYL_SELF_RECVR].recombination[10]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar11+ ions
  struct gkyl_gyrokinetic_species Ar11 = {
    .name = "Ar11",
    .charge = 11*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar11,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar11,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar11, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0, ctx.massAr, ctx.massElc, 11*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0, ctx.massAr, ctx.massIon, 11*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar1, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar2, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar3, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar4, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar5, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar6, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar7, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar8, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar9, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar10, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar12, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar13, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar14, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar15, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar16, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar17, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar11, ctx.n0Ar18, ctx.massAr, ctx.massAr, 11*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[10], 
		     reactions[GKYL_SELF_ION].recombination[10] ,
                     reactions[GKYL_SELF_DONOR].ionization[11],
		     reactions[GKYL_SELF_RECVR].recombination[11]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar12+ ions
  struct gkyl_gyrokinetic_species Ar12 = {
    .name = "Ar12",
    .charge = 12*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar12,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar12,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar12, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0, ctx.massAr, ctx.massElc, 12*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0, ctx.massAr, ctx.massIon, 12*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar1, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar2, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar3, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar4, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar5, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar6, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar7, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar8, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar9, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar10, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar11, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar13, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar14, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar15, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar16, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar17, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar12, ctx.n0Ar18, ctx.massAr, ctx.massAr, 12*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar13", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[11], 
		     reactions[GKYL_SELF_ION].recombination[11] ,
                     reactions[GKYL_SELF_DONOR].ionization[12],
		     reactions[GKYL_SELF_RECVR].recombination[12]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar13+ ions
  struct gkyl_gyrokinetic_species Ar13 = {
    .name = "Ar13",
    .charge = 13*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar13,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar13,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar13, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0, ctx.massAr, ctx.massElc, 13*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0, ctx.massAr, ctx.massIon, 13*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar1, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar2, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar3, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar4, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar5, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar6, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar7, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar8, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar9, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar10, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar11, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar12, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar14, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar15, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar16, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar17, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar13, ctx.n0Ar18, ctx.massAr, ctx.massAr, 13*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar14", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[12], 
		     reactions[GKYL_SELF_ION].recombination[12] ,
                     reactions[GKYL_SELF_DONOR].ionization[13],
		     reactions[GKYL_SELF_RECVR].recombination[13]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar14+ ions
  struct gkyl_gyrokinetic_species Ar14 = {
    .name = "Ar14",
    .charge = 14*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar14,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar14,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar14, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0, ctx.massAr, ctx.massElc, 14*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0, ctx.massAr, ctx.massIon, 14*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar1, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar2, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar3, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar4, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar5, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar6, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar7, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar8, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar9, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar10, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar11, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar12, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar13, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar15, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar16, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar17, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar14, ctx.n0Ar18, ctx.massAr, ctx.massAr, 14*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar15", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[13], 
		     reactions[GKYL_SELF_ION].recombination[13] ,
                     reactions[GKYL_SELF_DONOR].ionization[14],
		     reactions[GKYL_SELF_RECVR].recombination[14]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar15+ ions
  struct gkyl_gyrokinetic_species Ar15 = {
    .name = "Ar15",
    .charge = 15*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar15,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar15,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar15, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0, ctx.massAr, ctx.massElc, 15*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0, ctx.massAr, ctx.massIon, 15*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar1, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar2, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar3, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar4, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar5, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar6, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar7, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar8, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar9, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar10, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar11, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar12, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar13, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar14, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar16, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar17, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar15, ctx.n0Ar18, ctx.massAr, ctx.massAr, 15*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar16", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[14], 
		     reactions[GKYL_SELF_ION].recombination[14] ,
                     reactions[GKYL_SELF_DONOR].ionization[15],
		     reactions[GKYL_SELF_RECVR].recombination[15]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar16+ ions
  struct gkyl_gyrokinetic_species Ar16 = {
    .name = "Ar16",
    .charge = 16*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar16,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar16,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar16, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0, ctx.massAr, ctx.massElc, 16*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0, ctx.massAr, ctx.massIon, 16*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar1, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar2, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar3, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar4, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar5, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar6, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar7, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar8, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar9, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar10, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar11, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar12, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar13, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar14, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar15, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar17, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar16, ctx.n0Ar18, ctx.massAr, ctx.massAr, 16*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar17", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[15], 
		     reactions[GKYL_SELF_ION].recombination[15] ,
                     reactions[GKYL_SELF_DONOR].ionization[16],
		     reactions[GKYL_SELF_RECVR].recombination[16]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar17+ ions
  struct gkyl_gyrokinetic_species Ar17 = {
    .name = "Ar17",
    .charge = 17*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar17,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar17,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar17, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0, ctx.massAr, ctx.massElc, 17*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0, ctx.massAr, ctx.massIon, 17*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar1, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar2, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar3, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar4, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar5, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar6, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar7, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar8, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar9, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar10, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar11, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar12, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar13, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar14, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar15, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar16, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar17, ctx.n0Ar18, ctx.massAr, ctx.massAr, 17*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar12", "Ar13", "Ar14", "Ar15", "Ar16", "Ar18" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[16], 
		     reactions[GKYL_SELF_ION].recombination[16] ,
                     reactions[GKYL_SELF_DONOR].ionization[17],
		     reactions[GKYL_SELF_RECVR].recombination[17]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };

  // Ar18+ ions
  struct gkyl_gyrokinetic_species Ar18 = {
    .name = "Ar18",
    .charge = 18*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar18,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar18,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar18, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 18*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0, ctx.massAr, ctx.massElc, 18*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0, ctx.massAr, ctx.massIon, 18*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar1, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar2, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar3, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar4, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar5, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 5*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar6, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 6*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar7, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 7*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar8, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 8*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar9, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 9*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar10, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 10*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar11, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 11*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar12, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 12*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar13, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 13*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar14, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 14*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar15, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 15*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar16, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 16*ctx.chargeIon, ctx.TAr, ctx.TAr),
	norm_nu_func(ctx.nuFrac, ctx.n0Ar18, ctx.n0Ar17, ctx.massAr, ctx.massAr, 18*ctx.chargeIon, 17*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 19,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3", "Ar4", "Ar5", "Ar6", "Ar7", "Ar8", "Ar9", "Ar10", "Ar11", "Ar13", "Ar14", "Ar14", "Ar15", "Ar16", "Ar17" },
    },

    .react = {
      .num_react = 4,
      .react_type = {reactions[GKYL_SELF_ION].ionization[18], 
		     reactions[GKYL_SELF_ION].recombination[18]},
    },
    
    .bcx = {
      .lower = { .type = GKYL_SPECIES_ZERO_FLUX, },
      .upper = { .type = GKYL_SPECIES_ZERO_FLUX, },
    },
    /*
    .bcx = {
       .lower = {.type = GKYL_SPECIES_GK_SHEATH, },
       .upper = {.type = GKYL_SPECIES_GK_SHEATH, },
       },*/
    
    .num_diag_moments = 5,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp" },
  };


  // neutral Ar
  struct gkyl_gyrokinetic_neut_species Ar0 = {
    .name = "Ar0", .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, -ctx.vpar_max_Ar, -ctx.vpar_max_Ar},
    .upper = { ctx.vpar_max_Ar, ctx.vpar_max_Ar, ctx.vpar_max_Ar },
    .cells = { NV, NV, NV},
    .is_static = true,
    
    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ar0,
      .ctx_upar = &ctx,
      .udrift= eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    // can we add these diagnostics??
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"}, //, "M2par", "M2perp" },
  }; 
  
  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_PERIODIC, 
    .kperpSq = pow(ctx.kperp, 2.),
  };

  /*struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_ADIABATIC,
    .electron_mass = ctx.massElc,
    .electron_charge = ctx.chargeElc,
    .electron_temp = ctx.Te,
    .bmag_fac = ctx.B0, // Issue here. B0 from soloviev, so not sure what to do. Ours is not constant
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    };*/

  // GK app
  struct gkyl_gk gk = {
    .name = "c1x",

    //.cfl_frac = 0.1,

    .cdim = 1, .vdim = 2,
    .lower = { -ctx.Lz/2.0 },
    .upper = { ctx.Lz/2.0 },
    .cells = { NX },
    .poly_order = 1,
    .basis_type = app_args.basis_type,

    .geometry = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = &ctx,
      .bmag_func = bmag_func, // mapping of computational to physical space
      .bmag_ctx = &ctx
    },

    .num_periodic_dir = 1,
    .periodic_dirs = {0},

    .num_species = 8,
    .species = { elc, ion, Ar1, Ar2, Ar3, Ar4, Ar5, Ar6},
    //.num_species = 4,
    //.species = { elc, ion, Ar1, Ar2 },

    .num_neut_species = 1,
    .neut_species = {Ar0},
    .field = field,



    .use_gpu = app_args.use_gpu,
    .skip_field=true,
  };
  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.finalTime;
  double dt = tend-tcurr;
  int nframe = ctx.numFrames;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    //gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    //gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 1 == 0) {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    }
    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
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

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  
  return 0;
}
