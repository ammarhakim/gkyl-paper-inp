#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct ic_ctx {
double epsilon0;
  double mu0;
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double n0;
  double vAe;
  double B0;
  double beta;
  double vtElc;
  double vtIon;
  double nuElc;
  double nuIon;
  double delta_u0;
  double delta_B0;
  double Bx;
  double By;
  double Bz;
  double uxi;
  double uyi;
  double uzi;
  double uxe;
  double uye;
  double uze;
  double BxPhi;
  double ByPhi;
  double BzPhi;
  double uxiPhi;
  double uyiPhi;
  double uziPhi;
  double uxePhi;
  double uyePhi;
  double uzePhi;
  double L;
  double tend;
  double nframes;
  double min_dt;
  double th;
  double k;
  bool use_gpu;

};

static inline double sq(double x) { return x*x; }

static inline double
maxwellian(double n, double vx, double ux, double vy, double uy, double vz, double uz, double vth)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  return n/pow(sqrt(2*M_PI*vth*vth), 3)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ic_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double n0 = app->n0;
  double k = app->k;
  double th = app->th;
  double vtElc = app->vtElc;
  double uxe = app->uxe;
  double uye = app->uye;
  double uze = app->uze;
  double uxePhi = app->uxePhi;
  double uyePhi = app->uyePhi;
  double uzePhi = app->uzePhi;

  double u_xe = uxe*cos(k*x+uxePhi)*cos(th) - uze*cos(k*x+uzePhi)*sin(th);
  double u_ye = uye*cos(k*x+uyePhi);
  double u_ze = uxe*cos(k*x+uxePhi)*sin(th) + uze*cos(k*x+uzePhi)*cos(th);


  
  fout[0] = maxwellian(n0, vx, u_xe, vy, u_ye, vz, u_ze, vtElc);
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ic_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double n0 = app->n0;
  double k = app->k;
  double th = app->th;
  double vtIon = app->vtIon;
  double uxi = app->uxi;
  double uyi = app->uyi;
  double uzi = app->uzi;
  double uxiPhi = app->uxiPhi;
  double uyiPhi = app->uyiPhi;
  double uziPhi = app->uziPhi;

  double u_xi = uxi*cos(k*x+uxiPhi)*cos(th) - uzi*cos(k*x+uziPhi)*sin(th);
  double u_yi = uyi*cos(k*x+uyiPhi);
  double u_zi = uxi*cos(k*x+uxiPhi)*sin(th) + uzi*cos(k*x+uziPhi)*cos(th);
  
  fout[0] = maxwellian(n0, vx, u_xi, vy, u_yi, vz, u_zi, vtIon);
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ic_ctx *app = ctx;
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double L = app->L;
  double k = app->k;
  double th = app->th;
  double B0 = app->B0;
  double Bx = app->Bx;
  double By = app->By;
  double Bz = app->Bz;
  double BxPhi = app->BxPhi;
  double ByPhi = app->ByPhi;
  double BzPhi = app->BzPhi;
  double uxe = app->uxe;
  double uye = app->uye;
  double uze = app->uze;
  double uxePhi = app->uxePhi;
  double uyePhi = app->uyePhi;
  double uzePhi = app->uzePhi;

  double B_x = -B0*sin(th);
  double B_y = By*cos(k*x+ByPhi);
  double B_z = Bx*cos(k*x+BxPhi)*sin(th) +  (B0+Bz*cos(k*x+BzPhi))*cos(th);

  double u_xe = uxe*cos(k*x+uxePhi)*cos(th) - uze*cos(k*x+uzePhi)*sin(th);
  double u_ye = uye*cos(k*x+uyePhi);
  double u_ze = uxe*cos(k*x+uxePhi)*sin(th) + uze*cos(k*x+uzePhi)*cos(th);

  // E = - v_e x B ~  (J - u) x B
  double E_x = - (u_ye*B_z - u_ze*B_y);
  double E_y = - (u_ze*B_x - u_xe*B_z);
  double E_z = - (u_xe*B_y - u_ye*B_x);
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ic_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct ic_ctx *app = ctx;
  fout[0] = app->nuIon;
}


struct ic_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 100.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 1.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density
  double vAe = 0.1;
  double beta = 1.;

  double B0 = vAe*sqrt(mu0*n0*massElc);
  double vtElc = vAe*sqrt(beta/2.0); //no root(2)!
  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon); //Ti/Te = 1.0

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;
  double rhoi = sqrt(2.)*vtIon/omegaCi;

  // collision frequencies
  double nuElc = 0.001*omegaCi;
  double nuIon = 0.001*omegaCi/sqrt(massIon);

  // initial conditions
  double a = 1.e-3;
  double delta_B0 = a*B0;
  double delta_u0 = a*vAi;
/*  // kpar rhoi = 1
  double kperp = 0.01 / rhoi; 
  double kpar = 1. / rhoi; // Theta_kB = 0.6 degrees
  double Bx = 1.793626*delta_B0;
  double By = 1.793703*delta_B0;
  double Bz = 0.017936*delta_B0;
  double BxPhi = 2.493608;
  double ByPhi = 0.922717;
  double BzPhi = -0.647984;

  double uxi = 2.530374*delta_u0;
  double uyi = 2.530394*delta_u0;
  double uzi = 0.031262*delta_u0;
  double uxiPhi = -0.968017;
  double uyiPhi = -2.538784;
  double uziPhi = 2.347077;

  double uxe = 1.001658*delta_u0;
  double uye = 1.001707*delta_u0;
  double uze = 0.013768*delta_u0;
  double uxePhi = -1.566268;
  double uyePhi = -3.137179;
  double uzePhi = 2.155813;
*/

  // kpar rhoi = 0.4
  double kperp = 0.01 / rhoi; 
  double kpar = 0.4 / rhoi; // Theta_kB = 1.4 degrees
  double Bx = 1.639402*delta_B0;
  double By = 1.641197*delta_B0;
  double Bz = 0.040709*delta_B0;
  double BxPhi = 1.665294;
  double ByPhi = 0.093522;
  double BzPhi = -1.476299;

  double uxi = 1.656575*delta_u0;
  double uyi = 1.657299*delta_u0;
  double uzi = 0.039239*delta_u0;
  double uxiPhi = -1.532881;
  double uyiPhi = -3.103938;
  double uziPhi = 2.115794;

  double uxe = 0.997335*delta_u0;
  double uye = 0.998449*delta_u0;
  double uze = 0.025498*delta_u0;
  double uxePhi = -1.569722;
  double uyePhi = -3.141536;
  double uzePhi = 2.400447;


  double k = sqrt(kpar*kpar + kperp*kperp);
  double th = -atan(kpar/kperp); //Angle to rotate k to be along x
  // domain size and simulation time
  double L = 2.*M_PI/k;
  double tend = 100.0/omegaCi;
  double nframes = 100;
  
  struct ic_ctx ctx = {
  .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .Te_Ti = Te_Ti,
    .n0 = n0,
    .vAe = vAe,
    .B0 = B0,
    .beta = beta,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .Bx = Bx,
    .By = By,
    .Bz = Bz,
    .uxi = uxi,
    .uyi = uyi,
    .uzi = uzi,
    .uxe = uxe,
    .uye = uye,
    .uze = uze,
    .BxPhi = BxPhi,
    .ByPhi = ByPhi,
    .BzPhi = BzPhi,
    .uxiPhi = uxiPhi,
    .uyiPhi = uyiPhi,
    .uziPhi = uziPhi,
    .uxePhi = uxePhi,
    .uyePhi = uyePhi,
    .uzePhi = uzePhi,
    .delta_u0 = delta_u0,
    .delta_B0 = delta_B0,


    .L = L,
    .k = k,
    .th = th,
    .tend = tend,
    .nframes = nframes,
    .min_dt = 1.0e-2, 

  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
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
  struct ic_ctx ctx = create_ctx(); // context for init functions

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0*ctx.vtElc, -6.0*ctx.vtElc, -6.0*ctx.vtElc },
    .upper = { 6.0*ctx.vtElc, 6.0*ctx.vtElc, 6.0*ctx.vtElc }, 
    .cells = { 16, 16, 16 },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2ij" },
  };

    // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0*ctx.vtIon, -6.0*ctx.vtIon, -6.0*ctx.vtIon },
    .upper = { 6.0*ctx.vtIon, 6.0*ctx.vtIon, 6.0*ctx.vtIon }, 
    .cells = { 16, 16, 16 },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2ij" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "ic_damping_3v",

    .cdim = 1, .vdim = 3,
    .lower = { 0 },
    .upper = { ctx.L },
    .cells = { 16 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = ctx.nframes;
  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };
  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);

  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    write_data(&io_trig, app, tcurr);
    step += 1;
  }

  //gkyl_vlasov_app_write(app, tcurr, 1);
  //gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Field RHS calc took %g secs\n", stat.field_rhs_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
