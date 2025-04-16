#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_pkpm.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

#include <rt_arg_parse.h>

struct pkpm_kalf_ctx {
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
  int num_frames;
  double min_dt;
  double th;
  double k;
  bool use_gpu;
};

static inline double
maxwellian(double n, double v, double vth)
{
  double v2 = v*v;
  return n/sqrt(2*M_PI*vth*vth)*exp(-v2/(2*vth*vth));
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;
  
  double x = xn[0], vx = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  
  double fv = maxwellian(app->n0, vx, app->vtElc);
    
  fout[0] = fv;
  fout[1] = app->vtElc*app->vtElc*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;
  
  double x = xn[0], vx = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;

  double fv = maxwellian(app->n0, vx, app->vtIon);
    
  fout[0] = fv;
  fout[1] = app->vtIon*app->vtIon*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;
  
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double L = app->L;
  double k = app->k;
  double th = app->th;   
  double uxe = app->uxe;
  double uye = app->uye;
  double uze = app->uze;
  double uxePhi = app->uxePhi;
  double uyePhi = app->uyePhi;
  double uzePhi = app->uzePhi;   

  double u_xe = uxe*cos(k*x+uxePhi)*cos(th) - uze*cos(k*x+uzePhi)*sin(th);
  double u_ye = uye*cos(k*x+uyePhi);
  double u_ze = uxe*cos(k*x+uxePhi)*sin(th) + uze*cos(k*x+uzePhi)*cos(th);

  fout[0] = me*u_xe;
  fout[1] = me*u_ye;
  fout[2] = me*u_ze;  
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;
  
  double x = xn[0];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double L = app->L;
  double k = app->k;
  double th = app->th;    
  double uxi = app->uxi;
  double uyi = app->uyi;
  double uzi = app->uzi;
  double uxiPhi = app->uxiPhi;
  double uyiPhi = app->uyiPhi;
  double uziPhi = app->uziPhi;

  double u_xi = uxi*cos(k*x+uxiPhi)*cos(th) - uzi*cos(k*x+uziPhi)*sin(th);
  double u_yi = uyi*cos(k*x+uyiPhi);
  double u_zi = uxi*cos(k*x+uxiPhi)*sin(th) + uzi*cos(k*x+uziPhi)*cos(th);

  fout[0] = mi*u_xi;
  fout[1] = mi*u_yi;
  fout[2] = mi*u_zi;    // end new
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;

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
  struct pkpm_kalf_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_kalf_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_kalf_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // pemiability of free space

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.0; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 2.0; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density
  double vAe = 0.8;
  double beta = 0.01;

  double B0 = vAe*sqrt(mu0*n0*massElc);
  double vtElc = vAe*sqrt(Te_Ti*beta/2.0);
  // ion velocities
  double vAi = vAe/sqrt(massIon);
  double vtIon = vtElc/sqrt(massIon * Te_Ti); 

  // ion cyclotron frequency and gyroradius
  double omegaCi = chargeIon*B0/massIon;
  double di = vAi/omegaCi;
  double rhoi = sqrt(2.)*vtIon/omegaCi;

  // collision frequencies
  double nuElc = 1.0e-4*omegaCi;
  double nuIon = nuElc/sqrt(massIon);

  // initial conditions
  double a = 1.e-2; // modified from 1./3.
  double delta_B0 = a*B0;
  double delta_u0 = a*vAi;

  // kperp rhoi = 0.1
  double kperp = 0.1037786 / rhoi;
  double kpar = 0.01 / rhoi; // Theta = 84.5 degrees

  double Bx = 8.137698648752702e-05*delta_B0;  // new
  double By = 0.1008674311076785*delta_B0;
  double Bz = 0.0008445189786624518*delta_B0;
  double BxPhi = -1.56984935634485;
  double ByPhi = -0.0016123288625174861;
  double BzPhi = 1.5717432972475844;

  double uxi = 0.009998765125735562*delta_u0;
  double uyi = 0.10020010001266386*delta_u0;
  double uzi = 0.00105318096901496*delta_u0;
  double uxiPhi = -1.5724166575943115;
  double uyiPhi = 3.1415767548032663;
  double uziPhi = -1.8712400742477489;

  double uxe = 8.698077086693888e-05*delta_u0;
  double uye = 0.10108470000282252*delta_u0;
  double uze = 0.10568579211483767*delta_u0;
  double uxePhi = 1.5701036466287228;
  double uyePhi = 3.1415851806626343;
  double uzePhi = -1.5753424502174473; 

  double k = sqrt(kpar*kpar + kperp*kperp);
  double th = -atan(kpar/kperp); //Angle to rotate k to be along x
  // domain size and simulation time
  double L = 2.*M_PI/k;

  double tend = 320.0/omegaCi; 
  int num_frames = 400;

  struct pkpm_kalf_ctx ctx = {
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
    .Bx = Bx,   // new
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
    .uzePhi = uzePhi,   // end new
    .delta_u0 = delta_u0,
    .delta_B0 = delta_B0,
    .L = L,
    .k = k,
    .th = th,
    .tend = tend,
    .num_frames = num_frames, 
    .min_dt = 1.0e-2, 
  };
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_pkpm_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) 
    gkyl_pkpm_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 112);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 128);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_kalf_ctx ctx = create_ctx(); // context for init functions

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int cells[] = { NX };
  struct gkyl_range globalr;
  gkyl_create_global_range(1, cells, &globalr);
  
  // create decomposition
  int cuts[] = { 1 };
#ifdef GKYL_HAVE_MPI  
  if (app_args.use_mpi) {
    cuts[0] = app_args.cuts[0];
  }
#endif 
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(1, cuts, &globalr);

  // construct communcator for use in app
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  }
  else
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu        
      }
    );
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu      
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_sz;
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = cuts[0];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }
  
  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vtElc},
    .upper = { 6.0 * ctx.vtElc}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncElc,
    .init_fluid = evalFluidElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
    },    

    //.diffusion = {.D = 1.0e-5, .order=4},
  };
  
  // ions
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -6.0 * ctx.vtIon},
    .upper = { 6.0 * ctx.vtIon}, 
    .cells = { VX },

    .ctx_dist = &ctx,
    .ctx_fluid = &ctx,
    .init_dist = evalDistFuncIon,
    .init_fluid = evalFluidIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
    },    

    //.diffusion = {.D = 1.0e-5, .order=4},
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = 1.0, .mu0 = 1.0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_aw_1x_p2_beta001_kperp01_mime1836_TiTe05_A1_i",

    .cdim = 1, .vdim = 1,
    .lower = { 0.0},
    .upper = { ctx.L },
    .cells = { NX},
    .poly_order = 2,
    .basis_type = app_args.basis_type,
    
    .use_explicit_source=true,
    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .use_gpu = app_args.use_gpu,

    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = ctx.tend;
  double dt = tend-tcurr;
  int nframe = ctx.num_frames;
  // create trigger for IO 
  struct gkyl_tm_trigger io_trig = { .dt = tend/nframe };

  // initialize simulation
  gkyl_pkpm_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    if (step % 1000 == 0) {
      gkyl_pkpm_app_calc_field_energy(app, tcurr);
      gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
      gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
    }
    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    if (status.dt_actual < ctx.min_dt) {
      gkyl_pkpm_app_cout(app, stdout, "** Time step crashing! Aborting simulation and writing out last output ....\n");
      gkyl_pkpm_app_write(app, tcurr, 1000);
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_pkpm_app_calc_field_energy(app, tcurr);
  gkyl_pkpm_app_calc_integrated_L2_f(app, tcurr);
  gkyl_pkpm_app_calc_integrated_mom(app, tcurr);
  gkyl_pkpm_app_write_field_energy(app);
  gkyl_pkpm_app_write_integrated_L2_f(app);
  gkyl_pkpm_app_write_integrated_mom(app);
  gkyl_pkpm_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species RHS calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM Vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM Variables (bvar) calculation took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);

  gkyl_pkpm_app_cout(app, stdout, "Species BCs took %g secs\n", stat.species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid Species BCs took %g secs\n", stat.fluid_species_bc_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field BCs took %g secs\n", stat.field_bc_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);
  
  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);

  // simulation complete, free app
  gkyl_pkpm_app_release(app);

  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif  
  
  return 0;
}
