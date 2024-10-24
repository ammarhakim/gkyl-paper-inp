#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <rt_arg_parse.h>
#include <gkyl_tok_geo.h>
#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

struct gk_step_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double massAr; // Argon mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TAr; // Argon temperature
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac;
  double B0; // reference magnetic field
  double n0; // reference density
  double n0Ar; // Argon reference density
  double nsource;
  // Source parameters
  double T_source; // Source electron temperature
  double cx;
  double cz;
  // Simulation parameters
  double Lz; // Box size in z
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double vpar_max_Ar; // Velocity space extents in vparallel for Li ions
  double mu_max_Ar; // Velocity space extents in mu for Li ions    

  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

struct gkyl_tok_geo_efit_inp inp = {
  .filepath = "./data/eqdsk/step.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};


struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .rright= 6.2,
    .rleft= 2.0,
    .rmin = 1.1,
    .rmax = 6.2,
    .zmin = -8.4,
    .zmax = 8.4,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_fixed_z_nodes.gkyl"
  };

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/n0 < 1.0e-1)
    n = n0*1.0e-1;
  fout[0] = n;
}

void
eval_density_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = n0*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/n0 < 1.0e-1)
    n = n0*1.0e-1;

  fout[0] = n ;
}

// This is the plain exponential version
//void
//eval_density_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
//{
//  struct gk_step_ctx *app = ctx;
//  double x = xn[0], z = xn[1];
//  double n0 = app->n0Ar;
//  double cz = app->cz/1.4/6.0;
//  double zcenter = 3.14;
//  double n = 0.0;
//  if (z>0)
//    n = n0 * exp((z-zcenter)/cz);
//  else
//    n = n0 * exp(-(z+zcenter)/cz);
//  if (n < 1.0e8)
//    n = 1.0e8;
//
//  fout[0] = n;
//}
void
eval_density_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double n0 = app->n0Ar;
  double cz = app->cz/1.4/2.0;
  double zcenter = 3.14;
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
eval_density_arion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0e5;
}


void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
}

void
eval_udrift(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0; 
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Te;
  fout[0] = T;
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Ti;
  fout[0] = T;
}

void
eval_temp_ar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->TAr;
  fout[0] = T;
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double x = xn[0], z = xn[1];
  double nsource = app->nsource;
  double cx = app->cx;
  double cz = app->cz;
  double xcenter = 1.2014;
  double n = nsource*exp(-(x-xcenter)*(x-xcenter)/(2.0*cx*cx)) * exp(-z*z/(2.0*cz*cz));
  if (n/nsource < 1.0e-5)
    n = nsource*1.0e-5;
  fout[0] = n;
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double T = app->T_source;
  fout[0] = T;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  fout[0] = app->nuIon;
}


void mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
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

struct gk_step_ctx
create_ctx(void)
{
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  double me = GKYL_ELECTRON_MASS;
  double mAr = 39.95*GKYL_PROTON_MASS; // Ar ion mass
  double qi = eV; // ion charge
  double qe = -eV; // electron charge

  double tempfac = 1.0;
  double Te = 364.0*tempfac*eV;
  double Ti = 546.0*tempfac*eV;
  double TAr = 10.0*eV; 
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/tempfac; // Particle density in 1/m^3
  double n0Ar = n0*0.0001/3.0 * 10.0 * 1.5; // Particle density in 1/m^3

  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtAr = sqrt(TAr/mAr);

  // Source parameters.
  double nsource = 3.9e23/tempfac; // peak source rate in particles/m^3/s 
  double T_source = 1037.0*eV*tempfac;
  double cx = 0.0065612;
  double cz = 0.4916200*1.4;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).
  double Lz = 3.14*2;

  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = 1.5*12.*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = 1.5*12.*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_Ar = 4.0*vtAr;
  double mu_max_Ar = 1.5*12.*mAr*vtAr*vtAr/(2.0*B0);

  double t_end = 16.0e-3; 
  double num_frames = 160 * 10;

  struct gk_step_ctx ctx = {
    .chargeElc = qe, 
    .massElc = me,
    .chargeIon = qi, 
    .massIon = mi,
    .massAr = mAr,
    .Te = Te, 
    .Ti = Ti,
    .TAr = TAr,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0Ar = n0Ar, 
    .T_source = T_source, 
    .nsource = nsource,
    .cx = cx,
    .cz = cz,
    .Lz = Lz, 
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion,
    .vpar_max_Ar = vpar_max_Ar, 
    .mu_max_Ar = mu_max_Ar,
    .t_end = t_end, 
    .num_frames = num_frames, 
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
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_source_mom(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
    gkyl_gyrokinetic_app_write_integrated_source_mom(app);
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

  struct gk_step_ctx ctx = create_ctx(); // context for init functions
  //if (app_args.is_restart) {
  //  ctx.n0Ar = ctx.n0Ar * (2.0 + (double) ((app_args.restart_frame - 0)/ 20) );
  //}

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 40);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[2], 96);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], 16);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], 12);

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int ccells[] = { NX, NZ };
  int cdim = sizeof(ccells)/sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // create decomposition
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d=0; d<cdim; d++)
    cuts[d] = app_args.use_mpi? app_args.cuts[d] : 1;
#else
  for (int d=0; d<cdim; d++) cuts[d] = 1;
#endif
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(cdim, cuts, &cglobal_r);

  // construct communcator for use in app
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
#else
    printf("Using -g and -M together requires NCCL.\n");
    assert( 0 == 1);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  } else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank, comm_sz;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = 1;
  for (int d=0; d<cdim; d++) ncuts *= cuts[d];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }

  for (int d=0; d<cdim-1; d++) {
    if (cuts[d] > 1) {
      if (my_rank == 0)
        fprintf(stderr, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts[d], d);
      goto mpifinalize;
    }
  }

  // electrons
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -ctx.vpar_max_elc, 0.0},
    .upper = { ctx.vpar_max_elc, ctx.mu_max_elc}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,
    //.enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density,
      .ctx_upar = &ctx,
      .upar = eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_elc,      
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massElc, ctx.chargeElc, ctx.chargeElc, ctx.Te, ctx.Te),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massElc, ctx.massIon, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massElc, ctx.massAr, ctx.chargeElc, ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massElc, ctx.massAr, ctx.chargeElc, 2*ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massElc, ctx.massAr, ctx.chargeElc, 3*ctx.chargeIon, ctx.Te, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massElc, ctx.massAr, ctx.chargeElc, 4*ctx.chargeIon, ctx.Te, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 5,
      .collide_with = { "ion", "Ar1", "Ar2", "Ar3", "Ar4" },
    },

    .radiation = {
        .radiation_id = GKYL_GK_RADIATION,
        .num_cross_collisions = 5,
        .collide_with = { "Ar0", "Ar1", "Ar2", "Ar3", "Ar4"},
        .z = {18, 18, 18, 18, 18},
        .charge_state = { 0, 1, 2, 3, 4},
        .num_of_densities = {1, 1, 1, 1, 1},
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_source,      
      }, 
    },
    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1", // ion is always the higher charge state
          .donor_nm = "Ar0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    }, 

    .react= { 
      .num_react = 6, 
      .react_type = { 
    	  { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", // ion is always the higher charge state  
          .donor_nm = "Ar1", // interacts with elc to give up charge 
          .charge_state = 1, // corresponds to lower charge state (donor) 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .recvr_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
    	  { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar3",   
          .donor_nm = "Ar2",  
          .charge_state = 2, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar3", 
          .recvr_nm = "Ar2", 
          .charge_state = 2, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 

    	  { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar4",   
          .donor_nm = "Ar3", 
          .charge_state = 3,  
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB, 
          .type_self = GKYL_SELF_ELC, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar4", 
          .recvr_nm = "Ar3", 
          .charge_state = 3, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 

      },
    },  
    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // ions
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -ctx.vpar_max_ion, 0.0},
    .upper = { ctx.vpar_max_ion, ctx.mu_max_ion}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0,
    //.enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar = eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ion,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massIon, ctx.chargeIon, ctx.chargeIon, ctx.Ti, ctx.Ti),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0, ctx.massIon, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.Ti, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massIon, ctx.massAr, ctx.chargeIon, ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massIon, ctx.massAr, ctx.chargeIon, 2*ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massIon, ctx.massAr, ctx.chargeIon, 3*ctx.chargeIon, ctx.Ti, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0, ctx.n0Ar, ctx.massIon, ctx.massAr, ctx.chargeIon, 4*ctx.chargeIon, ctx.Ti, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 5,
      .collide_with = { "elc", "Ar1", "Ar2", "Ar3", "Ar4" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar = eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_source,      
      }, 
    },
    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
  };

  // Ar1+ ions
  struct gkyl_gyrokinetic_species Ar1 = {
    .name = "Ar1",
    .charge = ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -ctx.vpar_max_Ar, 0.0},
    .upper = { ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV, NMU },
    .polarization_density = ctx.n0Ar*2.0,
    .enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_arion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massElc, ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massIon, ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 5,
      .collide_with = { "elc", "ion", "Ar2", "Ar3", "Ar4" },
    },

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .donor_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar1",
          .recvr_nm = "Ar0",
          .charge_state = 0,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"  },
  };

  // Ar2+ ions
  struct gkyl_gyrokinetic_species Ar2 = {
    .name = "Ar2",
    .charge = 2*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -2*ctx.vpar_max_Ar, 0.0},
    .upper = { 2*ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV*2, NMU },
    .polarization_density = ctx.n0Ar*4.0,
    .enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_arion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massElc, 2*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massIon, 2*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 2*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 5,
      .collide_with = { "elc", "ion", "Ar1", "Ar3", "Ar4" },
    },

    .react = {
      .num_react = 4,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar2", 
          .donor_nm = "Ar1", 
          .charge_state = 1, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar2",
          .recvr_nm = "Ar1",
          .charge_state = 1,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar3", 
          .donor_nm = "Ar2", 
          .charge_state = 2, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar3",
          .recvr_nm = "Ar2",
          .charge_state = 2,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },

      },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"  },
  };
  
  // Ar3+ ions
  struct gkyl_gyrokinetic_species Ar3 = {
    .name = "Ar3",
    .charge = 3*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -3*ctx.vpar_max_Ar, 0.0},
    .upper = { 3*ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV*3, NMU },
    .polarization_density = ctx.n0Ar*8.0,
    .enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_arion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massElc, 3*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massIon, 3*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 3*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 5,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar4" },
    },

        .react = {
      .num_react = 4,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar3", 
          .donor_nm = "Ar2", 
          .charge_state = 2, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar3",
          .recvr_nm = "Ar2",
          .charge_state = 2,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_DONOR, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar4", 
          .donor_nm = "Ar3", 
          .charge_state = 3, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar4",
          .recvr_nm = "Ar3",
          .charge_state = 3,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"  },
  };

  // Ar4+ ions
  struct gkyl_gyrokinetic_species Ar4 = {
    .name = "Ar4",
    .charge = 4*ctx.chargeIon, .mass = ctx.massAr,
    .lower = { -6*ctx.vpar_max_Ar, 0.0},
    .upper = { 6*ctx.vpar_max_Ar, ctx.mu_max_Ar}, 
    .cells = { NV*6, NMU },
    .polarization_density = ctx.n0Ar*16.0,
    .enforce_positivity=true,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
      .ctx_density = &ctx,
      .density = eval_density_arion,
      .ctx_upar = &ctx,
      .upar= eval_upar,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .self_nu_fac = norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 4*ctx.chargeIon, ctx.TAr, ctx.TAr),
      .cross_nu_fac = {
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massElc, 4*ctx.chargeIon, ctx.chargeElc, ctx.TAr, ctx.Te), 
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0, ctx.massAr, ctx.massIon, 4*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.Ti),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 2*ctx.chargeIon, ctx.TAr, ctx.TAr),
        norm_nu_func(ctx.nuFrac, ctx.n0Ar, ctx.n0Ar, ctx.massAr, ctx.massAr, 4*ctx.chargeIon, 3*ctx.chargeIon, ctx.TAr, ctx.TAr)
      },
      .bmag_mid = 2.51,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 5,
      .collide_with = { "elc", "ion", "Ar1", "Ar2", "Ar3" },
    },

    .react = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ, 
          .type_self = GKYL_SELF_ION, 
          .ion_id = GKYL_ION_AR, 
          .elc_nm = "elc", 
          .ion_nm = "Ar4", 
          .donor_nm = "Ar3", 
          .charge_state = 3, 
          .ion_mass = ctx.massAr, 
          .elc_mass = ctx.massElc, 
        }, 
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_AR,
          .elc_nm = "elc",
          .ion_nm = "Ar4",
          .recvr_nm = "Ar3",
          .charge_state = 3,
          .ion_mass = ctx.massAr,
          .elc_mass = ctx.massElc,
        },
      },
    },


    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.03 }, 
      .order = 2, 
    }, 

    .bcx = {
      .lower={.type = GKYL_SPECIES_ABSORB,},
      .upper={.type = GKYL_SPECIES_ABSORB,},
    },
    .bcy = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"  },
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
      .density = eval_density_ar,
      .ctx_upar = &ctx,
      .udrift= eval_udrift,
      .ctx_temp = &ctx,
      .temp = eval_temp_ar,      
    },

    .bcx = { GKYL_SPECIES_ABSORB, GKYL_SPECIES_ZERO_FLUX },
    .bcy = { GKYL_SPECIES_ZERO_FLUX, GKYL_SPECIES_ZERO_FLUX },
    
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"}, //, "M2par", "M2perp" },
  };

  // field
  struct gkyl_gyrokinetic_field field = {
    .bmag_fac = ctx.B0, 
    .fem_parbc = GKYL_FEM_PARPROJ_NONE, 
    .poisson_bcs = {.lo_type = {GKYL_POISSON_DIRICHLET}, 
                    .up_type = {GKYL_POISSON_DIRICHLET}, 
                    .lo_value = {0.0}, .up_value = {0.0}}, 
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "sn13",

    .cdim = 2, .vdim = 2,
    .lower = { 1.0677, -ctx.Lz/2.0 },
    .upper = { 1.3351, ctx.Lz/2.0 },
    .cells = { NX, NZ },
    .poly_order = 1,
    .basis_type = app_args.basis_type,


    .geometry = {
      .geometry_id = GKYL_GEOMETRY_FROMFILE,
      //.world = {0.0},
      //.geometry_id = GKYL_TOKAMAK,
      //.tok_efit_info = &inp,
      //.tok_grid_info = &ginp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 6,
    .species = { elc, ion, Ar1, Ar2, Ar3, Ar4},
    .num_neut_species = 1,
    .neut_species = {Ar0},
    
    .field = field,

    .use_gpu = app_args.use_gpu,
    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // create app object
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);
  // Initial and final simulation times.
  int frame_curr = 0;
  double tcurr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }


    frame_curr = status.frame;
    tcurr = status.stime;

    // project Ar0 because we want to be able to change density
    gkyl_gyrokinetic_app_apply_ic_neut_species(app, 0, 0.0);

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", tcurr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, tcurr);
  }  

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = tcurr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = tcurr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, tcurr, false);
  write_data(&trig_write, app, tcurr, false);

  double dt = t_end-tcurr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((tcurr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, tcurr, tcurr > t_end);
    write_data(&trig_write, app, tcurr, tcurr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, tcurr, true);
        write_data(&trig_write, app, tcurr, true);
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
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);

  mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
