#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_vlasov.h>
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

struct sheath_ctx
{
 int cdim, vdim; // Dimensionality.

 double ep0; // Permittivity of free space.
 double mu0; // Vacuum permeability

 double qe; // Electron charge (C).
 double qi; // Ion charge (C).
 double me; // Electron mass (kg).
 double mi; // Ion mass (kg).

 double ne; // Electron density (m^-3).
 double ni; // Ion density (m^-3).

 double r; // Radius we are looking at (m)
 double I; // Pinch current (A)
 double a; // Pinch radius (m)

 double ue_x; // electron flow speed in x
 double ue_y; // electron flow speed in y
 double ui_x; // ion flow speed in x
 double ui_y; // ion flow speed in y

 double vte; // Electron thermal velocity (m/s).
 double vti; // Ion thermal velocity (m/s).

 double debye; // Debye length (m)

 double Lsrc; // Source length (m).

 double nu_ee; // electron-electron collision frequency
 double nu_ii; // ion-ion collision frequency

 double phi_bias; // bias potential
 double B;       // external magnetic field

 double x_min, x_max;           // Extents of the x grid (m).
 double vx_min_elc, vx_max_elc; // Extents of the electron vx grid (m/s).
 double vy_min_elc, vy_max_elc; // Extents of the electron vy grid (m/s).
 double vx_min_ion, vx_max_ion; // Extents of the ion vx grid (m/s).
 double vy_min_ion, vy_max_ion; // Extents of the ion vy grid (m/s).

 int Nx;                  // Number of cells along x.
 int Nvx;                 // Number of cells along vx.
 int Nvy;                 // Number of cells along vy.
 int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
 int poly_order;          // Polynomial order of the basis.

 double t_end;             // Final simulation time (s).
 int num_frames;           // Number of output frames.
 int field_energy_calcs; // Number of times to calculate field energy.
 int integrated_mom_calcs; // Number of times to calculate integrated moments.
 int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
 double dt_failure_tol;    // Minimum allowable fraction of initial time-step.
 int num_failures_max;     // Maximum allowable number of consecutive small time-steps.
};

static inline double sq(double x) { return x * x; }

// These next few functions are for defining hyperbolic trig functions that will be used in the source term
static inline double sech2(double x)  // Hyperbolic secant squared
{
    double sech = 1.0/cosh(x);
    return sq(sech);
}

static inline double asech(double x) // Inverse hyperbolic secant
{
    return log((1 + sqrt(1-sq(x)))/x);
}

static inline double maxwellian2D(double n, double vx, double vy, double ux, double uy, double vt)
{
 double v2 = sq(vx - ux) + sq(vy - uy);
 return n/(2.0*M_PI*sq(vt)) * exp(-(v2) / (2 * sq(vt)));
}

// The following four functions return certain parameters needed to calculate the density and magnetic field strength for different radii based on the Bennett pinch relation
// Equations come from Allen and Simons (2018):  http://doi.org/10.1017/S0022377818001137

// This function calculates the electron drift velocity.
// It assumes all of the current is carried by the electrons and that this drift does not change with radius
// We know these are not true but they allow for a good first approximation
static inline double uBennett(double e, double n0, double I, double a)  // I is current in A, a is pinch radius in m
{
    return I/e/n0/M_PI/sq(a);
}

// This function calculates a characteristic length called b
static inline double bBennett(double mu0, double e, double u, double Te, double Ti)
{
    return mu0*sq(e)*sq(u)/8/(Te+Ti);
}

// This function calculates the density as a function of radius
static inline double nBennett(double n0, double b, double r)
{
    return n0/sq(1 + b*n0*sq(r));
}

// This function calculates the magnetic field strength as a function of radius
static inline double BmagBennett(double mu0, double e, double u, double b, double n0, double r)
{
    return mu0*e*u*n0*r/2.0/(1+b*n0*sq(r));
}

// Electron intialization function
void evalDistFuncElc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0], vx = xn[1], vy = xn[2];
 fout[0] = maxwellian2D(app->ne, vx, vy, app->ue_x, app->ue_y, app->vte);
}

// Ion intialization function.
void evalDistFuncIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0], vx = xn[1], vy = xn[2];
 fout[0] = maxwellian2D(app->ni, vx, vy, app->ui_x, app->ui_y, app->vti);
}

// Electron source function.
void evalDistFuncElcSource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0], vx = xn[1], vy = xn[2];
 double fv = maxwellian2D(1.0, vx, vy, app->ue_x, app->ue_y, app->vte);

 double perc_value = 0.01;  // This is the percentage that we define the source term to be at 800 Debye lengths from the wall
 double b_source_coeff = asech(sqrt(perc_value))/((app->x_max) - 800.0*(app->debye));
 double c_source_coeff = b_source_coeff/(2.0*tanh(b_source_coeff*(app->x_max)));

 fout[0] = c_source_coeff * sech2(b_source_coeff*x)*fv;
}

// Ion source function.
void evalDistFuncIonSource(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0], vx = xn[1], vy = xn[2];
 double fv = maxwellian2D(1.0, vx, vy, app->ui_x, app->ui_y, app->vti);

 double perc_value = 0.01;  // This is the percentage that we define the source term to be at 800 Debye lengths from the wall
 double b_source_coeff = asech(sqrt(perc_value))/((app->x_max) - 800.0*(app->debye));
 double c_source_coeff = b_source_coeff/(2.0*tanh(b_source_coeff*(app->x_max)));

 fout[0] = c_source_coeff * sech2(b_source_coeff*x)*fv;
}

// Collision functions.
void evalNuEE(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0];
 double nu_ee = app->nu_ee;
 fout[0] = nu_ee;
}

void evalNuII(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 struct sheath_ctx *app = ctx;
 double x = xn[0];
 double nu_ii = app->nu_ii;
 fout[0] = nu_ii;
}

// External fields.
void external_fields(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
 double x = xn[0];
 struct sheath_ctx *app = ctx;

 double Ex = 0.0;
 double Ey = 0.0;
 double Ez = 0.0;
 double Bx = 0.0;
 double By = 0.0;
 double Bz = app->B;

 fout[0] = Ex;
 fout[1] = Ey;
 fout[2] = Ez;
 fout[3] = Bx;
 fout[4] = By;
 fout[5] = Bz;
}

struct sheath_ctx
create_ctx(void)
{
 int cdim = 1, vdim = 2; // Dimensionality.

 double r = 0.003;     // Radius at which we want to do the simulation. Density and magnetic field strength are scaled based on Bennet pinch relations

 // Set background Zpinch parameters for the FuZE regime
 // Values from Zhang et al. (2019): https://doi.org/10.1103/PhysRevLett.122.135001
 double n0 = 1.1e23;   // Base density (m^-3)
 double T0 = 2000*GKYL_ELEMENTARY_CHARGE;  // Base temperature (J)
 double a = 0.003;    // Pinch radius (m)
 double I = 200.e3;   // Pinch current (A) 

 double ep0 = GKYL_EPSILON0;          // Permittivity of free space.
 double mu0 = GKYL_MU0;
 double qe = -GKYL_ELEMENTARY_CHARGE; // Electron charge (C).
 double qi = GKYL_ELEMENTARY_CHARGE;  // Ion charge (C).
 double me = GKYL_ELECTRON_MASS;      // Electron mass (kg).
 double mi = GKYL_PROTON_MASS;    // Ion mass (kg).

 double Te = T0; // Electron temperature (J).
 double Ti = T0; // Ion temperature (J).

 // Calculate the electron drift velocity assuming they carry all the current uniformly in r (we know this isn't true but makes for a good first approximation)
 double u = uBennett(GKYL_ELEMENTARY_CHARGE, n0, I, a);

 // Calculate characteristic Bennett pinch length
 double b = bBennett(mu0, GKYL_ELEMENTARY_CHARGE, u, Te, Ti);

 // Calculate densities based on Bennett pinch profiles
 double ne = nBennett(n0, b, r);  // Density (m^-3)
 double ni = nBennett(n0, b, r);  // Density (m^-3)

 // Calculate magnetic field strength based on Bennett pinch profiles
 double B = BmagBennett(mu0, GKYL_ELEMENTARY_CHARGE, u, b, n0, r);  // Magnetic field (T)

 double ue_x = 0; // electron flow speed in x
 double ue_y = 0; // electron flow speed in y
 double ui_x = 0; // ion flow speed in x
 double ui_y = 0; // ion flow speed in y

 double vte = sqrt(Te / me); // Electron thermal velocity (m/s).
 double vti = sqrt(Ti / mi); // Ion thermal velocity (m/s).

 double debye = sqrt((ep0 * Te) / (ne * sq(qe))); // Debye length (m).
 double wpe = sqrt((sq(qe) * ne) / (ep0 * me));   // Electron plasma frequency (rad/s).
 double wpi = sqrt((sq(qi) * ni) / (ep0 * mi));   // Ion plasma frequency (rad/s).
 double Oci = qi*B/mi;

 double Lx = 2048.0 * debye;  // Domain length (m).
 double Lsrc = 1.0; // Source length (m).

 double mfp = 50.0 * debye; // mean free path
 double nu_ee = vte / mfp;   // electron-electron collision frequency
 double nu_ii = vti / mfp;   // ion-ion collision frequency

 double phi_bias = 10.0e3;            // bias potential (V)

 double x_min = -Lx, x_max = Lx;                        // Extents of the x grid (m).
 double vx_min_elc = -4.0 * vte, vx_max_elc = 4.0 * vte; // Extents of the electron vx grid (m/s).
 double vy_min_elc = -4.0 * vte, vy_max_elc = 4.0 * vte; // Extents of the electron vy grid (m/s).
 double vx_min_ion = -6.0 * vti, vx_max_ion = 6.0 * vti; // Extents of the ion vx grid (m/s).
 double vy_min_ion = -4.0 * vti, vy_max_ion = 4.0 * vti; // Extents of the ion vy grid (m/s).

 int Nx = 8192;       // Number of cells along x.
 int Nvx = 16;       // Number of cells along vx.
 int Nvy = 8;        // Number of cells along vy.
 int poly_order = 2; // Polynomial order of the basis.

 double t_end = 5.0*2.0*M_PI/Oci;          // Final simulation time (s).
 int num_frames = 500;                // Number of output frames.
 int field_energy_calcs = INT_MAX;   // Number of times to compute field energy.
 int integrated_mom_calcs = INT_MAX; // Number of times to compute integrated diagnostics.
 int integrated_L2_f_calcs = INT_MAX; // Number of times to compute L2 norm of distribution function.
 double dt_failure_tol = 1.0e-4;     // Minimum allowable fraction of initial time-step.
 int num_failures_max = 20;          // Maximum allowable number of consecutive small time-steps.

 struct sheath_ctx ctx = {
     .cdim = cdim,
     .vdim = vdim,

     .ep0 = ep0,
     .mu0 = mu0,    

     .qe = qe,
     .qi = qi,
     .me = me,
     .mi = mi,

     .ne = ne,
     .ni = ni,

     .r = r,
     .I = I,
     .a = a,

     .ue_x = ue_x,
     .ue_y = ue_y,
     .ui_x = ui_x,
     .ui_y = ui_y,

     .vte = vte,
     .vti = vti,

     .debye = debye,

     .Lsrc = Lsrc,

     .nu_ee = nu_ee,
     .nu_ii = nu_ii,

     .phi_bias = phi_bias,
     .B = B,

     .x_min = x_min,
     .x_max = x_max,
     .vx_min_elc = vx_min_elc,
     .vx_max_elc = vx_max_elc,
     .vy_min_elc = vy_min_elc,
     .vy_max_elc = vy_max_elc,
     .vx_min_ion = vx_min_ion,
     .vx_max_ion = vx_max_ion,
     .vy_min_ion = vy_min_ion,
     .vy_max_ion = vy_max_ion,
     .Nx = Nx,
     .Nvx = Nvx,
     .Nvy = Nvy,
     .cells = {Nx, Nvx, Nvy},
     .poly_order = poly_order,

     .t_end = t_end,
     .num_frames = num_frames,
     .field_energy_calcs = field_energy_calcs,
     .integrated_mom_calcs = integrated_mom_calcs,
     .integrated_L2_f_calcs = integrated_L2_f_calcs,
     .dt_failure_tol = dt_failure_tol,
     .num_failures_max = num_failures_max,
 };
 return ctx;
}

void calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_field_energy(app, t_curr);
  }
}

void calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_mom(app, t_curr);
  }
}

void calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_vlasov_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr) || force_calc) {
    gkyl_vlasov_app_calc_integrated_L2_f(app, t_curr);
  }
}

void write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double t_curr, bool force_write)
{
 if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_write)
 {
  int frame = iot->curr - 1;
  if (force_write)
  {
   frame = iot->curr;
  }

  gkyl_vlasov_app_write(app, t_curr, frame);
  gkyl_vlasov_app_write_field_energy(app);
  gkyl_vlasov_app_write_integrated_mom(app);
  gkyl_vlasov_app_write_integrated_L2_f(app);

  gkyl_vlasov_app_calc_mom(app);
  gkyl_vlasov_app_write_mom(app, t_curr, frame);
 }
}

int main(int argc, char **argv)
{
 struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
 if (app_args.use_mpi)
 {
  MPI_Init(&argc, &argv);
 }
#endif

 if (app_args.trace_mem)
 {
  gkyl_cu_dev_mem_debug_set(true);
  gkyl_mem_debug_set(true);
 }

 struct sheath_ctx ctx = create_ctx(); // Context for initialization functions.

 int NX = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.Nx);
 int NVX = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.Nvx);
 int NVY = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.Nvy);

 int nrank = 512; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
 if (app_args.use_mpi)
 {
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);
 }
#endif

 int ccells[] = {NX};
 int cdim = sizeof(ccells) / sizeof(ccells[0]);

 int cuts[cdim];
#ifdef GKYL_HAVE_MPI
 for (int d = 0; d < cdim; d++)
 {
  if (app_args.use_mpi)
  {
   cuts[d] = app_args.cuts[d];
  }
  else
  {
   cuts[d] = 1;
  }
 }
#else
 for (int d = 0; d < cdim; d++)
 {
  cuts[d] = 1;
 }
#endif

 // Construct communicator for use in app.
 struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
 if (app_args.use_gpu && app_args.use_mpi)
 {
#ifdef GKYL_HAVE_NCCL
  comm = gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){
      .mpi_comm = MPI_COMM_WORLD,
  });
#else
  printf(" Using -g and -M together requires NCCL.\n");
  assert(0 == 1);
#endif
 }
 else if (app_args.use_mpi)
 {
  comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){
      .mpi_comm = MPI_COMM_WORLD,
  });
 }
 else
 {
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
      .use_gpu = app_args.use_gpu});
 }
#else
 comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){
     .use_gpu = app_args.use_gpu});
#endif

 int my_rank;
 gkyl_comm_get_rank(comm, &my_rank);
 int comm_size;
 gkyl_comm_get_size(comm, &comm_size);

 int ncuts = 1;
 for (int d = 0; d < cdim; d++)
 {
  ncuts *= cuts[d];
 }

 if (ncuts != comm_size)
 {
  if (my_rank == 0)
  {
   fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
  }
  goto mpifinalize;
 }

 // Electrons.
 struct gkyl_vlasov_species elc = {
     .name = "elc",
     .charge = ctx.qe,
     .mass = ctx.me,
     .lower = {ctx.vx_min_elc, ctx.vy_min_elc},
     .upper = {ctx.vx_max_elc, ctx.vy_max_elc},
     .cells = {NVX, NVY},

     // Initial conditions.
     .num_init = 1,
     .projection[0] = {
         .proj_id = GKYL_PROJ_FUNC,
         .func = evalDistFuncElc,
         .ctx_func = &ctx,
     },

     // Source.
     .source = {
         .source_id = GKYL_BFLUX_SOURCE,
         .source_length = ctx.Lsrc,
         .source_species = "ion",
         .num_sources = 1,
         .projection[0] = {
             .proj_id = GKYL_PROJ_FUNC,
             .func = evalDistFuncElcSource,
             .ctx_func = &ctx,
         },
     },

     // Collisions.
     .collisions = {
         .collision_id = GKYL_LBO_COLLISIONS,
         .self_nu = evalNuEE,
         .ctx = &ctx,
         .num_cross_collisions = 1,
         .collide_with = {"ion"},
     },

     .bcx = {
         .lower = {
             .type = GKYL_SPECIES_ABSORB,
         },
         .upper = {
             .type = GKYL_SPECIES_ABSORB,
         },
     },

     .num_diag_moments = 4,
     .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2IJ, GKYL_F_MOMENT_M3, GKYL_F_MOMENT_LTE},
 };

 // Ions.
 struct gkyl_vlasov_species ion = {
     .name = "ion",
     .charge = ctx.qi,
     .mass = ctx.mi,
     .lower = {ctx.vx_min_ion, ctx.vy_min_ion},
     .upper = {ctx.vx_max_ion, ctx.vy_max_ion},
     .cells = {NVX, NVY},

     .num_init = 1,
     .projection[0] = {
         .proj_id = GKYL_PROJ_FUNC,
         .func = evalDistFuncIon,
         .ctx_func = &ctx,
     },

     .source = {
         .source_id = GKYL_BFLUX_SOURCE,
         .source_length = ctx.Lsrc,
         .source_species = "ion",
         .num_sources = 1,
         .projection[0] = {
             .proj_id = GKYL_PROJ_FUNC,
             .func = evalDistFuncIonSource,
             .ctx_func = &ctx,
         },
     },

     .collisions = {
         .collision_id = GKYL_LBO_COLLISIONS,
         .self_nu = evalNuII,
         .ctx = &ctx,
         .num_cross_collisions = 1,
         .collide_with = {"elc"},
     },

     .bcx = {
         .lower = {
             .type = GKYL_SPECIES_ABSORB,
         },
         .upper = {
             .type = GKYL_SPECIES_ABSORB,
         },
     },

     .num_diag_moments = 4,
     .diag_moments = {GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2IJ, GKYL_F_MOMENT_M3, GKYL_F_MOMENT_LTE},
 };

 // Field.
 struct gkyl_vlasov_field field = {
     .epsilon0 = ctx.ep0,
     .poisson_bcs = {
         .lo_type = {GKYL_POISSON_DIRICHLET},
         .up_type = {GKYL_POISSON_DIRICHLET},
         .lo_value = {ctx.phi_bias},
         .up_value = {0.0}},
     .ext_em = external_fields,
     .ext_em_ctx = &ctx,
     .ext_em_evolve = false,
 };

 // VP app
 struct gkyl_vm app_inp = {
     .name = "vp_magSheath_phiB_10kV_1x2v",

     .cdim = ctx.cdim,
     .vdim = ctx.vdim,
     .lower = {ctx.x_min},
     .upper = {ctx.x_max},
     .cells = {NX},
     .poly_order = ctx.poly_order,
     .basis_type = app_args.basis_type,

     .num_periodic_dir = 0,
     .periodic_dirs = {},

     .num_species = 2,
     .species = {elc, ion},

     .field = field,
     .is_electrostatic = true,

     .parallelism = {
         .use_gpu = app_args.use_gpu,
         .cuts = {app_args.cuts[0]},
         .comm = comm,
     }};

 // Create app object.
 gkyl_vlasov_app *app = gkyl_vlasov_app_new(&app_inp);

 // Initial and final simulation times.
 double t_curr = 0.0, t_end = ctx.t_end;

 // Initialize simulation.
 int frame_curr = 0;
 if (app_args.is_restart)
 {
  struct gkyl_app_restart_status status = gkyl_vlasov_app_read_from_frame(app, app_args.restart_frame);

  if (status.io_status != GKYL_ARRAY_RIO_SUCCESS)
  {
   gkyl_vlasov_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
   goto freeresources;
  }

  frame_curr = status.frame;
  t_curr = status.stime;

  gkyl_vlasov_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
  gkyl_vlasov_app_cout(app, stdout, " at time = %g\n", t_curr);
 }
 else
 {
  gkyl_vlasov_app_apply_ic(app, t_curr);
 }

 // Create trigger for field energy.
 int field_energy_calcs = ctx.field_energy_calcs;
 struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

 calc_field_energy(&fe_trig, app, t_curr, false);

 // Create trigger for integrated moments.
 int integrated_mom_calcs = ctx.integrated_mom_calcs;
 struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

 calc_integrated_mom(&im_trig, app, t_curr, false);

 // Create trigger for integrated L2 norm of the distribution function.
 int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
 struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

 calc_integrated_L2_f(&l2f_trig, app, t_curr, false);

 // Create trigger for IO.
 int num_frames = ctx.num_frames;
 struct gkyl_tm_trigger io_trig = {.dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr};

 write_data(&io_trig, app, t_curr, false);

 // Compute initial guess of maximum stable time-step.
 double dt = t_end - t_curr;

 // Initialize small time-step check.
 double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
 int num_failures = 0, num_failures_max = ctx.num_failures_max;

 long step = 1;
 while ((t_curr < t_end) && (step <= app_args.num_steps))
 {
  gkyl_vlasov_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
  struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
  gkyl_vlasov_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

  if (!status.success)
  {
   gkyl_vlasov_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
   break;
  }

  t_curr += status.dt_actual;
  dt = status.dt_suggested;

  calc_field_energy(&fe_trig, app, t_curr, false);
  calc_integrated_mom(&im_trig, app, t_curr, false);
  calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
  write_data(&io_trig, app, t_curr, false);

  if (dt_init < 0.0)
  {
   dt_init = status.dt_actual;
  }
  else if (status.dt_actual < dt_failure_tol * dt_init)
  {
   num_failures += 1;

   gkyl_vlasov_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
   gkyl_vlasov_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
   gkyl_vlasov_app_cout(app, stdout, " num_failures = %d\n", num_failures);
   if (num_failures >= num_failures_max)
   {
    gkyl_vlasov_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
    gkyl_vlasov_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);

    calc_field_energy(&fe_trig, app, t_curr, true);
    calc_integrated_mom(&im_trig, app, t_curr, true);
    calc_integrated_L2_f(&l2f_trig, app, t_curr, true);
    write_data(&io_trig, app, t_curr, true);

    break;
   }
  }
  else
  {
   num_failures = 0;
  }

  step += 1;
 }

 calc_field_energy(&fe_trig, app, t_curr, false);
 calc_integrated_mom(&im_trig, app, t_curr, false);
 calc_integrated_L2_f(&l2f_trig, app, t_curr, false);
 write_data(&io_trig, app, t_curr, false);
 gkyl_vlasov_app_stat_write(app);

 struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

 gkyl_vlasov_app_cout(app, stdout, "\n");
 gkyl_vlasov_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
 gkyl_vlasov_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
 gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
 if (stat.nstage_2_fail > 0)
 {
  gkyl_vlasov_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
  gkyl_vlasov_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
 }
 gkyl_vlasov_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
 gkyl_vlasov_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
 gkyl_vlasov_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
 gkyl_vlasov_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
 gkyl_vlasov_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
 gkyl_vlasov_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

 gkyl_vlasov_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
 gkyl_vlasov_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
 // Free resources after simulation completion.
 gkyl_comm_release(comm);
 gkyl_vlasov_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
 if (app_args.use_mpi)
 {
  MPI_Finalize();
 }
#endif

 return 0;
}
