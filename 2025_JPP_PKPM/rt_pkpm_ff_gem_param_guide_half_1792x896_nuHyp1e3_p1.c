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

struct pkpm_ff_ctx {
  // Fundamental constants
  double epsilon0; // permittivity of free space
  double mu0; // permeability of free space
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  // Plasma parameters
  double T_e; // electron temperature
  double T_i; // ion temperature
  double vtElc; // electron thermal velocity sqrt(2 T_e/massElc)
  double vtIon; // ion thermal velocity sqrt(2 T_i/massIon)
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double n0; // reference density
  double B0; // reference magnetic field strength
  // Reconnection parameters
  double guide; // guide field strength
  double w0; // layer width
  double psi0; // size of perturbation
  double noise_amp; // amplitude of noise perturbation
  double noise_index; // spectral index of noise perturbation
  int k_init; // initial wavenumber to perturb
  int k_final; // final wavenumber to perturb
  // Simulation parameters
  double Lx; // Box size in x
  double Ly; // Box size in y
  int poly_order; // Polynomial order.
  double cfl_frac; // CFL coefficient.

  double t_end; // Final simulation time.
  int num_frames; // Number of output frames.
  int field_energy_calcs; // Number of times to calculate field energy.
  int integrated_mom_calcs; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

static inline double
maxwellian(double n, double v, double temp, double mass)
{
  double v2 = v*v;
  return n/sqrt(2.0*M_PI*temp/mass)*exp(-v2/(2.0*temp/mass));
}

static inline double
sech2(double x)
{
  return 1.0/(cosh(x)*cosh(x));
}

static inline void
noise_init(double noise_amp, double noise_index, int k_init, int k_final, double Lx, double Ly, double x, double y, double noise[3])
{
  pcg64_random_t rng = gkyl_pcg64_init(0);
  double kindex = (noise_index + 1.0) / 2.0;
  double B_amp = 0.0; 
  double B_phase = 0.0;
  for (int i = k_init; i < k_final; ++i) {
    B_amp = gkyl_pcg64_rand_double(&rng);
    B_phase = gkyl_pcg64_rand_double(&rng);

    noise[0] -= 2.0*(2.0*M_PI/Ly)*(Lx/(i*2.0*M_PI))*B_amp*sin(2.0*M_PI*y/Ly)*(cos(2.0*M_PI*y/Ly)+1)*cos(i*2.0*M_PI*x/Lx +  2.0*M_PI*B_phase)*pow(i,kindex);
    noise[1] += B_amp*(cos(2.0*M_PI*y/Ly) + 1.0)*(cos(2.0*M_PI*y/Ly) + 1.0)*sin(i*2.0*M_PI*x/Lx + 2.0*M_PI*B_phase)*pow(i,kindex);
    noise[2] += (2.0*M_PI*i/Lx)*B_amp*(cos(2.0*M_PI*y/Ly) + 1.0)*(cos(2.0*M_PI*y/Ly) + 1.0)*cos(i*2.0*M_PI*x/Lx + 2*M_PI*B_phase)*pow(i,kindex) + 
                 2.0*(2.0*M_PI/Ly)*(2.0*M_PI/Ly)*(Lx/(i*2.0*M_PI))*B_amp*(sin(2.0*M_PI*y/Ly)*sin(2.0*M_PI*y/Ly) - cos(2.0*M_PI*y/Ly)*(cos(2.0*M_PI*y/Ly)+1.0))*cos(i*2.0*M_PI*x/Lx +  2.*M_PI*B_phase)*pow(i,kindex);
  }
  double kdiff = floor(k_final) - floor(k_init) + 1.0;
  noise[0] = noise_amp*noise[0]/sqrt(2.0*kdiff*kdiff/3.0);
  noise[1] = noise_amp*noise[1]/sqrt(2.0*kdiff*kdiff/3.0);
  noise[2] = noise_amp*noise[2]/sqrt(2.0*kdiff*kdiff/3.0);
}

void
evalDistFuncElc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double me = app->massElc;
  double T_e = app->T_e;
  double n0 = app->n0;
  
  double fv = maxwellian(n0, vx, T_e, me);
    
  fout[0] = fv;
  fout[1] = T_e/me*fv;
}
void
evalDistFuncIon(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  
  double x = xn[0], y = xn[1], vx = xn[2];

  double mi = app->massIon;
  double T_i = app->T_i;
  double n0 = app->n0;
  
  double fv = maxwellian(n0, vx, T_i, mi);
    
  fout[0] = fv;
  fout[1] = T_i/mi*fv;
}

void
evalFluidElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double guide = app->guide;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double noise_amp = app->noise_amp;
  double noise_index = app->noise_index;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double noise[3] = {0.0};
  noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y, noise);

  double Jx = -B0/w0*tanh(y/w0)*sech2(y/w0)/(sqrt(sech2(y/w0) + guide*guide));
  double Jy = 0.0;
  double Jz = B0/w0*sech2(y/w0) - psi0*cos(pi_2*x/Lx)*cos(M_PI*y/Ly)*((pi_2/Lx)*(pi_2/Lx) + (M_PI/Ly)*(M_PI/Ly)) + noise[2];

  fout[0] = me*Jx/qe;
  fout[1] = 0.0;
  fout[2] = me*Jz/qe;
}

void
evalFluidIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  
  double x = xn[0], y = xn[1];

  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;

  double x = xn[0], y = xn[1];

  double qe = app->chargeElc;
  double qi = app->chargeIon;
  double me = app->massElc;
  double mi = app->massIon;
  double Lx = app->Lx;
  double Ly = app->Ly;
  double B0 = app->B0;
  double guide = app->guide;
  double w0 = app->w0;
  double psi0 = app->psi0;  
  double noise_amp = app->noise_amp;
  double noise_index = app->noise_index;
  int k_init = app->k_init;
  int k_final = app->k_final;

  double pi_2 = 2.0*M_PI;
  double pi_4 = 4.0*M_PI;

  double noise[3] = {0.0};
  noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y, noise);

  double b1x = -B0*tanh(y/w0);
  double b1y = 0.0;
  double b1z = B0*sqrt(guide*guide + sech2(y/w0));

  double E_x = 0.0;
  double E_y = 0.0;
  double E_z = 0.0;
  double B_x = b1x + psi0 * (M_PI / Ly) * cos(2 * M_PI * x / Lx) * sin(M_PI * y / Ly) + noise[0];
  double B_y = b1y - psi0 * (2 * M_PI / Lx) * sin(2 * M_PI * x / Lx) * cos(M_PI * y / Ly) + noise[1];
  double B_z = b1z;
  
  fout[0] = E_x; fout[1] = E_y, fout[2] = E_z;
  fout[3] = B_x; fout[4] = B_y; fout[5] = B_z;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  fout[0] = app->nuElc;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct pkpm_ff_ctx *app = ctx;
  fout[0] = app->nuIon;
}

struct pkpm_ff_ctx
create_ctx(void)
{
  double epsilon0 = 1.0; // permittivity of free space
  double mu0 = 1.0; // permiability of free space
  double lightSpeed = 1.0/sqrt(epsilon0*mu0);

  double massElc = 1.0; // electron mass
  double chargeElc = -1.0; // electron charge
  double massIon = 1836.153; // ion mass
  double chargeIon = 1.0; // ion charge

  double Te_Ti = 0.2; // ratio of electron to ion temperature
  double n0 = 1.0; // initial number density

  // how non-relativistic the simulation is vtElc/c = 1/8 -> 4 kev electrons
  double vtElc = 1.0/16.0;
  double T_e = vtElc*vtElc*massElc/2.0;
  double T_i = T_e/Te_Ti;
  double vtIon = sqrt(2.0*T_i/massIon);

  // plasma beta and B0 based on *in-plane* field
  double guide = 0.5;
  double beta_elc = 1.0/6.0;
  double B0 = sqrt((2.0*mu0*n0*T_e)/beta_elc);

  // frequencies and skin depths
  double omegaCi = chargeIon*B0/massIon;
  double wpe = sqrt(n0*chargeElc*chargeElc/(epsilon0*massElc));
  double de = lightSpeed/wpe;
  double wpi = sqrt(n0*chargeIon*chargeIon/(epsilon0*massIon));
  double di = lightSpeed/wpi;

  // Layer width and perturbation
  double w0 = 0.5*di;
  double psi0 = 0.1*B0*di;

  // noise levels for perturbation
  double noise_amp = 0.01*B0;
  int k_init = 1;            // first wave mode to perturb with noise, 1.0 correspond to box size
  int k_final = 20;          // last wave mode to perturb with noise
  double noise_index = -1.0; // spectral index of the noise

  // collision frequencies
  double nuElc = 0.01*omegaCi;
  double nuIon = 0.01*omegaCi/sqrt(massIon);

  // domain size and simulation time
  double Lx = 8.0*M_PI*di;
  double Ly = 4.0*M_PI*di;
  int poly_order = 1; // Polynomial order.
  double cfl_frac = 1.0; // CFL coefficient.

  double t_end = 50.0/omegaCi; // Final simulation time.
  int num_frames = 100; // Number of output frames.
  int field_energy_calcs = 5000; // Number of times to calculate field energy.
  int integrated_mom_calcs = 5000; // Number of times to calculate integrated moments.
  int integrated_L2_f_calcs = 5000; // Number of times to calculate integrated L2 norm of distribution function.
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  
  struct pkpm_ff_ctx ctx = {
    .epsilon0 = epsilon0,
    .mu0 = mu0,
    .chargeElc = chargeElc,
    .massElc = massElc,
    .chargeIon = chargeIon,
    .massIon = massIon,
    .T_e = T_e,
    .T_i = T_i,
    .vtElc = vtElc,
    .vtIon = vtIon,
    .nuElc = nuElc,
    .nuIon = nuIon,
    .n0 = n0,
    .B0 = B0,
    .guide = guide,
    .w0 = w0,
    .psi0 = psi0,
    .noise_amp = noise_amp,
    .noise_index = noise_index,
    .k_init = k_init,
    .k_final = k_final,
    .Lx = Lx,
    .Ly = Ly,
    .poly_order = poly_order,
    .cfl_frac = cfl_frac,
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

void
write_data(struct gkyl_tm_trigger* iot, gkyl_pkpm_app* app, double t_curr, bool force_write)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    if (force_write) {
      frame = iot->curr;
    }

    gkyl_pkpm_app_write(app, t_curr, frame);
    gkyl_pkpm_app_write_field_energy(app);
    gkyl_pkpm_app_write_integrated_mom(app);
    gkyl_pkpm_app_write_integrated_L2_f(app);
  }
}

void
calc_field_energy(struct gkyl_tm_trigger* fet, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(fet, t_curr)) {
    gkyl_pkpm_app_calc_field_energy(app, t_curr);
  }
}

void
calc_integrated_mom(struct gkyl_tm_trigger* imt, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(imt, t_curr)) {
    gkyl_pkpm_app_calc_integrated_mom(app, t_curr);
  }
}

void
calc_integrated_L2_f(struct gkyl_tm_trigger* l2t, gkyl_pkpm_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(l2t, t_curr)) {
    gkyl_pkpm_app_calc_integrated_L2_f(app, t_curr);
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

  int NX = APP_ARGS_CHOOSE(app_args.xcells[0], 1792);
  int NY = APP_ARGS_CHOOSE(app_args.xcells[1], 896);
  int VX = APP_ARGS_CHOOSE(app_args.vcells[0], 32);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
     
  struct pkpm_ff_ctx ctx = create_ctx(); // context for init functions

  int nrank = 1; // Number of processors in simulation.
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
#endif  

  int ccells[] = { NX, NY };
  int cdim = sizeof(ccells) / sizeof(ccells[0]);

  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d = 0; d < cdim; d++) {
    if (app_args.use_mpi) {
      cuts[d] = app_args.cuts[d];
    }
    else {
      cuts[d] = 1;
    }
  }
#else
  for (int d = 0; d < cdim; d++) {
    cuts[d] = 1;
  }
#endif  
    
  // Construct communicator for use in app.
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    printf(" Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);
  int comm_size;
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) {
    ncuts *= cuts[d];
  }

  if (ncuts != comm_size) {
    if (my_rank == 0) {
      fprintf(stderr, "*** Number of ranks, %d, does not match total cuts, %d!\n", comm_size, ncuts);
    }
    goto mpifinalize;
  }

  // electrons
  struct gkyl_pkpm_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -8.0 * ctx.vtElc},
    .upper = { 8.0 * ctx.vtElc}, 
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

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .diffusion = {.D = 1.0e-3, .order=4},
  };
  
  // ions
  struct gkyl_pkpm_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -8.0 * ctx.vtIon},
    .upper = { 8.0 * ctx.vtIon}, 
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

    .bcy = { GKYL_SPECIES_REFLECT, GKYL_SPECIES_REFLECT },
    .diffusion = {.D = 1.0e-3, .order=4},
  };

  // field
  struct gkyl_pkpm_field field = {
    .epsilon0 = ctx.epsilon0, .mu0 = ctx.mu0,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc,
    .bcy = { GKYL_FIELD_PEC_WALL, GKYL_FIELD_PEC_WALL },
  };

  // pkpm app
  struct gkyl_pkpm pkpm = {
    .name = "pkpm_ff_gem_param_1792x896_nuHyp1e3_p1",

    .cdim = 2, .vdim = 1,
    .lower = { -ctx.Lx/2.0, -ctx.Ly/2.0 },
    .upper = { ctx.Lx/2.0, ctx.Ly/2.0 },
    .cells = { NX, NY },

    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .cfl_frac = ctx.cfl_frac,

    .use_explicit_source = true, 

    .num_periodic_dir = 1,
    .periodic_dirs = { 0 },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // create app object
  gkyl_pkpm_app *app = gkyl_pkpm_app_new(&pkpm);

  // Initial and final simulation times.
  double t_curr = 0.0, t_end = ctx.t_end;

  // Initialize simulation.
  int frame_curr = 0;
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_pkpm_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_pkpm_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_pkpm_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_pkpm_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_pkpm_app_apply_ic(app, t_curr);
  }

  // Create trigger for field energy.
  int field_energy_calcs = ctx.field_energy_calcs;
  struct gkyl_tm_trigger fe_trig = { .dt = t_end / field_energy_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_field_energy(&fe_trig, app, t_curr);

  // Create trigger for integrated moments.
  int integrated_mom_calcs = ctx.integrated_mom_calcs;
  struct gkyl_tm_trigger im_trig = { .dt = t_end / integrated_mom_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_mom(&im_trig, app, t_curr);

  // Create trigger for integrated L2 norm of the distribution function.
  int integrated_L2_f_calcs = ctx.integrated_L2_f_calcs;
  struct gkyl_tm_trigger l2f_trig = { .dt = t_end / integrated_L2_f_calcs, .tcurr = t_curr, .curr = frame_curr };

  calc_integrated_L2_f(&l2f_trig, app, t_curr);

  // Create trigger for IO.
  int num_frames = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = { .dt = t_end / num_frames, .tcurr = t_curr, .curr = frame_curr };

  write_data(&io_trig, app, t_curr, false);

  // Compute initial guess of maximum stable time-step.
  double dt = (ctx.Lx/NX)/6.0; 

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_pkpm_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_pkpm_update(app, dt);
    gkyl_pkpm_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      gkyl_pkpm_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_field_energy(&fe_trig, app, t_curr);
    calc_integrated_mom(&im_trig, app, t_curr);
    calc_integrated_L2_f(&l2f_trig, app, t_curr);
    write_data(&io_trig, app, t_curr, false);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_pkpm_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_pkpm_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_pkpm_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_pkpm_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_pkpm_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  calc_field_energy(&fe_trig, app, t_curr);
  calc_integrated_mom(&im_trig, app, t_curr);
  calc_integrated_L2_f(&l2f_trig, app, t_curr);
  write_data(&io_trig, app, t_curr, false);
  gkyl_pkpm_app_stat_write(app);

  struct gkyl_pkpm_stat stat = gkyl_pkpm_app_stat(app);

  gkyl_pkpm_app_cout(app, stdout, "\n");
  gkyl_pkpm_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_pkpm_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_pkpm_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_pkpm_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_pkpm_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_pkpm_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_pkpm_app_cout(app, stdout, "Fluid species RHD calc took %g secs\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species PKPM vars took %g secs\n", stat.species_pkpm_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, stdout, "EM variables (bvar) calc took %g secs\n", stat.field_em_vars_tm);
  gkyl_pkpm_app_cout(app, stdout, "Current evaluation and accumulate took %g secs\n", stat.current_tm);
  gkyl_pkpm_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_pkpm_app_cout(app, stdout, "Number of write calls %ld\n", stat.nio);
  gkyl_pkpm_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_pkpm_app_release(app);

mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif

  return 0;
}
