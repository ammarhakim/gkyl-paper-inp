#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>

struct jets_ctx {
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te_Ti; // electron to ion temperature ratio
  double vte; // electron thermal velocity
  double vti; // ion thermal velocity
  double cs; // sound speed
  double uShock; // in-flow velocity
  double Lx; // size of the box
  double n0; // initial number density
};

static inline double sq(double x) { return x*x; }

void
evalDistFuncElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vte, vdrift = app->uShock, n0 = app->n0;
  double fv = 0.0;
  if (x < 0)
    fv = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt))));
  else
    fv = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalDistFuncIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  double x = xn[0], v = xn[1];
  double vt = app->vti, vdrift = app->uShock, n0 = app->n0;
  double fv = 0.0;
  if (x < 0)
    fv = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v-vdrift)/(2*sq(vt))));
  else
    fv = n0/sqrt(2.0*M_PI*sq(vt))*(exp(-sq(v+vdrift)/(2*sq(vt))));
  fout[0] = fv;
}

void
evalTperpElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vte;
  double n = app->n0;
  fout[0] = n*vt*vt;
}

void
evalTperpIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  double x = xn[0];
  double vt = app->vti;
  double n = app->n0;
  fout[0] = n*vt*vt;
}

void
evalFieldFunc(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  double x = xn[0];
  
  fout[0] = 0.0; fout[1] = 0.0, fout[2] = 0.0;
  fout[3] = 0.0; fout[4] = 0.0; fout[5] = 0.0;
  fout[6] = 0.0; fout[7] = 0.0;
}

void
evalNuElc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  fout[0] = 1.038885e+9;
}

void
evalNuIon(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct jets_ctx *app = ctx;
  fout[0] = 3.516241e+06;
}

struct jets_ctx
create_ctx(void)
{
  struct jets_ctx ctx = {
    .n0 = 1e20,
    .chargeElc = -1.602177e-19,
    .massElc = 9.109384e-31,
    .chargeIon = 1.602177e-19,
    .massIon = 6.633521e-26,
    .Te_Ti = 1.0,
    .vte = 5.136370e+05,
    .vti = 1.903394e+03,
    .uShock = 15e3,
    .Lx = 2e-2
  };
  return ctx;
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct jets_ctx ctx = create_ctx(); // context for init functions

  // electron Tperp  
  struct gkyl_vlasov_fluid_species elc_Tperp = {
    .name = "elc_Tperp",

    .ctx = &ctx,
    .init = evalTperpElc,

    .advection = {.advect_with = "elc", .collision_id = GKYL_LBO_COLLISIONS},
  };  

  // electrons
  struct gkyl_vlasov_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -6.0 * ctx.vte},
    .upper = { 6.0 * ctx.vte}, 
    .cells = { 64 },

    .ctx = &ctx,
    .init = evalDistFuncElc,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuElc,
      .collide_with_fluid = "elc_Tperp",
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },    

    .num_diag_moments = 4,
    .diag_moments = { "M0", "M1i", "M2", "M3i" },
  };

  // ion Tperp                                                                                              
  struct gkyl_vlasov_fluid_species ion_Tperp = {
    .name = "ion_Tperp",

    .ctx = &ctx,
    .init = evalTperpIon,

    .advection = {.advect_with = "ion", .collision_id = GKYL_LBO_COLLISIONS},
  };  

  // ions
  struct gkyl_vlasov_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -3 * ctx.uShock},
    .upper = { 3 * ctx.uShock}, 
    .cells = { 64 },

    .ctx = &ctx,
    .init = evalDistFuncIon,

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,

      .ctx = &ctx,
      .self_nu = evalNuIon,
      .collide_with_fluid = "ion_Tperp" ,
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },    

    .num_diag_moments = 4,
    .diag_moments = { "M0", "M1i", "M2", "M3i" },
  };

  // field
  struct gkyl_vlasov_field field = {
    .epsilon0 = 8.8541878128e-12, .mu0 = 0.004280934795130751,
    .elcErrorSpeedFactor = 0.0,
    .mgnErrorSpeedFactor = 0.0,

    .ctx = &ctx,
    .init = evalFieldFunc
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "jets",

    .cdim = 1, .vdim = 1,
    .lower = { -ctx.Lx },
    .upper = { ctx.Lx },
    .cells = { 768 },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = { },

    .num_species = 2,
    .species = { elc, ion },
    .num_fluid_species = 2,
    .fluid_species = { elc_Tperp, ion_Tperp },
    .field = field,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(&vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = 2e-6;
  double dt = tend-tcurr;
  int nframes = 40, frame_count = 1;
  double frame_inc = tend / nframes;

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  
  gkyl_vlasov_app_write(app, tcurr, 0);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 0);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    if (tcurr > frame_count * frame_inc) {
      gkyl_vlasov_app_write(app, tcurr, frame_count);
      gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, frame_count);
      frame_count += 1;
    }

    tcurr += status.dt_actual;
    dt = status.dt_suggested;
    step += 1;
  }

  gkyl_vlasov_app_write(app, tcurr, nframes);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, nframes);
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
  printf("Current evaluation and accumulate took %g secs\n", stat.current_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
