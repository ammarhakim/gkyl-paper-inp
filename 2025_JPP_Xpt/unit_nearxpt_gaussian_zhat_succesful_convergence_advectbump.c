#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_run.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

struct gkyl_gk_block_geom*
create_gk_block_geom(void)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 6);

  struct gkyl_efit_inp efit_inp = {
    .filepath = "./cubic_null.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = false,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  //printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g, psisep = %g\n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry, efit->psisep);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "cubic_null_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "cubic_null_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "cubic_null_q.gkyl");
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  // psisep = 1.5093065418975686; //This is the value
  double wfac = 6;
  double wsol = 0.05*wfac;
  double wcore = 0.05*wfac;
  double wpf = 0.05*wfac;

  double psi_lo_sol = psisep;
  double psi_up_sol = psisep + wsol;

  double psi_lo_core = psisep - wcore;
  double psi_up_core = psisep;

  double psi_lo_pf = psisep - wpf;
  double psi_up_pf = psisep;


  double rcell_fac = 30.0;
  double pcell_fac = 30.0;

  int npsi_sol = 4*rcell_fac;
  int npsi_core = 4*rcell_fac;
  int npsi_pf = 4*rcell_fac;

  double ntheta_lower = 2*pcell_fac;
  double ntheta_middle = 8*pcell_fac;

  //int ncuts_lower = 1;//ntheta_lower/8;
  //int ncuts_middle = 1;//ntheta_middle/8;
  //enum gkyl_geometry_id id = GKYL_TOKAMAK;
  int ncuts_lower = 1;
  int ncuts_middle = 1;
  enum gkyl_geometry_id id = GKYL_GEOMETRY_FROMFILE;

  double rmax = 3.5;
  double rmin = 0.0;
  double zmax = 1.5;
  double zmin = -3.0;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;


  double compression_factor = 0.0;
  double radial_compression_factor = 0.0;

  // block 0. Lower outer PF region.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower},
      .cuts = { 1, ncuts_lower },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_R,
          .rright = rmax,
          .rleft = rmin,
          .rmin = rmin,
          .rmax = rmax,
          .zmin_right = zmin,
          .zmin_left = zmin,
          .plate_spec = false,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // Block 1: lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo },
      .upper = { psi_up_sol, theta_up },
      .cells = { npsi_sol, ntheta_lower },
      .cuts = { 1, ncuts_lower },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_LO,
          .rleft = rmin,
          .rright = rmax,
          .rmin = rmin,
          .rmax = rmax,
          .zmax = zmax,
          .zmin_right = zmin,
          .zmin_left = zmin,
          .plate_spec = false,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 2: mid SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo },
      .upper = { psi_up_sol, theta_up },
      .cells = { npsi_sol, ntheta_middle },
      .cuts = { 1, ncuts_middle },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_MID,
          .rleft = rmin,
          .rright = rmax,
          .rmin = rmin,
          .rmax = rmax,
          .zmin_right = zmin,
          .zmin_left = zmin,
          .zmax = zmax,
          .plate_spec = false,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 3: lower inner SOL.
  gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
      .lower = { psisep, theta_lo },
      .upper = { psi_up_sol, theta_up },
      .cells = { npsi_sol, ntheta_lower },
      .cuts = { 1, ncuts_lower },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_UP,
          .rleft = rmin,
          .rright = rmax,
          .rmin = rmin,
          .rmax = rmax,
          .zmin_right = zmin,
          .zmin_left = zmin,
          .zmax = zmax,
          .plate_spec = false,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
      },
      .connections[1] = { // z-direction.
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
      }
    }
  );




  // block 4. Lower inner PF region.
  gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower },
      .cuts = { 1, ncuts_lower },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_L,
          .rright = rmax,
          .rleft = rmin,
          .rmin = rmin,
          .rmax = rmax,
          .zmin_right = zmin,
          .zmin_left = zmin,
          .plate_spec = false,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );

  // Block 5: core region.
  gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_core, theta_lo },
      .upper = { psisep, theta_up },
      .cells = { npsi_core, ntheta_middle},
      .cuts = { 1, ncuts_middle },
      .geometry = {
        .world = {0.0},
        .geometry_id = id,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE,
          .rleft = rmin,
          .rright = rmax,
          .rmin = rmin,
          .rmax = rmax,
          .zmax = zmax,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
      },
      .connections[1] = { // z-direction.
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );
 
  return bgeom;
}

struct gk_step_ctx {
  int cdim, vdim; // Dimensionality.
  double chargeElc; // electron charge
  double massElc; // electron mass
  double chargeIon; // ion charge
  double massIon; // ion mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double vtIon;
  double vtElc;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nu_frac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  // Source parameters
  double nsource;
  double Tsource;
  double temp_fac;
  double n_fac;
  double cx;
  double cz;
  double xcenter;
  // Simulation parameters
  int Nx; // Cell count (configuration space: x-direction).
  int Nz; // Cell count (configuration space: z-direction).
  int Nvpar; // Cell count (velocity space: parallel velocity direction).
  int Nmu; // Cell count (velocity space: magnetic moment direction).
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  double vpar_max_elc; // Velocity space extents in vparallel for electrons
  double mu_max_elc; // Velocity space extents in mu for electrons
  double vpar_max_ion; // Velocity space extents in vparallel for ions
  double mu_max_ion; // Velocity space extents in mu for ions
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
  double write_phase_freq;
};


struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge


  double n_fac = 32.0;
  double temp_fac = 9.0/n_fac;

  double Te = 1500.0/temp_fac*eV;
  double Ti = 1500.0/temp_fac*eV;
  double B0 = 1.0; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);


  double nsource = 1.675e22*2.5*temp_fac/6.0;
  double Tsource = 3000.0*eV/temp_fac;
  double cx = 0.0065612;
  double cz = 0.4916200*1.4;
  double psisep = 1.5093065418975686; //This is the value from efit
  double xcenter = psisep;

  // Collision parameters.
  double nu_frac = 1.0;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nu_frac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nu_frac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = 36*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 36*mi*vtIon*vtIon/(2.0*B0);


  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 4;
  int Nmu = 4;

  double t_end = 2.0e-6; 
  double num_frames = 10;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-2; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.
  double write_phase_freq = 1.0;

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .Te = Te, 
    .Ti = Ti, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nu_frac = nu_frac,
    .B0 = B0, 
    .n0 = n0, 
    .nsource = nsource,
    .Tsource = Tsource,
    .temp_fac = temp_fac,
    .n_fac = n_fac,
    .cx = cx,
    .xcenter = xcenter,
    .cz = cz,
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .write_phase_freq = write_phase_freq,
  };
  return ctx;
}



void
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

void 
physical_phi_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx){
    double R = xn[0], Z = xn[1];
    double n0 = 1.0;
    double background_n = 1.0e-12;

    // Center at t = 1e-6 (shifted up by 0.1)
    double r3 = 2.2;
    double z3 = -1.9; 
    double sigR = 0.1, sigZ = 0.1;

    double dR = R - r3;
    double dZ = Z - z3;
    
    // Calculate normalized radius squared
    double rho2 = (dR * dR) / (sigR * sigR) + (dZ * dZ) / (sigZ * sigZ);

    // Smooth Gaussian: Peak + Background
    // Using exp(-rho2) means sigR is the 1/e width
    fout[0] = (n0 - background_n) * exp(-rho2) + background_n;
}

// INITIAL CONDITION (t = 0)
// Blob starts at -1.7
void 
physical_density_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx){
    double R = xn[0], Z = xn[1];
    double n0 = 1.0;
    double background_n = 1.0e-12;

    // Initial Center
    double r3 = 2.2;
    double z3 = -2.1; 
    double sigR = 0.1, sigZ = 0.1;

    double dR = R - r3;
    double dZ = Z - z3;

    double rho2 = (dR * dR) / (sigR * sigR) + (dZ * dZ) / (sigZ * sigZ);

    fout[0] = (n0 - background_n) * exp(-rho2) + background_n;
}

void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Te;
  fout[0] = T;
}

void
init_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = 2.0*app->Ti;
  fout[0] = T;
}

static inline void
mapc2p_vel_elc(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_elc = app->mu_max_elc;
  double vpar_max_elc = app->vpar_max_elc;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_elc*2.0*pow(cvpar,2);

  mu = mu_max_elc * (cmu * cmu);

  // Set rescaled electron velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar; vp[1] = mu;
}

static inline void
mapc2p_vel_ion(double t, const double* GKYL_RESTRICT vc, double* GKYL_RESTRICT vp, void* ctx)
{
  struct gk_step_ctx *app = ctx;
  double cvpar = vc[0], cmu = vc[1];

  double mu_max_ion = app->mu_max_ion;
  double vpar_max_ion = app->vpar_max_ion;

  double mu = 0.0;
  double vpar = 0.0;

  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vpar = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vpar = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vpar =  vpar_max_ion*2.0*pow(cvpar,2);

  mu = mu_max_ion * (cmu * cmu);

  // Set rescaled ion velocity space coordinates (vpar, mu) from old velocity space coordinates (cvpar, cmu):
  vp[0] = vpar ; vp[1] = mu;
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

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Construct block geometry.
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Ion Species
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[1];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM_PHYS,
      .ctx_density = &ctx,
      .density = init_density,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  struct gkyl_gyrokinetic_bc ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .vdim = ctx.vdim,
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M2},
    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES_NO_BY,
    },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,
    .boundary_flux_diagnostics = {
      .num_diag_moments = 1,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },

    
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 10,
    .bcs = ion_phys_bcs,
  };

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = ctx.B0,
    .time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_bc field_phys_bcs[] = {
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 6, 
    .bcs = field_phys_bcs,
    .time_rate_diagnostics = true,
    .gkfield_id = GKYL_GK_FIELD_ES,
    //.gkfield_id = GKYL_GK_FIELD_FULL_2X,
    //.half_domain=true,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "bump",

    .cdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,
    .gk_block_geom = bgeom,
    .num_species = 1,
    .species = {ion},
    .comm = comm,
    .field = field,
    .skip_field = true,

    .phys_density_func = physical_density_func,
    .phys_phi_func = physical_phi_func,
  };

  struct gkyl_gyrokinetic_run_inp run_inp = {
    .app_type = GKYL_GK_MULTIB,
    .multib_app_inp = app_inp,
    .time_stepping = {
      .t_end = ctx.t_end,
      .num_frames = ctx.num_frames,
      .write_phase_freq = ctx.write_phase_freq,
      .int_diag_calc_num = ctx.int_diag_calc_num,
      .dt_failure_tol = ctx.dt_failure_tol,
      .num_failures_max = ctx.num_failures_max,
      .is_restart = app_args.is_restart,
      .restart_frame = app_args.restart_frame,
      .num_steps = app_args.num_steps,
    },
    .print_verbosity = {
      .enabled = true,
      .frequency = 0.1,
      .estimate_completion_time = true,
    },
  };
  gkyl_gyrokinetic_run_simulation(&run_inp);

  // Create app object.
  //struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new_geom(&app_inp);


  //freeresources:
  // Free resources after simulation completion.
  //gkyl_gyrokinetic_multib_app_release_geom(app);
  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
