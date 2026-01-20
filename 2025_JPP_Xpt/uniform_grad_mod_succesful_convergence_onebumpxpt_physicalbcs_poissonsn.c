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

  //int ncuts_lower = 10;//ntheta_lower/8;
  //int ncuts_middle = 10;//ntheta_middle/8;
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
        { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
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
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
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
        { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL},
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
        { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL},
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
        { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
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
        { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL},
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

  double t_end = 1.0e-10; 
  double num_frames = 1;
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
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_step_ctx *app = ctx;

  fout[0] = 0.5; // Diffusivity [m^2/s].
}

void
init_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  fout[0] = n0;
}

// 1. Helper: Calculates the analytic source 'S' for a single bump
//    S = -Div(w Grad phi_bump)
double eval_bump_source(double R, double Z, double R0, double Z0, double sigma, double w) {
    double r2 = (R-R0)*(R-R0) + (Z-Z0)*(Z-Z0);
    if (r2 >= sigma*sigma) return 0.0;
    return ((pow(sigma, 2) > pow(R - R0, 2) + pow(Z - Z0, 2)) ? (
   -2*pow(sigma, 2)*w*(2*R*(pow(sigma, 2)*pow(R - R0, 2) + pow(sigma, 2)*pow(Z - Z0, 2) + 2*pow(R - R0, 2)*(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2)) + 2*pow(Z - Z0, 2)*(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2)) - pow(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2), 2)) - (R - R0)*pow(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2), 2))*exp((pow(R - R0, 2) + pow(Z - Z0, 2))/(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2)))/(R*pow(-pow(sigma, 2) + pow(R - R0, 2) + pow(Z - Z0, 2), 4))
)
: (
   0
));
}

void
physical_density_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
    double R = xn[0], Z = xn[1];
    double n_fac = 32.0;
    double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
    double eV=1.602e-19;
    double w = n0*3.34e-27;
    // --- Configuration: Two Bumps (Lower Single Null) ---
    double r1=2.2, z1=-2.2, sig1=0.2;
    double r2=2.0, z2=-1.7, sig2=0.2;


    double r3=2.0, z3=-2.0, sig3=0.2;

    // Calculate Analytic Source (Total Rho required for Poisson)
    //double rho_ana = eval_bump_source(R, Z, r1, z1, sig1, w) + eval_bump_source(R, Z, r2, z2, sig2, w);
    double rho_ana =  eval_bump_source(R, Z, r3, z3, sig3, w);

    // Try a gaussian
    //double s = (R-r2)*(R-r2) + (Z-z2)*(Z-z2);
    //double rho_ana = w*4/sig2/sig2 * exp(-s/sig2/sig2) * (1-s/sig2/sig2);

    double rho_electron_offset = eV * n0;
    double rho_ion = rho_ana + rho_electron_offset;
    fout[0] = rho_ion/eV;
}

double eval_bump_phi(double R, double Z, double R0, double Z0, double sigma)
{
    double r2 = (R-R0)*(R-R0) + (Z-Z0)*(Z-Z0);
    if (r2 >= sigma*sigma) return 0.0;
    return ((pow(sigma, 2) > pow(R - R0, 2) + pow(Z - Z0, 2)) ? ( exp(-pow(sigma, 2)/(pow(sigma, 2) - pow(R - R0, 2) - pow(Z - Z0, 2)) + 1)) : ( 0));
}
void
physical_phi_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
    double R = xn[0], Z = xn[1];
    double r1=2.2, z1=-2.2, sig1=0.2;
    double r2=2.0, z2=-1.7, sig2=0.2;

    double r3=2.0, z3=-2.0, sig3=0.2;


    //double phi_ana = eval_bump_phi(R, Z, r1, z1, sig1) + eval_bump_phi(R, Z, r2, z2, sig2);
    double phi_ana = eval_bump_phi(R, Z, r3, z3, sig3);


    // Try a gaussian
    //double s = (R-r2)*(R-r2) + (Z-z2)*(Z-z2);
    //double phi_ana = exp(-s/sig2/sig2);

    fout[0] = phi_ana;

}

void 
new_physical_phi_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx){

    
    double R = xn[0], Z = xn[1];
    double r3=2.0, z3=-2.0, sigR=0.2, sigZ = 0.2;

    double dR = R - r3;
    double dZ = Z - z3;
    double rho2 = (dR * dR) / (sigR * sigR) + (dZ * dZ) / (sigZ * sigZ);

    if (rho2 >= 1.0) {
      fout[0] = 0.0;
      return;
    }

    double u = 1.0 - rho2;
    double inv_u = 1.0 / u;
    double phi = exp(1.0 - inv_u);
    fout[0] = phi;

}
void 
new_physical_density_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx){

    double R = xn[0], Z = xn[1];
    double n_fac = 32.0;
    double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
    double eV= 1.602176634e-19;
    double w = n0*2.014*1.67262192e-27;

    double r3=2.0, z3=-2.0, sigR=0.2, sigZ = 0.2;

    double dPsi_dR = 2.0 * (R - 2.0);
    double dPsi_dZ = 2.0 * Z + (Z * Z);
    
    // B_R = -(1/R) * dPsi/dZ
    // B_Z =  (1/R) * dPsi/dR
    double B_R = -(1.0 / R) * dPsi_dZ;
    double B_Z =  (1.0 / R) * dPsi_dR;
    double F0=2.0;
    double B_phi = F0 / R;
    
    // Magnitude |B|
    double B_mag = sqrt(B_R*B_R + B_Z*B_Z + B_phi*B_phi);
    
    // Unit vectors
    double bR = B_R / B_mag;
    double bZ = B_Z / B_mag;


    double dR = R - r3;
    double dZ = Z - z3;
    double rho2 = (dR * dR) / (sigR * sigR) + (dZ * dZ) / (sigZ * sigZ);

    if (rho2 >= 1.0) {
      double rho_ana = 0.0;
      double rho_electron_offset = eV * n0;
      double rho_ion = rho_ana + rho_electron_offset;
      fout[0] = rho_ion/eV;
      return;
    }

    double u = 1.0 - rho2;
    double inv_u = 1.0 / u;
    double phi = exp(1.0 - inv_u);

    // Derivatives of the exponent function E(u) = 1 - 1/u
    // Let E = 1 - u^(-1)
    // dE/du = u^(-2) = 1/u^2
    // d2E/du2 = -2u^(-3) = -2/u^3
    double dE_du = inv_u * inv_u;
    double d2E_du2 = -2.0 * (inv_u * inv_u * inv_u);

    // Derivatives of u = 1 - ( (R-R0)^2/a^2 + ... ) with respect to R
    // du/dR = -2(R-R0)/a^2
    // d2u/dR2 = -2/a^2
    double du_dR = -2.0 * dR / (sigR * sigR);
    double d2u_dR2 = -2.0 / (sigR * sigR);
    double du_dZ = -2.0 * dZ / (sigZ * sigZ);
    double d2u_dZ2 = -2.0 / (sigZ * sigZ);

    // 4. Compute First Derivatives of Phi
    // dPhi/dR = Phi * (dE/du) * (du/dR)
    double dphi_dR = phi * dE_du * du_dR;
    double dphi_dZ = phi * dE_du * du_dZ;

    // 5. Compute Second Derivatives of Phi
    // d2Phi/dR2 = d/dR [ Phi * E' * u' ]
    // Product rule: (Phi') * (E' * u') + Phi * (E'' * (u')^2 + E' * u'')
    // Note: Phi' = Phi * E' * u'
    // So: Phi * [ (E' u')^2 + E'' (u')^2 + E' u'' ]
    
    double term_R = (dE_du * du_dR) * (dE_du * du_dR)  + // (E' u')^2
                    d2E_du2 * (du_dR * du_dR)          + // E'' (u')^2
                    dE_du * d2u_dR2;                     // E' u''
    double d2phi_dR2 = phi * term_R;

    double term_Z = (dE_du * du_dZ) * (dE_du * du_dZ)  + 
                    d2E_du2 * (du_dZ * du_dZ)          + 
                    dE_du * d2u_dZ2;
    double d2phi_dZ2 = phi * term_Z;

    // 6. Assemble Laplacian (Cylindrical: d2/dR2 + (1/R)d/dR + d2/dZ2)
    double laplacian_full = d2phi_dR2 + (1.0 / R) * dphi_dR + d2phi_dZ2;

    // Cross derivative d2Phi/dRdZ
    double d2phi_dRdZ = (dphi_dZ * dE_du * du_dR) + (phi * d2E_du2 * du_dZ * du_dR);

    // --- 6. The Parallel Laplacian Correction ---
    // grad_par^2 approx (b.grad)^2
    double par_term_2nd = (bR * bR * d2phi_dR2) + 
                          (2.0 * bR * bZ * d2phi_dRdZ) + 
                          (bZ * bZ * d2phi_dZ2);
                          
    double par_term_1st = (bR / R) * (bR * dphi_dR + bZ * dphi_dZ);
    
    double laplacian_par = par_term_2nd + par_term_1st;

    // --- 7. The Perpendicular Source ---
    double laplacian_perp = laplacian_full - laplacian_par;

    // C. Parallel Gradient
    double grad_par = bR * dphi_dR + bZ * dphi_dZ;

    // D. Divergence of b (The New Term)
    // Identity: div(b) = - (b . grad(B_mag)) / B_mag
    // Derivatives of B components for Cubic Single Null:
    // dB_R/dR = -B_R/R;  dB_Z/dR = 4/R^2;  dB_phi/dR = -B_phi/R
    // dB_R/dZ = -(2+2Z)/R; dB_Z/dZ = 0;    dB_phi/dZ = 0
    
    double dBR_dR = -B_R / R;
    double dBZ_dR = 4.0 / (R * R);
    double dBphi_dR = -B_phi / R;
    
    double dBR_dZ = -(2.0 + 2.0*Z) / R;
    
    // Grad |B|
    double dBmag_dR = (B_R*dBR_dR + B_Z*dBZ_dR + B_phi*dBphi_dR) / B_mag;
    double dBmag_dZ = (B_R*dBR_dZ) / B_mag; // Only B_R depends on Z here

    double div_b = - (bR * dBmag_dR + bZ * dBmag_dZ) / B_mag;

    // --- Final Result ---
    double lap_mod = laplacian_perp - (div_b * grad_par);
    
    double rho_ana = -w * lap_mod;

    // Combine with background
    double rho_electron_offset = eV * n0;
    double rho_ion = rho_ana + rho_electron_offset;
    fout[0] = rho_ion/eV;
}

// 2. Correct GK Source
//void new_physical_density_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
//{
//    double R = xn[0], Z = xn[1];
//    double n_fac = 32.0;
//    double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
//    double eV=1.602e-19;
//    double rho_electron_offset = eV * n0;
//
//    double w = n0*3.34e-27;
//    double R0=2.0, Z0=-2.0, sigma=0.2;
//    double r2 = (R-R0)*(R-R0) + (Z-Z0)*(Z-Z0);
//
//    if (r2 >= sigma*sigma) {
//      fout[0] =  rho_electron_offset/eV;
//      return;
//    }
//    double x0 = pow(sigma, 2);
//    double x1 = x0 - pow(R - R0, 2) - pow(Z - Z0, 2);
//    double x2 = pow(x1, -2);
//    double x3 = exp(-x0/x1 + 1);
//    double x4 = x0*x3;
//    double x5 = x2*x4;
//    double x6 = 2*x5;
//    double x7 = pow(sigma, 4);
//    double x8 = 2*Z;
//    double x9 = -2*Z0 + x8;
//    double x10 = pow(x9, 2);
//    double x11 = pow(x1, -4);
//    double x12 = x4/pow(x1, 3);
//    double x13 = x12*(4*Z - 4*Z0);
//    double x14 = x13*x9;
//    double x15 = pow(R, -3);
//    double x16 = 2.0*R - 4.0;
//    double x17 = 2.0*Z + 2;
//    double x18 = pow(Z, 2);
//    double x19 = Z + 0.5*x18;
//    double x20 = pow(R, -2);
//    double x21 = 0.25*x20;
//    double x22 = 0.5*R;
//    double x23 = pow(x22 - 1, 2);
//    double x24 = pow(x19, 2);
//    double x25 = x20*x23 + x21*x24 + x21;
//    double x26 = pow(x25, -3.0/2.0);
//    double x27 = 1.0/R;
//    double x28 = 0.25*x27;
//    double x29 = pow(x25, -1.0/2.0);
//    double x30 = x16*x29;
//    double x31 = x28*x30;
//    double x32 = x31*x9;
//    double x33 = 2*R - 2*R0;
//    double x34 = 1.0*x18 + x8;
//    double x35 = 0.25*x0*x2*x27*x29*x3*x33*x34 - x32*x5;
//    double x36 = x30*x5;
//    double x37 = 0.5*x27;
//    double x38 = x33*x5;
//    double x39 = x11*x3*x7;
//    double x40 = 0.03125*x15*x17*x19;
//    double x41 = x26*x5*x9;
//    double x42 = x29*x34;
//    double x43 = x28*x42;
//    double x44 = x33*x43;
//    double x45 = x26*x38;
//    double x46 = pow(x33, 2)*x39;
//    double x47 = x12*(4*R - 4*R0);
//    double x48 = x33*x47;
//    double x49 = x21*x42;
//    double x50 = 0.25*x15;
//    double x51 = x28*(x15*x23 - 1.0/2.0*x20*(x22 - 1.0) + x24*x50 + x50);
//    double x52 = x34*x51;
//    double x53 = x37*x5;
//
//    // --- Final Result ---
//    double rho_ana = -w*(x10*x11*x3*x7 - x14 + 0.03125*x15*x16*x17*x19*x26*x35 - x31*(x10*x31*x39 + x13*x44 - x14*x31 + x16*x40*x41 + x17*x28*x29*x38 - x34*x40*x45 - x36*x37 - x39*x44*x9) - x6) - x27*(R*w*(x26*x35*x52 - x35*x49 + x43*(-x16*x41*x51 + x21*x36*x9 - x29*x53*x9 + x32*x33*x39 - x32*x47 - x38*x49 + x42*x53 - x43*x46 + x43*x48 + x45*x52) + x46 - x48 - x6) + w*(0.25*x27*x29*x34*x35 - x38));
//
//
//    double rho_ion = rho_ana + rho_electron_offset;
//    fout[0] = rho_ion/eV;
//
//}


//void
//phi_func(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
//{
//    double R = xn[0], Z = xn[1];
//    double n_fac = 32.0;
//    double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
//    double eV=1.602e-19;
//    double w = n0*3.34e-27/2.51/2.51;
//    // --- Configuration: Two Bumps (Lower Single Null) ---
//    // Bump 1: Core O-point
//    double r1=2.0, z1=0.0, sig1=0.4;
//    // Bump 2: X-point Singularity
//    double r2=2.0, z2=-2.0, sig2=0.3;
//
//    // Calculate Analytic Source (Total Rho required for Poisson)
//    //double rho_ana = eval_bump_source(R, Z, r1, z1, sig1, w) + eval_bump_source(R, Z, r2, z2, sig2, w);
//    double phi_ana = eval_bump_phi(R, Z, r2, z2, sig2);
//
//    fout[0] = phi_ana;
//}

// ------------------------------------------------------------------
// AUTO-GENERATED SINUSOIDAL TEST FUNCTIONS
// Phi = A * cos(kR*(R-Rc)) * cos(kZ*(Z-Zc))
// ------------------------------------------------------------------
#include <math.h>

// 1. Target Potential Phi
double eval_phi_sine(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx) {
    // --- Settings ---
    double R = xn[0], Z = xn[1];
    double A = 1.0;
    // Wavenumbers: Adjust these to fit waves in your domain
    // Example: pi / 1.0 means half-wavelength is 1.0m
    double kR = 3.1415926535 / 1.0;
    double kZ = 3.1415926535 / 2.0;
    // Phase Shift: Center a peak at the X-point (2.0, -2.0)
    double Rc = 2.0;
    double Zc = -2.0;

    return 1.0*cos(kR*(R - Rc))*cos(kZ*(Z - Zc));
}

// 2. Analytic Source Rho (Charge Density)
double eval_rho_sine_source(double R, double Z) {
    // --- Same Settings (Must Match Phi!) ---
    double A = 1.0;
    double kR = 3.1415926535 / 1.0;
    double kZ = 3.1415926535 / 2.0;
    double Rc = 2.0;
    double Zc = -2.0;

    // Analytic Laplacian in (R,Z) Cylindrical
    return 1.0*(R*pow(kR, 2)*cos(kR*(R - Rc)) + R*pow(kZ, 2)*cos(kR*(R - Rc)) + kR*sin(kR*(R - Rc)))*cos(kZ*(Z - Zc))/R;
}

// 3. Initialization Function for Ion Particle Density
void
eval_ion_density_sine(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx) {
    double n_fac = 32.0;
    double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
    double eV=1.602e-19;
    double R = xn[0], Z = xn[1];
    // Get the analytic source required
    double rho_ana = eval_rho_sine_source(R, Z);

    // Convert to ion particle density:
    // ni = (rho_ana + e*n0) / q_ion
    fout[0] = (rho_ana + eV * n0) / eV;
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

  // Elc Species
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[1];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  struct gkyl_gyrokinetic_bc elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    //block 10 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 6, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
    //block 11 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 7, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .vdim = ctx.vdim,
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0 },
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
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nu_frac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .bmag_ref = ctx.B0,
      .collide_with = { "ion" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
      .write_diagnostics=true,
    },

    .duplicate_across_blocks = true,
    .blocks = elc_blocks,
    .num_physical_bcs = 16,
    .bcs = elc_phys_bcs,
  };

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
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH },

    //block 10 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 6, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
    //block 11 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 7, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_REFLECT},
  };

  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .vdim = ctx.vdim,
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { GKYL_F_MOMENT_M0},
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

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .nu_frac = ctx.nu_frac,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .bmag_ref = ctx.B0,
      .collide_with = { "elc" },
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
      .write_diagnostics=true,
    },

    .duplicate_across_blocks = true,
    .blocks = ion_blocks,
    .num_physical_bcs = 16,
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


    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },

    //{ .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 4, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    //{ .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },

    //{ .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    //{ .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },

    // Custom BCs for xpt bump (natural bcs)
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },

    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    { .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },

    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    { .bidx = 4, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0} },
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    { .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },

    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },

    // Custom BCs for two bumps
    //{ .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 4, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },
    //{ .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0} },

    //{ .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
    //{ .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_PERIODIC, .value = {0.0} },
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 6+6+12, 
    .bcs = field_phys_bcs,
    .time_rate_diagnostics = true,
    //.gkfield_id = GKYL_GK_FIELD_ES,
    .gkfield_id = GKYL_GK_FIELD_FULL_2X,
    //.half_domain=true,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "poissonsn",

    .cdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,
    .gk_block_geom = bgeom,
    .num_species = 2,
    .species = { elc, ion},
    .field = field,
    .comm = comm,

    .phys_density_func = new_physical_density_func,
    .phys_phi_func = new_physical_phi_func,
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
