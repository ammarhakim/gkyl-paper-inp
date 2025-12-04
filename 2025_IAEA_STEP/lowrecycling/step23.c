#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

// Old outer plates
//void shaped_pfunc_lower_outer(double s, double* RZ){
//  RZ[0] = 3.5+2.0*s;
//  RZ[1] = -8.29;
//}
//
//void shaped_pfunc_upper_outer(double s, double* RZ){
//  RZ[0] = 3.5+2.0*s;
//  RZ[1] = 8.29;
//}
// New outer plates
void shaped_pfunc_lower_outer(double s, double* RZ){
  double p0[2] = {5.488-0.6,-8.600};
  double p1[2] = {5.855-0.6,-8.52318};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  double p0[2] = {5.488-0.6,8.600};
  double p1[2] = {5.855-0.6,8.52318};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}

//old inner plates
void shaped_pfunc_upper_inner(double s, double* RZ){
    RZ[0] = 1.651 + (1.8 - 1.651)*s;
    RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
    RZ[0] = 1.651 + (1.8 - 1.651)*s;
    RZ[1] = -(6.331 + (6.777 - 6.331)*s);
}
//new inner plates
//void shaped_pfunc_upper_inner(double s, double* RZ){
//  double p0[2] = {1.65,6.1};
//  double p1[2] = {1.95,7.2};
//  p1[0] = (p1[0] - p0[0])*2 + p1[0];
//  p1[1] = (p1[1] - p0[1])*2 + p1[1];
//  RZ[0] = (1-s)*p0[0]+s*p1[0];
//  RZ[1] = (1-s)*p0[1]+s*p1[1];
//}
//
//void shaped_pfunc_lower_inner(double s, double* RZ){
//  double p0[2] = {1.65,-6.1};
//  double p1[2] = {1.95,-7.2};
//  p1[0] = (p1[0] - p0[0])*2 + p1[0];
//  p1[1] = (p1[1] - p0[1])*2 + p1[1];
//  RZ[0] = (1-s)*p0[0]+s*p1[0];
//  RZ[1] = (1-s)*p0[1]+s*p1[1];
//}

struct gkyl_block_geom*
create_block_geom(void)
{
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 8);

  /* Block layout and coordinates

   x  
   ^  
   |
   4  +------------------+------------------+------------------+
   |  |b1                |b2                |b3                |
   |  |lower outer SOL   |middle outer sol  |upper outer sol   |
   |  |                  |                  |                  |
   3  +------------------+------------------+------------------+
   |  |b0               x|o b10            %|$ b4              |
   |  |lower outer PF   x|o outer core     %|$ upper outer PF  |
   |  |                 x|o                %|$                 |
   |  +------------------+------------------+------------------+
   2  +------------------+------------------+------------------+
   |  |b9               x|o b11            %|$ b5              |
   |  |lower inner PF   x|o inner core     %|$ upper inner PF  |
   |  |                 x|o                %|$                 |
   1  +------------------+------------------+------------------+
   |  |b8                |b7                |b6                |
   |  |lower inner SOL   |middle inner SOL  |upper inner SOL   |
   |  |                  |                  |                  |
   0  +------------------+------------------+------------------+

      0 -----------1------------2------------3 -> z

      Edges that touch coincide are physically connected unless
      otherwise indicated by a special symbol. Edges with a special
      symbol such as o,x,%, or % are instead connected to the other
      edge with the same symbol. Edges that do not coincide with
      another edge are a physical boundary.
  */  


  struct gkyl_efit_inp efit_inp = {
    .filepath = "./step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  // psisep = 1.5093065418975686; //This is the value
  double dsep = 0.000;
  double wout = 0.1;
  double win = 0.05;
  double wcore = 0.1;
  double wpf = 0.05;

  double psi_lo_outer_sol = psisep - wout;
  double psi_up_outer_sol = psisep;

  double psi_lo_core = psisep;
  double psi_up_core = psisep + wcore;

  double psi_lo_pf = psisep;
  double psi_up_pf = psisep + wpf;

  double psi_lo_inner_sol = psisep - win ;
  double psi_up_inner_sol = psisep;

  int npsi_outer_sol = 8;
  int npsi_core = 8;
  int npsi_pf = 4;
  int npsi_inner_sol = 4;

  double ntheta_lower_inner  = 8;
  double ntheta_middle_inner = 8*4;
  double ntheta_upper_inner  = 8;

  double ntheta_lower_outer = 12;
  double ntheta_middle_outer = 12*4;
  double ntheta_upper_outer = 12;

  double ncuts_lower_inner  = 1;//ntheta_lower_inner;
  double ncuts_middle_inner = 1;//ntheta_middle_inner;
  double ncuts_upper_inner  = 1;//ntheta_upper_inner;

  double ncuts_lower_outer  = 1;//ntheta_lower_outer;
  double ncuts_middle_outer = 1;//ntheta_middle_outer;
  double ncuts_upper_outer  = 1;//ntheta_upper_outer;

  double zinner = 6.34;
  double zouter = 8.29;
  double rright_out = 5.2;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

  // block 0. Lower outer PF region.
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_outer },
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_R,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Lower outer SOL.
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_LO,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 2. Middle outer SOL.
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psi_up_outer_sol, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle_outer},
      .cuts = { 1, ncuts_middle_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_MID,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 0.0,        // closest R to inboard SOL
          .rmin = 0.7,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
          .zmin = -zouter,
          .zmax = zouter,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 6, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

   
  // block 7. Middle inner SOL.
  gkyl_block_geom_set_block(bgeom, 3, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_middle_inner},
      .cuts = { 1, ncuts_middle_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_MID,
          .rleft = 2.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 8. Lower inner SOL.
  gkyl_block_geom_set_block(bgeom, 4, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_lower_inner},
      .cuts = { 1, ncuts_lower_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_LO,
          .rleft = 2.0,
          .rright= rright_out,
          .rmin = 0.0,
          .rmax = 6.2,
          .zmin = -zinner,  
          .zmax = zinner,  
          .plate_spec = true,
          .plate_func_upper = shaped_pfunc_upper_inner,
          .plate_func_lower= shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 3, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_block_geom_set_block(bgeom, 5, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_inner },
      .cuts = { 1, ncuts_lower_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_L,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmin_right = -zouter,
          .zmin_left = -zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_lower_outer,
          .plate_func_upper = shaped_pfunc_lower_inner,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );


  // block 10. outer core.
  gkyl_block_geom_set_block(bgeom, 6, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_core, theta_lo},
      .upper = { psi_up_core,  theta_up},
      .cells = { npsi_core, ntheta_middle_outer},
      .cuts = { 1, ncuts_middle_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_R,
          .rclose = 6.2,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 2.0,        // closest R to inboard SOL
          .rmin = 1.58,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

  // block 11. Inner Core.
  gkyl_block_geom_set_block(bgeom, 7, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_core, theta_lo },
      .upper = { psi_up_core, theta_up },
      .cells = { npsi_core, ntheta_middle_inner},
      .cuts = { 1, ncuts_middle_inner },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE_L,
          .rclose = 0.0,       // Closest R to region of interest
          .rright = rright_out,       // Closest R to outboard SOL
          .rleft = 2.0,        // closest R to inboard SOL
          .rmin = 1.58,         // smallest R in machine
          .rmax = 6.2,         // largest R in machine
          .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 6, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
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
  double massD0; // Hydrogen mass
  double Te; // electron temperature
  double Ti; // ion temperature
  double TD0; // neutral hydrogen temperature
  double vtIon;
  double vtElc;
  double vtD0;
  double nuElc; // electron collision frequency
  double nuIon; // ion collision frequency
  double nuFrac; // Factor to multiply collision frequencies
  double B0; // reference magnetic field
  double n0; // reference density
  double n0D0; // neutral hydrogen reference density
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
  double vpar_max_D0; // Velocity space extents in vparallel for H0
  double t_end; // end time
  int num_frames; // number of output frames
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};



struct gk_step_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // ion mass
  double mD0 = 2.014*GKYL_PROTON_MASS; // D0 mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV; // ion charge
  double qe = -eV; // electron charge


  double n_fac = 32.0;
  double temp_fac = 9.0/n_fac;

  double Te = 1500.0/temp_fac*eV;
  double Ti = 1500.0/temp_fac*eV;
  double TD0 = 1600.0*eV; 
  double B0 = 2.51; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19/n_fac; // Particle density in 1/m^3
  double n0D0 = n0*1.0e-2; // Particle density in 1/m^3
                             
  // Derived parameters.
  double vtIon = sqrt(Ti/mi);
  double vtElc = sqrt(Te/me);
  double vtD0 = sqrt(TD0/mD0);


  double nsource = 1.675e22*2.5*temp_fac/6.0;
  double Tsource = 3000.0*eV/temp_fac;
  double cx = 0.0065612;
  double cz = 0.4916200*1.4;
  double psisep = 1.5093065418975686; //This is the value from efit
  double xcenter = psisep;

  // Collision parameters.
  double nuFrac = 0.25;
  double logLambdaElc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nuElc = nuFrac*logLambdaElc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(me)*(Te*sqrt(Te)));  // collision freq

  double logLambdaIon = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nuIon = nuFrac*logLambdaIon*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*eps0*eps0*sqrt(mi)*(Ti*sqrt(Ti)));

  // Simulation box size (m).

  double vpar_max_elc = 6.0*vtElc;
  double mu_max_elc = 12*me*vtElc*vtElc/(2.0*B0);

  double vpar_max_ion = 8.0*vtIon;
  double mu_max_ion = 12*mi*vtIon*vtIon/(2.0*B0);

  double vpar_max_D0 = 6.0*vtD0;

  // Number of cells.
  int Nx = 4;
  int Nz = 8;
  int Nvpar = 16;
  int Nmu = 12;

  double t_end = 8.0e-3; 
  double num_frames = 800;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-2; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_step_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .chargeElc = qe, 
    .massElc = me, 
    .chargeIon = qi, 
    .massIon = mi,
    .massD0 = mD0,
    .Te = Te, 
    .Ti = Ti, 
    .TD0 = TD0, 
    .vtIon = vtIon,
    .vtElc = vtElc,
    .vtD0 = vtD0,
    .nuElc = nuElc, 
    .nuIon = nuIon, 
    .nuFrac = nuFrac,
    .B0 = B0, 
    .n0 = n0, 
    .n0D0 = n0D0,
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
    .vpar_max_D0 = vpar_max_D0, 
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
init_density_core(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //1.5e20 at inner core and 2e19 at sep
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 3.0e20;
  double intercept = -4.22792e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_outboard_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double cz = 5.640389838180728;
  double n = app->n0D0;
  if (z < 0.0)
    n = n*exp(-cz*(z+M_PI));
  else 
    n = n*exp(-cz*(M_PI-z));
  fout[0] = n;
}

void
init_density_inboard_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double cz = 10.166078433144026;
  double n = app->n0D0;
  if (z < 0.0)
    n = n*exp(-cz*(z+M_PI));
  else 
    n = n*exp(-cz*(M_PI-z));
  fout[0] = n;
}

void
init_density_pfup_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double n = app->n0D0;
  if (z < -1.467030) {
    double cz = 4.125110078248697;
    n = n*exp(-cz*(z+M_PI));
  }
  else {
    double cz = 1.4988763017083815;
    n = n*exp(-cz*(M_PI - z));
  }
  fout[0] = n;
}

void
init_density_pflo_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  struct gk_step_ctx *app = ctx;
  double n = app->n0D0;
  if (z < 1.467030) {
    double cz = 1.4988763017083815;
    n = n*exp(-cz*(z+M_PI));
  }
  else {
    double cz = 4.125110078248697;
    n = n*exp(-cz*(M_PI - z));
  }
  fout[0] = n;
}




void
init_density_empty_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double n = app->n0D0*1e-5;
  fout[0] = n;
}


void
init_density_outer(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 2.7e20;
  double intercept = -3.77513e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_inner(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = 5.4e20;
  double intercept = -7.85026e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}

void
init_density_pf(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  //2e19 at sep 2e17 at outer boundary
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double n0 = app->n0;
  double n_fac = app->n_fac;
  double slope = -5.4e20;
  double intercept = 8.45026e+20;
  double n = slope*x + intercept;
  fout[0] = n/n_fac;
}


void
source_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];

  struct gk_step_ctx *app = ctx;
  double nsource = app->nsource;
  if(x >= 1.5593065418975687)
    fout[0] = nsource;
  else 
    fout[0] = nsource*1.0e-5;
}


void
init_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
init_udrift_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
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
void
source_temp(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->Tsource;
  fout[0] = T;
}


void
init_temp_D0(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *app = ctx;
  double T = app->TD0;
  fout[0] = T;
}

void
init_nu_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuElc;
}

void
init_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_step_ctx *input = ctx;
  fout[0] = input->nuIon;
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



void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_multib_app* app, double t_curr, double dt, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_multib_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_multib_app_calc_integrated_mom(app, t_curr);
    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_multib_app_save_dt(app, t_curr, dt);
  }
}

static int
write_data(struct gkyl_tm_trigger* iot, struct gkyl_tm_trigger* iot_phase, gkyl_gyrokinetic_multib_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  int status = 0;
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Writing data for frame %d\n", frame);

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Applying neutral profile for frame %d\n", frame);
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, 0, 0.0);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Applied neutral profile for frame %d\n", frame);

    gkyl_gyrokinetic_multib_app_write_conf(app, t_curr, frame);
    gkyl_gyrokinetic_multib_app_write_field_energy(app);
    gkyl_gyrokinetic_multib_app_write_integrated_mom(app);
    gkyl_gyrokinetic_multib_app_write_dt(app);
    status = gkyl_gyrokinetic_multib_app_write_flag(app, frame);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Wrote data and flag for frame %d | status = %d\n", frame, status);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_multib_app_write_phase(app, t_curr, frame);
  }

  return status;
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
  struct gkyl_block_geom *bgeom = create_block_geom();
  int nblocks = gkyl_block_geom_num_blocks(bgeom);

  struct gk_step_ctx ctx = create_ctx(); // Context for init functions.

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Elc Species
  // all data is common across blocks
  struct gkyl_gyrokinetic_multib_species_pb elc_blocks[8];
  elc_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };
  elc_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

  };

  elc_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      //.diagnostics = {
      //  .num_diag_moments = 5,
      //  .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
      //  .num_integrated_diag_moments = 1,
      //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      //}
    },

  };

  elc_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_elc,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      //.diagnostics = {
      //  .num_diag_moments = 5,
      //  .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
      //  .num_integrated_diag_moments = 1,
      //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      //}
    },

  };

  struct gkyl_gyrokinetic_block_physical_bcs elc_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },

    //block 10 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 6, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
    //block 11 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 7, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
  };

  struct gkyl_gyrokinetic_multib_species elc = {
    .name = "elc",
    .charge = ctx.chargeElc, .mass = ctx.massElc,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
    .no_by = true,
    //.num_integrated_diag_moments = 1,
    //.integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //.time_rate_diagnostics = true,
    //.boundary_flux_diagnostics = {
    //  .num_diag_moments = 1,
    //  .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //  .num_integrated_diag_moments = 1,
    //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //},

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_elc,
      .num_cross_collisions = 1,
      .bmag_mid = ctx.B0,
      .collide_with = { "ion" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.22 }, 
      .order = 2, 
    }, 

    .react_neut = {
      .num_react = 2,
      .react_type = {
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion", // ion is always the higher charge state
          .donor_nm = "D0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ELC,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "D0",
          .charge_state = 0,
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
      },
    }, 


    .duplicate_across_blocks = false,
    .blocks = elc_blocks,
    .num_physical_bcs = 16,
    .bcs = elc_phys_bcs,
  };


  // Ion Species
  struct gkyl_gyrokinetic_multib_species_pb ion_blocks[8];
  ion_blocks[0] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 0,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[1] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 1,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[2] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 2,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outer,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

  ion_blocks[3] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 3,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[4] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 4,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inner,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };
  ion_blocks[5] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 5,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pf,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

  };

 ion_blocks[6] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 6,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      //.diagnostics = {
      //  .num_diag_moments = 5,
      //  .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
      //  .num_integrated_diag_moments = 1,
      //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      //}
    },

  };

 ion_blocks[7] = (struct gkyl_gyrokinetic_multib_species_pb) {

    .block_id = 7,

    .polarization_density = ctx.n0,

    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_core,
      .ctx_upar = &ctx,
      .upar = init_upar,
      .ctx_temp = &ctx,
      .temp = init_temp_ion,
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = source_density,
        .ctx_upar = &ctx,
        .upar= init_upar,
        .ctx_temp = &ctx,
        .temp = source_temp, 
      }, 
      //.diagnostics = {
      //  .num_diag_moments = 5,
      //  .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
      //  .num_integrated_diag_moments = 1,
      //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      //}
    },

  };

  struct gkyl_gyrokinetic_block_physical_bcs ion_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB },
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},

    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_GK_SHEATH },

    //block 10 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 6, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
    //block 11 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ZERO_FLUX},
    { .bidx = 7, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_REFLECT},
  };


  struct gkyl_gyrokinetic_multib_species ion = {
    .name = "ion",
    .charge = ctx.chargeIon, .mass = ctx.massIon,
    .lower = { -1.0/sqrt(2.0), 0.0},
    .upper = {  1.0/sqrt(2.0), 1.0}, 
    .cells = { cells_v[0], cells_v[1] },
    .num_diag_moments = 7,
    .diag_moments = { "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
    .no_by = true,
    //.num_integrated_diag_moments = 1,
    //.integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //.time_rate_diagnostics = true,
    //.boundary_flux_diagnostics = {
    //  .num_diag_moments = 1,
    //  .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //  .num_integrated_diag_moments = 1,
    //  .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    //},

    
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .nuFrac = ctx.nuFrac,
      .n_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .T_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .ctx = &ctx,
      .self_nu = init_nu_ion,
      .num_cross_collisions = 1,
      .bmag_mid = ctx.B0,
      .collide_with = { "elc" },
    },

    .diffusion = {
      .num_diff_dir = 1, 
      .diff_dirs = { 0 },
      .D = { 0.22 }, 
      .order = 2, 
    }, 

    .react_neut = {
      .num_react = 3,
      .react_type = {
        { .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_D,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .partner_nm = "D0",
          .ion_mass = ctx.massIon,
          .partner_mass = ctx.massIon,
        },
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion", // ion is always the higher charge state
          .donor_nm = "D0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_ION,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "D0",
          .charge_state = 0,
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
      },
    },

  
    .duplicate_across_blocks = false,
    .blocks = ion_blocks,
    .num_physical_bcs = 16,
    .bcs = ion_phys_bcs,
  };

  struct gkyl_gyrokinetic_block_physical_bcs D0_phys_bcs[] = {
    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 0, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 1, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 2, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 2, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 2, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 3, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 4, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 5, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 5, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 6, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 6, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 6, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    { .bidx = 7, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 7, .dir = 1, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},
    { .bidx = 7, .dir = 1, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_SPECIES_ABSORB},

    };

  struct gkyl_gyrokinetic_multib_neut_species_pb D0_blocks[8];
  D0_blocks[0] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 0,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pflo_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };

  D0_blocks[1] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 1,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_outboard_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };

  D0_blocks[2] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 2,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_empty_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };

  D0_blocks[3] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 3,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_empty_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };
  D0_blocks[4] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 4,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_inboard_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };
  D0_blocks[5] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 5,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_pflo_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },

  };

 D0_blocks[6] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 6,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_empty_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },


  };

 D0_blocks[7] = (struct gkyl_gyrokinetic_multib_neut_species_pb) {

    .block_id = 7,


    .projection = {
      .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
      .ctx_density = &ctx,
      .density = init_density_empty_D0,
      .ctx_udrift = &ctx,
      .udrift = init_udrift_D0,
      .ctx_temp = &ctx,
      .temp = init_temp_D0,
    },


  };



  struct gkyl_gyrokinetic_multib_neut_species D0 = {
    .name = "D0",
    .mass = ctx.massD0,
    .lower = { -ctx.vpar_max_D0, -ctx.vpar_max_D0, -ctx.vpar_max_D0},
    .upper = {  ctx.vpar_max_D0, ctx.vpar_max_D0, ctx.vpar_max_D0}, 
    .cells = { 8, 8, 8 },
    .is_static = true,
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2"},

    .react_neut = {
      .num_react = 3,
      .react_type = {
        { .react_id = GKYL_REACT_CX,
          .type_self = GKYL_SELF_PARTNER,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .partner_nm = "D0",
          .ion_mass = ctx.massIon,
          .partner_mass = ctx.massIon,
        },
        { .react_id = GKYL_REACT_IZ,
          .type_self = GKYL_SELF_DONOR,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion", // ion is always the higher charge state
          .donor_nm = "D0", // interacts with elc to give up charge
          .charge_state = 0, // corresponds to lower charge state (donor)
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
        { .react_id = GKYL_REACT_RECOMB,
          .type_self = GKYL_SELF_RECVR,
          .ion_id = GKYL_ION_H,
          .elc_nm = "elc",
          .ion_nm = "ion",
          .recvr_nm = "D0",
          .charge_state = 0,
          .ion_mass = ctx.massIon,
          .elc_mass = ctx.massElc,
        },
      },
    },

    .duplicate_across_blocks = false,
    .blocks = D0_blocks,
    .num_physical_bcs = 32,
    .bcs = D0_phys_bcs,
  };



 

  // Field object
  struct gkyl_gyrokinetic_multib_field_pb field_blocks[1];
  field_blocks[0] = (struct gkyl_gyrokinetic_multib_field_pb) {
    .polarization_bmag = ctx.B0,
  };

  struct gkyl_gyrokinetic_block_physical_bcs field_phys_bcs[] = {
    // block 1 BCs
    { .bidx = 1, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    // block 2 BCs
    { .bidx = 2, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},

    // block 3 BCs
    { .bidx = 3, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    // block 4 BCs
    { .bidx = 4, .dir = 0, .edge = GKYL_LOWER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},

    // block 0 BCs
    { .bidx = 0, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},
    // block 5 BCs
    { .bidx = 5, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_DIRICHLET},

    // block 6 BCs
    { .bidx = 6, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NEUMANN},
    // block 7 BCs
    { .bidx = 7, .dir = 0, .edge = GKYL_UPPER_EDGE, .bc_type = GKYL_BC_GK_FIELD_NEUMANN},
  };

  struct gkyl_gyrokinetic_multib_field field = {
    .duplicate_across_blocks = true,
    .blocks = field_blocks, 
    .num_physical_bcs = 8, 
    .bcs = field_phys_bcs,
    //.time_rate_diagnostics = true,
  };

  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "step23",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,
    .cfl_frac_omegaH = 1.7,

    .block_geom = bgeom,
    
    .num_species = 2,
    .species = { elc, ion},

    .num_neut_species = 1,
    .neut_species = { D0 },

    .field = field,
    //.skip_field=true,

    .comm = comm
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new(&app_inp);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_multib_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_multib_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_multib_app_apply_ic(app, t_curr);
  }

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, -1.0, false);
  write_data(&trig_write, &trig_write_phase, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  // Flag for data writing
  int wrote = 0;
  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if(step%1000 == 0)
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = gkyl_gyrokinetic_multib_update(app, dt);
    if(step%1000 == 0)
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_multib_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, status.dt_actual, t_curr > t_end);
    wrote = write_data(&trig_write, &trig_write_phase, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_multib_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_multib_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_multib_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, status.dt_actual, t_curr, true);
        write_data(&trig_write, &trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_multib_app_stat_write(app);

  // Fetch simulation statistics.
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_multib_app_stat(app);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_multib_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_multib_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_multib_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release(app);
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
