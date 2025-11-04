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
  //RZ[0]=5.0;
  //RZ[1] = -(8.3+s/5.0);
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  double p0[2] = {5.488-0.6,8.600};
  double p1[2] = {5.855-0.6,8.52318};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
  //RZ[0]=5.0;
  //RZ[1] = (8.3+s/5.0);
}

//Old inner plates
//void shaped_pfunc_upper_inner(double s, double* RZ){
//    RZ[0] = 1.651 + (1.8 - 1.651)*s;
//    RZ[1] = 6.331 + (6.777 - 6.331)*s;
//}
//
//void shaped_pfunc_lower_inner(double s, double* RZ){
//    RZ[0] = 1.651 + (1.8 - 1.651)*s;
//    RZ[1] = -(6.331 + (6.777 - 6.331)*s);
//}
//new inner plates
void shaped_pfunc_upper_inner(double s, double* RZ){
  double p0[2] = {1.65,6.1};
  double p1[2] = {1.95,7.2};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}

void shaped_pfunc_lower_inner(double s, double* RZ){
  double p0[2] = {1.65,-6.1};
  double p1[2] = {1.95,-7.2};
  p1[0] = (p1[0] - p0[0])*2 + p1[0];
  p1[1] = (p1[1] - p0[1])*2 + p1[1];
  RZ[0] = (1-s)*p0[0]+s*p1[0];
  RZ[1] = (1-s)*p0[1]+s*p1[1];
}


struct gkyl_block_geom*
create_block_geom(void)
{
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 12);

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
  double wout = 0.125;
  double win = 0.05;
  double wcore = 0.1*2;
  double wpf = 0.05;

  double psi_lo_outer_sol = psisep - wout;
  double psi_up_outer_sol = psisep;

  double psi_lo_core = psisep;
  double psi_up_core = psisep + wcore;

  double psi_lo_pf = psisep;
  double psi_up_pf = psisep + wpf;

  double psi_lo_inner_sol = psisep - win ;
  double psi_up_inner_sol = psisep;

  int npsi_outer_sol = 5;
  int npsi_core = 4*2;
  int npsi_pf = 2;
  int npsi_inner_sol = 2;

  double ntheta_lower_inner  = 4;
  double ntheta_middle_inner = 16;
  double ntheta_upper_inner  = 4;

  double ntheta_lower_outer = 8;
  double ntheta_middle_outer = 24;
  double ntheta_upper_outer = 8;

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
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 1. Lower outer SOL.
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
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
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
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
        { .bid = 10, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 3. Upper outer SOL.
  gkyl_block_geom_set_block(bgeom, 3, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol, theta_up},
      .cells = { npsi_outer_sol, ntheta_upper_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        //.geometry_id = GKYL_GEOMETRY_FROMFILE,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_UP,
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
        { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
      }
    }
  );

  // block 4. Upper outer PF region.
  gkyl_block_geom_set_block(bgeom, 4, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_R,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );
  
  // block 5. Upper inner PF region.
  gkyl_block_geom_set_block(bgeom, 5, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_upper_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_UP_L,
          .rright = rright_out,
          .rleft = 0.0,
          .rmin = 1.7,
          .rmax = 6.2,
          .zmax_right = zouter,
          .zmax_left = zinner,
          .plate_spec = true,
          .plate_func_lower = shaped_pfunc_upper_inner,
          .plate_func_upper = shaped_pfunc_upper_outer,
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 6. Upper inner SOL.
  gkyl_block_geom_set_block(bgeom, 6, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_upper_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_IN_UP,
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
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 7. Middle inner SOL.
  gkyl_block_geom_set_block(bgeom, 7, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_middle_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 11, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 8. Lower inner SOL.
  gkyl_block_geom_set_block(bgeom, 8, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_inner_sol, theta_lo },
      .upper = { psi_up_inner_sol, theta_up },
      .cells = { npsi_inner_sol, ntheta_lower_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
      },
      .connections[1] = { // z-direction connections
        { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
      }
    }
  );

  // block 9. Lower inner PF region.
  gkyl_block_geom_set_block(bgeom, 9, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_pf, theta_lo},
      .upper = { psi_up_pf, theta_up},
      .cells = { npsi_pf, ntheta_lower_inner },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 8, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
      }
    }
  );


  // block 10. outer core.
  gkyl_block_geom_set_block(bgeom, 10, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_core, theta_lo},
      .upper = { psi_up_core,  theta_up},
      .cells = { npsi_core, ntheta_middle_outer},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 11. Inner Core.
  gkyl_block_geom_set_block(bgeom, 11, &(struct gkyl_block_geom_info) {
      .lower = { psi_lo_core, theta_lo },
      .upper = { psi_up_core, theta_up },
      .cells = { npsi_core, ntheta_middle_inner},
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
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
        { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
      },
      .connections[1] = { // z-direction connections
        { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
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


  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "gridstep2x",

    .cdim = 2, .vdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,

    .block_geom = bgeom,
    .comm = comm
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new_geom(&app_inp);


  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release_geom(app);
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
