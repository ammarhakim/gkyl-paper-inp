#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#endif

struct gkyl_gk_block_geom*
create_gk_block_geom(void)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_gk_block_geom_new(2, 3);

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
    .filepath = "./double_null.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g, psisep = %g\n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry, efit->psisep);
  gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, "double_null_psi.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->fpolflux, "double_null_fpol.gkyl");
  gkyl_grid_sub_array_write(&efit->fluxgrid, &efit->fluxlocal, 0, efit->qflux, "double_null_q.gkyl");
  double psisep = efit->psisep;
  gkyl_efit_release(efit);
  // psisep = 1.5093065418975686; //This is the value
  double wfac = 4.0;
  double wout = 0.1*wfac;
  double win = 0.1*wfac;
  double wcore = 0.1*wfac;
  double wpf = 0.1*wfac;

  double psi_lo_outer_sol = psisep;
  //double psi_lo_outer_sol = psisep + wout/4.0;
  double psi_up_outer_sol = psisep+wout;

  double psi_lo_core = psisep;
  double psi_up_core = psisep + wcore;

  double psi_lo_pf = psisep - wpf;
  double psi_up_pf = psisep;

  double psi_lo_inner_sol = psisep - win ;
  double psi_up_inner_sol = psisep;

  double cellfac = 16.0;

  int npsi_outer_sol = 1*cellfac;
  int npsi_core = 1*cellfac;
  int npsi_pf = 1*cellfac;
  int npsi_inner_sol = 1*cellfac;

  double ntheta_lower_inner  = 1*cellfac;
  double ntheta_middle_inner = 2*cellfac;
  double ntheta_upper_inner  = 1*cellfac;

  double ntheta_lower_outer = 1*cellfac;
  double ntheta_middle_outer = 2*cellfac;
  double ntheta_upper_outer = 1*cellfac;

  enum gkyl_geometry_id geo_type = GKYL_TOKAMAK;
  int ncuts_lower_outer = 1;
  int ncuts_middle_outer = 1;
  int ncuts_upper_outer = 1;

  double zinner = 2.75;
  double zouter = 2.75;
  double rmax = 4.0;
  double rmin = 0.0;

  double Lz = (M_PI-1e-14)*2.0;
  double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

  double compression_factor = 0.0;
  double radial_compression_factor = 0.0;

 // block 0. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = geo_type,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_LO,
          .rclose = rmax,       // Closest R to region of interest
          .rright = rmax,       // Closest R to outboard SOL
          .rleft = rmin,        // closest R to inboard SOL
          .rmin = rmin,         // smallest R in machine
          .rmax = rmax,         // largest R in machine
          .zmin = -zouter,
          .zmax = zouter,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // block 1. Middle outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo },
      .upper = { psi_up_outer_sol, theta_up },
      .cells = { npsi_outer_sol, ntheta_middle_outer},
      .cuts = { 1, ncuts_middle_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = geo_type,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_MID,
          .rclose = rmax,       // Closest R to region of interest
          .rright = rmax,       // Closest R to outboard SOL
          .rleft = rmin,        // closest R to inboard SOL
          .rmin = rmin,         // smallest R in machine
          .rmax = rmax,         // largest R in machine
          .zmin = -zouter,
          .zmax = zouter,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // z-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

 // block 2. Lower outer SOL.
  gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
      .lower = { psi_lo_outer_sol, theta_lo},
      .upper = { psi_up_outer_sol,  theta_up},
      .cells = { npsi_outer_sol, ntheta_lower_outer},
      .cuts = { 1, ncuts_lower_outer },
      .geometry = {
        .world = {0.0},
        .geometry_id = geo_type,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_DN_SOL_OUT_UP,
          .rclose = rmax,       // Closest R to region of interest
          .rright = rmax,       // Closest R to outboard SOL
          .rleft = rmin,        // closest R to inboard SOL
          .rmin = rmin,         // smallest R in machine
          .rmax = rmax,         // largest R in machine
          .zmin = -zouter,
          .zmax = zouter,
        },
        .position_map_info = {
          .id = GKYL_PMAP_XPT_COMPRESSION,
          .radial_compression_factor = radial_compression_factor,
          .compression_factor = compression_factor
        }
      },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }
      },
      .connections[1] = { // z-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
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
  struct gkyl_gk_block_geom *bgeom = create_gk_block_geom();
  int nblocks = gkyl_gk_block_geom_num_blocks(bgeom);


  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "griddn2x",

    .cdim = 2,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .use_gpu = app_args.use_gpu,
    .cfl_frac = 1.0,

    .gk_block_geom = bgeom,
    .comm = comm
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new_geom(&app_inp);


  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_multib_app_release_geom(app);
  gkyl_gk_block_geom_release(bgeom);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif

  return 0;
}
