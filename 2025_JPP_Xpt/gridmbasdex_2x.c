#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>
#include <gkyl_tok_geo.h>

#include <rt_arg_parse.h>

struct gk_asdex_ctx {
  int cdim, vdim; // Dimensionality.

  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass.
  double Te; // Electron temperature
  double Ti; // Ion temperature
  double c_s; // sound speed
  double nu_elc; // electron collision frequency
  double nu_ion; // ion collision frequency
  double B0; // Magnetic field.
  double n0; // Density.

  // Source parameters
  double psi_src; // Location.
  double lambda_src; // Width.
  double ndot_src; // Particle source rate.
  double Te_src; // Electron source temperature.
  double Ti_src; // Ion source temperature.

  // Domain parameters.            
  char geqdsk_file[128]; // File with equilibrium.
  double psi_axis; // Psi at the magnetic axis.
  double psi_sep; // Psi at the separatrix.
  double psi_min_core, psi_max_sol, psi_min_pf; // Psi extents.
  double Lx_core; // Box size in x
  // Z location of the X-point on psi=psi_min.
  double z_xpt_psi_sep_lo, z_xpt_psi_sep_up;

  // Grid.
  int Npsi_sol; // Number of cells in psi in the SOL.
  int Npsi_pf; // Number of cells in psi in the private flux.
  int Npsi_core; // Number of cells in psi in the core.
  int Ntheta_divertor; // Number of cells in theta in the divertor.
  int Ntheta_sol; // Number of cells in theta in the (upper) SOL (and core).
  int num_blocks; // Number of blocks.

  // Physical velocity space limits
  double vpar_max_elc; // Parallel velocity extents for electrons.
  double mu_max_elc; // Maximum magnetic moment for electrons.
  double vpar_max_ion; // Parallel velocity extents for ions.
  double mu_max_ion; // Maximum magnetic moment for ions.

  // Computational velocity space limits
  double vpar_min_elc_c, vpar_max_elc_c;
  double mu_min_elc_c, mu_max_elc_c;
  double vpar_min_ion_c, vpar_max_ion_c;
  double mu_min_ion_c, mu_max_ion_c;

};



struct gkyl_block_geom*
create_asdex_lsn_block_geom(void *ctx)
{
  struct gk_asdex_ctx *params = ctx;

  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(params->cdim, params->num_blocks);

  /* Block layout and coordinates.

    z  
    ^  
    |

    |                     +------------------+------------------+
    |                     |                  |                  | 
    |                     |b4                |b3                |
    |                     |inner PF          |lower inner sol   |
    |                     |                  |                  | 
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |  +------------------+------------------+------------------+
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  | b5              $|                  |$ b2              |
    |  + core            $+                  +$ upper sol       |
    |  |                 $|                  |$                 |
    |  |                 $|                  |$                 |
    |  +------------------+------------------+------------------+
    |                     |%%%%%%%%%%%%%%%%%%|                  |
    |                     |                  |                  | 
    |                     | b0               |b1                |
    |                     | outer PF         |lower outer sol   |
    |                     |                  |                  |
    0                     +------------------+------------------+

       0 -------------------------------------------------------- -> x

    Edges that touch coincide are physically connected unless
    otherwise indicated by a special symbol. Edges with a special
    symbol such as o,x,%, or % are instead connected to the other
    edge with the same symbol. Edges that do not coincide with
    another edge are a physical boundary.
  */  

  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, params->geqdsk_file, sizeof(params->geqdsk_file));

  // Theta limits are actually set internally by the code.
  double theta_min = -1.0, theta_max = 1.0;

  double psi_sep  = params->psi_sep; // Psi at the separatrix.
  double psi_axis = params->psi_axis; // Psi at the magnetic axis.
  double psi_min_core = params->psi_min_core; // Minimum psi the core.
  double psi_max_sol  = params->psi_max_sol ; // Maximum psi the SOL.
  double psi_min_pf   = params->psi_min_pf  ; // Minimum psi the private flux.

  // Number of cells.
  int Npsi_sol        = params->Npsi_sol       ;
  int Npsi_pf         = params->Npsi_pf        ;
  int Npsi_core       = params->Npsi_core      ;
  int Ntheta_divertor = params->Ntheta_divertor;
  int Ntheta_sol      = params->Ntheta_sol     ;

  // Block 0: outer private flux (PF) region.
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { psi_min_pf, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_pf, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_R,
          .rleft = 1.1,
          .rright = 1.7,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.25,
          .zmax = -0.9,
          .zmin_left = -1.25,
          .zmin_right = -1.25,
          .plate_spec = false,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
        { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  // Block 1: lower outer SOL.
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_LO,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmax = 1.0,
          .zmin = -1.25,
          .zmin_left = -1.25,
          .zmin_right = -1.25,
          .plate_spec = false,
          //.plate_func_lower = divertor_plate_func_out,
          //.plate_func_upper = divertor_plate_func_in,
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
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_sol },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_MID,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmax = 1.0,
          .zmin = -1.25,
          .zmin_left = -1.25,
          .zmin_right = -1.25,
          .plate_spec = false,
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
  gkyl_block_geom_set_block(bgeom, 3, &(struct gkyl_block_geom_info) {
      .lower = { psi_sep, theta_min },
      .upper = { psi_max_sol, theta_max },
      .cells = { Npsi_sol, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_LSN_SOL_UP,
          .rclose = 2.5,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax= 2.5,
          .zmax = 1.0,
          .zmin = -1.25,
          .zmin_left = -1.25,
          .zmin_right = -1.25,
          .plate_spec = false,
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

  // Block 4: inner private flux (PF) region.
  gkyl_block_geom_set_block(bgeom, 4, &(struct gkyl_block_geom_info) {
      .lower = { psi_min_pf, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_pf, Ntheta_divertor },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_PF_LO_L,
          .rclose = 1.1,
          .rleft = 1.1,
          .rright = 1.7,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmax = -0.9,
          .zmin = -1.25,
          .zmin_left = -1.25,
          .zmin_right = -1.25,
          .plate_spec = false,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // Physical boundary.
        { .bid = 4, .dir = 1, .edge = GKYL_PHYSICAL},
      }
    }
  );

  // Block 5: core region.
  gkyl_block_geom_set_block(bgeom, 5, &(struct gkyl_block_geom_info) {
      .lower = { psi_min_core, theta_min },
      .upper = { psi_sep, theta_max },
      .cells = { Npsi_core, Ntheta_sol },
      .cuts = { 1, 1 },
      .geometry = {
        .world = {0.0},
        .geometry_id = GKYL_TOKAMAK,
        .efit_info = efit_inp,
        .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
          .ftype = GKYL_CORE,
          .rclose = 2.0,
          .rleft = 0.8,
          .rright = 2.5,
          .rmin = 0.8,
          .rmax = 2.5,
          .zmin = -1.0,
          .zmax = 1.0,
        }
      },

      .connections[0] = { // x-direction.
        { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
      },
      .connections[1] = { // z-direction.
        { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
      }
    }
  );

  return bgeom;
}


struct gk_asdex_ctx
create_ctx(void)
{
  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double mi = 2.014*GKYL_PROTON_MASS; // Ion mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.
  double qi = eV; // Ion charge.
  double qe = -eV; // Electron charge.

  double Te = 150.0*eV; // Electron temperature.
  double Ti = 150.0*eV; // Ion temperature.
  double B0 = (1.937854e+00+3.937710e+00)/2.0; // B field amplitude.
  double n0 = 1.0e19; // Particle density.

  // Derived parameters.
  double vt_ion = sqrt(Ti/mi);
  double vt_elc = sqrt(Te/me);
  double c_s = sqrt(Te/mi);
  double omega_ci = fabs(qi*B0/mi);
  double rho_s = c_s/omega_ci;

  // Collision parameters.
  double nu_frac = 1.0;  
  double logLambda_elc = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Te/eV);
  double nu_elc = nu_frac*logLambda_elc*pow(eV, 4.0)*n0/(6.0*sqrt(2.0)*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(me)*pow(Te,3.0/2.0));

  double logLambda_ion = 6.6 - 0.5*log(n0/1e20) + 1.5*log(Ti/eV);
  double nu_ion = nu_frac*logLambda_ion*pow(eV, 4.0)*n0/(12.0*M_PI*sqrt(M_PI)*pow(eps0,2)*sqrt(mi)*pow(Ti,3.0/2.0));

  // Location of the numerical equilibrium.
  char geqdsk_file[128] = "./asdex.geqdsk";

  // Position space parameters.
  double num_blocks = 6;
  double R_axis = (1.61640+1.70022)/2.0; // R of the magnetic axis.
  double Z_axis = (-0.0013+0.1001)/2.0; // Z of the magnetic axis.
  double R_sep_OZA = 2.1389435; // Separatrix major at outboard Z axis.
  double R_sep_omp = 2.1334876; // Separatrix major at the OMP.
  double psi_axis = -9.276977e-02; // Psi at the magnetic axis.
  // Get the separatrix psi.
  struct gkyl_efit_inp efit_inp = {
    // psiRZ and related inputs
    .rz_poly_order = 2,
    .flux_poly_order = 1,
  };
  // Copy eqdsk file into efit_inp.
  memcpy(efit_inp.filepath, geqdsk_file, sizeof(geqdsk_file));
  struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
  double psi_sep = efit->psisep;
  double Rxpt = efit->Rxpt[0], Zxpt = efit->Zxpt[0];
  gkyl_efit_release(efit);
  double psi_min_core =  0.14;
  double psi_max_sol = 0.17;
  double psi_min_pf = 0.14;

  // Number of cells.
  int Npsi_sol = 4;
  int Npsi_pf = 4;
  int Npsi_core = 4;
  int Ntheta_divertor = 4;
  int Ntheta_sol = 24;

  printf("  X-point @ (R,Z) = (%.9e,%9e)\n",Rxpt,Zxpt);
  printf("  psi_axis = %.13e\n",psi_axis);
  printf("  psi_sep = %.13e\n",psi_sep);
  printf("  psi_min_core = %.13e\n",psi_min_core);
  printf("  psi_max_sol = %.13e\n",psi_max_sol);
  printf("  psi_min_pf = %.13e\n",psi_min_pf);
  printf("  Npsi_sol        = %d\n",Npsi_sol       );
  printf("  Npsi_pf         = %d\n",Npsi_pf        );
  printf("  Npsi_core       = %d\n",Npsi_core      );
  printf("  Ntheta_divertor = %d\n",Ntheta_divertor);
  printf("  Ntheta_sol      = %d\n",Ntheta_sol     );

  struct gk_asdex_ctx ctx = {
    .cdim = cdim,
    .B0 = B0, 
    .num_blocks = num_blocks,
    .psi_axis = psi_axis,
    .psi_sep = psi_sep,
    .psi_min_core = psi_min_core,
    .psi_max_sol = psi_max_sol,
    .psi_min_pf = psi_min_pf,
    .Npsi_sol        = Npsi_sol       ,
    .Npsi_pf         = Npsi_pf        ,
    .Npsi_core       = Npsi_core      ,
    .Ntheta_divertor = Ntheta_divertor,
    .Ntheta_sol      = Ntheta_sol     ,
  };

  // Copy eqdsk file into ctx.
  memcpy(ctx.geqdsk_file, geqdsk_file, sizeof(geqdsk_file));
  return ctx;
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

  struct gk_asdex_ctx ctx = create_ctx(); // Context for init functions.
                    
  // Construct block geometry
  struct gkyl_block_geom *bgeom = create_asdex_lsn_block_geom(&ctx);


  struct gkyl_gyrokinetic_multib app_inp = {
    .name = "gridmbasdex_2x",

    .cdim = ctx.cdim, .vdim = ctx.vdim,
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,

    .block_geom = bgeom,

    .comm = comm,
    .use_gpu = app_args.use_gpu,
  };

  // Create app object.
  struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new_geom(&app_inp);

 
  freeresources:
  // Free resources after simulation completion.
  gkyl_comm_release(comm);
  gkyl_block_geom_release(bgeom);
  gkyl_gyrokinetic_multib_app_release_geom(app);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  
  return 0;
}
