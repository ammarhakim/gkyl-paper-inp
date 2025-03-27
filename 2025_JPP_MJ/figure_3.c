#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>

// Make the context
struct ctx_mj_int {
  int j_T; // size of the box
};

// Make the context
struct ctx_mj_int_cold {
  double temp; // size of the box
};

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Helper functions to compute T scan on a log grid
double compute_T_log(int j, int j_max) {
    double T_min = 1e-5;
    double T_max = 1.0;
    return T_min * pow(T_max / T_min, (double)j / (j_max - 1));
}

// Helper functions to compute du scan on a log grid
double compute_du_log(int i, int i_max) {
    double du_min = 0.01;
    double du_max = 1.0;
    return du_min * pow(du_max / du_min, (double)i / (i_max - 1));
}

void 
eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void 
eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0; // Was 0.5 for the original test
}

void 
eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0; // Was 0.5 for the original test
  fout[1] = 0.0; // Was 0.5 for the original test
}

void 
eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0; // Was 0.5 for the original test
  fout[1] = 0.0; // Was 0.5 for the original test
  fout[2] = 0.0; // Was 0.5 for the original test
}


void 
eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  struct ctx_mj_int *app = ctx;
  int j_T = app->j_T;
  fout[0] = T + j_T;
}

void 
eval_M2_non_rel(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  struct ctx_mj_int_cold *app = ctx;
  double temp = app->temp;
  fout[0] = temp;
}

static bool
test_1x1v_rel(int poly_order, int i_vmax, int j_T)
{
  double vmax = 10.0 + i_vmax;
  double lower[] = {0.1, -vmax}, upper[] = {1.0, vmax};
  int cells[] = {2, 32}; 
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int ctx = {.j_T = j_T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, 1.0 + j_T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}



static bool
test_1x1v_non_rel(int poly_order, double du, double T, int i_vmax, int j_T)
{
  double vth = sqrt(T);
  double vmax = 10.0;
  int NV = 2*(vmax*vth)/du + 2;
  double eps_fac = du*NV/2.0 - vmax*vth - du;
  vmax = vmax*vth + du + eps_fac;
  double lower[] = {0.1, -vmax}, upper[] = {1.0, vmax};
  int cells[] = {2, NV}; 
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int_cold ctx = {.temp = T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_non_rel, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}


static bool
test_1x2v_rel(int poly_order, int i_vmax, int j_T)
{
  double vmax = 10.0 + i_vmax;
  double lower[] = {0.1, -vmax, -vmax}, upper[] = {1.0, vmax, vmax};
  int cells[] = {2, 32, 32}; 
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int ctx = {.j_T = j_T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, 1.0 + j_T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}



static bool
test_1x2v_non_rel(int poly_order, double du, double T, int i_vmax, int j_T)
{
  double vth = sqrt(T);
  double vmax = 10.0;
  int NV = 2*(vmax*vth)/du + 2;
  double eps_fac = du*NV/2.0 - vmax*vth - du;
  vmax = vmax*vth + du + eps_fac;
  double lower[] = {0.1, -vmax, -vmax}, upper[] = {1.0, vmax, vmax};
  int cells[] = {2, NV, NV}; 
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int_cold ctx = {.temp = T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_non_rel, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}

static bool
test_1x3v_rel(int poly_order, int i_vmax, int j_T)
{
  double vmax = 10.0 + i_vmax;
  double lower[] = {0.1, -vmax, -vmax, -vmax}, upper[] = {1.0, vmax, vmax, vmax};
  int cells[] = {1, 32, 32, 32}; 
  int vdim = 3, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int ctx = {.j_T = j_T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);
  
  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);


  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, 1.0 + j_T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  //char fname[1024];
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}



static bool
test_1x3v_non_rel(int poly_order, double du, double T, int i_vmax, int j_T)
{
  double vth = sqrt(T);
  double vmax = 5.0; //10.0;
  int NV = 2*(vmax*vth)/du + 2;
  double eps_fac = du*NV/2.0 - vmax*vth - du;
  vmax = vmax*vth + du + eps_fac;
  double lower[] = {0.1, -vmax, -vmax, -vmax}, upper[] = {1.0, vmax, vmax, vmax};
  int cells[] = {1, NV, NV, NV}; 
  int vdim = 3, cdim = 1;
  int ndim = cdim + vdim;

  if (NV > 200) {
    // Simulating how your output might look
    printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
    printf("dvx(%d) = %.16e;\n", i_vmax+1, du);
    printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, T);
    printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, -1);
    printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, 20);
    return -1;
  }

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext; 
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {0};
  struct gkyl_range confLocal, confLocal_ext; 
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext; 
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  struct ctx_mj_int_cold ctx = {.temp = T};
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_non_rel, &ctx);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  // build gamma and gamma_inv
  struct gkyl_array *gamma = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_array *gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  struct gkyl_dg_calc_sr_vars *sr_vars = gkyl_dg_calc_sr_vars_new(&grid, &vel_grid,
      &confBasis,  &velBasis, &confLocal, &velLocal, false);
  // Project gamma and its inverse
  gkyl_calc_sr_vars_init_p_vars(sr_vars, gamma, gamma_inv);
  // Free SR variable computation
  gkyl_dg_calc_sr_vars_release(sr_vars);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
    .phase_grid = &grid,
    .vel_grid = &vel_grid, 
    .conf_basis = &confBasis,
    .vel_basis = &velBasis, 
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .phase_range = &local,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .model_id = GKYL_MODEL_SR,
    .use_gpu = false,
    .max_iter = 20,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  //char fname[1024];
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Simulating how your output might look
  printf("vmax(%d) = %.16e;\n", i_vmax+1, vmax);
  printf("dvx(%d) = %.16e;\n", i_vmax+1, vel_grid.dx[0]);
  printf("T(%d,%d) = %.16e;\n",i_vmax+1, j_T+1, T);
  printf("status(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.iter_converged);
  printf("niter(%d,%d) = %d;\n",i_vmax+1, j_T+1, status_corr.num_iter);

  // release memory for moment data object
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);

  return status_corr.iter_converged;
}

// special note, the p1 basis does not function
void test_1x1v_rel_limits_p2() { 
  bool stat;
  printf("function [T,vmax,niter,status] = load_data_rel_1x1v()\n");
  for (int i=0; i<300; ++i){ //vmax, 300 originally
    stat = 0;
    for (int j=0; j<50; ++j){ //T, 50 originally
      stat = test_1x1v_rel(2,i,j); 
    }
  }
  printf("end\n");
}


// special note, the p1 basis does not function
// Test du range: 0.01 to 1.0
// Test temp range: 1.0 to 1e-5
void test_1x1v_non_rel_limits_p2() {
  int i_max = 50;
  int j_max = 50;
  printf("function [T,dvx,niter,status] = load_data_non_rel_1x1v()\n");
  for (int i = 0; i < i_max; ++i) { // du index
    for (int j = 0; j < j_max; ++j) { // T index
      double du = compute_du_log(i, i_max);
      double T = compute_T_log(j, j_max);
      test_1x1v_non_rel(2, du, T, i, j);
    }
  }
  printf("end\n");
}

// special note, the p1 basis does not function
void test_1x2v_rel_limits_p2() { 
  printf("function [T,vmax,niter,status] = load_data_rel_1x2v()\n");
  for (int i=0; i<300; ++i){ //vmax, 300 originally
    for (int j=0; j<50; ++j){ //T, 50 originally
      test_1x2v_rel(2,i,j); 
    }
  }
  printf("end\n");
}


// special note, the p1 basis does not function
// Test du range: 0.01 to 1.0
// Test temp range: 1.0 to 1e-5
void test_1x2v_non_rel_limits_p2() {
  int i_max = 50;
  int j_max = 50;
  printf("function [T,dvx,niter,status] = load_data_non_rel_1x2v()\n");
  for (int i = 0; i < i_max; ++i) { // du index
    for (int j = 0; j < j_max; ++j) { // T index
      double du = compute_du_log(i, i_max);
      double T = compute_T_log(j, j_max);
      test_1x2v_non_rel(2, du, T, i, j);
    }
  }
  printf("end\n");
}


// special note, the p1 basis does not function
void test_1x3v_rel_limits_p2() { 
  printf("function [T,vmax,niter,status] = load_data_rel_1x3v()\n");
  for (int i=0; i<300; ++i){ //vmax, 300 originally
   for (int j=0; j<50; ++j){ //T, 50 originally
     test_1x3v_rel(2,i,j); 
   }
  }
  //test_1x3v_rel(2,150,10); 
  printf("end\n");
}


// special note, the p1 basis does not function
// Test du range: 0.01 to 1.0
// Test temp range: 1.0 to 1e-5
void test_1x3v_non_rel_limits_p2() {
  int i_max = 50;
  int j_max = 50;
  printf("function [T,dvx,niter,status] = load_data_non_rel_1x3v()\n");
  for (int i = 0; i < i_max; ++i) { // du index
    for (int j = 0; j < j_max; ++j) { // T index
      double du = compute_du_log(i, i_max);
      double T = compute_T_log(j, j_max);
      test_1x3v_non_rel(2, du, T, i, j);
    }
  }
  printf("end\n");
}

TEST_LIST = {
  {"test_1x1v_rel_limits_p2", test_1x1v_rel_limits_p2},
  {"test_1x1v_non_rel_limits_p2", test_1x1v_non_rel_limits_p2},
  {"test_1x2v_rel_limits_p2", test_1x2v_rel_limits_p2},
  {"test_1x2v_non_rel_limits_p2", test_1x2v_non_rel_limits_p2},
  {"test_1x3v_rel_limits_p2", test_1x3v_rel_limits_p2},
  {"test_1x3v_non_rel_limits_p2", test_1x3v_non_rel_limits_p2},
  {NULL, NULL},
};
