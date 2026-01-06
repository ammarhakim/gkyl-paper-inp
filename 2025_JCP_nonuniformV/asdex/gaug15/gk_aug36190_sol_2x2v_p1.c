#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
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
  double diffusionD; // Diffusion coefficient.

  // Source parameters
  double psi_src; // Location.
  double lambda_src; // Width.
  double ndot_src; // Particle source rate.
  double Te_src; // Electron source temperature.
  double Ti_src; // Ion source temperature.

  // Domain parameters.            
  double psi_min, psi_max; // Psi extents.
  double Lx; // Box size in x
  double Lz; // Box size in z
  // Z location of the X-point on psi=psi_min.
  double z_xpt_psi_min_lo, z_xpt_psi_min_up;

  // Grid.
  int Nx, Nz; // Number of cells in x,z.
  int Nvpar, Nmu; // Number of cells in vpar,mu.
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.

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

  double t_end; // End time.
  int num_frames; // Number of output frames.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.
};

void divertor_plate_func_out(double s, double* RZ)
{
//  // Straight plate from (1.5966,-1.1421) to (1.6888,-0.8781).
//  RZ[0] = 1.5966 + (1.6888 - 1.5966)*s;
//  RZ[1] = -1.1421 + (-0.8781 - (-1.1421))*s;
  // Cubic spline approximation to the plate.
  const int npts = 100;
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        
  };
  double R[] = {
    1.58211701, 1.58373013, 1.58522876, 1.58662408, 1.58792727, 1.58914951,
    1.59030198, 1.59139587, 1.59244236, 1.59345262, 1.59443785, 1.59540921,
    1.59637789, 1.59735404, 1.59833796, 1.59932451, 1.60030849, 1.60128472,
    1.60224816, 1.60320156, 1.60415862, 1.60513384, 1.6061417 , 1.60719629,
    1.608296  , 1.6094194 , 1.61054379, 1.61164649, 1.61270478, 1.61369639,
    1.61461937, 1.61550099, 1.6163708 , 1.61725834, 1.61819317, 1.6192048 ,
    1.62030643, 1.6214671 , 1.62264849, 1.62381228, 1.62492016, 1.62595178,
    1.62693809, 1.62791909, 1.62893481, 1.63002528, 1.63123012, 1.63256489,
    1.63400806, 1.63553491, 1.63712072, 1.63874077, 1.64037203, 1.642006  ,
    1.64364166, 1.64527805, 1.64691421, 1.64854896, 1.65017582, 1.65178259,
    1.65335681, 1.65488603, 1.65635783, 1.65777956, 1.65920642, 1.66070072,
    1.66232475, 1.66412714, 1.66609031, 1.66817673, 1.67034882, 1.67256901,
    1.67480116, 1.67704013, 1.67931122, 1.681641  , 1.68405605, 1.68658294,
    1.68924066, 1.69201572, 1.69488593, 1.69782911, 1.70082309, 1.70384595,
    1.70688438, 1.70993515, 1.7129956 , 1.71606307, 1.71913489, 1.72220864,
    1.72528362, 1.72836005, 1.73143813, 1.73451804, 1.7376    , 1.7406842 ,
    1.74377084, 1.74686013, 1.74995225, 1.75304743
  };
  double Z[] = {
    -1.1797302 , -1.17584427, -1.17192351, -1.16797136, -1.16399128, -1.15998669,
    -1.15596104, -1.15191777, -1.14786031, -1.14379211, -1.13971661, -1.13563724,
    -1.13155745, -1.12748036, -1.12340608, -1.11933308, -1.11525979, -1.11118464,
    -1.10710612, -1.10302482, -1.09894425, -1.09486817, -1.09080032, -1.08674434,
    -1.08269988, -1.0786615 , -1.07462347, -1.07058004, -1.06652547, -1.06245411,
    -1.05836552, -1.05426666, -1.05016509, -1.04606836, -1.04198404, -1.03791967,
    -1.03387851, -1.02985215, -1.02583029, -1.02180258, -1.01775872, -1.01369345,
    -1.00961605, -1.00553836, -1.00147223, -0.99742949, -0.99342188, -0.98945418,
    -0.9855205 , -0.981614  , -0.97772786, -0.97385526, -0.96998977, -0.96612844,
    -0.96227006, -0.95841347, -0.95455747, -0.95070083, -0.94684033, -0.94297071,
    -0.93908658, -0.93518255, -0.93125327, -0.92730239, -0.92335531, -0.91944071,
    -0.91558724, -0.91181723, -0.90812254, -0.90448578, -0.90088956, -0.8973165 ,
    -0.89375003, -0.8901916 , -0.88666037, -0.88317622, -0.87975903, -0.87642869,
    -0.87319966, -0.87006321, -0.86700441, -0.8640083 , -0.86105996, -0.85814459,
    -0.85525204, -0.85237755, -0.84951673, -0.84666513, -0.84381835, -0.84097233,
    -0.83812611, -0.83528023, -0.83243524, -0.82959169, -0.82675013, -0.82391111,
    -0.82107517, -0.81824287, -0.81541475, -0.81259136,
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

void divertor_plate_func_in(double s, double* RZ)
{
//  // Straight plate from (1.2686,-1.0520) to (1.1886,-0.7294).
//  RZ[0] = 1.2686 + (1.1886 - 1.2686)*s;
//  RZ[1] = -1.0520 + (-0.7294 - (-1.0520))*s;
//  printf("s = %g\n",s);
  // Cubic spline approximation to the plate.
  const int npts = 100;
  double t[] = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        ,
  };
  double R[] = {
    1.24928187, 1.25037791, 1.25164662, 1.25307592, 1.25465372, 1.25636793,
    1.25820646, 1.26015722, 1.26220783, 1.26432432, 1.26643709, 1.2684732 ,
    1.27035972, 1.27204779, 1.27356159, 1.27493914, 1.27621846, 1.27743755,
    1.27861728, 1.27972182, 1.28070361, 1.28151512, 1.2821169 , 1.28252687,
    1.28278764, 1.28294196, 1.28303254, 1.28308818, 1.28309865, 1.2830469 ,
    1.28291589, 1.28268856, 1.28234789, 1.28187859, 1.2812746 , 1.28053274,
    1.27964986, 1.2786228 , 1.27744867, 1.276135  , 1.27470289, 1.27317436,
    1.27157144, 1.26991054, 1.26818273, 1.26637194, 1.26446206, 1.26243818,
    1.26029643, 1.25803897, 1.25566806, 1.2531884 , 1.25063723, 1.24807477,
    1.24556173, 1.2431587 , 1.24088867, 1.23868208, 1.23645539, 1.23412506,
    1.23161335, 1.22890566, 1.22602622, 1.22299979, 1.21985393, 1.21665898,
    1.21351899, 1.21053889, 1.20780586, 1.20524866, 1.20271309, 1.20004418,
    1.19709218, 1.19385834, 1.1905151 , 1.18724412, 1.18422655, 1.18154187,
    1.17904582, 1.1765643 , 1.17392753, 1.17107626, 1.16806833, 1.16496727,
    1.16183612, 1.15872996, 1.15569703, 1.15278535, 1.15004111, 1.14745534,
    1.14495522, 1.14246444, 1.13992525, 1.13736378, 1.13482978, 1.132373  ,
    1.13004319, 1.12789008, 1.12596345, 1.12431302,
  };
  double Z[] = {
   -1.08698226, -1.08196971, -1.07700781, -1.07209476, -1.06722871, -1.06240785,
   -1.05763035, -1.05289438, -1.048198  , -1.0435299 , -1.03886333, -1.0341701 ,
   -1.02942201, -1.02460028, -1.01971465, -1.01478027, -1.00981228, -1.00482584,
   -0.99983182, -0.99482694, -0.98980503, -0.98475992, -0.97968661, -0.97458858,
   -0.96947299, -0.96434696, -0.95921765, -0.95409065, -0.94896713, -0.9438475 ,
   -0.93873217, -0.93362155, -0.92851607, -0.92341684, -0.91832872, -0.91325773,
   -0.90820991, -0.90319127, -0.89820775, -0.89326079, -0.88834607, -0.88345886,
   -0.87859442, -0.87374979, -0.86893001, -0.86414237, -0.85939418, -0.85469247,
   -0.85004179, -0.84544538, -0.84090641, -0.83642654, -0.83198715, -0.82755527,
   -0.82309763, -0.81858101, -0.81399271, -0.80937046, -0.80475963, -0.80020557,
   -0.79575083, -0.79140744, -0.78716867, -0.78302749, -0.77897453, -0.77496393,
   -0.77092106, -0.76677053, -0.76244948, -0.75800658, -0.75354893, -0.74918416,
   -0.7450163 , -0.7410448 , -0.73715056, -0.73320813, -0.72909236, -0.72474766,
   -0.72027146, -0.71578157, -0.71139291, -0.70714603, -0.70300264, -0.69892065,
   -0.69485832, -0.69077974, -0.68665405, -0.68245053, -0.67813961, -0.67372661,
   -0.66925716, -0.66477913, -0.66032943, -0.6558955 , -0.65145088, -0.6469691 ,
   -0.64242369, -0.63778818, -0.63303609, -0.62814097,
  };
  // Find indices in t that bound s.
  int idx_tlo, idx_tup;
  if (s < 1e-8) {
    idx_tlo = 0;
    idx_tup = 0;
  }
  else if (fabs(s-1.0) < 1e-8) {
    idx_tlo = npts-1;
    idx_tup = npts-1;
  }
  else {
    for (int i=0; i<npts-1; i++) {
      if (t[i] <= s && s < t[i+1]) {
        idx_tlo = i;
        idx_tup = i+1;
        break;
      }
    }
  }
  // Interpolate the value of R and Z.
  double Dt = t[idx_tup]-t[idx_tlo];
  if (idx_tlo == idx_tup) {
    RZ[0] = R[idx_tlo];
    RZ[1] = Z[idx_tlo];
  }
  else {
    RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo];
    RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo];
  }
}

double rho_psi(double psi, double psi_axis, double psi_sep)
{
  // Normalized radial coordinate.
  return sqrt((psi-psi_axis) / (psi_sep - psi_axis));
}

double psi_rho(double rho, double psi_axis, double psi_sep)
{
  // Poloidal flux given the normalized radial coordinate.
  return pow(rho,2) * (psi_sep - psi_axis) + psi_axis;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double psi_min = app->psi_min;

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 6 experimental.
  double den0 = 0.7098096769979423e19;
  double lambda_psi = 0.009;
  fout[0] = den0 * exp(-(x - psi_min)/lambda_psi);
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double psi_min = app->psi_min;
  double eV = GKYL_ELEMENTARY_CHARGE;

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 7 experimental.
  double Te0 = 59.52408228927226*eV;
  double lambda_psi = 0.013;
  fout[0] = Te0 * exp(-(x - psi_min)/lambda_psi);
}

void
eval_temp_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double psi_min = app->psi_min;
  double eV = GKYL_ELEMENTARY_CHARGE;

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 7 experimental.
  double Ti0 = 64.57831221968155*eV;
  double lambda_psi = 0.01;
  fout[0] = Ti0 * exp(-(x - psi_min)/lambda_psi);
}

double main_chamber_parallel_mask(double z, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double z_xpt_psi_min_lo = app->z_xpt_psi_min_lo;

  double mu = 0.0; //z_xpt_psi_min_lo/2.0;
  double sig = 2.6;
  int p = 24;
  double floor = 0.001;

  return fmax( exp(-0.5*pow((z - mu)/sig,p)), floor );
}

void
eval_density_bcx_lo(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double psi_min = app->psi_min;

  // Mimics the SOL profile in D. Michels, et al. Phys. Plasmas 29, 032307
  // (2022), figure 6 experimental.
  double den0 = 0.7098096769979423e19;
  double lambda_psi = 0.009;
  fout[0] = den0 * exp(-(x - psi_min)/lambda_psi) * main_chamber_parallel_mask(z, ctx);
}

double
parallel_profile_source(double z, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double z_xpt_psi_min_lo = app->z_xpt_psi_min_lo;

  double mu = z_xpt_psi_min_lo/2.0;
  double sig = 1.0;
  int p = 6;
  double floor = 0.001;

  return fmax( exp(-0.5*pow((z - mu)/sig,p)), floor );
}

void
eval_density_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = app->lambda_src;
  double psi_src = app->psi_src;
  double ndot_src = app->ndot_src;

  double source_floor = 1e-10;
  if (x < psi_src + 3*lambda_src)
    source_floor = 1e-2;

  double src_prof = exp(-pow(x-psi_src,2)/(2*pow(lambda_src,2))) * parallel_profile_source(z, ctx);
  fout[0] = ndot_src * fmax(src_prof, source_floor);
}

void
eval_upar_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = app->lambda_src;
  double psi_src = app->psi_src;
  double Te_src = app->Te_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (x < psi_src + 3*lambda_src)
    fout[0] = Te_src * parallel_profile_source(z, ctx);
  else
    fout[0] = 2.0*eV;

}

void
eval_temp_ion_source(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double x = xn[0], z = xn[1];

  double lambda_src = app->lambda_src;
  double psi_src = app->psi_src;
  double Ti_src = app->Ti_src;
  double eV = GKYL_ELEMENTARY_CHARGE;

  if (x < psi_src + 3*lambda_src)
    fout[0] = Ti_src * parallel_profile_source(z, ctx);
  else
    fout[0] = 2.0*eV;
}

void
eval_nu_elc(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  fout[0] = app->nu_elc;
}

void
eval_nu_ion(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  fout[0] = app->nu_ion;
}

void
diffusion_D_func(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  struct gk_asdex_ctx *app = ctx;

  fout[0] = app->diffusionD; // Diffusivity [m^2/s].
}

// Velocity space mappings.
void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_asdex_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

struct gk_asdex_ctx
create_ctx(struct gkyl_app_args *app_args)
{
  // Reference simulation was run with.
  //  Nx = 16
  //  Nz = 16
  //  Nvpar = 16
  //  Nmu = 8
  int Nx; // Number of cells in x.
  int Nz; // Number of cells in z, originally 8.
  int Nvpar; // Number of cells in vpar.
  int Nmu; // Number of cells in mu.

  // Extract variables from command line arguments.
  sscanf(app_args->opt_args, "Nx=%d,Nz=%d,Nvpar=%d,Nmu=%d", &Nx, &Nz, &Nvpar, &Nmu);
  printf("Command line arguments:\n");
  printf("  Nx    = %d\n", Nx   );
  printf("  Nz    = %d\n", Nz   );
  printf("  Nvpar = %d\n", Nvpar);
  printf("  Nmu   = %d\n", Nmu  );

  int cdim = 2, vdim = 2; // Dimensionality.

  double eps0 = GKYL_EPSILON0; // Permittivity of free space.
  double eV = GKYL_ELEMENTARY_CHARGE; // Elementary charge.
  double mi = 2.014*GKYL_PROTON_MASS; // Ion mass.
  double me = GKYL_ELECTRON_MASS; // Electron mass.
  double qi = eV; // Ion charge.
  double qe = -eV; // Electron charge.

  double Te = 37.5*eV; // Electron temperature.
  double Ti = 38.0*eV; // Ion temperature.
  double B0 = (1.937830e+00+3.930574e+00)/2.0; // B field amplitude.
  double n0 = 0.4e19; // Particle density.

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

  double diffusionD = 0.3; // Particle diffusion coefficient.

  // Position space parameters.
  double R_axis = (1.61640+1.70022)/2.0; // R of the magnetic axis.
  double Z_axis = (-0.0013+0.1001)/2.0; // Z of the magnetic axis.
  double R_sep_OZA = 2.1389435; // Separatrix major at outboard Z axis.
  double R_sep_omp = 2.1334876; // Separatrix major at the OMP.
  double psi_sep = 0.1497542827844; // Psi of the separatrix.
  double psi_axis = -9.276977e-02; // Psi at the magnetic axis.
  // Here rho = sqrt((psi-psi_axis) / (psi_sep - psi_axis)).
  double rho_min = 1.005;
  double rho_max = 1.04;
  double psi_min = psi_rho(rho_min, psi_axis, psi_sep);
  double psi_max = psi_rho(rho_max, psi_axis, psi_sep);
  double Lx = psi_max - psi_min;
  double Lz = (M_PI-1.0e-14)*2;
  // z location of the X-point on psi=psi_min.
  double z_xpt_psi_min_lo = -2.8469;
  double z_xpt_psi_min_up =  2.8486;
  printf("  psi_axis = %.13e\n",psi_axis);
  printf("  psi_sep = %.13e\n",psi_sep);
  printf("  psi_min = %.13e\n",psi_min);
  printf("  psi_max = %.13e\n",psi_max);

  // Source parameters.
  double psi_src = psi_rho(rho_min, psi_axis, psi_sep);
  double Lc_src = 67.0; // Connection length where source is localized.
  double lambda_src = 0.14*Lx;
  double ndot_src = 1.2*2.0*0.7e19*c_s/(Lc_src/4.0);
  double Te_src = 2.0*1.15*1.3*1.5*Te;
  double Ti_src = 2.0*1.15*1.3*1.5*Ti;

  // Physical velocity space limits
  double vpar_max_elc = 6.0*vt_elc;
  double mu_max_elc = me*pow(4.0*vt_elc,2)/(2.0*B0);

  double vpar_max_ion = 6.0*vt_ion;
  double mu_max_ion = mi*pow(4.0*vt_ion,2)/(2.0*B0);

  // Computational velocity space limits.
  double vpar_min_ion_c = -1.0/sqrt(2.0);
  double vpar_max_ion_c = 1.0/sqrt(2.0);
  double mu_min_ion_c = 0.;
  double mu_max_ion_c = 1.;
  // Computational velocity space limits.
  double vpar_min_elc_c = -1.0/sqrt(2.0);
  double vpar_max_elc_c = 1.0/sqrt(2.0);
  double mu_min_elc_c = 0.;
  double mu_max_elc_c = 1.;

  double t_end = 2*4*5*100.0e-6;
  double num_frames = 2*4*5*50;
  double write_phase_freq = 0.2; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  struct gk_asdex_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .charge_elc = qe, 
    .charge_ion = qi, 
    .mass_elc = me, 
    .mass_ion = mi,
    .Te = Te, 
    .Ti = Ti, 
    .c_s = c_s, 
    .nu_elc = nu_elc, 
    .nu_ion = nu_ion, 
    .B0 = B0, 
    .n0 = n0, 
    .diffusionD = diffusionD,
    .psi_min = psi_min,
    .psi_max = psi_max,
    .Lx = Lx, 
    .Lz = Lz, 
    .z_xpt_psi_min_lo = z_xpt_psi_min_lo,
    .z_xpt_psi_min_up = z_xpt_psi_min_up,
    .lambda_src = lambda_src,
    .psi_src = psi_src,
    .ndot_src = ndot_src,
    .Te_src = Te_src,
    .Ti_src = Ti_src,
    // Physical velocity space limits
    .vpar_max_elc = vpar_max_elc, 
    .mu_max_elc = mu_max_elc, 
    .vpar_max_ion = vpar_max_ion, 
    .mu_max_ion = mu_max_ion, 
    // Computational velocity space limits
    .vpar_min_elc_c = vpar_min_elc_c,
    .vpar_max_elc_c = vpar_max_elc_c,
    .mu_min_elc_c = mu_min_elc_c,
    .mu_max_elc_c = mu_max_elc_c,
    .vpar_min_ion_c = vpar_min_ion_c,
    .vpar_max_ion_c = vpar_max_ion_c,
    .mu_min_ion_c = mu_min_ion_c,
    .mu_max_ion_c = mu_max_ion_c,
    .Nx = Nx,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Nz, Nvpar, Nmu},
    .t_end = t_end, 
    .num_frames = num_frames, 
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_field_energy(app);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
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

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  struct gk_asdex_ctx ctx = create_ctx(&app_args); // context for init functions

  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);

  // Electrons.
  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_elc,      
  };

  struct gkyl_gyrokinetic_projection elc_bcx_lo = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density_bcx_lo,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_elc, 
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.charge_elc, .mass = ctx.mass_elc,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_elc_c, ctx.mu_min_elc_c},
    .upper = { ctx.vpar_max_elc_c, ctx.mu_max_elc_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },

    .projection = elc_ic,
  
    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Te, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar = eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2 },
      }
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },
      
    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = elc_bcx_lo, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1,  GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN, },
    .boundary_flux_diagnostics = {
      .num_diag_moments = 3,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP, },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Ions.
  struct gkyl_gyrokinetic_projection ion_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_ion,      
  };

  struct gkyl_gyrokinetic_projection ion_bcx_lo = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
    .ctx_density = &ctx,
    .density = eval_density_bcx_lo,
    .ctx_upar = &ctx,
    .upar = eval_upar,
    .ctx_temp = &ctx,
    .temp = eval_temp_ion,      
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.charge_ion, .mass = ctx.mass_ion,
    .vdim = ctx.vdim,
    .lower = { ctx.vpar_min_ion_c, ctx.mu_min_ion_c},
    .upper = { ctx.vpar_max_ion_c, ctx.mu_max_ion_c},
    .cells = { cells_v[0], cells_v[1] },
    .polarization_density = ctx.n0,

    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .projection = ion_ic,

    .collisionless = {
      .type = GKYL_GK_COLLISIONLESS_ES,
    },

    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .den_ref = ctx.n0, // Density used to calculate coulomb logarithm
      .temp_ref = ctx.Ti, // Temperature used to calculate coulomb logarithm
      .num_cross_collisions = 1,
      .collide_with = { "elc" },
    },

    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_source,
        .ctx_upar = &ctx,
        .upar = eval_upar_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
      .diagnostics = {
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2 },
      }
    },

    .anomalous_diffusion = {
      .anomalous_diff_id = GKYL_GK_ANOMALOUS_DIFF_D,
      .D_profile = diffusion_D_func,
      .D_profile_ctx = &ctx,
//      .write_diagnostics = true,
    },
    
    .bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_FIXED_FUNC, .projection = ion_bcx_lo, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_ABSORB, },
      { .dir = 1, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
      { .dir = 1, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_SPECIES_SHEATH, },
    },

    .num_diag_moments = 4,
    .diag_moments = { GKYL_F_MOMENT_M1,  GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_BIMAXWELLIAN, },
    .boundary_flux_diagnostics = {
      .num_diag_moments = 3,
      .diag_moments = { GKYL_F_MOMENT_HAMILTONIAN, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP, },
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
//      .time_integrated = true,
    },
  };

  // Field
  struct gkyl_gyrokinetic_field field = {
    .poisson_bcs = {
      { .dir = 0, .edge = GKYL_LOWER_EDGE, .type = GKYL_BC_GK_FIELD_NEUMANN, .value = {0.0}, },
      { .dir = 0, .edge = GKYL_UPPER_EDGE, .type = GKYL_BC_GK_FIELD_DIRICHLET, .value = {0.0}, },
    },
  };

  struct gkyl_efit_inp efit_inp = {
    // Numerical equilibrium inputs.
    .filepath = "/pscratch/sd/m/mana/gkeyll/nonuniv_post_quadgeo/asdex/aug36190.geqdsk", // Equilibrium to use.
    .rz_poly_order = 2,              // Polynomial order for psi(R,Z) used for field line tracing.
    .flux_poly_order = 1,            // Polynomial order for fpol(psi).
  };

  struct gkyl_tok_geo_grid_inp grid_inp = {
    .ftype = GKYL_LSN_SOL, // Type of geometry.
    .rclose = 2.5,           // Closest R to region of interest.
    .rright = 2.5,           // Closest R to outboard SOL.
    .rleft = 0.7,            // Closest R to inboard SOL.
    .rmax = 2.5,             // Largest R in machine.
    .rmin = 0.7,             // Smallest R in machine.
    .zmin = -1.3,            // Lower Z boundary.
    .zmax = 1.0,             // Upper Z boundary.
    .zmin_left = -1.15,       // Z of inboard divertor plate.
    .zmin_right = -1.15,      // Z of outboard divertor plate.
    .plate_spec = true,
    .plate_func_lower = divertor_plate_func_out,
    .plate_func_upper = divertor_plate_func_in,
  };

  // GK app
  struct gkyl_gk gk = {
    .name = "gk_aug36190_sol_2x2v_p1",

    .cdim = ctx.cdim,
    .lower = { ctx.psi_min, -ctx.Lz/2.0 },
    .upper = { ctx.psi_max,  ctx.Lz/2.0 },
    .cells = { cells_x[0], cells_x[1] },
    .poly_order = 1,
    .basis_type = app_args.basis_type,
    .cfl_frac = 1.0,
//    .cfl_frac_omegaH = 1.0e9,

//    .enforce_positivity = true,

    .geometry = {
      .geometry_id = GKYL_TOKAMAK,
      .efit_info = efit_inp,
      .tok_grid_info = grid_inp,
    },

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 2,
    .species = { elc, ion },
    .field = field,

    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0], app_args.cuts[1] },
      .comm = comm,
    },
  };

  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }  

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false);

  // start, end and initial time-step
  double dt = t_end-t_curr;
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    if (step % 100 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", t_curr);

    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    if (step % 100 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);
  
  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
