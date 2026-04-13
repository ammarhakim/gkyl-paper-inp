#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// State of the pseudo orbit-averaged integrator.
enum gk_poa_state {
  GK_POA_NONE = 0, // Haven't started.
  GK_POA_OAP, // Orbit averaged phase.
  GK_POA_FDP, // Full dynamics phase.
  GK_POA_COMPLETED, // Finished simulation.
};

struct gk_poa_phase_params {
  enum gk_poa_state phase; // Type of phase.
  int num_frames; // Number of frames.
  double duration; // Duration.
  double alpha, alpha_ion, alpha_elc; // Factor multiplying collisionless terms.
  bool is_static_field; // Whether to evolve the field.
  bool is_positivity_enabled, is_positivity_enabled_ion, is_positivity_enabled_elc; // Whether positivity is enabled.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type, fdot_mult_type_ion, fdot_mult_type_elc; // Type of df/dt multipler.
  double f_threshold, f_threshold_ion, f_threshold_elc; // Threshold for df/dt multiplier.
  double cfl_factor_times_omega_max, cfl_factor_times_omega_max_ion, cfl_factor_times_omega_max_elc; // CFL factor for fixed factor times omega max multiplier.
  double time_dilation_scale_const, time_dilation_scale_const_ion, time_dilation_scale_const_elc; // Constant time dilation scale factor (if fdot_mult_type is GKYL_GK_FDOT_MULTIPLIER_CONSTANT).
  enum gkyl_gyrokinetic_damping_type damping_type; // Type of damping to apply.
  double damping_rate_const;
};

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p; // Hardcoded magnetic field at midplane 0.53 for double lorentzian
  double Bmag_midp; // Magnetic field magnitude at midplane (Z=0)
  double beta;
  double tau;
  double Ti0;
  double nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  double vti, vte;
  double RatZeq0; // Radius of the field line at Z=0.
  double kperp; // Perpendicular wavenumber for ES gyrokinetics.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_eval;
  double psi_max;
  double psi_min;
  double theta_eval;
  double theta_min;
  double theta_max;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double mu_max_ion;
  double vpar_max_elc;
  double mu_max_elc;

  int Npsi;
  int Ntheta;
  int Nz;
  int Nvpar;
  int Nmu;
  int Nvpar_elc;
  int Nmu_elc;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;

  // Source parameters
  double ion_source_amplitude;
  double ion_source_sigma;
  double ion_source_temp;

  double t_end; // End time.
  int num_frames; // Number of output frames.
  int num_phases; // Number of phases.
  struct gk_poa_phase_params *poa_phases; // Phases to run.
  double write_phase_freq; // Frequency of writing phase-space diagnostics (as a fraction of num_frames).
  double int_diag_calc_freq; // Frequency of calculating integrated diagnostics (as a factor of num_frames).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Geometry parameters for Lorentzian mirror
  double mcB;     // Magnetic field parameter
  double gamma;   // Width parameter for Lorentzian profile
  double Z_m;     // Mirror throat location
  double Z_min;   // Minimum Z coordinate
  double Z_max;   // Maximum Z coordinate
  double psi_in;  // Working variable for psi integration
  double z_in;    // Working variable for z integration
};

static inline const char*
fdot_multiplier_type_to_str(enum gkyl_gyrokinetic_fdot_multiplier_type type)
{
  switch (type) {
    case GKYL_GK_FDOT_MULTIPLIER_NONE:
      return "None";
    case GKYL_GK_FDOT_MULTIPLIER_USER_INPUT:
      return "User input";
    case GKYL_GK_FDOT_MULTIPLIER_LOSS_CONE:
      return "Loss cone";
    case GKYL_GK_FDOT_MULTIPLIER_CONSTANT:
      return "Constant";
    case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT:
      return "Fixed dt";
    case GKYL_GK_FDOT_MULTIPLIER_FIXED_FACTOR_TIMES_OMEGA_MAX:
      return "Fixed factor times omega max";
    case GKYL_GK_FDOT_MULTIPLIER_FIXED_DT_OMEGAH:
      return "Fixed dt omega_H";
    case GKYL_GK_FDOT_MULTIPLIER_DT_SET_BY_SPECIES:
      return "dt set by species";
    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_THRESHOLD:
      return "Mask f threshold";
    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_LOCAL:
      return "Mask f frac local";
    case GKYL_GK_FDOT_MULTIPLIER_MASK_F_FRAC_GLOBAL:
      return "Mask f frac global";
    default:
      return "Unknown";
  }
}

void
print_ctx(struct gk_mirror_ctx *ctx)
{
  printf("Plasma parameters:\n");
  printf("  Ion mass (mi) = %g\n", ctx->mi);
  printf("  Ion charge (qi) = %g\n", ctx->qi);
  printf("  Electron mass (me) = %g\n", ctx->me);
  printf("  Electron charge (qe) = %g\n", ctx->qe);
  printf("  Electron temperature (Te0) = %g eV\n", ctx->Te0/GKYL_ELEMENTARY_CHARGE);
  printf("  Ion temperature (Ti0) = %g eV\n", ctx->Ti0/GKYL_ELEMENTARY_CHARGE);
  printf("  Density (n0) = %g m^-3\n", ctx->n0);
  printf("  Ion thermal speed (vti) = %g m/s\n", ctx->vti);
  printf("  Electron thermal speed (vte) = %g m/s\n", ctx->vte);
  printf("  Magnetic field (B_p) = %g T\n", ctx->B_p);
  printf("  Beta = %g\n", ctx->beta);
  printf("  Tau = Ti0/Te0 = %g\n", ctx->tau);
  printf("  Ion-ion collision frequency factor (nuFrac) = %g\n", ctx->nuFrac);
  printf("  Ion-ion Coulomb logarithm (logLambdaIon) = %g\n", ctx->logLambdaIon);
  printf("  Ion-ion collision frequency (nuIon) = %g Hz\n", ctx->nuIon);
  printf("  Magnetic field magnitude at midplane (Bmag_midp) = %g T\n", ctx->Bmag_midp);
  printf("  Perpendicular wavenumber (kperp) = %g 1/m\n", ctx->kperp);
  
  printf("\nGeometry parameters:\n");
  printf("  Mirror throat radius (RatZeq0) = %g m\n", ctx->RatZeq0);
  printf("  Psi evaluated (psi_eval) = %g Wb\n", ctx->psi_eval);
  printf("  Psi extents: [%g, %g] Wb\n", ctx->psi_min, ctx->psi_max);
  printf("  Z extents: [%g, %g] m\n", ctx->Z_min, ctx->Z_max);
  printf("  z extents: [%g, %g] m\n", ctx->z_min, ctx->z_max);
  printf("  Theta extents: [%g, %g] rad\n", ctx->theta_min, ctx->theta_max);
  printf("  Mirror throat Z location (Z_m) = %g m\n", ctx->Z_m);
  printf("  Magnetic field parameter (mcB) = %g\n", ctx->mcB);
  printf("  Lorentzian width parameter (gamma) = %g\n", ctx->gamma);
  
  printf("\nGrid parameters:\n");
  printf("  Configuration space dimensions (cdim) = %d\n", ctx->cdim);
  printf("  Velocity space dimensions (vdim) = %d\n", ctx->vdim);
  printf("  Number of cells (Npsi, Ntheta, Nz, Nvpar, Nmu) = (%d, %d, %d, %d, %d)\n", ctx->Npsi, ctx->Ntheta, ctx->Nz, ctx->Nvpar, ctx->Nmu);
  printf("  Number of cells for electrons (Nvpar_elc, Nmu_elc) = (%d, %d)\n", ctx->Nvpar_elc, ctx->Nmu_elc);
  printf("  Polynomial order = %d\n", ctx->poly_order);
  printf("  Max ion parallel velocity (vpar_max_ion) = %g m/s (%.2f vti)\n", 
         ctx->vpar_max_ion, ctx->vpar_max_ion/ctx->vti);
  printf("  Max ion magnetic moment (mu_max_ion) = %g J/T\n", ctx->mu_max_ion);
  printf("  Max electron parallel velocity (vpar_max_elc) = %g m/s (%.2f vte)\n", 
         ctx->vpar_max_elc, ctx->vpar_max_elc/ctx->vte);
  printf("  Max electron magnetic moment (mu_max_elc) = %g J/T\n", ctx->mu_max_elc);

  
  printf("\nSource parameters:\n");
  printf("  Ion source amplitude = %g m^-3/s\n", ctx->ion_source_amplitude);
  printf("  Ion source sigma = %g\n", ctx->ion_source_sigma);
  printf("  Ion source temperature = %g eV\n", ctx->ion_source_temp/GKYL_ELEMENTARY_CHARGE);
  
  printf("\nSimulation parameters:\n");
  printf("  Total simulation time (t_end) = %g s\n", ctx->t_end);
  printf("  Total number of frames = %d\n", ctx->num_frames);
  printf("  Number of POA phases = %d\n", ctx->num_phases);
  printf("  Write phase-space frequency = %g\n", ctx->write_phase_freq);
  printf("  Integrated diagnostics calc frequency = %g\n", ctx->int_diag_calc_freq);
  printf("  Time-step failure tolerance = %g\n", ctx->dt_failure_tol);
  printf("  Maximum consecutive failures = %d\n", ctx->num_failures_max);
  
  int oap_count = 0, fdp_count = 0;
  double total_oap_time = 0.0, total_fdp_time = 0.0;
  for (int i = 0; i < ctx->num_phases; i++) {
    struct gk_poa_phase_params *p = &ctx->poa_phases[i];
    if (p->phase == GK_POA_OAP) {
      oap_count++;
      total_oap_time += p->duration;
    } else if (p->phase == GK_POA_FDP) {
      fdp_count++;
      total_fdp_time += p->duration;
    }
  }
  printf("\nPOA Summary:\n");
  printf("  Total OAP phases: %d (total time: %.3e s, %.1f%%)\n", 
         oap_count, total_oap_time, 100.0*total_oap_time/ctx->t_end);
  printf("  Total FDP phases: %d (total time: %.3e s, %.1f%%)\n", 
         fdp_count, total_fdp_time, 100.0*total_fdp_time/ctx->t_end);
  printf("  OAP/FDP time ratio: %.2f\n", total_oap_time/total_fdp_time);
  printf("\n");
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool force_calc, double dt)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
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
    gkyl_gyrokinetic_app_write_dt(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

struct time_frame_state {
  double t_curr; // Current simulation time.
  double t_end; // End time of current phase.
  int frame_curr; // Current frame.
  int num_frames; // Number of frames at the end of current phase.
};

void reset_io_triggers(struct gk_mirror_ctx *ctx, struct time_frame_state *tfs,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag)
{
  // Reset I/O triggers:
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  int frame_curr = tfs->frame_curr;
  int num_frames = tfs->num_frames;
  int num_int_diag_calc = ctx->int_diag_calc_freq*num_frames;

  // Prevent division by zero when frame_curr equals num_frames
  int frames_remaining = num_frames - frame_curr;
  double time_remaining = t_end - t_curr;

  trig_write_conf->dt = time_remaining / frames_remaining;
  trig_write_conf->tcurr = t_curr;
  trig_write_conf->curr = frame_curr;

  trig_write_phase->dt = time_remaining / (ctx->write_phase_freq * frames_remaining);
  trig_write_phase->tcurr = t_curr;
  trig_write_phase->curr = frame_curr;

  int diag_frames = GKYL_MAX2(frames_remaining, (num_int_diag_calc/num_frames) * frames_remaining);
  trig_calc_intdiag->dt = time_remaining / diag_frames;
  trig_calc_intdiag->tcurr = t_curr;
  trig_calc_intdiag->curr = frame_curr;
}

/////////// Boltzmann Electron POA phases

void run_phase(gkyl_gyrokinetic_app* app, struct gk_mirror_ctx *ctx, double num_steps,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag, struct time_frame_state *tfs,
  struct gk_poa_phase_params *pparams, int phase_idx)
{
  tfs->t_end = tfs->t_curr + pparams->duration;
  tfs->num_frames = tfs->frame_curr + pparams->num_frames;

  gkyl_gyrokinetic_app_cout(app, stdout, "----------------------------------------------\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Starting phase %d of type %s with parameters:\n",
    phase_idx + 1,
    (pparams->phase == GK_POA_OAP)?"OAP":"FDP");
  gkyl_gyrokinetic_app_cout(app, stdout, "  Duration = %g\n", pparams->duration);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Number of frames = %d\n", pparams->num_frames);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Alpha = %g\n", pparams->alpha);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Evolve field = %s\n", pparams->is_static_field?"No":"Yes");
  gkyl_gyrokinetic_app_cout(app, stdout, "  Positivity = %s\n", pparams->is_positivity_enabled?"Yes":"No");
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt multiplier = %s\n",
    fdot_multiplier_type_to_str(pparams->fdot_mult_type));
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt threshold = %g\n", pparams->f_threshold);
  gkyl_gyrokinetic_app_cout(app, stdout, "  CFL factor times omega max = %g\n", pparams->cfl_factor_times_omega_max);
  gkyl_gyrokinetic_app_cout(app, stdout, "----------------------------------------------\n");

  // Run an OAP or FDP.
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  
  // Reset I/O triggers:
  reset_io_triggers(ctx, tfs, trig_write_conf, trig_write_phase, trig_calc_intdiag);

  // Reset simulation parameters and function pointers.
  struct gkyl_gyrokinetic_collisionless collisionless_inp = {
    .type = GKYL_GK_COLLISIONLESS_ES,
    .scale_factor = pparams->alpha,
  };

  struct gkyl_gyrokinetic_fdot_multipliers fdot_mult_inp = {
    .num_multipliers = 1,
    .multiplier[0] = {
      .type = pparams->fdot_mult_type,
      .f_threshold = pparams->f_threshold,
      .cfl_factor_times_omega_max = pparams->cfl_factor_times_omega_max,
      .time_dilation_scale_const = pparams->time_dilation_scale_const,
      .cellwise_const = true,
      .write_diagnostics = true,
    },
  };
  struct gkyl_gyrokinetic_field field_inp = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx->me,
    .electron_charge = ctx->qe,
    .electron_temp = ctx->Te0,
    .polarization_bmag = ctx->B_p,
    .is_static = pparams->is_static_field,
  };
  struct gkyl_gyrokinetic_positivity positivity_inp = {
    .type = pparams->is_positivity_enabled? GKYL_GK_POSITIVITY_SHIFT : GKYL_GK_POSITIVITY_NONE,
    .write_diagnostics = pparams->is_positivity_enabled,
  };
  struct gkyl_gyrokinetic_damping damping_inp = {
    .type = pparams->damping_type,
    .rate_const = pparams->damping_rate_const,
    .write_rate = false,
    .write_fbar = true,
    .cellwise_const = false,
  };

  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "ion", collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "ion", fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "ion", positivity_inp);
  gkyl_gyrokinetic_app_reset_species_damping(app, t_curr, "ion", damping_inp);
  gkyl_gyrokinetic_app_reset_field(app, t_curr, field_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;

  // Initialize wall-time tracking for phase completion estimate
  struct timespec phase_start_time;
  clock_gettime(CLOCK_MONOTONIC, &phase_start_time);
  double t_start = tfs->t_curr;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps))
  {
    if (step%1000 == 1 || step==1)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    dt = t_end - t_curr; // Ensure we don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    if (step%1000 == 1 || step==1) {
      // Calculate elapsed wall time and estimated time remaining
      struct timespec current_time;
      clock_gettime(CLOCK_MONOTONIC, &current_time);
      double wall_time_elapsed = (current_time.tv_sec - phase_start_time.tv_sec) + 
                                  (current_time.tv_nsec - phase_start_time.tv_nsec) / 1e9;
      double sim_time_elapsed = t_curr - t_start;
      double sim_time_remaining = t_end - t_curr;
      
      double wall_time_per_sim_time = wall_time_elapsed / sim_time_elapsed;
      double wall_time_remaining = wall_time_per_sim_time * sim_time_remaining;
      int hours = (int)(wall_time_remaining / 3600);
      int minutes = (int)((wall_time_remaining - hours*3600) / 60);
      int seconds = (int)(wall_time_remaining - hours*3600 - minutes*60);
      
      double progress_pct = 100.0 * sim_time_elapsed / (sim_time_elapsed + sim_time_remaining);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g (phase %.1f%% complete, est. %dh %dm %ds remaining)\n", 
                                status.dt_actual, progress_pct, hours, minutes, seconds);
    }

    if (!status.success)
    {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, t_curr >= t_end, status.dt_actual);
    write_data(trig_write_conf, trig_write_phase, app, t_curr, t_curr >= t_end);

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
        calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, true, status.dt_actual);
        write_data(trig_write_conf, trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  tfs->t_curr = t_curr;
  tfs->frame_curr = tfs->frame_curr+pparams->num_frames;
}


//////////// Kinetic Electron POA ///////////////

void run_phase_kinetic_elc(gkyl_gyrokinetic_app* app, struct gk_mirror_ctx *ctx, double num_steps,
  struct gkyl_tm_trigger *trig_write_conf, struct gkyl_tm_trigger *trig_write_phase,
  struct gkyl_tm_trigger *trig_calc_intdiag, struct time_frame_state *tfs,
  struct gk_poa_phase_params *pparams, int phase_idx)
{
  tfs->t_end = tfs->t_curr + pparams->duration;
  tfs->num_frames = tfs->frame_curr + pparams->num_frames;

  gkyl_gyrokinetic_app_cout(app, stdout, "----------------------------------------------\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Starting phase %d of type %s with parameters:\n",
    phase_idx + 1,
    (pparams->phase == GK_POA_OAP)?"OAP":"FDP");
  gkyl_gyrokinetic_app_cout(app, stdout, "  Duration = %g\n", pparams->duration);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Number of frames = %d\n", pparams->num_frames);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Alpha ion = %g\n", pparams->alpha_ion);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Alpha elc = %g\n", pparams->alpha_elc);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Evolve field = %s\n", pparams->is_static_field?"No":"Yes");
  gkyl_gyrokinetic_app_cout(app, stdout, "  Positivity ion = %s\n", pparams->is_positivity_enabled_ion?"Yes":"No");
  gkyl_gyrokinetic_app_cout(app, stdout, "  Positivity elc = %s\n", pparams->is_positivity_enabled_elc?"Yes":"No");
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt multiplier ion = %s\n",
    fdot_multiplier_type_to_str(pparams->fdot_mult_type_ion));
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt multiplier elc = %s\n",
    fdot_multiplier_type_to_str(pparams->fdot_mult_type_elc));
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt threshold ion = %g\n", pparams->f_threshold_ion);
  gkyl_gyrokinetic_app_cout(app, stdout, "  df/dt threshold elc = %g\n", pparams->f_threshold_elc);
  gkyl_gyrokinetic_app_cout(app, stdout, "  CFL factor times omega max ion = %g\n", pparams->cfl_factor_times_omega_max_ion);
  gkyl_gyrokinetic_app_cout(app, stdout, "  CFL factor times omega max elc = %g\n", pparams->cfl_factor_times_omega_max_elc);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Time dilation scale const ion = %g\n", pparams->time_dilation_scale_const_ion);
  gkyl_gyrokinetic_app_cout(app, stdout, "  Time dilation scale const elc = %g\n", pparams->time_dilation_scale_const_elc);
  
  gkyl_gyrokinetic_app_cout(app, stdout, "----------------------------------------------\n");

  // Run an OAP or FDP.
  double t_curr = tfs->t_curr;
  double t_end = tfs->t_end;
  
  // Reset I/O triggers:
  reset_io_triggers(ctx, tfs, trig_write_conf, trig_write_phase, trig_calc_intdiag);

  // Reset simulation parameters and function pointers.
  struct gkyl_gyrokinetic_collisionless ion_collisionless_inp = {
    .type = GKYL_GK_COLLISIONLESS_ES,
    .scale_factor = pparams->alpha_ion,
  };
  struct gkyl_gyrokinetic_collisionless elc_collisionless_inp = {
    .type = GKYL_GK_COLLISIONLESS_ES,
    .scale_factor = pparams->alpha_elc,
  };
  struct gkyl_gyrokinetic_fdot_multipliers ion_fdot_mult_inp = {
    .num_multipliers = 1,
    .multiplier[0] = {
      .type = pparams->fdot_mult_type_ion,
      .f_threshold = pparams->f_threshold_ion,
      .cfl_factor_times_omega_max = pparams->cfl_factor_times_omega_max_ion,
      .time_dilation_scale_const = pparams->time_dilation_scale_const_ion,
      .dt_set_by_species = "elc",
      .cellwise_const = true,
      .write_diagnostics = true,
    },
  };
  struct gkyl_gyrokinetic_fdot_multipliers elc_fdot_mult_inp = {
    .num_multipliers = 1,
    .multiplier[0] = {
      .type = pparams->fdot_mult_type_elc,
      .f_threshold = pparams->f_threshold_elc,
      .cfl_factor_times_omega_max = pparams->cfl_factor_times_omega_max_elc,
      .time_dilation_scale_const = pparams->time_dilation_scale_const_elc,
      .dt_set_by_species = "ion",
      .cellwise_const = true,
      .write_diagnostics = true,
    },
  };
  struct gkyl_gyrokinetic_positivity ion_positivity_inp = {
    .type = pparams->is_positivity_enabled_ion? GKYL_GK_POSITIVITY_SHIFT : GKYL_GK_POSITIVITY_NONE,
    .write_diagnostics = pparams->is_positivity_enabled_ion,
  };
  struct gkyl_gyrokinetic_positivity elc_positivity_inp = {
    .type = pparams->is_positivity_enabled_elc? GKYL_GK_POSITIVITY_SHIFT : GKYL_GK_POSITIVITY_NONE,
    .write_diagnostics = pparams->is_positivity_enabled_elc,
  };
  struct gkyl_gyrokinetic_field field_inp = {
    .polarization_bmag = ctx->B_p,
    .kperpSq = pow(ctx->kperp, 2.),
    .is_static = pparams->is_static_field,
  };
  struct gkyl_gyrokinetic_damping damping_inp = {
    .type = pparams->damping_type,
    .rate_const = pparams->damping_rate_const,
    .write_rate = false,
    .write_fbar = true,
    .cellwise_const = false,
  };
  
  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "ion", ion_collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "ion", ion_fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "ion", ion_positivity_inp);
  gkyl_gyrokinetic_app_reset_species_damping(app, t_curr, "ion", damping_inp);
  
  gkyl_gyrokinetic_app_reset_species_collisionless(app, t_curr, "elc", elc_collisionless_inp);
  gkyl_gyrokinetic_app_reset_species_fdot_multiplier(app, t_curr, "elc", elc_fdot_mult_inp);
  gkyl_gyrokinetic_app_reset_species_positivity(app, t_curr, "elc", elc_positivity_inp);
  gkyl_gyrokinetic_app_reset_species_damping(app, t_curr, "elc", damping_inp);

  gkyl_gyrokinetic_app_reset_field(app, t_curr, field_inp);

  // Compute initial guess of maximum stable time-step.
  double dt = t_end - t_curr;

  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx->dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx->num_failures_max;

  // Initialize wall-time tracking for phase completion estimate
  struct timespec phase_start_time;
  clock_gettime(CLOCK_MONOTONIC, &phase_start_time);
  double t_start = tfs->t_curr;

  long step = 1;
  while ((t_curr < t_end) && (step <= num_steps))
  {
    if (step%1000 == 1 || step==1)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    dt = t_end - t_curr; // Ensure we don't step beyond t_end.
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
    if (step%1000 == 1 || step==1) {
      // Calculate elapsed wall time and estimated time remaining
      struct timespec current_time;
      clock_gettime(CLOCK_MONOTONIC, &current_time);
      double wall_time_elapsed = (current_time.tv_sec - phase_start_time.tv_sec) + 
                                  (current_time.tv_nsec - phase_start_time.tv_nsec) / 1e9;
      double sim_time_elapsed = t_curr - t_start;
      double sim_time_remaining = t_end - t_curr;
      
      double wall_time_per_sim_time = wall_time_elapsed / sim_time_elapsed;
      double wall_time_remaining = wall_time_per_sim_time * sim_time_remaining;
      int hours = (int)(wall_time_remaining / 3600);
      int minutes = (int)((wall_time_remaining - hours*3600) / 60);
      int seconds = (int)(wall_time_remaining - hours*3600 - minutes*60);
      
      double progress_pct = 100.0 * sim_time_elapsed / (sim_time_elapsed + sim_time_remaining);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g (phase %.1f%% complete, est. %dh %dm %ds remaining)\n", 
                                status.dt_actual, progress_pct, hours, minutes, seconds);
    }

    if (!status.success)
    {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, t_curr >= t_end, status.dt_actual);
    write_data(trig_write_conf, trig_write_phase, app, t_curr, t_curr >= t_end);

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
        calc_integrated_diagnostics(trig_calc_intdiag, app, t_curr, true, status.dt_actual);
        write_data(trig_write_conf, trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

    step += 1;
  }

  tfs->t_curr = t_curr;
  tfs->frame_curr = tfs->frame_curr+pparams->num_frames;
}


void
run_poa_simulation(struct gkyl_gk app_inp, struct gk_mirror_ctx ctx, struct gkyl_app_args app_args, bool is_kinetic_elc){
  // Create app object.

  struct timespec app_init_start;
  clock_gettime(CLOCK_MONOTONIC, &app_init_start);

  printf("App initialization started ...\n");
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  struct timespec app_init_end;
  clock_gettime(CLOCK_MONOTONIC, &app_init_end);
  double app_init_time = (app_init_end.tv_sec - app_init_start.tv_sec) + 
                         (app_init_end.tv_nsec - app_init_start.tv_nsec) / 1e9;
  gkyl_gyrokinetic_app_cout(app, stdout, "App initialization time: %.3f seconds\n\n", app_init_time);

  // Triggers for IO.
  struct gkyl_tm_trigger trig_write_conf, trig_write_phase, trig_calc_intdiag;

  struct time_frame_state tfs = {
    .t_curr = 0.0, // Initial simulation time.
    .frame_curr = 0, // Initial frame.
    .t_end = ctx.poa_phases[0].duration, // Final time of 1st phase.
    .num_frames = ctx.poa_phases[0].num_frames, // Number of frames in 1st phase.
  };

  int phase_idx_init = 0, phase_idx_end = ctx.num_phases; // Initial and final phase index.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n", gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    tfs.frame_curr = status.frame;
    tfs.t_curr = status.stime;

    // Find out what phase we are in.
    double time_count = 0.0;
    int frame_count = 0;
    int pit_curr = 0;
    for (int pit=0; pit<ctx.num_phases; pit++) {
      time_count += ctx.poa_phases[pit].duration;
      frame_count += ctx.poa_phases[pit].num_frames;
      if ((tfs.t_curr <= time_count) && (tfs.frame_curr <= frame_count)) {
        pit_curr = pit;
        break;
      }
    };
    phase_idx_init = pit_curr;

    // Calculate time and frames at the START of the current phase
    double phase_start_time = time_count - ctx.poa_phases[phase_idx_init].duration;
    int phase_start_frame = frame_count - ctx.poa_phases[phase_idx_init].num_frames;
    
    // Calculate how much of the phase has elapsed
    double time_elapsed_in_phase = tfs.t_curr - phase_start_time;
    int frames_elapsed_in_phase = tfs.frame_curr - phase_start_frame;

    // Change the duration and number frames so this phase reaches the expected
    // time and number of frames and not beyond.
    struct gk_poa_phase_params *pparams = &ctx.poa_phases[phase_idx_init];
    int original_frames = pparams->num_frames;
    double original_duration = pparams->duration;
    pparams->num_frames = frame_count - tfs.frame_curr;
    pparams->duration = time_count - tfs.t_curr;

    gkyl_gyrokinetic_app_cout(app, stdout, "\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "==============================================\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "          RESTART INFORMATION\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "==============================================\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d at time = %.6e s\n", tfs.frame_curr, tfs.t_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, "Restart phase index: %d (of %d total phases)\n", phase_idx_init+1, ctx.num_phases);
    gkyl_gyrokinetic_app_cout(app, stdout, "Restart phase type: %s\n", 
      (pparams->phase == GK_POA_OAP) ? "OAP (Orbit Averaging Phase)" : "FDP (Full Dynamics Phase)");
    gkyl_gyrokinetic_app_cout(app, stdout, "\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "Phase completion status:\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "  Time elapsed in current phase: %.6e s (%.1f%% of phase)\n", 
      time_elapsed_in_phase, 
      100.0*time_elapsed_in_phase/original_duration);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Frames completed in current phase: %d (of %d total)\n", 
      frames_elapsed_in_phase, original_frames);
    gkyl_gyrokinetic_app_cout(app, stdout, "\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "Remaining simulation:\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "  Remaining time in current phase: %.6e s\n", pparams->duration);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Remaining frames in current phase: %d\n", pparams->num_frames);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Total remaining phases: %d\n", ctx.num_phases - phase_idx_init);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Total remaining time: %.6e s (%.1f%% of simulation)\n", 
      ctx.t_end - tfs.t_curr, 100.0*(ctx.t_end - tfs.t_curr)/ctx.t_end);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Total remaining frames: %d (of %d total)\n", 
      ctx.num_frames - tfs.frame_curr, ctx.num_frames);
    gkyl_gyrokinetic_app_cout(app, stdout, "==============================================\n");
    gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  }
  else {

    clock_gettime(CLOCK_MONOTONIC, &app_init_start);
    gkyl_gyrokinetic_app_cout(app, stdout, "Applying initial conditions ...\n");
    gkyl_gyrokinetic_app_apply_ic(app, tfs.t_curr);
    clock_gettime(CLOCK_MONOTONIC, &app_init_end);
    double ic_apply_time = (app_init_end.tv_sec - app_init_start.tv_sec) + (app_init_end.tv_nsec - app_init_start.tv_nsec) / 1e9;
    gkyl_gyrokinetic_app_cout(app, stdout, "Initial condition application time: %.3f seconds\n\n", ic_apply_time);

    // Write out ICs.
    reset_io_triggers(&ctx, &tfs, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag);

    calc_integrated_diagnostics(&trig_calc_intdiag, app, tfs.t_curr, true, -1.0);
    write_data(&trig_write_conf, &trig_write_phase, app, tfs.t_curr, true);
  }

  if (app_args.num_steps != INT_MAX)
    phase_idx_end = 1;

  // Loop over number of number of phases;
  for (int pit=phase_idx_init; pit<phase_idx_end; pit++) {
    struct gk_poa_phase_params *phase_params = &ctx.poa_phases[pit];
    if (is_kinetic_elc) {
      run_phase_kinetic_elc(app, &ctx, app_args.num_steps, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag, &tfs, phase_params, pit);
    } else {
      run_phase(app, &ctx, app_args.num_steps, &trig_write_conf, &trig_write_phase, &trig_calc_intdiag, &tfs, phase_params, pit);
    }
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0)
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

  freeresources:
  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
}



//////////// MAPC2P functions ////////////////


double
psi_RZ(double RIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double psi = 0.5 * pow(RIn, 2.) * mcB *
               (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
                1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))));
  return psi;
}

double
R_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rout = sqrt(2. * psiIn / (mcB * 
    (1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.))) +
     1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))))));
  return Rout;
}

void
Bfield_psiZ(double psiIn, double ZIn, void *ctx, double *BRad, double *BZ, double *Bmag)
{
  struct gk_mirror_ctx *app = ctx;
  double mcB = app->mcB;
  double gamma = app->gamma;
  double Z_m = app->Z_m;

  double Rcoord = R_psiZ(psiIn, ZIn, ctx);

  BRad[0] = -(1. / 2.) * Rcoord * mcB *
          (-2. * (ZIn - Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn - Z_m) / gamma, 2.), 2.)))
           -2. * (ZIn + Z_m) / (M_PI * pow(gamma, 3.) * (pow(1.0 + pow((ZIn + Z_m) / gamma, 2.), 2.))));

  BZ[0] = mcB *
        ( 1. / (M_PI * gamma * (1. + pow((ZIn - Z_m) / gamma, 2.)))
         +1. / (M_PI * gamma * (1. + pow((ZIn + Z_m) / gamma, 2.))) );

  Bmag[0] = sqrt(pow(BRad[0], 2) + pow(BZ[0], 2));
}

double
integrand_z_psiZ(double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_in;
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, ZIn, ctx, &BRad, &BZ, &Bmag);
  return Bmag / BZ;
}

double
z_psiZ(double psiIn, double ZIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double eps = 0.0;
  app->psi_in = psiIn;
  struct gkyl_qr_res integral;
  if (eps <= ZIn)
  {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, eps, ZIn, 7, 1e-14);
  }
  else
  {
    integral = gkyl_dbl_exp(integrand_z_psiZ, ctx, ZIn, eps, 7, 1e-14); 
    integral.res = -integral.res;
  }
  return integral.res;
}

// Invert z(Z) via root-finding.
double
root_Z_psiz(double Z, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  return app->z_in - z_psiZ(app->psi_in, Z, ctx);
}

double
Z_psiz(double psiIn, double zIn, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double maxL = app->Z_max - app->Z_min;
  double eps = maxL / app->Nz;   // Interestingly using a smaller eps yields larger errors in some geo quantities.
  app->psi_in = psiIn;
  app->z_in = zIn;
  struct gkyl_qr_res Zout;
  if (0.0 <= zIn)
  {
    double fl = root_Z_psiz(-eps, ctx);
    double fr = root_Z_psiz(app->Z_max + eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, -eps, app->Z_max + eps, fl, fr, 1000, 1e-14);
  }
  else
  {
    double fl = root_Z_psiz(app->Z_min - eps, ctx);
    double fr = root_Z_psiz(eps, ctx);
    Zout = gkyl_ridders(root_Z_psiz, ctx, app->Z_min - eps, eps, fl, fr, 1000, 1e-14);
  }
  return Zout.res;
}

// Geometry evaluation functions for the gk app
void
mapc2p(double t, const double *xc, double *GKYL_RESTRICT xp, void *ctx)
{
  double psi = xc[0], theta = xc[1], z = xc[2];

  double Z = Z_psiz(psi, z, ctx);
  double R = R_psiZ(psi, Z, ctx);

  // Cartesian coordinates on plane perpendicular to Z axis.
  double x = R * cos(theta);
  double y = R * sin(theta);

  xp[0] = x;  xp[1] = y;  xp[2] = Z;
}

void
bfield_func(double t, const double *xc, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xc[2];
  double psi = psi_RZ(app->RatZeq0, 0.0, ctx); // Magnetic flux function psi of field line.
  double Z = Z_psiz(psi, z, ctx);
  double BRad, BZ, Bmag;
  Bfield_psiZ(psi, Z, ctx, &BRad, &BZ, &Bmag);

  double phi = xc[1];
  // zc are computational coords. 
  // Set Cartesian components of magnetic field.
  fout[0] = BRad*cos(phi);
  fout[1] = BRad*sin(phi);
  fout[2] = BZ;
}

void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}