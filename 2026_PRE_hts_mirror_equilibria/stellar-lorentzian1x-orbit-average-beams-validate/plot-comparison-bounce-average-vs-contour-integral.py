"""
Beam source comparison: spatial vs orbit-averaged vs contour-integrated sources.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize
from multiprocessing import Pool, cpu_count
import matplotlib
import postgkyl as pg

matplotlib.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 12,
    'axes.titlesize': 18,
    'axes.labelsize': 22,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
})

def _bounce_avg_worker(args):
  """Worker function for parallel contour integration."""
  z, vpar, mu, obj = args
  return obj.bounce_average_beam_source(z, vpar, mu)

class BeamSourceComparison:
  """
  Class to compare different beam source models in a mirror machine.
  
  Models:
  1. Spatial beam source (Dorf et al. 2025): Gaussian in z, beam distribution in v
  2. Orbit-averaged source: Maps velocities to midplane using energy conservation
  3. Contour-integrated source: Averages spatial source over energy contours
  """
  
  # Physical constants
  eV = 1.602176634e-19
  mp = 1.67262192369e-27
  mi = 2.014 * mp  # Deuteron mass
  
  def __init__(self, 
         # Geometry parameters
         mcB=6.51292,
         gamma_geo=0.124904,
         Z_m=0.98,
         B_p=0.53,
         # Beam parameters
         gamma_spatial=4198.25,
         gamma_midmap=1110.13,
         T_beam_eV_spatial=200,
         T_beam_eV_midmap=200,
         E_beam_eV_spatial=25000,
         E_beam_eV_midmap=25000,
         Lb=0.2,
         use_potential=True):
    """
    Initialize beam source comparison.
    
    Parameters:
      mcB : float - Magnetic field parameter
      gamma_geo : float - Lorentzian width parameter
      Z_m : float - Mirror throat location (m)
      B_p : float - Reference magnetic field (T)
      gamma_spatial : float - Spatial source normalization
      gamma_midmap : float - Midplane source normalization
      T_beam_eV_spatial : float - Beam temperature for spatial source (eV)
      T_beam_eV_midmap : float - Beam temperature for midplane source (eV)
      E_beam_eV_spatial : float - Beam energy for spatial source (eV)
      E_beam_eV_midmap : float - Beam energy for midplane source (eV)
      Lb : float - Spatial source width (m)
    """
    # Geometry
    self.mcB = mcB
    self.gamma_geo = gamma_geo # Gamma used in Bmag calculation
    self.Z_m = Z_m
    self.B_p = B_p
    self.qi = self.eV
    
    # Beam parameters
    self.gamma_spatial = gamma_spatial  # Same normalization for spatial source
    self.gamma_midmap = gamma_midmap
    self.T_beam_spatial = T_beam_eV_spatial * self.eV
    self.T_beam_midmap = T_beam_eV_midmap * self.eV
    self.E_beam_spatial = E_beam_eV_spatial * self.eV
    self.E_beam_midmap = E_beam_eV_midmap * self.eV
    self.v_beam_spatial = np.sqrt(self.E_beam_spatial / self.mi)
    self.v_beam_midmap = np.sqrt(self.E_beam_midmap / self.mi)
    self.sigma_beam_spatial = 2 * self.T_beam_spatial / self.mi
    self.sigma_beam_midmap = 2 * self.T_beam_midmap / self.mi
    self.Lb = Lb
    
    # Precompute bmag spline for fast evaluation in contour integrals
    self._z_grid = np.linspace(-1.0, 1.0, 2001)
    self._bmag_grid = self._bmag_raw(self._z_grid)
    self._bmag_spline = CubicSpline(self._z_grid, self._bmag_grid)
    
    self._use_potential = use_potential
    if use_potential:
      data_phi = pg.GData('python-data/gk_lorentzian_mirror-field_332.gkyl')
      interp_phi = pg.GInterpModal(data_phi, 1, 'ms')
      x, phi_data = interp_phi.interpolate()
      phi_data = np.squeeze(phi_data)

      data_mc2p = pg.GData('python-data/gk_lorentzian_mirror-mapc2p.gkyl')
      interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
      nodes_Z = interp_mc2p.interpolate(2)[1]
      nodes_Z = np.squeeze(nodes_Z)
      # Store raw data for pickling compatibility
      self._phi_z_grid = nodes_Z
      self._phi_grid = phi_data
      self._phi_spline = CubicSpline(nodes_Z, phi_data)
    else:
      self._phi_z_grid = None
      self._phi_grid = None
      self._phi_spline = None
  
  def __getstate__(self):
    """Prepare state for pickling (needed for multiprocessing)."""
    state = self.__dict__.copy()
    # Remove spline objects - they will be recreated
    state['_bmag_spline'] = None
    state['_phi_spline'] = None
    return state
  
  def __setstate__(self, state):
    """Restore state after unpickling."""
    self.__dict__.update(state)
    # Recreate splines from stored data
    self._bmag_spline = CubicSpline(self._z_grid, self._bmag_grid)
    if self._use_potential and self._phi_z_grid is not None:
      self._phi_spline = CubicSpline(self._phi_z_grid, self._phi_grid)
  
  def phi(self, z):
    """Evaluate electrostatic potential at position z."""
    if not self._use_potential:
      return np.zeros_like(np.atleast_1d(z)) if hasattr(z, '__len__') else 0.0
    return self._phi_spline(z)
    
  def _bmag_raw(self, Z):
    """Raw magnetic field calculation (used for spline construction)."""
    Z = np.atleast_1d(Z)
    L1 = 1.0 / (1 + ((Z - self.Z_m) / self.gamma_geo)**2)
    L2 = 1.0 / (1 + ((Z + self.Z_m) / self.gamma_geo)**2)
    Bmag = self.mcB / (np.pi * self.gamma_geo) * (L1 + L2)
    return Bmag.item() if Bmag.size == 1 else Bmag
  
  def bmag(self, Z, use_spline=False):
    """Compute magnetic field magnitude at position Z.
    
    Args:
      Z: Position(s) to evaluate
      use_spline: If True, use precomputed spline (faster for many evaluations)
    """
    if use_spline:
      return self._bmag_spline(Z)
    return self._bmag_raw(Z)
  
  def dbmag_dz(self, Z):
    """Compute derivative of magnetic field dB/dZ."""
    Z = np.atleast_1d(Z)
    # d/dZ of 1/(1 + ((Z - Z_m)/gamma)^2) = -2(Z - Z_m) / (gamma^2 * (1 + ((Z-Z_m)/gamma)^2)^2)
    dL1 = -2 * (Z - self.Z_m) / (self.gamma_geo**2 * (1 + ((Z - self.Z_m) / self.gamma_geo)**2)**2)
    dL2 = -2 * (Z + self.Z_m) / (self.gamma_geo**2 * (1 + ((Z + self.Z_m) / self.gamma_geo)**2)**2)
    dBdZ = self.mcB / (np.pi * self.gamma_geo) * (dL1 + dL2)
    return dBdZ.item() if dBdZ.size == 1 else dBdZ
  
  def midplane_beam_source(self, vpar, mu, gamma_beam, v_beam, sigma_beam):
    """Beam source evaluated at midplane velocities."""
    vperp = np.sqrt(2.0 * mu * self.B_p / self.mi)
    return gamma_beam * np.exp(-((np.abs(vpar) - v_beam)**2 + 
                      (vperp - v_beam)**2) / sigma_beam)
  
  def spatial_beam_source(self, z, vpar, mu):
    """Spatial beam source with Gaussian z-dependence (Dorf et al. 2025)."""
    zdep = np.exp(-(z / self.Lb)**2)
    vdep = self.midplane_beam_source(vpar, mu, self.gamma_spatial, self.v_beam_spatial, self.sigma_beam_spatial)
    return zdep * vdep
  
  def midplane_mapped_beam_source(self, z, vpar, mu):
    """Orbit-averaged source: map vpar to midplane using energy conservation.
    
    Uses E = 0.5*m*vpar² + mu*B + q*phi = const along orbit.
    Maps vpar(z) -> vpar(0) at midplane.
    
    Returns 0 for particles that cannot reach the midplane (trapped by potential).
    """
    z = np.atleast_1d(z)
    vpar = np.atleast_1d(vpar)
    mu = np.atleast_1d(mu)
    
    # vpar_midplane² = vpar² + 2*mu*(B(z)-B(0))/m + 2*q*(phi(z)-phi(0))/m
    vpar_sq_midplane = (vpar**2 
                        + 2 * mu * (self.bmag(z) - self.bmag(0)) / self.mi 
                        + 2 * self.qi / self.mi * (self.phi(z) - self.phi(0)))
    
    # If vpar_sq_midplane < 0, particle cannot reach midplane - set source to 0
    valid = vpar_sq_midplane >= 0
    vpar_midplane = np.where(valid, np.sqrt(np.maximum(vpar_sq_midplane, 0)), 0.0)
    
    result = self.midplane_beam_source(vpar_midplane, mu, 
                                       self.gamma_midmap, self.v_beam_midmap, self.sigma_beam_midmap)
    # Zero out invalid regions and outside mirror throat
    result = np.where(valid & (np.abs(z) <= self.Z_m), result, 0.0)
    return result.item() if result.size == 1 else result
  
  def plot_energy_diagnostics(self, z_cut=0.0):
    """Plot diagnostic information about the energy terms."""
    z_coords = np.linspace(-0.95, 0.95, 200)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Panel 1: B(z) and phi(z)
    ax1 = axes[0, 0]
    ax1.plot(z_coords, self.bmag(z_coords), 'b-', label=r'$B(z)$ [T]')
    ax1.set_xlabel('z (m)')
    ax1.set_ylabel('B (T)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1_twin = ax1.twinx()
    phi_vals = self.phi(z_coords)
    ax1_twin.plot(z_coords, phi_vals, 'r-', label=r'$\phi(z)$ [V]')
    ax1_twin.set_ylabel(r'$\phi$ (V)', color='r')
    ax1_twin.tick_params(axis='y', labelcolor='r')
    ax1.set_title('Magnetic field and electrostatic potential')
    ax1.axhline(self.bmag(0), color='b', linestyle='--', alpha=0.5, label=f'B(0)={self.bmag(0):.2f} T')
    ax1_twin.axhline(self.phi(0), color='r', linestyle='--', alpha=0.5)
    ax1.legend(loc='upper left')
    
    # Panel 2: Energy terms vs z for a typical particle
    ax2 = axes[0, 1]
    mu_typical = 0.5 * self.mi * self.v_beam_spatial**2 / self.B_p
    vpar_typical = self.v_beam_spatial
    
    mu_B = mu_typical * self.bmag(z_coords)
    q_phi = self.qi * self.phi(z_coords)
    kinetic_par = 0.5 * self.mi * vpar_typical**2 * np.ones_like(z_coords)
    total_E = kinetic_par + mu_B + q_phi
    
    ax2.plot(z_coords, mu_B / (1000 * self.eV), 'g-', label=r'$\mu B$')
    ax2.plot(z_coords, q_phi / (1000 * self.eV), 'r-', label=r'$q\phi$')
    ax2.plot(z_coords, kinetic_par / (1000 * self.eV), 'b--', label=r'$\frac{1}{2}m v_\parallel^2$ (at z)')
    ax2.plot(z_coords, total_E / (1000 * self.eV), 'k-', linewidth=2, label='Total E')
    ax2.set_xlabel('z (m)')
    ax2.set_ylabel('Energy (keV)')
    ax2.set_title(f'Energy components\n(mu={mu_typical:.2e}, vpar={vpar_typical/1e6:.2f} Mm/s)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: vpar_midplane² for different starting z
    ax3 = axes[1, 0]
    for vpar_frac in [0.5, 1.0, 1.5]:
      vpar_test = vpar_frac * self.v_beam_spatial
      vpar_sq_mid = (vpar_test**2 
                     + 2 * mu_typical * (self.bmag(z_coords) - self.bmag(0)) / self.mi 
                     + 2 * self.qi / self.mi * (self.phi(z_coords) - self.phi(0)))
      ax3.plot(z_coords, vpar_sq_mid / 1e12, label=f'vpar = {vpar_frac:.1f} v_beam')
    ax3.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax3.set_xlabel('z (m)')
    ax3.set_ylabel(r'$v_{\parallel,midplane}^2$ ($10^{12}$ m²/s²)')
    ax3.set_title(r'Midplane mapping: $v_{\parallel,0}^2 = v_\parallel^2 + \frac{2\mu}{m}(B-B_0) + \frac{2q}{m}(\phi-\phi_0)$')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(bottom=-0.5)
    
    # Panel 4: phi(z) - phi(0) contribution
    ax4 = axes[1, 1]
    delta_phi = self.phi(z_coords) - self.phi(0)
    delta_B = self.bmag(z_coords) - self.bmag(0)
    
    # Convert to equivalent vpar² change
    delta_vpar_sq_from_B = 2 * mu_typical * delta_B / self.mi
    delta_vpar_sq_from_phi = 2 * self.qi * delta_phi / self.mi
    
    ax4.plot(z_coords, delta_vpar_sq_from_B / 1e12, 'g-', label=r'From $\Delta B$: $\frac{2\mu}{m}(B-B_0)$')
    ax4.plot(z_coords, delta_vpar_sq_from_phi / 1e12, 'r-', label=r'From $\Delta\phi$: $\frac{2q}{m}(\phi-\phi_0)$')
    ax4.plot(z_coords, (delta_vpar_sq_from_B + delta_vpar_sq_from_phi) / 1e12, 'k-', 
             linewidth=2, label='Total')
    ax4.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax4.set_xlabel('z (m)')
    ax4.set_ylabel(r'$\Delta v_\parallel^2$ ($10^{12}$ m²/s²)')
    ax4.set_title(r'Contributions to midplane velocity mapping')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

  def find_bounce_points(self, E, mu, z_min=-0.98, z_max=0.98):
    """Find bounce points where E = mu * B(z) + q * phi(z).
    
    At bounce points, vpar = 0, so all kinetic energy is in perpendicular motion.
    """
    if mu <= 0:
      return z_min, z_max
    
    # Function whose roots give bounce points: E - mu*B(z) - q*phi(z) = 0
    def f(zp):
      return E - mu * self.bmag(zp) - self.qi * self.phi(zp)
    
    # Check if particle can access any region (check at midplane)
    f_mid = f(0)
    if f_mid < 0:
      return None  # Particle cannot exist at midplane - energy too low
    
    # Find left bounce point (between z_min and 0)
    f_left = f(z_min)
    if f_left * f_mid < 0:
      z_left = brentq(f, z_min, 0, xtol=1e-12)
    elif f_left >= 0:
      z_left = z_min  # Accessible all the way to z_min
    else:
      return None  # No accessible region on left
    
    # Find right bounce point (between 0 and z_max)
    f_right = f(z_max)
    if f_mid * f_right < 0:
      z_right = brentq(f, 0, z_max, xtol=1e-12)
    elif f_right >= 0:
      z_right = z_max  # Accessible all the way to z_max
    else:
      return None  # No accessible region on right
    
    return z_left, z_right
  
  def bounce_average_beam_source(self, z, vpar, mu, z_min=-0.98, z_max=0.98, use_fast=True):
    """
    Compute orbit-averaged source by integrating along energy contours.
    Uses time weighting: <S> = (∫ S dt) / (∫ dt) = (∫ S dz/|v_par|) / (∫ dz/|v_par|)
    
    Args:
      use_fast: If True, use fixed-order Gaussian quadrature with spline bmag (much faster)
    """
    E = 0.5 * self.mi * vpar**2 + mu * self.bmag(z) + self.qi * self.phi(z)
    
    # Find bounce points
    bounce_pts = self.find_bounce_points(E, mu, z_min, z_max)
    if bounce_pts is None:
      return 0.0
    z_left, z_right = bounce_pts
    
    if use_fast:
      result = self._bounce_average_integral_fast(E, mu, z_left, z_right)
    else:
      result = self._bounce_average_integral_adaptive(E, mu, z_left, z_right)

    result = np.where(np.abs(z) > self.Z_m, 0.0, result)
    return result
  
  def _bounce_average_integral_fast(self, E, mu, z_left, z_right, n_quad=32):
    """
    Fast contour integral using fixed-order Gaussian quadrature.
    
    Uses change of variables to handle sqrt singularity at bounce points:
    z = z_mid + (z_half) * sin(theta), theta in [-pi/2, pi/2]
    dz = z_half * cos(theta) d(theta)
    
    This removes the 1/sqrt singularity since vpar ~ sqrt(cos(theta)) near boundaries.
    """
    z_mid = 0.5 * (z_left + z_right)
    z_half = 0.5 * (z_right - z_left)
    
    two_over_mi = 2.0 / self.mi
    
    def integrand_transformed(theta):
      """Integrand in transformed coordinates (vectorized)."""
      theta = np.atleast_1d(theta)
      cos_th = np.cos(theta)
      z = z_mid + z_half * np.sin(theta)
      
      # Use spline for fast bmag evaluation
      B = self._bmag_spline(z)
      phi = self.phi(z)
      kinetic = E - mu * B - self.qi * phi
      vpar_sq = two_over_mi * np.maximum(kinetic, 1e-30)
      vpar = np.sqrt(vpar_sq)
      jacobian = z_half * cos_th  # dz/dtheta
      source = self.spatial_beam_source(z, vpar, mu)
      weight = jacobian / vpar
      return source * weight, weight
    
    # Use fixed_quad which is much faster than adaptive quad
    # Integrate from -pi/2 to pi/2
    def num_integrand(theta):
      num, _ = integrand_transformed(theta)
      return num
    
    def den_integrand(theta):
      _, den = integrand_transformed(theta)
      return den
    
    numerator, _ = fixed_quad(num_integrand, -np.pi/2, np.pi/2, n=n_quad)
    denominator, _ = fixed_quad(den_integrand, -np.pi/2, np.pi/2, n=n_quad)
    result = numerator / denominator if denominator > 1e-30 else 0.0
    return result
  
  def _bounce_average_integral_adaptive(self, E, mu, z_left, z_right):
    """Original adaptive quadrature method (slower but more accurate)."""
    def vpar_at_z(zp):
      return np.sqrt(max((2.0 / self.mi) * (E - mu * self.bmag(zp) - self.qi * self.phi(zp)), 1e-30))
    
    def integrand_num(zp):
      vp = vpar_at_z(zp)
      return self.spatial_beam_source(zp, vp, mu) / vp
    
    def integrand_den(zp):
      return 1.0 / vpar_at_z(zp)
    
    # Suppress integration warnings - near bounce points 1/vpar diverges but integral is still finite
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      numerator, _ = quad(integrand_num, z_left, z_right, limit=200, epsabs=1e-8, epsrel=1e-6)
      denominator, _ = quad(integrand_den, z_left, z_right, limit=200, epsabs=1e-8, epsrel=1e-6)
    
    return numerator / denominator if denominator > 1e-30 else 0.0
  
  def _compute_bounce_average_parallel(self, Z, VPAR, MU, n_workers=None):
    """Compute contour integral in parallel using multiprocessing."""
    shape = Z.shape
    args_list = [(z, vpar, mu, self) for z, vpar, mu in 
           zip(Z.flatten(), VPAR.flatten(), MU.flatten())]
    
    if n_workers is None:
      n_workers = cpu_count()
    print(f"Using multiprocessing with {n_workers} workers for {len(args_list)} integrals...")
    with Pool(n_workers) as pool:
      results = pool.map(_bounce_avg_worker, args_list)
    print("Done.")
    
    return np.array(results).reshape(shape)
  
  def compute_all_sources_on_grid(self, z_coords, vpar_coords, mu_coords, n_workers=None):
    """
    Compute all three source models on a grid.
    
    Returns:
      dict with keys 'spatial', 'midplane_mapped', 'bounce_average' containing 3D arrays
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    
    spatial = self.spatial_beam_source(Z, VPAR, MU)
    midmap = self.midplane_mapped_beam_source(Z, VPAR, MU)
    
    # Parallel computation for contour integral
    bounce_avg = self._compute_bounce_average_parallel(Z, VPAR, MU, n_workers)
    
    return {
      'Z': Z, 'VPAR': VPAR, 'MU': MU,
      'spatial': spatial,
      'midplane_mapped': midmap,
      'bounce_average': bounce_avg
    }

  def _integrate_M0_from_3d(self, source_3d, z_coords, vpar_coords, mu_coords):
    """Helper: compute M0 flux from precomputed 3D source array.
    
    M0 = ∫∫∫ S(z, vpar, mu) * jacobian dz dvpar dmu
    """
    jacobian = 2 * np.pi / self.mi
    vals = source_3d * jacobian
    # Integrate over vpar, then z, then mu
    vals = np.trapz(vals, vpar_coords, axis=1)
    vals = np.trapz(vals, z_coords, axis=0)
    return np.trapz(vals, mu_coords, axis=0)
  
  def _integrate_power_from_3d(self, source_3d, Z_grid, VPAR_grid, MU_grid, 
                               z_coords, vpar_coords, mu_coords):
    """Helper: compute power from precomputed 3D source array.
    
    Power = ∫∫∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    """
    energy = 0.5 * self.mi * VPAR_grid**2 + MU_grid * self.bmag(Z_grid) + self.qi * self.phi(Z_grid)
    jacobian = 2 * np.pi / self.mi
    vals = source_3d * energy * jacobian
    # Integrate over vpar, then z, then mu
    vals = np.trapz(vals, vpar_coords, axis=1)
    vals = np.trapz(vals, z_coords, axis=0)
    return np.trapz(vals, mu_coords, axis=0)

  def compute_M0_midplane_mapped(self, z_coords, vpar_coords, mu_coords):
    """Compute integrated M0 flux for orbit-averaged source only (3D integral).
    
    M0 = ∫∫∫ S(z, vpar, mu) * jacobian dz dvpar dmu
    
    This is a lightweight version that only computes the orbit-averaged source,
    useful for optimization loops where we don't need all three sources.
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    midmap = self.midplane_mapped_beam_source(Z, VPAR, MU)
    return self._integrate_M0_from_3d(midmap, z_coords, vpar_coords, mu_coords)
  
  def compute_power_midplane_mapped(self, z_coords, vpar_coords, mu_coords):
    """Compute total power for orbit-averaged source only (3D integral).
    
    Power = ∫∫∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    
    This is a lightweight version that only computes the orbit-averaged source,
    useful for optimization loops where we don't need all three sources.
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    midmap = self.midplane_mapped_beam_source(Z, VPAR, MU)
    return self._integrate_power_from_3d(midmap, Z, VPAR, MU, z_coords, vpar_coords, mu_coords)

  def compute_M0_fluxes(self, z_coords, vpar_coords, mu_coords, results=None, n_workers=None):
    """Compute integrated M0 fluxes for each source model (3D integral over z, vpar, mu).
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    fluxes = {}
    for key in ['spatial', 'midplane_mapped', 'bounce_average']:
      fluxes[key] = self._integrate_M0_from_3d(results[key], z_coords, vpar_coords, mu_coords)
    
    print("\nTotal M0 fluxes (particles per second per m^3):")
    print(f"  Spatial source:    {fluxes['spatial']:.3e} s^-1 m^-3")
    print(f"  Contour integral:  {fluxes['bounce_average']:.3e} s^-1 m^-3")
    print(f"  Orbit-averaged:    {fluxes['midplane_mapped']:.3e} s^-1 m^-3")

    # pgkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 0 info
    # pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel -c0 info

    return fluxes
  
  def compute_power(self, z_coords, vpar_coords, mu_coords, results=None, n_workers=None):
    """Compute total power (energy per second) deposited by each source model.
    
    Power = ∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    
    Returns:
      dict with keys 'spatial', 'midplane_mapped', 'bounce_average' containing power in Watts/m^3
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    Z = results['Z']
    VPAR = results['VPAR']
    MU = results['MU']
    
    power = {}
    for key in ['spatial', 'midplane_mapped', 'bounce_average']:
      power[key] = self._integrate_power_from_3d(results[key], Z, VPAR, MU, 
                                                  z_coords, vpar_coords, mu_coords)
    
    print("\nTotal power deposited (energy per second per m^3):")
    print(f"  Spatial source:    {power['spatial']:.3e} W/m^3 = {power['spatial']/self.eV:.3e} eV/s/m^3")
    print(f"  Bounce-averaged:  {power['bounce_average']:.3e} W/m^3 = {power['bounce_average']/self.eV:.3e} eV/s/m^3")
    print(f"  Midplane-mapped:    {power['midplane_mapped']:.3e} W/m^3 = {power['midplane_mapped']/self.eV:.3e} eV/s/m^3")
    
    # Also compute average energy per particle
    fluxes = self.compute_M0_fluxes(z_coords, vpar_coords, mu_coords, results, n_workers)
    print("\nAverage energy per injected particle:")
    print(f"  Spatial source:    {power['spatial']/fluxes['spatial']:.3e} J = {(power['spatial']/fluxes['spatial'])/self.eV:.3e} eV")
    print(f"  Bounce-averaged:  {power['bounce_average']/fluxes['bounce_average']:.3e} J = {(power['bounce_average']/fluxes['bounce_average'])/self.eV:.3e} eV")
    print(f"  Midplane-mapped:    {power['midplane_mapped']/fluxes['midplane_mapped']:.3e} J = {(power['midplane_mapped']/fluxes['midplane_mapped'])/self.eV:.3e} eV")
    
    return power
  
  def plot_comparison_2d(self, z_coords, vpar_coords, mu_coords, results=None, mu_idx=0, n_workers=None):
    """Plot 2D comparison of all three source models.
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))
    
    Z = results['Z'][:,:,mu_idx]
    VPAR = results['VPAR'][:,:,mu_idx] / 1e6
    
    im1 = ax1.pcolormesh(Z, VPAR, results['spatial'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im1, ax=ax1)
    ax1.set_xlabel('z (m)')
    ax1.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax1.set_title('Spatial beam source')
    
    im2 = ax2.pcolormesh(Z, VPAR, results['bounce_average'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im2, ax=ax2)
    ax2.set_xlabel('z (m)')
    ax2.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax2.set_title('Bounce-averaged spatial source')
    
    im3 = ax3.pcolormesh(Z, VPAR, results['midplane_mapped'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im3, ax=ax3)
    ax3.set_xlabel('z (m)')
    ax3.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax3.set_title('Midplane-mapped source')
    
    plt.tight_layout()
    return fig
  
  def optimize_midplane_map_to_match_bounce_avg(self, z_coords, vpar_coords, mu_coords,
                                          n_workers=None, verbose=True):
    """
    Optimize orbit_avg_beam_source parameters to match bounce_average_beam_source.
    
    Creates a NEW BeamSourceComparison object with optimized parameters (E_beam, 
    gamma_midmap, T_beam) such that midplane_mapped_beam_source best matches the 
    bounce_average_beam_source of the ORIGINAL (self) spatial source.
    
    This optimizes for the FULL 3D integrals over z, vpar, and mu, matching
    the total M0 flux and power computed by compute_M0_fluxes() and compute_power().
    
    The original object (self) is NOT modified.
    
    Args:
      z_coords: Array of z coordinates for the 3D grid
      vpar_coords: Array of vpar coordinates for the 3D grid
      mu_coords: Array of mu coordinates for the 3D grid
      n_workers: Number of parallel workers for bounce average
      verbose: Print optimization progress
    
    Returns:
      BeamSourceComparison: New object with optimized parameters
    """
    n_z = len(z_coords)
    n_vpar = len(vpar_coords)
    n_mu = len(mu_coords)
    
    if verbose:
      print(f"Computing bounce average on 3D grid ({n_z} x {n_vpar} x {n_mu} = {n_z*n_vpar*n_mu} points)...")
    
    # Create 3D grids
    Z_grid, VPAR_grid, MU_grid = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    
    # Compute bounce average for all points (this is expensive but only done once)
    bnc_avg_target_3d = self._compute_bounce_average_parallel(
      Z_grid.reshape(-1, 1, 1), 
      VPAR_grid.reshape(-1, 1, 1), 
      MU_grid.reshape(-1, 1, 1), 
      n_workers
    ).reshape(n_z, n_vpar, n_mu)
    
    # Precompute target M0 and power using class helper methods (consistent with compute_M0_fluxes/compute_power)
    M0_bnc_avg_total = self._integrate_M0_from_3d(bnc_avg_target_3d, z_coords, vpar_coords, mu_coords)
    power_bnc_avg_total = self._integrate_power_from_3d(bnc_avg_target_3d, Z_grid, VPAR_grid, MU_grid,
                                                        z_coords, vpar_coords, mu_coords)
    
    if verbose:
      print(f"  Target M0 flux (3D integral): {M0_bnc_avg_total:.3e}")
      print(f"  Target power (3D integral):   {power_bnc_avg_total:.3e} W = {power_bnc_avg_total/self.eV:.3e} eV/s")
    
    # Initial parameter guesses in physical units (eV)
    E_beam_init_eV = self.E_beam_midmap / self.eV
    T_beam_init_eV = self.T_beam_midmap / self.eV
    x0 = np.array([E_beam_init_eV, self.gamma_midmap, T_beam_init_eV])
    
    # Parameter bounds
    bounds = [
      (0.5 * E_beam_init_eV, 2.0 * E_beam_init_eV),
      (0.1 * self.gamma_midmap, 10.0 * self.gamma_midmap),
      (0.1 * T_beam_init_eV, 10.0 * T_beam_init_eV),
    ]
    
    def objective(params):
      """Objective function: combined residual of shape, M0 flux, and power (3D integrals)."""
      E_beam_test_eV, gamma_midmap_test, T_beam_test_eV = params
      
      # Create temporary object with test parameters
      temp_obj = BeamSourceComparison(
        mcB=self.mcB, gamma_geo=self.gamma_geo, Z_m=self.Z_m, B_p=self.B_p,
        gamma_midmap=gamma_midmap_test,
        T_beam_eV_midmap=T_beam_test_eV,
        E_beam_eV_midmap=E_beam_test_eV,
        use_potential=self._use_potential,
      )
      
      # Compute orbit_avg on 3D grid (this is fast - no contour integration needed)
      midmap_3d = temp_obj.midplane_mapped_beam_source(Z_grid, VPAR_grid, MU_grid)
      
      # 1. Shape residual (normalized by max value)
      max_val = max(np.max(np.abs(bnc_avg_target_3d)), 1e-30)
      shape_residual = np.sum(((midmap_3d - bnc_avg_target_3d) / max_val)**2) / (n_z * n_vpar * n_mu)
      
      # 2. M0 flux residual using compute_M0_orbit_avg (consistent with compute_M0_fluxes)
      M0_midmap_total = temp_obj.compute_M0_midplane_mapped(z_coords, vpar_coords, mu_coords)
      M0_ref = max(abs(M0_bnc_avg_total), 1e-30)
      M0_residual = ((M0_midmap_total - M0_bnc_avg_total) / M0_ref)**2
      
      # 3. Power residual using compute_power_orbit_avg (consistent with compute_power)
      power_midmap_total = temp_obj.compute_power_midplane_mapped(z_coords, vpar_coords, mu_coords)
      power_ref = max(abs(power_bnc_avg_total), 1e-30)
      power_residual = ((power_midmap_total - power_bnc_avg_total) / power_ref)**2
      
      total_residual = shape_residual + M0_residual + power_residual
      
      return total_residual, shape_residual, M0_residual, power_residual
    
    def objective_scalar(params):
      return objective(params)[0]
    
    if verbose:
      print("Optimizing orbit-averaged parameters...")
      print(f"  Initial: E_beam_midmap={x0[0]:.1f} eV, gamma_midmap={x0[1]:.3e}, T_beam_midmap={x0[2]:.1f} eV")
      init_total, init_shape, init_M0, init_power = objective(x0)
      print(f"  Initial residuals: shape={init_shape:.3e}, M0={init_M0:.3e}, power={init_power:.3e}, total={init_total:.3e}")
    
    result = minimize(objective_scalar, x0, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': 1000, 'ftol': 1e-16, 'gtol': 1e-16})
    
    E_beam_opt_eV, gamma_midmap_opt, T_beam_opt_eV = result.x
    
    if verbose:
      print(f"  Optimized: E_beam_midmap={E_beam_opt_eV:.1f} eV, gamma_midmap={gamma_midmap_opt:.3e}, T_beam_midmap={T_beam_opt_eV:.1f} eV")
      print(f"  Optimization converged: {result.success}, iterations: {result.nit}, message: {result.message}")
      final_total, final_shape, final_M0, final_power = objective(result.x)
      print(f"  Final residuals: shape={final_shape:.3e}, M0={final_M0:.3e}, power={final_power:.3e}, total={final_total:.3e}")
    
    if verbose:
      print(f"\nOptimized orbit-averaged beam parameters:")
      print(f"  E_beam_midmap: {self.E_beam_midmap/self.eV:.1f} eV -> {E_beam_opt_eV:.1f} eV")
      print(f"  T_beam_midmap: {self.T_beam_midmap/self.eV:.1f} eV -> {T_beam_opt_eV:.1f} eV")
      print(f"  gamma_midmap: {self.gamma_midmap:.3e} -> {gamma_midmap_opt:.3e}")
      print(f"\nSpatial source parameters (unchanged):")
      print(f"  E_beam_spatial: {self.E_beam_spatial/self.eV:.1f} eV")
      print(f"  T_beam_spatial: {self.T_beam_spatial/self.eV:.1f} eV")
      print(f"  gamma_spatial: {self.gamma_spatial:.3e}")
    
    # Create new object with optimized parameters
    optimized = BeamSourceComparison(
      mcB=self.mcB,
      gamma_geo=self.gamma_geo,
      Z_m=self.Z_m,
      B_p=self.B_p,
      gamma_spatial=self.gamma_spatial,
      gamma_midmap=gamma_midmap_opt,
      T_beam_eV_spatial=self.T_beam_spatial / self.eV,
      T_beam_eV_midmap=T_beam_opt_eV,
      E_beam_eV_spatial=self.E_beam_spatial / self.eV,
      E_beam_eV_midmap=E_beam_opt_eV,
      Lb=self.Lb,
      use_potential=self._use_potential
    )
    
    return optimized
  
  def plot_vpar_cut(self, mu, z_cut=0.0, n_workers=None):
    """Plot 1D comparison at fixed z, with x-axis in energy (keV)."""
    vpar_coords = np.linspace(0, 2e6, 500)
    midmap = self.midplane_mapped_beam_source(z_cut, vpar_coords, mu)
    
    # Use parallel computation for contour integral
    Z = np.full_like(vpar_coords, z_cut)
    MU = np.full_like(vpar_coords, mu)
    bnc_avg = self._compute_bounce_average_parallel(
      Z.reshape(-1, 1, 1), vpar_coords.reshape(-1, 1, 1), 
      MU.reshape(-1, 1, 1), n_workers).flatten()

    # Compute the integral for both sources
    bnc_avg_integral = np.trapz(bnc_avg, vpar_coords)
    midmap_integral = np.trapz(midmap, vpar_coords)
    print(f"At z={z_cut} m, mu={mu:.3e}:")
    print(f"  Bounce averaged M0 flux: {bnc_avg_integral:.3e}")
    print(f"  Midplane-mapped M0 flux:   {midmap_integral:.3e}")
    
    # Convert vpar to energy in keV: E = 0.5 * m * vpar^2
    # Use signed energy to preserve direction information
    energy_keV = np.sign(vpar_coords) * 0.5 * self.mi * vpar_coords**2 / (1000 * self.eV) \
                                      + mu * self.bmag(z_cut) / (1000 * self.eV) \
                                      + self.qi * self.phi(z_cut) / (1000 * self.eV)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(energy_keV, bnc_avg, label='Bounce Averaged', linestyle='-')
    ax.plot(energy_keV, midmap, label='Midplane Mapped', linestyle='--')
    ax.set_xlabel(r'$E =\frac{1}{2} m v_\parallel^2 + \mu B $ (keV)')
    ax.set_ylabel('Source')
    ax.set_title(f'Source comparison at z = {z_cut} m, mu = {mu:.3e} J/T')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return fig

  def plot_z_slice(self, z_cut=0.0, n_vpar=100, n_mu=100, n_workers=None):
    """Plot 2D comparison at fixed z: bounce averaged,midplane-mapped, and difference.
    
    Args:
      z_cut: z position for the slice
      n_vpar: Number of vpar grid points
      n_mu: Number of mu grid points
      n_workers: Number of parallel workers for contour integral
    
    Returns:
      fig: Matplotlib figure with 3 panels
    """
    # Create vpar and mu grids
    vpar_coords = np.linspace(0, 2e6, n_vpar)
    mu_max = 2.0 * 0.5 * self.mi * self.v_beam_spatial**2 / self.B_p
    mu_coords = np.linspace(0, mu_max, n_mu)
    
    VPAR, MU = np.meshgrid(vpar_coords, mu_coords, indexing='ij')
    Z = np.full_like(VPAR, z_cut)
    
    # Compute orbit-averaged source (fast, vectorized)
    midmap = self.midplane_mapped_beam_source(Z, VPAR, MU)
    
    # Compute contour integral (parallel)
    orbt_avg = self._compute_bounce_average_parallel(
      Z.reshape(-1, 1, 1), 
      VPAR.reshape(-1, 1, 1), 
      MU.reshape(-1, 1, 1), 
      n_workers
    ).reshape(n_vpar, n_mu)
    
    # Convert vpar to 10^6 m/s for plotting
    vpar_plot = VPAR / 1e6
    # Use mu directly for y-axis (in units of J/T)
    mu_plot = MU
    
    # Compute difference
    diff = midmap - orbt_avg
    
    # Create figure with 3 panels
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))
    
    # Common color scale for contour and orbit_avg
    vmax = max(np.max(orbt_avg), np.max(midmap))
    vmin = 0
    
    # Left panel: Contour integral
    im1 = ax1.pcolormesh(vpar_plot, mu_plot, orbt_avg, shading='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(im1, ax=ax1, label='Source')
    ax1.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax1.set_ylabel(r'$\mu$ (J/T)')
    ax1.set_title(f'Bounce Averaged at z = {z_cut} m')
    
    # Middle panel: Orbit-averaged
    im2 = ax2.pcolormesh(vpar_plot, mu_plot, midmap, shading='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(im2, ax=ax2, label='Source')
    ax2.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax2.set_ylabel(r'$\mu$ (J/T)')
    ax2.set_title(f'Midplane Mapped at z = {z_cut} m')
    
    # Right panel: Difference (midmap - orbt_avg)
    diff_max = max(abs(np.min(diff)), abs(np.max(diff)))
    im3 = ax3.pcolormesh(vpar_plot, mu_plot, diff, shading='auto', cmap='RdBu_r', 
                         vmin=-diff_max, vmax=diff_max)
    plt.colorbar(im3, ax=ax3, label='Difference')
    ax3.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax3.set_ylabel(r'$\mu$ (J/T)')
    ax3.set_title(f'Midplane map - Bounce Averaged at z = {z_cut} m')
    
    plt.tight_layout()
    return fig

# ============================================================================
# Script execution
# ============================================================================

if __name__ == "__main__":
  bsc = BeamSourceComparison(use_potential=True)
  
  # z_coords = np.linspace(-0.6, 0.6, 100)
  # vpar_coords = np.linspace(-1.5e6, 1.5e6, 100)
  z_coords = np.linspace(-0.8, 0.8, 50)
  vpar_coords = np.linspace(-2e6, 2e6, 50)
  mu_val = 0.5 * bsc.mi * bsc.v_beam_spatial**2 / bsc.B_p
  mu_coords = np.linspace(0.0, 2.0 * mu_val, 50)
  
  # Plot energy diagnostics to understand the potential effects
  # fig_diag = bsc.plot_energy_diagnostics()
  
  # results = bsc.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords)
  # power = bsc.compute_power(z_coords, vpar_coords, mu_coords, results)
  # fluxes = bsc.compute_M0_fluxes(z_coords, vpar_coords, mu_coords, results)
  
  # # Plot 2D comparison
  # fig1 = bsc.plot_comparison_2d(z_coords, vpar_coords, mu_val)
  # plt.show()
  
  # # Plot 1D vpar cut at z=0 (original)
  # fig2 = bsc.plot_vpar_cut(mu_val, z_cut=0.0)
  # fig2.suptitle('Before optimization')
  
  # # Plot 2D slice at z=0 (before optimization)
  # fig_slice = bsc.plot_z_slice(z_cut=0.0, n_vpar=100, n_mu=100)
  # fig_slice.suptitle('Before optimization')
  
  # Optimize using full 3D integrals
  bsc_optimized = bsc.optimize_midplane_map_to_match_bounce_avg(z_coords, vpar_coords, mu_coords)
  
  # mu_val = 4e-15
  # fig3 = bsc_optimized.plot_vpar_cut(mu_val, z_cut=0.0)

  bsc_optimized.compute_power(z_coords, vpar_coords, mu_coords)
  # fig3 = bsc_optimized.plot_z_slice(z_cut=0.0, n_vpar=100, n_mu=100)
  # fig3.suptitle('After optimization')

  # bsc_optimized.plot_comparison_2d(z_coords, vpar_coords, mu_val)
  plt.show()