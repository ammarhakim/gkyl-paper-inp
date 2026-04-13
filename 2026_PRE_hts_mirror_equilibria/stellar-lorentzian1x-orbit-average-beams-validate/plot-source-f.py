import postgkyl as pg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as multip
from concurrent.futures import ProcessPoolExecutor
import time
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.patches import FancyBboxPatch

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

read_frame = 100

directory_POA_simulation = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average'

# Hardcode a few tricky to read values from the data
nu_ii = 6.706930e+00 # Read from the simulation at midplane

# Universal constants (matching input_file.c)
eps0 = 8.8541878128e-12  # F/m (permittivity of free space)
mu0 = 4 * np.pi * 1e-7  # H/m (permeability of free space)
eV = 1.602176634e-19  # J (elementary charge / electron volt)
mp = 1.67262192369e-27  # kg (proton mass)
me = 9.1093837015e-31  # kg (electron mass)

# Plasma parameters (from input_file.c)
mi = 2.014 * mp  # kg (deuteron mass)
qi = eV  # C (ion charge)
qe = -eV  # C (electron charge)
Te0 = 940 * eV  # J (electron temperature)
n0 = 3e+19 # Read from the simulation at midplane
B_p = 0.53  # T (magnetic field strength)
beta = 0.4  # plasma beta
tau = (B_p**2) * beta / (2.0 * mu0 * n0 * Te0) - 1.0  # temperature ratio
Ti0 = tau * Te0  # J (ion temperature)
kperpRhos = 0.1  # normalized perpendicular wavenumber

# Derived parameters
vti = np.sqrt(Ti0 / mi)  # m/s (ion thermal velocity)
vte = np.sqrt(Te0 / me)  # m/s (electron thermal velocity)
c_s = np.sqrt(Te0 / mi)  # m/s (sound speed)
omega_ci = eV * B_p / mi  # rad/s (ion cyclotron frequency)
rho_s = c_s / omega_ci  # m (sound gyroradius)
kperp = kperpRhos / rho_s  # m^-1 (perpendicular wavenumber)

# Collision frequencies
nuFrac = 1.0
elc_nuFrac = 1/5.489216862238348
logLambdaElc = 6.6 - 0.5 * np.log(n0 / 1e20) + 1.5 * np.log(Te0 / eV)
nuElc = elc_nuFrac * nuFrac * logLambdaElc * (eV**4) * n0 / (6 * np.sqrt(2) * \
                        (np.pi**(3/2)) * (eps0**2) * np.sqrt(me) * (Te0**(3/2)))

# Geometry parameters
z_min = -2.5
z_max = 2.5
psi_min = 1e-6
psi_eval = 1e-3
psi_max = 3e-3

# Velocity space parameters
vpar_max_elc = 16 * vte
mu_max_elc = me * (3 * vte)**2 / (2 * B_p)
vpar_max_ion = 16 * vti
mu_max_ion = mi * (3 * vti)**2 / (2 * B_p)

# Source parameters
ion_source_amplitude = 1e20  # m^-3/s
ion_source_sigma = 0.5
ion_source_temp = 5000 * eV  # J

def plot_biMax():
  # Plot BiMaxwellian Moments
  data_BiMax = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_BiMaxwellianMoments_'+str(read_frame)+'.gkyl')
  data_BiMax.info()
  interp_BiMax = pg.GInterpModal(data_BiMax, 1, 'ms')
  x, dens = interp_BiMax.interpolate(0)
  x, upar = interp_BiMax.interpolate(1)
  x, Tpar_div_m = interp_BiMax.interpolate(2)
  x, Tperp_div_m = interp_BiMax.interpolate(3)

  x = np.squeeze(x)
  dens = np.squeeze(dens)
  upar = np.squeeze(upar)
  Tpar_div_m = np.squeeze(Tpar_div_m)
  Tperp_div_m = np.squeeze(Tperp_div_m)

  Tpar_ev = Tpar_div_m * mi / eV
  Tperp_ev = Tperp_div_m * mi / eV

  data_mc2p = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-mapc2p.gkyl')
  interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
  nodes_Z = interp_mc2p.interpolate(2)[1]
  nodes_Z = np.squeeze(nodes_Z)

  # Calculate global x-limits from all position maps
  x_min = nodes_Z.min()
  x_max = nodes_Z.max()

  # Function to split data at x=0
  def split_at_zero(x, y):
      """Split data into left (x<0) and right (x>=0) portions"""
      # Ensure x is 1D for comparison
      if x.ndim > 1:
          x_1d = x.flatten()
      else:
          x_1d = x
      
      # Ensure y is 1D for indexing
      if y.ndim > 1:
          y_1d = y.flatten()
      else:
          y_1d = y
      
      left_mask = x_1d < 0
      right_mask = x_1d >= 0
      return x_1d[left_mask], y_1d[left_mask], x_1d[right_mask], y_1d[right_mask]

  # Split data at zero
  dens_x_left, dens_y_left, dens_x_right, dens_y_right = split_at_zero(nodes_Z, dens)
  upar_x_left, upar_y_left, upar_x_right, upar_y_right = split_at_zero(nodes_Z, upar/c_s)
  Tpar_x_left, Tpar_y_left, Tpar_x_right, Tpar_y_right = split_at_zero(nodes_Z, Tpar_ev)
  Tperp_x_left, Tperp_y_left, Tperp_x_right, Tperp_y_right = split_at_zero(nodes_Z, Tperp_ev)

  fig, ax = plt.subplots(2, 4, figsize=(10, 6))

  # Density - Left panel (linear scale)
  ax[0, 0].plot(dens_x_left, dens_y_left, color='blue')
  ax[0, 0].plot([-0.98, -0.98], [0, 4e19], color='grey', linestyle='--')
  ax[0, 0].text(-2.3, 2e19, r'$\mathbf{(a)}$')
  ax[0, 0].set_ylabel(r'$n_i$ (m$^{-3}$)')
  ax[0, 0].set_xlim(x_min, 0)
  ax[0, 0].set_ylim(0, 2.2e19)

  # Density - Right panel (log scale)
  ax[0, 1].plot(dens_x_right, dens_y_right, color='blue')
  ax[0, 1].plot([0.98, 0.98], [1e13, 1e20], color='grey', linestyle='--')
  ax[0, 1].set_yscale('log')
  ax[0, 1].yaxis.tick_right()
  ax[0, 1].yaxis.set_label_position("right")
  ax[0, 1].set_xlim(0, x_max)
  ax[0, 1].set_ylim(1e13, 4e19)

  # Upar - Left panel
  ax[0, 2].plot([-0.98, -0.98], [-10, 10], color='grey', linestyle='--')
  ax[0, 2].plot(upar_x_left, upar_y_left, color='blue')
  ax[0, 2].set_ylabel(r'$u_{||} / c_s$')
  ax[0, 2].text(-2.3, 6, r'$\mathbf{(b)}$')
  ax[0, 2].set_ylim(-7, 7)
  ax[0, 2].set_xlim(x_min, 0)

  # Upar - Right panel
  ax[0, 3].plot([0.98, 0.98], [-10, 10], color='grey', linestyle='--')
  ax[0, 3].plot(upar_x_right, upar_y_right, color='blue')
  ax[0, 3].yaxis.tick_right()
  ax[0, 3].set_yscale('log')
  ax[0, 3].yaxis.set_label_position("right")
  ax[0, 3].set_ylim(0.001, 10)
  ax[0, 3].set_xlim(0, x_max)

  # Tpar - Left panel
  ax[1, 0].plot([-0.98, -0.98], [0, 15000], color='grey', linestyle='--')
  ax[1, 0].plot(Tpar_x_left, Tpar_y_left, color='blue')
  ax[1, 0].set_ylabel(r'$T_{||}$ (eV)')
  ax[1, 0].text(-2.3, 10000, r'$\mathbf{(c)}$')
  ax[1, 0].set_ylim(0, 11000)
  ax[1, 0].set_xlim(x_min, 0)

  # Tpar - Right panel
  ax[1, 1].plot([0.98, 0.98], [0, 150000], color='grey', linestyle='--')
  ax[1, 1].plot(Tpar_x_right, Tpar_y_right, color='blue')
  ax[1, 1].yaxis.tick_right()
  ax[1, 1].set_yscale('log')
  ax[1, 1].yaxis.set_label_position("right")
  ax[1, 1].set_ylim(0, 13000)
  ax[1, 1].set_xlim(0, x_max)

  # Tperp - Left panel
  ax[1, 2].plot([-0.98, -0.98], [0, 35000], color='grey', linestyle='--')
  ax[1, 2].plot(Tperp_x_left, Tperp_y_left, color='blue')
  ax[1, 2].set_ylabel(r'$T_{\perp}$ (eV)')
  ax[1, 2].set_ylim(0, 24000)
  ax[1, 2].text(-2.3, 22000, r'$\mathbf{(d)}$')
  ax[1, 2].set_xlim(x_min, 0)

  # Tperp - Right panel
  ax[1, 3].plot([0.98, 0.98], [0, 350000], color='grey', linestyle='--')
  ax[1, 3].plot(Tperp_x_right, Tperp_y_right, color='blue')
  ax[1, 3].yaxis.tick_right()
  ax[1, 3].set_yscale('log')
  ax[1, 3].yaxis.set_label_position("right")
  ax[1, 3].set_ylim(0, 100000)
  ax[1, 3].set_xlim(0, x_max)

  # # Add a grid to each subplot
  # for a in ax.flat:
  #   a.grid(True, alpha=0.3)

  xlabel_1 = ax[1, 1].set_xlabel(r'$Z$ (m)')
  ax[1, 1].xaxis.set_label_coords(-0, -0.15)
  # For the right pair (columns 2 and 3), place xlabel on column 3
  xlabel_3 = ax[1, 3].set_xlabel(r'$Z$ (m)')
  ax[1, 3].xaxis.set_label_coords(-0, -0.15)

  plt.tight_layout()
    
  # Double all subplot widths
  for i in range(2):  # rows
      for j in range(4):  # columns
          pos = ax[i, j].get_position()
          new_width = pos.width * 1.5
          new_pos = [pos.x0, pos.y0, new_width, pos.height]
          ax[i, j].set_position(new_pos)

  # Manually adjust subplot positions to control spacing
  # Get current positions
  for i in range(2):  # rows
      for j in range(4):  # columns
          pos = ax[i, j].get_position()
          
          # Adjust horizontal positions
          if j == 0:  # First column - shift slightly right
              new_pos = [pos.x0 + 0.02, pos.y0, pos.width, pos.height]
          elif j == 1:  # Second column - no left margin (merge with first)
              new_pos = [ax[i, 0].get_position().x0 + ax[i, 0].get_position().width, pos.y0, pos.width, pos.height]
          elif j == 2:  # Third column - normal position with gap
              new_pos = [pos.x0, pos.y0, pos.width, pos.height]
          elif j == 3:  # Fourth column - no left margin (merge with third)
              new_pos = [ax[i, 2].get_position().x0 + ax[i, 2].get_position().width, pos.y0, pos.width, pos.height]
          
          ax[i, j].set_position(new_pos)

  plt.savefig('bimaxwellian_moments_lorentzian.pdf')
  # plt.show()
  plt.close()  # Close figure to free memory

def plot_field():
  data_0 = pg.GData(directory_POA_simulation + '/field_time_trace_z0_eq_0.gkyl')
  data_098 = pg.GData(directory_POA_simulation + '/field_time_trace_z0_eq_0,98.gkyl')
  data_25 = pg.GData(directory_POA_simulation + '/field_time_trace_z0_eq_2,5.gkyl')

  t = data_0.get_grid()
  t = np.array(t[0])[:-1]

  phi_0 = np.array(data_0.get_values())[:,0,0] * -qe / Te0
  phi_098 = np.array(data_098.get_values())[:,0,0] * -qe / Te0
  phi_25 = np.array(data_25.get_values())[:,0,0] * -qe / Te0

  # print('Potential drop from midplane to mirror: ' + str(phi_0[-1] - phi_098[-1]))
  # print('Potential at z=0.98 to edge: ' + str(phi_098[-1]))
  # Get the final field
  data_phi = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-field_'+str(read_frame)+'.gkyl')
  interp_phi = pg.GInterpModal(data_phi, 1, 'ms')
  x, phi_final = interp_phi.interpolate(0)

  data_mc2p = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-mapc2p.gkyl')
  interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
  nodes_Z = interp_mc2p.interpolate(2)[1]
  nodes_Z = np.squeeze(nodes_Z)

  fig, ax = plt.subplots( 1, 2, figsize=(10,4))

  ax[0].plot(t*nu_ii, phi_0, label=r'$Z=0.0$ m')
  # ax[0].plot(t, phi_05, label='z=0.5')
  ax[0].plot(t*nu_ii, phi_098, label=r'$Z=0.98$ m')
  # ax[0].plot(t, phi_15, label='z=1.5')
  ax[0].plot(t*nu_ii, phi_25, label=r'$Z=2.5$ m')

  # Make a brace betwen the final value of phi_0 and phi_1
  # Draw a vertical bracket between phi_0[-1] and phi_1[-1] at t[-1]
  y0, y1 = phi_098[-1], phi_0[-1]
  x = t[-1] * nu_ii
  # Draw the bracket using annotate with arrowstyle '-['
  ax[0].annotate(
      '', 
      xy=((x-3), (y0 + y1)/2), 
      xytext=((x - 0.1*(t[-1] - t[0]) - 15), (y0 + y1)/2), 
      arrowprops=dict(arrowstyle='<-[, lengthB=0.5, widthB=3.9', lw=1, color='black'), 
      annotation_clip=False
  )
  # Draw a grey rectangle (box) around the text

  # Box position and size
  box_x = (x - 0.1*(t[-1] - t[0]) - 0.2 - 10*nu_ii)  # shift left a bit more for padding
  box_y = (y0 + y1)/2 - 1.1  # lower edge of box
  box_width = 7.5*nu_ii  # enough for both lines
  box_height = 2.2

  # Add the rectangle patch
  rect = FancyBboxPatch(
      (box_x, box_y),
      box_width,
      box_height,
      boxstyle="round,pad=0.2",
      linewidth=1,
      edgecolor='grey',
      facecolor='lightgrey',
      alpha=0.5,
      zorder=2
  )
  ax[0].add_patch(rect)

  # Place the text to the left of the bracket, vertically centered, inside the box
  ax[0].text(
      (x - 0.1*(t[-1] - t[0]) - 18), 
      (y0 + y1)/2 + 0.5, 
      r'$\Delta e \phi / T_e = {:.2f}$'.format(phi_0[-1] - phi_098[-1]), 
      fontsize=12, 
      ha='right', 
      va='center',
      zorder=3
  )
  ax[0].text(
      (x - 0.1*(t[-1] - t[0]) - 18), 
      (y0 + y1)/2 - 0.7, 
      r'$\Delta e\phi_{{\mathrm{{theory}}}} / T_e = {:.2f}$'.format(7.31), 
      fontsize=12, 
      ha='right', 
      va='center',
      zorder=3
  )

  # Make another bracket between the final value of phi_1 and phi_19

  y0, y1 = phi_25[-1], phi_098[-1]
  ax[0].annotate(
    '',
    xy=((x-3), (y0 + y1)/2),
    xytext=((x - 0.1*(t[-1] - t[0]) - 15), (y0 + y1)/2),
    arrowprops=dict(arrowstyle='<-[, lengthB=0.5, widthB=3.5', lw=1, color='black'),
    annotation_clip=False
  )
  # Draw a grey rectangle (box) around the text
  # Box position and size
  box_x = (x - 0.1*(t[-1] - t[0]) - 0.2 - 10*nu_ii)  # shift left a bit more for padding
  box_y = (y0 + y1)/2 - 1.1  # lower edge of box
  box_width = 7.5*nu_ii  # enough for both lines
  box_height = 2.2

  # Add the rectangle patch
  rect = FancyBboxPatch(
      (box_x, box_y),
      box_width,
      box_height,
      boxstyle="round,pad=0.2",
      linewidth=1,
      edgecolor='grey',
      facecolor='lightgrey',
      alpha=0.5,
      zorder=2
  )
  ax[0].add_patch(rect)

  # Place the text to the left of the bracket, vertically centered, inside the box
  ax[0].text(
      (x - 0.1*(t[-1] - t[0]) - 18), 
      (y0 + y1)/2 + 0.5, 
      r'$\Delta e \phi / T_e = {:.2f}$'.format(phi_098[-1] - phi_25[-1]), 
      fontsize=12, 
      ha='right', 
      va='center',
      zorder=3
  )
  ax[0].text(
      (x - 0.1*(t[-1] - t[0]) - 18), 
      (y0 + y1)/2 - 0.7, 
      r'$\Delta e\phi_{{\mathrm{{theory}}}} / T_e = {:.2f}$'.format(6.32), 
      fontsize=12, 
      ha='right', 
      va='center',
      zorder=3
  )
  

  ax[0].set_xlabel(r'$t \nu_{ii}$')
  ax[0].set_ylabel(r'$e\phi / T_e$')
  # ax[0].set_title(r'Time trace of $\phi$ at various $z$ locations')
  ax[0].set_xlim(-0.5*nu_ii, t[-1]*nu_ii)
  # ax[0].set_ylim(0, 1.1*np.max(phi_0))
  ax[0].legend()
  # ax[0].grid()


  ax[1].plot([0.98, 0.98], [0, 16], color='grey', linestyle='--', label=r'$Z=0.98 \mathrm{m}$')
  ax[1].plot([-0.98, -0.98], [0, 16], color='grey', linestyle='--', label=r'$Z=-0.98 \mathrm{m}$')
  ax[1].plot(nodes_Z, -phi_final * qe / Te0, label='Final Field', color='blue')
  ax[1].set_xlabel(r'$Z$, (m)')
  # ax[1].set_title('Potential at t = {:.2f} s'.format(t[-1]))
  ax[1].set_xlim(z_min, z_max)
  ax[1].set_ylim(0, 15)

  ax[0].text(45,12, r'$\mathbf{(a)}$')
  ax[1].text(-2.3,14, r'$\mathbf{(b)}$')

  plt.tight_layout()
  plt.savefig('field_fig_lorentzian.pdf')
  # plt.show()
  plt.close()  # Close figure to free memory

def plot_integrated_moments():
  data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_integrated_moms.gkyl')

  t = data.get_grid()
  t = np.array(t[0])[:-1] * nu_ii

  integ_moms = np.array(data.get_values())[:-1,:]

  M0 = integ_moms[:,0]  # Total ion density
  M1 = integ_moms[:,1]  # Parallel ion momentum
  M2par = integ_moms[:,2]  # Parallel ion temperature
  M2perp = integ_moms[:,3]  # Perpendicular ion temperature
  
  fig, ax = plt.subplots(2, 2, figsize=(10, 6))

  ax[0, 0].plot(t, M0, label='Total Ion Density', color='blue')
  ax[0, 0].set_ylabel(r'$\left< M_0 \right>$')
  # ax[0, 0].set_yscale('log')

  # Inset zoom for t=5 to t=10 on ax[0, 0]
  axins = inset_axes(ax[0, 0], width="40%", height="40%", loc='lower left', borderpad=2)
  axins.plot(t, M0, color='blue')
  axins.set_xlim(10, 20)
  axins.set_ylim(np.min(M0[(t >= 10) & (t <=20)])*0.995, np.max(M0[(t >= 10) & (t <= 20)]))
  axins.set_xticks([])
  axins.set_yticks([])
  mark_inset(ax[0, 0], axins, loc1=2, loc2=1, fc="none", ec="0.5")

  ax[0, 1].plot(t, M1, label='Parallel Ion Momentum', color='blue')
  ax[0, 1].set_ylabel(r'$\left< M_1 \right>$')

  ax[1, 0].plot(t, M2par, label='Parallel Ion Temperature', color='blue')
  ax[1, 0].set_ylabel(r'$\left< M_{2,\parallel} \right>$')

  ax[1, 1].plot(t, M2perp, label='Perpendicular Ion Temperature', color='blue')
  ax[1, 1].set_ylabel(r'$\left< M_{2,\perp} \right>$')


  # Set the x-axis label for the bottom row
  ax[1, 0].set_xlabel(r'$t \nu_{ii}$')
  ax[1, 1].set_xlabel(r'$t \nu_{ii}$')

  ax[0, 0].text(130,0.1e19, r'$\mathbf{(a)}$')
  ax[0, 1].text(130,-1.9e19, r'$\mathbf{(b)}$')
  ax[1, 0].text(130,0.1e31, r'$\mathbf{(c)}$')
  ax[1, 1].text(130,0.1e31, r'$\mathbf{(d)}$')

  plt.tight_layout()
  plt.savefig('integrated_moments_lorentzian.pdf')
  plt.close()  # Close figure to free memory

def plot_f_at_0():
  # Plot f at z=0
  f_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_'+str(read_frame)+'.gkyl',\
                     mapc2p_vel_name=directory_POA_simulation + '/gk_lorentzian_mirror-ion_mapc2p_vel.gkyl')
  Jv_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_jacobvel.gkyl')
  jb_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-jacobtot.gkyl')
  bmag_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-bmag.gkyl')
  phi_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-field_100.gkyl')
  f_c = f_data.get_values()
  Jv_c = Jv_data.get_values()
  f_data._values = f_c/Jv_c
  dg = pg.GInterpModal(f_data, 1, 'gkhyb')
  bmag_dg = pg.GInterpModal(bmag_data, 1, 'ms')
  jb_dg = pg.GInterpModal(jb_data, 1, 'ms')
  phi_dg = pg.GInterpModal(phi_data, 1, 'ms')
  xInt_i, fIon = dg.interpolate()
  xInt, bmag = bmag_dg.interpolate(0)
  xInt, phi = phi_dg.interpolate(0)
  XInt, jb = jb_dg.interpolate(0)
  fIon = np.squeeze(fIon)
  bmag = np.squeeze(bmag)
  phi = np.squeeze(phi)
  jb = np.squeeze(jb)

  data_mc2p = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-mapc2p.gkyl')
  interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
  nodes_Z = interp_mc2p.interpolate(2)[1]
  nodes_Z = np.squeeze(nodes_Z)

  # Set fIon = 1e-60 where it is <0
  fIon[fIon < 1e-16] = 1e-16

  ndim = len(xInt_i)
  nxInt_i = [np.size(xInt_i[d]) for d in range(ndim)]

  #[ Cell center coordinates
  xIntC_i = [np.zeros(np.size(xInt_i[d])) for d in range(ndim)]
  for d in range(len(xIntC_i)):
    xIntC_i[d] = 0.5*(xInt_i[d][:-1]+xInt_i[d][1:])

  nxIntC_i = [np.size(xIntC_i[d]) for d in range(ndim)]

  xInt_i[0] = nodes_Z

  #[ Get indices along z of slices we wish to plot:
  plot_z_locations = [0., 0.75, 0.9, 0.98, 1.1, 2.5]

  plotzIdx = [np.argmin(np.abs(nodes_Z-val)) for val in plot_z_locations]

  Bmin = bmag[plotzIdx[0]]
  B075 = bmag[plotzIdx[1]]
  B09 = bmag[plotzIdx[2]]
  Bmax = bmag[plotzIdx[3]]

  phi_center = phi[plotzIdx[0]]
  phi075 = phi[plotzIdx[1]]
  phi09 = phi[plotzIdx[2]]
  phi_mirror = phi[plotzIdx[3]]

  mu_loss = (( 1/2 * mi * xInt_i[1]**2 ) + qi * (phi_center - phi_mirror)) / (Bmax - Bmin)
  mu_loss   = mu_loss  /(0.5 * mi*(vti**2)/B_p)
  mu_loss_075 = (( 1/2 * mi * xInt_i[1]**2 ) + qi * (phi075 - phi_mirror)) / (Bmax - B075)
  mu_loss_075 = mu_loss_075 / (0.5 * mi*(vti**2)/B_p)
  mu_loss_09 = (( 1/2 * mi * xInt_i[1]**2 ) + qi * (phi09 - phi_mirror)) / (Bmax - B09)
  mu_loss_09 = mu_loss_09 / (0.5 * mi*(vti**2)/B_p)

  #[ Normalize velocity space
  xInt_i[1] = xInt_i[1]/vti
  xInt_i[2] = xInt_i[2]/(0.5*mi*(vti**2)/B_p)
  xIntC_i[1] = xIntC_i[1]/vti
  xIntC_i[2] = xIntC_i[2]/(0.5*mi*(vti**2)/B_p)


  #[ Create colorplot grid. Recall coordinates have to be nodal.
  Xnodal_i = [np.outer(xInt_i[1], np.ones(nxInt_i[2])),
              np.outer(np.ones(nxInt_i[1]), xInt_i[2])]
  
  extreme_vals = [0., np.amax(fIon[plotzIdx[0],:,:]/jb[plotzIdx[0]])]

  fig, ax = plt.subplots(2,3, figsize=(10,6))

  colorbar_zmin = 1e-9
  xlim = (-10, 10)

  pcm = ax[0,0].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[0],:,:]/jb[plotzIdx[0]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  ax[0,0].plot(xInt_i[1], mu_loss, color='white', linestyle='--', label='Loss Cone')
  # ax[0,0].plot(xInt_i[1], mu_loss*2, color='grey', linestyle='--', label='Loss Cone')
  ax[0,0].set_ylim(0, 9)
  ax[0,0].set_ylabel(r'$\mu / \mu_{ti,0}$')
  ax[0,0].set_title(r'Z = {:.1f} m'.format(np.abs(nodes_Z[plotzIdx[0]])))
  ax[0,0].set_xlim(xlim)


  pcm = ax[0,1].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[1],:,:]/jb[plotzIdx[1]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  ax[0,1].plot(xInt_i[1], mu_loss_075, color='white', linestyle='--', label='Loss Cone')
  # ax[0,1].plot(xInt_i[1], mu_loss_075*2, color='grey', linestyle='--', label='Loss Cone')
  ax[0,1].set_ylim(0, 9)
  ax[0,1].set_title(r'Z = {:.2f} m'.format(nodes_Z[plotzIdx[1]]))
  ax[0,1].set_xlim(xlim)

  pcm = ax[0,2].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[2],:,:]/jb[plotzIdx[2]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  ax[0,2].plot(xInt_i[1], mu_loss_09, color='white', linestyle='--', label='Loss Cone')
  # ax[0,2].plot(xInt_i[1], mu_loss_09*2, color='grey', linestyle='--', label='Loss Cone')
  ax[0,2].set_ylim(0, 9)
  cbar = plt.colorbar(pcm, ax=ax[0,2], orientation='vertical')
  cbar.set_label(r'$f_i$'.format(nodes_Z[plotzIdx[2]]), fontsize=22)
  cbar.ax.tick_params(labelsize=10)
  ax[0,2].set_title(r'Z = {:.1f} m'.format(nodes_Z[plotzIdx[2]]))
  ax[0,2].set_xlim(xlim)

  pcm = ax[1,0].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[3],:,:]/jb[plotzIdx[3]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  ax[1,0].set_ylabel(r'$\mu / \mu_{ti,0}$')
  ax[1,0].set_title(r'Z = {:.2f} m'.format(nodes_Z[plotzIdx[3]]))
  ax[1,0].set_xlim(xlim)

  pcm = ax[1,1].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[4],:,:]/jb[plotzIdx[4]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  ax[1,1].set_title(r'Z = {:.1f} m'.format(nodes_Z[plotzIdx[4]]))
  ax[1,1].set_xlim(xlim)

  pcm = ax[1,2].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[5],:,:]/jb[plotzIdx[5]],
                         cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
  cbar = plt.colorbar(pcm, ax=ax[1,2], orientation='vertical')
  cbar.set_label(r'$f_i$'.format(nodes_Z[plotzIdx[5]]), fontsize=22)
  cbar.ax.tick_params(labelsize=10)
  ax[1,2].set_title(r'Z = {:.1f} m'.format(nodes_Z[plotzIdx[5]]))
  ax[1,2].set_xlim(xlim)


  ax[1,1].set_xlabel(r'$v_{||}/v_{ti,0}$')

  ax[0,0].text(-9, 8 , r'$\mathbf{(a)}$',color='white')
  ax[0,1].text(-9, 8 , r'$\mathbf{(b)}$',color='white')
  ax[0,2].text(-9, 8 , r'$\mathbf{(c)}$',color='white')
  ax[1,0].text(-9, 8 , r'$\mathbf{(d)}$',color='white')
  ax[1,1].text(-9, 8 , r'$\mathbf{(e)}$',color='white')
  ax[1,2].text(-9, 8 , r'$\mathbf{(f)}$',color='white')

  # plt.show()
  plt.tight_layout()
  plt.savefig('f_distributions_lorentzian.pdf')
  plt.close()  # Close figure to free memory

def run_plotting_functions_parallel(functions):
  """
  Run a list of plotting functions in parallel.
  
  Parameters:
  functions: List of function objects to execute in parallel
  """
  print(f"Running {len(functions)} plotting functions in parallel...")
  start_time = time.time()

  # Ensure the number of workers does not exceed available CPU cores
  max_workers = min(len(functions), multip.cpu_count())

  # Use ProcessPoolExecutor for true parallelism (avoids GIL issues)
  with ProcessPoolExecutor(max_workers=max_workers) as executor:
    # Submit all functions to the executor
    futures = [executor.submit(func) for func in functions]
    
    # Wait for all functions to complete and collect results
    for i, future in enumerate(futures):
      try:
        future.result()  # This will raise an exception if the function failed
        print(f"  ✓ Function {i+1}/{len(functions)} completed successfully")
      except Exception as e:
        print(f"  ✗ Function {i+1}/{len(functions)} failed with error: {e}")
  
  end_time = time.time()
  print(f"All plotting functions completed in {end_time - start_time:.2f} seconds")

# Array of plotting functions to run in parallel
plotting_functions = [
  plot_biMax,
  plot_field,
  plot_integrated_moments,
  plot_f_at_0,
]

if __name__ == "__main__":
  # Run all plotting functions in parallel
  run_plotting_functions_parallel(plotting_functions)