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
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})

read_frame = 100

directory_POA_simulation = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average-exploration/stellar-lorentzian1x-orbit-average-time-dilation-unif-Nz800'
directory_hires_simulation = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average-high-res-no-positivity'

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

# Plot BiMaxwellian Moments
data_BiMax = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_BiMaxwellianMoments_'+str(read_frame)+'.gkyl')
interp_BiMax = pg.GInterpModal(data_BiMax, 1, 'ms')
x, dens = interp_BiMax.interpolate(0)
x, upar = interp_BiMax.interpolate(1)
x, Tpar_div_m = interp_BiMax.interpolate(2)
x, Tperp_div_m = interp_BiMax.interpolate(3)

data_BiMax_hires = pg.GData(directory_hires_simulation + '/gk_lorentzian_mirror-ion_BiMaxwellianMoments_'+str(read_frame)+'.gkyl')
interp_BiMax_hires = pg.GInterpModal(data_BiMax_hires, 1, 'ms')
x_hires, dens_hires = interp_BiMax_hires.interpolate(0)
x_hires, upar_hires = interp_BiMax_hires.interpolate(1)
x_hires, Tpar_div_m_hires = interp_BiMax_hires.interpolate(2)
x_hires, Tperp_div_m_hires = interp_BiMax_hires.interpolate(3)

x = np.squeeze(x)
dens = np.squeeze(dens)
upar = np.squeeze(upar)
Tpar_div_m = np.squeeze(Tpar_div_m)
Tperp_div_m = np.squeeze(Tperp_div_m)

dens_hires = np.squeeze(dens_hires)
upar_hires = np.squeeze(upar_hires)
Tpar_div_m_hires = np.squeeze(Tpar_div_m_hires)
Tperp_div_m_hires = np.squeeze(Tperp_div_m_hires)

Tpar_ev = Tpar_div_m * mi / eV
Tperp_ev = Tperp_div_m * mi / eV

Tpar_ev_hires = Tpar_div_m_hires * mi / eV
Tperp_ev_hires = Tperp_div_m_hires * mi / eV

data_mc2p = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-mapc2p.gkyl')
interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
nodes_Z = interp_mc2p.interpolate(2)[1]
nodes_Z = np.squeeze(nodes_Z)

data_mc2p_hires = pg.GData(directory_hires_simulation + '/gk_lorentzian_mirror-mapc2p.gkyl')
interp_mc2p_hires = pg.GInterpModal(data_mc2p_hires, 1, 'ms')
nodes_Z_hires = interp_mc2p_hires.interpolate(2)[1]
nodes_Z_hires = np.squeeze(nodes_Z_hires)


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

# Plot deflated position map
# Handle different possible shapes of the deflated position map
# Since nodes_Z is already squeezed to 1D, we can use it directly
u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(nodes_Z, nodes_Z)
u_x_left_hires, u_y_left_hires, u_x_right_hires, u_y_right_hires = split_at_zero(nodes_Z_hires, nodes_Z_hires)

# Split density data at zero
dens_x_left, dens_y_left, dens_x_right, dens_y_right = split_at_zero(nodes_Z, dens)
dens_x_left_hires, dens_y_left_hires, dens_x_right_hires, dens_y_right_hires = split_at_zero(nodes_Z_hires, dens_hires)

fig, ax = plt.subplots(2, 4, figsize=(10, 6))

# Density - Left panel (linear scale)
ax[0, 0].plot(dens_x_left, dens_y_left, label='Nonuniform', color='blue')
ax[0, 0].plot(dens_x_left_hires, dens_y_left_hires, label='Uniform', color='red', linestyle='--')
ax[0, 0].plot([-0.98, -0.98], [0, 4e19], color='grey', linestyle='--')
# ax[0, 0].legend(fontsize=9)
ax[0, 0].text(-2.3, 2e19, r'$\mathbf{(a)}$', fontsize=18)
ax[0, 0].set_ylabel(r'$n$ (m$^{-3}$)')
ax[0, 0].set_xlim(x_min, 0)
ax[0, 0].set_ylim(0, 2.2e19)

# Density - Right panel (log scale)
ax[0, 1].plot(dens_x_right, dens_y_right, color='blue')
ax[0, 1].plot(dens_x_right_hires, dens_y_right_hires, color='red', linestyle='--')
ax[0, 1].plot([0.98, 0.98], [1e13, 1e20], color='grey', linestyle='--')
ax[0, 1].set_yscale('log')
# ax[0, 1].set_ylabel(r'$n_i$ (m$^{-3}$)')
ax[0, 1].yaxis.tick_right()
ax[0, 1].yaxis.set_label_position("right")
ax[0, 1].set_xlim(0, x_max)
ax[0, 1].set_ylim(1e13, 4e19)

# Split upar data at zero
upar_x_left, upar_y_left, upar_x_right, upar_y_right = split_at_zero(nodes_Z, upar/c_s)
upar_x_left_hires, upar_y_left_hires, upar_x_right_hires, upar_y_right_hires = split_at_zero(nodes_Z_hires, upar_hires/c_s)

# Upar - Left panel
ax[0, 2].plot(upar_x_left, upar_y_left, color='blue')
ax[0, 2].plot(upar_x_left_hires, upar_y_left_hires, color='red', linestyle='--')
ax[0, 2].plot([-0.98, -0.98], [-10, 10], color='grey', linestyle='--')
ax[0, 2].text(-2.3, 2, r'$\mathbf{(b)}$', fontsize=18)
ax[0, 2].set_ylabel(r'$u_{||} / c_s$')
ax[0, 2].set_ylim(-7, 7)
ax[0, 2].legend(['Time dilation', 'Without'])
ax[0, 2].set_xlim(x_min, 0)

# Upar - Right panel
ax[0, 3].plot([0.98, 0.98], [-10, 10], color='grey', linestyle='--')
ax[0, 3].plot(upar_x_right, upar_y_right, color='blue')
ax[0, 3].plot(upar_x_right_hires, upar_y_right_hires, color='red', linestyle='--')
# ax[0, 3].set_ylabel(r'$u_{||} / c_s$')
ax[0, 3].yaxis.tick_right()
ax[0, 3].yaxis.set_label_position("right")
ax[0, 3].set_ylim(.001, 10)
ax[0, 3].set_yscale('log')
ax[0, 3].set_xlim(0, x_max)

# Split Tpar data at zero
Tpar_x_left, Tpar_y_left, Tpar_x_right, Tpar_y_right = split_at_zero(nodes_Z, Tpar_ev)
Tpar_x_left_hires, Tpar_y_left_hires, Tpar_x_right_hires, Tpar_y_right_hires = split_at_zero(nodes_Z_hires, Tpar_ev_hires)

# Tpar - Left panel
ax[1, 0].plot([-0.98, -0.98], [0, 15000], color='grey', linestyle='--')
ax[1, 0].plot(Tpar_x_left, Tpar_y_left, color='blue')
ax[1, 0].text(-2.3, 10000, r'$\mathbf{(c)}$')
ax[1, 0].plot(Tpar_x_left_hires, Tpar_y_left_hires, color='red', linestyle='--')
ax[1, 0].set_ylabel(r'$T_{||}$ (eV)')
ax[1, 0].set_ylim(0, 11000)
ax[1, 0].set_xlim(x_min, 0)

# Tpar - Right panel
ax[1, 1].plot([0.98, 0.98], [0, 15000], color='grey', linestyle='--')
ax[1, 1].plot(Tpar_x_right, Tpar_y_right, color='blue')
ax[1, 1].plot(Tpar_x_right_hires, Tpar_y_right_hires, color='red', linestyle='--')
# ax[1, 1].set_ylabel(r'$T_{||}$ (eV)')
ax[1, 1].yaxis.tick_right()
ax[1, 1].yaxis.set_label_position("right")
ax[1, 1].set_yscale('log')
ax[1, 1].set_ylim(0, 14000)
ax[1, 1].set_xlim(0, x_max)

# Split Tperp data at zero
Tperp_x_left, Tperp_y_left, Tperp_x_right, Tperp_y_right = split_at_zero(nodes_Z, Tperp_ev)
Tperp_x_left_hires, Tperp_y_left_hires, Tperp_x_right_hires, Tperp_y_right_hires = split_at_zero(nodes_Z_hires, Tperp_ev_hires)

# Tperp - Left panel
ax[1, 2].plot([-0.98, -0.98], [0, 35000], color='grey', linestyle='--')
ax[1, 2].plot(Tperp_x_left, Tperp_y_left, color='blue')
ax[1, 2].text(-2.3, 22000, r'$\mathbf{(d)}$')
ax[1, 2].plot(Tperp_x_left_hires, Tperp_y_left_hires, color='red', linestyle='--')
ax[1, 2].set_ylabel(r'$T_{\perp}$ (eV)')
ax[1, 2].set_ylim(0, 24000)
ax[1, 2].set_xlim(x_min, 0)

# Tperp - Right panel
ax[1, 3].plot([0.98, 0.98], [0, 3500000], color='grey', linestyle='--')
ax[1, 3].plot(Tperp_x_right, Tperp_y_right, color='blue')
ax[1, 3].plot(Tperp_x_right_hires, Tperp_y_right_hires, color='red', linestyle='--')
# ax[1, 3].set_ylabel(r'$T_{\perp}$ (eV)')
ax[1, 3].yaxis.tick_right()
ax[1, 3].yaxis.set_label_position("right")
ax[1, 3].set_yscale('log')
ax[1, 3].set_ylim(0, 100000)
ax[1, 3].set_xlim(0, x_max)

# Add a grid to each subplot
# for a in ax.flat:
#     a.grid(True, alpha=0.3)

# Set the x-axis label for the bottom row at the center of each pair of subplots
# For the left pair (columns 0 and 1), place xlabel on column 1
# Shift xlabel position to the left for both bottom subplots
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

# plt.tight_layout()
plt.savefig('compare-Nz800-2sOAP-time-dilation.pdf')
plt.show()
plt.close()  # Close figure to free memory