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

read_frame_maxwellian = 65
read_frame_beam = 65

directory_POA_simulation = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average'
directory_POA_simulation_beam = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average-beams'

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

f_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_'+str(read_frame_maxwellian)+'.gkyl',\
                    mapc2p_vel_name=directory_POA_simulation + '/gk_lorentzian_mirror-ion_mapc2p_vel.gkyl')
Jv_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-ion_jacobvel.gkyl')
jb_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-jacobtot.gkyl')
bmag_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-bmag.gkyl')
phi_data = pg.GData(directory_POA_simulation + '/gk_lorentzian_mirror-field_' + str(read_frame_maxwellian) + '.gkyl')
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
_, nodes_Z = interp_mc2p.interpolate(1)
nodes_Z = np.squeeze(nodes_Z)

# Set fIon = 1e-60 where it is <0
fIon = np.abs(fIon)
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
plot_z_locations = [0., 0.98, 2.5]

plotzIdx = [np.argmin(np.abs(nodes_Z-val)) for val in plot_z_locations]

Bmin = bmag[plotzIdx[0]]
Bmax = bmag[plotzIdx[1]]

phi_center = phi[plotzIdx[0]]
phi_mirror = phi[plotzIdx[1]]

mu_loss = (( 1/2 * mi * xInt_i[1]**2 ) + qi * (phi_center - phi_mirror)) / (Bmax - Bmin)
mu_loss = mu_loss  /(0.5 * mi*(vti**2)/B_p)

#[ Normalize velocity space
xInt_i[1] = xInt_i[1]/vti
xInt_i[2] = xInt_i[2]/(0.5*mi*(vti**2)/B_p)
xIntC_i[1] = xIntC_i[1]/vti
xIntC_i[2] = xIntC_i[2]/(0.5*mi*(vti**2)/B_p)

#[ Create colorplot grid. Recall coordinates have to be nodal.
Xnodal_i = [np.outer(xInt_i[1], np.ones(nxInt_i[2])),
            np.outer(np.ones(nxInt_i[1]), xInt_i[2])]

extreme_vals = [0., np.amax(fIon[plotzIdx[0],:,:]/jb[plotzIdx[0]])]

#[ ---- Beams simulation data ----
f_data_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-ion_'+str(read_frame_beam)+'.gkyl',
                    mapc2p_vel_name=directory_POA_simulation_beam + '/gk_lorentzian_mirror-ion_mapc2p_vel.gkyl')
Jv_data_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-ion_jacobvel.gkyl')
jb_data_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-jacobtot.gkyl')
bmag_data_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-bmag.gkyl')
phi_data_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-field_'+str(read_frame_beam)+'.gkyl')
f_c_b = f_data_b.get_values()
Jv_c_b = Jv_data_b.get_values()
f_data_b._values = f_c_b/Jv_c_b
dg_b = pg.GInterpModal(f_data_b, 1, 'gkhyb')
bmag_dg_b = pg.GInterpModal(bmag_data_b, 1, 'ms')
jb_dg_b = pg.GInterpModal(jb_data_b, 1, 'ms')
phi_dg_b = pg.GInterpModal(phi_data_b, 1, 'ms')
xInt_i_b, fIon_b = dg_b.interpolate()
_, bmag_b = bmag_dg_b.interpolate(0)
_, phi_b = phi_dg_b.interpolate(0)
_, jb_b = jb_dg_b.interpolate(0)
fIon_b = np.squeeze(fIon_b)
bmag_b = np.squeeze(bmag_b)
phi_b = np.squeeze(phi_b)
jb_b = np.squeeze(jb_b)

data_mc2p_b = pg.GData(directory_POA_simulation_beam + '/gk_lorentzian_mirror-mapc2p.gkyl')
interp_mc2p_b = pg.GInterpModal(data_mc2p_b, 1, 'ms')
_, nodes_Z_b = interp_mc2p_b.interpolate(1)
nodes_Z_b = np.squeeze(nodes_Z_b)

fIon_b = np.abs(fIon_b)
fIon_b[fIon_b < 1e-16] = 1e-16

ndim_b = len(xInt_i_b)
nxInt_i_b = [np.size(xInt_i_b[d]) for d in range(ndim_b)]

xInt_i_b[0] = nodes_Z_b

plotzIdx_b = [np.argmin(np.abs(nodes_Z_b-val)) for val in plot_z_locations]

Bmin_b = bmag_b[plotzIdx_b[0]]
Bmax_b = bmag_b[plotzIdx_b[1]]

phi_center_b = phi_b[plotzIdx_b[0]]
phi_mirror_b = phi_b[plotzIdx_b[1]]

mu_loss_b = ((1/2 * mi * xInt_i_b[1]**2) + qi * (phi_center_b - phi_mirror_b)) / (Bmax_b - Bmin_b)
mu_loss_b = mu_loss_b / (0.5 * mi*(vti**2)/B_p)

#[ Normalize beam velocity space
xInt_i_b[1] = xInt_i_b[1]/vti
xInt_i_b[2] = xInt_i_b[2]/(0.5*mi*(vti**2)/B_p)

#[ Create colorplot grid for beam simulation
Xnodal_i_b = [np.outer(xInt_i_b[1], np.ones(nxInt_i_b[2])),
              np.outer(np.ones(nxInt_i_b[1]), xInt_i_b[2])]

extreme_vals_b = [0., np.amax(fIon_b[plotzIdx_b[0],:,:]/jb_b[plotzIdx_b[0]])]

fig, ax = plt.subplots(2, 3, figsize=(10, 6), gridspec_kw={'width_ratios': [1.2, 1, 1.2]})

colorbar_zmin = 1e-9
xlim = (-8, 8)
inset_xlim = (0, 5)
inset_ylim = (0, 1)

pcm = ax[0,0].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[0],:,:]/jb[plotzIdx[0]],
                        cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
ax[0,0].plot(xInt_i[1], mu_loss, color='white', linestyle='--', label='Loss Cone')
# ax[0,0].plot(xInt_i[1], mu_loss*2, color='grey', linestyle='--', label='Loss Cone')
ax[0,0].set_ylim(0, 8)
ax[0,0].set_ylabel(r'Maxwellian' + '\n' + r'$\mu / \mu_{ti,0}$')
ax[0,0].set_title(r'z = {:.1f} m'.format(np.abs(nodes_Z[plotzIdx[0]])))
ax[0,0].set_xlim(xlim)

# Add a linear colorbar for the leftmost column
from matplotlib.colors import Normalize
linear_norm = Normalize(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1])
pcm_linear = ax[0,0].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[0],:,:]/jb[plotzIdx[0]],
                        cmap='inferno', norm=linear_norm, rasterized=True)
cb_linear = plt.colorbar(pcm_linear, ax=ax[0,0], orientation='vertical', pad=0.02)
# cb_linear.set_label(r'$ |f_i | $ (linear)', fontsize=18)
cb_linear.ax.tick_params(labelsize=10)

pcm = ax[0,1].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[1],:,:]/jb[plotzIdx[1]],
                        cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
ax[0,1].set_title(r'z = {:.2f} m'.format(nodes_Z[plotzIdx[1]]))
ax[0,1].set_xlim(xlim)
ax[0,1].set_ylim(0, 8)

axins_01 = inset_axes(ax[0,1], width="60%", height="42%", loc='upper center', borderpad=1.0)
axins_01.pcolormesh(
  Xnodal_i[0],
  Xnodal_i[1],
  fIon[plotzIdx[1],:,:]/jb[plotzIdx[1]],
  cmap='inferno',
  norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]),
  rasterized=True,
)
axins_01.set_xlim(0, 2)
axins_01.set_ylim(0, 0.5)
axins_01.set_xticks([0, 2])
axins_01.set_yticks([0, 0.5])
axins_01.tick_params(labelsize=8)
for spine in axins_01.spines.values():
  spine.set_edgecolor('0.6')
  spine.set_linewidth(1.1)
mark_inset(ax[0,1], axins_01, loc1=3, loc2=4, fc='none', ec='white', lw=0.8)

pcm = ax[0,2].pcolormesh(Xnodal_i[0], Xnodal_i[1], fIon[plotzIdx[2],:,:]/jb[plotzIdx[2]],
                        cmap='inferno', norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]), rasterized=True)
cbar = plt.colorbar(pcm, ax=ax[0,2], orientation='vertical')
cbar.set_label(r'$ f_i $'.format(nodes_Z[plotzIdx[2]]), fontsize=22)
cbar.ax.tick_params(labelsize=10)
ax[0,2].set_title(r'z = {:.1f} m'.format(nodes_Z[plotzIdx[2]]))
ax[0,2].set_xlim(xlim)
ax[0,2].set_ylim(0, 8)

axins_02 = inset_axes(ax[0,2], width="60%", height="42%", loc='upper center', borderpad=1.0)
axins_02.pcolormesh(
  Xnodal_i[0],
  Xnodal_i[1],
  fIon[plotzIdx[2],:,:]/jb[plotzIdx[2]],
  cmap='inferno',
  norm=LogNorm(vmin=max(extreme_vals[0], colorbar_zmin), vmax=extreme_vals[1]),
  rasterized=True,
)
axins_02.set_xlim(inset_xlim)
axins_02.set_ylim(inset_ylim)
axins_02.set_xticks([0, 5])
axins_02.set_yticks([0, 1])
axins_02.tick_params(labelsize=8)
for spine in axins_02.spines.values():
  spine.set_edgecolor('0.6')
  spine.set_linewidth(1.1)
mark_inset(ax[0,2], axins_02, loc1=3, loc2=4, fc='none', ec='white', lw=0.8)



pcm_b = ax[1,0].pcolormesh(Xnodal_i_b[0], Xnodal_i_b[1], fIon_b[plotzIdx_b[0],:,:]/jb_b[plotzIdx_b[0]],
                        cmap='inferno', rasterized=True)
ax[1,0].plot(xInt_i_b[1], mu_loss_b, color='white', linestyle='--', label='Loss Cone')
ax[1,0].set_ylim(0, 8)
ax[1,0].set_ylabel(r'Beams' + '\n' + r'$\mu / \mu_{ti,0}$')
# ax[1,0].set_xlabel(r'$v_{||}/v_{ti,0}$')
ax[1,0].set_xlim(xlim)
cb_linear = plt.colorbar(pcm_b, ax=ax[1,0], orientation='vertical', pad=0.02)
cb_linear.ax.tick_params(labelsize=10)

pcm_b = ax[1,1].pcolormesh(Xnodal_i_b[0], Xnodal_i_b[1], fIon_b[plotzIdx_b[1],:,:]/jb_b[plotzIdx_b[1]],
                        cmap='inferno', norm=LogNorm(vmin=max(extreme_vals_b[0], colorbar_zmin), vmax=extreme_vals_b[1]), rasterized=True)
ax[1,1].set_xlabel(r'$v_{||}/v_{ti,0}$')
ax[1,1].set_xlim(xlim)

axins_11 = inset_axes(ax[1,1], width="60%", height="42%", loc='upper center', borderpad=1.0)
axins_11.pcolormesh(
  Xnodal_i_b[0],
  Xnodal_i_b[1],
  fIon_b[plotzIdx_b[1],:,:]/jb_b[plotzIdx_b[1]],
  cmap='inferno',
  norm=LogNorm(vmin=max(extreme_vals_b[0], colorbar_zmin), vmax=extreme_vals_b[1]),
  rasterized=True,
)
axins_11.set_xlim(0, 2)
axins_11.set_ylim(0, 0.5)
axins_11.set_xticks([0, 2])
axins_11.set_yticks([0, 0.5])
axins_11.tick_params(labelsize=8)
for spine in axins_11.spines.values():
  spine.set_edgecolor('0.6')
  spine.set_linewidth(1.1)
mark_inset(ax[1,1], axins_11, loc1=3, loc2=4, fc='none', ec='white', lw=0.8)

ax[1,1].set_ylim(0, 8)
ax[1,2].set_ylim(0, 8)


pcm_b = ax[1,2].pcolormesh(Xnodal_i_b[0], Xnodal_i_b[1], fIon_b[plotzIdx_b[2],:,:]/jb_b[plotzIdx_b[2]],
                        cmap='inferno', norm=LogNorm(vmin=max(extreme_vals_b[0], colorbar_zmin), vmax=extreme_vals_b[1]), rasterized=True)
cbar_b = plt.colorbar(pcm_b, ax=ax[1,2], orientation='vertical')
cbar_b.set_label(r'$ f_i $', fontsize=22)
cbar_b.ax.tick_params(labelsize=10)
# ax[1,2].set_xlabel(r'$v_{||}/v_{ti,0}$')
ax[1,2].set_xlim(xlim)
ax[1,2].set_ylim(0, 8)

axins_12 = inset_axes(ax[1,2], width="60%", height="42%", loc='upper center', borderpad=1.0)
axins_12.pcolormesh(
  Xnodal_i_b[0],
  Xnodal_i_b[1],
  fIon_b[plotzIdx_b[2],:,:]/jb_b[plotzIdx_b[2]],
  cmap='inferno',
  norm=LogNorm(vmin=max(extreme_vals_b[0], colorbar_zmin), vmax=extreme_vals_b[1]),
  rasterized=True,
)
axins_12.set_xlim(inset_xlim)
axins_12.set_ylim(inset_ylim)
axins_12.set_xticks([0, 5])
axins_12.set_yticks([0, 1])
axins_12.tick_params(labelsize=8)
for spine in axins_12.spines.values():
  spine.set_edgecolor('0.6')
  spine.set_linewidth(1.1)
mark_inset(ax[1,2], axins_12, loc1=3, loc2=4, fc='none', ec='white', lw=0.8)


ax[0,0].text(-7, 7 , r'$\mathbf{(a)}$',color='white', fontsize=18)
ax[0,1].text(-7, 7 , r'$\mathbf{(b)}$',color='white', fontsize=18)
ax[0,2].text(-7, 7 , r'$\mathbf{(c)}$',color='white', fontsize=18)
ax[1,0].text(-7, 7 , r'$\mathbf{(d)}$',color='white', fontsize=18)
ax[1,1].text(-7, 7 , r'$\mathbf{(e)}$',color='white', fontsize=18)
ax[1,2].text(-7, 7 , r'$\mathbf{(f)}$',color='white', fontsize=18)

# plt.show()
plt.tight_layout()
plt.savefig('f_distributions_lorentzian.pdf')
plt.show()  # Close figure to free memory