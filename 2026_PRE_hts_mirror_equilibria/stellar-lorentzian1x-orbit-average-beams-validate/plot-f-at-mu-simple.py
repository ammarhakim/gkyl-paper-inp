"""
Simple script to plot the distribution function f(z, vpar) at a single mu index.
"""
import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ============================================================================
# INPUT PARAMETERS
# ============================================================================
directory = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-exploration/stellar-lorentzian1x-orbit-average-beams'
dist_func_file = 'gk_lorentzian_mirror-ion_source_0.gkyl'
mapc2p_vel_file = 'gk_lorentzian_mirror-ion_mapc2p_vel.gkyl'
jacobvel_file = 'gk_lorentzian_mirror-ion_jacobvel.gkyl'
jacobtot_file = 'gk_lorentzian_mirror-jacobtot.gkyl'

mu_index = int(0.7*64)  # Index into mu grid to plot

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================
eV = 1.602176634e-19
mp = 1.67262192369e-27
mi = 2.014 * mp
Te0 = 940 * eV
B_p = 0.53
mu0 = 4 * np.pi * 1e-7
n0 = 3e19
beta = 0.4
tau = (B_p**2) * beta / (2.0 * mu0 * n0 * Te0) - 1.0
Ti0 = tau * Te0
vti = np.sqrt(Ti0 / mi)
mu_ti = 0.5 * mi * vti**2 / B_p  # Thermal mu for normalization

# ============================================================================
# LOAD DATA
# ============================================================================
# Load distribution function with velocity mapping
f_data = pg.GData(f'{directory}/{dist_func_file}',
                  mapc2p_vel_name=f'{directory}/{mapc2p_vel_file}')

# Load Jacobians
Jv_data = pg.GData(f'{directory}/{jacobvel_file}')
jb_data = pg.GData(f'{directory}/{jacobtot_file}')

# Divide by velocity Jacobian
f_data._values = f_data.get_values() / Jv_data.get_values()

# Interpolate to grid
dg = pg.GInterpModal(f_data, 1, 'gkhyb')
jb_dg = pg.GInterpModal(jb_data, 1, 'ms')

coords, f = dg.interpolate()
_, jb = jb_dg.interpolate(0)

f = np.squeeze(f)
jb = np.squeeze(jb)

# coords[0] = z, coords[1] = vpar, coords[2] = mu (all nodal)
z = coords[0]
vpar = coords[1] / vti  # Normalize by thermal velocity
mu = coords[2] / mu_ti  # Normalize by thermal mu

# Get cell-center mu values for labeling
mu_centers = 0.5 * (mu[:-1] + mu[1:])

print(f"Grid: z={len(z)-1} cells, vpar={len(vpar)-1} cells, mu={len(mu)-1} cells")
print(f"Plotting at mu index {mu_index}, mu/mu_ti = {mu_centers[mu_index]:.3f}")

# ============================================================================
# EXTRACT SLICE AND PLOT
# ============================================================================
# f has shape (nz, nvpar, nmu) - extract f(z, vpar) at fixed mu
f_slice = f[:, :, mu_index] / jb[:, np.newaxis]
# f_slice[f_slice < 1e-40] = 1e-40  # Floor for log scale

# Create meshgrid for pcolormesh (nodal coordinates)
Z, VPAR = np.meshgrid(z, vpar, indexing='ij')

# Plot
fig, ax = plt.subplots(figsize=(8, 6))

pcm = ax.pcolormesh(Z, VPAR, f_slice, cmap='inferno')
                    # norm=LogNorm(vmin=1e-16, vmax=f_slice.max()),
                    # shading='flat')

ax.set_xlabel('z (m)')
ax.set_ylabel(r'$v_{||} / v_{ti}$')
ax.set_title(rf'$f(z, v_{{||}})$ at $\mu/\mu_{{ti}} = {mu_centers[mu_index]:.2f}$')

cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label(r'$f / (J_v \cdot J_{tot})$')

plt.tight_layout()
plt.show()
