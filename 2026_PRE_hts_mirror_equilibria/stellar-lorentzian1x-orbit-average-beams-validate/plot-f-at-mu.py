import postgkyl as pg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

# ============================================================================
# INPUT PARAMETERS - Modify these as needed
# ============================================================================
# Directory containing simulation data
directory = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-exploration/stellar-lorentzian1x-orbit-average-beams'

# Frame number to read
read_frame = 0

# Distribution function file (relative to directory)
dist_func_file = f'gk_lorentzian_mirror-ion_source_{read_frame}.gkyl'

# Velocity space mapping file (for non-uniform grids)
mapc2p_vel_file = 'gk_lorentzian_mirror-ion_mapc2p_vel.gkyl'

# Jacobian files
jacobvel_file = 'gk_lorentzian_mirror-ion_jacobvel.gkyl'
jacobtot_file = 'gk_lorentzian_mirror-jacobtot.gkyl'

# Position mapping file
mapc2p_file = 'gk_lorentzian_mirror-mapc2p.gkyl'

# Mu indices to plot (indices into the mu grid)
# These are the array indices for mu at which to slice the distribution
# The simulation has 64 mu cells (indices 0-63), with mu=0 at index 0
plot_mu_indices = [0, 5, 10, 20, 30, 40]

# Colorbar limits (set to None for auto-scaling)
colorbar_vmin = 1e-16  # Minimum value for colorbar
colorbar_vmax = None   # Maximum value for colorbar (None = auto)

# Output filename
output_filename = 'f_at_mu_slices.pdf'
# ============================================================================

# Universal constants
eV = 1.602176634e-19  # J (elementary charge / electron volt)
mp = 1.67262192369e-27  # kg (proton mass)

# Plasma parameters (from input_file.c)
mi = 2.014 * mp  # kg (deuteron mass)
Te0 = 940 * eV  # J (electron temperature)
n0 = 3e+19
B_p = 0.53  # T (magnetic field strength)
mu0 = 4 * np.pi * 1e-7  # H/m
beta = 0.4
tau = (B_p**2) * beta / (2.0 * mu0 * n0 * Te0) - 1.0
Ti0 = tau * Te0  # J (ion temperature)

# Derived parameters
vti = np.sqrt(Ti0 / mi)  # m/s (ion thermal velocity)


def plot_f_at_mu():
    """
    Plot the distribution function at specific mu coordinates,
    showing f(z, vpar) slices.
    """
    print(f"Loading distribution function from: {directory}/{dist_func_file}")
    
    # Load distribution function data
    f_data = pg.GData(
        f'{directory}/{dist_func_file}',
        mapc2p_vel_name=f'{directory}/{mapc2p_vel_file}'
    )
    
    # Load Jacobians
    Jv_data = pg.GData(f'{directory}/{jacobvel_file}')
    jb_data = pg.GData(f'{directory}/{jacobtot_file}')
    
    # Get values and divide by velocity Jacobian
    f_c = f_data.get_values()
    Jv_c = Jv_data.get_values()
    f_data._values = f_c / Jv_c
    
    # Interpolate
    dg = pg.GInterpModal(f_data, 1, 'gkhyb')
    jb_dg = pg.GInterpModal(jb_data, 1, 'ms')
    
    xInt_i, fIon = dg.interpolate()
    _, jb = jb_dg.interpolate(0)
    
    fIon = np.squeeze(fIon)
    jb = np.squeeze(jb)
    
    # Load position mapping to get physical Z coordinates
    data_mc2p = pg.GData(f'{directory}/{mapc2p_file}')
    interp_mc2p = pg.GInterpModal(data_mc2p, 1, 'ms')
    nodes_Z = interp_mc2p.interpolate(2)[1]
    nodes_Z = np.squeeze(nodes_Z)
    
    # Set small/negative values to minimum for log plot
    fIon[fIon < 1e-40] = 1e-40
    
    # Get grid dimensions
    ndim = len(xInt_i)
    nxInt_i = [np.size(xInt_i[d]) for d in range(ndim)]
    
    # Get cell center coordinates
    xIntC_i = [0.5 * (xInt_i[d][:-1] + xInt_i[d][1:]) for d in range(ndim)]
    
    # Print grid info
    print(f"\nGrid dimensions:")
    print(f"  z (dim 0): {nxInt_i[0]} nodes, range [{xInt_i[0].min():.3f}, {xInt_i[0].max():.3f}]")
    print(f"  vpar (dim 1): {nxInt_i[1]} nodes, range [{xInt_i[1].min():.3e}, {xInt_i[1].max():.3e}] m/s")
    print(f"  vpar/vti range: [{xInt_i[1].min()/vti:.2f}, {xInt_i[1].max()/vti:.2f}]")
    print(f"  mu (dim 2): {nxInt_i[2]} nodes, range [{xInt_i[2].min():.3e}, {xInt_i[2].max():.3e}]")
    
    # Use physical Z coordinates for plotting
    # nodes_Z from mapc2p gives the physical Z coordinates corresponding to computational z
    # xInt_i[0] has the computational z nodal coordinates (201 points for 200 cells)
    # We need to map computational z to physical Z for plotting
    
    # The interpolated grid xInt_i[0] has nodal coordinates
    # fIon has cell-centered values, so fIon.shape[0] = len(xInt_i[0]) - 1
    z_nodal = xInt_i[0]  # Use computational z coordinates (nodal)
    
    # Normalize velocity coordinates for labeling
    vpar_normalized = xInt_i[1] / vti
    mu_normalized = xInt_i[2] / (0.5 * mi * (vti**2) / B_p)
    
    # Get mu values at the specified indices (using cell centers)
    mu_center_normalized = xIntC_i[2] / (0.5 * mi * (vti**2) / B_p)
    
    # Validate mu indices
    n_mu_cells = len(xIntC_i[2])
    actual_mu_indices = []
    actual_mu_values = []
    for mu_idx in plot_mu_indices:
        if mu_idx < 0 or mu_idx >= n_mu_cells:
            print(f"  WARNING: mu index {mu_idx} out of range [0, {n_mu_cells-1}], skipping")
            continue
        actual_mu_indices.append(mu_idx)
        actual_mu_values.append(mu_center_normalized[mu_idx])
        print(f"  Using mu index {mu_idx}, mu/mu_ti = {actual_mu_values[-1]:.3f}")
    
    if len(actual_mu_indices) == 0:
        print("ERROR: No valid mu indices specified!")
        return
    
    # Determine subplot layout
    n_plots = len(actual_mu_indices)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    
    # Handle case of single plot
    if n_plots == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Create meshgrid for plotting (z, vpar)
    # xInt_i contains nodal coordinates, fIon has cell values
    # For pcolormesh with shading='flat', we need nodal coords with shape (N+1, M+1) for data (N, M)
    # z_nodal has 201 points, vpar_normalized has 193 points
    # f_slice will have shape (200, 192) after slicing at a mu index
    Z_nodal, VPAR_nodal = np.meshgrid(z_nodal, vpar_normalized, indexing='ij')
    
    print(f"\nMeshgrid shapes: Z_nodal={Z_nodal.shape}, VPAR_nodal={VPAR_nodal.shape}")
    print(f"fIon shape: {fIon.shape}")
    
    # Find global colorbar limits if not specified
    if colorbar_vmax is None:
        vmax = 0
        for idx in actual_mu_indices:
            # Divide by jacobtot (which only depends on z)
            f_slice = fIon[:, :, idx] / jb[:, np.newaxis]
            vmax = max(vmax, np.max(f_slice))
    else:
        vmax = colorbar_vmax
    
    vmin = colorbar_vmin if colorbar_vmin is not None else 1e-16
    
    print(f"\nColorbar range: [{vmin:.2e}, {vmax:.2e}]")
    
    # Plot each mu slice
    pcm = None
    for i, (mu_idx, actual_mu) in enumerate(zip(actual_mu_indices, actual_mu_values)):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # Extract f(z, vpar) at this mu index
        # fIon shape is (nz, nvpar, nmu), we want fIon[:, :, mu_idx]
        f_slice = fIon[:, :, mu_idx] / jb[:, np.newaxis]
        
        # Set small values for log scale
        f_slice[f_slice < vmin] = vmin
        
        # pcolormesh with shading='flat' (default): X,Y should have shape (N+1, M+1) for C shape (N, M)
        pcm = ax.pcolormesh(
            Z_nodal, VPAR_nodal, f_slice,
            cmap='inferno',
            norm=LogNorm(vmin=vmin, vmax=vmax),
            rasterized=True
        )
        
        ax.set_title(rf'$\mu/\mu_{{ti}} = {actual_mu:.2f}$ (idx={mu_idx})')
        
        if col == 0:
            ax.set_ylabel(r'$v_{||} / v_{ti}$')
        if row == n_rows - 1:
            ax.set_xlabel(r'$z$ (m)')
        
        # # Add panel label
        # label = chr(ord('a') + i)
        # ax.text(0.05, 0.95, rf'$\mathbf{{({label})}}$', 
        #         transform=ax.transAxes, color='white', 
        #         fontsize=14, va='top', ha='left')
    
    # Remove empty subplots
    for i in range(n_plots, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    # Add colorbar
    if pcm is not None:
        cbar = fig.colorbar(pcm, ax=axes.ravel().tolist(), orientation='vertical', 
                           fraction=0.02, pad=0.02)
        cbar.set_label(r'$f_i / (J_v \cdot J_{\mathrm{tot}})$', fontsize=18)
    
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{directory}/{output_filename}', dpi=150, bbox_inches='tight')
    print(f"\nSaved plot to: {directory}/{output_filename}")
    plt.close()


if __name__ == "__main__":
    plot_f_at_mu()
