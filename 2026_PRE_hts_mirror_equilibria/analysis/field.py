import postgkyl as pg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import os

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

# ── Simulation settings ────────────────────────────────────────────────────────
read_frame = 65
directory = '/home/mr1884/scratch/gkylmax/stellar-lorentzian1x-orbit-average'
# directory = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-exploration/stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil20'
sim_name = 'gk_lorentzian_mirror'
species = 'ion'
z_mirror_throat = 0.98
simulation_edge = 2.5

# Ion-ion collision frequency from the simulation profile at the midplane
nusum_data = pg.GData(
    directory + f'/{sim_name}-{species}_lbo_nu_sum_{read_frame}.gkyl',
    mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl',
)
pg.GInterpModal(nusum_data).interpolate(comp=0, overwrite=True)
_, nu_sum_sel = pg.data.select(nusum_data, z0=0.0)
nu_ii = float(np.squeeze(nu_sum_sel))

# Universal constants (matching input_file.c)
eV  = 1.602176634e-19   # J  (elementary charge)
mu0 = 4 * np.pi * 1e-7  # H/m (permeability of free space)
eps0 = 8.8541878128e-12  # F/m (permittivity of free space)
mp  = 1.67262192369e-27  # kg (proton mass)
me  = 9.1093837015e-31   # kg (electron mass)

# Plasma parameters (from input_file.c)
mi  = 2.014 * mp   # kg (deuteron mass)
qe  = -eV          # C  (electron charge)
Te0 = 940 * eV     # J  (electron temperature)
n0  = 3e+19        # m^-3 (density at midplane)
B_p = 0.53         # T  (magnetic field strength)
beta = 0.4         # plasma beta
tau  = (B_p**2) * beta / (2.0 * mu0 * n0 * Te0) - 1.0  # temperature ratio
Ti0  = tau * Te0   # J  (ion temperature)

# Derived parameters
vti = np.sqrt(Ti0 / mi)   # m/s (ion thermal velocity)
vte = np.sqrt(Te0 / me)   # m/s (electron thermal velocity)
c_s = np.sqrt(Te0 / mi)   # m/s (sound speed)
omega_ci = eV * B_p / mi  # rad/s (ion cyclotron frequency)
rho_s = c_s / omega_ci    # m  (sound gyroradius)

# ── Helper ─────────────────────────────────────────────────────────────────────
def draw_bracket(ax, x_bracket, t_span, y_lo, y_hi, box_width, delta_sim, delta_theory,
                 x_offset=0):
    """Draw a labelled comparison bracket on *ax*.

    Parameters
    ----------
    x_tip       : x-coordinate of the bracket tip (data units)
    t_span      : total time span (t[-1] - t[0]) in nu_ii units
    y_lo, y_hi  : lower / upper y-values spanned by the bracket
    box_width   : width of the box containing the labels (data units)
    delta_sim   : simulated Δeφ/Te value
    delta_theory: theoretical Δeφ/Te value
    x_offset    : shift the box, arrow, and labels rightward by this amount
    """
    y_mid = (y_lo + y_hi) / 2
    y_range = y_hi - y_lo
    x_text = x_bracket - x_offset
    x_arrow = x_text + box_width * 0.55
    box_x = x_text - 0.5 * box_width

    ax.annotate(
        '',
        xy=(x_bracket*0.98, y_mid),
        xytext=(x_arrow, y_mid),
        arrowprops=dict(
            arrowstyle=f'<-[, lengthB=0.2, widthB={y_range*0.5}',
            lw=1, color='black'
        ),
        annotation_clip=False,
    )

    box = FancyBboxPatch(
        (box_x, y_mid - 1.1), box_width, 2.2,
        boxstyle='round,pad=0.2',
        linewidth=1, edgecolor='grey', facecolor='lightgrey',
        alpha=0.5, zorder=2,
    )
    ax.add_patch(box)

    ax.text(x_text, y_mid + 0.5,
            r'$\Delta e \phi / T_e = {:.2f}$'.format(delta_sim),
            ha='center', va='center', zorder=3)
    ax.text(x_text, y_mid - 0.7,
            r'$\Delta e\phi_{{\mathrm{{theory}}}} / T_e = {:.2f}$'.format(delta_theory),
            ha='center', va='center', zorder=3)


# ── Load field data ────────────────────────────────────────────────────────────
t_list, phi_0_list, phi_098_list, phi_25_list = [], [], [], []

for iframe in range(read_frame + 1):
    fpath = directory + f'/gk_lorentzian_mirror-field_{iframe}.gkyl'
    data_phi = pg.GData(fpath, mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
    # CLI equivalent: ... interp sel --z0 <value>
    pg.GInterpModal(data_phi).interpolate(comp=0, overwrite=True)
    _, phi_sel = pg.data.select(data_phi, z0='0.0')
    phi_0_list.append(float(np.squeeze(phi_sel)))
    _, phi_sel = pg.data.select(data_phi, z0=f'{z_mirror_throat}')
    phi_098_list.append(float(np.squeeze(phi_sel)))
    _, phi_sel = pg.data.select(data_phi, z0=f'{simulation_edge}')
    phi_25_list.append(float(np.squeeze(phi_sel)))
    t_list.append(data_phi.ctx['time'])

# Prepend a duplicate of the first point so the trace starts as a flat line
t_list    = [t_list[0],    t_list[1]]    + t_list[1:]
phi_0_list   = [phi_0_list[0]]   * 2 + phi_0_list[1:]
phi_098_list = [phi_098_list[0]] * 2 + phi_098_list[1:]
phi_25_list  = [phi_25_list[0]]  * 2 + phi_25_list[1:]

t      = np.array(t_list)
phi_0   = np.array(phi_0_list)   * -qe / Te0
phi_098 = np.array(phi_098_list) * -qe / Te0
phi_25  = np.array(phi_25_list)  * -qe / Te0

# ── Plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(1, 1, figsize=(5, 4))

cmap   = plt.cm.inferno
colors = [cmap(i / 3) for i in range(4)]

ax.plot(t * nu_ii, phi_0,   label=r'$z=0.0$ m',  color=colors[0])
ax.plot(t * nu_ii, phi_098, label=r'$z=0.98$ m', color=colors[1])
ax.plot(t * nu_ii, phi_25,  label=r'$z=2.5$ m',  color=colors[2])

t_span = (t[-1] - t[0]) * nu_ii
x_tip  = t[-1] * nu_ii

print(f'Delta phi (midplane → mirror throat): {phi_0[-1] - phi_098[-1]:.4f}')
print(f'percent error in delta phi (midplane → mirror throat): {100 * (phi_0[-1] - phi_098[-1] - 6.93) / 6.93:.2f}%')
draw_bracket(ax, x_tip, t_span, phi_098[-1], phi_0[-1], box_width=2.6,
             delta_sim=phi_0[-1] - phi_098[-1],   delta_theory=6.93, x_offset=2.2)
draw_bracket(ax, x_tip, t_span, phi_25[-1],  phi_098[-1], box_width=2.6,
             delta_sim=phi_098[-1] - phi_25[-1],  delta_theory=6.49, x_offset=2.2)

ax.set_xlabel(r'$t \nu_{ii}$')
ax.set_ylabel(r'$e\phi / T_e$')
ax.set_xlim(0, t[-1] * nu_ii)
ax.legend(loc='upper left', bbox_to_anchor=(0.05, 0.80), framealpha=1.0)

plt.tight_layout()
plt.savefig('field_fig_lorentzian.pdf')
plt.show()