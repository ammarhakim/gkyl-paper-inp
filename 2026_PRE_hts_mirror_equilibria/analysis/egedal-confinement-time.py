import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 12,
        "axes.titlesize": 18,
        "axes.labelsize": 22,
        "legend.fontsize": 12,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    }
)

### Beam sources
beams_intM0dx_R3 = 1.653159e+19
beams_intM0dx_R5 = 2.044288e+19
beams_intM0dx_R10 = 2.579418e+19
beams_intM0dx_R15 = 2.860074e+19
beams_intM0dx_R22 = 3.147811e+19
beams_intM0dx_R32 = 3.419167e+19
beams_intM0dx_R50 = 3.754514e+19

beams_intM0Sdx_R3 = 3.171798e+20
beams_intM0Sdx_R5 = 3.172110e+20
beams_intM0Sdx_R10 = 3.172137e+20
beams_intM0Sdx_R15 = 3.172138e+20
beams_intM0Sdx_R22 = 3.172137e+20
beams_intM0Sdx_R32 = 3.172138e+20
beams_intM0Sdx_R50 = 3.166954e+20

beams_ion_M0mid_R3 = 1.277128e+19
beams_ion_M0mid_R5 = 1.503446e+19
beams_ion_M0mid_R10 = 1.921113e+19
beams_ion_M0mid_R15 = 2.107043e+19
beams_ion_M0mid_R22 = 2.282627e+19
beams_ion_M0mid_R32 = 2.440953e+19
beams_ion_M0mid_R50 = 2.612092e+19

beams_tau_pi_R3 = beams_intM0dx_R3 / beams_intM0Sdx_R3
beams_tau_pi_R5 = beams_intM0dx_R5 / beams_intM0Sdx_R5
beams_tau_pi_R10 = beams_intM0dx_R10 / beams_intM0Sdx_R10
beams_tau_pi_R15 = beams_intM0dx_R15 / beams_intM0Sdx_R15
beams_tau_pi_R22 = beams_intM0dx_R22 / beams_intM0Sdx_R22
beams_tau_pi_R32 = beams_intM0dx_R32 / beams_intM0Sdx_R32
beams_tau_pi_R50 = beams_intM0dx_R50 / beams_intM0Sdx_R50

beams_dens_vec = np.array([beams_ion_M0mid_R3, beams_ion_M0mid_R5, beams_ion_M0mid_R10, \
                           beams_ion_M0mid_R15, beams_ion_M0mid_R22, beams_ion_M0mid_R32, beams_ion_M0mid_R50])
beams_dens_vec /= 1e20  # Convert to 10^19 m^-3 for better scaling in fits and plots.

R_values = np.array([3, 5, 10, 15, 22, 32, 50])


def egedal_fit_func(R, k, c): # 25 kev beam
    return k * 25 ** (3 / 2) * np.log10(R) + c

# --- Beam fit and analysis ---
beams_tau_values = np.array(
    [
        beams_tau_pi_R3,
        beams_tau_pi_R5,
        beams_tau_pi_R10,
        beams_tau_pi_R15,
        beams_tau_pi_R22,
        beams_tau_pi_R32,
        beams_tau_pi_R50,
    ]
)

# Fit to Egedal scaling
# Recast as y = n_i * tau_pi and x = 25^(3/2) * log10(R).
x_values = 25 ** (3 / 2) * np.log10(R_values)
time_unit = 1e-4  # seconds
y_values = beams_dens_vec * beams_tau_values / time_unit

# Least-squares fit and covariance estimate for [k, c].
(k_fit_egedal, c_fit_egedal), pcov_egedal = np.polyfit(x_values, y_values, deg=1, cov=True)
k_err_egedal, c_err_egedal = np.sqrt(np.diag(pcov_egedal))
y_fit = egedal_fit_func(R_values, k_fit_egedal, c_fit_egedal)

print(f"[Egedal] Fitted parameter k (10^-4 s units): {k_fit_egedal}")
print(f"[Egedal] Fitted parameter C (10^-4 s units): {c_fit_egedal}")
print(f"[Egedal] k standard error (10^-4 s units): {k_err_egedal}")
print(f"[Egedal] C standard error (10^-4 s units): {c_err_egedal}")
print(
    f"[Egedal] R^2: {1 - np.sum((y_values - y_fit) ** 2) / np.sum((y_values - np.mean(y_values)) ** 2)}"
)

plt.figure(figsize=(5, 4))
plt.plot(
    R_values,
    y_values,
    "o",
    label="Beam data",
)
plt.plot(
    R_values,
    y_fit,
    label="Egedal fit",
)
plt.xlabel("R (m)")
plt.ylabel(r"$n_i\,\tau_i\;[10^{-4}\,\mathrm{s}]$")
plt.title(r"$n_i\,\tau_i$ vs R with Egedal Fit ($10^{-4}$ s units)")
plt.xscale("log")
plt.legend()
plt.tight_layout()
plt.show()