import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, curve_fit
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


def B_mag(z, Bhat, gamma):
    z_m = 0.98
    lor_1 = 1 / (1 + ((z - z_m) / (gamma)) ** 2)
    lor_2 = 1 / (1 + ((z + z_m) / (gamma)) ** 2)
    Bmag = Bhat / np.pi / gamma * (lor_1 + lor_2)
    return Bmag


Bhat = 6.5
gamma = 0.124

z = np.linspace(-2.0, 2.0, 1000)
B = B_mag(z, Bhat, gamma)

R = B_mag(0.98, Bhat, gamma) / B_mag(0.0, Bhat, gamma)
print("Initial R =", R)
print("Initial B(0) =", B_mag(0.0, Bhat, gamma))

# Optimize Bhat and gamma for a desired R, holding B_mag(0) fixed
B0_target = B_mag(0.0, Bhat, gamma)  # Fix B(0) at current value
desired_R_values = np.array([3, 5, 10, 15, 22, 31.86])  # Array of desired mirror ratios


def constraints(params, desired_R):
    """
    System of equations to solve:
    1. B_mag(0, Bhat, gamma) = B0_target
    2. B_mag(0.98, Bhat, gamma) / B_mag(0, Bhat, gamma) = desired_R
    """
    Bhat_new, gamma_new = params

    # Constraint 1: B(0) should equal B0_target
    B0_new = B_mag(0.0, Bhat_new, gamma_new)
    eq1 = B0_new - B0_target

    # Constraint 2: Mirror ratio should equal desired_R
    B_mirror = B_mag(0.98, Bhat_new, gamma_new)
    R_new = B_mirror / B0_new
    eq2 = R_new - desired_R

    return [eq1, eq2]


# Store results
results = []

print("\n--- Optimized Parameters for Different Mirror Ratios ---")
print(f"{'R':>6s}  {'Bhat':>10s}  {'gamma':>10s}  {'B(0)':>10s}  {'R_actual':>10s}")
print("-" * 60)

for desired_R in desired_R_values:
    # Initial guess for optimization
    initial_guess = [Bhat, gamma]

    # Solve the system of equations
    solution = fsolve(constraints, initial_guess, args=(desired_R,))
    Bhat_opt, gamma_opt = solution

    # Verify results
    B0_result = B_mag(0.0, Bhat_opt, gamma_opt)
    R_result = B_mag(0.98, Bhat_opt, gamma_opt) / B0_result

    results.append(
        {
            "R_target": desired_R,
            "Bhat": Bhat_opt,
            "gamma": gamma_opt,
            "B0": B0_result,
            "R_actual": R_result,
        }
    )

    print(
        f"{desired_R:6.1f}  {Bhat_opt:10.6f}  {gamma_opt:10.6f}  {B0_result:10.6f}  {R_result:10.6f}"
    )



### Maxwellian source
intM0dx_R3 = 1.121816e+19
intM0dx_R5 = 1.625025e+19
intM0dx_R10 = 2.158982e+19
intM0dx_R15 = 2.411287e+19
intM0dx_R22 = 2.671260e+19
intM0dx_R32 = 2.898487e+19
intM0dx_R50 = 3.145267e+19

intM0Sdx_R3 = 3.172150e+20
intM0Sdx_R5 = 3.172143e+20
intM0Sdx_R10 = 3.172140e+20
intM0Sdx_R15 = 3.172140e+20
intM0Sdx_R22 = 3.172140e+20
intM0Sdx_R32 = 3.172138e+20
intM0Sdx_R50 = 3.172138e+20

nu_ii_R3 = 6.241776e+00
nu_ii_R5 = 8.557117e+00
nu_ii_R10 = 1.116917e+01
nu_ii_R15 = 1.220401e+01
nu_ii_R22 = 1.336318e+01
nu_ii_R32 = 1.413121e+01
nu_ii_R50 = 1.512440e+01

tau_pi_R3 = intM0dx_R3 / intM0Sdx_R3 * nu_ii_R3
tau_pi_R5 = intM0dx_R5 / intM0Sdx_R5 * nu_ii_R5
tau_pi_R10 = intM0dx_R10 / intM0Sdx_R10 * nu_ii_R10
tau_pi_R15 = intM0dx_R15 / intM0Sdx_R15 * nu_ii_R15
tau_pi_R22 = intM0dx_R22 / intM0Sdx_R22 * nu_ii_R22
tau_pi_R32 = intM0dx_R32 / intM0Sdx_R32 * nu_ii_R32
tau_pi_R50 = intM0dx_R50 / intM0Sdx_R50 * nu_ii_R50

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

beams_nu_ii_R3 = 6.687136e+00
beams_nu_ii_R5 = 9.616921e+00
beams_nu_ii_R10 = 1.262483e+01
beams_nu_ii_R15 = 1.383968e+01
beams_nu_ii_R22 = 1.511641e+01
beams_nu_ii_R32 = 1.624971e+01
beams_nu_ii_R50 = 1.760742e+01

beams_ion_M0mid_R3 = 1.277128e+19
beams_ion_M0mid_R5 = 1.503446e+19
beams_ion_M0mid_R10 = 1.921113e+19
beams_ion_M0mid_R15 = 2.107043e+19
beams_ion_M0mid_R22 = 2.282627e+19
beams_ion_M0mid_R32 = 2.440953e+19
beams_ion_M0mid_R50 = 2.612092e+19

beams_tau_pi_R3 = beams_intM0dx_R3 / beams_intM0Sdx_R3 * beams_nu_ii_R3
beams_tau_pi_R5 = beams_intM0dx_R5 / beams_intM0Sdx_R5 * beams_nu_ii_R5
beams_tau_pi_R10 = beams_intM0dx_R10 / beams_intM0Sdx_R10 * beams_nu_ii_R10
beams_tau_pi_R15 = beams_intM0dx_R15 / beams_intM0Sdx_R15 * beams_nu_ii_R15
beams_tau_pi_R22 = beams_intM0dx_R22 / beams_intM0Sdx_R22 * beams_nu_ii_R22
beams_tau_pi_R32 = beams_intM0dx_R32 / beams_intM0Sdx_R32 * beams_nu_ii_R32
beams_tau_pi_R50 = beams_intM0dx_R50 / beams_intM0Sdx_R50 * beams_nu_ii_R50

R_values = np.array([3, 5, 10, 15, 22, 32, 50])
# --- Maxwellian fit and analysis ---
def fit_func(R, k, c):
    return k * np.log10(R) + c


tau_values = np.array(
    [tau_pi_R3, tau_pi_R5, tau_pi_R10, tau_pi_R15, tau_pi_R22, tau_pi_R32, tau_pi_R50]
)
popt, pcov = curve_fit(fit_func, R_values, tau_values)
k_fit = popt[0]
c_fit = popt[1]
print(f"[Maxwellian] Fitted parameter k: {k_fit}")
print(f"[Maxwellian] Fitted parameter c: {c_fit}")

# Print the R^2 value
residuals = tau_values - fit_func(R_values, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((tau_values - np.mean(tau_values)) ** 2)
r_squared = 1 - (ss_res / ss_tot)
print(f"[Maxwellian] R^2: {r_squared}")

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
R_values_beam = np.array([3, 5, 10, 15, 22, 32, 50])
popt_beam, pcov_beam = curve_fit(fit_func, R_values_beam, beams_tau_values)
k_fit_beam = popt_beam[0]
c_fit_beam = popt_beam[1]
print(f"[Beam] Fitted parameter k: {k_fit_beam}")
print(f"[Beam] Fitted parameter c: {c_fit_beam}")

residuals_beam = beams_tau_values - fit_func(R_values_beam, *popt_beam)
ss_res_beam = np.sum(residuals_beam**2)
ss_tot_beam = np.sum((beams_tau_values - np.mean(beams_tau_values)) ** 2)
r_squared_beam = 1 - (ss_res_beam / ss_tot_beam)
print(f"[Beam] R^2: {r_squared_beam}")


plt.figure(figsize=(5, 4))

colors = plt.cm.inferno(np.linspace(0, 1, len(desired_R_values)))


# Draw an elipce centred around x=32, y=1.2 with width 10 and height 0.5
# from matplotlib.patches import Ellipse
# ellipse = Ellipse((32, 1.5), width=10, height=1.0, edgecolor='grey', facecolor='none', linestyle='--')
# plt.gca().add_patch(ellipse)
# plt.text(32, 0.8, 'WHAM', ha='center', va='bottom', fontsize=12)

plt.vlines(32, -1, 2.5, color='grey', linestyle=":", linewidth=2, alpha=0.7)
plt.text(36, 0.5, '$R_m=32$', ha='center', va='bottom', fontsize=12, rotation=90)

plt.plot(
    R_values,
    tau_values,
    marker="o",
    color=colors[1],
    markersize=8,
    linestyle="None",
    label="Maxwellian",
)

# Maxwellian fit
R_fit = np.linspace(2.5, 60, 100)
tau_fit = fit_func(R_fit, k_fit, c_fit)
plt.plot(
    R_fit,
    tau_fit,
    label=r'$ {:.2f} \log_{{10}} (R_m) $'.format(k_fit),
    color=colors[2],
    linestyle="--",
    linewidth=2,
    alpha=0.7,
)


# Beam points
plt.plot(
    R_values_beam,
    beams_tau_values,
    marker="s",
    color=colors[4],
    markersize=8,
    linestyle="None",
    label="Beam",
)
# Beam fit
R_fit_beam = np.linspace(2.5, 60, 100)
tau_fit_beam = fit_func(R_fit_beam, k_fit_beam, c_fit_beam)
plt.plot(
    R_fit_beam,
    tau_fit_beam,
    label=r'$ {:.2f} \log_{{10}} (R_m)$'.format(k_fit_beam),
    color=colors[3],
    linestyle="-.",
    linewidth=2,
    alpha=0.7,
)

# Kinetic electron
intM0dx_kin = 3.724339e+19
intM0Sdx_kin = 3.172138e+20
nu_ii_kin = 1.543262e+01
tau_pi_kin = intM0dx_kin / intM0Sdx_kin * nu_ii_kin
plt.plot(32, tau_pi_kin, marker="x", linestyle="None", color=colors[0], markersize=10, label="Beam (kinetic $e^{-}$)")

# # Add an arrow at R=32 pointing up with the text WHAM
# plt.arrow(32, 0.8, 0, 0.3, fc='black', ec='black', head_width=2.5, head_length=0.1)
# plt.text(32, 0.6, 'WHAM', ha='center', va='bottom', fontsize=12)

plt.ylim(0, 2.2)

plt.legend(loc="upper left")
plt.xlabel(r"$R_m$")
plt.ylabel(r"$ \tau_{p} \nu_{ii} $")
plt.xscale("log")
# plt.text(2.5, 0.7, r"$\mathbf{(b)}$", fontsize=18)
# plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("confinement-time.pdf")
plt.show()
