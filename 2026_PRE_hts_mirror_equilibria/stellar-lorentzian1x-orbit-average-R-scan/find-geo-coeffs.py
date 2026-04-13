import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def B_mag(z, Bhat, gamma):
  z_m = 0.98
  lor_1 = 1/(1 + ((z - z_m)/(gamma))**2)
  lor_2 = 1/(1 + ((z + z_m)/(gamma))**2)
  Bmag = Bhat/np.pi/gamma*(lor_1 + lor_2)
  return Bmag

Bhat=6.5
gamma=0.124

z = np.linspace(-2.0, 2.0, 1000)
B = B_mag(z, Bhat, gamma)

R = B_mag(0.98, Bhat, gamma)/B_mag(0.0, Bhat, gamma)
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
    
    results.append({
        'R_target': desired_R,
        'Bhat': Bhat_opt,
        'gamma': gamma_opt,
        'B0': B0_result,
        'R_actual': R_result
    })
    
    print(f"{desired_R:6.1f}  {Bhat_opt:10.6f}  {gamma_opt:10.6f}  {B0_result:10.6f}  {R_result:10.6f}")

# Plot comparison
plt.figure(figsize=(12, 8))

# Plot initial profile
# B_initial = B_mag(z, Bhat, gamma)
# plt.plot(z, B_initial, label=f'Initial: R={R:.2f}', linewidth=2.5, color='black', linestyle='-')

# Plot optimized profiles
colors = plt.cm.Set2(np.linspace(0, 1, len(desired_R_values)))
for i, result in enumerate(results):
    B_optimized = B_mag(z, result['Bhat'], result['gamma'])
    plt.plot(z, B_optimized, 
             label=f"R={result['R_target']:.0f}", 
             linewidth=2, color=colors[i])

plt.xlabel("z", fontsize=12)
plt.ylabel("B", fontsize=12)
plt.title(f"Magnetic Field Profiles for Different Mirror Ratios (B(0) = {B0_target:.2f} fixed)", fontsize=14)
plt.legend(fontsize=9, loc='best')
# plt.grid(alpha=0.3)
plt.yscale('log')
plt.tight_layout()
plt.show()

# Print arrays for easy copying
print("\n--- Arrays for copy-paste ---")
print("Bhat_values = [", end="")
print(", ".join([f"{r['Bhat']:.6f}" for r in results]), end="")
print("]")
print("gamma_values = [", end="")
print(", ".join([f"{r['gamma']:.6f}" for r in results]), end="")
print("]")

