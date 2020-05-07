from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

# Normalization parameters for checking Gauss' Law in code units.
epsilon0 = 1.0
elcCharge = -1.0
ionCharge = 1.0

def calcDivError(pre,p,numInterp,fr):
    # Inputs:
    # pre - prefix for file name.
    # p - polynomial order (used for assisting in distinguishing between different prefixes).
    # numInterp - number of points on which to interpolate DG basis functions.
    # fr - frame number to read-in.
    # Outputs:
    # Xc - cell center coordinates for interpolated fields.
    # E_dx - x-derivative of the x-electric field.
    # rho_c - charge density divide by permittivity of free space (epsilon0).
    
    # Read in electric field data of specified frame.
    d = postgkyl.GData("%s-conservation-test-%s_field_%d.bp" % (pre,p,fr))
    dg = postgkyl.GInterpModal(d, 2, "ms", numInterp)
    # Interpolate the zeroth component (Ex).
    X, Ex = dg.interpolate(0)

    # Read in and interpolate electron number density of specified frame.
    d = postgkyl.GData("%s-conservation-test-%s_elc_M0_%d.bp" % (pre,p,fr))
    dg = postgkyl.GInterpModal(d, 2, "ms", numInterp)
    X, n_e = dg.interpolate()

    # Read in and interpolate proton number density of specified frame.
    d = postgkyl.GData("%s-conservation-test-%s_ion_M0_%d.bp" % (pre,p,fr))
    dg = postgkyl.GInterpModal(d, 2, "ms", numInterp)
    X, n_i = dg.interpolate()

    # Dividing charge density by epsilon0 here for comparison to divergence of the electric field.
    rho_c = (elcCharge*n_e + ionCharge*n_i)/epsilon0

    Xn = X[0]; dx = Xn[1]-Xn[0]
    # Cell-center coordinates.
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)
    # Compute derivative of Ex in x.
    E_dx = gradient(Ex,  dx , edge_order=2, axis=0)
    return Xc, E_dx[:,0], rho_c[:,0]

# p = 2, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}.
X_c14, E_dx_c14, rho_c_c14 = calcDivError("c14/c14","p2", 9, 1)
# p = 3, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}.
X_c15, E_dx_c15, rho_c_c15 = calcDivError("c15/c15","p3", 9, 1)
# p = 2, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
X_c16, E_dx_c16, rho_c_c16 = calcDivError("c16/c16","p2", 9, 1)
# p = 3, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
X_c17, E_dx_c17, rho_c_c17 = calcDivError("c17/c17","p3", 9, 1)
# p = 2, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
X_c20, E_dx_c20, rho_c_c20 = calcDivError("c20/c20","p2", 9, 1)
# p = 3, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
X_c21, E_dx_c21, rho_c_c21 = calcDivError("c21/c21","p3", 9, 1)

figure(1)
plot(X_c14, E_dx_c14, "b--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_t$")
plot(X_c14, rho_c_c14, "b*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_t$")
plot(X_c16, E_dx_c16, "k--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_t$")
plot(X_c16, rho_c_c16, "k*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_t$")
plot(X_c20, E_dx_c20, "r--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_t$")
plot(X_c20, rho_c_c20, "r*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_t$")
xlabel(r"$X (\lambda_D)$")
xlim(0, 96)
title(r"Polynomial Order 2, $t=1000 \omega_{pe}^{-1}$")
legend(loc="upper right", ncol=1, prop={"size": 11})
savefig("div-E-vs-rhoc-p2.png", dpi=300)
close()

figure(2)
plot(X_c15, E_dx_c15, "b--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_t$")
plot(X_c15, rho_c_c15, "b*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_t$")
plot(X_c17, E_dx_c17, "k--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_t$")
plot(X_c17, rho_c_c17, "k*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_t$")
plot(X_c21, E_dx_c21, "r--", label=r"$\frac{\partial E_x}{\partial x}, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_t$")
plot(X_c21, rho_c_c21, "r*", label=r"$\frac{\rho_c}{\epsilon_0}, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_t$")
xlabel(r"$X (\lambda_D)$")
xlim(0, 96)
ylim(-0.006, 0.01)
title(r"Polynomial Order 3, $t=1000 \omega_{pe}^{-1}$")
legend(loc="upper right", ncol=2, mode="expand", prop={"size": 11})
savefig("div-E-vs-rhoc-p3.png", dpi=300)
close()
