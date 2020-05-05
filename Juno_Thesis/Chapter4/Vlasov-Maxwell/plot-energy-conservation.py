from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

# Normalization parameters for computing energies in code units.
epsilon0 = 1.0
elcMass = 1.0
ionMass = 1836.153

def calcNormError(M):
    err = abs((M-M[0])/M[0])
    return err

def calcEnergyError(pre,p):
    # Inputs:
    # pre - prefix for file name.
    # p - polynomial order (used for assisting in distinguishing between different prefixes).
    # Outputs:
    # t - array of time data.
    # err - normalized error ( abs(err - err[0])/err[0] ) in the total energy.
    
    # Read in electron data.
    d = postgkyl.GData("%s-conservation-test-%s_elc_intM2Thermal.bp" % (pre,p))
    elcM2Thermal =  0.5*elcMass*d.getValues()
    d = postgkyl.GData("%s-conservation-test-%s_elc_intM2Flow.bp" % (pre,p))
    elcM2Flow =  0.5*elcMass*d.getValues()
    # Read in proton data.
    d = postgkyl.GData("%s-conservation-test-%s_ion_intM2Thermal.bp" % (pre,p))
    ionM2Thermal =  0.5*ionMass*d.getValues()
    d = postgkyl.GData("%s-conservation-test-%s_ion_intM2Flow.bp" % (pre,p))
    ionM2Flow =  0.5*ionMass*d.getValues()
    # Read in electromagnetic field data.
    d = postgkyl.GData("%s-conservation-test-%s_fieldEnergy.bp" % (pre,p))
    fieldEnergy =  0.5*epsilon0*d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(elcM2Flow[:,0]+elcM2Thermal[:,0]+ionM2Flow[:,0]+ionM2Thermal[:,0]+fieldEnergy[:,0])

    return t, err

# p = 1, dx = 24 lambda_D, dv = 1 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c1, err_c1 = calcEnergyError("c1/c1","p1")

# p = 1, dx = 24 lambda_D, dv = 1/2 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c2, err_c2 = calcEnergyError("c2/c2","p1")

# p = 1, dx = 24 lambda_D, dv = 1/4 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c3, err_c3 = calcEnergyError("c3/c3","p1")

# p = 2, dx = 24 lambda_D, dv = 1 v_{th}, dt = 2/5 omega_{pe}^{-1}.
t_c4, err_c4 = calcEnergyError("c4/c4","p2")

# p = 2, dx = 24 lambda_D, dv = 1 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c5, err_c5 = calcEnergyError("c5/c5","p2")

# p = 2, dx = 24 lambda_D, dv = 1 v_{th}, dt = 1/10 omega_{pe}^{-1}.
t_c6, err_c6 = calcEnergyError("c6/c6","p2")

# p = 2, dx = 12 lambda_D, dv = 1 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c7, err_c7 = calcEnergyError("c7/c7","p2")

# p = 2, dx = 6 lambda_D, dv = 1 v_{th}, dt = 1/10 omega_{pe}^{-1}.
t_c8, err_c8 = calcEnergyError("c8/c8","p2")

# p = 3, dx = 24 lambda_D, dv = 1 v_{th}, dt = 2/5 omega_{pe}^{-1}.
t_c9, err_c9 = calcEnergyError("c9/c9","p3")

# p = 3, dx = 24 lambda_D, dv = 1 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c10, err_c10 = calcEnergyError("c10/c10","p3")

# p = 3, dx = 24 lambda_D, dv = 1 v_{th}, dt = 1/10 omega_{pe}^{-1}.
t_c11, err_c11 = calcEnergyError("c11/c11","p3")

# p = 3, dx = 12 lambda_D, dv = 1 v_{th}, dt = 1/5 omega_{pe}^{-1}.
t_c12, err_c12 = calcEnergyError("c12/c12","p3")

# p = 3, dx = 6 lambda_D, dv = 1 v_{th}, dt = 1/10 omega_{pe}^{-1}.
t_c13, err_c13 = calcEnergyError("c13/c13","p3")

figure(1)
semilogy(t_c4, err_c4, "bo--", label=r"$p=2, \Delta t = 0.4 \omega_{pe}^{-1}$")
semilogy(t_c5, err_c5, "ko--", label=r"$p=2, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c6, err_c6, "ro--", label=r"$p=2, \Delta t = 0.1 \omega_{pe}^{-1}$")
semilogy(t_c9, err_c9, "b*-", label=r"$p=3, \Delta t = 0.4 \omega_{pe}^{-1}$")
semilogy(t_c10, err_c10, "k*-", label=r"$p=3, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c11, err_c11, "r*-", label=r"$p=3, \Delta t = 0.1 \omega_{pe}^{-1}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{E}$")
xlim(0, 1000)
ylim(5e-10, 5e-7)
title("Polynomial Order Comparison")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("energy-conservation-polyOrder-comparison.png", dpi=300)
close()

figure(2)
semilogy(t_c1, err_c1, "bs:", label=r"$p=1, \Delta v = 1 v_t, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c2, err_c2, "ks:", label=r"$p=1, \Delta v = 1/2 v_t, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c3, err_c3, "rs:", label=r"$p=1, \Delta v = 1/4 v_t, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c5, err_c5, "ro--", label=r"$p=2, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c10, err_c10, "r*-", label=r"$p=3, \Delta t = 0.2 \omega_{pe}^{-1}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{E}$")
xlim(0, 1000)
ylim(1e-10, 1e-5)
title("Polynomial Order Comparison")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("energy-conservation-polyOrder-1.png", dpi=300)
close()

figure(3)
semilogy(t_c4, err_c4, "bo--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.4 \omega_{pe}^{-1}$")
semilogy(t_c5, err_c5, "ko--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c6, err_c6, "ro--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.1 \omega_{pe}^{-1}$")
semilogy(t_c7, err_c7, "k*--", label=r"$\Delta x = 12 \lambda_D, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c8, err_c8, "r*--", label=r"$\Delta x = 6 \lambda_D, \Delta t = 0.1 \omega_{pe}^{-1}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{E}$")
xlim(0, 1000)
ylim(5e-10, 5e-7)
title("Polynomial Order 2")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("energy-conservation-polyOrder-2.png", dpi=300)
close()

figure(4)
semilogy(t_c9, err_c9, "bo--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.4 \omega_{pe}^{-1}$")
semilogy(t_c10, err_c10, "ko--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c11, err_c11, "ro--", label=r"$\Delta x = 24 \lambda_D, \Delta t = 0.1 \omega_{pe}^{-1}$")
semilogy(t_c12, err_c12, "k*--", label=r"$\Delta x = 12 \lambda_D, \Delta t = 0.2 \omega_{pe}^{-1}$")
semilogy(t_c13, err_c13, "r*--", label=r"$\Delta x = 6 \lambda_D, \Delta t = 0.1 \omega_{pe}^{-1}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{E}$")
xlim(0, 1000)
ylim(5e-10, 5e-7)
title("Polynomial Order 3")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("energy-conservation-polyOrder-3.png", dpi=300)
close()

print("Total Energy Orders of Convergence")
print("=============================================================================")
print("Comparing $p=2, dt = 0.4 \omega_{pe}^{-1}$ to $p=2, dt = 0.2 \omega_{pe}^{-1}$")
print(err_c4[-1]/err_c5[-1])
print("$\log_2$ of the comparison")
print(log(err_c4[-1]/err_c5[-1])/log(2))
print(r"Comparing $p=2, dt = 0.2 \omega_{pe}^{-1}$ to $p=2, dt = 0.1 \omega_{pe}^{-1}$")
print(err_c5[-1]/err_c6[-1])
print("$\log_2$ of the comparison")
print(log(err_c5[-1]/err_c6[-1])/log(2))
print("=============================================================================")
print("Comparing $p=3, dt = 0.4 \omega_{pe}^{-1}$ to $p=3, dt = 0.2 \omega_{pe}^{-1}$")
print(err_c9[-1]/err_c10[-1])
print("$\log_2$ of the comparison")
print(log(err_c9[-1]/err_c10[-1])/log(2))
print("Comparing $p=3, dt = 0.4 \omega_{pe}^{-1}$ to $p=3, dt = 0.2 \omega_{pe}^{-1}$")
print(err_c10[-1]/err_c11[-1])
print("$\log_2$ of the comparison")
print(log(err_c10[-1]/err_c11[-1])/log(2))
print("=============================================================================")
