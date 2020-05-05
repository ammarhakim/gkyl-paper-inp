from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

# Normalization parameters for computing momenta in code units.
elcMass = 1.0
ionMass = 1836.153

def calcNormError(M):
    err = abs((M-M[0])/M[0])
    return err

def calcMomentumError(pre,p):
    # Inputs:
    # pre - prefix for file name.
    # p - polynomial order (used for assisting in distinguishing between different prefixes).
    # Outputs:
    # t - array of time data.
    # err - normalized error ( abs(err - err[0])/err[0] ) in the total momentum.
    
    # Read in electron data.
    d = postgkyl.GData("%s-conservation-test-%s_elc_intM1i.bp" % (pre,p))
    elcM1i =  elcMass*d.getValues()
    # Read in proton data.
    d = postgkyl.GData("%s-conservation-test-%s_ion_intM1i.bp" % (pre,p))
    ionM1i =  ionMass*d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(elcM1i[:,0]+ionM1i[:,0])

    return t, err

# p = 2, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
t_c16, err_c16 = calcMomentumError("c16/c16","p2")
# p = 3, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
t_c17, err_c17 = calcMomentumError("c17/c17","p3")
# p = 2, dx = 3 lambda_{D}, dv = 1/4 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c18, err_c18 = calcMomentumError("c20/c20","p2")
# p = 3, dx = 3 lambda_{D}, dv = 1/4 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c19, err_c19 = calcMomentumError("c21/c21","p3")
# p = 2, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c20, err_c20 = calcMomentumError("c20/c20","p2")
# p = 3, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c21, err_c21 = calcMomentumError("c21/c21","p3")
# p = 2, dx = 1.5 lambda_{D}, dv = 1/4 v_{th}, dt = 0.025 omega_{pe}^{-1}.
t_c22, err_c22 = calcMomentumError("c22/c22","p2")
# p = 3, dx = 1.5 lambda_{D}, dv = 1/4 v_{th}, dt = 0.025 omega_{pe}^{-1}.
t_c23, err_c23 = calcMomentumError("c23/c23","p3")
# p = 2, dx = 1.5 lambda_{D}, dv = 1/8 v_{th}, dt = 0.025 omega_{pe}^{-1}.
t_c24, err_c24 = calcMomentumError("c24/c24","p2")
# p = 3, dx = 1.5 lambda_{D}, dv = 1/8 v_{th}, dt = 0.025 omega_{pe}^{-1}.
t_c25, err_c25 = calcMomentumError("c25/c25","p3")
# p = 2, dx = 0.75 lambda_{D}, dv = 1/4 v_{th}, dt = 0.0125 omega_{pe}^{-1}.
t_c26, err_c26 = calcMomentumError("c26/c26","p2")
# p = 3, dx = 0.75 lambda_{D}, dv = 1/4 v_{th}, dt = 0.0125 omega_{pe}^{-1}.
t_c27, err_c27 = calcMomentumError("c27/c27","p3")
# p = 2, dx = 0.75 lambda_{D}, dv = 1/8 v_{th}, dt = 0.0125 omega_{pe}^{-1}.
t_c28, err_c28 = calcMomentumError("c28/c28","p2")
# p = 3, dx = 0.75 lambda_{D}, dv = 1/8 v_{th}, dt = 0.0125 omega_{pe}^{-1}.
t_c29, err_c29 = calcMomentumError("c29/c29","p3")
# p = 2, dx = 0.375 lambda_{D}, dv = 1/4 v_{th}, dt = 0.00625 omega_{pe}^{-1}.
t_c30, err_c30 = calcMomentumError("c30/c30","p2")
# p = 3, dx = 0.375 lambda_{D}, dv = 1/4 v_{th}, dt = 0.00625 omega_{pe}^{-1}.
t_c31, err_c31 = calcMomentumError("c31/c31","p3")
# p = 2, dx = 0.375 lambda_{D}, dv = 1/8 v_{th}, dt = 0.00625 omega_{pe}^{-1}.
t_c32, err_c32 = calcMomentumError("c32/c32","p2")
# p = 3, dx = 0.375 lambda_{D}, dv = 1/8 v_{th}, dt = 0.00625 omega_{pe}^{-1}.
t_c33, err_c33 = calcMomentumError("c33/c3","p3")

figure(1)
semilogy(t_c16, err_c16, "bo--", label=r"$\Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c18, err_c18, "ko-", label=r"$\Delta x = 3 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c20, err_c20, "k*--", label=r"$\Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c22, err_c22, "ro-", label=r"$\Delta x = 1.5 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c24, err_c24, "r*--", label=r"$\Delta x = 1.5 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c26, err_c26, "go-", label=r"$\Delta x = 0.75 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c28, err_c28, "g*--", label=r"$\Delta x = 0.75 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c30, err_c30, "mo-", label=r"$\Delta x = 0.375 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c32, err_c32, "m*--", label=r"$\Delta x = 0.375 \lambda_D, \Delta v = 1/8 v_{th}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{M}$")
xlim(0, 1000)
ylim(1e-12, 2e-6)
title("Polynomial Order 2")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("momentum-conservation-polyOrder-2.png", dpi=300)
close()

figure(2)
semilogy(t_c17, err_c17, "bo--", label=r"$\Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c19, err_c19, "ko-", label=r"$\Delta x = 3 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c21, err_c21, "k*--", label=r"$\Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c23, err_c23, "ro-", label=r"$\Delta x = 1.5 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c25, err_c25, "r*--", label=r"$\Delta x = 1.5 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c27, err_c27, "go-", label=r"$\Delta x = 0.75 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c29, err_c29, "g*--", label=r"$\Delta x = 0.75 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c31, err_c31, "mo-", label=r"$\Delta x = 0.375 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c33, err_c33, "m*--", label=r"$\Delta x = 0.375 \lambda_D, \Delta v = 1/8 v_{th}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$\Delta \mathcal{M}$")
xlim(0, 1000)
ylim(2e-14, 1e-6)
title("Polynomial Order 3")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 11})
savefig("momentum-conservation-polyOrder-3.png", dpi=300)
close()

print("Total Momentum Orders of Convergence")
print("=============================================================================")
print("Comparing $p = 2, dx = 6 \lambda_{D}, dv = 1/4 v_{th}$ to $p = 2, dx = 3 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c16[-1]/err_c20[-1])
print("$\log_2$ of the comparison")
print(log(err_c16[-1]/err_c20[-1])/log(2))
print(r"Comparing $p = 2, dx = 3 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 2, dx = 1.5 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c20[-1]/err_c24[-1])
print("$\log_2$ of the comparison")
print(log(err_c20[-1]/err_c24[-1])/log(2))
print("Comparing $p = 2, dx = 1.5 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 2, dx = 0.75 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c24[-1]/err_c28[-1])
print("$\log_2$ of the comparison")
print(log(err_c24[-1]/err_c28[-1])/log(2))
print(r"Comparing $p = 2, dx = 0.75 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 2, dx = 0.375 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c28[-1]/err_c32[-1])
print("$\log_2$ of the comparison")
print(log(err_c28[-1]/err_c32[-1])/log(2))
print("=============================================================================")
print("Comparing $p = 3, dx = 6 \lambda_{D}, dv = 1/4 v_{th}$ to $p = 3, dx = 3 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c17[-1]/err_c21[-1])
print("$\log_2$ of the comparison")
print(log(err_c17[-1]/err_c21[-1])/log(2))
print(r"Comparing $p = 3, dx = 3 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 3, dx = 1.5 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c21[-1]/err_c25[-1])
print("$\log_2$ of the comparison")
print(log(err_c21[-1]/err_c25[-1])/log(2))
print("Comparing $p = 3, dx = 1.5 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 3, dx = 0.75 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c25[-1]/err_c29[-1])
print("$\log_2$ of the comparison")
print(log(err_c15[-1]/err_c29[-1])/log(2))
print(r"Comparing $p = 3, dx = 0.75 \lambda_{D}, dv = 1/8 v_{th}$ to $p = 3, dx = 0.375 \lambda_{D}, dv = 1/8 v_{th}$")
print(err_c29[-1]/err_c33[-1])
print("$\log_2$ of the comparison")
print(log(err_c29[-1]/err_c33[-1])/log(2))
print("=============================================================================")
