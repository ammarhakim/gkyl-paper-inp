from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

def calcNormError(M):
    err = abs((M-M[0])/M[0])
    return err

def calcL2Error(pre,p):
    # Inputs:
    # pre - prefix for file name.
    # p - polynomial order (used for assisting in distinguishing between different prefixes).
    # Outputs:
    # t - array of time data.
    # err_elcL2 - normalized error ( abs(err - err[0])/err[0] ) in the electron L^2 norm (f_e^2).
    # ion_elcL2 - normalized error ( abs(err - err[0])/err[0] ) in the proton L^2 norm (f_p^2).
    
    # Read in electron data.
    d = postgkyl.GData("%s-conservation-test-%s_elc_intL2.bp" % (pre,p))
    elcL2 = d.getValues()
    # Read in proton data.
    d = postgkyl.GData("%s-conservation-test-%s_ion_intL2.bp" % (pre,p))
    ionL2 = d.getValues()

    t = d.getGrid()[0]
    err_elcL2 = calcNormError(elcL2[:,0])
    err_ionL2 = calcNormError(ionL2[:,0])

    return t, err_elcL2, err_ionL2

# p = 2, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.4 omega_{pe}^{-1}.
t_c4, err_elcL2_c4, err_ionL2_c4 = calcL2Error("c4/c4","p2")
# p = 3, dx = 24 lambda_{D}, dv = 1 v_{th}, dt = 0.4 omega_{pe}^{-1}.
t_c9, err_elcL2_c9, err_ionL2_c9 = calcL2Error("c9/c9","p3")
# p = 2, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}.
t_c14, err_elcL2_c14, err_ionL2_c14 = calcL2Error("c14/c14","p2")
# p = 3, dx = 12 lambda_{D}, dv = 1/2 v_{th}, dt = 0.2 omega_{pe}^{-1}.
t_c15, err_elcL2_c15, err_ionL2_c15 = calcL2Error("c15/c15","p3")
# p = 2, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
t_c16, err_elcL2_c16, err_ionL2_c16 = calcL2Error("c16/c16","p2")
# p = 3, dx = 6 lambda_{D}, dv = 1/4 v_{th}, dt = 0.1 omega_{pe}^{-1}.
t_c17, err_elcL2_c17, err_ionL2_c17 = calcL2Error("c17/c17","p3")
# p = 2, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c20, err_elcL2_c20, err_ionL2_c20 = calcL2Error("c20/c20","p2")
# p = 3, dx = 3 lambda_{D}, dv = 1/8 v_{th}, dt = 0.05 omega_{pe}^{-1}.
t_c21, err_elcL2_c21, err_ionL2_c21 = calcL2Error("c21/c21","p3")

figure(1)
semilogy(t_c4, err_elcL2_c4, "bo-", label=r"$p=2, \Delta x = 24 \lambda_D, \Delta v = 1 v_{th}$")
semilogy(t_c14, err_elcL2_c14, "ko-", label=r"$p=2, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_{th}$")
semilogy(t_c16, err_elcL2_c16, "ro-", label=r"$p=2, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c16, err_elcL2_c20, "go-", label=r"$p=2, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c9, err_elcL2_c9, "bx-", label=r"$p=3, \Delta x = 24 \lambda_D, \Delta v = 1 v_{th}$")
semilogy(t_c15, err_elcL2_c15, "kx-", label=r"$p=3, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_{th}$")
semilogy(t_c17, err_elcL2_c17, "rx-", label=r"$p=3, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c21, err_elcL2_c21, "gx-", label=r"$p=3, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$L^2_{electron}$")
xlim(0, 1000)
ylim(8e-4, 3e-1)
title("Polynomial Order Comparison Electrons")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 10})
savefig("elc-L2-change.png", dpi=300)
close()

figure(2)
semilogy(t_c4, err_ionL2_c4, "bo-", label=r"$p=2, \Delta x = 24 \lambda_D, \Delta v = 1 v_{th}$")
semilogy(t_c14, err_ionL2_c14, "ko-", label=r"$p=2, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_{th}$")
semilogy(t_c16, err_ionL2_c16, "ro-", label=r"$p=2, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c16, err_ionL2_c20, "go-", label=r"$p=2, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
semilogy(t_c9, err_ionL2_c9, "bx-", label=r"$p=3, \Delta x = 24 \lambda_D, \Delta v = 1 v_{th}$")
semilogy(t_c15, err_ionL2_c15, "kx-", label=r"$p=3, \Delta x = 12 \lambda_D, \Delta v = 1/2 v_{th}$")
semilogy(t_c17, err_ionL2_c17, "rx-", label=r"$p=3, \Delta x = 6 \lambda_D, \Delta v = 1/4 v_{th}$")
semilogy(t_c21, err_ionL2_c21, "gx-", label=r"$p=3, \Delta x = 3 \lambda_D, \Delta v = 1/8 v_{th}$")
xlabel(r"$t (\omega_{pe}^{-1})$")
ylabel(r"$L^2_{proton}$")
xlim(0, 1000)
ylim(8e-4, 3e-1)
title("Polynomial Order Comparison Protons")
legend(loc="lower right", ncol=2, mode="expand", prop={"size": 10})
savefig("ion-L2-change.png", dpi=300)
close()
