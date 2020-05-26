from pylab import *
import postgkyl as pg
import math
import numpy as np
import scipy.optimize as opt
style.use("../postgkyl.mplstyle")

# --------------------------------------------------------------------
# Growth rate fitting stuff ------------------------------------------
def exp2(x, a, b):
    """Define custom exponential a*exp(2b*x)

    Parameters:
    x -- independent variable
    a -- scaling parameter
    b -- growth rate

    Notes:
    Energy (quantity^2) is often used for the growth-rate study,
    therefore the factor 2
    """
    return a*np.exp(2*b*x)

def fitGrowth(x, y, function=exp2, minN=100, maxN=None, p0=(1, 0.1)):
    """Fit function to continuously increasing region of data

    Parameters:
    x -- independet variable
    y -- dependent variable
    minN -- minimal number of fitted points (default: 100)
    maxN -- maximal number of fitted points (default: full length)
    function -- function to fit (default: exp2)
    p0 -- initial guess (default: 1, 0.1)

    Notes:
    The best is determined based on the coeficient of determination,
      R^2 https://en.wikipedia.org/wiki/Coefficient_of_determination
    """
    bestR2 = 0
    if maxN is None:
        maxN = len(x)
    bestParams = p0

    print("fitGrowth: fitting region {:d} -> {:d}".format(minN, maxN))
    for n in np.linspace(minN, maxN-1, maxN-minN):
        n = int(n)
        xn = x[0 : n]  # continuously increasing fitting region
        yn = y[0 : n]
        try:
            params, cov = opt.curve_fit(function, xn, yn, bestParams)
            residual = yn - function(xn, *params)
            ssRes = np.sum(residual**2)
            ssTot = np.sum((yn - np.mean(yn))**2)
            R2 = 1 - ssRes/ssTot
            if R2 > bestR2:
                bestR2 = R2
                bestParams = params
                bestN = n
            percent = float(n-minN)/(maxN-minN)*100
            progress = "[" + int(percent/10)*"=" + (10-int(percent/10))*" " + "]"
            sys.stdout.write(
                "\rgamma = {:+.4e} (best {:+.5e}) R^2 = {:.3e}   {:6.2f}% done {}".format(params[1], bestParams[1], R2, percent, progress))
            sys.stdout.flush()
        except RuntimeError:
            print("fitGrowth: curve_fit failed for N = {}".format(n))

    print("\rgamma = {:+.4e} (best {:+.5e}) R^2 = {:.3e}   {:6.2f}% done {}".
          format(params[1], bestParams[1], R2, 100, "[==========]"))
    return bestParams, bestR2, bestN

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i, figsize=(14,8))
    data_fieldEnergy = pg.GData("kink_fieldEnergy_")
    t = data_fieldEnergy.getGrid()[0]
    fieldEnergy = data_fieldEnergy.getValues()
    ExEnergy = 0.5*fieldEnergy[:,0]

    omegaCi = 0.0055277079839257
    t = t*0.0055277079839257
    #param, R2, N = fitGrowth(t[40000:50000], ExEnergy[40000:50000])
    subplot(1, 2, 1)
    semilogy(t, ExEnergy, "k")
    xlabel("$t (\Omega_{ci}^{-1})$")
    ylabel("$\int E_x^2$")
    xlim(0, 8)
    ylim(1e-12, 1e-2)


    data_field = pg.GData("kink_field_%d.bp" % fr)
    dg_field = pg.data.GInterpModal(data_field, 2, "ms")
    XX, Ex = dg_field.interpolate(0)
    #center the grid values
    for d in range(2):
        XX[d] = 0.5*(XX[d][:-1] + XX[d][1:])
    #normalization for electric field
    ionMass = 36
    norm = (omegaCi*ionMass)**2/math.sqrt(ionMass)
    #computing ion larmor radius
    vtIon = 0.0316227833333333
    rhoi = vtIon/omegaCi
    subplot(1, 2, 2)
    pcolormesh(XX[0]/rhoi, XX[1]/rhoi, Ex[:,:,0].transpose()/norm, cmap = "seismic", shading="gouraud")
    colorbar(label="$E_x$")
    xlabel(r"$X(\rho_i)$")
    ylabel(r"$Y(\rho_i)$")
    title(r"$t=6 \Omega_{ci}^{-1}$")

    tight_layout()
    savefig("lhdi-ex-energy-and-ex.pdf")
    
for i in range(30,31):
    plotFig(i,i)
    
show()

