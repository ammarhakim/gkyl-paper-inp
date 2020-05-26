from pylab import *
import postgkyl as pg
import math
import numpy as np
style.use("../postgkyl.mplstyle")

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i, figsize=(14,8))

    #Selective read-in
    #plotting at vy = 0, coord3 = 12, y = -1.7 rhoi, coord1 = 94
    data_distf = pg.GData("kink_ion_%d.bp" % fr, coord1 = 94, coord3 = 12)
    dg_distf = pg.data.GInterpModal(data_distf, 2, "ms", 9)
    XX, distfIon = dg_distf.interpolate(0)
    #center the grid values
    for d in range(4):
        XX[d] = 0.5*(XX[d][:-1] + XX[d][1:])

    #computing ion larmor radius
    omegaCi = 0.0055277079839257
    vtIon = 0.0316227833333333
    rhoi = vtIon/omegaCi
    #want to plot at x = 2.3 rhoi, y = -1.7 rhoi at edge of current sheet
    #corresponds to x = 330, y = 0 after interpolation (and selective read-in)
    #also plotting at vy = 0, vy = 0 after interpolation (and selective read-in)
    subplot(1, 2, 1)
    pcolormesh(XX[0]/rhoi, XX[2]/vtIon, distfIon[:,0,:,0,0].transpose(), shading="gouraud")
    colorbar()
    xlabel(r"$X(\rho_i)$")
    ylabel(r"$V_X(v_{th_p})$")
    ylim(-4, 4)
    title(r"$t=6 \Omega_{ci}^{-1}$")

    subplot(1, 2, 2)
    plot(XX[2]/vtIon, distfIon[990,0,:,0,0],"k")
    axvline(x=-1.0, color = "r", linestyle = "-", label="Initial drift velocity")
    axvline(x=-0.55, color = "g", linestyle = "--", label="Phase velocity")
    xlim(-4, 4)
    xlabel(r"$V_X(v_{th_p})$")
    ylabel(r"$f_p(v_x)$")
    title(r"$t=6 \Omega_{ci}^{-1}$")
    legend(loc="best")
    
    tight_layout()
    savefig("proton-distribution-function-lhdi.pdf")
    
for i in range(30,31):
    plotFig(i,i)
    
show()
