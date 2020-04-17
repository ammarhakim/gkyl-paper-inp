#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import postgkyl as pg
import scipy.optimize as opt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#tIdxs = (117, 142)
tIdxs = (118,)
q = -1
x = 80

def doubleMaxwell(v, n1, u1, t1, n2, u2, t2):
    return n1*np.exp(-(v-u1)**2/(2*t2**2)) + n2*np.exp(-(v-u2)**2/(2*t2**2))

def fit(x, y, p0=(2, 0.3, 0.1, 2, -0.3, 0.1)):
     return opt.curve_fit(doubleMaxwell, x, y, p0)

def plotFrame(tIdx):
    dD = pg.GData('ck_1X2V_distfElc_{:03d}.h5'.format(int(tIdx)))
    time = dD.time
    dg = pg.data.GInterpNodal(dD, 2, 'ns')
    c3D, distf = dg.interpolate(0)

    VX, VY = np.meshgrid(c3D[1], c3D[2], indexing='ij')

    dx = c3D[0][1] - c3D[0][0]
    dvx = c3D[1][1] - c3D[1][0]
    dvy = c3D[2][1] - c3D[2][0]

    boundaryIdx = np.zeros(distf.shape[0])
    vBoundary = np.zeros(distf.shape[0])
    numDensity1 = np.zeros(distf.shape[0])
    numDensity2 = np.zeros(distf.shape[0])
    ux1 = np.zeros(distf.shape[0])
    ux2 = np.zeros(distf.shape[0])
    uy1 = np.zeros(distf.shape[0])
    uy2 = np.zeros(distf.shape[0])

    distfXVX = np.sum(distf, axis=2)*dvy
    distfXVY = np.sum(distf, axis=1)*dvx

    params = (2, 0.3, 0.1, 2, -0.3, 0.1)
    for i in range(distf.shape[0]):
        tmp = np.sum(distf[i, :, :], axis=0)*dvx
        params, c = fit(c3D[2], tmp, p0=params)
        idx1 = pg.tools.findNearestIdx(c3D[2], params[1])
        idx2 = pg.tools.findNearestIdx(c3D[2], params[4])
        boundary = np.argmin(doubleMaxwell(c3D[2], *params)[idx2:idx1])+idx2
        boundaryIdx[i] = boundary
        vBoundary[i] = c3D[2][int(boundary)]

        numDensity1[i] = np.sum(distf[i, :, slice(boundary, 120)])*dvx*dvy
        numDensity2[i] = np.sum(distf[i, :, slice(0, boundary)])*dvx*dvy

        ux1[i] = np.sum(distf[i, :, slice(boundary, 120)]*
                        VX[:, slice(boundary, 120)])*dvx*dvy/numDensity1[i]
        ux2[i] = np.sum(distf[i, :, slice(0, boundary)]*
                        VX[:, slice(0, boundary)])*dvx*dvy/numDensity2[i]
        uy1[i] = np.sum(distf[i, :, slice(boundary, 120)]*
                        VY[:, slice(boundary, 120)])*dvx*dvy/numDensity1[i]
        uy2[i] = np.sum(distf[i, :, slice(0, boundary)]*
                         VY[:, slice(0, boundary)])*dvx*dvy/numDensity2[i]

    dEM = pg.data.GData('ck_1X2V_em_{:03d}.h5'.format(int(tIdx)))
    dg = pg.data.GInterpNodal(dEM, 2, 'ns')
    x, Ex = dg.interpolate(0)
    x, Ey = dg.interpolate(1)
    x, Bz = dg.interpolate(5)

    dND = pg.data.GData('ck_1X2V_numDensityElc_{:03d}.h5'.
                        format(int(tIdx)))
    dg = pg.data.GInterpNodal(dND, 2, 'ns')
    xK, numDensity = dg.interpolate(0)

    ExEnergy = pg.data.GHistoryData('ck_1X2V_ExEnergy')
    BzEnergy = pg.data.GHistoryData('ck_1X2V_BzEnergy')

    x[0] = x[0]*0.4/2/np.pi

    ff1 = -uy1*Bz
    ff2 = -uy2*Bz

    phiE = np.zeros(len(Ex))
    phiFF1 = np.zeros(len(Ex))
    phiFF2 = np.zeros(len(Ex))
    for i in np.arange(len(Ex)-1)+1:
        phiE[i] = phiE[i-1] + Ex[i]*dx
        phiFF1[i] = phiFF1[i-1] - ff1[i]*dx
        phiFF2[i] = phiFF2[i-1] - ff2[i]*dx

    

    # plotting
    fig, ax = plt.subplots(5, 2, figsize=(12, 16))

    ax[0, 0].plot(x[0], numDensity, color='C3', label='$n$')
    ax[0, 0].plot(x[0], numDensity1, color='C0', label='$n^+$')
    ax[0, 0].plot(x[0], numDensity2, color='C1', label='$n^-$')
    ax[0, 0].set_ylim((0, 1.4))
    ax[0, 0].xaxis.set_ticklabels([])
    ax[0, 0].set_ylabel('$n$')
    ax[0, 0].text(0.04, 0.96, 'a)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[0, 0].transAxes)

    ax[0, 1].plot(ExEnergy.time, ExEnergy.values, color='C2', label='$E_x^2$')
    ax[0, 1].plot(BzEnergy.time, BzEnergy.values, color='C4', label='$B_z^2$')
    ax[0, 1].set_autoscale_on(False)
    ax[0, 1].plot((tIdx*0.5, tIdx*0.5), (-1, 1),
                  '--', color='C9', zorder=0, linewidth=1.5)
    ax[0, 1].xaxis.set_ticklabels([])
    ax[0, 1].set_ylabel('$energies$')
    ax[0, 1].text(0.04, 0.96, 'b)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[0, 1].transAxes)

    ax[1, 0].plot(x[0], ux1, color='C0', label='$u_x^+$')
    ax[1, 0].plot(x[0], ux2, color='C1', label='$u_x^-$')
    ax[1, 0].set_ylim((-0.5, 0.5))
    ax[1, 0].xaxis.set_ticklabels([])
    ax[1, 0].set_yticks((-0.25, 0.0, 0.25))
    ax[1, 0].set_ylabel('$u_x$')
    ax[1, 0].text(0.04, 0.96, 'c)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[1, 0].transAxes)

    ax[1, 1].plot(x[0], uy1, color='C0', label='$u_y^+$')
    ax[1, 1].plot(x[0], uy2, color='C1', label='$u_y^-$')
    ax[1, 1].set_ylim((-0.6, 0.6))
    ax[1, 1].xaxis.set_ticklabels([])
    ax[1, 1].set_ylabel('$u_y$')
    ax[1, 1].text(0.04, 0.96, 'd)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[1, 1].transAxes)

    ax[2, 0].pcolormesh(x[0], c3D[1], distfXVX.transpose())
    ax[2, 0].xaxis.set_ticklabels([])
    ax[2, 0].set_ylabel('$v_x$')
    ax[2, 0].text(0.04, 0.96, 'e)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[2, 0].transAxes, color='w')

    ax[2, 1].pcolormesh(x[0], c3D[2], distfXVY.transpose())
    ax[2, 1].plot(x[0], vBoundary, '--', color='C7', linewidth=1.5)
    ax[2, 1].xaxis.set_ticklabels([])
    ax[2, 1].set_ylabel('$v_y$')
    ax[2, 1].text(0.04, 0.96, 'f)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[2, 1].transAxes, color='w')

    ax[3, 0].plot(x[0], -Ex, color='C2', label='$E_x$')
    ax[3, 0].plot(x[0], ff1, color='C0', label='$u_y^+B_z$')
    ax[3, 0].plot(x[0], ff2, color='C1', label='$u_y^-B_z$')
    ax[3, 0].set_ylim((-0.1, 0.1))
    ax[3, 0].xaxis.set_ticklabels([])
    ax[3, 0].set_yticks((-0.05, 0.0, 0.05))
    ax[3, 0].set_ylabel('$F_x$')
    ax[3, 0].text(0.04, 0.96, 'g)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[3, 0].transAxes)

    ax[3, 1].plot(x[0], ff1-Ex, color='C0', label='$E_x+u_y^+B_z$')
    ax[3, 1].plot(x[0], ff2-Ex, color='C1', label='$E_x+u_y^-B_z$')
    ax[3, 1].set_ylim((-0.15, 0.15))
    ax[3, 1].xaxis.set_ticklabels([])
    ax[3, 1].set_yticks((-0.05, -0.1, 0.0, 0.1, 0.05))
    ax[3, 1].set_ylabel('$F_x$')
    ax[3, 1].text(0.04, 0.96, 'h)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[3, 1].transAxes)
    ax[3, 1].legend(loc=5)

    ax[4, 0].plot(x[0], phiE, color='C2', label='$-\phi_E$')
    ax[4, 0].plot(x[0], phiFF1, color='C4', label='$-\phi_{ff}$')
    ax[4, 0].plot(x[0], phiE+phiFF1, color='C0', label='$-\phi_E-\phi_{ff}$')
    ax[4, 0].set_ylim((-0.10, 0.20))
    ax[4, 0].set_xlabel('$x\cdot k_0/2\pi$')
    ax[4, 0].set_ylabel('$\phi$')

    ax[4, 1].plot(x[0], Bz, color='C4', label='$B_z$')
    ax[4, 1].set_ylim((-0.2, 0.2))
    ax[4, 1].set_xlabel('$x\cdot k_0/2\pi$')
    ax[4, 1].set_ylabel('$B_z$')
    ax[4, 1].set_yticks((-0.1, 0.0, 0.1))
    ax[4, 1].text(0.04, 0.96, 'j)',
                  verticalalignment='top', horizontalalignment='left',
                  transform=ax[4, 1].transAxes)

    for i in range(10):
        ax.flatten()[i].set_autoscale_on(False)
        ax.flatten()[i].grid(True)
        if i != 1:
            ax.flatten()[i].set_xticks((0.25, 0.5, 0.75))
        ax.flatten()[i].legend(loc=5)
        #ax.flatten()[i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        ax.flatten()[i].autoscale(enable=True, axis='x', tight=True)

    fig.suptitle('t = {:.1f}'.format(tIdx*0.5))

    plt.tight_layout(pad=0.3)
    fig.subplots_adjust(top=0.94, bottom=0.06)

    return fig


plt.style.use('../../../paper.mplstyle')

for i in range(301):
    print(i)
    try:
        fig = plotFrame(i)
        plt.savefig('fig/f_{:03d}.png'.format(int(i)))
        plt.close(fig)
    except:
        print('fail')
   
