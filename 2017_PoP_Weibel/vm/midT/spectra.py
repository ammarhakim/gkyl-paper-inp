#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg
from scipy import fftpack
from glob import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable

def getSpectraAndDisp(k0, nK, nW, root, comp=0):
    files = glob(root)
    numFiles = len(files)

    d = pg.GData(files[0])
    dg = pg.GInterpNodal(d, 2, 'ns')
    x, val = dg.interpolate(0)
    evolution = np.zeros((numFiles, len(x[0])))
    spectrum = np.zeros((numFiles, nK))
    t = np.zeros(numFiles)
    for i, f in enumerate(files):
        d = pg.GData(f)
        t[i] = d.time 
        dg = pg.GInterpNodal(d, 2, 'ns')
        x, val = dg.interpolate(comp)
    
        evolution[i, :] = val

        valf = fftpack.fft(val)
        spectrum[i, :] = np.abs(valf[1:(nK+1)])**2
    k = fftpack.fftfreq(len(val), x[0][1]-x[0][0])
    k = 2*np.pi*k[1:(nK+1)]/k0
    x = x[0]*k0/2/np.pi
        
    disp = np.zeros((nW, nK))
    for i in range(nK):
        ft = fftpack.fft(spectrum[:, i])
        disp[:, i] = np.abs(ft[1:(nW+1)])**2
    omega = fftpack.fftfreq(numFiles, t[1]-t[0])
    omega = 2*np.pi*omega[1:(nW+1)]

    return x, k, t, omega, evolution, spectrum, disp

def cb(fig, ax, im):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='2%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

plt.style.use('../../../paper.mplstyle')

#---------------------------------------------------------------------
#-- Number density --------------------------------------------------- 
fig, ax = plt.subplots(1,3, figsize=(16,4))
x, k, t, w, ev, sp, dp = getSpectraAndDisp(0.4, 15, 75,
                                           'ck_1X2V_numDensityElc_*.h5')
im0 = ax[0].pcolormesh(t, x, ev.transpose())
cb(fig, ax[0], im0)
im1 = ax[1].pcolormesh(t, k-0.5, sp.transpose())
cb(fig, ax[1], im1)
im2 = ax[2].pcolormesh(w, k-0.5, dp.transpose())
cb(fig, ax[2], im2)

ax[0].set_ylabel(r'$x\cdot k_0/2\pi$')
ax[0].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[0].set_yticks((0.25, 0.5, 0.75))
ax[0].grid()
ax[1].set_ylabel(r'$k/k_0$')
ax[1].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[1].set_yticks(range(1,15))
ax[1].grid()
ax[2].set_ylabel(r'$k/k_0$')
ax[2].set_xlabel(r'$\omega/\omega_{pe}$')
ax[2].set_yticks(range(1,15))
ax[2].grid()
fig.suptitle('Number density')
plt.tight_layout()
plt.savefig('fig/fft_n.png', dpi=150)

#---------------------------------------------------------------------
#-- Bz --------------------------------------------------------------- 
fig, ax = plt.subplots(1,3, figsize=(16,4))
x, k, t, w, ev, sp, dp = getSpectraAndDisp(0.4, 15, 75,
                                           'ck_1X2V_em_*.h5', 5)
im0 = ax[0].pcolormesh(t, x, ev.transpose())
cb(fig, ax[0], im0)
im1 = ax[1].pcolormesh(t, k-0.5, sp.transpose())
cb(fig, ax[1], im1)
im2 = ax[2].pcolormesh(w, k-0.5, dp.transpose())
cb(fig, ax[2], im2)

ax[0].set_ylabel(r'$x\cdot k_0/2\pi$')
ax[0].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[0].set_yticks((0.25, 0.5, 0.75))
ax[0].grid()
ax[1].set_ylabel(r'$k/k_0$')
ax[1].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[1].set_yticks(range(1,15))
ax[1].grid()
ax[2].set_ylabel(r'$k/k_0$')
ax[2].set_xlabel(r'$\omega/\omega_{pe}$')
ax[2].set_yticks(range(1,15))
ax[2].grid()
fig.suptitle('$B_z$')
plt.tight_layout()
plt.savefig('fig/fft_Bz.png', dpi=150)

#---------------------------------------------------------------------
#-- Ex---------------------------------------------------------------- 
fig, ax = plt.subplots(1,3, figsize=(16,4))
x, k, t, w, ev, sp, dp = getSpectraAndDisp(0.4, 15, 75,
                                           'ck_1X2V_em_*.h5')
im0 = ax[0].pcolormesh(t, x, ev.transpose())
cb(fig, ax[0], im0)
im1 = ax[1].pcolormesh(t, k-0.5, sp.transpose())
cb(fig, ax[1], im1)
im2 = ax[2].pcolormesh(w, k-0.5, dp.transpose())
cb(fig, ax[2], im2)

ax[0].set_ylabel(r'$x\cdot k_0/2\pi$')
ax[0].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[0].set_yticks((0.25, 0.5, 0.75))
ax[0].grid()
ax[1].set_ylabel(r'$k/k_0$')
ax[1].set_xlabel(r'$t\cdot\omega_{pe}$')
ax[1].set_yticks(range(1,15))
ax[1].grid()
ax[2].set_ylabel(r'$k/k_0$')
ax[2].set_xlabel(r'$\omega/\omega_{pe}$')
ax[2].set_yticks(range(1,15))
ax[2].grid()
fig.suptitle('$E_x$')
plt.tight_layout()
plt.savefig('fig/fft_Ex.png', dpi=150)


plt.show()
