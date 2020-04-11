from pylab import *
import postgkyl

import numpy

style.use('../code/postgkyl.mplstyle')

def _colorbar(obj, fig, ax, label=""):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(obj, cax=cax, label=label)

def calcNormError(M):
    err = (M-M[0])/max(M)
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-bi-maxwellian-relax_neut_intM2Thermal_" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-bi-maxwellian-relax_neut_intM2Flow_" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(m2Flow+m2Thermal)

    return t, err

def calcMomentumError(pre):
    d = postgkyl.GData("%s-bi-maxwellian-relax_neut_intM1i_" % pre)
    m1i = d.getValues()

    t = d.getGrid()[0]
    err1 = calcNormError(m1i[:,0])
    err2 = calcNormError(m1i[:,1])

    return t, err1, err2

def getDist(pre, fr):
    d = postgkyl.GData("%s-bi-maxwellian-relax_neut_%d.bp" % (pre, fr))
    dg = postgkyl.GInterpModal(d, 2, "ms")
    XX, fv = dg.interpolate()

    X, V = meshgrid(XX[1], XX[2])
    return X, V, fv

def calcEntropy(X, V, fv):
    dx = X[1]-X[0]
    dv = V[1]-V[0]
    S = -dx*dv*sum(fv*log(abs(fv)))
    return S

def getEntropy(polyOrder, pre):
    svals = zeros((100,), float)
    for i in range(0,100):
        d = postgkyl.GData("%s-bi-maxwellian-relax_neut_%d.bp" % (pre, i))
        dg = postgkyl.GInterpModal(d, polyOrder, "ms")
        XX, fv = dg.interpolate()
        svals[i] = calcEntropy(XX[0], XX[1], fv)

    return svals

fig = figure(1)
X, V, f0 = getDist("r5", 0)
X, V, f100 = getDist("r5", 100)

ax = subplot(1,2,1)
im = pcolormesh(X, V, transpose(f0[0,:,:,0]))
axis ('image')
xlabel('$v_x$')
ylabel('$v_y$')
_colorbar(im, fig, ax)

ax = subplot(1,2,2)
im = pcolormesh(X, V, transpose(f100[0,:,:,0]))
ax.set_yticklabels("")
xlabel('$v_x$')
axis ('image')
_colorbar(im, fig, ax)

tight_layout()

savefig('bi-maxwellian-dist.png', dpi=150)

figure(2)

# plot of energy and momentum error v/s time
t, m2_err = calcEnergyError("r5")
t, m1i_0_err, m1i_1_err = calcMomentumError("r5")

plot(t, m2_err, 'b-', label='$\Delta M_2/M_2(0)$')
plot(t, m1i_0_err, 'm-', label='$\Delta M_{1,x}/M_{1,x}(0)$')
plot(t, m1i_1_err, 'r-', label='$\Delta M_{1,y}/M_{1,y}(0)$')
#legend(loc='best')
xlabel(r'$t\nu$')
grid()
gca().set_xlim([0, 5])

savefig('bi-maxwellian-er.png', dpi=150)

figure(3)
T = linspace(0, 5, 100)
s = getEntropy(2, "r5")
semilogx(T[1:], s[1:]/s[1]-1, '-b', label='$\Delta S/S(0)$')
xlabel(r'$t\nu$')
grid()
gca().set_xlim([T[1], 5])

savefig('bi-maxwellian-entropy.png', dpi=150)

show()
