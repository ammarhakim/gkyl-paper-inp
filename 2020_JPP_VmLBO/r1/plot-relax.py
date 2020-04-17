from pylab import *
import postgkyl

import numpy

style.use('../code/postgkyl.mplstyle')


def calcNormError(M):
    err = (M-M[0])/max(M)
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-relax_neut_intM2Thermal_" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-relax_neut_intM2Flow_" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(m2Flow+m2Thermal)

    return t, err

# p=1 case
p4_t, p4_err = calcEnergyError("../r4/r4")
# p=2 case
p2_t, p2_err = calcEnergyError("../r2/r2")

r = 1.0*p2_t.shape[0]/p4_t.shape[0]

# plot of energy error v/s time
figure(1)
plot(p4_t, p4_err, 'r-', label='$p=1$')
plot(p2_t, p2_err/r, 'k-', label='$p=2$')
legend(loc='best')
xlabel(r'$t\nu$')
title(r'$\Delta M_2/M_2(0)$')
grid()
gca().set_xlim([0, 5])
savefig('square-relax-er.png', dpi=150)

def calcEntropy(X, V, fv):
    dx = X[1]-X[0]
    dv = V[1]-V[0]
    S = -dx*dv*sum(fv*log(abs(fv)))
    return S

def getEntropy(polyOrder, pre):
    svals = zeros((100,), float)
    for i in range(0,100):
        d = postgkyl.GData("%s-relax_neut_%d.bp" % (pre, i))
        dg = postgkyl.GInterpModal(d, polyOrder, "ms")
        XX, fv = dg.interpolate()
        svals[i] = calcEntropy(XX[0], XX[1], fv)

    return svals

# plot of entropy v/s time
figure(2)
T = linspace(0, 5, 100)
#p1_s = getEntropy(1, "r1")
#semilogx(T[1:], p1_s[1:]/p1_s[1], 'r-', label='$p=1$')

p4_s = getEntropy(1, "../r4/r4")
semilogx(T[1:], p4_s[1:]/p4_s[1]-1, 'r-', label='$p=1$')

p2_s = getEntropy(2, "../r2/r2")
semilogx(T[1:], p2_s[1:]/p2_s[1]-1, 'k-', label='$p=2$')

legend(loc='best')
xlabel(r'$t\nu$')
title(r'$\Delta S/S(0)$')
gca().set_xlim([T[1], 5])
grid()
savefig('square-relax-entropy.png', dpi=150)

show()
