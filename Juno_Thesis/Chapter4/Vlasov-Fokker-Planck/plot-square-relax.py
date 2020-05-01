from pylab import *
import postgkyl

import numpy

style.use('postgkyl.mplstyle')

def calcNormError(M):
    err = abs((M-M[0])/M[0])
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-relax_neut_intM2Thermal.bp" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-relax_neut_intM2Flow.bp" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(m2Flow+m2Thermal)

    return t, err

# p=1 case
r1_t, r1_err = calcEnergyError("r1/r1")
# p=2 case
r2_t, r2_err = calcEnergyError("r2/r2")

# plot of energy error v/s time
figure(1)
semilogy(r1_t, r1_err, 'r-', label='$p=1$')
semilogy(r2_t, r2_err, 'k-', label='$p=2$')
legend(loc='best')
xlabel(r'$t\nu$')
ylabel(r'$|\Delta M_2|/M_2(0)$')
xlim(0,5)
savefig('square-relax-er.png', dpi=300)

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

r1_s = getEntropy(1, "r1/r1")
semilogx(T[1:], r1_s[1:]/r1_s[1]-1, 'r-', label='$p=1$')

r2_s = getEntropy(2, "r2/r2")
semilogx(T[1:], r2_s[1:]/r2_s[1]-1, 'k-', label='$p=2$')

legend(loc='best')
xlabel(r'$t\nu$')
ylabel(r'$\Delta S/S(0)$')
xlim(T[1], 5)
savefig('square-relax-entropy.png', dpi=300)

show()
