from pylab import *
import postgkyl

import numpy

style.use('../code/postgkyl.mplstyle')

def calcNormError(M):
    err = (M-M[0])/max(M)
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-sod-shock_neut_intM2Thermal_" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-sod-shock_neut_intM2Flow_" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    eErr = calcNormError(m2Flow+m2Thermal)

    d = postgkyl.GData("%s-sod-shock_neut_intM1i_" % pre)
    m1i = d.getValues()
    mErr = calcNormError(m1i)

    return t, eErr, mErr

# p=2 case
n2_t, n2_eErr, n2_mErr = calcEnergyError("n2")
# p=1 case
n3_t, n3_eErr, n3_mErr = calcEnergyError("../n3/n3")

r = 1.0*n2_t.shape[0]/n3_t.shape[0]

# plot of energy error v/s time
figure(1)

plot(n3_t, n3_eErr, 'r-', label='$p=1$')
plot(n2_t, n2_eErr/r, 'k-', label='$p=2$')
legend(loc='best')
xlabel(r'T')
title(r'$\Delta M_2/M_2(0)$')
grid()
gca().set_xlim([0, 0.1])

savefig('sod-shock-er-conservation.png', dpi=150)

# plot of momentum error v/s time
figure(2)
plot(n3_t, n3_mErr, 'r-', label='$p=1$')
plot(n2_t, n2_mErr/r, 'k-', label='$p=2$')
legend(loc='best')
xlabel(r'T')
title(r'$\Delta M_1/M_1(0)$')
grid()
gca().set_xlim([0, 0.1])

savefig('sod-shock-mom-conservation.png', dpi=150)

show()
