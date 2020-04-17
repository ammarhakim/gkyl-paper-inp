from pylab import *
import postgkyl

import numpy

style.use('../code/postgkyl.mplstyle')

def calcGrad(x, ff):
    f = ff
    if len(ff.shape) > 1:
        f = ff[:,0]
        
    dx = x[1]-x[0]
    N = f.shape[0]
    df = zeros((f.shape[0]+1,), float)
    df[0] = (f[1]-f[0])/dx
    for i in range(1,N):
        df[i] = (f[i]-f[i-1])/dx
    df[N] = (f[N-1]-f[N-2])/dx

    return df

def calcMidVals(f):
    return 0.5*(f[0:-1]+f[1:])

def getXc(Xn):
    dx = Xn[1]-Xn[0]
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)
    return Xc

def getMoments(polyOrder, pre, numInterp):
    d = postgkyl.GData("%s-sod-shock_neut_M0_1.bp" % pre)
    dg = postgkyl.GInterpModal(d, polyOrder, "ms", numInterp=numInterp)
    X, m0 = dg.interpolate()

    d = postgkyl.GData("%s-sod-shock_neut_M1i_1.bp" % pre)
    dg = postgkyl.GInterpModal(d, polyOrder, "ms", numInterp=numInterp)
    X, m1i = dg.interpolate()

    d = postgkyl.GData("%s-sod-shock_neut_M2_1.bp" % pre)
    dg = postgkyl.GInterpModal(d, polyOrder, "ms", numInterp=numInterp)
    X, m2 = dg.interpolate()

    d = postgkyl.GData("%s-sod-shock_neut_M3i_1.bp" % pre)
    dg = postgkyl.GInterpModal(d, polyOrder, "ms", numInterp=numInterp)
    X, m3i = dg.interpolate()

    u = m1i/m0 # velocity
    nvt2 = m2 - m0*u**2 # ptcl internal energy
    q = m3i - (3*u*m2 - 3*u**2*m1i + u**3*m0) # heat-flux

    Xn = X[0]; dx = Xn[1]-Xn[0]
    # cell-center coordinates
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)

    return Xc, m0, u, nvt2, q

X, s1_n0, s1_u, s1_ie, s1_q = getMoments(2, "s1", 3)
X, s2_n0, s2_u, s2_ie, s2_q = getMoments(2, "../s2/s2", 3)
X, s3_n0, s3_u, s3_ie, s3_q = getMoments(2, "../s3/s3", 3)

# exact Euler
eu_n0 = loadtxt("../s4/s4-sod-shock-exact-density.txt")
eu_u = loadtxt("../s4/s4-sod-shock-exact-velocity.txt")
eu_ie = loadtxt("../s4/s4-sod-shock-exact-internal-energy.txt")

figure(1)
ax = subplot(2,2,1)
plot(X, s1_n0, '-r')
plot(X, s2_n0, '-m')
plot(X, s3_n0, '-b')
plot(eu_n0[:,0], eu_n0[:,1], 'k--')
grid()
title('Density')
gca().set_xlim([0,1])
ax.set_xticklabels([""])

ax = subplot(2,2,2)
plot(X, s1_u, '-r')
plot(X, s2_u, '-m')
plot(X, s3_u, '-b')
plot(eu_n0[:,0], eu_u[:,1], 'k--')
grid()
title('Velocity')
gca().set_xlim([0,1])
ax.set_xticklabels([""])

subplot(2,2,3)
# plot(X, s1_ie/2, '-r')
# plot(X, s2_ie/2, '-m')
# plot(X, s3_ie/2, '-b')
# plot(eu_n0[:,0], eu_n0[:,1]*eu_ie[:,1], 'k--')
# title('Internal Energy')
# xlabel('X')
# grid()
# gca().set_xlim([0,1])

plot(X, s1_ie/s1_n0, '-r')
plot(X, s2_ie/s2_n0, '-m')
plot(X, s3_ie/s3_n0, '-b')
plot(eu_n0[:,0], 2*eu_ie[:,1], 'k--')
title('Temperature')
xlabel('X')
grid()
gca().set_xlim([0,1])

subplot(2,2,4)
plot(X, s1_q, '-r')
plot(X, s2_q, '-m')
plot(X, s3_q, '-b')
plot(eu_n0[:,0], 0.0*eu_n0[:,0], 'k--') # no heat-flux in Euler
grid()
title('Heat flux')
xlabel('X')
gca().set_xlim([0,1])

#tight_layout()

savefig('sod-shock-moments.png', dpi=150)

show()



    
