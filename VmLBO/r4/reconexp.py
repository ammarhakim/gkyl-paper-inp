from numpy import *
from pylab import *

import postgkyl
import scipy.optimize

style.use('../code/postgkyl.mplstyle')

def func(an, f0, f1):
    g0 = an[0]; g1 = an[1]
    rhs0 = ((exp(g1+g0))/(g1)-(exp(g0-g1))/(g1))/(sqrt(2.0))
    rhs1 = (sqrt(3)*(((exp(g0)*g1-exp(g0))*exp(g1))/(g1**2)+((exp(g0)*g1+exp(g0))*exp(-g1))/(g1**2)))/(sqrt(2))

    return rhs0-f0, rhs1-f1

def fitExp(f0, f1):
    aout = scipy.optimize.fsolve(func, [1.0, 0.01], args=(f0, f1))
    return aout[0], aout[1]

def getExpRecon(pre, fr):
    d = postgkyl.GData("%s-relax_neut_%d.bp" % (pre, fr))
    q = d.getValues()
    f0 = q[0,:,0]
    f1 = q[0,:,2]

    gcoeff = zeros((f0.shape[0],2), float)
    for i in range(f0.shape[0]):
        gcoeff[i][0], gcoeff[i][1] = fitExp(f0[i], f1[i])

    return gcoeff

def evalLin(f0, f1, X):
    return f0/sqrt(2.0) + sqrt(3.0/2.0)*f1*X

def evalExp(g0, g1, X):
    return exp(g0+g1*X)

vlo = -6.0
vup = 6.0

d = postgkyl.GData("r4-relax_neut_100.bp")
q = d.getValues()
f0 = q[0,:,0]
f1 = q[0,:,2]

X = linspace(-1, 1, 20)  # for interpolation
dx = 12/f0.shape[0]

figure(1)
for i in range(f0.shape[0]):
    v = linspace(-6.0+i*dx, -6.0+(i+1)*dx, X.shape[0])
    plot(v, evalLin(f0[i], f1[i], X), 'r-')
    if f0[i] < 0:
        plot(0.5*(v[0]+v[-1]), f0[0], 'mo')

for i in range(f0.shape[0]):
    if f0[i] > 0:
        g0, g1 = fitExp(f0[i], f1[i])
        v = linspace(-6.0+i*dx, -6.0+(i+1)*dx, X.shape[0])
        plot(v, evalExp(g0, g1, X), 'k-')

grid()    
show()

