from numpy import *
import scipy.optimize

def func(an, f0, f1):
    a0 = an[0]; a1 = an[1]
    rhs0 = ((exp(g1+g0))/(g1)-(exp(g0-g1))/(g1))/(sqrt(2.0))
    rhs1 = (sqrt(3)*(((exp(g0)*g1-exp(g0))*exp(g1))/(g1^2)+((exp(g0)*g1+exp(g0))*exp(-g1))/(g1^2)))/(sqrt(2))

    return rhs0-f0, rhs1-f1

# compute g0 and g1 for f0=1, f1=1.0 with initial guess 1.0, 0.01
aout = scipy.optimize.fsolve(func, [1.0, 0.01], args=(1.0, 1.0))
