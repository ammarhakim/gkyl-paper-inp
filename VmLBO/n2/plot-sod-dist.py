from pylab import *
import postgkyl
import numpy

style.use('../code/postgkyl.mplstyle')

def getXc(Xn):
    dx = Xn[1]-Xn[0]
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)
    return Xc

# density and velocity
d = postgkyl.GData("n2-sod-shock_neut_M0_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, m0 = dg.interpolate()
Xc = getXc(X[0])

d = postgkyl.GData("n2-sod-shock_neut_M1i_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, m1i = dg.interpolate()

# final solution
d = postgkyl.GData("n2-sod-shock_neut_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, fv1 = dg.interpolate()

XX, VV = meshgrid(X[0], X[1])

figure(1)

ax = subplot2grid((2,2),(0,0))
plot(Xc, m0)
gca().set_xlim([-1,1])
grid()
title('Density')

ax = subplot2grid((2,2),(0,1))
plot(Xc, m1i/m0)
gca().set_xlim([-1,1])
grid()
title('Velocity')

ax = subplot2grid((2,2),(1,0), colspan=2)
pcolormesh(XX, VV, transpose(fv1[:,:,0]))
title("Distribution Function")
ylabel('V')
xlabel('X')

tight_layout()

savefig('sod-shock-distf.png', dpi=150)

show()

