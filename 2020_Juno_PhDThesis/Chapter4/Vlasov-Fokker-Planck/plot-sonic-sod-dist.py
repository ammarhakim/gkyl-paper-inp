from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

def getXc(Xn):
    dx = Xn[1]-Xn[0]
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)
    return Xc

# Density and flow.
d = postgkyl.GData("n2/n2-sonic-sod-shock_neut_M0_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, m0 = dg.interpolate()

d = postgkyl.GData("n2/n2-sonic-sod-shock_neut_M1i_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, m1i = dg.interpolate()

Xn = X[0]; dx = Xn[1]-Xn[0]
# Cell-center configuration space coordinates.
Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)

# Distribution function.
d = postgkyl.GData("n2/n2-sonic-sod-shock_neut_1.bp")
dg = postgkyl.GInterpModal(d, 2, "ms")
X, fv1 = dg.interpolate()

# Center velocity space coordinates and construct meshgrid
Xn = X[0]; dx = Xn[1]-Xn[0]
Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)
Vn = X[1]; dv = Vn[1]-Vn[0]
Vc = linspace(Vn[0]+0.5*dv, Vn[-1]-0.5*dv, Vn.shape[0]-1)
XX, VV = meshgrid(Xc, Vc)

figure(1, figsize=(14,8))

subplot2grid((2,2),(0,0))
plot(Xc, m0)
xlim(-1,1)
ylabel("Density")

subplot2grid((2,2),(0,1))
plot(Xc, m1i/m0)
xlim(-1,1)
ylabel("Velocity")

subplot2grid((2,2),(1,0), colspan=2)
pcolormesh(XX, VV, transpose(fv1[:,:,0]))
title("Distribution Function",fontsize=20)
xlim(-1,1)
ylim(-6,6)
xticks(linspace(-1, 1, 5))
yticks(linspace(-6, 6, 5))
ylabel(r"$V (V_{th})$")
xlabel(r"$X$")

tight_layout()

savefig("sod-shock-distf.png", dpi=300)

show()

