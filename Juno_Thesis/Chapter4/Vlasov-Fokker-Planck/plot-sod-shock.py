from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

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

    u = m1i/m0 # Velocity.
    nvt2 = m2 - m0*u**2 # Pressure (density*temperature).
    q = m3i - (3*u*m2 - 3*u**2*m1i + u**3*m0) # Plasma-frame heat-flux.

    Xn = X[0]; dx = Xn[1]-Xn[0]
    # Cell-center coordinates.
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)

    return Xc, m0, u, nvt2, q

X, s1_n0, s1_u, s1_ie, s1_q = getMoments(2, "s1/s1", 3)
X, s2_n0, s2_u, s2_ie, s2_q = getMoments(2, "s2/s2", 3)
X, s3_n0, s3_u, s3_ie, s3_q = getMoments(2, "s3/s3", 3)

# Exact Euler solution.
eu_n0 = loadtxt("s4/s4-sod-shock-exact-density.txt")
eu_u = loadtxt("s4/s4-sod-shock-exact-velocity.txt")
eu_ie = loadtxt("s4/s4-sod-shock-exact-internal-energy.txt")

figure(1, figsize=(14,8))
subplot(2,2,1)
plot(X, s1_n0, "-r", label="Kn = 1/10")
plot(X, s2_n0, "-m", label="Kn = 1/100")
plot(X, s3_n0, "-b", label="Kn = 1/500")
plot(eu_n0[:,0], eu_n0[:,1], "k--", label="Inviscid Euler")
ylabel("Density")
xlim(0,1)
legend(loc="best",fontsize=16)

subplot(2,2,2)
plot(X, s1_u, "-r")
plot(X, s2_u, "-m")
plot(X, s3_u, "-b")
plot(eu_n0[:,0], eu_u[:,1], "k--")
ylabel("Velocity")
xlim(0,1)

subplot(2,2,3)
plot(X, s1_ie/s1_n0, "-r")
plot(X, s2_ie/s2_n0, "-m")
plot(X, s3_ie/s3_n0, "-b")
plot(eu_n0[:,0], 2*eu_ie[:,1], "k--")
ylabel("Temperature")
xlabel("X")
xlim(0,1)

subplot(2,2,4)
plot(X, s1_q, "-r")
plot(X, s2_q, "-m")
plot(X, s3_q, "-b")
plot(eu_n0[:,0], 0.0*eu_n0[:,0], "k--") # No heat-flux in Euler solution.
ylabel("Heat flux")
xlabel("X")
xlim(0,1)

tight_layout()

savefig("sod-shock-moments.png", dpi=300)

show()



    
