from pylab import *
import postgkyl
style.use("postgkyl.mplstyle")

# Normalizations from code.
E0 = 0.5
w = 1.0

# Resonant analytic solution. 
def wR(t):
    wx = E0/2*(t*cos(t)+sin(t))
    wy = -E0/2*t*sin(t)
    return wx, wy

# Compute analytic solution in time for comparison.
Tex = linspace(0, 20, 500)
wx, wy = wR(Tex)

def plotFig():
    ux = zeros((21,), float)
    uy = zeros((21,), float)

    wx = zeros((21,), float)
    wy = zeros((21,), float)
    
    for i in range(21):
        print("Working on %d ..." % i)
        # Note that this read-in assumes polynomial order 2 and default interpolation.
        data = postgkyl.GData("a2/a2-oscillating-E_elc_M1i_%d.bp" % i)
        dg = postgkyl.data.GInterpModal(data, 2, "ms")
        XX, Ux = dg.interpolate(0)
        XX, Uy = dg.interpolate(1)

        ux[i] = Ux[0]
        uy[i] = Uy[0]        
        
        tm = data.time
        wx1, wy1 = wR(tm)
        wx[i] = wx1
        wy[i] = wy1

    return ux, uy

figure(1)
T = linspace(0, 20, 21)
ux, uy = plotFig()
subplot(2,1,1)
plot(T, ux, "ro", Tex, -wx, "k-")
ylabel(r"$u_x$")
xlim(0, 20)
ylim(-5, 5)

subplot(2,1,2)
plot(T, uy, "ro", Tex, wy, "k-")
xlabel(r"Time $(\Omega_c^{-1})$")
ylabel(r"$u_y$")
xlim(0, 20)
ylim(-5, 5)

savefig("a2-oscillating-E-c-cmp.png", dpi=300)

show()    

