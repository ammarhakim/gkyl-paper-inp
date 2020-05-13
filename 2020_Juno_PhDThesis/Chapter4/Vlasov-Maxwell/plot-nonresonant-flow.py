from pylab import *
import postgkyl
style.use("postgkyl.mplstyle")

# Normalizations from code.
E0 = 1.0
w = 0.5

# Non-resonant analytic solution. 
def wNR(t):
    wx = E0/(1-w**2)*(sin(t)-w*sin(w*t))
    wy = E0/(1-w**2)*(cos(t)-cos(w*t))
    return wx, wy

# Compute analytic solution in time for comparison.
Tex = linspace(0, 100, 5000)
wx, wy = wNR(Tex)

def plotFig():
    ux = zeros((101,), float)
    uy = zeros((101,), float)

    wx = zeros((101,), float)
    wy = zeros((101,), float)
    
    for i in range(101):
        print("Working on %d ..." % i)
        # Note that this read-in assumes polynomial order 2 and default interpolation.
        data = postgkyl.GData("a1/a1-oscillating-E_elc_M1i_%d.bp" % i)
        dg = postgkyl.data.GInterpModal(data, 2, "ms")
        XX, Ux = dg.interpolate(0)
        XX, Uy = dg.interpolate(1)

        ux[i] = Ux[0]
        uy[i] = Uy[0]        
        
        tm = data.time
        wx1, wy1 = wNR(tm)
        wx[i] = wx1
        wy[i] = wy1

    return ux, uy

figure(1)
T = linspace(0, 100, 101)
ux, uy = plotFig()
subplot(2,1,1)
plot(T, ux, "ro", Tex, -wx, "k-")
ylabel(r"$u_x$")
xlim(0, 100)
ylim(-3, 3)

subplot(2,1,2)
plot(T, uy, "ro", Tex, wy, "k-")
xlabel(r"Time $(\Omega_c^{-1})$")
ylabel(r"$u_y$")
xlim(0, 100)
ylim(-3, 3)

savefig("a1-oscillating-E-c-cmp.png", dpi=300)

show()    

