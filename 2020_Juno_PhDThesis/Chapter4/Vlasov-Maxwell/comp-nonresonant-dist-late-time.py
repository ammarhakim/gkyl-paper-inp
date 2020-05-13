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
T = linspace(0, 100, 1000)
WX, WY = wNR(T)

def plotFig(nr,nc,i,fr):
    # Inputs:
    # nr - number of subplot rows.
    # nc - number of subplot columns.
    # i - location in subplot.
    # fr - frame number to be read in.
    # Output: subplot for figure.

    # Note that this read-in assumes default interpolation.
    data = postgkyl.GData("a3/a3-oscillating-E_elc_%d.bp" % fr)
    dg = postgkyl.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()

    # Center the grid values.
    for d in range(3):
        XX[d] = 0.5*(XX[d][:-1] + XX[d][1:])

    f = subplot(nr,nc,i)
    # Since solution is uniform in x, choose single x-point to visualize.
    # In this case we choose the x = Lx/2.
    pcolormesh(XX[1], XX[2], transpose(q[3,:,:,0]), shading="gouraud")
    plot(WX, WY, "w", linewidth=1.0)
    xlim(-6,6)
    ylim(-6,6)
    xlabel(r"$v_x (v_{th})$")
    ylabel(r"$v_y (v_{th})$")
    xticks(linspace(-6,6,5))
    yticks(linspace(-6,6,5))
    clim(0.0, 0.16)
    colorbar(format='%.2f', ticks=np.linspace(0.0, 0.16, 5), fraction=0.046, pad=0.04)
    # Make aspect ratio equal.
    setp(f, aspect=1.0)
    title("Polynomial Order 2")

    # Note that this read-in assumes default interpolation.
    data = postgkyl.GData("a4/a4-oscillating-E_elc_%d.bp" % fr)
    dg = postgkyl.data.GInterpModal(data, 3, "ms")
    XX, q = dg.interpolate()

    # Center the grid values.
    for d in range(3):
        XX[d] = 0.5*(XX[d][:-1] + XX[d][1:])

    f = subplot(nr,nc,i+1)
    # Since solution is uniform in x, choose single x-point to visualize.
    # In this case we choose the x = Lx/2.
    pcolormesh(XX[1], XX[2], transpose(q[3,:,:,0]), shading="gouraud")
    plot(WX, WY, "w", linewidth=1.0)
    xlim(-6,6)
    ylim(-6,6)
    xticks(linspace(-6,6,5))
    yticks(linspace(-6,6,5))
    xlabel(r"$v_x (v_{th})$")
    clim(0.0, 0.16)
    colorbar(format='%.2f', ticks=np.linspace(0.0, 0.16, 5), fraction=0.046, pad=0.04)
    # Make aspect ratio equal.
    setp(f, aspect=1.0)
    title("Polynomial Order 3")

figure(1, figsize=(14,7))
plotFig(1, 2, 1, 1)
suptitle(r"$t = 1000 \Omega_c^{-1}$")
tight_layout()
savefig("comp-long-time-oscc-E-cmp.png", dpi=300)

show()

