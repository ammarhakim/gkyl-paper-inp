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
T = linspace(0, 20, 200)
WX, WY = wR(T)

def plotFig(nr,nc,i,fr,showX=False,showY=False,showColorBar=False):
    # Inputs:
    # nr - number of subplot rows.
    # nc - number of subplot columns.
    # i - location in subplot.
    # fr - frame number to be read in.
    # showX, showY, showColorbar - booleans for axis/colorbar showing.
    # Output: subplot for figure.
    print("Working on %d ..." % i)
    # Note that this read-in assumes polynomial order 2 and default interpolation.
    data = postgkyl.GData("a2/a2-oscillating-E_elc_%d.bp" % fr)
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
    xticks(linspace(-6,6,5))
    yticks(linspace(-6,6,5))
    if showX:
        xlabel(r"$v_x (v_{th})$")
    if showY:
        ylabel(r"$v_y (v_{th})$")
    if showColorBar:
        colorbar(format='%.2f', ticks=np.linspace(0.0, 0.16, 5), fraction=0.046, pad=0.04)
        clim(0.0, 0.16)
    # Make aspect ratio equal.
    setp(f, aspect=1.0)
    title(r"$t = %d \Omega_{c}^{-1}$" % fr)

figure(1, figsize=(9,8))
plotFig(2, 2, 1, 4, False, True, False)
plotFig(2, 2, 2, 8, False, False, True)
plotFig(2, 2, 3, 12, True, True, False)
plotFig(2, 2, 4, 16, True, False, True)
tight_layout()
savefig("a2-oscillating-E-cmp.png", dpi=300)

show()    

