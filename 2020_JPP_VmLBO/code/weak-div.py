from pylab import *

style.use('postgkyl.mplstyle')

def calcN(X, n0, n1):
    return (n0/sqrt(2) + sqrt(3.0/2.0)*n1*X)

def calcU(X, n0, n1):
    return sqrt(2)/(n0**2-n1**2)*(n0-sqrt(3)*n1*X)

def calcUex(X, n0, n1):
    return 1/(n0/sqrt(2) + sqrt(3.0/2.0)*n1*X)

X = linspace(-1, 1, 10)

# for computing ranges
n1Max = 0.85
nMax = calcN(X, 1.0, n1Max)
uMax = calcU(X, 1.0, n1Max)
nMax_lims = [nMax.min(), nMax.max()]
uMax_lims = [uMax.min(), uMax.max()]

fig = figure(1)

# Case 1
ax = subplot2grid((2,3), (0,0))
ax.plot(X, calcN(X, 1.0, 0.25), 'k-')
ax.set_ylabel(r'$M_0$')
ax.set_xticklabels("")
gca().set_xlim([-1,1])
gca().set_ylim(nMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], nMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

ax = subplot2grid((2,3), (1,0))
ax.plot(X, calcU(X, 1.0, 0.25), 'r')
ax.set_ylabel('U')
ax.set_xlabel('X')
gca().set_xlim([-1,1])
gca().set_ylim(uMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], uMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

# Case 2
ax = subplot2grid((2,3), (0,1))
ax.plot(X, calcN(X, 1.0, 1/sqrt(3.0)), 'k-')
ax.set_xticklabels("")
ax.set_yticklabels("")
gca().set_xlim([-1,1])
gca().set_ylim(nMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], nMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

ax = subplot2grid((2,3), (1,1))
ax.plot(X, calcU(X, 1.0, 1/sqrt(3.0)), 'r')
ax.set_yticklabels("")
ax.set_xlabel('X')
gca().set_xlim([-1,1])
gca().set_ylim(uMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], uMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

# Case 3
ax = subplot2grid((2,3), (0,2))
ax.plot(X, calcN(X, 1.0, n1Max), 'k-')
ax.set_xticklabels("")
ax.set_yticklabels("")
gca().set_xlim([-1,1])
gca().set_ylim(nMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], nMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

ax = subplot2grid((2,3), (1,2))
ax.plot(X, calcU(X, 1.0, n1Max), 'r')
ax.set_yticklabels("")
ax.set_xlabel('X')
gca().set_xlim([-1,1])
gca().set_ylim(uMax_lims)

plot([-1,1], [0,0], 'm-', linewidth=0.5)
plot([0,0], uMax_lims, 'm-', linewidth=0.5)
plot([1/sqrt(3.0)], [0.0], 'mx')
plot([-1/sqrt(3.0)], [0.0], 'mx')

#grid()

#tight_layout()
savefig('weak-divide-p1.png', dpi=150)
show()

