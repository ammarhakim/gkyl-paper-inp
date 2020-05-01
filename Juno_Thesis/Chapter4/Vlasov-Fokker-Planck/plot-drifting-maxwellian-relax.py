from pylab import *
import postgkyl

style.use('postgkyl.mplstyle')

def calcNormError(M):
    err = abs((M-M[0])/M[0])
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-drifting-maxwellian-relax_neut_intM2Thermal.bp" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-drifting-maxwellian-relax_neut_intM2Flow.bp" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    err = calcNormError(m2Flow+m2Thermal)

    return t, err

def calcMomentumError(pre):
    d = postgkyl.GData("%s-drifting-maxwellian-relax_neut_intM1i.bp" % pre)
    m1i = d.getValues()

    t = d.getGrid()[0]
    err1 = calcNormError(m1i[:,0])
    err2 = calcNormError(m1i[:,1])

    return t, err1, err2

def getDist(pre, fr):
    d = postgkyl.GData("%s-drifting-maxwellian-relax_neut_%d.bp" % (pre, fr))
    dg = postgkyl.GInterpModal(d, 2, "ms")
    XX, fv = dg.interpolate()

    X, V = meshgrid(XX[1], XX[2])
    return X, V, fv

def calcEntropy(X, V, fv):
    dx = X[1]-X[0]
    dv = V[1]-V[0]
    S = -dx*dv*sum(fv*log(abs(fv)))
    return S

def getEntropy(polyOrder, pre):
    svals = zeros((100,), float)
    for i in range(0,100):
        d = postgkyl.GData("%s-drifting-maxwellian-relax_neut_%d.bp" % (pre, i))
        dg = postgkyl.GInterpModal(d, polyOrder, "ms")
        XX, fv = dg.interpolate()
        svals[i] = calcEntropy(XX[0], XX[1], fv)

    return svals

X, V, f0 = getDist("r3/r3", 0)
X, V, f100 = getDist("r3/r3", 100)
figure(1)
pcolormesh(X, V, transpose(f0[0,:,:,0]))
setp(gca(), aspect=1.0)
xlabel('$V_X (V_{th})$')
ylabel('$V_Y (V_{th})$')
colorbar(format='%.2f', fraction=0.046, pad=0.04)
savefig('drifting-maxwellian-dist-init.png', dpi=300)

figure(2)
pcolormesh(X, V, transpose(f100[0,:,:,0]))
setp(gca(), aspect=1.0)
xlabel('$V_X (V_{th})$')
ylabel('$V_Y (V_{th})$')
colorbar(format='%.2f', fraction=0.046, pad=0.04)
savefig('drifting-maxwellian-dist-final.png', dpi=300)

figure(3)
# plot of energy and momentum error v/s time
t, m2_err = calcEnergyError("r3/r3")
t, m1i_0_err, m1i_1_err = calcMomentumError("r3/r3")

semilogy(t, m2_err, 'b-', label='$\Delta M_2/M_2(0)$')
semilogy(t, m1i_0_err, 'm-', label='$\Delta M_{1,x}/M_{1,x}(0)$')
semilogy(t, m1i_1_err, 'r-', label='$\Delta M_{1,y}/M_{1,y}(0)$')
legend(loc='best')
xlabel(r'$t\nu$')
ylabel(r'$|\Delta M|/M(0)$')
xlim(0,5)

savefig('drifting-maxwellian-er.png', dpi=300)

figure(4)
T = linspace(0, 5, 100)
s = getEntropy(2, "r3/r3")
semilogx(T[1:], s[1:]/s[1]-1, '-b', label='$\Delta S/S(0)$')
xlabel(r'$t\nu$')
ylabel(r'$\Delta S/S(0)$')
xlim([T[1], 5])
ylim(0.0, 0.4)

savefig('drifting-maxwellian-entropy.png', dpi=300)

show()
