from pylab import *
import postgkyl

style.use("postgkyl.mplstyle")

def calcNormError(M):
    err = abs(M-M[0])/max(M)
    return err

def calcEnergyError(pre):
    d = postgkyl.GData("%s-sonic-sod-shock_neut_intM2Thermal.bp" % pre)
    m2Thermal = d.getValues()
    d = postgkyl.GData("%s-sonic-sod-shock_neut_intM2Flow.bp" % pre)
    m2Flow = d.getValues()

    t = d.getGrid()[0]
    eErr = calcNormError(m2Flow+m2Thermal)

    d = postgkyl.GData("%s-sonic-sod-shock_neut_intM1i.bp" % pre)
    m1i = d.getValues()
    mErr = calcNormError(m1i)

    return t, eErr, mErr

# p=1 case.
n1_t, n1_eErr, n1_mErr = calcEnergyError("n1/n1")
# p=2 case.
n2_t, n2_eErr, n2_mErr = calcEnergyError("n2/n2")

# Plot of energy error versus time.
figure(1)

semilogy(n1_t, n1_eErr, "r-", label=r"$p=1$")
semilogy(n2_t, n2_eErr, "k-", label=r"$p=2$")
legend(loc="lower right")
xlabel(r"Normalized time, $t$")
ylabel(r"$|\Delta M_2|/M_2(0)$")
xlim(0.0, 0.1)

savefig("sod-shock-er-conservation.png", dpi=300)

# Plot of momentum error versus time.
figure(2)
semilogy(n1_t, n1_mErr, "r-", label=r"$p=1$")
semilogy(n2_t, n2_mErr, "k-", label=r"$p=2$")
legend(loc="lower right")
xlabel(r"Normalized time, $t$")
ylabel(r"$|\Delta M_1|/M_1(0)$")
xlim(0.0, 0.1)

savefig("sod-shock-mom-conservation.png", dpi=300)

show()
