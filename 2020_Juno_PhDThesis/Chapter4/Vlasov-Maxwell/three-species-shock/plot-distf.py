from pylab import *
import postgkyl as pg
import numpy as np
style.use("../postgkyl.mplstyle")

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i, figsize=(14,8))
    dataAl = pg.GData("vlasov-riemann-3-species-1x1v_Al_%d.bp" % fr)
    dgAl = pg.data.GInterpModal(dataAl, 2, "ms")
    XX_Al, q_Al = dgAl.interpolate()
    lambdaD = 0.1
    vthAl = 0.00034788351
    #center the grid values
    for d in range(2):
        XX_Al[d] = 0.5*(XX_Al[d][:-1] + XX_Al[d][1:])

    q_Al_rescaled = np.zeros((q_Al.shape[0], q_Al.shape[1]))
    max_q_Al = np.max(q_Al)
    for j in range(0, q_Al.shape[0]):
        for k in range(0, q_Al.shape[1]):
            if q_Al[j,k,0] > 0.0:
                q_Al_rescaled[j,k] = sqrt(q_Al[j,k]/max_q_Al)
                
    dataIon = pg.GData("vlasov-riemann-3-species-1x1v_ion_%d.bp" % fr)
    dgIon = pg.data.GInterpModal(dataIon, 2, "ms")
    XX_ion, q_ion = dgIon.interpolate()
    vthIon = 0.00034788344
    #center the grid values
    for d in range(2):
        XX_ion[d] = 0.5*(XX_ion[d][:-1] + XX_ion[d][1:])

    q_ion_rescaled = np.zeros((q_ion.shape[0], q_ion.shape[1]))
    max_q_ion = np.max(q_ion)
    for j in range(0, q_ion.shape[0]):
        for k in range(0, q_ion.shape[1]):
            if q_ion[j,k,0] > 0.0:
                q_ion_rescaled[j,k] = sqrt(q_ion[j,k]/max_q_ion)

    subplot(1, 2, 1)
    pcolormesh(XX_Al[0]/lambdaD, XX_Al[1]/vthAl-5.66, transpose(q_Al_rescaled), shading="gouraud")
    colorbar()
    axis("tight")
    xlim(60, 90)
    ylim(-8, 0)
    xlabel("$X (\lambda_D)$")
    ylabel("$V_x (v_{th_{p}})$")
            
    subplot(1, 2, 2)
    pcolormesh(XX_ion[0]/lambdaD, XX_ion[1]/vthIon-5.66, transpose(q_ion_rescaled), shading="gouraud")
    colorbar()
    axis("tight")
    xlim(60, 90)
    ylim(-8, 8)
    xlabel("$X (\lambda_D)$")
    ylabel("$V_x (v_{th_{p}})$")
    
    tight_layout()
    savefig("vlasov-riemann-3-species-1x1v-distf-%05d.png" % i)
    
    #Find reflected density
    shock_frame_vt = XX_ion[1]/vthIon-5.435
    dv_shock_frame = shock_frame_vt[1] - shock_frame_vt[0]
    n_incoming = np.zeros(q_ion.shape[0])
    n_reflected = np.zeros(q_ion.shape[0])
    n_incoming[:] = np.sum(q_ion[:,0:136,0], axis=1)*dv_shock_frame
    n_reflected[:] = np.sum(q_ion[:,136:,0], axis=1)*dv_shock_frame
    figure(i+1)
    plot(XX_ion[0]/lambdaD, n_reflected/n_incoming)

    dataphi = pg.GData("vlasov-riemann-3-species-1x1v_field_%d.bp" % fr)
    dgphi = pg.data.GInterpModal(dataphi, 2, "ms")
    XX_phi, Ex = dgphi.interpolate()
    #center the grid values
    for d in range(1):
        XX_phi[d] = 0.5*(XX_phi[d][:-1] + XX_phi[d][1:])
        
    dx = XX_phi[0][1] - XX_phi[0][0]
    phi = np.zeros(Ex.shape[0])
    for j in range(0, Ex.shape[0]):
        phi[j] = -np.sum(Ex[0:j, 0])*dx

    AlCharge = 13.0
    AlMass = 49577.0
    HamilAl = np.zeros((q_Al.shape[0], q_Al.shape[1]))
    for j in range(0, q_Al.shape[0]):
        for k in range(0, q_Al.shape[1]):
            HamilAl[j,k] = 0.5*AlMass*(XX_Al[1][k]-5.66*vthIon)**2 + AlCharge*phi[j]
    
    figure(i+2, figsize=(14,8))
    subplot(1, 2, 1)
    pcolormesh(XX_Al[0]/lambdaD, XX_Al[1]/vthAl-5.66, transpose(q_Al_rescaled), shading="gouraud")
    colorbar(label=r"$\sqrt{\frac{f_i}{max(f_i)}}$")
    contour(XX_Al[0]/lambdaD, XX_Al[1]/vthAl-5.66, transpose(HamilAl), 1, colors="w")
    axis("tight")
    xlim(65, 90)
    ylim(-8, 0)
    xlabel("$X (\lambda_D)$")
    ylabel("$V_x (v_{th_{p}}) - V_{shock}$")
    title("Aluminum")

    ionCharge = 1.0
    ionMass = 1836.2
    HamilIon = np.zeros((q_ion.shape[0], q_ion.shape[1]))
    HamilIonConstant = np.zeros((q_ion.shape[0], q_ion.shape[1]))
    for j in range(0, q_ion.shape[0]):
        for k in range(0, q_ion.shape[1]):
            HamilIon[j,k] = 0.5*ionMass*(XX_ion[1][k]-5.66*vthIon)**2 + ionCharge*phi[j]
            
    subplot(1, 2, 2)
    pcolormesh(XX_ion[0]/lambdaD, XX_ion[1]/vthIon-5.66, transpose(q_ion_rescaled), shading="gouraud")
    colorbar(label=r"$\sqrt{\frac{f_p}{max(f_p)}}$")
    contour(XX_ion[0]/lambdaD, XX_ion[1]/vthIon-5.66, transpose(HamilIon), 8, colors="w")
    axis("tight")
    xlim(65, 90)
    ylim(-8, 8)    
    xlabel("$X (\lambda_D)$")
    ylabel("$V_x (v_{th_{p}}) - V_{shock}$")
    title("Proton impurity")
    
    tight_layout()
    savefig("vlasov-riemann-3-species-1x1v-distf-contours.pdf")
    
for i in range(35,36):
    plotFig(i,i)
    
show()

