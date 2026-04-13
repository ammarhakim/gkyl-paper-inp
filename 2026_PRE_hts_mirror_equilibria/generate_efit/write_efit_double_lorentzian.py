import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union

# Remember to delete the top line and that Gkeyll takes the _psi file which is made from the efit file reader

# Parameters for different mirror ratios R
mcB_values   = [2.130115, 2.665626, 3.691260, 4.490901, 5.416264, 6.51292, 7.274615, 8.125522]
gamma_values = [0.451454, 0.331696, 0.226381, 0.182792, 0.149893, 0.124904, 0.110435, 0.098619]
R_values     = [3, 5, 10, 15, 22, 32, 40, 50]

B0 = 0.0
R0 = .02
Z_m = 0.98  # Same for all mirror ratios

def psi_f(R, Z, mcB, gamma):
  return 0.5*np.power(R,2)*mcB*( 1./(np.pi*gamma*(1.+((Z-Z_m)/gamma)**2)) \
                                +1./(np.pi*gamma*(1.+((Z+Z_m)/gamma)**2)) )

#RZ box
NW = 257
NH = 257
RMIN,RMAX = 1e-3, .3 # The _zero is for Rmin=0, but 1e-3 was used for 1x cases
ZMIN,ZMAX = -2.5, 2.5
RDIM = RMAX - RMIN
ZDIM = ZMAX - ZMIN
RLEFT = RMIN
ZMID = (ZMAX+ZMIN)/2.0
RMAXIS = 0.0
ZMAXIS = 0.0
NPSI = NW #don't write
RCENTR = 0.0
BCENTR = B0
CURRENT = 0

#Solve GS in RZ coords
Rgrid = np.linspace(RMIN,RMAX,NW)
Zgrid = np.linspace(ZMIN,ZMAX,NH)

#Header stuff
header_fmt = "(a48,3i4)"
label = 'FREEGS'
creation_date = date.today().strftime("%d/%m/%Y")
shot = int(0)
time = int(0)
shot_str = f"# {shot:d}"
time_str = f"  {time:d}ms"
comment = f"{label:11}{creation_date:10s}   {shot_str:>8s}{time_str:16s}"

def write_line(data: Iterable[Any], fh: TextIO, fmt: str) -> None:
    r"""
    Writes to a Fortran formatted ASCII data file. The file handle will be left on a
    newline.

    Parameters
    ---------
    data:
        The data to write.
    fh:
        File handle. Should be in a text write mode, i.e. ``open(filename, "w")``.
    fmt:
        A Fortran IO format string, such as ``'(6a8,3i3)'``.
    """
    fh.write(ff.FortranRecordWriter(fmt).write(data))
    fh.write("\n")

# Loop over all mirror ratios
for idx, (mcB, gamma, R_mirror) in enumerate(zip(mcB_values, gamma_values, R_values)):
    outFileName = f'lorentzian_R{R_mirror}.geqdsk'
    
    print(f"\n--- Processing R = {R_mirror} ---")
    print(f"mcB = {mcB}, gamma = {gamma}")
    
    SIMAG = psi_f(2.0, 0, mcB, gamma) 
    SIBRY = psi_f(4.0, 0, mcB, gamma)
    
    print("simag = %g"%SIMAG)
    print("sibry = %g"%SIBRY)
    
    #rthetagrid = np.zeros((len(Rgrid),len(Zgrid),2))
    psiRZ = np.zeros((len(Rgrid),len(Zgrid)))
    for i,Ri in enumerate(Rgrid):
        for j,Zj in enumerate(Zgrid):
            psiRZ[i,j] = psi_f(Ri, Zj, mcB, gamma)
    
    plt.figure()
    plt.contour(Rgrid,Zgrid, psiRZ.T, levels = (np.linspace(0.001,0.004,12)))
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.colorbar()
    plt.title(f"Psi calculated in RZ coords (R={R_mirror})")
    
    #PSI quantities
    PSIGRID = np.linspace(SIMAG, SIBRY,NPSI)
    FPOL = (B0*R0/Rgrid)*Rgrid # F = RBphi
    FFPRIM = np.repeat(0.0, NPSI)
    PPRIME = np.repeat(-1e-6,NPSI)
    PRES = integrate.cumulative_trapezoid(PPRIME,PSIGRID,initial=0)
    PSIZR = psiRZ.T
    
    QPSI = np.zeros(NPSI)
    for i in range(NPSI):
        QPSI[i] = 0
    
    writeList = [NW, NH,                                        #3i4
                 RDIM, ZDIM, RCENTR, RLEFT, ZMID,               #5E16.9
                 RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR,          #5e16.9
                 CURRENT, SIMAG, 0, RMAXIS, 0,                  #5E16.9
                 ZMAXIS, 0, SIBRY, 0, 0,                        #5E16.9
                 FPOL, PRES, FFPRIM, PPRIME,                    #5E16.9
                 PSIZR, QPSI]                                   #5E16.9
    
    #Now write the EFIT FILE
    with open(outFileName,'w',newline='') as f:
        write_line((comment, 3, NW, NH), f, header_fmt) #3 is idum
        # rdim,zdim,rcentr,rleft,zmid
        for i in range(2,7):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        # rmaxis,zmaxis,simag,sibry,bcentr
        for i in range(7,12):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        # current, simag, xdum, rmaxis, xdum
        for i in range(12,17):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        #zmaxis,xdum,sibry,xdum,xdum
        for i in range(17,22):
             f.write('%16.9E'%writeList[i])
        f.write('\n')
        #FPOL,PRES,FFPRIM,PPRIME
        for i in range(22,26):
            count=0
            for j in range(0,NW):
                f.write('%16.9E'%writeList[i][j])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0
        #PSIZR
        for i in range(26,27):
            count = 0
            for j in range(0,NH):
                for k in range(0,NW):
                    f.write('%16.9E'%writeList[i][j][k])
                    count = count+1
                    if count==5:
                        f.write('\n')
                        count=0
        #QPSI
        for i in range(27,28):
            count=0
            for j in range(0,NW):
                f.write('%16.9E'%writeList[i][j])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0
    
    print(f"Written: {outFileName}")

plt.show()
