#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# customization of plotting
plt.rcParams['lines.linewidth']            = 4
plt.rcParams['font.size']                  = 18
plt.rcParams['font.weight']                = 'bold'
plt.rcParams['axes.labelsize']             = 'large'
plt.rcParams['figure.facecolor']           = 'white'
plt.rcParams['image.interpolation']        = 'none'
plt.rcParams['image.origin']               = 'lower'
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams['mathtext.default'] = 'regular' # match the font used for
                                             # regular text

# Defining exponential
def exponential(x, a, b):
    return a*np.exp(2*b*x)

files = ["k_0_2/5m_bz2.txt"]

# Loop over files
for i, fileName in enumerate(files):
    minN = 100
    maxR_squared = 0
    data = np.loadtxt(fileName)
    for n in np.linspace(minN, data.shape[0]-1, data.shape[0]-minN):
        x = data[0:n,0]
        y = data[0:n,1]
        try:
            # Main fitting
            popt, pcov = curve_fit(exponential, x, y, p0=(1,0.1))
            # Calculation of R^2
            residuals = y-exponential(x, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y-np.mean(y))**2)
            R_squared = 1-ss_res/ss_tot
            if R_squared > maxR_squared:
                maxR_squared = R_squared
                bestFit = popt
            percentage = float(n-minN)/(data.shape[0]-minN)*100
            print('{} {:5.2f} %: gamma = {:6.4f} ({:6.4f}); R2 = {:6.4f} ({:6.4f})'.format(fileName, percentage, popt[1], bestFit[1], R_squared, maxR_squared))
        except RuntimeError:
            print("Error - curve_fit failed")
    print('{}: gamma = {}'.format(fileName, bestFit[1]))
    plt.figure(i)
    plt.plot(data[:,0], data[:,1], '.')
    plt.plot(data[:,0], exponential(data[:,0], *bestFit))
    plt.title(fileName)
    plt.xlabel('Time [Ut]')
    plt.ylabel('B_z^2')
    plt.ylim([0,np.max(data[:,1])])
    #plt.draw()
plt.show()

    
        
    
