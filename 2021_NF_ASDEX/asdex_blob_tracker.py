##### Movie of Electric field and density contour:
#https://jychoi-hpc.github.io/adios-python-docs/quick.html
import adios as ad
import postgkyl as pg
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('seaborn-white')
import matplotlib.animation as animation
from random import seed
from random import random
from time import sleep
from tqdm import tqdm
import os


seed(99999999)

show_anim = True
save_anim = False

interval = 0.1#in seconds
fstart = 900
fend = fstart+300#3510
dt = 1
DIR ="Run06/"
outDir = "blobDataRun06"

thresholdDensity = 3.1e18
#========== Blob Data Directory Setup =============
if os.path.exists(outDir):
    os.system('rm -rf '+outDir)
    os.system('mkdir '+outDir)
else:
    os.system('mkdir '+outDir)
############################################
data_num = np.arange(start=fstart, stop=fend, step=dt, dtype=int)
f = ad.file(DIR+'asdex_phi_%d'%data_num[0]+'.bp')

blob_size_file = open(outDir+"/blob_size.txt", "w")

Nx = f['numCells'][0]
Ny = f['numCells'][1]
Nz = f['numCells'][2]

Xmin = f['lowerBounds'][0]
Ymin = f['lowerBounds'][1]
Zmin = f['lowerBounds'][2]

Xmax = f['upperBounds'][0]
Ymax = f['upperBounds'][1]
Zmax = f['upperBounds'][2]

dx = (Xmax - Xmin) / Nx
dy = (Ymax - Ymin) / Ny

z_slice = 10
cnum = 100
cnumout = 30
color = 'jet'


################### INTERPOLATION ###########################

def interpTestPoint(xWeight,yWeight,dx,dy,probeDensity):
    testDensity00 = probeDensity[0,0] * (dx-xWeight) * (dy-yWeight)
    testDensity01 = probeDensity[0,1] * xWeight * (dy-yWeight)
    testDensity10 = probeDensity[1,0] * (dx-xWeight) * yWeight
    testDensity11 = probeDensity[1,1] * xWeight * yWeight
    testDensity = ( testDensity00 + testDensity01 + testDensity10 + testDensity11 ) / (dx*dy)
    return testDensity

################### Shoelace formula to find polygon Area ###########################

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
####################################################################################
#################### RAY TRACING ALGORITHM #########################################
####################################################################################

# A Python3 program to check if a given point lies inside a given polygon
# Refer https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
# for explanation of functions onSegment(), orientation() and doIntersect()

# Define Infinite (Using INT_MAX caused overflow problems)
INF = 10000

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

# Given three colinear points p, q, r, the function checks if
# point q lies on line segment 'pr'
def onSegment(p, q, r):
    if ( (q.x <= max(p.x, r.x)) and (q.x >= min(p.x, r.x)) and
           (q.y <= max(p.y, r.y)) and (q.y >= min(p.y, r.y))):
        return True
    return False

def orientation(p, q, r):
    # to find the orientation of an ordered triplet (p,q,r)
    # function returns the following values:
    # 0 : Colinear points
    # 1 : Clockwise points
    # 2 : Counterclockwise

    # See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/
    # for details of below formula.

    val = (float(q.y - p.y) * (r.x - q.x)) - (float(q.x - p.x) * (r.y - q.y))
    if (val > 0):
        # Clockwise orientation
        return 1
    elif (val < 0):
        # Counterclockwise orientation
        return 2
    else:
        # Colinear orientation
        return 0

# The main function that returns true if
# the line segment 'p1q1' and 'p2q2' intersect.
def doIntersect(p1,q1,p2,q2):

    # Find the 4 orientations required for
    # the general and special cases
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if ((o1 != o2) and (o3 != o4)):
        return True

    # Special Cases

    # p1 , q1 and p2 are colinear and p2 lies on segment p1q1
    if ((o1 == 0) and onSegment(p1, p2, q1)):
        return True

    # p1 , q1 and q2 are colinear and q2 lies on segment p1q1
    if ((o2 == 0) and onSegment(p1, q2, q1)):
        return True

    # p2 , q2 and p1 are colinear and p1 lies on segment p2q2
    if ((o3 == 0) and onSegment(p2, p1, q2)):
        return True

    # p2 , q2 and q1 are colinear and q1 lies on segment p2q2
    if ((o4 == 0) and onSegment(p2, q1, q2)):
        return True

    # If none of the cases
    return False

# Returns true if the point p lies inside the polygon[] with n vertices
def isInside(polygon, n, p):
    # There must be at least 3 vertices in polygon[]
    if (n < 3):
        return False

    # Create a point for line segment from p to infinite
    extreme = Point(INF, p.y)

    # Count intersections of the above line with sides of polygon
    count = 0
    i = 0

    # To initialize i for the first iteration of do-while loop of C++ type
    next = (i+1)%n
    # Check if the line segment from 'p' to 'extreme' intersects
    # with the line segment from 'polygon[i]' to 'polygon[next]'
    if (doIntersect(polygon[i], polygon[next], p, extreme)):
        # If the point 'p' is colinear with line segment 'i-next',
        # then check if it lies on segment. If it lies, return true,
        # otherwise false
        if (orientation(polygon[i], p, polygon[next]) == 0):
            return onSegment(polygon[i], p, polygon[next])
        count = count + 1
    i = next

    while (i != 0):
        next = (i+1)%n
        # Check if the line segment from 'p' to 'extreme' intersects
        # with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], p, extreme)):
            # If the point 'p' is colinear with line segment 'i-next',
            # then check if it lies on segment. If it lies, return true,
            # otherwise false
            if (orientation(polygon[i], p, polygon[next]) == 0):
                return onSegment(polygon[i], p, polygon[next])
            count = count + 1
        i = next
        if (i == 0):
            break
    # Return true if count is odd, false otherwise
    if (count%2 == 1):
        return True
    else:
        return False

####################################################################################
####################################################################################
####################################################################################

def func_data(ionDensityData,phiData):
	ionDensityInterp = pg.data.GInterpModal(ionDensityData, 1, 'ms')
	phiInterp = pg.data.GInterpModal(phiData, 1, 'ms')
	interpGrid, ionDensityValues = ionDensityInterp.interpolate()
	interpGrid, phiValues = phiInterp.interpolate()

	#exValues = - np.gradient(phiValues,dx,axis = 0)
	#dexdxValues = np.gradient(exValues,dx,axis = 0)
	eyValues = - np.gradient(phiValues,dy,axis = 1)

	# get cell center coordinates
	CCC = []
	for j in range(0,len(interpGrid)):
	    CCC.append((interpGrid[j][1:] + interpGrid[j][:-1])/2)

	x_vals = CCC[0]
	y_vals = CCC[1]
	z_vals = CCC[2]
	X, Y = np.meshgrid(x_vals, y_vals)
	ionDensityGrid = np.transpose(ionDensityValues[:,:,z_slice,0])
	eyGrid = np.transpose(eyValues[:,:,z_slice,0])
	return x_vals,y_vals,X,Y,ionDensityGrid,eyGrid

def animate(i):
        blob_counter = 0
        ionDensity=DIR+'asdex_ion_GkM0_%d'%data_num[i]+'.bp'
        phi=DIR+'asdex_phi_%d'%data_num[i]+'.bp'
        ionDensityData = pg.data.GData(ionDensity)
        phiData = pg.data.GData(phi)

        x_vals,y_vals,X,Y,ionDensityGrid,eyGrid = func_data(ionDensityData,phiData)

        Nx = len(x_vals)
        Ny = len(y_vals)

        ax1.cla()
        ax1.set_title('Time = %d'%i+' $\\mu$s')

        cp1 = ax1.contourf(X, Y, ionDensityGrid, cnum, cmap=color)
        #cp2 = ax1.contour(X, Y, eyGrid, cnum, linewidths=0.1, colors='black', linestyles='solid')
        cp3 = ax1.contour(X, Y, ionDensityGrid, cnumout, linewidths=0.1, colors='black', linestyles='solid')
        #cp3 = ax1.contour(X, Y, ionDensityGrid, cnumout, linewidths=1, cmap=color)
        cp4 = ax1.contour(X, Y, ionDensityGrid, [thresholdDensity], linewidths=1, colors='black', linestyles='solid')
        # #plt.grid()
        # ax1.set_xticks(x_vals , minor=True)
        # ax1.set_yticks(y_vals , minor=True)
        # #ax1.grid(which='both')
        # ax1.grid(which='minor', alpha=0.9, color='k', linestyle='-')

        p = cp4.collections[0].get_paths()
        contour_number = len(p)
        imageCounter = 0
        for j in range(contour_number):
            p_new = cp4.collections[0].get_paths()[j]
            v = p_new.vertices
            x = v[:,0]
            y = v[:,1]
            x_min = np.min(x)
            x_max = np.max(x)
            y_min = np.min(y)
            y_max = np.max(y)
            blobMidX = (x_min + x_max)/2
            blobMidY = (y_min + y_max)/2
            blobLimX = abs(x_max - x_min)
            blobLimY = abs(y_max - y_min)
            if (abs(x[0]-x[len(x)-1]) <= 1e-10) and blobLimX > 2*dx and blobLimY > 2*dy:
                polygon = []
                for plgn in range(len(x)):
                    polygon.append(Point(x[plgn],y[plgn]))
                npoly = len(polygon)
                numTrial = 100
                blobConfidence = 0
                insideTrialPoints = 0
                for numT in range(numTrial):
                    xT = 0.5*(x_max+x_min) - 0.5*(x_max-x_min)*(random()-0.5)
                    yT = 0.5*(y_max+y_min) - 0.5*(y_max-y_min)*(random()-0.5)
                    #print("Trial point",numT,"with",round(xT,4),round(yT,4),'for contour number %d'%j)
                    trialPoint = Point(xT,yT)
                    if isInside(polygon, npoly, trialPoint):
                        insideTrialPoints = insideTrialPoints + 1
                        #print("Trial point", numT, "is INSIDE for contour number %d"%j)
                        xd = abs(x_vals-xT)
                        yd = abs(y_vals-yT)
                        idx = np.where(xd <= 0.5*dx)
                        idy = np.where(yd <= 0.5*dy)
                        ionDensityFind = np.reshape(ionDensityGrid,Nx*Ny)
                        probeDensity = np.zeros((2,2))
                        for id in range(len(idx[0])):
                            for jd in range(len(idy[0])):
                                probeDensity[id,jd] = ionDensityFind[(idy[0][jd] * Nx) + (idx[0][id] + 1)]

                        xGrid = np.zeros(2)
                        yGrid = np.zeros(2)
                        for id in range(len(idx[0])):
                            xGrid[id] = x_vals[idx[0][id]]
                        for jd in range(len(idy[0])):
                            yGrid[jd] = y_vals[idy[0][jd]]

                        xWeight = abs(xGrid[0]-xT)
                        yWeight = abs(yGrid[0]-yT)
                        testDensity = interpTestPoint(xWeight,yWeight,dx,dy,probeDensity)
                        if (testDensity >= thresholdDensity):
                                #print("Interpolated point",numT,"with",round(xInterp,4),round(yInterp,4)," for Contour number %d"%j+" is INSIDE & truly a BLOB! Yeyy...")
                                blobConfidence = blobConfidence + 1

                        else:
                                None
                    else:
                        None
                        #print("Trial point", numT, " lies Outside before interpolation")

                confidence = blobConfidence/insideTrialPoints
                #print("Confidence = ",confidence*100,"%")
                if (confidence > 0.80):
                    blob_counter = blob_counter + 1
                    polyArea = PolyArea(x,y)
                    #print(polyArea)
                    # print('File number = %d'%data_num[i]+', contour number %d'%j+' = It is TRULY a blob with confidence',confidence*100,"%")
                    blob_size_file.write('%d'%data_num[i]+'\t%d'%j+'\t%.8f'%blobLimX+'\t%.8f'%blobLimY+'\t%.8f'%blobMidX+'\t%.8f'%blobMidY+'\t%.8f'%polyArea+'\n')
                    if imageCounter == 0:
                        plt.savefig(outDir+"/file_number%d"%data_num[i]+"_blob_snap.png")   # save the figure to file
                    imageCounter = imageCounter + 1
                    #print("blobConfidence=",blobConfidence,"insideTrialPoints=",insideTrialPoints)
                    blob_file = open(outDir+"/file_number%d"%data_num[i]+"_contour_number_%d"%j+".txt", "w")
                    for k in range(len(x)):
                        blob_file.write('%.8f'%x[k]+'\t%.8f'%y[k]+'\n')
                    blob_file.close()
            elif (abs(x[0]-x[len(x)-1]) <= 1e-10):
                None
                # print('File number = %d'%data_num[i]+', contour number %d'%j+' = It is a sub-grid-sized closed contour')
            else:
                None
                # print('File number = %d'%data_num[i]+', contour number %d'%j+' = It is open line & NOT a blob')
                #print(x,y)
                #for k in range(len(x)):
                    #print(round(x[k],7),round(y[k],7))

        if blob_counter == 0:
            None
            # print("No blob found for file number = %d"%data_num[i])
        sleep(0.1)
        pbar.update(pstep)
        #plt.grid(True)
        ax1.set_xlabel("X",fontsize=14)
        ax1.set_ylabel("Y",fontsize=14)
        #ax1.tick_params(axis='both', which='major', labelsize=12)
        del ionDensityData
        del phiData

if (show_anim == True):
    fig,ax1 = plt.subplots(1,1,figsize=(8,5),dpi=150)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["font.family"] = "Times New Roman"
    #To keep the colorbar static:
    ionDensity=DIR+'asdex_ion_GkM0_%d'%data_num[0]+'.bp'
    phi=DIR+'asdex_phi_%d'%data_num[0]+'.bp'
    ionDensityData = pg.data.GData(ionDensity)
    phiData = pg.data.GData(phi)
    x_vals,y_vals,X,Y,ionDensityGrid,eyGrid = func_data(ionDensityData,phiData)
    cp1 = ax1.contourf(X, Y, ionDensityGrid, cnum, cmap=color)
    fig.colorbar(cp1)
    #TColorbar fixing completed:
    pstep = 1#len(data_num)/100
    pbar = tqdm(total=len(data_num))
    ani = animation.FuncAnimation(fig,animate,frames=len(data_num),interval=interval*1e+3,blit=False,repeat=False)
    ax1.set_xticks(x_vals , minor=True)
    ax1.set_yticks(y_vals , minor=True)
    ax1.grid(which='both')
    ax1.grid(which='minor', alpha=0.2, color='b', linestyle='--')
    #ax1.grid(b=True, which='major', color='b', linestyle='-')
    plt.show()
    if(save_anim == True):
        try:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=(1/interval), metadata=dict(artist='Me'), bitrate=1800)
        except RuntimeError:
            print("ffmpeg not available trying ImageMagickWriter")
            writer = animation.ImageMagickWriter(fps=(1/interval))
        ani.save('animation.mp4')
pbar.close()
blob_size_file.close()
