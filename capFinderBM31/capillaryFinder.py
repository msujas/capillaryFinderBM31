import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit



file = r'zscan.dat'

def readZscan(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    header = 'zpos i1 monitor'
    for c,line in enumerate(lines[::-1]):
        if header in line:
            break
    skipline = len(lines)-c
    zpos,i1s, mons = np.loadtxt(filename,skiprows = skipline,unpack=True)
    return zpos, i1s, mons

def deriv(filename):
    zpos,i1s,mons = readZscan(filename)
    i1norm = i1s/mons
    return np.gradient(i1norm)

def peakFind(x,y, capsize = 1):
    ymean = np.average(y)
    ystdev = np.std(y)
    peaks = []
    for i in range(len(y[:-1])):
        if (y[i]-ymean)**2  > (y[i-1]-ymean)**2 and (y[i]-ymean)**2 > (y[i+1]-ymean)**2 and np.abs(y[i]-ymean) > ystdev:
            if not peaks or np.abs(x[i]-peaks[-1]) > capsize:
                peaks.append(x[i])

    return np.array(peaks)

def dpeakPair(peaks):
    capPositions = []
    i = 1
    while i < len(peaks)-1:
        if (peaks[i] - peaks[i-1])**2 < (peaks[i] - peaks[i+1])**2 and np.abs(peaks[i] - peaks[i-1]) < 3:
            capPos = (peaks[i] + peaks[i-1])/2
            capPositions.append(capPos)
            i += 1
        elif (peaks[i] - peaks[i+1])**2 < (peaks[i] - peaks[i-1])**2 and np.abs(peaks[i] - peaks[i+1]) < 3:
            capPos = (peaks[i] + peaks[i+1])/2
            capPositions.append(capPos)
            i += 2
        else:
            i+=1
            
    return np.array(capPositions)

def gauss(x,a,x0,c,d):
    return (a*np.exp(-(x-x0)**2/(2*c**2))) + d

def refineCapPositions(capPositions,x,y, capsize = 1):
    ymean = np.mean(y)
    capPosRefined = []
    xarray = np.array([])
    yfitArray = np.array([])
    minindex = np.abs(x-capPositions[0]+capsize).argmin()
    maxindex = np.abs(x-capPositions[1]-capsize).argmin()
    baseline = np.mean(y[minindex:maxindex])
    yround = np.round(y,3)
    ymode = stats.mode(yround)[0]
    for cap in capPositions:
        minindex = np.abs(x-(cap-capsize*2)).argmin()
        maxindex = np.abs(x-(cap+capsize*2)).argmin()
        peakindex = np.abs(x-cap).argmin()
        yval = y[peakindex] - baseline
        pguess = [yval,cap,capsize/2.35,ymode]

        limits = ([-np.inf, cap-capsize/2, 0, ymode - 0.002], [np.inf,cap+capsize/2, (capsize/2.35)*1.1, ymode+0.002])
        popt,pcov = curve_fit(gauss,x[minindex:maxindex],y[minindex:maxindex], p0=pguess, maxfev = 100000, bounds=limits)
        capPosRefined.append(popt[1])
        xfit = x[minindex:maxindex]
        yfit = gauss(xfit,*popt)
        xarray = np.append(xarray,xfit)
        yfitArray = np.append(yfitArray,yfit)
    return np.array(capPosRefined), xarray, yfitArray

def run(filename, capsize = 1):
    z,i1,mon = readZscan(filename)
    i1norm = i1/mon
    #d = deriv(file)
    peaks = peakFind(z,i1norm)
    capPositions = peaks #dpeakPair(peaks)
    capRefined, xfit,yfit = refineCapPositions(capPositions,z,i1norm, capsize)
    return capRefined, xfit, yfit

def getPositions(filename, capsize = 1):
    capRefined,x,y = run(filename, capsize)
    return capRefined

def plotResults(filename, capsize = 1, block = True):
    z, i1, mon = readZscan(filename)
    i1norm = i1/mon
    capRefined, xfit,yfit = run(filename,capsize)
    plt.figure()
    plt.plot(z,i1norm)
    plt.plot(xfit,yfit)
    plt.plot(capRefined, np.zeros(len(capRefined))+np.average(i1norm),'o')
    plt.bar(capRefined,np.zeros(len(capRefined))+np.max(i1norm), width=0.02)
    plt.ylim(np.min(i1norm) - (np.max(i1norm)-np.min(i1norm))*0.1, np.max(i1norm) + (np.max(i1norm)-np.min(i1norm))*0.1)
    plt.ylabel('i1/monitor')
    plt.xlabel('z position (mm)')
    plt.show(block=block)
    return capRefined


if __name__ == '__main__':
    capPositions = plotResults(file, 1.1)
    print(capPositions)