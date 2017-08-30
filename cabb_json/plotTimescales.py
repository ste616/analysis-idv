import json
import matplotlib.pyplot as plt
import numpy as np
import sys
import math

# Load the JSON file.
with open(sys.argv[1], 'r') as fp:
    data = json.load(fp)

# We're going to make a plot with DOY on the x-axis, and hours
# on the y-axis.

doys = []
years = []
timescales = []
timescaleErrors = []
freqs = []
limits = []

for src in data:
    print src
    srcData = data[src]
    for i in xrange(0, len(srcData)):
        scaleFactor = 1.0
        if srcData[i]['timescaleUnits'] == "minutes":
            scaleFactor = 1.0 / 60.0
        doys.append(srcData[i]['doy']['doy'])
        years.append(srcData[i]['doy']['year'])
        timescales.append(srcData[i]['timescaleValue'] * scaleFactor)
        if srcData[i]['timescaleType'] == "lower_limit":
            timescaleErrors.append(0.5)
            limits.append(1)
        else:
            timescaleErrors.append(srcData[i]['timescaleError'] * scaleFactor)
            limits.append(0)
        freqs.append(srcData[i]['frequency'])

# Change arrays to numpy.
doysnp = np.array(doys)
timescalesnp = np.array(timescales)
timescaleErrorsnp = np.array(timescaleErrors)
errorRatiosnp = (timescaleErrorsnp / timescalesnp)**2
freqsnp = np.array(freqs)
limitsnp = np.array(limits)

# Get the average time scale per epoch.
avDoys = []
avYears = []
avTimescales = []
avTimescaleErrors = []
for i in xrange(0, len(doys)):
    ddoy = int(doys[i])
    if ddoy in avDoys:
        next
    idx = np.where(np.abs(doysnp - ddoy) <= 1.)
    avgDoy = np.mean(doysnp[idx])
    avgTimescale = np.mean(timescalesnp[idx])
    avgTimescaleError = avgTimescale * math.sqrt(np.sum(errorRatiosnp[idx]))
    avDoys.append(avgDoy)
    avTimescales.append(avgTimescale)
    avTimescaleErrors.append(avgTimescaleError)

#pFreqs = list(set(freqs))
#pYears = list(set(years))
#for fq in xrange(0, len(pFreqs)):
#    idx = np.where(freqsnp == pFreqs[fq])
#    x = doysnp[idx]
#    y = timescalesnp[idx]
#    ye = timescaleErrorsnp[idx]
#    yl = limitsnp[idx]
#    plt.errorbar(x, y, yerr=ye, fmt='o', lolims=yl, label="%d MHz" % pFreqs[fq])
plt.errorbar(avDoys, avTimescales, yerr=avTimescaleErrors, fmt='o')
plt.xlim(0, 365)
plt.xlabel("DOY")
plt.ylim(0, 20)
plt.ylabel("Timescale (hours)")
#plt.legend()
plt.savefig("testts.png")
plt.close()
