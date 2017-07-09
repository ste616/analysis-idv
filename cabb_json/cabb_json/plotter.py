# This is analysis-idv/cabb_json/plotter.py
# Jamie Stevens 2017
# Licensed under GPL 3.0 or later. A copy of the license can
# be found at https://www.gnu.org/licenses/gpl.html

# This module contains routines to take spectra and time series and make
# useful and good-looking plots.

# Our necessary imports.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata

# Plot a set of spectra from a getSpecta method, in some ways.
def spectraPlot(spectra, timeRange=None, plotType="dynamic", frequencyRange=None,
                ampRange=None, includeZero=False, frequencyResolution=None,
                timeResolution=None):
    # Work out some ranges for the data.
    # Minimum and maximum available times.
    minTime = min(spectra['mjd'])
    maxTime = max(spectra['mjd'])
    # Minimum and maximum available frequencies.
    minFreqs = []
    maxFreqs = []
    for i in xrange(0, len(spectra['spectra'])):
        for j in xrange(0, len(spectra['spectra'][i]['freq'])):
            minFreqs.append(min(spectra['spectra'][i]['freq'][j]))
            maxFreqs.append(max(spectra['spectra'][i]['freq'][j]))
    minFreq = min(minFreqs)
    maxFreq = max(maxFreqs)
    print "minimum / maximum frequencies %.3f / %.3f" % (minFreq, maxFreq)
    # Minimum and maximum available flux densities.
    minAmps = []
    maxAmps = []
    for i in xrange(0, len(spectra['spectra'])):
        for j in xrange(0, len(spectra['spectra'][i]['amp'])):
            minAmps.append(min(spectra['spectra'][i]['amp'][j]))
            maxAmps.append(max(spectra['spectra'][i]['amp'][j]))
    minAmp = min(minAmps)
    maxAmp = max(maxAmps)

    # Choose some sensible defaults.
    if frequencyResolution is None and 'frequencyResolution' in spectra:
        frequencyResolution = spectra['frequencyResolution']
        print "found frequency resolution of %.3f" % frequencyResolution
    if timeResolution is None:
        # We default to 2 minutes.
        timeResolution = 2. / 1440.
        
    if plotType == "dynamic":
        # We are making a dynamic spectrum.
        ## Assemble all the points into three areas, mjd, freq and amp.
        #mjd = []
        #freq = []
        #amp = []

        #for i in xrange(0, len(spectra['spectra'])):
        #    for j in xrange(0, len(spectra['spectra'][i]['freq'])):
        #        for k in xrange(0, len(spectra['spectra'][i]['freq'][j])):
        #            mjd.append(spectra['mjd'][i])
        #            freq.append(spectra['spectra'][i]['freq'][j][k])
        #            amp.append(spectra['spectra'][i]['amp'][j][k])
        ## Now mesh the frequency and mjd into a grid.
        #mjdMesh, freqMesh = np.meshgrid(mjd, freq)
        ## And grid the amps onto this grid.
        #ampMesh = griddata(freq, mjd, amp, freqMesh, mjdMesh, interp='nn')
        ## And plot.
        #plt.pcolor(ampMesh)
        #plt.colorbar()
        #plt.savefig("test_grid.png")
        #plt.close()
        
        # We first form a regular grid.
        # The frequency grid has 1 MHz increments.
        freqGrid = [ minFreq ]
        ff = minFreq
        while ff <= maxFreq:
            if 'frequencyUnits' in spectra and spectra['frequencyUnits'].lower() == "ghz":
                ff += 0.001
            else:
                ff += 1.
            freqGrid.append(ff)
        #print freqGrid
        timeGrid = [ minTime ]
        tt = minTime
        while tt <= maxTime:
            tt += timeResolution
            timeGrid.append(tt)
        # The array that will hold the values.
        imageData = []
        for i in xrange(0, len(freqGrid)):
            # Accumulate all the data for this frequency.
            #print "frequency %.3f" % freqGrid[i]
            ftimes = []
            famps = []
            for j in xrange(0, len(spectra['spectra'])):
                #print " searching spectra %d" % j
                for k in xrange(0, len(spectra['spectra'][j]['freq'])):
                    #print "  searching frequency spectrum %d" % k
                    idx = -1
                    #print " min max %.3f %.3f" % (min(spectra['spectra'][j]['freq'][k]),
                    #                              max(spectra['spectra'][j]['freq'][k]))
                    if (freqGrid[i] >= min(spectra['spectra'][j]['freq'][k]) - frequencyResolution and
                        freqGrid[i] <= max(spectra['spectra'][j]['freq'][k]) + frequencyResolution):
                        # There is a chance we'll find this frequency.
                        print "   continuing search"
                        for l in xrange(0, len(spectra['spectra'][j]['freq'][k])):
                            if (abs(freqGrid[i] - spectra['spectra'][j]['freq'][k][l]) < (frequencyResolution / 2.)):
                                idx = l
                                break
                    if idx > -1:
                        print "   found index %d" % idx
                        ftimes.append(spectra['mjd'][j])
                        famps.append(spectra['spectra'][j]['amp'][k][idx])
            print "arrays follow"
            print ftimes
            print famps
            if len(ftimes) > 0:
                finter = np.interp(timeGrid, ftimes, famps, left=float('nan'), right=float('nan'))
            else:
                finter = [ float('nan') for x in timeGrid ]
            imageData.append(finter)
        #print imageData
        print len(freqGrid)
        print len(timeGrid)
        print len(imageData)
        print len(imageData[0])
        print "image data follows"
        print imageData
        #plt.imshow(imageData, extent=(minTime, maxTime, minFreq, maxFreq), origin='lower', aspect='auto',
        #           vmin=minAmp, vmax=maxAmp)
        plt.pcolormesh(timeGrid, freqGrid, imageData, vmin=minAmp, vmax=maxAmp)
        plt.colorbar()
        plt.savefig("test_grid.png")
        plt.close()
