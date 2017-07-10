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
                timeResolution=None, outputName='test_spectraPlot.png',
                colourMap='gist_ncar', animationMax=5):
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
    # Minimum and maximum available flux densities.
    minAmps = []
    maxAmps = []
    for i in xrange(0, len(spectra['spectra'])):
        for j in xrange(0, len(spectra['spectra'][i]['amp'])):
            minAmps.append(min(spectra['spectra'][i]['amp'][j]))
            maxAmps.append(max(spectra['spectra'][i]['amp'][j]))
    minAmp = min(minAmps)
    maxAmp = max(maxAmps)
    if includeZero == True:
        minAmp = 0.

    # Choose some sensible defaults.
    if frequencyResolution is None and 'frequencyResolution' in spectra:
        frequencyResolution = spectra['frequencyResolution']
    if timeResolution is None:
        # We default to 2 minutes.
        timeResolution = 2. / 1440.
        
    if plotType == "dynamic":
        # We are making a dynamic spectrum.
        # We first form a regular grid.
        # The frequency grid has 1 MHz increments.
        freqGrid = [ minFreq ]
        ff = minFreq
        fincr = 1.
        if 'frequencyUnits' in spectra and spectra['frequencyUnits'].lower() == "ghz":
            fincr = fincr / 1000.

        while ff <= maxFreq:
            ff += fincr
            freqGrid.append(ff)
        timeGrid = [ minTime ]
        tt = minTime
        while tt <= maxTime:
            tt += timeResolution
            timeGrid.append(tt)
        # The array that will hold the values.
        imageData = []
        for i in xrange(0, len(freqGrid)):
            # Accumulate all the data for this frequency.
            ftimes = []
            famps = []
            for j in xrange(0, len(spectra['spectra'])):
                for k in xrange(0, len(spectra['spectra'][j]['freq'])):
                    idx = -1
                    if (freqGrid[i] >= min(spectra['spectra'][j]['freq'][k]) - frequencyResolution and
                        freqGrid[i] <= max(spectra['spectra'][j]['freq'][k]) + frequencyResolution):
                        # There is a chance we'll find this frequency.
                        for l in xrange(0, len(spectra['spectra'][j]['freq'][k])):
                            if (abs(freqGrid[i] - spectra['spectra'][j]['freq'][k][l]) < (frequencyResolution / 2.)):
                                idx = l
                                break
                    if idx > -1:
                        ftimes.append(spectra['mjd'][j])
                        famps.append(spectra['spectra'][j]['amp'][k][idx])
            if len(ftimes) > 0:
                finter = np.interp(timeGrid, ftimes, famps, left=float('nan'), right=float('nan'))
            else:
                finter = [ float('nan') for x in timeGrid ]
            imageData.append(finter)
        # Mask the nans.
        maskedImageData = np.ma.masked_where(np.isnan(imageData), imageData)
        plt.pcolormesh(timeGrid, freqGrid, maskedImageData, vmin=minAmp, vmax=maxAmp, cmap=colourMap)
        plt.colorbar()
        plt.xlabel('MJD')
        plt.ylabel('Frequency [%s]' % spectra['frequencyUnits'])
        plt.savefig(outputName)
        plt.close()
    elif plotType == "animation":
        # We will plot the time evolution of the spectra as a series of images,
        # where each plot has the latest spectrum is coloured, and the previous are
        # coloured as lighter shades of grey.

        # Go through each of the measurements in the time series.
        recentColour = 'red'
        for i in xrange(0, len(spectra['spectra'])):
            # Plot only the most recent N spectra, and no more than animationMax.
            startIndex = i - animationMax
            if startIndex < 0:
                startIndex = 0
            minLinewidth = 1.0
            maxLinewidth = 3.0
            incrLinewidth = (maxLinewidth - minLinewidth) / float(animationMax)
            for j in xrange(startIndex, i + 1):
                lw = maxLinewidth - float(i - j) * incrLinewidth
                for k in xrange(0, len(spectra['spectra'][j]['freq'])):
                    colour = recentColour
                    if j != i:
                        colour = 'grey'
                    plt.plot(spectra['spectra'][j]['freq'][k], spectra['spectra'][j]['amp'][k], '-',
                             color=colour, lw=lw)
            plt.ylim((minAmp, maxAmp))
            plt.xlim((minFreq, maxFreq))
            plt.xlabel('Frequency [%s]' % spectra['frequencyUnits'])
            plt.ylabel('Flux Density [Jy]')
            timeLabel = "MJD = %.4f" % spectra['mjd'][i]
            plt.text(minFreq, maxAmp + ((maxAmp - minAmp) / 40.), timeLabel, color=recentColour)
            plt.savefig(outputName % i)
            plt.close()
