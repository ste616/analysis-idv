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
import datetime
from astropy.time import Time
from analysis import userGauss
import math

# Plot a set of spectra from a getSpecta method, in some ways.
def spectraPlot(spectra, timeRange=None, plotType="dynamic", frequencyRange=None,
                ampRange=None, includeZero=False, frequencyResolution=None,
                timeResolution=None, outputName='test_spectraPlot.png',
                colourMap='gist_ncar', animationMax=5, title=None):
    # Work out some ranges for the data.
    minTime = 0.
    maxTime = 0.
    # Minimum and maximum available times.
    if (spectra['timeUnits'].lower() == "mjd" or spectra['timeUnits'].lower() == "astro" or
        spectra['timeUnits'].lower() == "datetime"):
        minTime = min(spectra['time'])
        maxTime = max(spectra['time'])
    elif spectra['timeUnits'].lower() == "doy":
        minTime = spectra['time'][0]
        maxTime = spectra['time'][0]
        for i in xrange(0, len(spectra['time'])):
            cmin = minTime['year'] + minTime['doy'] / 365.
            cmax = maxTime['year'] + maxTime['doy'] / 365.
            t = spectra['time'][i]['year'] + spectra['time'][i]['doy'] / 365.
            if t < cmin:
                cmin = spectra['time'][i]
            if t > cmax:
                cmax = spectra['time'][i]
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
    if includeZero == True or minAmp < 0:
        minAmp = 0.
    if ampRange is not None:
        minAmp = ampRange[0]
        maxAmp = ampRange[1]

    # Choose some sensible defaults.
    if frequencyResolution is None and 'frequencyResolution' in spectra:
        frequencyResolution = spectra['frequencyResolution']
    if timeResolution is None:
        # We default to 2 minutes.
        if spectra['timeUnits'].lower() == "mjd":
            timeResolution = 2. / 1440.
        elif spectra['timeUnits'].lower() == "datetime":
            timeResolution = datetime.timedelta(minutes=2)
        
    if plotType == "dynamic":
        # We are making a dynamic spectrum.
        # We only support one type of time labelling.
        if spectra['timeUnits'].lower() != 'mjd':
            print "Dynamic spectrum plotter supports only MJD labelling."
            return
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
        mjdDay = math.floor(minTime)
        timeGridHours = [ (minTime - mjdDay) * 24. ]
        tt = minTime
        while tt <= maxTime:
            tt += timeResolution
            timeGrid.append(tt)
            timeGridHours.append((tt - mjdDay) * 24.)
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
                        try:
                            ftimes.append(spectra['time'][j].value)
                        except AttributeError:
                            ftimes.append(spectra['time'][j])
                        famps.append(spectra['spectra'][j]['amp'][k][idx])
            if len(ftimes) > 0:
                finter = np.interp(timeGrid, ftimes, famps, left=float('nan'), right=float('nan'))
            else:
                finter = [ float('nan') for x in timeGrid ]
            imageData.append(finter)
        # Mask the nans.
        maskedImageData = np.ma.masked_where(np.isnan(imageData), imageData)
        plt.pcolormesh(np.array(timeGridHours), np.array(freqGrid), maskedImageData, vmin=minAmp, vmax=maxAmp, cmap=colourMap)
        plt.colorbar()
        if spectra['timeUnits'].lower() == "mjd":
            plt.xlabel('Hours')
        elif spectra['timeUnits'].lower() == "datetime":
            plt.xlabel('Date')
            plt.gcf().autofmt_xdate()
        plt.ylabel('Frequency [%s]' % spectra['frequencyUnits'])
        if title is None and 'sourceName' in spectra:
            # Default to just putting the name of the source on the top.
            plt.title(spectra['sourceName'])
        elif title is not None:
            plt.title(title)
        plt.savefig(outputName)
        plt.close()
        return outputName
    elif plotType == "animation":
        # We will plot the time evolution of the spectra as a series of images,
        # where each plot has the latest spectrum is coloured, and the previous are
        # coloured as lighter shades of grey.
        # We support three types of time labelling.
        if (spectra['timeUnits'].lower() != "mjd" and spectra['timeUnits'].lower() != "datetime" and
            spectra['timeUnits'].lower() != "doy"):
            print "Animation plotter supports MJD, DateTime and DOY plotting only."
            return
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
            timeLabel = ""
            if spectra['timeUnits'].lower() == "mjd":
                timeLabel = "MJD = %.4f" % spectra['time'][i]
            elif spectra['timeUnits'].lower() == "datetime":
                timeLabel = "Date = %s" % spectra['time'][i].strftime("%Y-%m-%d %H:%M:%S")
            elif spectra['timeUnits'].lower() == "doy":
                timeLabel = "DOY = %d / %.4f" % (spectra['time'][i]['year'], spectra['time'][i]['doy'])
            plt.text(minFreq, maxAmp + ((maxAmp - minAmp) / 40.), timeLabel, color=recentColour)
            plt.savefig(outputName % i)
            plt.close()

# A routine to plot up time series.
def timeSeriesPlot(timeSeries, timeRange=None, plotType="fluxDensity",
                   includeZero=False, ampRange=None,
                   outputName='test_timeSeriesPlot.png', title=None,
                   plotLegend=True):

    # Get some details about the data.
    # Minimum and maximum available times.
    minTime = 0.
    maxTime = 0.
    if timeSeries['timeUnits'].lower() == "doy":
        # The minimum and maximum are set by the number of days in a year!
        minTime = 0.
        maxTime = 366.
    else:
        # Transfer times into hours, if we need to.
        minTime = min(timeSeries['times'][0])
        maxTime = max(timeSeries['times'][0])
        for i in xrange(1, len(timeSeries['times'])):
            tmin = min(timeSeries['times'][i])
            tmax = max(timeSeries['times'][i])
            if tmin < minTime:
                minTime = tmin
            if tmax > maxTime:
                maxTime = tmax
        mjdDay = math.floor(minTime)
        minTime = (minTime - mjdDay) * 24.
        maxTime = (maxTime - mjdDay) * 24.
    # Minimum and maximum flux densities.
    minFlux = min(timeSeries['fluxDensities'][0])
    maxFlux = max(timeSeries['fluxDensities'][0])
    for i in xrange(1, len(timeSeries['fluxDensities'])):
        fmin = min(timeSeries['fluxDensities'][i])
        fmax = max(timeSeries['fluxDensities'][i])
        if fmin < minFlux:
            minFlux = fmin
        if fmax > maxFlux:
            maxFlux = fmax
    if includeZero == True or minFlux < 0:
        minFlux = 0.
    if ampRange is not None:
        minFlux = ampRange[0]
        maxFlux = ampRange[1]
    # Make the plot.
    # Set up a colour map.
    ax = plt.axes()
    ax.set_prop_cycle('color', plt.cm.gist_rainbow(
        np.linspace(0, 1, len(timeSeries['frequencies']))))
    for i in xrange(0, len(timeSeries['frequencies'])):
        sLabel = ""
        if timeSeries['frequencyUnits'].lower() == "mhz":
            sLabel = "%d MHz" % int(timeSeries['frequencies'][i])
        elif timeSeries['frequencyUnits'].lower() == "ghz":
            sLabel = "%.3f GHz" % timeSeries['frequencies'][i]
        timeArray = (np.array(timeSeries['times'][i]) - mjdDay) * 24.
        #plt.plot(timeSeries['times'][i], timeSeries['fluxDensities'][i], 'o-', label=sLabel)
        plt.plot(timeArray, timeSeries['fluxDensities'][i], 'o-', label=sLabel)
    if plotLegend == True:
        plt.legend()
    plt.ylabel("Flux Density (Jy)")
    if timeSeries['timeUnits'].lower() == "doy":
        plt.xlabel("DOY")
    elif timeSeries['timeUnits'].lower() == "mjd":
        plt.xlabel("Hours")
    elif timeSeries['timeUnits'].lower() == "datetime":
        plt.gcf().autofmt_xdate()
        plt.xlabel("Time")
    #print "minTime / maxTime"
    #print minTime
    #print maxTime
    if title is None and 'sourceName' in timeSeries:
        # Default to just putting the name of the source on the top.
        plt.title(timeSeries['sourceName'])
    elif title is not None:
        plt.title(title)
    plt.xlim((minTime, maxTime))
    plt.ylim((minFlux, maxFlux))
    plt.savefig(outputName)
    plt.close()
    return outputName

def acfPlot(acfResults, idx=[], outputName='test_acfPlot.png', plotErrors=False, title=None,
            separatePlots=False):
    # Make a plot of the autocorrelation spectrum, and optionally overlay the
    # calculated timescale Gaussian and values.
    
    # First thing to do is plot the autocorrelation spectra as a function of lap.
    if 'cor' not in acfResults:
        # Something is wrong with the thing we were passed.
        return None
    plotMade = False
    plots = []
    plotNames = []
    # Plot the requested autocorrelation spectra.
    for i in xrange(0, len(acfResults['cor'])):
        if len(idx) == 0 or i in idx:
            if separatePlots == True:
                plt.clf()
            pd = plt.errorbar(acfResults['lag'][i], acfResults['cor'][i], fmt='o-', yerr=acfResults['corError'][i],
                              xerr=acfResults['lagError'][i], label='%d %s' % (int(acfResults['frequencies'][i]),
                                                                               acfResults['frequencyUnits']))
            plotMade = True
            plots.append(pd[0].get_color())
        else:
            plots.append(None)

        if plotMade == True and ((separatePlots == True) or (i == len(acfResults['cor']) - 1)):
            # Then plot the timescale plots if they exist. First, grab the x plot range.
            xrge = plt.xlim()
            yrge = plt.ylim()
            if 'timescale' in acfResults and len(acfResults['timescale']) == len(acfResults['cor']):
                # Make a smooth enough plot of the Gaussian.
                xg = np.arange(xrge[0], xrge[1], 0.1)
                for j in xrange(0, len(acfResults['cor'])):
                    if separatePlots == True:
                        j = i
                    if acfResults['timescale'][j]['randomInterval'] is None:
                        if separatePlots == True:
                            break
                        else:
                            next
                    # Shade in the random confidence interval.
                    plt.axhspan(acfResults['timescale'][j]['randomInterval'][0],
                                acfResults['timescale'][j]['randomInterval'][1], facecolor='r', alpha=0.1)
                    vg = userGauss(xg, 1.0, 0.0, acfResults['timescale'][j]['sigma'])
                    vgp = userGauss(xg, 1.0, 0.0, (acfResults['timescale'][j]['sigma'] +
                                                   acfResults['timescale'][j]['sigmaError']))
                    vgm = userGauss(xg, 1.0, 0.0, (acfResults['timescale'][j]['sigma'] -
                                                   acfResults['timescale'][j]['sigmaError']))
                    # Plot the Gaussian line.
                    plt.plot(xg, vg, '-.', color=plots[j])
                    # We plot the error lines if requested.
                    if plotErrors == True:
                        plt.plot(xg, vgp, ':', color=plots[j])
                        plt.plot(xg, vgm, ':', color=plots[j])
                    # And then the line the user cares about.
                    ypoints = []
                    xpoints = []
                    if acfResults['timescale'][j]['mode'] == "fwhm":
                        ypoints = [0.5, 0.5]
                        xpoints = [-1 * acfResults['timescale'][j]['value'] / 2.,
                                    acfResults['timescale'][j]['value'] / 2.]
                    elif acfResults['timescale'][j]['mode'] == 'hwhm':
                        ypoints = [0.5, 0.5]
                        xpoints = [0, acfResults['timescale'][j]['value']]
                    elif acfResults['timescale'][j]['mode'] == 'fwhme':
                        ypoints = [ 1. / math.exp(1.), 1. / math.exp(1.) ]
                        xpoints = [-1 * acfResults['timescale'][j]['value'] / 2.,
                                    acfResults['timescale'][j]['value'] / 2.]
                    elif acfResults['timescale'][j]['mode'] == 'hwhme':
                        ypoints = [ 1. / math.exp(1.), 1. / math.exp(1.) ]
                        xpoints = [0, acfResults['timescale'][j]['value'] ]
                    plt.plot(xpoints, ypoints, '--', color=plots[j])
                    # Plot the range used for the Gaussian fitting.
                    xpoints = [ acfResults['timescale'][j]['fitRegion'][0],
                                acfResults['timescale'][j]['fitRegion'][0] ]
                    ypoints = [ yrge[0], yrge[1] ]
                    plt.plot(xpoints, ypoints, '-', color=plots[j])
                    xpoints = [ acfResults['timescale'][j]['fitRegion'][1],
                                acfResults['timescale'][j]['fitRegion'][1] ]
                    ypoints = [ yrge[0], yrge[1] ]
                    plt.plot(xpoints, ypoints, '-', color=plots[j])
                    if separatePlots == True:
                        break
                # Put the labels on.
                plt.xlabel("Lag [%s]" % acfResults['timeUnits'])
                plt.ylabel("ACF Strength")
                plt.legend()
                if title is None and 'sourceName' in acfResults:
                    # Default to just putting the name of the source on the top.
                    plt.title(acfResults['sourceName'])
                elif title is not None:
                    plt.title(title)
                if (separatePlots == True):
                    sOutputName = outputName.replace(".png", "_f%d.png" % int(acfResults['frequencies'][i]))
                    plt.savefig(sOutputName)
                    plotNames.append({ 'frequency': int(acfResults['frequencies'][i]),
                                       'fileName': sOutputName })
                else:
                    plt.savefig(outputName)
                    plotNames.append({ 'frequency': 0, 'fileName': outputName })
                plt.close()
    return plotNames
        
def epochTimescalePlot(epochData=None, outputName='test_epochts.png', title=None,
                       maxTimescale=None):
    # Plot the timescales measured for a single epoch as a function of frequency.
    if epochData is None:
        return

    frequencies = np.array(epochData['acf']['frequencies'])
    timescaleValues = np.array([ ((aa['value'] / 60.) if aa['value'] is not None else 0.) for aa in epochData['acf']['timescale'] ])
    timescaleValueErrors = np.array([ ((aa['valueError'] / 60.) if aa['valueError'] is not None else 0. ) for aa in epochData['acf']['timescale'] ])
    timescaleLimit = np.array([ aa['lowerLimit'] for aa in epochData['acf']['timescale'] ])
    timescaleInterval = np.array([ (aa['lagInterval'][1] - aa['lagInterval'][0]) for aa in epochData['acf']['timescale'] ])
    
    rIdx = np.where(timescaleLimit==True)
    timescaleValueErrors[rIdx] = 2.
    timescaleValues[rIdx] = timescaleInterval[rIdx]
    
    plt.errorbar(frequencies, timescaleValues, yerr=timescaleValueErrors, fmt='o', lolims=timescaleLimit)
    if title is not None:
        plt.title(title)
    plt.xlabel("Frequency [%s]" % epochData['acf']['frequencyUnits'])
    plt.ylabel("Timescale [hours]")
    if maxTimescale is not None:
        plt.ylim(0, maxTimescale)
    else:
        plt.ylim(ymin=0)
    plt.savefig(outputName)
    plt.close()
    
    return outputName

def timescaleVariationPlot(allEpochData=None, outputName='test_timescalevar.png', title=None,
                           mode='frequency', frequency=None, maxTimescale=None):
    # Plot how the timescale varies over many epochs.
    if allEpochData is None:
        return

    # Assemble the time and timescale arrays that we'll plot.
    mjds = []
    doys = []
    timescales = []
    timescaleErrors = []
    timescaleLimit = []
    
    if mode == "frequency" and frequency is not None:
        # We'll plot a single frequency's timescale variation.
        for i in xrange(0, len(allEpochData)):
            ed = allEpochData[i]
            for j in xrange(0, len(ed['acf']['frequencies'])):
                if ed['acf']['frequencies'][j] == frequency:
                    # This is the frequency we want.
                    if ed['acf']['timescale'][j]['value'] is not None:
                        timescales.append(ed['acf']['timescale'][j]['value'] / 60.)
                        timescaleErrors.append(ed['acf']['timescale'][j]['valueError'] / 60.)
                        timescaleLimit.append(ed['acf']['timescale'][j]['lowerLimit'])
                    else:
                        timescales.append(0.)
                        timescaleErrors.append(0.)
                        timescaleLimit.append(True)
                    doys.append(ed['timeTuple'][3])
                    
        ndoys = np.array(doys)
        ntimescales = np.array(timescales)
        ntimescaleErrors = np.array(timescaleErrors)
        ntimescaleLimit = np.array(timescaleLimit)
        rIdx = np.where(ntimescaleLimit==True)
        ntimescaleErrors[rIdx] = 2.
        plt.errorbar(ndoys, ntimescales, yerr=ntimescaleErrors, 
                     fmt='o', lolims=ntimescaleLimit)
        if title is not None:
            plt.title(title)
        plt.xlabel("Day of Year")
        plt.ylabel("Timescale [hours]")
        plt.xlim(0, 365)
        if maxTimescale is not None:
            plt.ylim(0, maxTimescale)
        else:
            plt.ylim(ymin=0)
        plt.savefig(outputName)
        plt.close()

    return outputName
