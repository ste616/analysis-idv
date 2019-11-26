#!/usr/bin/env python
"""Time-variable analysis for ATCA

Usage:
  source_investigation.py <sourcename> [--plot-dir=<dir>] [--root-dir=<dir>] [--web-dir=<dir>] [--spectral-averaging=<width>] [--frequency=<freq>]... [--low-frequency=<freq>] [--high-frequency=<freq>] [--frequency-interval=<intv>]

-h --help                      show this
-p --plot-dir DIR              output all the plots and other outputs in the specified DIR [default: .]
-r --root-dir DIR              find all the input data files below this DIR [default: .]
-w --web-dir DIR               output all the web page data in the specified DIR [default: .]
-A --spectral-averaging WIDTH  how many MHz to average into each data point [default: 1.0]
-f --frequency FREQ            include this frequency (MHZ) in the plots (no default)
-L --low-frequency FREQ        lowest frequency (MHz) to plot [default: 4500]
-H --high-frequency FREQ       highest-frequency (MHz) to plot [default: 10600]
-I --frequency-interval INTV   the interval (MHz) between each plotted frequency [default: 100]
"""

from docopt import docopt
import matplotlib
matplotlib.use("Agg")

import analysis_idv as ihv
import sys
import os
import math
from astropy.time import Time
import json
import copy
import shutil

# This script looks for all files relating to the named source, then tries
# to form several things:
# 1. The time series over the entire available time range.
# 2. Dynamic spectra for each individual epoch.
# 3. Timescale calculations for all data as a function of DOY.

def find_source_files(source_name, data_directory):
    sourceFiles = []
    for root, dirs, files in os.walk(data_directory, followlinks=True):
        for file in files:
            if file.endswith(".json") and file.startswith(source_name):
                sourceFiles.append(os.path.join(root, file))
                #print "    found new file %s" % os.path.join(root, file)
    return sourceFiles

def read_source_file(fname):
    print "  reading file %s" % (fname)
    data = ihv.readJson(fname)
    # Collect some metadata for ease of use.
    metadata = { 'data': data, 'useable': True }
    metadata['nTimeIntervals'] = len(data.timeSeries['I'].measurements)
    try:
        metadata['mjdLow'] = data.timeSeries['I'].measurements[0].mjd['low']
        metadata['mjdHigh'] = data.timeSeries['I'].measurements[-1].mjd['high']
        metadata['mjdMid'] = metadata['mjdLow'] + (metadata['mjdHigh'] -
                                                   metadata['mjdLow']) / 2.
        metadata['timeInterval'] = (metadata['mjdHigh'] - metadata['mjdLow']) * 1440.
        epochTime = Time(metadata['mjdLow'], format='mjd')
        metadata['timeLow'] = Time(metadata['mjdLow'], format='mjd')
        metadata['timeHigh'] = Time(metadata['mjdHigh'], format='mjd')
        metadata['timeMid'] = Time(metadata['mjdMid'], format='mjd')
        bandsPresent = {}
        for j in xrange(0, len(data.timeSeries['I'].measurements)):
            band = ihv.frequency2Band(data.timeSeries['I'].measurements[j].bandRanges[0][0])
            bandsPresent[band] = True
        metadata['bandsPresent'] = bandsPresent.keys()
        epochTime.format = "iso"
        epochTime.out_subfmt = "date"
        metadata['epochName'] = epochTime.value
    except IndexError:
        # No data.
        print " encountered index error!"
        sys.exit(0)
        metadata['useable'] = False
        metadata['mjdLow'] = 0.
    return metadata

def analyse(args):
    print "Searching for data for source %s" % args['<sourcename>']
    source_files = find_source_files(args['<sourcename>'], args['--root-dir'])
    print "  found %d files" % len(source_files)

    # We'll collect data on the flux density ranges in a file that
    # we can continually append to.
    fdFile = "%s/fluxdensities.csv" % (args['--plot-dir'])
    fdOutputLines = []

    # Read in the data for each of these files.
    epochData = []
    print "READING FILES"
    for i in xrange(0, len(source_files)):
        fdata = read_source_file(source_files[i])
        if fdata['useable'] == True:
            epochData.append(fdata)
        else:
            print "FILE %d is not useable" % i

    # Sort the data into ascending time order.
    epochData = sorted(epochData, key=lambda data: data['mjdLow'])
    # Remove this if it doesn't have enough data.
    repochs = []
    for i in xrange(0, len(epochData)):
        if (epochData[i]['nTimeIntervals'] < 3):
            repochs.append(i)

    for i in xrange(0, len(repochs)):
        print "  Removing epoch %d" % (repochs[i] - i)
        del epochData[repochs[i] - i]

    # Now exit if we have nothing.
    print "  %d EPOCHS REMAINING" % (len(epochData))
    if (len(epochData) == 0):
        sys.exit()

    allBands = []
    for i in xrange(0, len(epochData)):
        if epochData[i]['useable'] == False:
            continue
        # Print out information about each file.
        print "FILE %d INFORMATION:" % (i + 1)
        print " Source Name (RA, Dec): %s (%s, %s)" % (
            epochData[i]['data'].name,
            epochData[i]['data'].rightAscension,
            epochData[i]['data'].declination)
        print "      # time intervals: %d" % epochData[i]['nTimeIntervals']
        print "  MJD range low - high: %.4f - %.4f" % ( epochData[i]['mjdLow'],
                                                        epochData[i]['mjdHigh'] )
        print "         time interval: %.2f minutes" % epochData[i]['timeInterval']
        print "         bands present: %s" % ", ".join(epochData[i]['bandsPresent'])
        for j in xrange(0, len(epochData[i]['bandsPresent'])):
            if epochData[i]['bandsPresent'][j] not in allBands:
                allBands.append(epochData[i]['bandsPresent'][j])
        fdOutputLines.append([ args['<sourcename>'],
                               epochData[i]['data'].rightAscension,
                               epochData[i]['data'].declination,
                               epochData[i]['nTimeIntervals'],
                               epochData[i]['mjdLow'], epochData[i]['mjdHigh'],
                               epochData[i]['timeInterval'] ])

    # Some settings we need to use later.
    frequencyUnits = "MHz"
    spectralAveraging = float(args['--spectral-averaging'])
    if (len(args['--frequency']) > 0):
        bandFrequencies = []
        for i in xrange(0, len(args['--frequency'])):
            bandFrequencies.append(int(args['--frequency'][i]))
    else:
        bandFrequencies = range(int(args['--low-frequency']),
                                int(args['--high-frequency']),
                                int(args['--frequency-interval']))
    maxDiff = spectralAveraging

    print ""
    print "STAGE 1: COLLATE TIME SERIES"

    # Start by collating a complete time-series and plotting it.
    completeTimeSeries = {}
    for i in xrange(0, len(epochData)):
        if epochData[i]['useable'] == False:
            continue
        print "  GETTING TIME SERIES FOR FILE %d" % (i + 1)
        timeSeries = epochData[i]['data'].getTimeSeries({
            'spectralAveraging': spectralAveraging,
            'frequencyUnits': frequencyUnits,
            'frequencies': bandFrequencies,
            'maxDiff': maxDiff,
            'alwaysPresent': False,
            'timeUnits': 'mjd' })
        if len(timeSeries['frequencies']) == 0:
            print "  DISCARDING TIME SERIES"
            continue
        print "  STORING TIME SERIES"
        epochData[i]['timeSeries'] = timeSeries
        if 'timeUnits' not in completeTimeSeries:
            # Start storing the time series.
            completeTimeSeries = copy.deepcopy(timeSeries)
        else:
            # Append this new data to the old one.
            for j in xrange(0, len(timeSeries['frequencies'])):
                pos = -1
                #print "matching frequency %.1f" % timeSeries['frequencies'][j]
                for k in xrange(0, len(completeTimeSeries['frequencies'])):
                    #print "  looking at frequency %.1f" % completeTimeSeries['frequencies'][k]
                    if (completeTimeSeries['frequencies'][k] ==
                        timeSeries['frequencies'][j]):
                        #print "  found match at %d" % k
                        pos = k
                        break
                #print "pos now %d" % pos
                if pos > -1:
                    # Append this new information.
                    completeTimeSeries['times'][pos] += timeSeries['times'][j]
                    completeTimeSeries['fluxDensities'][pos] += timeSeries['fluxDensities'][j]
                else:
                    # Make a new frequency object.
                    completeTimeSeries['frequencies'].append(
                        timeSeries['frequencies'][j])
                    completeTimeSeries['times'].append(timeSeries['times'][j])
                    completeTimeSeries['fluxDensities'].append(
                        timeSeries['fluxDensities'][j])

        # Make a plot of the individual day's time series.
        roughMJD = math.floor(epochData[i]['mjdLow'])
        tt = epochData[i]['timeLow'].datetime.timetuple()
        year = tt.tm_year
        roughDOY = tt.tm_yday
        epochData[i]['timeTuple'] = ( args['<sourcename>'], roughMJD,
                                      year, roughDOY )
        out_plot = "%s/daily_%s_timeSeries_mjd%d_doy%04d-%03d.png" % (
            (args['--plot-dir'],) + epochData[i]['timeTuple'])
        epochData[i]['timeSeriesImage'] = ihv.timeSeriesPlot(
            timeSeries, plotLegend=False,
            title='%s (MJD %d DOY %04d-%03d)' % epochData[i]['timeTuple'],
            outputName=out_plot)
        modIndex = ihv.calculateModulationIndex(timeSeries)
        epochData[i]['modulationIndex'] = modIndex
        for j in xrange(0, len(modIndex['frequencies'])):
            #print fdOutputLines[i]
            fdOutputLines[i] += [ modIndex['frequencies'][j],
                                  modIndex['averageFlux'][j],
                                  modIndex['modulationIndex'][j] ]
            #print fdOutputLines[i]

    # Output the file now.
    with open(fdFile, 'a') as fd:
        for i in xrange(0, len(fdOutputLines)):
            outLine = ""
            for j in xrange(0, len(fdOutputLines[i])):
                if j > 0:
                    outLine += ","
                outLine += str(fdOutputLines[i][j])
            outLine += "\n"
            fd.write(outLine)

    # We can now plot the complete time series.
    try:
        tsplot = "%s/full_%s_timeSeries.png" % (args['--plot-dir'],
                                                args['<sourcename>'])
        ihv.timeSeriesPlot(completeTimeSeries, outputName=tsplot)
    except KeyError as e:
        print "KeyError caught"

    print ""
    print "STAGE 2: CALCULATE TIMESCALES"

    # We store our timescale results into an object that we will turn into a JSON
    # output file at the end.
    timescaleResults = { args['<sourcename>']: [] }
    # Keep track of the maximum time scale that we find.
    maxTimescale = 0.
    for i in xrange(0, len(epochData)):
        if epochData[i]['useable'] == False:
            continue
        if 'timeSeries' not in epochData[i]:
            epochData[i]['useable'] = False
            continue
        print '  CALCULATING AUTO-CORRELATION FUNCTION FOR EPOCH %d' % (i + 1)
        acf = ihv.calculateACF(epochData[i]['timeSeries'])
        acf['timescale'] = []
        epochData[i]['acf'] = acf
        for j in xrange(0, len(acf['cor'])):
            print '    CALCULATING TIMESCALE FOR FREQUENCY %d' % (
                int(acf['frequencies'][j]))
            timescale = ihv.calculateTimescale(acf['lag'][j], acf['cor'][j],
                                               acf['corError'][j][0], mode='fwhme')
            if timescale is not None and timescale['value'] is not None:
                timescaleUnits = timescale['timeUnits']
                if timescale['lowerLimit'] == False:
                    print '      timescale is %.1f +/- %.1f %s' % (
                        timescale['value'], timescale['valueError'],
                        timescale['timeUnits'])
                    timescaleValue = timescale['value']
                    timescaleError = timescale['valueError']
                    timescaleType = timescale['mode']
                else:
                    print '      timescale is > %.1f %s' % (
                        timescale['lagInterval'][1], timescale['timeUnits'])
                    timescaleValue = timescale['lagInterval'][1]
                    timescaleError = -1.
                    timescaleType = 'lower_limit'
            else:
                timescaleValue = -1
                timescaleError = -1
                timescaleType = 'error'
                timescaleUnits = 'error'
            if timescaleValue > maxTimescale:
                maxTimescale = timescaleValue
            epochData[i]['acf']['timescale'].append(timescale)
            # Put this in our timescale results object.
            timescaleResults[args['<sourcename>']].append({
                'frequency': acf['frequencies'][j],
                'frequencyUnits': acf['frequencyUnits'],
                'mjd': epochData[i]['mjdMid'],
                'year': epochData[i]['timeMid'].datetime.timetuple().tm_year,
                'doy': ihv.astro2doy(epochData[i]['timeMid']),
                'timescaleType': timescaleType,
                'timescaleValue': timescaleValue,
                'timescaleError': timescaleError,
                'timescaleUnits': timescaleUnits
                })
        try:
            acfout = "%s/daily_%s_acf_mjd%d_doy%04d-%03d.png" % (
                (args['--plot-dir'],) + epochData[i]['timeTuple'])
            acfPlots = ihv.acfPlot(
                epochData[i]['acf'], plotErrors=True, separatePlots=True,
                title='%s (MJD %d DOY %04d-%03d)' % epochData[i]['timeTuple'],
                outputName=acfout)
        except ValueError:
            acfPlots = None
        epochData[i]['acfPlotNames'] = acfPlots
    # Write out our JSON file.
    outputJsonFile = '%s/timescales_%s.json' % (
        args['--plot-dir'], args['<sourcename>'])
    with open(outputJsonFile, 'w') as oj:
        json.dump(timescaleResults, oj)

    print ""
    print "STAGE 3: PRODUCE DYNAMIC SPECTRA"

    for i in xrange(0, len(epochData)):
        if epochData[i]['useable'] == False:
            continue
        print "  GETTING SPECTRA FOR FILE %d" % (i + 1)
        spectra = epochData[i]['data'].getSpectra({
            'splitBand': True,
            'spectralAveraging': spectralAveraging,
            'frequencyUnits': frequencyUnits,
            'timeUnits': 'mjd' })
        print "  STORING SPECTRA"
        epochData[i]['spectra'] = spectra
        # Make a dynamic spectrum plot.
        dsplot = "%s/daily_%s_dynamic_mjd%d_doy%04d-%03d.png" % (
            (args['--plot-dir'],) + epochData[i]['timeTuple'])
        epochData[i]['dynamicImage'] = ihv.spectraPlot(
            spectra, outputName=dsplot,
            title='%s (MJD %d DOY %04d-%03d)' % epochData[i]['timeTuple'])
        animationDir = "%s/animation_frames" % args['--plot-dir']
        if (not os.path.isdir(animationDir)):
            os.makedirs(animationDir)
        animationName = '%s/daily_%s_mjd%d_doy%04d-%03d_animation' % (
            (animationDir,) + epochData[i]['timeTuple'])
        animationName = animationName + "_frame%d.png"
        ihv.spectraPlot(spectra, plotType="animation", outputName=animationName,
                        title='%s (MJD %d DOY %04d-%03d)' % epochData[i]['timeTuple'])

    print ""
    print "STAGE 4: PRODUCE WEB PAGES"

    # To make nice next/previous links, we need to construct the names of the
    # pages first.
    for i in xrange(0, len(epochData)):
        if epochData[i]['useable'] == False:
            continue
        epochData[i]['htmlTimescalesFile'] = "%s_doy%04d-%03d_epochinfo.html" % (
            epochData[i]['timeTuple'][0], epochData[i]['timeTuple'][2],
            epochData[i]['timeTuple'][3])

    # Make the plot for the entire time period for each frequency.
    timescalesDoyPlots = []
    outputDirectory = "%s/%s" % (args['--web-dir'], args['<sourcename>'])
    if (not os.path.isdir(outputDirectory)):
        # Make this directory.
        os.makedirs(outputDirectory)
    # Increase the maximum time scale plotted by 1 hour.
    maxTimescale /= 60.
    maxTimescale += 1.
    nf = 0
    fi = -1
    for i in xrange(0, len(epochData)):
        if 'acf' in epochData[i]:
            nf = len(epochData[i]['acf']['frequencies'])
            fi = i
            break
    if (nf == 0 or fi == -1):
        print "  WARNING: Unable to find frequencies"
        #sys.exit()
    for i in xrange(0, nf):
        f = int(epochData[fi]['acf']['frequencies'][i])
        atsplot = "%s/%s_%d_timescales.png" % (args['--plot-dir'],
                                               epochData[fi]['timeTuple'][0], f)
        allTimescalePlot = ihv.timescaleVariationPlot(
            epochData, frequency=f, maxTimescale=maxTimescale,
            outputName=atsplot)
        timescalesDoyPlots.append({ 'frequency': f, 'plotFile': allTimescalePlot })
        if os.path.isfile(allTimescalePlot):
            shutil.move(allTimescalePlot,
                        "%s/%s" % (outputDirectory,
                                   os.path.basename(allTimescalePlot)))

    # Make the index file.
    if (nf > 0):
        ihv.outputIndex(epochData, outputDirectory, "index.html", timescalesDoyPlots)

    for i in xrange(0, len(epochData)):
        if (not os.path.isdir(outputDirectory)):
            # Make this directory.
            os.makedirs(outputDirectory)
        if ("dynamicImage" in epochData[i] and
            os.path.isfile(epochData[i]['dynamicImage'])):
            # Move this image into the directory.
            shutil.move(epochData[i]['dynamicImage'], "%s/%s" % (
                outputDirectory,
                os.path.basename(epochData[i]['dynamicImage'])))
        if ("timeSeriesImage" in epochData[i] and
            os.path.isfile(epochData[i]['timeSeriesImage'])):
            shutil.move(epochData[i]['timeSeriesImage'], "%s/%s" % (
                outputDirectory,
                os.path.basename(epochData[i]['timeSeriesImage'])))
        if "acfPlotNames" in epochData[i]:
            for j in xrange(0, len(epochData[i]['acfPlotNames'])):
                if (epochData[i]['acfPlotNames'] is not None and
                    os.path.isfile(epochData[i]['acfPlotNames'][j]['fileName'])):
                    shutil.move(epochData[i]['acfPlotNames'][j]['fileName'], "%s/%s" % (
                        outputDirectory,
                        os.path.basename(epochData[i]['acfPlotNames'][j]['fileName'])))
        linkNext = None
        if i < (len(epochData) - 1):
            linkNext = epochData[i + 1]['htmlTimescalesFile']
        linkPrevious = None
        if i > 0:
            linkPrevious = epochData[i - 1]['htmlTimescalesFile']
        linkUp = "index.html"
        # Make a plot of all the timescales.
        tspname = "%s/%s_doy%04d-%03d_timescaleFrequency.png" % (
            args['--plot-dir'], epochData[i]['timeTuple'][0],
            epochData[i]['timeTuple'][2], epochData[i]['timeTuple'][3])
        tsptitle = "%s (DOY %04d-%03d)" % (epochData[i]['timeTuple'][0],
                                           epochData[i]['timeTuple'][2],
                                           epochData[i]['timeTuple'][3])
        tscalePlot = ihv.epochTimescalePlot(
            epochData[i], maxTimescale=maxTimescale,
            outputName=tspname, title=tsptitle)
        if os.path.isfile(tscalePlot):
            shutil.move(tscalePlot, "%s/%s" % (outputDirectory,
                                               os.path.basename(tscalePlot)))
            epochData[i]['timescaleFrequencyImage'] = tscalePlot

        ihv.outputEpoch(epochData=epochData[i], outputDirectory=outputDirectory, 
                        outputName=epochData[i]['htmlTimescalesFile'],
                        linkNext=linkNext,
                        linkPrevious=linkPrevious, linkUp=linkUp)


if __name__ == '__main__':
    arguments = docopt(__doc__, version="Time-variable analyser 1.0")
    valid = True
    if (('<sourcename>' not in arguments) or
        (arguments['<sourcename>'] == "")):
        print "No source specified!"
        valid = False
    if (not os.path.isdir(arguments['--plot-dir'])):
        print "Output directory not found, making it."
        os.makedirs(arguments['--plot-dir'])
        if (not os.path.isdir(arguments['--plot-dir'])):
            print "Output directory unable to be made!"
            valid = False
    if (not os.path.isdir(arguments['--web-dir'])):
        print "Web output directory not found, making it."
        os.makedirs(arguments['--web-dir'])
        if (not os.path.isdir(arguments['--web-dir'])):
            print "Web output directory unable to be made!"
            valid = False
    if (not os.path.isdir(arguments['--root-dir'])):
        print "Root directory not found."
        valid = False
    if (valid == True):
        analyse(arguments)
