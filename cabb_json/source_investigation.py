import cabb_json as ihv
import sys
import os

# This script looks for all files relating to the named source, then tries
# to form several things:
# 1. The time series over the entire available time range.
# 2. Dynamic spectra for each individual epoch.
# 3. Timescale calculations for all data as a function of DOY.

# The argument is the name of the source to look for.
sourceName = sys.argv[1]
# The second optional argument is the directory under which all data is found.
dataDirectory = "."
if len(sys.argv) == 3:
    dataDirectory = sys.argv[2]

print "Searching for data for source %s" % sourceName
sourceFiles = []
for root, dirs, files in os.walk(dataDirectory, followlinks=True):
    for file in files:
        if file.endswith(".json") and file.startswith(sourceName):
            sourceFiles.append(os.path.join(root, file))
print "  found %d files" % len(sourceFiles)

# Read in the data for each of these files.
epochData = []
for i in xrange(0, len(sourceFiles)):
    data = ihv.readJson(sourceFiles[i])
    # Collect some metadata for ease of use.
    metadata = { 'data': data }
    metadata['nTimeIntervals'] = len(data.timeSeries['I'].measurements)
    metadata['mjdLow'] = data.timeSeries['I'].measurements[0].mjd['low']
    metadata['mjdHigh'] = data.timeSeries['I'].measurements[-1].mjd['high']
    metadata['timeInterval'] = (metadata['mjdHigh'] - metadata['mjdLow']) * 1440.
    bandsPresent = {}
    for j in xrange(0, len(data.timeSeries['I'].measurements)):
        band = ihv.frequency2Band(data.timeSeries['I'].measurements[j].bandRanges[0][0])
        bandsPresent[band] = True
    metadata['bandsPresent'] = bandsPresent.keys()
    epochData.append(metadata)

# Sort the data into ascending time order.
epochData = sorted(epochData, key=lambda data: data['mjdLow'])
for i in xrange(0, len(epochData)):
    # Print out information about each file.
    print "FILE %d INFORMATION:" % (i + 1)
    print " Source Name (RA, Dec): %s (%s, %s)" % (epochData[i]['data'].name,
                                                   epochData[i]['data'].rightAscension,
                                                   epochData[i]['data'].declination)
    print "      # time intervals: %d" % epochData[i]['nTimeIntervals']
    print "  MJD range low - high: %.4f - %.4f" % ( epochData[i]['mjdLow'],
                                                    epochData[i]['mjdHigh'] )
    print "         time interval: %.2f minutes" % epochData[i]['timeInterval']
    print "         bands present: %s" % ", ".join(epochData[i]['bandsPresent'])

# Some settings we need to use later.
frequencyUnits = "MHz"
spectralAveraging = 1
bandFrequencies = [ 4800, 8400, 17400 ]

# Start by collating a complete time-series and plotting it.
completeTimeSeries = {}
for i in xrange(0, len(epochData)):
    print "  GETTING TIME SERIES FOR FILE %d" % (i + 1)
    timeSeries = epochData[i]['data'].getTimeSeries({ 'spectralAveraging': spectralAveraging,
                                                      'frequencyUnits': frequencyUnits,
                                                      'frequencies': bandFrequencies,
                                                      'alwaysPresent': False,
                                                      'timeUnits': 'mjd' })
    print "  STORING TIME SERIES"
    if 'timeUnits' not in completeTimeSeries:
        # Start storing the time series.
        completeTimeSeries = timeSeries
    else:
        # Append this new data to the old one.
        for j in xrange(0, len(timeSeries['frequencies'])):
            pos = -1
            for k in xrange(0, len(completeTimeSeries['frequencies'])):
                if completeTimeSeries['frequencies'][k] == timeSeries['frequencies'][j]:
                    pos = k
                    break
            if pos > -1:
                # Append this new information.
                completeTimeSeries['times'][k] += timeSeries['times'][j]
                completeTimeSeries['fluxDensities'][k] += timeSeries['fluxDensities'][j]
            else:
                # Make a new frequency object.
                completeTimeSeries['frequencies'].append(timeSeries['frequencies'][j])
                completeTimeSeries['times'].append(timeSeries['times'][j])
                completeTimeSeries['frequencies'].append(timeSeries['frequencies'][j])
# We can now plot the complete time series.
ihv.timeSeriesPlot(completeTimeSeries, outputName='full_%s_timeSeries.png' % sourceName)
