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

epochData = sorted(epochData, key=lambda data: data['mjdLow'])
for i in xrange(0, len(epochData)):
    print "FILE %d INFORMATION:" % (i + 1)
    print " Source Name (RA, Dec): %s (%s, %s)" % (epochData[i]['data'].name,
                                                   epochData[i]['data'].rightAscension,
                                                   epochData[i]['data'].declination)
    print "      # time intervals: %d" % epochData[i]['nTimeIntervals']
    print "  MJD range low - high: %.4f - %.4f" % ( epochData[i]['mjdLow'],
                                                    epochData[i]['mjdHigh'] )
    print "         time interval: %.2f minutes" % epochData[i]['timeInterval']
    print "         bands present: %s" % ", ".join(epochData[i]['bandsPresent'])
