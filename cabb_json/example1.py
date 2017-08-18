import cabb_json as ihv
import sys
import matplotlib.pyplot as plt
import numpy as np

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)

print "There are %d time intervals" % (len(source.timeSeries['I'].measurements))
# Get all the spectra, averaged to 16 MHz resolution.
allSpectra = source.getSpectra({ 'splitBand': True, 'spectralAveraging': 0.064,
                                 'frequencyUnits': "GHz", 'timeUnits': "mjd" })

# Get a time series at two frequencies.
timeSeries = source.getTimeSeries({ 'spectralAveraging': 64, 'frequencyUnits': "MHz",
                                    'frequencies': [ 5700, 7000, 17400 ], 'alwaysPresent': False,
                                    'timeUnits': 'mjd' })
print timeSeries
ihv.timeSeriesPlot(timeSeries, outputName='test_spectralplot.png')

ihv.spectraPlot(allSpectra, outputName='test_gistncar.png', includeZero=True)
ihv.spectraPlot(allSpectra, outputName='test_animation_%04d.png', plotType='animation')

acfResults = ihv.calculateACF(timeSeries)

acfResults['timescale'] = []
for i in xrange(0, len(acfResults['cor'])):
    timescaleResults = ihv.calculateTimescale(acfResults['lag'][i], acfResults['cor'][i], acfResults['corError'][i][0],
                                              mode='hwhm')
    print "time scale is %.3f +/- %.3f %s" % (timescaleResults['value'], timescaleResults['valueError'],
                                              timescaleResults['timeUnits'])
    acfResults['timescale'].append(timescaleResults)
ihv.acfPlot(acfResults, plotErrors=True)
