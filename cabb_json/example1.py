import cabb_json as ihv
import sys
import matplotlib.pyplot as plt
import numpy as np

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)

print "There are %d time intervals" % (len(source.timeSeries['I'].measurements))
# Get all the spectra, averaged to 16 MHz resolution.
allSpectra = source.getSpectra({ 'splitBand': True, 'spectralAveraging': 0.016,
                                 'frequencyUnits': "GHz", 'timeUnits': "DOY" })

# Get a time series at two frequencies.
timeSeries = source.getTimeSeries({ 'spectralAveraging': 16, 'frequencyUnits': "MHz",
                                    'frequencies': [ 4800, 8400, 17400 ], 'alwaysPresent': False,
                                    'timeUnits': 'mjd' })
print timeSeries
ihv.timeSeriesPlot(timeSeries, outputName='test_spectralplot.png')

ihv.spectraPlot(allSpectra, outputName='test_gistncar.png', includeZero=True)
ihv.spectraPlot(allSpectra, outputName='test_animation_%04d.png', plotType='animation')

acfResults = ihv.calculateACF(timeSeries)
#for i in xrange(0, len(acfResults['cor'])):
#    print acfResults['frequencies'][i]
#    print acfResults['lag'][i]
#    print acfResults['cor'][i]
#    print acfResults['corError'][i]
#    plt.errorbar(acfResults['lag'][i], acfResults['cor'][i], fmt='o-', yerr=acfResults['corError'][i],
#                 xerr=acfResults['lagError'][i],
#                 label="%d MHz" % int(acfResults['frequencies'][i]))
#    timescaleResults = ihv.calculateTimescale(acfResults['lag'][i], acfResults['cor'][i])
#    print "the %s of the time scale is %f %s" % (timescaleResults['mode'], timescaleResults['value'],
#                                                 timescaleResults['timeUnits'])
#    
#    
#    plt.plot(acfResults['lag'][i],
#             ihv.userGauss(np.array(acfResults['lag'][i]), 1.0, 0.0, timescaleResults['sigma']), 'ro:', label='fit')
#plt.legend()
#plt.ylim((0., 1.2))
#plt.savefig('test_acfplot.png')
#plt.close()

acfResults['timescale'] = []
for i in xrange(0, len(acfResults['cor'])):
    timescaleResults = ihv.calculateTimescale(acfResults['lag'][i], acfResults['cor'][i], mode='fwhme')
    print "time scale is %f %s" % (timescaleResults['value'], timescaleResults['timeUnits'])
    acfResults['timescale'].append(timescaleResults)
ihv.acfPlot(acfResults)
