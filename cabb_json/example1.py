import cabb_json as ihv
import sys
import matplotlib.pyplot as plt

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)

print "There are %d time intervals" % (len(source.timeSeries['I'].measurements))
# Get all the spectra, averaged to 16 MHz resolution.
allSpectra = source.getSpectra({ 'splitBand': True, 'spectralAveraging': 0.016,
                                 'frequencyUnits': "GHz" })

# Get a time series at two frequencies.
timeSeries = source.getTimeSeries({ 'spectralAveraging': 16, 'frequencyUnits': "MHz",
                                    'frequencies': [ 4800, 8400 ], 'alwaysPresent': False })
for i in xrange(0, len(timeSeries['frequencies'])):
    plt.plot(timeSeries['times'][i], timeSeries['fluxDensities'][i], 'o-', label="%.3f GHz" % timeSeries['frequencies'][i])
plt.legend()
plt.savefig('test_ts.png')
plt.close()

ihv.spectraPlot(allSpectra, outputName='test_gistncar.png', includeZero=True)
ihv.spectraPlot(allSpectra, outputName='test_animation_%04d.png', plotType='animation')
