import cabb_json as ihv
import sys
import matplotlib.pyplot as plt

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)

print "There are %d time intervals" % (len(source.timeSeries['I'].measurements))
for i in xrange(0, len(source.timeSeries['I'].measurements)):
    a = source.timeSeries['I'].measurements[i]

# Get all the spectra, averaged to 4 MHz resolution.
allSpectra = source.getSpectra({ 'splitBand': True, 'spectralAveraging': 16.0,
                                 'frequencyUnits': "GHz" })
for i in xrange(0, len(allSpectra['spectra'])):
    s = allSpectra['spectra'][i]
    for j in xrange(0, len(s['freq'])):
        plt.plot(s['freq'][j], s['amp'][j])
plt.savefig('test.png')
plt.close()

# Get a time series at two frequencies.
timeSeries = source.getTimeSeries({ 'spectralAveraging': 16.0, 'frequencyUnits': "GHz",
                                    'frequencies': [ 4.8, 8.4 ], 'alwaysPresent': False })
print timeSeries
for i in xrange(0, len(timeSeries['frequencies'])):
    plt.plot(timeSeries['times'][i], timeSeries['fluxDensities'][i], 'o-', label="%.3f GHz" % timeSeries['frequencies'][i])
plt.legend()
plt.savefig('test_ts.png')
plt.close()

