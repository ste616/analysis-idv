import cabb_json as ihv
import sys

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)

print "There are %d time intervals" % (len(source.timeSeries['I'].measurements))
for i in xrange(0, len(source.timeSeries['I'].measurements)):
    a = source.timeSeries['I'].measurements[i]
    print a.centreFrequencies
    print a.bandRanges
    print a.channelWidth
    print "=="
                                                                   
