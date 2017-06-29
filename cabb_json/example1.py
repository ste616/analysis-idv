import cabb_json as ihv
import sys

# Read in a single file, which is specified on the command line.
source = ihv.readJson(sys.argv[1])

print "Read in file with source %s (%s, %s)" % (source.name, source.rightAscension, source.declination)
