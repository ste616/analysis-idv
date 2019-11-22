#!/usr/bin/python
"""Time-variable calibration for ATCA

Usage:
  timevariation_cal.py <dataset> [--calibrator=<cal>] [--segment=<min>] [--slop=<min>]

-h --help             show this
-c --calibrator CAL   the name of the calibrator [default: 1934-638]
-s --segment LEN      the length of each calibration segment, in minutes [default: 2]
-S --slop SLOP        a tolerance on the segment length, if it makes sense to extend it, in minutes [default: 1]
"""

from docopt import docopt
from mirpy import miriad
import ephem
from datetime import date, datetime, timedelta
import fnmatch
import os
import shutil
import numpy as np
import sys
import json

# A routine to turn a Miriad type time string into a ephem Date.
def mirtime_to_date(mt):
    year = 2000 + int(mt[0:2])
    monthShort = mt[2:5]
    date = int(mt[5:7])
    hour = int(mt[8:10])
    minute = int(mt[11:13])
    second = int(round(float(mt[14:18])))
    monthDict = { 'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
                  'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12 }
    dateString = "%4d/%02d/%02d %02d:%02d:%02d" % (year, monthDict[monthShort], date,
                                                   hour, minute, second)
    return ephem.Date(dateString)

def datetime_to_mirtime(dt):
    # Output a Miriad formatted date.
    rs = dt.strftime("%y%b%d:%H:%M:%S").lower()
    return rs

# Use a uvlist log to get the cycle time.
def filter_uvlist_variables(logfile_name):
    # We send back a dictionary.
    rd = { 'cycle_time': -1. }
    with open(logfile_name, "r") as fp:
        loglines = fp.readlines()
    for i in xrange(0, len(loglines)):
        index_elements = loglines[i].split()
        if ((len(index_elements) > 2) and
            (index_elements[0] == "inttime") and (index_elements[1] == ":")):
            rd['cycle_time'] = float(index_elements[2])
    return rd

# Use uvindex to work out the necessary parameters of this dataset.
def filter_uvindex(output):
    # We send back a dictionary.
    rd = { 'index': { 'time': [], 'source': [], 'calcode': [], 'antennas': [],
                      'spectral_channels': [], 'wideband_channels': [], 'freq_config': [],
                      'record_number': [] },
           'total_time': 0,
           'freq_configs': [], 'polarisations': [], 'sources': [] }
    outlines = output.split('\n')
    section = 0
    freqconfig_n = 0
    freqconfig_found = 0
    fc = None
    sourcearea = 0
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if ((section == 0) and (len(outlines[i]) >= 74)):
            indexTime = mirtime_to_date(outlines[i][0:18])
            if ((index_elements[1] != "Total") and (index_elements[2] != "number")):
                # This is a regular line.
                offset = 0
                rd['index']['time'].append(indexTime)
                rd['index']['source'].append(index_elements[1])
                # Check if we have a calibrator code.
                calcode = outlines[i][36:37]
                if (calcode == " "):
                    # No calibrator code.
                    offset = 1
                rd['index']['calcode'].append(calcode)
                rd['index']['antennas'].append(int(index_elements[3 - offset]))
                rd['index']['spectral_channels'].append(int(index_elements[4 - offset]))
                rd['index']['wideband_channels'].append(int(index_elements[5 - offset]))
                rd['index']['freq_config'].append(int(index_elements[6 - offset]))
                rd['index']['record_number'].append(int(index_elements[7 - offset]))
            else:
                # We've moved to the next section
                section = 1
        elif ((section == 1) and (len(index_elements) > 0) and (index_elements[0] == "Total") and
              (index_elements[1] == "observing")):
            # We've found the total amount of observing time.
            rd['total_time'] = float(index_elements[4])
            section = 2
        elif (section == 2):
            if ((len(index_elements) > 0) and (index_elements[0] == "Frequency")
                and (index_elements[1] == "Configuration")):
                freqconfig_n = int(index_elements[2])
                freqconfig_found = 1
                if (fc is not None):
                    rd['freq_configs'].append(fc)
                fc = { 'number': freqconfig_n, 'nchannels': [],
                       'frequency1': [], 'frequency_increment': [],
                       'rest_frequency': [], 'ifchain': [] }
            elif (freqconfig_found == 1):
                freqconfig_found = 2
            elif (freqconfig_found == 2):
                if (outlines[i] == ""):
                    freqconfig_found = 0
                else:
                    # This is the actual line.
                    fc['nchannels'].append(int(index_elements[0]))
                    fc['frequency1'].append(float(index_elements[1]))
                    fc['frequency_increment'].append(float(index_elements[2]))
                    fc['rest_frequency'].append(float(index_elements[3]))
                    fc['ifchain'].append(int(index_elements[5]))
            elif (outlines[i] == "------------------------------------------------"):
                if (fc is not None):
                    rd['freq_configs'].append(fc)
                section = 3
        elif (section == 3):
            if ((len(index_elements) > 0) and (index_elements[0] == "There") and
                (index_elements[3] == "records") and (index_elements[5] == "polarization")):
                rd['polarisations'].append(index_elements[6])
            elif (outlines[i] == "------------------------------------------------"):
                section = 4
        elif (section == 4):
            if ((len(index_elements) > 0) and (index_elements[0] == "Source")):
                sourcearea = 1
            elif ((len(index_elements) > 2) and (sourcearea == 1)):
                src = { 'name': index_elements[0], 'calcode': index_elements[1],
                        'right_ascension': index_elements[2], 'declination': index_elements[3],
                        'dra': index_elements[4], 'ddec': index_elements[5] }
                rd['sources'].append(src)

    # Convert things into numpy arrays for easy where-ing later.
    rd['index']['time'] = np.array(rd['index']['time'])
    rd['index']['source'] = np.array(rd['index']['source'])
    rd['index']['calcode'] = np.array(rd['index']['calcode'])
    rd['index']['antennas'] = np.array(rd['index']['antennas'])
    rd['index']['spectral_channels'] = np.array(rd['index']['spectral_channels'])
    rd['index']['wideband_channels'] = np.array(rd['index']['wideband_channels'])
    rd['index']['freq_config'] = np.array(rd['index']['freq_config'])
    rd['index']['record_number'] = np.array(rd['index']['record_number'])
    
    return rd

def split_into_segments(idx):
    # We go through a uvindex dictionary and return segments.
    # Each segment is a single source, at a single frequency,
    # with a start and end time.
    segs = []
    oldsrc = ""
    oldconfig = -1
    sseg = None
    for i in xrange(0, len(idx['index']['source'])):
        if ((idx['index']['source'][i] != oldsrc) or
            (idx['index']['freq_config'][i] != oldconfig)):
            if ((oldsrc != "") and (oldconfig != -1)):
                # Put the segment on the list.
                segs.append(sseg)
            oldsrc = idx['index']['source'][i]
            oldconfig = idx['index']['freq_config'][i]
            sseg = { 'source': idx['index']['source'][i],
                     'freq_config': idx['index']['freq_config'][i],
                     'start_time': ephem.Date(idx['index']['time'][i]),
                     'end_time': ephem.Date(idx['index']['time'][i]) }
        else:
            sseg['end_time'] = ephem.Date(idx['index']['time'][i])
    # Have to push the last segment on.
    segs.append(sseg)
    return segs

def dataset_find(srcname, freq=None):
    srcpat = "%s.*" % srcname
    if (freq is not None):
        srcpat = "%s.%d" % (srcname, freq)
    matches = []
    for root, dirnames, filenames in os.walk('.'):
        for filename in fnmatch.filter(dirnames, srcpat):
            # Check if this has data.
            fname = os.path.join(root, filename)
            cname = "%s/visdata" % fname
            if ((os.path.isdir(fname)) and (os.path.isfile(cname))):
                matches.append(fname)
    return matches

def filter_uvplt(output):
    outlines = output.split('\n')

    rd = { 'nvisibilities': 0 }
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if (len(index_elements) < 3):
            continue
        if ((index_elements[0] == "Read") and
            (index_elements[2] == "visibilities")):
            rd['nvisibilities'] = int(index_elements[1])
    return rd

def filter_closure(output):
    outlines = output.split('\n')

    rd = { 'theoretical_rms': 0, 'measured_rms': 0 }
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if (len(index_elements) < 1):
            continue
        if (index_elements[0] == "Actual"):
            rd['measured_rms'] = float(index_elements[-1])
        elif (index_elements[0] == "Theoretical"):
            rd['theoretical_rms'] = float(index_elements[-1])
    return rd

def filter_uvfstats(output):
    outlines = output.split('\n')

    rd = { 'flagged_fraction': 0 }
    strt = 0
    nchans = 0
    nflagged = 0
    for i in xrange(0, len(outlines)):
        index_elements = outlines[i].split()
        if (len(index_elements) < 1):
            continue
        if (strt == 1):
            nchans = nchans + 1
            if (index_elements[1] < 15):
                nflagged = nflagged + 1
        else:
            if (index_elements[0] == "-------"):
                strt = 1
    rd['flagged_fraction'] = "%.2f" % (nflagged / nchans)
    return rd

def calibrate(srcname, calname, fconfig, stime, etime):
    smtime = datetime_to_mirtime(stime)
    emtime = datetime_to_mirtime(etime)
    selstring = "time(%s,%s)" % (smtime, emtime)

    print "  calibrating source %s with selection %s" % (srcname, selstring)

    # Find all the relevant datasets.
    dsets = dataset_find(srcname)
    
    cfreqs = []
    rv = { 'code': 0, 'frequencies': [] }

    for i in xrange(0, len(dsets)):
        fname = dsets[i]
        # Check if this is one of the frequencies in this configuration.
        setf = int(fname.split(".")[-1])
        for j in xrange(0, len(fconfig['nchannels'])):
            c = (fconfig['frequency1'][j] + (fconfig['nchannels'][j] - 1) *
                 fconfig['frequency_increment'][j] / 2.) * 1000.
            fdiff = abs(c - setf)
            if (fdiff <= abs(fconfig['frequency_increment'][j] * 1000.)):
                cfreqs.append([ setf, j, fname ])
                rv['frequencies'].append(setf)

    # Calibrate per frequency.
    for i in xrange(0, len(cfreqs)):
        dset = cfreqs[i][2]
        csets = dataset_find(calname, cfreqs[i][0])
        cset = csets[0]
        # Check that we actually have data in this time range.
        miriad.set_filter('uvplt', filter_uvplt)
        
        uvout = miriad.uvplt(vis=dset, axis="time,amp", device="/null",
                             options="nopol,nocal,nopass", stokes="xx,yy",
                             select=selstring)
        if (uvout['nvisibilities'] > 0):
            # Do the calibration.
            print "  calibrating frequency %d MHz" % cfreqs[i][0]
            miriad.gpcal(vis=dset, interval="0.1", options="xyvary,nopol,qusolve",
                         nfbin="2", refant="3", select=selstring)
            miriad.gpboot(vis=dset, cal=cset, select=selstring)
            rv['code'] = 1
        else:
            print "  Unable to find any data for this time range!" % selstring
            rv['code'] = 0
    return rv


def measure_closure_phase(srcname, freq, stime, etime):
    print "    closure phase"
    closurelog = "closure_log.txt"
    selstring = "time(%s,%s)" % (datetime_to_mirtime(stime),
                                 datetime_to_mirtime(etime))
    if (os.path.isfile(closurelog)):
        os.remove(closurelog)
    # Find the data set.
    dsets = dataset_find(srcname, freq)
    miriad.set_filter('closure', filter_closure)
    cout = miriad.closure(vis=dsets[0], stokes="i", device="/null",
                          options="log", select=selstring)

    rv = { 'closure_phase': { 'theoretical_rms': cout['theoretical_rms'],
                              'measured_rms': cout['measured_rms'],
                              'average_value': -999 } }
    if (os.path.isfile(closurelog)):
        with open(closurelog, "r") as fp:
            loglines = fp.readlines()
        pvals = []
        for i in xrange(0, len(loglines)):
            lels = loglines[i].split()
            if (len(lels) < 1):
                continue
            if (lels[0] == "Antennas"):
                continue
            pvals.append(float(lels[-1]))
        rv['closure_phase']['average_value'] = np.average(pvals)
    return rv

def measure_flagging_statistic(srcname, freq, stime, etime):
    print "    flagging"
    miriad.set_filter('uvfstats', filter_uvfstats)
    selstring = "time(%s,%s)" % (datetime_to_mirtime(stime),
                                 datetime_to_mirtime(etime))

    # Find the data set.
    dsets = dataset_find(srcname, freq)
    uo = miriad.uvfstats(vis=dsets[0], mode="channel",
                         options="absolute,unflagged",
                         select=selstring)
    return uo['flagged_fraction']

def measure(srcname, datadir, coords, fconfig, stime, etime, calresult):
    print "  measuring source %s parameters" % srcname

    ro = { 'source': srcname, 'closurePhase': [], 'flaggedFraction': [] }

    stokes = [ 'i' ]
    orders = [ 3 ]
    options = [ ",log" ]

    plotdir = "plots"
    if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)

    for i in xrange(0, len(stokes)):
        p = "%s/%s" % (plotdir, srcname)
        o = "%s/%s" % (datadir, srcname)
        if (stokes[i] != "i"):
            o = "%s.%s" % (o, stokes[i])
        o = "%s.plotgen" % o

        if (i == 0):
            # Get the closure phase and flagging statistic.
            for k in xrange(0, len(calresult['frequencies'])):
                clop = measure_closure_phase(srcname, calresult['frequencies'][k],
                                             stime, etime)
                ro['closurePhase'].append({
                    'IF': calresult['frequencies'][k],
                    'average_value': clop['closure_phase']['average_value'],
                    'measured_rms': clop['closure_phase']['measured_rms'],
                    'theoretical_rms': clop['closure_phase']['theoretical_rms']
                    })
                ro['flaggedFraction'].append(
                    measure_flagging_statistic(srcname, calresult['frequencies'][k],
                                               stime, etime))
            # Calculate the hour angle ranges.
    return ro
                

def calibrate_and_measure(args):
    # Run uvindex on the dataset.
    dataset_name = args['<dataset>']
    miriad.set_filter('uvindex', filter_uvindex)
    index_data = miriad.uvindex(vis=dataset_name, interval='0.1')

    # Get the cycle time.
    uvlist_log_name = "uvlist.log"
    if (os.path.isfile(uvlist_log_name)):
        os.remove(uvlist_log_name)
    miriad.uvlist(vis=dataset_name, options="variables,full",
                  log=uvlist_log_name)
    telescope_variables = filter_uvlist_variables(uvlist_log_name)
    cycle_time = round(telescope_variables['cycle_time'])

    # Split the observations up into segments.
    segments = split_into_segments(index_data)

    # We will write out files per source.
    source_data = {}
    
    # Each segment is a source for some time, so we go through
    # them in order.
    for i in xrange(0, len(segments)):
        # Check it isn't the calibrator source.
        if (segments[i]['source'] != args['--calibrator']):
            source_name = segments[i]['source']
            print " Source %s is observed from %s to %s" % (source_name,
                                                            segments[i]['start_time'],
                                                            segments[i]['end_time'])
            segment_length = (segments[i]['end_time'].datetime() -
                              segments[i]['start_time'].datetime()).total_seconds() / 60.
            #print " segment length is %.2f min" % (segment_length)
            iter_starttime = segments[i]['start_time'].datetime()
            finished = False
            if source_name not in source_data:
                # Find the source in the index data.
                for j in xrange(0, len(index_data['sources'])):
                    if (index_data['sources'][j]['name'] == source_name):
                        source_data[source_name] = {
                            'sourceName': source_name,
                            'rightAscension': index_data['sources'][j]['right_ascension'],
                            'declination': index_data['sources'][j]['declination'],
                            'arrayConfiguration': "",
                            'timeSeries': []
                        }
                        break
            while not finished:
                iter_endtime = segments[i]['end_time'].datetime()
                segment_length = (iter_endtime - iter_starttime).total_seconds() / 60.
                if ((segment_length < args['--segment']) or
                    (segment_length < (args['--segment'] + args['--slop']))):
                    # This will be all we need to do.
                    finished = True
                else:
                    iter_endtime = iter_starttime + timedelta(minutes=args['--segment'])
                # Calibrate this chunk.
                fcn = segments[i]['freq_config'] - 1
                #print index_data['freq_configs'][fcn]
                calout = calibrate(segments[i]['source'], args['--calibrator'],
                                   index_data['freq_configs'][fcn],
                                   iter_starttime, iter_endtime)
                if (calout['code'] == 1):
                    # The calibration worked, we can now make a measurement.
                    mset = measure(segments[i]['source'], "data", None,
                                   index_data['freq_configs'][fcn],
                                   iter_starttime, iter_endtime, calout)
                    tso = { 'closurePhase': mset['closurePhase'],
                            'centreFreqs': calout['frequencies'] }
                    source_data[source_name]['timeSeries'].append(tso)
                iter_starttime = iter_endtime + timedelta(seconds=cycle_time)

    # Output the JSON files.
    for src in source_data:
        jfile = "data/%s.json" % src
        if (not os.path.isdir("data")):
            os.makedirs("data")
        with open(jfile, "w+") as fp:
            json.dump(source_data[src], fp)
            
if __name__ == '__main__':
    arguments = docopt(__doc__, version="Time-variable calibrator 1.0")
    valid = True
    # Do some existence checking.
    if (('<dataset>' not in arguments) or
        (not os.path.isdir(arguments['<dataset>']))):
        print "No valid dataset found."
        valid = False
    if (valid == True):
        calibrate_and_measure(arguments)
