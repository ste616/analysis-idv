# This is analysis-idv/cabb_json/reader.py
# Jamie Stevens 2017
# Licensed under GPL 3.0 or later. A copy of the license can
# be found at https://www.gnu.org/licenses/gpl.html

# This module contains classes and methods to store and read CABB flux density
# data written in a JSON format.

# Our necessary imports.
import json
import numpy as np

# Spectrum class: this stores a single spectrum.
class Spectrum:
    def __init__(self, dataObject=None):
        # We store the type of spectrum this is.
        self.averagingMode = ""
        # And then an array of N=2 lists, with frequency in N=0, amplitude in N=1.
        self.spectrum = []
        self.storeFromObject(dataObject)
        return None

    def storeFromObject(self, dataObject=None):
        if dataObject is not None:
            if 'stokes' in dataObject and 'mode' in dataObject and 'data' in dataObject:
                # This looks like the right type of object.
                self.averagingMode = dataObject['mode']
                self.spectrum = dataObject['data']
        return None

    def getAveragingMode(self):
        return self.averagingMode

    def getSpectrum(self):
        return self.spectrum

    def getSpectrumArrays(self):
        try:
            self.freq = self.npSpectrum[:,0]
        except AttributeError:
            self.npSpectrum = np.array(self.spectrum)
            self.freq = self.npSpectrum[:,0]
        self.amp = self.npSpectrum[:,1]
        return { 'freq': self.freq, 'amp': self.amp }

# Measurement class: this stores a spectrum along with some metadata.
class Measurement:
    def __init__(self, measureObject=None, stokes=None):
        # We keep one spectrum.
        self.spectrum = None
        # The metadata includes:
        # 1. The central frequencies measured at this time.
        self.centreFrequencies = []
        # 2. The fraction flagged for each band.
        self.flaggedFraction = []
        # 3. The closure phase measurements for each band.
        self.closurePhase = []
        # 4. The defect overall.
        self.defect = None
        # 5. The time quantities, including MJD and hour angles.
        self.mjd = None
        self.hourAngle = None
        # We only care about a single Stokes, which we need to be told about.
        if stokes is None:
            stokes = "I"
        self.stokes = stokes
        self.storeFromObject(measureObject)
        return None

    def storeFromObject(self, measureObject=None):
        if measureObject is not None:
            if 'defect' in measureObject and 'mjdRange' in measureObject and 'centreFreqs' in measureObject:
                # This looks right.
                # Store our own metadata first.
                self.centreFrequencies = measureObject['centreFreqs']
                self.flaggedFraction = measureObject['flaggedFraction']
                for i in xrange(0, len(measureObject['closurePhase'])):
                    self.closurePhase.append({
                        'averageValue': measureObject['closurePhase'][i]['average_value'],
                        'measuredRms': measureObject['closurePhase'][i]['measured_rms']
                    })
                for i in xrange(0, len(measureObject['defect'])):
                    if measureObject['defect'][i]['stokes'] == self.stokes:
                        self.defect = measureObject['defect'][i]
                self.mjd = measureObject['mjdRange']
                self.hourAngle = measureObject['hourAngleRange']
                # Now find the right spectrum and add it.
                for i in xrange(0, len(measureObject['fluxDensityData'])):
                    if measureObject['fluxDensityData'][i]['stokes'] == self.stokes:
                        self.spectrum = Spectrum(measureObject['fluxDensityData'][i])
        return None
    
# TimeSeries class: this stores a series of spectra.
class TimeSeries:
    def __init__(self, seriesObject=None, stokes=None):
        # We keep an array of Measurement objects.
        self.measurements = []
        # We have metadata about each measurement.
        # 1. The array configuration.
        self.arrayConfigurations = []
        # We only care about a single Stokes, which we need to be told about.
        if stokes is None:
            stokes = "I"
        self.storeFromObject(seriesObject, stokes)
        return None

    def storeFromObject(self, seriesObject=None, stokes=None):
        if seriesObject is not None:
            if 'arrayConfiguration' in seriesObject and 'timeSeries' in seriesObject:
                # This looks right.
                # Go through the time series.
                for i in xrange(0, len(seriesObject['timeSeries'])):
                    self.arrayConfigurations.append(seriesObject['arrayConfiguration'])
                    self.measurements.append(Measurement(seriesObject['timeSeries'][i], stokes))
        return None

# Source class: this stores information about a source, with time series of all
# the types.
class Source:
    def __init__(self, sourceObject=None):
        # We keep several TimeSeries, one for each Stokes.
        self.timeSeries = { 'I': None, 'Q': None, 'U': None, 'V': None }
        # Our metadata includes our position and name.
        self.rightAscension = ""
        self.declination = ""
        self.name = ""
        self.storeFromObject(sourceObject)
        return None

    def storeFromObject(self, sourceObject=None):
        if sourceObject is not None:
            if 'sourceName' in sourceObject and 'rightAscension' in sourceObject:
                # This looks right.
                self.rightAscension = sourceObject['rightAscension']
                self.declination = sourceObject['declination']
                self.name = sourceObject['sourceName']
                # Try to load each time series.
                if 'timeSeries' in sourceObject:
                    for stokes in self.timeSeries:
                        if self.timeSeries[stokes] is None:
                            self.timeSeries[stokes] = TimeSeries(sourceObject, stokes)
                        else:
                            self.timeSeries[stokes].storeFromObject(sourceObject, stokes)
        return None

# The routine which reads in a JSON file.
def readJson(fileName=None):
    if fileName is None:
        return None

    # Try to open the file.
    with open(fileName, 'r') as fp:
        fileContents = json.load(fp)

    # Act based on what type of file this is.
    if 'sourceName' in fileContents and 'timeSeries' in fileContents:
        # This is likely an IHV file.
        return Source(fileContents)
        
