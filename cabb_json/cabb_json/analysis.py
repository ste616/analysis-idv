# This is analysis-idv/cabb_json/analysis.py
# Jamie Stevens 2017
# Licensed under GPL 3.0 or later. A copy of the license can
# be found at https://www.gnu.org/licenses/gpl.html

# This module contains the routines that can do analysis on the time-series
# of the IHV sources.

# Our necessary imports.
import numpy as np
import math
from scipy.optimize import curve_fit

# Define model function to be used to fit the time-scale.
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def userGauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# Define model function where amplitude and mean is fixed, because ACF needs
# to have a peak strength of 1.0 at lag = 0.0
def acfGauss(x, *p):
    sigma, = p
    return 1.0 * np.exp(-1. * x**2 / (2. * sigma**2))

# The width of a Gaussian at half height.
def gaussHalfHeight(sigma):
    return 2. * math.sqrt(2. * math.log(2.)) * sigma

# The width of a Gaussian at 1/e height.
def gaussEHeight(sigma):
    return 2. * math.sqrt(2.) * sigma

# Return the width of a Gaussian with sigma, at one of two heights.
def gaussHeight(sigma, sigmaError=0., height="half"):
    heightFn = gaussHalfHeight
    if height == "e":
        heightFn = gaussEHeight
    width = heightFn(sigma)
    widthPlus = abs(width - heightFn(sigma + sigmaError))
    widthMinus = abs(width - heightFn(sigma - sigmaError))
    widthError = (widthPlus + widthMinus) / 2.
    return (width, widthError)

def calculateTimescale(acfLag, acfCor, acfCorError, mode='fwhm', tfitmin=0.1, tfitmax=0.8):
    # We take the ACF as a function of lag, and determine what the timescale is.

    # Our initial guess at the Gaussian parameters.
    # [ std dev ]
    p0 = [ 1. ]

    # We restrict the points that can be used to do the fitting. We do this
    # with numpy arrays because it's easier.
    aLag = np.array(acfLag)
    aCor = np.array(acfCor)
    aError = np.array(acfCorError)
    # Work out the region where the values have fallen below the cutoff.
    outsideNegative = np.where((aCor < tfitmin) & (aLag < 0))
    outsidePositive = np.where((aCor < tfitmin) & (aLag > 0))
    # We want to know the smallest negative number and smallest positive number as the
    # index bounds.
    if len(outsideNegative[0]) > 0:
        smallestNegative = np.max(outsideNegative)
    else:
        smallestNegative = 0
    if len(outsidePositive[0]) > 0:
        smallestPositive = np.min(outsidePositive)
    else:
        smallestPositive = None
    lagBoundNegative = aLag[smallestNegative + 1]
    if smallestPositive is not None:
        lagBoundPositive = aLag[smallestPositive]
    else:
        lagBoundPositive = aLag[-1]

    # Do the fit.
    popt, pcov = curve_fit(acfGauss, aLag[smallestNegative + 1:smallestPositive],
                           aCor[smallestNegative + 1:smallestPositive],
                           sigma=aError[smallestNegative + 1:smallestPositive], p0=p0)

    # We return the lag at some value on the Gaussian fit.
    rval = { 'value': None, 'mode': None, 'sigma': popt[0], 'timeUnits': 'minutes',
             'sigmaError': math.sqrt(pcov[0,0]), 'valueError': None,
             'fitRegion': [ lagBoundNegative, lagBoundPositive ] }
    if mode == 'fwhm' or mode == 'hwhm':
        rval['mode'] = mode
        # We want to know the width at the half amplitude.
        (fwhm, ewhm) = gaussHeight(popt[0], math.sqrt(pcov[0,0]))
        if mode == 'fwhm':
            rval['valueError'] = ewhm
            rval['value'] = fwhm
        elif mode == 'hwhm':
            rval['valueError'] = ewhm / 2.
            rval['value'] = fwhm / 2.
    elif mode == 'fwhme' or mode == 'hwhme':
        rval['mode'] = mode
        # We want to know the width at amplitude 1/e.
        (fwhme, ewhme) = gaussHeight(popt[0], math.sqrt(pcov[0,0]), height='e')
        if mode == 'fwhme':
            rval['valueError'] = ewhme
            rval['value'] = fwhme
        elif mode == 'hwhme':
            rval['valueError'] = ewhme / 2.
            rval['value'] = fwhme / 2.
    return rval

def calculateACF(timeSeries):
    if timeSeries['timeUnits'].lower() != "mjd":
        print "ACF analysis requires MJD times."
        return None

    # We do the ACF for each frequency separately.
    retVal = { 'cor': [], 'lag': [], 'frequencies': [], 'frequencyUnits': timeSeries['frequencyUnits'],
               'corError': [], 'lagError': [], 'timeUnits': "minutes" }
    for i in xrange(0, len(timeSeries['frequencies'])):
        times = np.array(timeSeries['times'][i]) - min(timeSeries['times'][i])
        #print "times"
        # Convert into minutes.
        times *= 1440.
        #print times
        fluxes = np.array(timeSeries['fluxDensities'][i])
        #print "fluxes"
        #print fluxes
        start = np.min(times)
        stop = np.max(times)
        #print "start / stop"
        #print start
        #print stop
        trange = stop - start
        #print "trange"
        #print trange
        minlag = -trange
        maxlag = trange
        # Check for some unusable conditions.
        if 0. in fluxes:
            continue
        if len(fluxes) < 3:
            continue
        acfData = zdcf(timeSeries['times'][i], fluxes, timeSeries['times'][i], fluxes, minlag=minlag, maxlag=maxlag)
        retVal['cor'].append(acfData['cor'])
        retVal['corError'].append([ acfData['corErrMinus'], acfData['corErrPlus'] ])
        retVal['lagError'].append(acfData['lagErr'])
        retVal['lag'].append(acfData['lag'])
        retVal['frequencies'].append(timeSeries['frequencies'][i])
    return retVal

class TimeLagPair:
    def __init__(self, aTime, aAmp, aIndex, bTime, bAmp, bIndex, timeUnits='minutes'):
        # The time should be supplied in MJD.
        #print "times suppled for pair = %f %f" % (aTime, bTime)
        self.times = [ aTime, bTime ]
        self.amps = [ aAmp, bAmp ]
        self.indices = [ aIndex, bIndex ]
        # We will calculate the time in the time unit.
        self.timeLag = aTime - bTime
        if timeUnits == 'minutes':
            # Converting from days to minutes.
            self.timeLag *= 1440.
        elif timeUnits == 'seconds':
            # Converting from days to seconds.
            self.timeLag *= 86400.
        elif timeUnits == 'hours':
            # Converting from days to hours.
            self.timeLag *= 24.
        self.used = False
        return None

    def __lt__(self, other):
        return self.timeLag < other.timeLag

def zdcf(timesA, ampsA, timesB, ampsB, minlag=0., maxlag=None, minpt=11, epsilon=0.5, timeUnits='minutes'):
    # Make the DCF using the z-transformed direct correlation function.
    # Described in the Alexander (2013) paper.

    # First step is to calculate the time delay for each pair of points.
    pairs = []
    for i in xrange(0, len(timesA)):
        for j in xrange(0, len(timesB)):
            pairs.append(TimeLagPair(timesA[i], ampsA[i], i, timesB[j], ampsB[j], j, timeUnits=timeUnits))
    # Now sort according to time lag.
    pairs.sort()

    # We exit now if we have less than the required number of points for a single bin.
    if len(pairs) < minpt:
        return None

    # Now go through from the minimum lag, adding the minimum number of points to begin with.
    lag = []
    lagErr = []
    cor = []
    corMinErr = []
    corMaxErr = []
    finished = False
    loopNumber = 0
    while finished == False and loopNumber < 100:
        # Go through each of the pairs in turn.
        binPoints = []
        loopNumber += 1
        #print "LL this is loop number %d" % loopNumber
        #print "minlag maxlag = %f %f" % (minlag, maxlag)
        for i in xrange(0, len(pairs)):
            # Has this pair already been used?
            #print "investigating pair %d" % i
            #print pairs[i].timeLag
            #print pairs[i].used
            if pairs[i].used == True:
                continue
            # Is this lag within the lag limits?
            if minlag is not None and pairs[i].timeLag < (minlag - epsilon):
                #print "lag too small"
                continue
            if maxlag is not None and pairs[i].timeLag > (maxlag + epsilon):
                #print "lag too large"
                continue
            # We now add think about adding this to our list.
            #print "could we add this point?"
            canAdd = False
            if len(binPoints) >= minpt:
                # We can still add this point, if it is close to the last added lag.
                lagdiff = pairs[i].timeLag - binPoints[-1].timeLag
                if lagdiff <= epsilon:
                    #print "may add small offset pair"
                    canAdd = True
                else:
                    # Nothing else can be added.
                    #print "cannot add this pair, stopping"
                    break
            else:
                # We don't yet have enough, so we add this point.
                #print "pair can be added"
                canAdd = True
            if canAdd == True:
                #print "point may be added, if possible"
                # We check if one of these values is already in the bin.
                #print "checking for interdependency"
                for j in xrange(0, len(binPoints)):
                    if (pairs[i].indices[0] == binPoints[j].indices[0] or
                        pairs[i].indices[1] == binPoints[j].indices[1]):
                        # We discard this point by marking it as used.
                        pairs[i].used = True
                        break
                if pairs[i].used == True:
                    #print "pair interdependent"
                    continue
                # We get here, we can add.
                #print "adding pair"
                pairs[i].used = True
                binPoints.append(pairs[i])
        # We're out of the loop now.
        #print "there are %d points in this bin" % len(binPoints)
        if len(binPoints) < minpt:
            # We can no longer make any valid bins.
            finished = True
            #print "we're done"
            break
        # Otherwise, we can calculate the correlation for this bin.
        # Calculate the average amplitude and time lag.
        #print "calculating r"
        meanAmpA = 0.
        meanAmpB = 0.
        meanLag = 0.
        nAmp = 0
        for i in xrange(0, len(binPoints)):
            meanAmpA += binPoints[i].amps[0]
            meanAmpB += binPoints[i].amps[1]
            meanLag += binPoints[i].timeLag
            nAmp += 1
        if nAmp > 0:
            meanAmpA /= float(nAmp)
            meanAmpB /= float(nAmp)
            meanLag /= float(nAmp)
        #print "mean amplitudes A B = %f %f" % (meanAmpA, meanAmpB)
        #print "mean lag = %f" % meanLag
        # Now calculate the standard deviation of the amps and lags.
        stdAmpA = 0.
        stdAmpB = 0.
        stdLag = 0.
        for i in xrange(0, len(binPoints)):
            stdAmpA += (binPoints[i].amps[0] - meanAmpA)**2
            stdAmpB += (binPoints[i].amps[1] - meanAmpB)**2
            stdLag += (binPoints[i].timeLag - meanLag)**2
        stdDevA = stdAmpA / float(nAmp - 1)
        stdDevB = stdAmpB / float(nAmp - 1)
        stdDevLag = stdLag / float(nAmp - 1)
        #print "standard deviation A B = %f %f" % (stdDevA, stdDevB)
        #print "standard deviation lag = %f" % math.sqrt(stdDevLag)
        # Now calculate r.
        r = 0.
        for i in xrange(0, len(binPoints)):
            r += (binPoints[i].amps[0] - meanAmpA) * (binPoints[i].amps[1] - meanAmpB) / float(nAmp - 1)
        #print "r numerator %f" % r
        r /= (math.sqrt(stdDevA) * math.sqrt(stdDevB))
        #print "r full %f" % r
        # From that we get z.
        if r < 1:
            z = 0.5 * math.log10( (1 + r) / (1 - r) )
            xi = z
            zbar = xi + ((r / (2. * (nAmp - 1))) * (1. + ((5 + r**2) / (4. * (nAmp - 1))) +
                                                    ((11. + 2. * r**2 + 3. * r**4) / (8. * (nAmp - 1)**2))))
            szsquared = (1 / (nAmp - 1)) * (1. + ((4 - r**2) / (2. * (nAmp - 1))) +
                                            ((22. - 6. * r**2 - 3. * r**4) / (6. * (nAmp - 1)**2)))
            deltarplus = abs(math.tanh(zbar + math.sqrt(szsquared)) - r)
            deltarminus = abs(math.tanh(zbar - math.sqrt(szsquared)) - r)
        else:
            # We can't use the z equation because 1 / (1 - r) is bad.
            # But essentially, there is no error on r when r = 1, because it is perfect
            # auto-correlation.
            deltarplus = 0.01
            deltarminus = 0.01
        #print "z now %f" % z
        # Form our arrays.
        lag.append(meanLag)
        lagErr.append(math.sqrt(stdDevLag))
        cor.append(r)
        corMinErr.append(deltarplus)
        corMaxErr.append(deltarminus)

    # Return the arrays.
    return { 'lag': lag, 'cor': cor, 'lagErr': lagErr,
             'corErrMinus': corMinErr, 'corErrPlus': corMaxErr }

# This following routine is from Andrew O'Brien.
'''
Discrete Correlation Function.
INPUTS:
	ta: time stamps of first light curve.
	ra: rate (flux densities) of first light curve.
	tb: time stamps of second light curve.
	rb: rate (flux densities) of second light curve.

OPTIONAL INPUTS:
	erra: rate errors of first light curve.
	errb: rate errors of second light curve.
	minlag: minimum lag value to consider.
	maxlag: maximum lag value to consider.
	numf: number of lag bins.
	minpt: minimum number of data points required per bin for DCF to be computed.
'''
def dcf(ta, ra, tb, rb, erra = np.empty(0), errb = np.empty(0), minlag = None, maxlag = None, numf = None, minpt = None):
	# characteristic error of time series A
	if(erra.size == 0):	erra = 0
	elif(erra.size > 1): erra = np.mean(erra) # should I use the average error?

	#characteristic error of time series B
	if(errb.size == 0): errb = 0
	elif(errb.size > 1): errb = np.mean(errb) # should I use the average error?

	# minimum number of points to consider a CCF point valid
	if(minpt is None): minpt = 10
	if(minpt < 0):
		print "ERROR: minpt must be positive."
		return

	if(ra.size != ta.size):
		print "ERROR: ta and ra must have same length."
		return
	if(rb.size != tb.size):
		print "ERROR: tb and rb must have same length."
		return

	# number of lags
	if(numf is None): numf = int(min(ra.size,rb.size)/10)

	udcf = np.zeros([ra.size,rb.size])
	dt = np.zeros([ra.size,rb.size])
	
	# This could be done with a Python nested loop, but using NumPy's broadcasting
	# will speed up computation
	for i in range(0, ra.size):
		udcf[i,:] = ((ra[i] - ra.mean()) * (rb - rb.mean()))
		dt[i,:] = ta[i]-tb

	udcf = udcf / np.sqrt((ra.std()**2 - erra**2) * (rb.std()**2 - errb**2))
	#udcf = udcf / np.sqrt((ra.std() - erra)**2 * (rb.std() - errb)**2)
        
	if(minlag is None): minlag = np.amin(dt)
	if(maxlag is None): maxlag = np.amax(dt)

	lag = minlag + (maxlag - minlag) * np.arange(numf+1)/numf
	
	# now sort all coefficients as a function of the lag
	# (makes finding the individual elements easier and
	# significantly speeds up the code)
	dt = dt.flatten()
	udcf = udcf.flatten()
	ndx = dt.argsort() # numpy magic to sort an array by another array
	dt = dt[ndx]
	udcf = udcf[ndx]
	#dt = dt[ndx]
	#udcf = udcf[ndx]
	
	
	npt1 = dt.size - 1
	cor = np.zeros(numf)
	cor.fill(np.nan)
	numpt = np.zeros(numf)
	sta = 0
	
	'''
	while(sta < npt1):
		for k in dt[sta]:
			#print k
			if k < lag[0]:
				sta += 1
	'''
	while(sta < npt1 and dt[sta] < lag[0]): sta += 1
	sto = sta

	for i in range(0, lag.size - 1):
		# search next element in dt that is greater than the current lag
		'''
		while(sto < npt1):
			for k in dt[sto]:
				if k < lag[i+1]: sto += 1
		'''
		while(sto < npt1 and dt[sto] < lag[i+1]): sto += 1
		# compute the mean correlation if we have any udcf elements
		if(sta < sto):
			if(sto == npt1 and dt[sto] < lag[i+1]):
				cor[i] = np.mean(udcf[sta:sto])
				numpt[i] = sto - sta
			else:
				cor[i] = np.mean(udcf[sta:sto-1])
				numpt[i] = sto - sta - 1
			sta = sto
	
	lag = (lag[:-1] + lag[1:])/2
	
	ndx = np.flatnonzero(numpt < minpt) # flatnonzero returns indexes of non-zero elements in a flattened array (np.nonzero() returns results in a tuple per dimension)
	if(ndx.size != 0): cor[ndx] = np.nan
	
	return (lag, cor)
