import matplotlib.pyplot as plt
import numpy as np
import cabb_json as ihv
from scipy.optimize import curve_fit

# Make a Gaussian with mean 0, peak of 1, and a standard deviation of some number.
gMean = 0.
gPeak = 1.
gStdDev = np.arange(1, 40000, 5)
gMeasured = []

gXvalues = np.arange(-180, 20, 1)

for i in xrange(0, len(gStdDev)):
    gYvalues = np.array([ ihv.userGauss(a, gPeak, gMean, gStdDev[i]) for a in gXvalues ])
    gEvalues = 0.001 * np.random.randn(len(gXvalues))
    gZvalues = gYvalues + gEvalues

    # Now do a fit to this data.
    popt, pcov = curve_fit(ihv.acfGauss, gXvalues, gZvalues, sigma=0.001, p0=[1.])
    print "standard deviation measured = %.3f" % popt[0]
    gMeasured.append(popt[0])


#plt.plot(gXvalues, gYvalues, 'o')
#plt.plot(gXvalues, gZvalues, 'o')
plt.plot(gStdDev, gMeasured, 'o')

plt.show()
