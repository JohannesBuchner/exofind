#

import numpy
import json
import sys

n = len(numpy.loadtxt(sys.argv[1], ndmin = 2))
data = [numpy.loadtxt(d, ndmin = 2) for d in sys.argv[2:]]
fits = [(d[:,1].min(), d.shape[1] - 2) for d in data]

def aicc(bestfit, nparams):
	""" bestfit is -2 * loglikelihood of best fit; nparams is number of model parameters """
	k = nparams
	return (2 * k + bestfit + 2 * k * (k + 1) / (n - k - 1.)) / 2.

for (bestfit, nparams), filename in zip(fits, sys.argv[2:]):
	print filename, nparams, bestfit, aicc(bestfit, nparams)


