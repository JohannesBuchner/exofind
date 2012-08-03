import pymultinest
import math
import os, sys
import numpy

from ctypes import cdll
#cdll.LoadLibrary('libnest3.so')
lib = cdll.LoadLibrary('exo.so')
from ctypes import *

# put data into double[n][3] array
data = numpy.loadtxt(sys.argv[1])
cdata = ((c_double * 3) * data.shape[0])()
for i in range(data.shape[0]):
	for j in range(3):
		cdata[i][j] = data[i,j]
lib.set_data(cdata, data.shape[0])

def plot(a, parameters):
	import matplotlib.pyplot as plt
	p = pymultinest.PlotMarginal(a)

	n_params = len(parameters)
	values = a.get_equal_weighted_posterior()
	s = a.get_stats()
	assert n_params == len(s['marginals'])
	prefix = a.outputfiles_basename

	plt.figure(figsize=(5*n_params, 5*n_params))
	for i in range(n_params):
		plt.subplot(n_params, n_params, n_params * i + i + 1)
		plt.xlabel(parameters[i])
	
		m = s['marginals'][i]
		plt.xlim(m['5sigma'])
	
		oldax = plt.gca()
		x,w,patches = oldax.hist(values[:,i], bins=20, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
		oldax.set_ylim(0, x.max())
	
		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
		newax.set_ylim(0, 1)
	
		ylim = newax.get_ylim()
		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
		center = m['median']
		low1, high1 = m['1sigma']
		newax.errorbar(x=center, y=y,
			xerr=numpy.transpose([[center - low1, high1 - center]]), 
			color='blue', linewidth=2, marker='s')
		oldax.set_yticks([])
		#newax.set_yticks([])
		newax.set_ylabel("Probability")
		ylim = oldax.get_ylim()
		newax.set_xlim(m['5sigma'])
		oldax.set_xlim(m['5sigma'])
		#plt.close()
	
		for j in range(i):
			plt.subplot(n_params, n_params, n_params * j + i + 1)
			p.plot_conditional(i, j, bins=20, cmap = plt.cm.gray_r)
			plt.xlabel(parameters[i])
			plt.ylabel(parameters[j])
			#plt.savefig('cond_%s_%s.pdf' % (params[i], params[j]), bbox_tight=True)
			#plt.close()

	plt.savefig(prefix + 'marg.pdf')
	plt.close()

def run(n_planets):
	# number of dimensions our problem has
	parameters = ['V'] 
	# plus possibly a noise term 
	# parameters += ['s']
	for i in range(n_planets):
		parameters += ['%s%d' % (var, i) for var in ['P', 'K', 'chi', 'e', 'omega']]
	n_params = len(parameters)
	params_low = (c_double * n_params)()
	params_high = (c_double * n_params)()
	for i, p in enumerate(parameters):
		if p.startswith('V'):
			params_low[i], params_high[i] = -2000., +2000.
		if p.startswith('P'):
			params_low[i], params_high[i] = 0.2, 10000.
		if p.startswith('K'):
			params_low[i], params_high[i] = 1., 20000
		if p.startswith('chi'):
			params_low[i], params_high[i] = 0., 1.
		if p.startswith('e'):
			params_low[i], params_high[i] = 0., 1.
		if p.startswith('omega'):
			params_low[i], params_high[i] = 0., 1.
		if p.startswith('s'):
			params_low[i], params_high[i] = 1., 2000
	lib.set_param_limits(params_low, params_high)
	#print 'n_params', n_params
	#lib.check(n_params)
	
	lib.LogLike.restype = c_double
	
	pymultinest.run(lib.LogLike, None, n_params, outputfiles_basename='%d-' % n_planets, resume = True, verbose = False, sampling_efficiency = 'model', n_live_points = 2000, max_modes=1000, evidence_tolerance = 2)
	
	# lets analyse the results
	a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='%d-' % n_planets)
	s = a.get_stats()

	import json
	json.dump(parameters, file('%sparams.json' % a.outputfiles_basename, 'w'), indent=2)
	json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
	#plot(a, parameters)
	return ( s['global evidence'], s['global evidence error'] )

def evidence(diff):
	if diff > math.log(100):
		return 'decisive'
	if diff > math.log(30):
		return 'very strong'
	if diff > math.log(10):
		return 'strong'
	return 'indifferent'

# number of planets
n = 0 # start value, we'll increase it later

ev, everr = run(0)
while True:
	n = n + 1
	evnext, evnexterr = run(n)
	print evnext, evnexterr, ev, everr
	if evnext - evnexterr > ev + everr + math.log(10):
		print "%d planets preferred over %d" % (n, n - 1)
		print "evidence difference: %s" % evidence((evnext - evnexterr - (ev + everr)))
	elif evnext + evnexterr > ev - everr + math.log(10):
		print "%d planned could be preferred over %d if accuracy was better" % (n, n - 1)
	else:
		print "stick with %d planets, %d are not preferred" % (n - 1, n)
		break
	ev, everr = evnext, evnexterr



