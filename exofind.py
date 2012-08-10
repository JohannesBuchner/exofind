import pymultinest
import math
import os, sys
import numpy
import json

from ctypes import cdll
#cdll.LoadLibrary('libnest3.so')
lib = cdll.LoadLibrary('exo.so')
from ctypes import *

# put data into double[n][3] array
filename = sys.argv[1]
data = numpy.loadtxt(filename)
cdata = ((c_double * 3) * data.shape[0])()
for i in range(data.shape[0]):
	for j in range(3):
		cdata[i][j] = data[i,j]
lib.set_data(cdata, data.shape[0])
if 'MODE' in os.environ:
	fastmode = 'fast' in os.environ['MODE']
	plotmode = 'noplot' not in os.environ['MODE']
else:
	fastmode = False
	plotmode = True

# simplify parameter space:
# strongest impact is at most the variation of the data
kmax = (data[:,1].max() - data[:,1].min())
# sampling theorem -- highest frequency = longest period = 1 / dataduration / 2.
pmax = (data[:,0].max() - data[:,0].min()) / 2

# exclude very bad fits (worse than 5 sigma on all datapoints)
log_bad_fit = (-0.5 * 3**2.) * len(data)
log_very_bad_fit = (-0.5 * 5**2.) * len(data)

def plot_marg(parameters, values, i, s, p):
	import matplotlib.pyplot as plt
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
	newax.set_yticks([])
	newax.set_ylabel("Probability")
	ylim = oldax.get_ylim()
	newax.set_xlim(m['5sigma'])
	oldax.set_xlim(m['5sigma'])
	#plt.close()


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
		plot_marg(parameters, values, i, s, p)
		for j in range(i):
			plt.subplot(n_params, n_params, n_params * j + i + 1)
			p.plot_conditional(i, j, bins=20, cmap = plt.cm.gray_r)
			plt.xlabel(parameters[i])
			plt.ylabel(parameters[j])
			#plt.savefig('cond_%s_%s.pdf' % (params[i], params[j]), bbox_tight=True)
			#plt.close()

	plt.savefig(prefix + 'cond.pdf')
	plt.close()
	plt.figure(figsize=(4*4, 2*math.ceil(n_params / 4.)))
	for i in range(n_params):
		plt.subplot(math.ceil(n_params / 4.), 4, i + 1)
		plot_marg(parameters, values, i, s, p)
	plt.savefig(prefix + 'marg.pdf', bbox_inches='tight')
	plt.close()

def mkdir(path):
	if not os.path.exists(path): os.mkdir(path)

def run(n_planets, previous_periods = [], log_zero = -1e90):
	# number of dimensions our problem has
	parameters = ['V'] 
	# plus possibly a noise term 
	parameters += ['s']
	for i in range(n_planets):
		parameters += ['%s%d' % (var, i) for var in ['P', 'K', 'chi', 'e', 'omega']]
	n_params = len(parameters)
	params_low = (c_double * n_params)()
	params_high = (c_double * n_params)()
	simplification = 0.
	wrapped_params = []
	for i, p in enumerate(parameters):
		if p.startswith('V'):
			params_low[i], params_high[i] = -2000., +2000.
			wrapped_params.append(0)
		if p.startswith('P'):
			if not fastmode or len(previous_periods) == 0:
				params_low[i], params_high[i] = 0.2, pmax
			else:
				params_low[i], params_high[i] = previous_periods.pop()
				if params_low[i] < 0.2: params_low[i] = 0.2
				if params_high[i] > pmax: params_high[i] = pmax
				print 'fixing period (parameter %s) to %f..%f' % (p, params_low[i], params_high[i])
				simplification += math.log((math.log(params_high[i]) - math.log(params_low[i])) / (math.log(pmax) - math.log(0.2)))
			wrapped_params.append(0)
		if p.startswith('K'):
			params_low[i], params_high[i] = 1., kmax
			wrapped_params.append(0)
		if p.startswith('chi'):
			params_low[i], params_high[i] = 0., 1.
			wrapped_params.append(1)
		if p.startswith('e'):
			params_low[i], params_high[i] = 0., 1.
			wrapped_params.append(0)
		if p.startswith('omega'):
			params_low[i], params_high[i] = 0., 2 * math.pi
			wrapped_params.append(1)
		if p.startswith('s'):
			params_low[i], params_high[i] = 0.1, 100
			wrapped_params.append(0)
	lib.set_param_limits(params_low, params_high)
	#print 'n_params', n_params
	#lib.check(n_params)
	
	lib.LogLike.restype = c_double
	
	mkdir(filename + '.out')
	basename = filename + '.out/%d/' % n_planets
	mkdir(basename)
	if simplification != 0:
		print 'simplification: %f' % simplification
	pymultinest.run(lib.LogLike, None, n_params, resume = True, verbose = True, 
		outputfiles_basename=basename, sampling_efficiency = 'model', 
		n_live_points = 2000, evidence_tolerance = math.log(5),
		wrapped_params = wrapped_params, const_efficiency_mode = fastmode,
		null_log_evidence = log_zero,
		max_iter = 1000000,
		) # max_modes=1000, 
	
	# lets analyse the results
	a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename=basename)
	s = a.get_stats()

	json.dump(parameters, file('%sparams.json' % a.outputfiles_basename, 'w'), indent=2)
	json.dump(s, file('%s.json' % a.outputfiles_basename, 'w'), indent=2)
	if plotmode:
		print 'plotting ...'
		plot(a, parameters)
		print 'plotting done.'
	print 'evidence: %.1f +- %.1f' % (s['global evidence'] + simplification, s['global evidence error'])
	return ( s['global evidence'] + simplification, s['global evidence error'] )

def get_previous_periods(n_planets):
	if n_planets == 0:
		return []
	basename = filename + '.out/%d/' % n_planets
	parameters = json.load(file('%sparams.json' % basename))
	a = pymultinest.Analyzer(n_params = len(parameters), outputfiles_basename=basename)
	s = a.get_stats()
	return map(lambda (p,m): (m['5sigma'][0] - (m['5sigma'][1] - m['5sigma'][0]), m['5sigma'][1] + (m['5sigma'][1] - m['5sigma'][0])),
		filter(lambda (p,m): p.startswith('P'), 
		zip(parameters, s['marginals'])))
def get_previous_best_fit(n_planets):
	basename = filename + '.out/%d/' % n_planets
	s = numpy.loadtxt("%s.txt" % basename)
	return (-0.5*s[:,1]).max()
	
	
def evidence(diff):
	if diff > math.log(100):
		return 'decisive'
	if diff > math.log(30):
		return 'very strong'
	if diff > math.log(10):
		return 'strong'
	return 'indifferent'


if len(sys.argv) == 3:
	n = int(sys.argv[2])
	run(n)
	sys.exit(0)

# number of planets
n = 0 # start value, we'll increase it later

ev, everr = run(0)
while True:
	previous_periods = get_previous_periods(n)
	previous_best_fit = get_previous_best_fit(n)
	log_zero = log_very_bad_fit
	print 'using log_zero %.1f' % log_zero
	n = n + 1
	evnext, evnexterr = run(n, previous_periods, log_zero = log_zero)
	if evnext - evnexterr > ev + everr + math.log(10):
		print "%d planets preferred over %d" % (n, n - 1)
		print "evidence difference: %s" % evidence((evnext - evnexterr - (ev + everr)))
	elif evnext + evnexterr > ev - everr + math.log(10):
		print "%d planned could be preferred over %d if accuracy was better" % (n, n - 1)
	else:
		print "stick with %d planets, %d are not preferred" % (n - 1, n)
		break
	ev, everr = evnext, evnexterr

# write out marginal posteriors in a nice way
#   calculate pdf on a_s * sin i = K*T*(1-e^2)**0.5 / 2 / pi
#   calculate pdf on m_p * sin i = K*mstar^(2/3)*T^(1/3)*(1-e^2)**0.5 / 2 / pi / G
#   calculate pdf on a = mstar^(1/3) * T^(2/3) * G

# ?  calculate pdf of m_p / mstar using an assumption on sin i

# write out nice html / latex table of discovered planet system

# go through all params with associated prob. p(params)
# go through all times when observation could be made
# go through all observations that could be made at that time with the given params (proportional to its probability): p(y|params, time)
# search time that maximizes 
#    \sum{ log( p(y|params, time) ) * p(y|params, time) * p(params) }
#     - \sum{  log( p(y|time) ) * p(y|time)  dy  }





