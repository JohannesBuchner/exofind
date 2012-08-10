import pymultinest
import math
import os, sys
import numpy
import json
from random import Random
from time import time
import inspyred
prng = Random()
prng.seed(123)
#prng.seed(time())

from ctypes import *
#cdll.LoadLibrary('libnest3.so')
lib = cdll.LoadLibrary('exo.so')
lib.LogLike.restype = c_double

# put data into double[n][3] array
filename = sys.argv[1]
data = numpy.loadtxt(filename)
cdata = ((c_double * 3) * data.shape[0])()
for i in range(data.shape[0]):
	for j in range(3):
		cdata[i][j] = data[i,j]
lib.set_data(cdata, data.shape[0])

# simplify parameter space:
# strongest impact is at most the variation of the data
kmax = (data[:,1].max() - data[:,1].min())
print 'kmax: %f' % (kmax)
pmax = (data[:,0].max() - data[:,0].min()) * 2
# sampling theorem -- highest frequency --> longest period = 1 / dataduration / 2.
fs = 1 / (data[1:,0] - data[:-1,0]).min()
pmin = 1 / (fs / 2.)
print 'P limits: %f..%f' % (pmin, pmax)
#pmax = 15000

# exclude very bad fits (worse than 5 sigma on all datapoints)
log_bad_fit = (-0.5 * 3**2.) * len(data)
log_very_bad_fit = (-0.5 * 5**2.) * len(data)


def mkdir(path):
	if not os.path.exists(path): os.mkdir(path)

n_params = None
params_low = None
params_high = None
cube = None

# number of dimensions our problem has
parameters = ['V'] 
# plus possibly a noise term 
parameters += ['s']

def set_param_limits():
	global n_params
	global params_low
	global params_high
	global cube

	n_params = len(parameters)
	params_low = (c_double * n_params)()
	params_high = (c_double * n_params)()
	cube = (c_double * n_params)()
	
	for i, p in enumerate(parameters):
		if p.startswith('V'):
			params_low[i], params_high[i] = -2000., +2000.
		if p.startswith('P'):
			if True or len(previous_periods) == 0:
				params_low[i], params_high[i] = pmin, pmax
			else:
				params_low[i], params_high[i] = previous_periods.pop()
				if params_low[i] < 0.2: params_low[i] = 0.2
				if params_high[i] > pmax: params_high[i] = pmax
				print 'fixing period (parameter %s) to %f..%f' % (p, params_low[i], params_high[i])
		if p.startswith('K'):
			params_low[i], params_high[i] = 1., kmax
		if p.startswith('chi'):
			params_low[i], params_high[i] = 0., 1.
		if p.startswith('e'):
			params_low[i], params_high[i] = 0., 1.
		if p.startswith('omega'):
			params_low[i], params_high[i] = 0., 2 * math.pi
		if p.startswith('s'):
			params_low[i], params_high[i] = 0.1, 100
	lib.set_param_limits(params_low, params_high)

set_param_limits()
# later add more params	
n_planets = 0
#for i in range(n_planets):
#	parameters += ['%s%d' % (var, i) for var in ['P', 'K', 'chi', 'e', 'omega']]

def likelihood(*args):
	if len(args) == 1:
		args = args[0]
	if numpy.any(numpy.isnan(args)):
		return 1e300
	for i,a in enumerate(args):
		cube[i] = a
	v = lib.LogLike(cube, n_params, 0)
	#print args, v
	return -float(v)
def transform(*args):
	for i,a in enumerate(args):
		cube[i] = a
	lib.set_params(cube, len(args))
	return [float(cube[i]) for i in range(len(args))]

fn = {}
fn[2] = lambda c0, c1: likelihood(c0, c1)
fn[7] = lambda c0, c1, c2, c3, c4, c5, c6: likelihood(c0, c1, c2, c3, c4, c5, c6)
fn[12] = lambda c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11: likelihood(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)
fn[17] = lambda c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16: likelihood(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16)
fn[22] = lambda c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21: likelihood(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21)
fn[25] = lambda c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26: likelihood(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26)

# do minuit search
def run_minuit(start, errors):
	import minuit
	m = minuit.Minuit(fn[n_params])
	minuitparams = ['c%d' % i for i in range(len(start))]
	for s, e, paramname in zip(start, errors, minuitparams):
		m.limits[paramname] = (0,1)
		m.errors[paramname] = e
		m.values[paramname] = s
	
	m.maxcalls = 100000
	# minimization succeeds when edm is less than 0.001*tol*up
	#m.tol = 50
	m.printMode = 0
	m.up = 0.5
	m.hesse()
	#m.migrad()
	for pminuit, paramname, val in zip(minuitparams, parameters, transform(*map(m.values.get, minuitparams))):
		print "\t%s: %f  [+-%f]" % (paramname, val, m.errors[pminuit])
	#print 'values', m.values
	#print 'errors', m.errors
	optparams = map(m.values.get, minuitparams)
	errors = map(m.errors.get, minuitparams)
	return {'opt': optparams, 'errors': errors, 'minuit':m}

def run_fmin(start, errors):
	import scipy.optimize
	n = len(start)
	cons = [lambda x: x[i] for i in range(n)] + [lambda x: 1 - x[i] for i in range(n)]
	start = scipy.optimize.fmin_cobyla(func=likelihood, x0=start, cons=cons, disp=1, rhobeg = numpy.mean(errors))
	start = scipy.optimize.fmin(func=likelihood, x0=start, disp=1)
	return start

import matplotlib.pyplot as plt
plt.ion()
plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], label='data', ls=' ')
plt.show()

def run_de(seeds):
	class GlobalOpt(inspyred.benchmarks.Benchmark):
		def __init__(self):
			inspyred.benchmarks.Benchmark.__init__(self, n_params)
			self.bounder = inspyred.ec.Bounder(0, 1)
			self.maximize = True
		
		def generator(self, random, args):
			return [random.uniform(0,1) for _ in range(n_params)]
		
		def evaluator(self, candidates, args):
			fitness = []
			for c in candidates:
				fitness.append(-likelihood(*c))
			return fitness
	
	problem = GlobalOpt()
	
	#ea.terminator = inspyred.ec.terminators.diversity_termination
	def stats_obs(population, num_generations, num_evaluations, args):
	    stats = inspyred.ec.analysis.fitness_statistics(population)
	    worst_fit = '{0:>10}'.format(stats['worst'])[:10]
	    best_fit = '{0:>10}'.format(stats['best'])[:10]
	    avg_fit = '{0:>10}'.format(stats['mean'])[:10]
	    med_fit = '{0:>10}'.format(stats['median'])[:10]
	    std_fit = '{0:>10}'.format(stats['std'])[:10]

	    sys.stdout.write("{0:>10} {1:>10} {2:>10} {3:>10} {4:>10} {5:>10} {6:>10}\r".format(
	    	num_generations, num_evaluations, worst_fit, best_fit, med_fit, avg_fit, std_fit))
	    sys.stdout.flush()

	
	if n_params <= 2:
		timelimit = 2
	else:
		timelimit = 120 + 10 * n_params

	"""print 'SA:'
	ea = inspyred.ec.SA(prng)
	ea.terminator = inspyred.ec.terminators.evaluation_termination
	ea.observer = stats_obs
	final_pop = ea.evolve(generator=problem.generator, 
		evaluator=problem.evaluator, 
		seeds = seeds,
		bounder=problem.bounder,
		maximize=problem.maximize,
		max_evaluations=30000)
	
	seeds = [c.candidate for c in final_pop]
	print """ 
	ea = inspyred.ec.DEA(prng)
	ea.observer = stats_obs # inspyred.ec.observers.stats_observer
	#ea.observer = [inspyred.ec.observers.file_observer, inspyred.ec.observers.stats_observer]
	#ea.observer = inspyred.ec.observers.plot_observer
	ea.terminator = inspyred.ec.terminators.time_termination
	#ea.selector = inspyred.ec.selectors.fitness_proportionate_selection
	ea.selector = inspyred.ec.selectors.truncation_selection
	
	print 'EA:'
	final_pop = ea.evolve(generator=problem.generator, 
		evaluator=problem.evaluator, 
		pop_size=200,
		seeds=seeds,
		bounder=problem.bounder,
		maximize=problem.maximize,
		mutation_rate = 0.1,
		gaussian_stdev = 0.04,
		min_diversity = 0.1,
		max_time = timelimit,
		num_selected = 20,
		)
	print 
	#plt.close()
	return final_pop

def plot(args, label = ''):
	x = (c_double * 1000)(*numpy.linspace(data[:,0].min(), data[:,0].max(), 1000))
	y = (c_double * 1000)()
	for i,a in enumerate(args):
		cube[i] = a
	lib.predict(cube, len(cube), 0, 1000, x, y)
	plt.plot(x, y, '-', label=label)
	plt.legend(loc = 'upper left')
	#plt.savefig(label + '.pdf')
	plt.savefig(label + '.png')

def aicc(res):
	n = 2 * len(data)
	k = len(res['opt'])
	return (k - 2 * res['optvalue'] + 2 * k * (k + 1) / (n - k - 1.)) / 2.


# optimize globally first
pop = run_de(seeds = [[prng.uniform(0,1) for _ in range(n_params)] for _ in range(10000)])
print 'optimization found highest loglikelihood %.1f' % -pop[-1].fitness, pop[-1].candidate, numpy.std([c.candidate for c in pop[-10:]], axis=0)
plot(pop[-1].candidate, '%d planets' % n_planets)

start = pop[-1].candidate
errors = numpy.std([c.candidate for c in pop[-10:]], axis=0) 
start = run_fmin(start, errors)
res = run_minuit(start = start, errors = errors)
res['optvalue'] = -likelihood(*res['opt'])
del res['minuit']
allresults = {0:res}

while True:
	# add a planet
	parameters += ['%s%d' % (var, n_planets) for var in ['P', 'K', 'chi', 'e', 'omega']]
	n_planets += 1
	print '---- %d planets -------------------------' % n_planets

	set_param_limits()

	print 'creating seed population...'
	# create seed population from last population:
	seeds = []
	for c in pop[-10:]:
		seeds += [c.candidate + [prng.uniform(0,1) for _ in range(5)] for i in range(100000)]
	
	print 'running DE...'
	pop = run_de(seeds = [[prng.uniform(0,1) for _ in range(n_params)] for _ in range(10000)])
	print 'optimization found highest loglikelihood %.1f' % -pop[-1].fitness, pop[-1].candidate, numpy.std([c.candidate for c in pop[-10:]], axis=0)
	plot(pop[-1].candidate, '%d planets' % n_planets)

	res = run_minuit(start = pop[-1].candidate, errors = numpy.std([c.candidate for c in pop[-10:]], axis=0) )
	del res['minuit']
	res['optvalue'] = -likelihood(*res['opt'])
	print res
	allresults[n_planets] = res
	json.dump([allresults[i] for i in range(n_planets + 1)], 
		file('result.json', 'w'), indent=4)
		
	
	print 'aicc of %d planets: %.1f' % (n_planets, aicc(allresults[n_planets]))
	for k,v in allresults.iteritems():
		print '  %d: %f' % (k,v['optvalue'])
	lndiff = aicc(allresults[n_planets]) - aicc(allresults[n_planets - 1])
	print 'ln diff: %.1f' % lndiff
	if lndiff > 0:
		print 'only taking %d planets' % (n_planets - 1)
		break
	#if abs(lndiff) < numpy.log(100):
	#	break
	
	


