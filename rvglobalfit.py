import numpy
import math
import json
import sys
from random import Random
import copy
import inspyred
import time
import scipy.optimize, scipy.stats
from numpy import pi, log, exp, nan, sin, cos, arccos

prng = Random()
prng.seed(time.time())
#prng.seed(123)
#prng.seed(1234)

from ctypes import *
lib = cdll.LoadLibrary('exo.so')
lib.LogLike.restype = c_double
lib.ModLogPrior.restype = c_double

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
kavg = data[:,1].mean()
print 'kmax: %f' % (kmax)
pmax = (data[:,0].max() - data[:,0].min()) * 10
# sampling theorem -- highest frequency --> longest period = 1 / dataduration / 2.
fs = 1 / (data[1:,0] - data[:-1,0]).min()
pmin = 1 / (fs / 2.)
print 'P limits: %f..%f' % (pmin, pmax)
#pmax = 15000

# exclude very bad fits (worse than 5 sigma on all datapoints)
log_bad_fit = (-0.5 * 3**2.) * len(data)
log_good_fit = (-0.5 * 3**2.) * len(data)
log_very_bad_fit = (-0.5 * 5**2.) * len(data)

n_params = None
params_low = None
params_high = None
cube = None

n_data = len(data)

def prepare(params):
	global n_params
	global cube
	global params_low
	global params_high
	assert not numpy.any(numpy.isinf(params))
	assert not numpy.any(numpy.isnan(params))
	n_params = len(params)
	
	# number of dimensions our problem has
	global parameters
	parameters = ['V'] 
	# plus possibly a noise term 
	parameters += ['s']
	
	n_planets = (n_params - 2) / 5
	
	for i in range(n_planets):
		parameters += ['%s%d' % (var, i) for var in ['P', 'K', 'chi', 'e', 'omega']]
	
	params_low = (c_double * n_params)()
	params_high = (c_double * n_params)()
	cube = (c_double * n_params)()
	
	for i, p in enumerate(parameters):
		if p.startswith('V'):
			params_low[i], params_high[i] = kavg - 2 * kmax, kavg + 2 * kmax, 
		if p.startswith('P'):
			params_low[i], params_high[i] = pmin, pmax
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
	
	if numpy.any(numpy.isnan(params)):
		return 1e300
	for i,a in enumerate(params):
		cube[i] = a

def predict(c, x = None):
	prepare(c)
	if x is None:
		n = 1000
		x = (c_double * 1000)(*numpy.linspace(data[:,0].min(), data[:,0].max(), 1000))
	else:
		n = len(x)
		x = (c_double * len(x))(*x)
	y = (c_double * len(x))()
	for i,a in enumerate(c):
		cube[i] = a
	n_params = len(c)
	lib.predict(cube, n_params, 0, n, x, y)
	return numpy.asarray([x[i] for i in range(n)]), numpy.asarray([y[i] for i in range(n)])

def plot_predict_resid(newcube, predicted, resid):
	import matplotlib.pyplot as plt
	plt.figure()
	plt.subplot(2,1,1)
	lastplanet = newcube[-5:]
	if len(lastplanet) == 5:
		lastplanet[0] = 1. / lastplanet[0]
		lastplanet[1] = lastplanet[1]**2
		lastplanet[2] = (numpy.fmod(lastplanet[2] + 0.5, 1) - 0.5)  *2*pi
		lastplanet[4] = (numpy.fmod(lastplanet[4] + 0.5, 1) - 0.5)
	
		#plt.title(", ".join(["%.2f" % c for c in lastplanet]))
		plt.title("freq=%.5f pow=%.1f ang1=%.3f e=%.2f ang2=%.3f" % tuple(lastplanet), fontsize='small')
	elif len(lastplanet) == 2:
		plt.title("V=%.3f s=%.3f" % tuple(lastplanet))
	else:
		plt.title("V=%.3f" % tuple(lastplanet))
	plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], marker='x', ls=' ')
	plt.plot(data[:,0], predicted, 's ')
	plt.subplot(2,1,2)
	plt.plot(data[:,0], resid, 's ')

def plot_fft_results(freq, pow, cum, angles, x, y, freqi, powi, angi):
	import matplotlib.pyplot as plt
	plt.figure(figsize=(7,10))
	plt.subplot(5,1,1)
	plt.plot(freq, pow, 'x')
	plt.xlabel('frequency')
	plt.ylabel('power')
	plt.hlines(powi, plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,2)
	plt.plot(freq, cum, 'x')
	plt.ylabel('cumulative probability')
	plt.xlabel('frequency')
	plt.vlines(freqi, plt.ylim()[0], plt.ylim()[1], alpha=0.3)
	plt.subplot(5,1,3)
	plt.plot(cum, pow**0.5, 'x')
	plt.xlabel('cumulative probability')
	plt.ylabel('amplitude')
	plt.hlines(powi**0.5, plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,4)
	plt.plot(cum, angles, 'x')
	plt.xlabel('cumulative probability')
	plt.ylabel('phase')
	plt.hlines(angi, plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,5)
	plt.plot(x, y, 'x')
	plt.xlabel('time')
	plt.ylabel('value')
	xnew = numpy.linspace(x.min(), x.max(), 200)
	ynew = powi**0.5 * numpy.cos(xnew * freqi * pi * 2 + angi)
	plt.plot(xnew, ynew, '-', alpha=0.5, lw=3, 
		label="freq=%.5f\npow=%.3f\nang=%.3f" % (freqi, powi, angi))
	plt.legend(loc='best', prop=dict(size=10))


def fft_prior(cube, plot = False):
	import matplotlib.pyplot as plt
	import dftransform
	# K0 is gaussian with average + stdev of data
	bias = 0.
	print 'FFT transforming', cube
	
	newcube = []
	
	Vrv = scipy.stats.norm(numpy.mean(data[:,1]), numpy.std(data[:,1]))
	V = Vrv.ppf(cube[0])
	newcube.append( V )
	bias += log(Vrv.pdf(V))
	# would have been uniform prior, so don't need anything too special here
	
	ndim = len(cube)
	n_planets = (ndim - 1) / 5
	i0 = 1
	if ndim > n_planets * 5 + 1:
		# using s
		newcube.append( float(lib.ModLogPrior(c_double(cube[1]), c_double(params_low[1]), c_double(params_high[1]) )))
		# no bias here, same as before
		i0 = i0 + 1
	
	# now P, K, chi, e, omega
	# K * (sin(true_anomaly(ti + chi * P) + omega) + e * sin(omega))
	# K * (sin(true_anomaly(ti, P, e, chi) + omega + PI) + e * sin(omega + PI))
	# omega = phase
	# P = period
	# ti = 
	# chi = 
	# e = eccentricity: uniform prior
	for i in range(n_planets):
		for k in range(2, len(cube)):
			cube[k] = 0.5
	
		j0 = i * 5 + i0
		# calculate psd of residuals
		print 'calculating psd of residuals'
		x, preddata = predict(newcube, data[:,0])
		resdata = data[:,1] - preddata
		if plot:
			print 'plotting prediction...'
			plot_predict_resid(newcube, preddata, resdata)
			plt.savefig("predict_%02d.pdf" % i)
			print 'plotting prediction ... done'
		
		print 'dft prepare...'
		dftransform.prepare(data[:,0], resdata)
		dftransform.updateNyquist()
		dftransform.updateZeropoint()
		print 'dft calc...'
		out = dftransform.calc()
		print 'dft calc... done'
		freq, pow, angles = out[1:,0], out[1:,1], out[1:,2]
		
		low = pow.sum() / len(pow) / 100
		pdf = numpy.where(pow < low, low, pow)
		pdf = pdf / freq
		cum = pdf.cumsum()
		cum /= cum[-1]
		
		freqi = numpy.interp(x=cube[j0], xp=cum, fp=freq)
		powi = numpy.interp(x=cube[j0], xp=cum, fp=pow)
		angi = numpy.interp(x=cube[j0], xp=cum, fp=angles) + 3*pi/2.
		if plot:
			print 'plotting results ...'
			plot_fft_results(freq, pow, cum, angles, data[:,0], resdata, freqi, powi, angi)
			plt.savefig("fft_%02d.pdf" % i)
			print 'plotting results ... done'
		
		# powi**0.5 * numpy.sin(xnew * freqi * pi * 2 + angi)
		
		# generate distributions
		# draw from distributions using cube[] ppf
		# also compute probability there, vs. usual draw --> bias
		
		period = 1. / (freqi)
		bias -= scipy.stats.uniform.logpdf(log(period), log(pmin), log(pmax))
		
		amps = pow**0.5
		amprv = scipy.stats.norm(powi**0.5, powi**0.25)
		amplitude = amprv.ppf(cube[j0 + 1])
		bias += log( amprv.pdf(amplitude) )
		bias -= scipy.stats.uniform.logpdf(log(amplitude), log(1.), log(kmax))
		
		# make angle 3 sigma correspond to 2 pi
		chirv = scipy.stats.norm(angi / (2*pi), 1./3.)
		chi = numpy.fmod(chirv.ppf(cube[j0 + 2]) + 3., 1)
		bias += log( chirv.pdf(chi) )
		#print 'new angle: ', angi, (numpy.fmod(chi + 0.5, 1) - 0.5) * 2*pi
		if plot:
			print 'plotting results ...'
			plot_fft_results(freq, pow, cum, angles, data[:,0], resdata, 
				1 / period, amplitude**2, 
				(numpy.fmod(chi + 0.5, 1) - 0.5) * 2*pi)
			plt.savefig("fftrand_%02d.pdf" % i)
			print 'plotting results ... done'
		# would have been uniform prior, so don't need anything too special here
		
		# draw e, no bias there
		e = cube[j0 + 3]**2
		
		# calculate omega from chi and angle
		cosE = cos(chi / (1 + e))
		f0 = arccos( (cosE - e) / ( 1 - e * cosE) )
		omega0 = numpy.fmod(f0 - (chi) + 4*pi, 2*pi)
		omegarv = scipy.stats.norm(omega0, 1./3. * (2*numpy.pi))
		omega = numpy.fmod(omegarv.ppf(cube[j0 + 4]) + 3., 1)
		bias += log( omegarv.pdf(omega) )
		#chi = numpy.fmod(chi + 0.25, 1)
		omega = numpy.fmod(omega + 3*pi/2., 2*pi)
		
		newcube += [period, amplitude, chi, e, omega]
		# next round
		sys.exit(0)
	
	assert not numpy.any(numpy.isinf(newcube)), newcube
	assert not numpy.any(numpy.isnan(newcube)), newcube
	# so finally we have parameters and a parameter space deformation correction bias
	return bias, newcube

def calc(params):
	prepare(params)
	# transform cube, calculate residuals, etc.
	bias, fftparams = fft_prior(params, plot=True)
	prepare(fftparams)
	
	# AICc as fitness to be maximized
	v = lib.LogLike(cube, len(params), 0) # - bias
	#print params, [cube[i] for i in range(n_params)], v
	
	n = n_data
	k = len(params)
	assert not numpy.isinf(v)
	assert not numpy.isnan(v)
	aicc = (k - 2 * v + 2 * k * (k + 1) / (n - k - 1.)) / 2.
	assert not numpy.isnan(aicc)
	assert not numpy.isinf(aicc)
	return -aicc

"""
import sys

cmax = [prng.uniform(0,1) for _ in range(2)]
vmax = calc(cmax)

for i in range(10):
	improved = False
	for j in range(1000):
		c = [prng.uniform(0,1) for _ in range(2 + 5*i)]
		v = calc(c)
		if v > vmax:
			cmax = c
			vmax = v
			improved = True
			print 'found improvement, %d planets' % i
	if not improved:
		break

prepare(cmax)
fft_prior(cmax, plot=True)
sys.exit(0)
"""

def transform(args):
	prepare(args)
	return fft_prior(args, plot=False)[1]
	#lib.set_params(cube, len(args))
	#return [float(cube[i]) for i in range(len(args))]

class GlobalOpt(inspyred.benchmarks.Benchmark):
	def __init__(self):
		inspyred.benchmarks.Benchmark.__init__(self, None)
		self.bounder = inspyred.ec.Bounder(1e-9, 1-1e-9)
		self.maximize = True
	
	def generator(self, random, args):
		return [random.uniform(0,1) for _ in range(2)]
	
	def evaluator(self, candidates, args):
		fitness = []
		for c in candidates:
			assert not numpy.any(numpy.isnan(c))
			fitness.append(calc(c))
		return fitness

best = {}

def local_opt(candidate):
	n = len(candidate)
	cons = [lambda x: x[i] for i in range(n)] + [lambda x: 1 - x[i] for i in range(n)]
	start = list(scipy.optimize.fmin_cobyla(func=lambda x: -calc(x), x0=candidate, cons=cons, disp=1, rhobeg = 0.01, maxfun=2000))
	return start


# mutation that adds or removes a planet
def dimension_variator(random, candidates, args):
	mutants = []
	#print 'dimension_variator'
	b = args.get('add_remove_rate')
	a = args.get('add_rate')
	
	for c in candidates:
		#mutants.append(c)
		d = random.uniform(0,1)
		n_planets = (len(c) - 2) / 5
		if d < a: # rarely add
			# add planet
			#mutants += [ c + [random.uniform(0,1) for j in range(5)] for i in range(1000)]
			#print 'adding  planet: ', c, '-->'
			k = random.randint(0, n_planets) if n_planets > 0 else 0
			if k == n_planets:
				c += [random.uniform(0,1) for j in range(5)]
			else:
				for i in range(5):
					c[k*5+2 + i] = random.uniform(0,1)
			#c += [0.5 for j in range(5)]
			#c = local_opt(c)
			#print c
		elif d < b:
			# remove planet
			if n_planets > 0:
				k = random.randint(0, n_planets - 1)
				start = 2 +  k      * 5
				end   = 2 + (k + 1) * 5
				
				#print 'removing planet: ', c, '-->'
				c = c[:start] + c[end:]
		# else: remain
				#c = local_opt(c)
				#print c
				# mutants += [ [random.gauss(ci, 0.02) for ci in c] for i in range(1000)]
		else:
			# allow cycling for wrapped parameters chi and omega
			for i in range(n_planets):
				j = 2 + i * 5 + 2
				k = 2 + i * 5 + 4
				if c[j] > 0.95 and random.uniform(0,1) < 0.05:
					c[j] = 0
				elif c[j] < 0.05 and random.uniform(0,1) < 0.05:
					c[j] = 0
				if c[k] > 0.95 and random.uniform(0,1) < 0.05:
					c[k] = 1
				if c[k] < 0.05 and random.uniform(0,1) < 0.05:
					c[k] = 1
		
		mutants.append(c)
	#print 'dimension_variator done: ', mutants
	
	# reintroduce previous best ones from lower dimensions
	#lowdim = min([len(c) for c in candidates])
	#for bv in best.values():
	#	if len(bv.candidate) < lowdim:
	#		mutants.append(bv.candidate)
	
	return mutants

@inspyred.ec.variators.crossover
def variable_dim_crossover(random, mom, dad, args):
	mutants = []
	#print 'dimension_variator'
	a = args.get('crossover_rate')
	crossover_rate = args.setdefault('crossover_rate', 1.0)
	num_crossover_points = args.setdefault('num_crossover_points', 1)
	children = []
	if random.random() < crossover_rate:
		num_cuts = min(len(mom)-1, num_crossover_points)
		n = max(len(mom), len(dad))
		n_planets = (n - 2) / 5
		cut_points = random.sample([1] + [2 + i * 5 for i in range(n_planets)], num_cuts)
		cut_points.sort()
		#print 'crossing over: '
		#print '  ', dad
		#print '  ', mom, '--> to:'
		bro = copy.copy(dad)
		sis = copy.copy(mom)
		normal = True
		for i in range(n):
			if i in cut_points:
				normal = not normal
			if not normal:
				if len(bro) > i and len(mom) > i:
					bro[i] = mom[i]
				if len(sis) > i and len(dad) > i:
					sis[i] = dad[i]
		
		#print '  ', bro
		#print '  ', sis
		assert not numpy.any(numpy.isnan(bro))
		assert not numpy.any(numpy.isnan(sis))
		assert not numpy.any(numpy.isinf(bro))
		assert not numpy.any(numpy.isinf(sis))
		children.append(bro)
		children.append(sis)
	else:
		children.append(mom)
		children.append(dad)
	return children

# gaussian mutation -- should be adaptive
## start with 0.1 as stdev

def stats_obs(population, num_generations, num_evaluations, args):
	stats = inspyred.ec.analysis.fitness_statistics(population)
	worst_fit = '{0:>10}'.format(stats['worst'])[:10]
	best_fit = '{0:>10}'.format(stats['best'])[:10]
	avg_fit = '{0:>10}'.format(stats['mean'])[:10]
	med_fit = '{0:>10}'.format(stats['median'])[:10]
	std_fit = '{0:>10}'.format(stats['std'])[:10]
	
	dims = [(len(c.candidate) - 2) / 5 for c in population]
	
	#sys.stdout.write("G{0:>10} E{1:>10} W{2:>10} B{3:>10} M{4:>10} A{5:>10} S{6:>10} lD{7:>2} hD{8:>2} bD{9:>2}\r".format(
	#num_generations, num_evaluations, worst_fit, best_fit, med_fit, avg_fit, std_fit, min(dims), max(dims), (len(max(population).candidate) - 2)/5))
	best = max(population)
	values = transform(best.candidate)
	
	sys.stdout.write("G{0:>8} E{1:>8} D{7:>1}[{8:>1}-{9:>1}] W{2:>6} S{6:>6} B{3:>6} {10} \r".format(
		num_generations, num_evaluations, worst_fit, best_fit, med_fit, 
		avg_fit, std_fit, min(dims), max(dims), 
		(len(best.candidate) - 2)/5,
		",".join(["%s:%2.2f" % (p,f) for p,f in zip(parameters, values)]),
	))
	sys.stdout.flush()

import matplotlib.pyplot as plt
#plt.ion()

def progress_plot_obs(population, num_generations, num_evaluations, args):
    if num_generations % 100 != 0:
        return
    import pylab

    stats = inspyred.ec.analysis.fitness_statistics(population)
    best_fitness = -stats['best']
    worst_fitness = -stats['worst']
    median_fitness = -stats['median']
    average_fitness = -stats['mean']
    colors = ['black', 'blue', 'green', 'red']
    labels = ['average', 'median', 'best', 'worst']
    data = []
    pylab.figure('progress')
    if num_generations == 0:
        pylab.ion()
        data = [[num_evaluations], [average_fitness], [median_fitness], [best_fitness], [worst_fitness]]
        lines = []
        for i in range(4):
            line, = pylab.plot(data[0], data[i+1], color=colors[i], label=labels[i])
            lines.append(line)
        # Add the legend when the first data is added.
        pylab.legend(loc='lower right')
        args['plot_data'] = data
        args['plot_lines'] = lines
        pylab.xlabel('Evaluations')
        pylab.ylabel('Fitness')
    else:
        data = args['plot_data']
        data[0].append(num_evaluations)
        data[1].append(average_fitness)
        data[2].append(median_fitness)
        data[3].append(best_fitness)
        data[4].append(worst_fitness)
        lines = args['plot_lines']
        for i, line in enumerate(lines):
            line.set_xdata(numpy.array(data[0]))
            line.set_ydata(numpy.array(data[i+1]))
        args['plot_data'] = data
        args['plot_lines'] = lines
    ymin = min([min(d) for d in data[1:]])
    ymax = max([max(d) for d in data[1:]])
    yrange = ymax - ymin
    pylab.xlim((0, num_evaluations))
    pylab.ylim((1, 10000))
    pylab.gca().set_yscale('log')
    pylab.draw()
    pylab.show()

def plot_obs(population, num_generations, num_evaluations, args):
	
	if num_generations % 20 != 0:
		return
	#if num_evaluations > 200000:
	#	args['mutation_rate'] = 0.1
	#	args['gaussian_stdev'] = 0.01

	changed = False
	for c in population:
		if len(c.candidate) not in best or c.fitness > best[len(c.candidate)].fitness:
			best[len(c.candidate)] = c
			changed = True
	if not changed:
		return
	plt.figure()
	plt.clf()
	plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], label='data', ls=' ')
	
	for i in sorted(best.keys()):
		c = best[i]
		n = n_data
		k = len(c.candidate)
		bad_aicc = - (k - 2 * log_bad_fit + 2 * k * (k + 1) / (n - k - 1.)) / 2.
		bad = '[<%.0f]' % bad_aicc if c.fitness < bad_aicc else ''
		
		values = transform(c.candidate)
		
		label = '%d planets%s: aicc %.1f %s' % (
			(len(c.candidate) - 2)/5, 
			bad,
			c.fitness, 
			",".join(["%s:%2.2f" % (p,f) for p,f in zip(parameters, values)]))
		x, y = predict(fft_prior(c.candidate)[1])
		plt.plot(x, y, '-', label=label)
	plt.legend(loc = 'upper center', prop={"size":8}, 
		bbox_to_anchor=(0., 1.02, 1., .102), mode = 'expand', borderaxespad=0.)
	plt.savefig(filename + '_current_%05d.png' % (num_generations / 100))
	plt.savefig(filename + '_current.png')
	plt.close()
	
problem = GlobalOpt()

def good_val_termination(population, num_generations, num_evaluations, args):
	f = max(population).fitness
	n = n_data
	k = 2
	bad_aicc = - (k - 2 * log_good_fit + 2 * k * (k + 1) / (n - k - 1.)) / 2.
	return f > bad_aicc

"""
import logging
logger = logging.getLogger('inspyred.ec')
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('inspyred.log', mode='w')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)"""

ea = inspyred.ec.DEA(prng)
#ea = inspyred.ec.GA(prng)
#ea.selector = inspyred.ec.selectors.default_selection
#ea.selector = inspyred.ec.selectors.truncation_selection
ea.selector = inspyred.ec.selectors.fitness_proportionate_selection
#ea.replacer = inspyred.ec.replacers.generational_replacement
#ea.replacer = inspyred.ec.replacers.crowding_replacement
#ea.selector = inspyred.ec.selectors.tournament_selection
ea.replacer = inspyred.ec.replacers.truncation_replacement
#ea.replacer = [inspyred.ec.replacers.steady_state_replacement, inspyred.ec.replacers.truncation_replacement]

ea.variator = [variable_dim_crossover, dimension_variator, inspyred.ec.variators.gaussian_mutation]
ea.terminator = [inspyred.ec.terminators.evaluation_termination, inspyred.ec.terminators.time_termination] #, good_val_termination]

ea.observer = [stats_obs, plot_obs] #, progress_plot_obs]

seeds = None
#seeds = [[prng.uniform(0,1) for _ in range(2)] for _ in range(10000)]
#seeds += [[prng.uniform(0,1) for _ in range(2 + 5)] for _ in range(100)]
last_best = -1e100

print 'starting ...'
i = 1
while i < 8:
	final_pop = ea.evolve(generator=problem.generator, 
		evaluator=problem.evaluator, 
		seeds = seeds,
		pop_size=40,
		bounder=problem.bounder,
		maximize=problem.maximize,
		mutation_rate = 0.1,
		gaussian_stdev = 0.1 / i**2,
		max_evaluations=1000 / i,
		add_rate = 0.04,
		add_remove_rate = 0.04,
		num_selected = 20,
	#	num_selected = 100,
	#	tournament_size = 100,
	#	num_elites = 1,
		max_time = 10*60 # 10 min
		)
	print
	if max(final_pop).fitness < last_best + 10:
		i = i + 1
		print 'no significant improvements made (%.1f-->%.1f); cooling to %d' % (last_best, max(final_pop).fitness, i)
	else:
		i = i - 2
		print '   significant improvements made (%.1f-->%.1f); heating to %d' % (last_best, max(final_pop).fitness, i)
		if i < 1: 
			i = 1
	last_best = max(final_pop).fitness
	seeds = [c.candidate for c in final_pop]

final = max(final_pop)

print final.fitness, transform(final.candidate)

best_results = {}
for i in sorted(best.keys()):
	c = best[i]
	best_results[(i - 2) / 5] = {'candidate': c.candidate, 'values':transform(c.candidate), 'AICc': c.fitness}
import os
for i in range(100):
	if os.path.exists(filename + '.opt%02d.json' % i):
		continue
	json.dump({'all': best_results, 'best': {'candidate': final.candidate, 'values':transform(final.candidate), 'AICc': final.fitness}}, 
		file(filename + '.opt%02d.json' % i, 'w'), 
		indent=4)
	break

plt.figure('progress')
plt.savefig('progress_%02d.png' % i)
plt.close()

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



