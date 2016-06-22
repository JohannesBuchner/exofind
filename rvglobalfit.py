import numpy
import math
import json
import sys
from random import Random
import copy
import evolve
import inspyred
import time
import numpy, scipy, scipy.interpolate, scipy.stats
from numpy import pi, log, exp, nan, sin, cos, arccos
import matplotlib.pyplot as plt
import dftransform
import shutil
plt.rcParams['pdf.compression'] = 0
prng = Random()
prng.seed(time.time())
#prng.seed(123)
#prng.seed(1234)

from ctypes import *
lib = cdll.LoadLibrary('exo_noprior.so')
lib.LogLike.restype = c_double
lib.ModLogPrior.restype = c_double

# put data into double[n][3] array
evolve.parse_args()
filename = evolve.args.data
prefix = filename + '_out_'

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
sturn = data[:,2].min()
smax = (data[:,1].max() - data[:,1].min())
print 's limits: %f..%f' % (sturn, smax)
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

def update_n_params(new_n_params):
	global n_params
	global cube
	global params_low
	global params_high
	n_params = new_n_params
	
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
			params_low[i], params_high[i] = sturn, smax
	lib.set_param_limits(params_low, params_high)

def prepare(params):
	global n_params
	global cube
	global params_low
	global params_high
	if len(params) != n_params:
		update_n_params(len(params))
	assert not numpy.any(numpy.isinf(params))
	assert not numpy.any(numpy.isnan(params))
	for i,a in enumerate(params):
		cube[i] = a

def predict(c, x = None):
	if x is None:
		n = 1000
		x = (c_double * 1000)(*numpy.linspace(data[:,0].min(), data[:,0].max(), 1000))
	else:
		n = len(x)
		x = (c_double * len(x))(*x)
	y = (c_double * len(x))()
	prepare(c)
	global cube
	global n_params
	lib.predict(cube, n_params, 0, n, x, y)
	return numpy.asarray([x[i] for i in range(n)]), numpy.asarray([y[i] for i in range(n)])

def plot_predict_resid(newcube, predicted, resid):
	plt.figure("predict", figsize=(7,10))
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
	#plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], marker='x', ls=' ')
	plt.plot(data[:,0], predicted, 's ')
	plt.subplot(2,1,2)
	plt.plot(data[:,0], resid, 's ')

def plot_fft_results(freq, pow, cum, angles, x, y, freqi, powi, angi):
	plt.figure(figsize=(10,15))
	plt.subplot(5,1,1)
	plt.plot(freq, pow, 'x')
	plt.gca().set_xscale('log')
	plt.xlabel('frequency')
	plt.ylabel('power')
	plt.hlines(powi, plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,2)
	plt.plot(freq, cum, 'x')
	plt.gca().set_xscale('log')
	plt.ylabel('cumulative probability')
	plt.xlabel('frequency')
	plt.vlines(freqi, 0, 1, alpha=0.3)
	plt.subplot(5,1,3)
	plt.plot(cum, pow**0.5, 'x')
	plt.xlabel('cumulative probability')
	plt.ylabel('amplitude')
	plt.hlines(powi**0.5, plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,4)
	plt.plot(cum, numpy.fmod(angles+10*numpy.pi, 2*numpy.pi), 'x')
	plt.xlabel('cumulative probability')
	plt.ylabel('phase')
	plt.hlines(numpy.fmod(angi+10*numpy.pi, 2*numpy.pi), plt.xlim()[0], plt.xlim()[1], alpha=0.3)
	plt.subplot(5,1,5)
	plt.plot(x, y, 'x')
	plt.xlabel('time')
	plt.ylabel('value')
	xnew = numpy.linspace(x.min(), x.max(), 200)
	ynew = powi**0.5 * numpy.cos(xnew * freqi * pi * 2 + angi)
	plt.plot(xnew, ynew, '-', alpha=0.5, lw=3, 
		label="freq=%.5f\npow=%.3f\nang=%.3f" % (freqi, powi, angi))
	plt.legend(loc='best', prop=dict(size=10))

Vrv = scipy.stats.norm(numpy.mean(data[:,1]), numpy.std(data[:,1]))

def fft_prior(cube, plot = False):
	global params_low
	global params_high
	global n_params
	if len(cube) != n_params:
		update_n_params(len(cube))
	assert len(cube) == n_params, [cube, n_params]
	# K0 is gaussian with average + stdev of data
	bias = 0.
	print 'FFT transforming', cube
	
	newcube = []
	lowprob = scipy.stats.norm().cdf(-5)
	V = Vrv.ppf(lowprob + (cube[0]) / (1 + 2 * lowprob))
	newcube.append(V)
	bias += log(Vrv.pdf(V))
	# would have been uniform prior, so don't need anything too special here
	
	ndim = len(cube)
	n_planets = (ndim - 1) / 5
	i0 = 1
	if ndim > n_planets * 5 + 1:
		# using s
		noise = float(lib.ModLogPrior(c_double(cube[1]), c_double(params_low[1]), c_double(params_high[1]) ))
		newcube.append(noise)
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
		#for k in range(2, len(cube)):
		#	cube[k] = 0.5
	
		j0 = i * 5 + i0
		# calculate psd of residuals
		#print 'calculating psd of residuals'
		x, preddata = predict(newcube, data[:,0])
		resdata = data[:,1] - preddata
		if plot:
			#print 'plotting prediction...'
			plt.figure("predict_%02d" % i)
			plot_predict_resid(newcube, preddata, resdata)
			plt.subplot(2,1,1)
			plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], label='data', marker='x', ls=' ', color='black')
			plt.savefig("predict_%02d_.pdf" % i)
			shutil.move("predict_%02d_.pdf" % i, "predict_%02d.pdf" % i)
			plt.close()
			#print 'plotting prediction ... done'
		
		#print 'dft prepare...'
		dftransform.prepare(data[:,0], resdata)
		dftransform.updateNyquist()
		dftransform.updateZeropoint()
		#print 'dft calc...'
		out = dftransform.calc()
		#print 'dft calc... done'
		freq, pow, angles = out[1:,0], out[1:,1], out[1:,2]
		
		low = pow.sum() / len(pow) / 20.
		quitelow = pow.max() / 4.
		pdf = numpy.where(pow < quitelow, low, pow)
		pdf = pdf / freq # weighting: prefer lower frequencies
		cum = pdf.cumsum()
		cum /= cum[-1]
		
		freqi = numpy.interp(x=cube[j0], xp=cum, fp=freq)
		powi = numpy.interp(x=cube[j0], xp=cum, fp=pow)
		angi = numpy.interp(x=cube[j0], xp=cum, fp=angles) + 3*pi/2.
		if plot:
			#print 'plotting results ...'
			plot_fft_results(freq, pow, cum, angles, data[:,0], resdata, freqi, powi, angi)
			plt.savefig("fft_%02d_.pdf" % i)
			shutil.move("fft_%02d_.pdf" % i, "fft_%02d.pdf" % i)
			plt.close()
			#print 'plotting results ... done'
		
		# powi**0.5 * numpy.sin(xnew * freqi * pi * 2 + angi)
		
		# generate distributions
		# draw from distributions using cube[] ppf
		# also compute probability there, vs. usual draw --> bias
		
		period = 1. / (freqi)
		bias -= scipy.stats.uniform.logpdf(log(period), log(pmin), log(pmax))
		
		amps = powi**0.5
		amprv = scipy.stats.norm(amps, amps / 5.) # so that at lowprob, amprv.ppf gives 0
		amplitude = amprv.ppf(lowprob + (cube[j0 + 1]) / (1 + 2 * lowprob))
		bias += log( amprv.pdf(amplitude) )
		bias -= scipy.stats.uniform.logpdf(log(amplitude), log(1.), log(kmax))
		
		# make angle 2 sigma correspond to 2 pi
		chirv = scipy.stats.norm(angi / (2*pi), 1./2.)
		chi = numpy.fmod(chirv.ppf(lowprob + (cube[j0 + 2]) / (1 + 2 * lowprob)) + 3., 1)
		bias += log( chirv.pdf(chi) )
		#print 'new angle: ', angi, (numpy.fmod(chi + 0.5, 1) - 0.5) * 2*pi
		if plot:
			#print 'plotting results ...'
			plot_fft_results(freq, pow, cum, angles, data[:,0], resdata, 
				1 / period, amplitude**2, 
				chi * 2*pi)
			#print 'plotting results ... saving'
			plt.savefig("fftrand_%02d_.pdf" % i)
			shutil.move("fftrand_%02d_.pdf" % i, "fftrand_%02d.pdf" % i)
			plt.close()
			#print 'plotting results ... done'
		# would have been uniform prior, so don't need anything too special here
		
		# draw e, no bias there
		e = cube[j0 + 3]**2
		
		# calculate omega from chi and angle
		cosE = cos(chi / (1 + e))
		f0 = arccos( (cosE - e) / ( 1 - e * cosE) )
		omega0 = numpy.fmod(f0 - (chi) + 4*pi, 2*pi)
		omegarv = scipy.stats.norm(omega0, 1./5. * (2*numpy.pi))
		omega = numpy.fmod(omegarv.ppf(lowprob + (cube[j0 + 4]) / (1 + 2 * lowprob)) + 5., 1)
		bias += log( omegarv.pdf(omega) )
		#chi = numpy.fmod(chi + 0.25, 1)
		omega = numpy.fmod(omega + 3*pi/2., 2*pi)
		
		newcube += [period, amplitude, chi, e, omega]
		print '  Period%2d   : %2.5f' % (i, period)
		print '  Amplitude%2d: %2.5f' % (i, amplitude)
		# next round
	print '  Velocity: %2.2f' % V
	print '  Noise   : %.3f' % noise
	#if plot:
	#	sys.stdout.write("press return: ")
	#	sys.stdout.flush()
	#	sys.stdin.readline()
	
	assert not numpy.any(numpy.isinf(newcube)), newcube
	assert not numpy.any(numpy.isnan(newcube)), newcube
	# so finally we have parameters and a parameter space deformation correction bias
	return bias, newcube

class Planets(evolve.AbstractComponents): 
	""" rv planet """
	def __init__(self, n):
		self.n = n
		self.component_length = 5 # ['P', 'K', 'chi', 'e', 'omega']
		self.head = 2 # V s
		self.length = self.head + self.component_length * n
		self.low_bounds = [0] * self.length
		self.high_bounds = [1] * self.length
	def decode(self, encoding, plot=False):
		bias, fftparams = fft_prior(encoding, plot=plot)
		return fftparams
	def visualize(self, encodings, limits):
		plt.figure("viz", figsize=(7,10))
		for i, encoding in enumerate(encodings[:len(encodings)/3+2]):
		#for i, encoding in enumerate(encodings[0:1]):
			alpha = 1 if i == 0 else 0.3
			fftparams = self.decode(encoding, plot = (i == 0))
			plt.figure("viz")
			x, predicted = predict(fftparams, data[:,0])
			resid = data[:,1] - predicted
			plt.subplot(3,1,1)
			plt.plot(data[:,0], predicted, 's ', alpha=alpha)
			plt.subplot(3,1,2)
			plt.plot(predicted, 's ', alpha=alpha)
			plt.subplot(3,1,3)
			plt.plot(resid, 's ', alpha=alpha)
		
		plt.subplot(3,1,1)
		plt.errorbar(x=data[:,0], y=data[:,1], yerr=data[:,2], label='data', marker='x', ls=' ', color='black')
		plt.subplot(3,1,2)
		plt.errorbar(x=range(len(data[:,1])), y=data[:,1], yerr=data[:,2], label='data', marker='x', ls=' ', color='black')
		plt.savefig(prefix + "_viz_.pdf")
		shutil.move(prefix + "_viz_.pdf", prefix + "_viz.pdf")
		plt.close()


class RVLikelihoodModel(object):
	def __init__(self, encoder, limits, events=[]):
		self.events = events
		self.encoder = encoder
		self.limits = limits
	
	def likelihood(self, candidate):
		fftparams = self.encoder.decode(candidate)
		prepare(fftparams)
		v = lib.LogLike(cube, len(cube), 0) # - bias
		#print params, [cube[i] for i in range(n_params)], v
		return v
	
	def generate_events(self, candidate, random, k):
		raise NotImplementedError

# mutation that adds or removes a planet
def rotation_variator(random, candidates, args):
	mutants = []
	#print 'dimension_variator'
	r = args.get('rotate_rate', 0.05)
	
	for c in candidates:
		n_planets = (len(c) - 2) / 5
		for i in range(n_planets):
			j = 2 + i * 5 + 2
			k = 2 + i * 5 + 4
			if c[j] > 0.95 and random.uniform(0,1) < r:
				print 'applying rotation', c[j]
				c[j] = 1 - c[j]
			elif c[j] < 0.05 and random.uniform(0,1) < r:
				print 'applying rotation', c[j]
				c[j] = 1 - c[j]
			if c[k] > 0.95 and random.uniform(0,1) < r:
				print 'applying rotation', c[k]
				c[k] = 1 - c[k]
			elif c[k] < 0.05 and random.uniform(0,1) < r:
				print 'applying rotation', c[k]
				c[k] = 1 - c[k]
		mutants.append(c)
	return mutants

if __name__ == '__main__':
	limits = [[data[:,0].min(), data[:,0].max()]]
	
	encoder = Planets(evolve.args.n_components)
	
	pm = RVLikelihoodModel(encoder=encoder, limits=limits)
	ea, evolve_args = evolve.setup(pm)
	# inject rotation variator
	ea.variator += [rotation_variator]
	
	gc_debug = False
	if gc_debug:
		from collections import defaultdict
		from gc import get_objects
		before = defaultdict(int)
		after = defaultdict(int)
		for i in get_objects():
			before[type(i)] += 1
		import cherrypy
		import dowser

		"""cherrypy.tree.mount(dowser.Root())
		cherrypy.config.update({
			'environment': 'embedded',
			'server.socket_port': 8081
			})
		cherrypy.engine.start()"""
	evolve.evolve(ea, evolve_args)
	if gc_debug:
		objs = get_objects()
		for i in objs:
			after[type(i)]+=1
		items = [(k,after[k]-before[k]) for k in after if after[k]-before[k]]
		for j, (k,v) in enumerate(sorted(items, key = lambda (k,v):v, reverse=True)):
			print '  ', k, v, [str(i) for i in evolve.random.sample(objs, 1000) if type(i) == k][:10]
			if j > 5: break
		del objs
		del items
	
		from guppy import hpy
		h = hpy()
		print h.heap()
		del h
	#cherrypy.engine.block()

