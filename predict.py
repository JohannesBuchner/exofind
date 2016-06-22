import numpy
from numpy import pi, log, exp, nan, sin, cos, arccos
from ctypes import *
lib = cdll.LoadLibrary('./exo_noprior.so')
lib.LogLike.restype = c_double
lib.ModLogPrior.restype = c_double

n_params = None
params_low = None
params_high = None
cube = None

kavg = 0
kmax = 1000
sturn = 0.1
smax = 10
pmax = 100
fs = 0.01
pmin = 0.01

# number of dimensions our problem has
parameters = ['V'] 
# plus possibly a noise term 
parameters += ['s']

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
			params_low[i], params_high[i] = 0., 2 * pi
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

if __name__ == '__main__':
	x = numpy.linspace(0, 40, 1000)
	# parameters: 
	# V: Velocity of the entire stellar system
	# s: noise term (not doing anything for predict)
	# parameters for each added planet:
	# P: Period (days)
	# K: Amplitude (m/s)
	# chi: anomaly parameter
	# e: eccentricity of orbit
	# omega: starting point (longitude of periastron)
	import matplotlib.pyplot as plt
	x, y = predict([0, 0.1, 10, 2, 0, 0, 0], x)
	plt.plot(x, y, ':')
	x, y = predict([0, 0.1, 10, 2, 0.2, 0.3, 0.1], x)
	plt.plot(x, y)
	plt.ylabel('RV Line Shift [m/s]')
	plt.xlabel('time [d]')
	plt.savefig('predict.pdf', bbox_inches='tight')
	plt.close()

