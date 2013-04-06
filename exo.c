#include <math.h>
#include <stdio.h>
#include "rungekutta.h"
#include "interpol.h"

#define MAX_NPOINTS 1024
#define MAX_PLANETS 20
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define NOISE

#define gsl_matrix_get(data, i, j) data[i*3 + j]
/*
#define USE_RK4
*/
struct planetorbit {
	double P;
	double e;
	double chi;
	double omega;
	double K;

	double B;
	double v_x[MAX_NPOINTS];
	double v_y[MAX_NPOINTS];
	double v_yy[MAX_NPOINTS];
	double minx;
	double maxx;
	double theta_0;
};

static struct planetorbit * current;

struct planetorbit planets[MAX_PLANETS]; 

double * params_low;
double * params_high;
double * data;
unsigned int n_data;

/* for RK4 */
double f(double t, double theta) {
	(void)t;
	return current->B * (1 - current->e * cos(theta + current->P * current->chi));
}

double UniformPrior(double c, double xmin, double xmax){
	return xmin + c * (xmax - xmin);
}
double LogPrior(double c, double xmin, double xmax){
	double log_min = log(xmin);
	double log_max = log(xmax);
	return exp(log_min + c * (log_max - log_min));
}
double ModLogPrior(double c, double xturn, double xmax) {
	return xturn * (exp( (c + 1e-10) * log(xmax/xturn + 1)) - 1);
}

void set_param_limits(double * new_params_low, double * new_params_high) {
	params_low = new_params_low;
	params_high = new_params_high;
}
void set_data(double * new_data, int new_n_data) {
	data = new_data;
	n_data = new_n_data;
}

/* 
 * estimate theta from t0 = gsl_matrix_get(m->data, 0, 0), P, e, chi
 */
static void theta_setup() {
	current->minx = 0;
	current->theta_0 = 0;
	current->maxx = current->minx + current->P;
	current->B = 2 * M_PI / (current->P * pow(1 - pow(current->e, 2), 3./2));
	rungekutta(current->v_x, current->v_y, current->v_yy, 
		n_data, (current->maxx - current->minx)/n_data, 
		current->minx, current->theta_0);
}

static double theta_eval(double x) {
	int j;
	double y;
	int n_points = n_data;
	
	if (x < current->minx || x > current->maxx)
		x = fmod(x, current->P);
	
	j = floor(n_points * (x - current->minx) / (current->maxx - current->minx));

	if (j + 1 > n_points - 1)
		j = n_points - 2;

	#if INTERPOL == 1
		interpol_lin_setup(
			current->v_x[j], current->v_x[j+1], 
			current->v_y[j], current->v_y[j+1]);
		y = interpol_lin_at(x);
	#elif INTERPOL == 2
		if (j+2 >= n_points)
			j--;
		interpol_sq_setup(
			current->v_x[j], current->v_x[j+1], current->v_x[j+2], 
			current->v_y[j], current->v_y[j+1], current->v_y[j+2]);
		y = interpol_sq_at(x);
	#elif INTERPOL == 3
		assert(j < n_points);
		assert(j+1 < n_points);
		interpol_diff_setup(
			current->v_x[j], current->v_x[j+1], 
			current->v_y[j],  current->v_y[j+1], 
			current->v_yy[j], current->v_yy[j+1]);
		y = interpol_diff_at(x);
	#else
		#error "INTERPOL not in [1 (lin), 2 (square) or 3 (diff)]"
	#endif
	return y;
}

void check(int ndim) {
	unsigned int i;
	printf(" PARAMETERS \n");
	for (i = 0; (signed)i < ndim; i++) {
		printf("  param %i/%i low, high: %f %f\n", i, ndim, params_low[i], params_high[i]);
	}
	printf(" DATA \n");
	for (i = 0; i < n_data; i++) {
		printf("  data point %i/%i: %f %f %f\n", i+1, n_data, gsl_matrix_get(data,i,0), gsl_matrix_get(data,i,1), gsl_matrix_get(data,i,2));
	}
}

const double D2S = 60*60*24;

// this function calculates eccentric anomaly
double ecc_anomaly(double time, double prd, double ecc, double peri_pass)
{
	time = time * D2S;
	double T = prd * D2S;
	double n = 2 * M_PI / T;
	double tau = peri_pass * T;
	double M = n * (time + tau);
	double E0; // start value
	double Mstar;// check equation  6.6.9 Danby
	double sigma;
	int  int_M2PI;
	// to get an interger value for M/(2*PI)
	// so that for M we have a value between 0 and 2PI
	int_M2PI = M / (2 * M_PI);
	Mstar = M - int_M2PI * 2 * M_PI;

	// define a SIGN function
	sigma = fabs(sin(Mstar));

	E0 = Mstar + sigma * 0.85 * ecc;
	// the value for k=0.85 is arbitrary
	// the only condition is 0<= k <=1 check, again Danby   

	double TINY = 1e-6;
	int count = 0;
	for(count = 0; count < 100; count++) {
		double eSinE = ecc * sin(E0);    // a dummy
		double f = E0 - eSinE - Mstar;
		if (fabs(f) < TINY)
			break;
		double eCosE = ecc * cos(E0);
		double f1 = 1 - eCosE;
		double f2 = eSinE;
		double f3 = eCosE;
		double dE0 = -f / f1;
		dE0 = -f / (f1 + 0.5 * dE0 * f2);
		dE0 = -f / (f1 + 0.5 * dE0 * f2 + dE0 * dE0 * f3 / 6);
		E0 = E0 + dE0;
	}

	return E0;
}

double true_anomaly(double time, double prd, double ecc, double peri_pass)
{
	double E = ecc_anomaly(time,prd,ecc,peri_pass);
	double f = acos((cos(E) - ecc) / (1 - ecc * cos(E)));
	//acos gives the principal values ie [0:PI]
	//when E goes above PI we need another condition
	if (E > M_PI)
		f = 2 * M_PI - f;
	return f;
}

double calc(unsigned int n_planets, double ti) {
	double fi = 0;
	unsigned int j;
	
	for(j = 0; j < n_planets; j++) {
		current = &planets[j];
#ifdef USE_RK4
		fi += current->K * (cos(theta_eval(ti + current->chi * current->P) 
			+ current->omega) + current->e * cos(current->omega));
#else
		fi += current->K * (sin(true_anomaly(ti, current->P, current->e, current->chi) 
			+ current->omega + M_PI) + current->e * sin(current->omega + M_PI));
#endif
	}
	
	return fi;
}

/* (9)
 * A not needed. just sum over the datapoints.
 */
unsigned int set_params(double *params, int ndim) {
	unsigned int i;
	unsigned int j;
	
	double s;
	
	unsigned int n_planets;
	i = 0;
	n_planets = (ndim - 1) / 5;

	/* V prior: uniform improper */
	/* (17) priors */
	#ifdef USE_PRIORS
	params[i] = UniformPrior(params[i], params_low[i], params_high[i]);
	#endif
	/*V     = params[0]; */
	i++;
	
	if (ndim > (signed) n_planets * 5 + 1) {
		#ifdef USE_PRIORS
		params[i] = ModLogPrior(params[i], params_low[i], params_high[i]);
		#endif
		s     = params[i++];
	} else {
		s     = 0;
	}
	
	for (j = 0; j < n_planets; j++) {
		current = &planets[j];

		/* P prior: jeffreys */
		#ifdef USE_PRIORS
		params[i] = LogPrior(params[i], params_low[i], params_high[i]);
		#endif
		current->P     = params[i++]; /* days */
		
		/* K prior: mod jeffreys */
		#ifdef USE_PRIORS
		params[i] = ModLogPrior(params[i], params_low[i], params_high[i]);
		#endif
		current->K     = params[i++];

		/* Chi prior */
		#ifdef USE_PRIORS
		params[i] = UniformPrior(params[i], params_low[i], params_high[i]);
		#endif
		current->chi   = params[i++];

		/* e */
		#ifdef USE_PRIORS
		params[i] = UniformPrior(params[i], params_low[i], params_high[i]);
		#endif
		current->e     = params[i++];
		/* omega */
		#ifdef USE_PRIORS
		params[i] = UniformPrior(params[i], params_low[i], params_high[i]);
		#endif
		current->omega = params[i++];

#ifdef USE_RK4
		theta_setup();
#endif
	}
	return n_planets;
}

int planets_are_sorted(unsigned int n_planets, int ndim) {
	unsigned int i;
	unsigned int j = 1;
	double lastp = planets[0].P;
	if (n_planets == 0)
		return 0;
	if (ndim > (signed) n_planets * 5 + 1) {
		j = 2;
	}
	for (i = 1; i < n_planets; i++) {
		/* avoid the planets that have been forced to a range */
		if (params_low[i * 5 + j] / params_high[i * 5 + j] < 100) {
			continue;
		}
		
		/* avoid double solutions, by making sure the periods 
		 * are sorted in increasing order */
		if (planets[i].P < lastp)
			return 1;
		lastp = planets[i].P;
	}
	return 0;
}

double LogLike(double *params, int ndim, int npars) {
	unsigned int n_planets = set_params(params, ndim);
	unsigned int i;
	
	double ti;
	double vi;
	double fi;
	double varsum;
	double vari;
	double prob = 0;
	double V = params[0];
	#ifdef NOISE
	double s = params[1];
	#endif
	
	/*if (planets_are_sorted(n_planets, ndim) != 0) 
		return -1e300;*/
	
	for (i = 0; i < n_data; i++) {
		ti   = gsl_matrix_get(data,i,0);
		vi   = gsl_matrix_get(data,i,1);
		vari = pow(gsl_matrix_get(data,i,2), 2);
		varsum = vari 
			#ifdef NOISE
				+ pow(s, 2)
			#endif
			;
		
		/* (5)  f_i = */
		fi = V + calc(n_planets, ti);
		
		/* (9)  2*p = */
		prob -= pow(vi - fi, 2) / varsum;
		
		/* (10) 2*A = */
		prob -= log(varsum);
	}
	
	prob += log(2 * M_PI) * n_data;
	/*for (i = 0; (signed)i < ndim; i++) {
		printf("  %.3f \t", params[i]);
	}
	printf(" --> %f\n", prob);*/
	
	return prob;
}

void predict(double *params, int ndim, int npars, int n_data, double * x, double * y) {
	unsigned int n_planets = set_params(params, ndim);
	unsigned int i;
	double V = params[0];
	
	for (i = 0; i < n_data; i++) {
		y[i] = V + calc(n_planets, x[i]);
	}
}

