
#include<stdio.h>
#include<stdlib.h>
#include <gsl/gsl_rng.h>

extern double LogLike(double *params, int ndim, int npars);

/*
int main(int argc, char ** argv) {
	int i;
	int j = 0;
	int n;
	double params[100];
	if (argc != 2) {
		fprintf(stderr, "SYNOPSIY: %s <nparams> < infile \n");
		return -127;
	}
	n = atoi(argv[1]);
	for (j = 0; !feof(stdin); j++) {
		for (i = 0; i < n + 3; i++) {
			if (fscanf(stdin, "%le", &params[i]) != 1) {
				printf("couldn't read token %i.\n" % i);
			}
		}
		
		LogLike(params, n, 0);
	}
	printf("handled %i lines\n" % j);
	
	return 0;
}
*/

extern double * params_low;
extern double * params_high;
extern double * data;
extern unsigned int n_data;

double low[100] = {-2000, 0.2, 1, 0, 0, 0, 0.2, 1, 0, 0, 0};
double high[100] = {2000, 10000, 20000, 1, 1, 1, 10000, 20000, 1, 1, 1};
double newdata[1000] = {
	1400.152000, -25.300000, 2.100000,
	1413.877490, -18.200000, 2.800000,
	1434.465730, -22.500000, 1.400000,
	1455.053960, -26.200000, 2.000000,
	1489.367690, -52.700000, 2.100000,
	1496.230430, -50.800000, 2.400000,
	1503.093180, -51.200000, 1.700000,
	1564.857880, -45.300000, 1.700000,
	1571.720630, -46.800000, 2.000000,
	1585.446120, -50.700000, 1.600000,
	1592.308860, -44.700000, 1.900000,
	1612.897100, -14.700000, 1.900000,
	1619.759840, -3.070000, 1.900000,
	1654.073570, 45.300000, 2.400000,
	1688.387290, 50.300000, 1.600000,
	1695.250040, 54.900000, 1.700000,
	1715.838270, 75.300000, 2.900000,
	1722.701020, 78.000000, 2.700000,
	1736.426510, 78.000000, 1.800000,
	1743.289250, 81.100000, 2.600000,
	1770.740240, 65.500000, 2.300000,
	1777.602980, 57.800000, 2.000000,
	1784.465730, 55.200000, 2.700000,
	1798.191220, 53.400000, 2.400000,
	1811.916710, 56.000000, 1.900000,
	1818.779450, 58.800000, 2.200000,
	1832.504940, 58.200000, 1.700000,
	1901.132390, 25.900000, 1.600000,
	1907.995140, 29.600000, 2.000000,
	1921.720630, 33.700000, 2.200000,
	1928.583370, 40.100000, 3.700000,
	1935.446120, 36.800000, 2.900000,
	1962.897100, 22.400000, 2.200000,
	1969.759840, 16.300000, 2.200000
};

void setup_data() {
	params_low = low;
	params_high = high;
	data = newdata;
	n_data = 34;
}

extern void check(int ndim);

#define IFDEBUG if(0)

int main(int argc, char ** argv) {
	const gsl_rng_type * T;
	gsl_rng * r;
	int n = 1000000;
	int i;
	int j;
	double params[100];
	double v;
	int n_params;
	if (argc != 2) {
		fprintf(stderr, "SYNOPSIY: %s <nparams> < infile \n", argv[0]);
		return -127;
	}
	n_params = atoi(argv[1]);

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	setup_data();

	for (i = 0; i < n; i++) 
	{
		IFDEBUG printf(" iteration %d\n", i);
		for (j = 0; j < n_params; j++) {
			params[j] = gsl_rng_uniform (r);
			IFDEBUG printf(" %f\t", params[j]);
		}
		IFDEBUG check(n_params);
		IFDEBUG printf("   calc ... \n");
		v = LogLike(params, n_params, 0);
		IFDEBUG printf("   calc ... %e\n", v);
	}

	gsl_rng_free (r);
}

