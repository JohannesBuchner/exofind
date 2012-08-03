#include<stdio.h>
#include<math.h>
#include<assert.h>

#include "interpol.h"
#include "rungekutta.h"

double p[] = {
#include "p.dat"
};

#ifndef NEVAL
#define NEVAL 10
#endif

double P;
double e;
double chi;
double A;

double f(double t, double theta) {
	/*printf("A = %f, e = %f, theta = %f, P = %f, chi = %f\n", A, e, theta, P, chi);*/
	return A * (1 - e * cos(theta + P * chi));
}

/*double function(double x) {
	return sin(x) * (pow(x, 2)/10 + 1);
}
double derivative(double x) {
	return (x * sin(x)) / 5 + (pow(x,2) / 10 + 1) * cos(x);
}*/

int main(void) {
	
	double minx;
	double maxx;
	double v_x[NPOINTS];
	double v_y[NPOINTS];
	double v_yy[NPOINTS];

	double x;
	double y;
	
	int i;
	int j;
	int k;
	double h;

	double theta_0 = 0.1;
	minx = theta_0;
	e = 0.5;
	chi = 0.1;
	P = 1;
	h = P / NPOINTS;
	maxx = minx + 2 * P;

	for(k = 0; k < 100000; k++) {
		e *= 0.9;
		chi *= 0.9;
		P *= 1.01;
		A = 2 * 3.14 / (P * pow(1 - pow(e, 2), 3./2));
		rungekutta(v_x, v_y, v_yy, NPOINTS, (maxx - minx)/NPOINTS, minx, theta_0);
		for(i = 0; i < NPOINTS; i++) {
			printf("%f\tnan\t%f\n", v_x[i], v_y[i]);
		}

		for(i = 0; i < NEVAL; i++) {
			x = p[i] / 20 + 1;
			if (x < minx || x > maxx) {
				fprintf(stderr, "skipped %f\n", x);
				continue;
			}
			j = floor(NPOINTS * (x - minx)/(maxx - minx));

			#if INTERPOL == 1
				interpol_lin_setup(v_x[j], v_x[j+1], v_y[j], v_y[j+1]);
				y = interpol_lin_at(x);
			#elif INTERPOL == 2
				if (j+2 >= NPOINTS)
					j--;
				interpol_sq_setup(v_x[j], v_x[j+1], v_x[j+2], v_y[j], v_y[j+1], v_y[j+2]);
				y = interpol_sq_at(x);
			#elif INTERPOL == 3
				interpol_diff_setup(v_x[j], v_x[j+1], v_y[j], v_y[j+1], v_yy[j], v_yy[j+1]);
				y = interpol_diff_at(x);
			#else
				#error "INTERPOL not in [1 (lin), 2 (square) or 3 (diff)]"
			#endif
		
			printf("%f\t%f\n", x, y);
		}
	}

	return 0;
}

