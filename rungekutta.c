#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double f(double t, double theta);


/*
#define RK4_DEBUG
*/
void rungekutta(double * t, double * thetas, double * diff, 
	int n, double h, double t0, double theta_0
) {
	double kn1, kn2, kn3, kn4;
	double theta_i;
	double t_i;

	int i;

	theta_i = thetas[0] = theta_0;
	#ifdef RK4_DEBUG
	printf("%d\t%f\t%f\n", 0, theta_i, t0);
	#endif

	for(i = 1; i < n; i++) {
		t_i = t[i - 1] = (i - 1) * h + t0;
		/* function evaluated at previous point */
		kn1 = f(t_i       , theta_i);
		kn2 = f(t_i + h/2., theta_i + h * kn1/2.);
		kn3 = f(t_i + h/2., theta_i + h * kn2/2.);
		kn4 = f(t_i + h   , theta_i + h * kn3);
		
		diff[i - 1] = kn1;
		
		theta_i += h/6. * (kn1 + 2*kn2 + 2*kn3 + kn4);
		thetas[i] = theta_i;
		#ifdef RK4_DEBUG
		printf("%d\t%f\t%f\t%f (above)\n", i, theta_i, i * h + t0, kn1);
		#endif

	}
	i = n - 1;
	t_i = t[i] = i * h + t0;
	diff[i] = f(t_i, theta_i);
	#ifdef RK4_DEBUG
	printf("%d\t%f\t%f\t%f\n", i, theta_i, t_i, diff[i]);
	printf("end of RK4\n");
	#endif
}


