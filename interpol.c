#include<stdio.h>
#include<math.h>

#include "interpol.h"

/* 
 * INTERPOLATION 
 * linear, square, and with derivatives
 */

static double a, b, c, d;
static double deltax;
static double xn;

/* derivative interpolation */

double interpol_diff_at(double x) {
	return a + 
		b * (x - xn) / deltax + 
		c * pow((x - xn) / deltax, 2) + 
		d * pow((x - xn) / deltax, 3);
}

double interpol_diff_at_normalized(double x) {
	return a + b * x + c * pow(x, 2) + d * pow(x, 3);
}

void interpol_diff_setup(double xn_in, double xm, double yn, double ym, double yyn, double yym) {
	double deltay  = ym  - yn;
	double deltayy = yym - yyn;
	xn = xn_in;
	deltax  = xm  - xn;
	
	a = yn;
	b = yyn * deltax;
	d = (yym + yyn) * deltax - 2 * deltay;
	c = (deltayy * deltax - 3 * d) / 2;
}

/* linear interpolation */

double interpol_lin_at(double x) {
	return a + 
		b * (x - xn) / deltax;
}

double interpol_lin_at_normalized(double x) {
	return a + b * x;
}


void interpol_lin_setup(double xn_in, double xm, double yn, double ym) {
	double deltay  = ym  - yn;
	xn = xn_in;
	deltax  = xm  - xn;
	
	a = yn;
	b = deltay;
}

/* square interpolation */

double interpol_sq_at(double x) {
	return a + 
		b * x + 
		c * pow(x, 2);
}

double interpol_sq_at_normalized(double x) {
	return a + 
		b * (x * deltax + xn) + 
		c * pow((x * deltax + xn), 2);
}

void interpol_sq_setup(double xn_in, double xm, double xo, double yn, double ym, double yo) {
	/*double deltay  = ym  - yn;*/
	double xn_in2 = xn = xn_in;
	double xo2 = pow(xo,2);
	double xm2 = pow(xm,2);
	double xn2 = pow(xn,2);
	
	double anom = xn*(xo2*ym-xm2*yo)+xn2*(xm*yo-xo*ym)+(xm2*xo-xm*xo2)*yn;
	double div = xn*(xo2-xm2)-xm*xo2+xm2*xo+xn2*(xm-xo);
	
	double bnom = xn2*(yo-ym)-xm2*yo+(xm2-xo2)*yn+xo2*ym;
	
	double cnom = xn*(yo-ym)-xm*yo+(xm-xo)*yn+xo*ym;

	(void)xn_in2;

	deltax = xo - xn;

	a = anom / div;
	b = -bnom / div;
	c = cnom / div;
}

/* END OF INTERPOLATION */

