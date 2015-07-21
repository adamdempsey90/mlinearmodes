#include "eigen.h"


const static double h_p = .05;
const static double eps_p = .1;
const static double delta_R = 5.0;
const static double R1 = 1;
const static double R2 = 2;

#define s_p  (1 + flare_index)
#define delta1 (delta_R * h_p * pow(R1,flare_index + 1))
#define  delta2  (delta_R * h_p * pow(R2,flare_index + 1))

double bump_function(double x);
double f1_func(double x);
double drf1_func(double x);
double d2rf1_func(double x);
double f2_func(double x);
double drf2_func(double x);
double d2rf2_func(double x);


double sech(double x) {
	return 1/cosh(x);
}

double sigma_func(double x) {
	double sigfac = h_p /(2 * M_PI * bump_function(2.0));
	return sigfac * bump_function(x) * pow(x,-s_p);
}

double dlogsigma_func(double x) {
	return -s_p + drf1_func(x)/f1_func(x) + drf2_func(x)/f2_func(x);
}


double d2logsigma_func(double x) {
	double d2lf1 = d2rf1_func(x) - drf1_func(x)*drf1_func(x)/f1_func(x);
	double d2lf2 = d2rf2_func(x) - drf2_func(x)*drf2_func(x)/f2_func(x);

	return d2lf1/f1_func(x) + d2lf2/f2_func(x);
}


double temp_func(double x) {
	return  h0*h0 * pow(x,temp_index);
}

double dlogtemp_func(double x) {
	return temp_index;
}

double d2logtemp_func(double x) {
	return 0;
}

double omk_func(double x) {
	return pow(x,-1.5);
}

double dlogomk_func(double x) {
	return -1.5;
}

double d2logomk_func(double x) {
	return 0;
}

double scaleH_func(double x) {
	return h0*x*pow(x,flare_index);
}


double f1_func(double x) {

	return .5*(1-eps_p)*(1+tanh((x-R1)/delta1))+eps_p;

}

double f2_func(double x) {

	return  .5*(1-eps_p)*(1-tanh((x-R2)/delta2))+eps_p;

}

double drf1_func(double x) {
	double arg = (x - R1)/delta1;

	return  -x * (eps_p - 1) * sech(arg)*sech(arg)/(2*delta1);

}
double drf2_func(double x) {
	double arg = (x - R2)/delta2;

	return  x * (eps_p - 1) * sech(arg)*sech(arg)/(2*delta2);

}

double d2rf1_func(double x) {
	double arg = (x - R1)/delta1;
	return x*(eps_p - 1)*sech(arg)*sech(arg)*(2*x*tanh(arg)-delta1)/(2*delta1*delta1);
}

double d2rf2_func(double x) {
	double arg = (x - R2)/delta2;
	return x*(eps_p - 1)*sech(arg)*sech(arg)*(-2*x*tanh(arg)+delta2)/(2*delta2*delta2);
}


double bump_function(double x) {

	return f1_func(x) * f2_func(x);
}

int analytic_potential(void) {
	return 0;
}

double omega_prec_grav_analytic(double x) {
	return 0;
}
