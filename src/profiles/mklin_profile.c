#include "eigen.h"

double bump_function(double rval);

double sigma_func(double x) {
	double sigfac = h0 /(2 * M_PI * bump_function(2.0));
	return sigfac * bump_function(x) * pow(x,-(1.5 + .5*temp_index));
}

double dlogsigma_func(double x) {
	return sigma_index;
}

double d2logsigma_func(double x) {
	return 0;
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


double bump_function(double rval) {
	
	
	double delta1 = 5.0 * h0 * pow(1.0,1.5-.5*temp_index);
	double delta2 = 5.0 * h0 * pow(2.0,1.5-.5*temp_index);

	double fac1 = .5*(1-.1)*(1+tanh((rval-1.0)/delta1))+.1;
	double fac2 = .5*(1-.1)*(1-tanh((rval-2.0)/delta2))+.1;
	
	return (fac1*fac2);

}

int analytic_potential(void) {
	return 0;
}

double omega_prec_grav_analytic(double x) {
	return 0;
}