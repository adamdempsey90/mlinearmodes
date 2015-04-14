#include "eigen.h"

double rdecay = 1;
double p = 2;
double sigma_func(double x) {
	return sigma0 * pow(x, sigma_index) * exp(-pow(x/rdecay,p));
}

double dlogsigma_func(double x) {
	return sigma_index - p * pow(x/rdecay,p);
}

double d2logsigma_func(double x) {
	return -p*p*pow(x/rdecay,p);
}


double temp_func(double x) {
	return h0*h0*pow(x,temp_index);
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

int analytic_potential(void) {
	return 0;
}

double omega_prec_grav_analytic(double x) {
	return 0;
}