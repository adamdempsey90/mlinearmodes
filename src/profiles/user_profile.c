#include "eigen.h"


// g = sigma0

double sigma_func(double x) {
	return pow(x, -1.5);
}

double dlogsigma_func(double x) {
	return -1.5;
}

double d2logsigma_func(double x) {
	return 0;
}


double temp_func(double x) {
	return h0*h0*pow(x,-.5);
}

double dlogtemp_func(double x) {
	return -.5;
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
	return h0*x*pow(x,.25);
}

int analytic_potential(void) {
	return 1;
}

double omega_prec_grav_analytic(double x) {
//	return -pow(x,1.5)*pow(x*x + 1,-2.5) ;
	return pow(scaleH_func(x)/x,2)*omk_func(x);
}
