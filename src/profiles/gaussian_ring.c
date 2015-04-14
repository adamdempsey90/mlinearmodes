#include "eigen.h"

double sig_param;

double sigma_func(double x) {
	sig_param = Mdisk/(pow(2*M_PI,1.5)*sigma_index*exp(-.5*sigma_index*sigma_index));
	
	return sig_param/x * exp(-log(x)*log(x)/(2*sigma_index));
}

double dlogsigma_func(double x) {
	return -1 - log(x)/(sigma_index*sigma_index);
}

double d2logsigma_func(double x) {
	return -1.0/(sigma_index*sigma_index);
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