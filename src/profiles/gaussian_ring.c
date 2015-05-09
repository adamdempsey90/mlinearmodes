#include "eigen.h"

double sig_param;

double sigma_func(double x) {
	sig_param = Mdisk/(pow(2*M_PI,1.5)*sigma_index*exp(-.5*sigma_index*sigma_index));
	return sig_param/x * exp(-log(x)*log(x)/(2*sigma_index*sigma_index));
}

double dlogsigma_func(double x) {
	return -1 - log(x)/(sigma_index*sigma_index);
}

double d2logsigma_func(double x) {
	return -1.0/(sigma_index*sigma_index);
}


double temp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = temp_index*pow(sigma_func(x),temp_index-1);
#else	
	res = h0*h0*pow(x,temp_index);
#endif
	return res;
}

double dlogtemp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = (temp_index - 1)*dlogsigma_func(x);
#else
	res = temp_index;
#endif
	return res;
}

double d2logtemp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = (temp_index - 1)*d2logsigma_func(x);
#else
	res = 0;
#endif
	return res;
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
	double res;
#ifdef POLYTROPE
	res = sqrt(temp_index)*pow(sigma_func(x),.5*(temp_index-1))*pow(omk_func(x),-1);
#else
	res =  h0*x*pow(x,flare_index);
#endif
	return res;
}

int analytic_potential(void) {
	return 0;
}

double omega_prec_grav_analytic(double x) {
	return 0;
}