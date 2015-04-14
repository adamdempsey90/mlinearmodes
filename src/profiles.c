#include "eigen.h"


double sig_param = 1; 

double sigma_potential_func(double x);
double dlogsigma_potential_func(double x);

double sigma_func(double x) {
	return Mdisk*sig_param*pow(x*x+sig_param*sig_param,-1.5) / (2*M_PI);
}

double dlogsigma_func(double x) {
	return -3*x*x/(x*x + sig_param*sig_param);
}

double d2logsigma_func(double x) {
	return -6*sig_param*sig_param*x*x * pow(sig_param*sig_param+x*x,-2);
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
	return 1;
}

double sigma_potential_func(double x) {
	return -Mdisk/sqrt(x*x + sig_param*sig_param);
}
double dlogsigma_potential_func(double x) {
	return Mdisk*x*x*pow(x*x + sig_param*sig_param,-1.5);
}
double omega_prec_grav_analytic(double x) {
	return -1.5*Mdisk*sig_param*sig_param*pow(x*x+sig_param*sig_param,-2.5)/omk_func(x);
}

