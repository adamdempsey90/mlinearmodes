#include "eigen.h"

//#define ANALYTICPOTENTIAL

double width = 1;
double center = 1;
double amplitude = 5;
double Qmin = 2;

double h_p = .05;



double sigma_func(double x) {
	double sig_ref = temp_func(x)* pow(center,-sigma_index-1.5)/( M_PI*Qmin*(1+amplitude));
	double exp_fac = exp(-(x-center)*(x-center)/(2*width*width));
	return sig_ref * pow(x, sigma_index) * (1  + amplitude*exp_fac);
}

double dlogsigma_func(double x) {
	double exp_fac = exp((x-center)*(x-center)/(2*width*width));
	return sigma_index + (center-x)*x*amplitude/(width*width*(amplitude + exp_fac)); 
}

double d2logsigma_func(double x) {
	double result;
	double exp_fac = exp((x-center)*(x-center)/(2*width*width));
	result = (center-2*x)*amplitude*width*width + exp_fac*((center-x)*(center-x)*x+(center-2*x)*width*width);
	return result * x*amplitude * pow(width,-4)*pow(exp_fac+amplitude,-2);
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
