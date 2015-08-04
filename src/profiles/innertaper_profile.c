#include "eigen.h"

//#define ANALYTICPOTENTIAL

static const double r_max = .5;
static const double inner_slope = 3;



double sigma_func(double x) {
	double outer_slope = sigma_index;
	return sigma0 /( pow(x/r_max,-inner_slope) + pow(x/r_max,-outer_slope));
}

double dlogsigma_func(double x) {
	double outer_slope = sigma_index;
	double denom = pow(x/r_max,inner_slope) + pow(x/r_max,outer_slope);
	return outer_slope + (inner_slope - outer_slope)*pow(x/r_max,outer_slope)/denom;
}

double d2logsigma_func(double x) {
	double outer_slope = sigma_index;
	double denom = pow(x/r_max,inner_slope) + pow(x/r_max,outer_slope);
	denom *= denom;
	return -(outer_slope-inner_slope)*(outer_slope-inner_slope)*pow(x/r_max,inner_slope+outer_slope)/denom;
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
