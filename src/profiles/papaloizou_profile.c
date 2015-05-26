#include "eigen.h"

//#define ANALYTICPOTENTIAL

static const double rin = 1;
static const double rout = 100;
static const double polyn = 1.5;

double sigma_func(double x) {
	return sigma0 *  pow(temp_func(x),polyn);
}

double dlogsigma_func(double x) {
	return polyn * dlogtemp_func(x);
}

double d2logsigma_func(double x) {
	return polyn * d2logtemp_func(x);
}


double temp_func(double x) {
	double fac_in = (1 - pow(rin/x,10));
	double fac_out = (1 - pow(x/rout,10));
	return scaleH_func(x)*scaleH_func(x) * fac_in * fac_out * pow(rin/x,3);
}

double dlogtemp_func(double x) {
	double fac_in = (1 - pow(x/rin,10));
	double fac_out = (1- pow(rout/x,10));
	return temp_index - 10/fac_in + 10/fac_out;
}

double d2logtemp_func(double x) {
	double fac_in = pow(rout/x,10);
	double fac_out = pow(rin/x,10);

	return -100*fac_in *pow(1-fac_in,-2) - 100*fac_out*pow(1-fac_out,-2);
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
