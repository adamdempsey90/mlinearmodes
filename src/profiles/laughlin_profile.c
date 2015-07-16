#include "eigen.h"

//#define ANALYTICPOTENTIAL

const double dens_width2 = .05;
const double dens_peak_rad = .45;


double sigma_func(double x) {
	return sigma0 * exp(-(x-dens_peak_rad)*(x-dens_peak_rad)/(dens_width2));
}

double dlogsigma_func(double x) {
	return 2*(dens_peak_rad - x)*x/dens_width2;
}

double d2logsigma_func(double x) {
	return 2*(dens_peak_rad - 2*x)*x/dens_width2;
}


double temp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = flare_index*pow(sigma_func(x),flare_index-1);
//	printf("%.2e\t%.2e\n",sigma_func(x),res);
#else
	res = h0*h0*pow(x,temp_index);
#endif
	return res;
}

double dlogtemp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = (flare_index - 1)*dlogsigma_func(x);
#else
	res = temp_index;
#endif
	return res;
}

double d2logtemp_func(double x) {
	double res;
#ifdef POLYTROPE
	res = (flare_index - 1)*d2logsigma_func(x);
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
	res = sqrt(flare_index)*pow(sigma_func(x),.5*(flare_index-1))*pow(omk_func(x),-1);
#else
	res =  h0*x*pow(x,flare_index);
#endif
	return res;
}


// int analytic_potential(void) {
// 	return 0;
// }
//
// double omega_prec_grav_analytic(double x) {
// 	return 0;
// }
