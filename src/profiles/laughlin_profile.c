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
/* Polytrope with gamma = flare_index
 * Treat this as an isothermal disk with
 * T = K \sigma^(\gamma -1), so that
 * P = K \sigma^\gamma
*/
//1.76398*sigma0*flare_index*
	return h0*h0*pow(x,temp_index);
//	return 0.74528*sigma0*flare_index*pow(sigma_func(x),flare_index-1);
}

double dlogtemp_func(double x) {
	return temp_index;
//	return (flare_index - 1)*dlogsigma_func(x);
}

double d2logtemp_func(double x) {
	return 0;
//	return  (flare_index - 1)*d2logsigma_func(x);
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
	return sqrt(temp_func(x))/omk_func(x);
}


// int analytic_potential(void) {
// 	return 0;
// }
//
// double omega_prec_grav_analytic(double x) {
// 	return 0;
// }
