#include "eigen.h"

double sigma_func(double x) {
	return sigma0 * pow(x, sigma_index);
}

double dlogsigma_func(double x) {
	return sigma_index;
}

double d2logsigma_func(double x) {
	return 0;
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