#include "eigen.h"


double sig_param = 1; 
double frozen_param = exp(-5.);

double sigma_potential_func(double x);
double dlogsigma_potential_func(double x);

double sigma_func(double x) {
	double live_fac =  Mdisk*sig_param*pow(x*x+sig_param*sig_param,-1.5) / (2*M_PI);
	double frozen_fac = 0; //Mdisk*frozen_param*pow(x*x+frozen_param*frozen_param,-1.5)/(2*M_PI);
	return  fmax(live_fac - frozen_fac,0) +1e-8;
}

double dlogsigma_func(double x) {
	double live_fac =  -3*x*x/(x*x + sig_param*sig_param);
	double frozen_fac =  0; //-3*x*x/(x*x + frozen_param*frozen_param);
	return live_fac - frozen_fac;
}

double d2logsigma_func(double x) {
	double live_fac = -6*sig_param*sig_param*x*x * pow(sig_param*sig_param+x*x,-2);
	double frozen_fac = 0; //-6*frozen_param*frozen_param*x*x * pow(frozen_param*frozen_param+x*x,-2);
	return live_fac - frozen_fac;
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

