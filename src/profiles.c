#include "eigen.h"



double sigma_function(double rval) {

	double ans = 0;

#ifdef EXPDECAY

	ans = pow(rval,sigma_index) * exp(-rval/RDECAY);

#else
#ifdef EXPDECAY2

	ans = pow(rval,sigma_index) * exp(-pow(rval/RDECAY,2));

#else
	
	ans = pow(rval,sigma_index);

#endif
#endif
	return ans;

}

double dlogsigma_dlogr(double rval) {
	double ans = 0;
#ifdef EXPDECAY

	ans = sigma_index - rval/RDECAY;
#else
#ifdef EXPDECAY2

	ans = sigma_index - 2*pow(rval/RDECAY,2);
#else

	ans = sigma_index;	

#endif
#endif
	
	return ans ;
	
}

double dlogsigma2_dlogr2(double rval) {
	double ans = 0;
#ifdef EXPDECAY
	ans = -rval/RDECAY;
#else
#ifdef EXPDECAY2

	ans = -4*pow(rval/RDECAY,2);
#else

	ans = 0;	
#endif
#endif
	return ans;
	
}


double integrand(double x, void *params) {
	double r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee; 
	double ans;
	int i=  *(int *)params;
	
	r1 = r[i];
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r2*r2;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp3 = rp2*rp1;
	rp4 = rp2*rp2;
	eps1 = eps * scaleH_func(r1);
	eps2 = eps1*eps1;
	eps4 = eps2*eps2;
	eps6 = eps4 * eps2;
	r_p_rp = r1 + rp1;
	r_p_rp *= r_p_rp;
	r_m_rp = r1 - rp1;
	r_m_rp *= r_m_rp;
	
	kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
		
	ans=2*(eps2 - 2 *eps1 *r1 - r2 + rp2)*(eps2 + 2 *eps1*r1 - r2 + rp2)*(eps2 + r2 + rp2)*ee; 
	ans -= 2*(eps2 + (r1 - rp1)*(r1-rp1))*(eps4 + r4 + 2*(eps1 - r1)*(eps1 + r1)*rp2 + rp4)*ek;
	ans /= (pow(eps2 + (r1 - rp1)*(r1-rp1),2) *pow(eps2 + (r1 + rp1)*(r1+rp1),1.5));
	
	ans *= -rp2*sigma_func(rp1);
	
	return ans;

}

double sigma_func(double x) {

	return sigma0 * pow(x, sigma_index);

}

double temp_func(double x) {
	return h0*h0*pow(x,2*flare_index-1);
}

double omk_func(double x) {
	return pow(x,-1.5);
}
double scaleH_func(double x) {

	return h0*x*pow(x,flare_index);

}