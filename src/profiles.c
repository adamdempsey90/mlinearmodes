#include "eigen.h"


void set_papaloizou_profile(void) {
	int i;
	
	for(i=0;i<N;i++) {
		c2[i] = scaleH[i] * scaleH[i] / (r[i]*r[i]*r[i])*(1 - pow(r[i],-10))*(1-pow(r[i]/rout,10));
		sigma[i] = pow(c2[i],1.5);
		temp[i] = c2[i];
	}

	return;

}


void set_heemskerk_profile(void) {
	int i;
	
	for(i=0;i<N;i++) {
		sigma[i] = pow(r[i]-ri + .05,2)*pow(ro - r[i] +.05,2.5);

#ifndef INPUTMASS
		sigma[i] *= sigma0;
#endif

		c2[i] = sigma[i]*h0*h0;
		temp[i] = c2[i];
	}

	return;

}





void set_mlin_profile(void) {
	int i;
	
	for(i=0;i<N;i++) {
		
		c2[i] = h0 * pow(r[i],-.5*flare_index);
	
		scaleH[i] = c2[i]/omega[i];
		c2[i] *= c2[i];

	
//		sigfac = .05 * pow(2.0,2 - .5*flare_index - 1.5)/(2 * M_PI * bump_function(2.0));
//		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-2.0);
		sigfac = h0 /(2 * M_PI * bump_function(2.0));
		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-(1.5 + .5*flare_index));
	
		
		temp[i] = c2[i];
	}

	return;

}





double bump_function(double rval) {
	
	
	double delta1 = 5.0 * .05 * pow(1.0,1.5-.5*flare_index);
	double delta2 = 5.0 * .05 * pow(2.0,1.5-.5*flare_index);

	double fac1 = .5*(1-.1)*(1+tanh((rval-1.0)/delta1))+.1;
	double fac2 = .5*(1-.1)*(1-tanh((rval-2.0)/delta2))+.1;
	
	return (fac1*fac2);

}

	
	
#ifdef TESTFUNCTION

double complex test_function(double rval) {
	double complex result;
	
	result = 1 / pow(rval,2);
		
	return result;
}

void fill_mat(double complex *mat, double complex *bcmat) {
	int i,j,indx;

#ifdef OPENMP
#pragma omp parallel private(indx,i,j) shared(D2,r,mat)
#pragma omp for schedule(static)
#endif	
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		mat[indx] = D2[indx] * test_function(r[i]);
		if (i==j) {
			bcmat[indx] = 1;
		}
		else {
			bcmat[indx] = 0;
		}
	}
	
	for(j=0;j<N;j++) {
		indx= j;
	
		mat[indx] = D[indx];
		bcmat[indx] = 0;
		indx +=  N*(N-1);
		mat[indx] = D[indx];
		bcmat[indx] = 0;
	}
	
	return;
}

#endif




double sigma_function(double rval) {

	ans = 0;

#ifdef EXPDECAY
	ans = pow(rval,sigma_index) * exp(-rval/rdecay);
#else

#ifdef EXPDECAY2
	ans = pow(rval,sigma_index) * exp(-pow(rval/rdecay,2));
#else
	ans = pow(rval,sigma_index);

#endif

#endif


	return ans;
}

double dlogsigma_dlogr(double rval) {
	
	ans = 0;
#ifdef EXPDECAY
	ans = sigma_index - rval/r_decay;
#else
#ifdef EXPDECAY2
	ans = sigma_index - 2*pow(rval/rdecay,2);
#else
	ans = sigma_index;	
#endif
#endif


	return ans;

}

double dlogsigma2_dlogr2(double rval) {

	ans = 0;
#ifdef EXPDECAY
	ans = -rval/r_decay;
#else
#ifdef EXPDECAY2
	ans = -4*pow(rval/rdecay,2);
#else
	ans =  0;	
#endif
#endif


	return ans;

}



