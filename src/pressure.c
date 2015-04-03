#include "eigen.h"


void pressure_omega_correction(void) {
	int i;
	
	double *omegap2 = (double *)malloc(sizeof(double)*N);
	double *kappap2 = (double *)malloc(sizeof(double)*N); 
	
	for(i=0;i<N;i++) {
	
#ifdef BAROTROPIC

		omegap2[i] = c2[i]*dlds[i]/(r[i]*r[i]);
		
		
#ifdef EXACTKAPPA
		kappap2[i] = c2[i]*((2 + dldc2[i]) * dlds[i] + d2lds[i])/(r[i]*r[i]);
#else
		kappap2[i] = 4*omegap2[i];
#endif


#else

		omegap2[i] = dldpres[i] * temp[i]/(r[i]*r[i]);
		
		
#ifdef EXACTKAPPA
		kappap2[i] = temp[i]*((2 + dldtemp[i]) * dldpres[i] + d2ldpres[i])/(r[i]*r[i]);
#else
		kappap2[i] = 4*omegap2[i];
#endif

#endif // BAROTROPIC
	}
	
#ifndef EXACTKAPPA
	matvec(D,omegap2,kappap2,1,1,N);
#endif

	for(i=0;i<N;i++) {
		omega2[i] += omegap2[i];
		kappa2[i] += kappap2[i];
	}


	free(omegap2);
	free(kappap2);
	return;

}