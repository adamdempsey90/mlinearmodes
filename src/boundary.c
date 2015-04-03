#include "eigen.h"


void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at inner boundary */

	int j,indx;
#ifdef COOLING
	double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
	eta_fac /= (1 - eta_fac);
#endif
	for(j=0;j<N;j++) {
		indx = j;
		mat[indx] = D[indx];
		bcmat[indx] = 0;		
		
		if (j==0) {
#ifdef COOLING
//			mat[indx] += I*beta_cool * dldpres[j] / adi_gam;
//			mat[indx] -= dldtemp[j]*(beta_cool*beta_cool - I*adi_gam*beta_cool)/(beta_cool*beta_cool + adi_gam*adi_gam);
			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
#endif	
#ifdef ISOTHERMAL
			mat[indx] -= dldc2[j];
#endif

		}	
	}
	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at outer boundary */	
#ifdef COOLING
	double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
	eta_fac /= (1 - eta_fac);
#endif
	int j,indx;
	for(j=0;j<N;j++) {
		indx= j + N*(N-1);
	
		mat[indx] = D[indx];
		bcmat[indx] = 0;
		
		if (j==(N-1)) {
#ifdef COOLING
//			mat[indx] += I*beta_cool * dldpres[j] / adi_gam;
//			mat[indx] -= dldtemp[j] * (beta_cool*beta_cool - I*adi_gam*beta_cool)/(beta_cool*beta_cool + adi_gam*adi_gam);
			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
#endif	
#ifdef ISOTHERMAL
			mat[indx] -= dldc2[j];
#endif

		}
		
	}
	return;
}