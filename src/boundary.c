#include "eigen.h"



void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at inner boundary */

	int j,indx;
#ifdef COOLING
// 	double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
// 	eta_fac /= (1 - eta_fac);
		double complex cool_fac = (adi_gam - 1)/(adi_gam*(1-I*beta_cool*(adi_gam-1)));
		double complex fac1, fac2;
#endif
	for(j=0;j<N;j++) {
		indx = j;
		mat[indx] = D[indx];
		bcmat[indx] = 0;

		if (j==0) {
#ifdef COOLING
//			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
				fac1 = dldtemp[j] - cool_fac*adi_gam*dldtemp[j];
				fac2 = 1 + cool_fac*adi_gam;
#ifdef VISCHEATING
				fac1 -= cool_fac*.75*alpha_s*I;
				fac2 += cool_fac*1.5*alpha_s*I;
#endif
			mat[indx] -= fac1/fac2;

#endif
#ifdef ISOTHERMAL
			mat[indx] -= dldtemp[j];
#endif

		}
	}
	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at outer boundary */
#ifdef COOLING
	// double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
	// eta_fac /= (1 - eta_fac);
	double complex cool_fac = (adi_gam - 1)/(adi_gam*(1-I*beta_cool*(adi_gam-1)));
	double complex fac1, fac2;
#endif
	int j,indx;
	for(j=0;j<N;j++) {
		indx= j + nrows*(N-1);

		mat[indx] = D[j+N*(N-1)];
		bcmat[indx] = 0;

		if (j==(N-1)) {
#ifdef COOLING
//			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
			fac1 = dldtemp[j] - cool_fac*adi_gam*dldtemp[j];
			fac2 = 1 + cool_fac*adi_gam;
#ifdef VISCHEATING
			fac1 -= cool_fac*.75*alpha_s*I;
			fac2 += cool_fac*1.5*alpha_s*I;
#endif
			mat[indx] -= fac1/fac2;
#endif
#ifdef ISOTHERMAL
			mat[indx] -= dldtemp[j];
#endif

		}

	}
	return;
}


void zero_e_bc_inner(double complex *mat, double complex *bcmat) {

	int j,indx;
	for(j=0;j<N;j++) {
		indx = j;
		mat[indx] = 0;
		bcmat[indx] = 0;
		if (j==0) {
			mat[indx] = 1;
		}

	}

	return;
}


void zero_e_bc_outer(double complex *mat, double complex *bcmat) {

	int j,indx;
	for(j=0;j<N;j++) {
		indx = j + nrows*(N-1);
		mat[indx] = 0;
		bcmat[indx] = 0;
		if (j==N-1) {
			mat[indx] = 1;
		}
	}

	return;
}
void user_gradient_bc_outer(double complex *mat, double complex *bcmat,double complex val) {

/* Set up zero Lagrangian Pressure B.C at outer boundary */
	int j,indx;
	for(j=0;j<N;j++) {
		indx= j + nrows*(N-1);

		mat[indx] = D[j+N*(N-1)];
		bcmat[indx] = 0;

		if (j==(N-1)) {
			mat[indx] -= val;
		}

	}
	return;
}
void set_bc(double complex *mat, double complex *bcmat) {

#ifndef NOPRESSURE

#ifdef ZEROBCIN
	zero_e_bc_inner(mat,bcmat);
#endif

#ifdef ZEROBCOUT
	zero_e_bc_outer(mat,bcmat);
#endif

#ifdef PRESBCIN
	lagrangian_pressure_bc_inner(mat, bcmat);
#endif

#ifdef PRESBCOUT
	lagrangian_pressure_bc_outer(mat, bcmat);
#endif

#ifdef GRADBCIN
	user_gradient_bc_outer(mat,bcmat,.5);
#endif


#endif
	return;
}
