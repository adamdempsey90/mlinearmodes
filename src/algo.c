#include "eigen.h"

void calc_coefficients(void);


int calc_matrices(double complex *mat, double complex *bcmat) {
	int i,j,indx,mindx;

/* Compute the matrix including all of its component matrices */



#ifndef NOPRESSURE
	calc_coefficients();
#endif
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			mindx=  j + ncols*i;
			indx = j+ N*i;

			mat[mindx] = ( omega_prec[i] + coeffs_A[i])*Identity[indx]
					    + coeffs_B[i]*D[indx]
					    + coeffs_C[i]*D2[indx];
			bcmat[mindx] = Identity[indx];
		}
	}


#ifdef SELFGRAVITY
	add_sg(mat,bcmat);
#endif




/* Set Boundary Conditions */
	set_bc(mat,bcmat);


#ifdef PLANETS
	add_planets(mat,bcmat);
#endif

	return 0;
}


void calc_coefficients(void) {
	int i;
	double complex norm;
	for(i=0;i<N;i++) {


#ifdef BAROTROPIC
		coeffs_C[i]  = c2[i]/(2*omega[i]*r[i]*r[i]);

		coeffs_A[i]  = coeffs_C[i]  * ( dlds[i]*(2 + dldc2[i]) + d2lds[i]);
		coeffs_B[i]  = coeffs_C[i]  * (2 + dldc2[i] + dlds[i]) ;



#endif


#ifdef ISOTHERMAL

		coeffs_C[i]  = c2[i]/(2*omega[i]*r[i]*r[i]);

		coeffs_A[i]  = coeffs_C[i]  *( 2 * dlds[i] + d2lds[i]) ;
		coeffs_B[i]  = coeffs_C[i]  *( 2 + dlds[i]);



#endif




#ifdef ADIABATIC

		coeffs_C[i]  = temp[i]/(2*omega[i]*r[i]*r[i]);

		coeffs_B[i]  = coeffs_C[i] * adi_gam*(2 + dldpres[i]);

		coeffs_A[i] = coeffs_C[i]  * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]);

		coeffs_C[i]  *= adi_gam;

#endif


#ifdef COOLING

//	double complex cool_fac = ( 1 + 1j*beta_cool)/(1 + beta_cool*beta_cool);
	double complex cool_fac = ( 1 + 1j*beta_cool*(adi_gam-1))/(1 + beta_cool*beta_cool*(adi_gam-1)*(adi_gam-1));


	coeffs_C[i] = temp[i]/(2 * omega[i] * r[i] * r[i]);

	coeffs_A[i]  = coeffs_C[i]  * ( 2 * dlds[i] + d2lds[i]
				  + cool_fac * (d2ldtemp[i] + dldtemp[i]*(2 + dlds[i] + dldtemp[i])));

	coeffs_B[i]  = coeffs_C[i]  * (2 + dlds[i] + cool_fac * ( adi_gam*dldtemp[i] + ( adi_gam - 1)*(dlds[i] + 2)));

	coeffs_C[i]  *= (1 + cool_fac * (adi_gam -1));

#endif






/*	Shear Viscosity	*/
	if (alpha_s !=0) {
		norm = - I * alpha_s * temp[i]/(24 * omega[i]*r[i]*r[i]);
		coeffs_A[i]  += norm*3*( 9 + 4*dlds[i]*(-2+3*dlds[i])+2*dldtemp[i]*(-1+6*dlds[i])-12*d2lds[i] );
		coeffs_B[i] += norm*(-82+64*dldtemp[i]+28*dlds[i]);
		coeffs_C[i] -= norm*8;
	}
/* Bulk Viscosity */
	if (alpha_b != 0) {

		norm = I*alpha_b * temp[i]/(2 * omega[i]*r[i]*r[i]);

		coeffs_B[i] += norm*(-1 + dldtemp[i] + dlds[i]);
		coeffs_C[i] += norm;
	}



	}

	return;

}
