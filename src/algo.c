#include "eigen.h"

void calc_coefficients(void);


int calc_matrices(double complex *mat, double complex *bcmat) {
	int i,j,indx;
	
/* Compute the matrix including all of its component matrices */



#ifndef NOPRESSURE
	calc_coefficients();
#endif	
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			indx=  j +N*i;
			
			mat[indx] = ( omega_prec[i] + coeffs_A[i])*Identity[indx] 
					    + coeffs_B[i]*D[indx] 
					    + coeffs_C[i]*D2[indx];
			bcmat[indx] = Identity[indx];
		}
	}	
			

#ifdef SELFGRAVITY		
	add_sg(mat,bcmat);
#endif	


	

/* Set Boundary Conditions */
	set_bc(mat,bcmat);
	

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




	

// 	norm = temp[i]/(2*omega[i]*r[i]*r[i]);
// /* Shear Viscosity */
// 
// 	coeffs_A[i]  -= I * (alpha_s / 6.) * norm * (33.5 + dldtemp[i] - 2*dlds[i] + 18*d2lds[i]);
// 	coeffs_B[i]  -= I * (alpha_s / 3.) * norm * (-9.5 - 7*dldtemp[i] + 2*dlds[i]);
// 	coeffs_C[i]  -= I * (2 *alpha_s / 3.) * norm;
	if (alpha_s !=0) {
		norm  = -I * alpha_s * temp[i]/(24 * omega[i] *r[i]*r[i]);
	
		coeffs_A[i] += 3*norm*(-3+2*dldtemp[i] - 4*dlds[i] + 12*d2lds[i]);
		coeffs_B[i] -= 2*norm *(25+18*dldtemp[i]);
		coeffs_C[i] += 8*norm; 
	}
/* Bulk Viscosity */
	if (alpha_b != 0) {
	
		norm = I*alpha_b * temp[i]/(2 * omega[i]*r[i]*r[i]);
	
		coeffs_B[i] += norm*(2 + dldtemp[i] + dlds[i]);
		coeffs_C[i] += norm;
	}
// 	coeffs_B[i]  += I * alpha_b * norm * (2.5 + dldtemp[i] + dlds[i]);
// 	coeffs_C[i] -= I * alpha_b * norm;



	}

	return; 

}

