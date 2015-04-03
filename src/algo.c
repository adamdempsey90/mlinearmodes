#include "eigen.h"


int calc_matrices(double complex *mat, double complex *bcmat) {
	int i,j,indx;
	
	double complex A, B, C, G;
/* Compute the matrix including all of its component matrices */

	
// #ifdef OPENMP
// #pragma omp parallel 
// #pragma omp for private(i,A,B,C,indx,j,
// #endif	

	for(i=0;i<N;i++ ) {
		calc_coefficients(i,&A,&B,&C,&G);

		for(j=0 ; j<N ; j++) {
			indx = j + N*i;
		
			
			KD[indx] = G*sigma[j]*D[indx];
//			HL[indx] = G*D[indx];
			HL[indx] = 0;
			mat[indx] = 0;
			bcmat[indx] = 0;
			
			
			if (i==j) {
				mat[indx] += A;
				bcmat[indx] += 1;
//				KD[indx] = -sigma[j]* ( dlds[j] + D[indx]);
//				HL[indx] = G*(2. + D[indx]);
				KD[indx] += G*sigma[j]*dlds[j];
//				HL[indx] += 2*G;
			}
// 			else {
// 				mat[indx] = 0;
// 				bcmat[indx] = 0;
// 				KD[indx] = -D[indx]*sigma[j];
// 				HL[indx] = G*D[indx];
// 			}
			
			mat[indx] += B*D[indx] + C*D2[indx];
			
			H[indx] = 0;
			DKD[indx] = 0;
			
		
		}
	}
			
			
			
//	matmat(kernel,D,KD,1,0,N);
//	matmat(D,KD,DKD,1,0,N);
	
//	cmatmat(HL,H,DKD,1,0,N);
	

#ifdef SELFGRAVITY	
	cmatmat(kernel,KD,DKD,1,0,N);
	for(i=0;i<N*N;i++) mat[i] += DKD[i];
#endif


/* Set Boundary Conditions */

	lagrangian_pressure_bc_inner(mat, bcmat);
	lagrangian_pressure_bc_outer(mat, bcmat);

	return 0;
}


void calc_coefficients(int i, double complex *A, double complex *B, double complex *C, double complex *G) {



#ifdef BAROTROPIC
	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) * ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega[i]-kappa[i];
	*B = (*C) * (2 + dldc2[i] + dlds[i]) ;
	


#endif


#ifdef ISOTHERMAL

	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) *( 2 * dlds[i] + d2lds[i]) + omega[i]-kappa[i];
	*B = (*C) *( 2 + dlds[i]);
	


#endif




#ifdef ADIABATIC

	*C = temp[i]/(2*omega[i]*r[i]*r[i]);

	*B = (*C) * adi_gam*(2 + dldpres[i]);
	
	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];

	*C *= adi_gam;

#endif


#ifdef COOLING

//	double complex cool_fac = ( 1 + 1j*beta_cool)/(1 + beta_cool*beta_cool);
	double complex cool_fac = ( 1 + 1j*beta_cool*(adi_gam-1))/(1 + beta_cool*beta_cool*(adi_gam-1)*(adi_gam-1));
	
// 	*C = cool_fac *  temp[i]/(2*omega[i]*r[i]*r[i]);
// 
// 	*B = (*C) * ( adi_gam*(2 + dldpres[i]) + 1j*beta_cool*dldpres[i]);
// 	
// 	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];
// 
// 	*C *= adi_gam;

	*C = temp[i]/(2 * omega[i] * r[i] * r[i]);
	
	*A = omega[i] - kappa[i] + (*C) * ( 2 * dlds[i] + d2lds[i] 
				  + cool_fac * (d2ldtemp[i] + dldtemp[i]*(2 + dlds[i] + dldtemp[i])));
				  
	*B = (*C) * (2 + dlds[i] + cool_fac * ( adi_gam*dldtemp[i] + ( adi_gam - 1)*(dlds[i] + 2)));
	
	*C *= 1 + cool_fac * (adi_gam -1);
	
#endif







// #ifdef BAROTROPIC
// 	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
// 	
// 	*A = ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega[i]-kappa[i];
// 	*B = (2 + dldc2[i] + dlds[i]) ;
// 	
// 	*A *= (*C);
// 	*B *= (*C);
// #else
// 
// #ifdef ISOTHERMAL
// 
// 	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
// 	
// 	*A = 2 * dlds[i] + d2lds[i] + omega[i]-kappa[i];
// 	*B = 2 + dlds[i];
// 	
// 	*A *= (*C);
// 	*B *= (*C);
// 	
// #else // ADIABATIC OR COOLING
// 	double complex cool_fac = ( 1 + 1j*beta_cool)/(1 + beta_cool*beta_cool);
// 	
// 	*C = cool_fac *  temp[i]/(2*omega[i]*r[i]*r[i]);
// 
// 	*B = (*C) * ( adi_gam*(2 + dldpres[i]) + 1j*beta_cool*dldpres[i]);
// 	
// 	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];
// 
// 	*C *= adi_gam;
// 
// #endif
// #endif


#ifdef SELFGRAVITY
	*G = -1.0/(2*omega[i]*r[i]*r[i]*r[i]);
#else
	*G = 0;
#endif


#ifdef VISCOSITY
	double complex norm;
// 	double q = dldom[i];
// 	double qp = d2dom[i];
// 	double beta = dlds[i];
// 	double gam = dldnu[i];
// 	double betap = d2lds[i];
	
// 	norm = I/(r[i]*omega[i]);
// 	
// 	temp = -2.*(12.+9.*beta + 9.*gam + 7.*qp) + q*(-31.+4.*beta*(-2.+3*beta)-14.*gam-12.*betap);
// 
// 	*A += temp*(nu[i]*omega[i]*norm/(12*r[i]));
// 	
// 	temp = 31. + 14.*beta + 2*q*(11. + 6*beta) + 14.*gam;
// 	
// 	*B += -temp*norm*(nu[i]*omega[i]/(12*r[i]));
// 
// 	*C += -(nu[i]*omega[i]*norm/(6*r[i]))*(7 + 6*q);

	norm = temp[i]/(2*omega[i]*r[i]*r[i]);
/* Shear Viscosity */

	*A -= I * (alpha_s / 6.) * norm * (33.5 + dldtemp[i] - 2*dlds[i] + 18*d2lds[i]);
	*B -= I * (alpha_s / 3.) * norm * (-9.5 - 7*dldtemp[i] + 2*dlds[i]);
	*C -= I * (2 *alpha_s / 3.) * norm;

/* Bulk Viscosity */

	*B += I * alpha_b * norm * (2.5 + dldtemp[i] + dlds[i]);
	*C -= I * alpha_b * norm;

#endif

	return; 

}

