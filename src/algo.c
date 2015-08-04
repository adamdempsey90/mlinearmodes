#include "eigen.h"

void calc_coefficients(void);
void calc_A_coeff(int i, double complex *Amat);
void calc_B_coeff(int i, double complex *Bmat);
void calc_C_coeff(int i, double complex *Cmat);

int calc_matrices(double complex *mat, double complex *bcmat) {
	int i,j,indx,mindx;
	int k,l,cindx,s,dindx;

/* Compute the matrix including all of its component matrices */

// OPENMP
	for(i=0;i<N*NF*NF;i++) {
		coeffs_A[i] = 0;
		coeffs_B[i] = 0;
		coeffs_C[i] = 0;
	}

	printf("Calculating coefficients...\n");
	calc_coefficients();

	printf("Filling matrix...\n\tPressure...\n\tViscosity...\n");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			indx = getindex2(i,j,N);
			for(k=0;k<NF;k++) {
				for(l=0;l<NF;l++) {
					mindx = getindex4(i,j,k,l,NF,N);
					bcmat[mindx] = Identity[mindx];
					if (i==j) {
						mat[mindx] = coeffs_A[getindex3(i,k,l,NF)];
					}
					else {
						mat[mindx] = 0;
					}
					for(s=0;s<NF;s++) {
						cindx = getindex3(i,k,s,NF);
						dindx = getindex4(i,j,s,l,NF,N);
						mat[mindx] += coeffs_B[cindx]*D[dindx] + coeffs_C[cindx]*D2[dindx];
					}

//					mat[mindx] = coeffs_A[cindx]*Identity[mindx] + coeffs_B[cindx]*D[mindx]
//											+ coeffs_C[cindx]*D2[mindx];


				}
			}
		}
	}


#ifdef SELFGRAVITY
	printf("\tSelfgravity...\n");
	add_sg(mat,bcmat);
#endif




/* Set Boundary Conditions */
	printf("Setting Boundary Conditions...\n");
	set_bc(mat,bcmat);


	return 0;
}


void calc_coefficients(void) {
	int i,k,l,indx;
	double rval,rval2;

// OPENMP
	for(i=0;i<N;i++) {
		calc_A_coeff(i,&coeffs_A[NF*NF*i]);
		calc_B_coeff(i,&coeffs_B[NF*NF*i]);
		calc_C_coeff(i,&coeffs_C[NF*NF*i]);
	}
	printf("\t\tConverting to lnr derivatives...\n");
	for(i=0;i<N;i++) {
		rval = r[i];
		rval2 = rval*rval;
		for(k=0;k<NF;k++) {
			for(l=0;l<NF;l++) {
				indx = getindex3(i,k,l,NF);
				coeffs_B[indx] = coeffs_B[indx]/rval - coeffs_C[indx]/rval2;
				coeffs_C[indx] = coeffs_C[indx]/rval2;
				coeffs_A[indx] /= (I*mval);
				coeffs_B[indx] /= (I*mval);
				coeffs_C[indx] /= (I*mval);
				// coeffs_A[indx] /= (I*mval);
				// coeffs_C[indx] /= (I*mval*rval2);
				// coeffs_B[indx] /= (I*mval*rval);
				// coeffs_B[indx] -= coeffs_C[indx];
			}
		}
	}

	return;
}


void calc_A_coeff(int i, double complex *Amat) {
/* Calculate the coefficients for the A matrix assuming we're given the correct
 offset in memory to start from ( NF^2 *i )
*/
double mu, delta, om, mup, deltap,kap2,rval,sigval,kapom,rval2,Tval,Pval,vr,dlvr;
mu = dlds[i];
mup = d2lds[i];
delta = dldtemp[i];
deltap = d2ldtemp[i];
om = omega[i];
rval = r[i];
rval2 = rval*rval;
sigval = sigma[i];
Tval  = temp[i];
kap2 = kappa2[i];
kapom = kap2/(2*om);
Pval = pres[i];
vr = vrbar[i];
dlvr = dlvrbar[i];


double nus = alpha_s * Tval/om;
double dldnu = delta + 1.5;


#ifdef BAROTROPIC
	Amat[getindex2(0,0,NF)] = I*mval*om;
	Amat[getindex2(0,1,NF)] = - 2 *om;
	Amat[getindex2(0,2,NF)] = Tval*(delta-mu)/(rval*sigval);

	Amat[getindex2(1,0,NF)] = kapom;
	Amat[getindex2(1,1,NF)] = I*mval*om;
	Amat[getindex2(1,2,NF)] = I*mval*Tval/(rval*sigval);

	Amat[getindex2(2,0,NF)] = ( 1 + mu)*sigval/rval;
	Amat[getindex2(2,1,NF)] = I*mval*sigval/rval;
	Amat[getindex2(2,2,NF)] = I*mval*om;

#endif

#ifdef ISOTHERMAL
	Amat[getindex2(0,0,NF)] = I*mval*om;
	Amat[getindex2(0,1,NF)] = - 2 *om;
	Amat[getindex2(0,2,NF)] = - Tval*mu/(rval*sigval);

	Amat[getindex2(1,0,NF)] = kapom;
	Amat[getindex2(1,1,NF)] = I*mval*om;
	Amat[getindex2(1,2,NF)] = I*mval*Tval/(rval*sigval);

	Amat[getindex2(2,0,NF)] = ( 1 + mu)*sigval/rval;
	Amat[getindex2(2,1,NF)] = I*mval*sigval/rval;
	Amat[getindex2(2,2,NF)] = I*mval*om;

#endif

#ifdef COOLING
	Amat[getindex2(0,0,NF)] = I*mval*om;
	Amat[getindex2(0,1,NF)] = - 2 *om;
	Amat[getindex2(0,2,NF)] =  - Tval*(mu+delta)/(rval*sigval)
	Amat[getindex2(0,3,NF)] = 0;

	Amat[getindex2(1,0,NF)] = kapom;
	Amat[getindex2(1,1,NF)] = I*mval*om;
	Amat[getindex2(1,2,NF)] = 0;
	Amat[getindex2(1,3,NF)] = I*mval/(rval*sigval);

	Amat[getindex2(2,0,NF)] = ( 1 + mu)*sigval/rval;
	Amat[getindex2(2,1,NF)] = I*mval*sigval/rval;
	Amat[getindex2(2,2,NF)] = I*mval*om;
	Amat[getindex2(2,3,NF)] = 0;

	Amat[getindex2(3,0,NF)] = ( adi_gam + mu +delta)*Pval/rval;
	Amat[getindex2(3,1,NF)] = I*mval*adi_gam*Pval/rval;
	Amat[getindex2(3,2,NF)] = 0;
	Amat[getindex2(3,3,NF)] = I*mval*om;
#endif

// Viscosity
	Amat[getindex2(0,0,NF)] += nus*(2+mval2)/rval2;
	Amat[getindex2(0,1,NF)] += nus*3*I*mval/rval2;
	Amat[getindex2(0,2,NF)] += nus*-I*mval*(kapom-2*om)/rval2;

	Amat[getindex2(1,0,NF)] += nus*-I*mval*(3+dldnu+mu)/rval2;
	Amat[getindex2(1,1,NF)] += nus*(1+2*mval2+dldnu + mu)/rval2;
	Amat[getindex2(1,2,NF)] += nus*mu/(sigval*rval)*(kapom-2*om);

// Radial Velocity

	Amat[getindex2(0,0,NF)] += dlvr/rval;
	Amat[getindex2(0,2,NF)] += nus*mu*dlvr/(rval2*sigval);
	Amat[getindex2(1,1,NF)] += vr/rval;
	Amat[getindex2(1,2,NF)] += -2*I*mval*nus*vr/(rval2*sigval);
	Amat[getindex2(2,2,NF)] += (vr + dlvr)/rval;

	return;

}

void calc_B_coeff(int i, double complex *Bmat) {
/* Calculate the coefficients for the A matrix assuming we're given the correct
 offset in memory to start from ( NF^2 *i )
*/
	double mu, delta, om,mup, deltap,kap2,rval,sigval,kapom,rval2,Tval,Pval,vr,dlvr;
	mu = dlds[i];
	mup = d2lds[i];
	delta = dldtemp[i];
	deltap = d2ldtemp[i];
	om = omega[i];
	rval = r[i];
	rval2 = rval*rval;
	sigval = sigma[i];
	Tval  = temp[i];
	kap2 = kappa2[i];
	kapom = kap2/(2*om);
	Pval = temp[i]*sigma[i];
	vr = vrbar[i];
	dlvr = dlvrbar[i];

	double nus = alpha_s * Tval/om;
	double dldnu = delta + 1.5;

#ifdef BAROTROPIC
	Bmat[getindex2(0,0,NF)] = 0;
	Bmat[getindex2(0,1,NF)] = 0;
	Bmat[getindex2(0,2,NF)] = Tval/sigval;

	Bmat[getindex2(1,0,NF)] = 0;
	Bmat[getindex2(1,1,NF)] = 0;
	Bmat[getindex2(1,2,NF)] = 0;

	Bmat[getindex2(2,0,NF)] = sigval;
	Bmat[getindex2(2,1,NF)] = 0;
	Bmat[getindex2(2,2,NF)] = 0;

#endif

#ifdef ISOTHERMAL
	Bmat[getindex2(0,0,NF)] = 0;
	Bmat[getindex2(0,1,NF)] = 0;
	Bmat[getindex2(0,2,NF)] = Tval/sigval;;

	Bmat[getindex2(1,0,NF)] = 0;
	Bmat[getindex2(1,1,NF)] = 0;
	Bmat[getindex2(1,2,NF)] = 0;

	Bmat[getindex2(2,0,NF)] = sigval;
	Bmat[getindex2(2,1,NF)] = 0;
	Bmat[getindex2(2,2,NF)] = 0;

#endif

#ifdef COOLING
	Bmat[getindex2(0,0,NF)] = 0;
	Bmat[getindex2(0,1,NF)] = 0;
	Bmat[getindex2(0,2,NF)] = 0;
	Bmat[getindex2(0,3,NF)] = 1./sigval;

	Bmat[getindex2(1,0,NF)] = 0;
	Bmat[getindex2(1,1,NF)] = 0;
	Bmat[getindex2(1,2,NF)] = 0;
	Bmat[getindex2(1,3,NF)] = 0;

	Bmat[getindex2(2,0,NF)] = sigval;
	Bmat[getindex2(2,1,NF)] = 0;
	Bmat[getindex2(2,2,NF)] = 0;
	Bmat[getindex2(2,3,NF)] = 0;

	Bmat[getindex2(3,0,NF)] = adi_gam*Pval;
	Bmat[getindex2(3,1,NF)] = 0;
	Bmat[getindex2(3,2,NF)] = 0;
	Bmat[getindex2(3,3,NF)] = 0;

#endif
//Viscosity
	Bmat[getindex2(0,0,NF)] += nus*-2*(1+dldnu+mu)/rval;
	Bmat[getindex2(0,1,NF)] += nus*-I*mval/rval;
	Bmat[getindex2(0,2,NF)] += 0;

	Bmat[getindex2(1,0,NF)] += nus*-I*mval/rval;
	Bmat[getindex2(1,1,NF)] += -nus*(1 + dldnu+mu)/rval;
	Bmat[getindex2(1,2,NF)] += -nus*(kapom-2*om)/sigval;

// Radial Velocity

	Bmat[getindex2(0,0,NF)] += vr;
	Bmat[getindex2(0,2,NF)] += -2*nus*dlvr/(rval*sigval);
	Bmat[getindex2(1,1,NF)] += vr;
	Bmat[getindex2(2,2,NF)] += vr;
	return;
}

void calc_C_coeff(int i, double complex *Cmat) {
/* Calculate the coefficients for the A matrix assuming we're given the correct
 offset in memory to start from ( NF^2 *i )
*/
	double mu, delta, om,mup, deltap,kap2,rval,sigval,kapom,rval2,Tval,Pval;
	mu = dlds[i];
	mup = d2lds[i];
	delta = dldtemp[i];
	deltap = d2ldtemp[i];
	om = omega[i];
	rval = r[i];
	rval2 = rval*rval;
	sigval = sigma[i];
	Tval  = temp[i];
	kap2 = kappa2[i];
	kapom = kap2/(2*om);
	Pval = temp[i]*sigma[i];

	double nus = alpha_s * Tval/om;
	double dldnu = delta + 1.5;

#ifdef BAROTROPIC
	Cmat[getindex2(0,0,NF)] = 0;
	Cmat[getindex2(0,1,NF)] = 0;
	Cmat[getindex2(0,2,NF)] = 0;

	Cmat[getindex2(1,0,NF)] = 0;
	Cmat[getindex2(1,1,NF)] = 0;
	Cmat[getindex2(1,2,NF)] = 0;

	Cmat[getindex2(2,0,NF)] = 0;
	Cmat[getindex2(2,1,NF)] = 0;
	Cmat[getindex2(2,2,NF)] = 0;

#endif

#ifdef ISOTHERMAL
	Cmat[getindex2(0,0,NF)] = 0;
	Cmat[getindex2(0,1,NF)] = 0;
	Cmat[getindex2(0,2,NF)] = 0;

	Cmat[getindex2(1,0,NF)] = 0;
	Cmat[getindex2(1,1,NF)] = 0;
	Cmat[getindex2(1,2,NF)] = 0;

	Cmat[getindex2(2,0,NF)] = 0;
	Cmat[getindex2(2,1,NF)] = 0;
	Cmat[getindex2(2,2,NF)] = 0;

#endif

#ifdef COOLING
	Cmat[getindex2(0,0,NF)] = 0;
	Cmat[getindex2(0,1,NF)] = 0;
	Cmat[getindex2(0,2,NF)] = 0;
	Cmat[getindex2(0,3,NF)] = 0;

	Cmat[getindex2(1,0,NF)] = 0;
	Cmat[getindex2(1,1,NF)] = 0;
	Cmat[getindex2(1,2,NF)] = 0;
	Cmat[getindex2(1,3,NF)] = 0;

	Cmat[getindex2(2,0,NF)] = 0;
	Cmat[getindex2(2,1,NF)] = 0;
	Cmat[getindex2(2,2,NF)] = 0;
	Cmat[getindex2(2,3,NF)] = 0;

	Cmat[getindex2(3,0,NF)] = 0;
	Cmat[getindex2(3,1,NF)] = 0;
	Cmat[getindex2(3,2,NF)] = 0;
	Cmat[getindex2(3,3,NF)] = 0;

#endif

// Viscosity
	Cmat[getindex2(0,0,NF)] += -2*nus;
	Cmat[getindex2(0,1,NF)] += 0;
	Cmat[getindex2(0,2,NF)] += 0;

	Cmat[getindex2(1,0,NF)] += 0;
	Cmat[getindex2(1,1,NF)] += -nus;
	Cmat[getindex2(1,2,NF)] += 0;

	return;
}
