#include "eigen.h"



#ifdef MLIN
double bump_function(double rval) {
	
	
	double delta1 = 5.0 * .05 * pow(1.0,1.5-.5*flare_index);
	double delta2 = 5.0 * .05 * pow(2.0,1.5-.5*flare_index);

	double fac1 = .5*(1-.1)*(1+tanh((rval-1.0)/delta1))+.1;
	double fac2 = .5*(1-.1)*(1-tanh((rval-2.0)/delta2))+.1;
	
	return (fac1*fac2);

}
#endif

	
	
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

#ifdef EXACTKAPPA 

double sigma_function(double rval) {
//	return pow(rval,sigma_index);
//	return pow(rval,sigma_index) * exp(-rval/rdecay);
	return pow(rval,sigma_index) * exp(-pow(rval/rdecay,2));
}

double dlogsigma_dlogr(double rval) {
	
//	return sigma_index;	
//	return sigma_index - rval/r_decay;
	return sigma_index - 2*pow(rval/rdecay,2);
}

double dlogsigma2_dlogr2(double rval) {

//	return 0;	
//	return -rval/r_decay;
	return -4*pow(rval/rdecay,2);
}


#endif

