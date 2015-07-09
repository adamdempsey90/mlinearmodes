#include "eigen.h"


void alloc_globals(void) {

	D = (double *)malloc(sizeof(double)*N*N*NF*NF);
	D2 = (double *)malloc(sizeof(double)*N*N*NF*NF);
	Identity = (double *)malloc(sizeof(double)*N*N*NF*NF);
	kernel = (double *)malloc(sizeof(double)*N*N);
	weights = (double *)malloc(sizeof(double)*N);


	cwork = (double complex *)malloc(sizeof(double complex)*nrows*ncols);
	work = (double *)malloc(sizeof(double)*nrows*ncols);

	r = (double *)malloc(sizeof(double)*N);
	lr = (double *)malloc(sizeof(double)*N);

	scaleH = (double *)malloc(sizeof(double)*N);
	aspect_ratio = (double *)malloc(sizeof(double)*N);

	omega = (double *)malloc(sizeof(double)*N);
	dldom = (double *)malloc(sizeof(double)*N);
	d2ldom = (double *)malloc(sizeof(double)*N);
	kappa2 = (double *)malloc(sizeof(double)*N);

	sigma = (double *)malloc(sizeof(double)*N);
	dlds = (double *)malloc(sizeof(double)*N);
	d2lds = (double *)malloc(sizeof(double)*N);

	pres = (double *)malloc(sizeof(double)*N);
	dldpres = (double *)malloc(sizeof(double)*N);
	d2ldpres = (double *)malloc(sizeof(double)*N);

	temp = (double *)malloc(sizeof(double)*N);
	dldtemp = (double *)malloc(sizeof(double)*N);
	d2ldtemp = (double *)malloc(sizeof(double)*N);

	c2 = (double *)malloc(sizeof(double)*N);
	dldc2 = (double *)malloc(sizeof(double)*N);
	d2ldc2 = (double *)malloc(sizeof(double)*N);

	coeffs_A =  (double complex *)malloc(sizeof(double complex)*N*NF*NF);
	coeffs_B =  (double complex *)malloc(sizeof(double complex)*N*NF*NF);
	coeffs_C =  (double complex *)malloc(sizeof(double complex)*N*NF*NF);

	return;

}
void free_globals(void) {

	free(D);
	free(D2);
	free(Identity);
	free(kernel);
	free(weights);


	free(cwork);
	free(work);

	free(r );
	free(lr);

	free(scaleH);
	free(aspect_ratio);

	free(omega );
	free(dldom);
	free(d2ldom );
	free(kappa2);

	free(sigma);
	free(dlds);
	free(d2lds);

	free(pres );
	free(dldpres);
	free(d2ldpres );

	free(temp);
	free(dldtemp );
	free(d2ldtemp);

	free(c2);
	free(dldc2);
	free(d2ldc2);


	free(coeffs_A);
	free(coeffs_B);
	free(coeffs_C);



	return;

}
