#include "eigen.h"

void alloc_globals(void) {
	
	D = (double *)malloc(sizeof(double)*N*N);
	D2 = (double *)malloc(sizeof(double)*N*N);
	
	weights = (double *)malloc(sizeof(double)*N);
	KD = (double complex*)malloc(sizeof(double complex)*N*N);
	DKD = (double complex*)malloc(sizeof(double complex)*N*N);
	H = (double complex *)malloc(sizeof(double complex)*N*N);
	HL = (double complex*)malloc(sizeof(double complex)*N*N);
	cwork = (double complex *)malloc(sizeof(double complex)*N*N);
	work = (double *)malloc(sizeof(double)*N*N);
	
	r = (double *)malloc(sizeof(double)*N);
	lr = (double *)malloc(sizeof(double)*N);
	c2 = (double *)malloc(sizeof(double)*N);
	sigma = (double *)malloc(sizeof(double)*N);
	omega = (double *)malloc(sizeof(double)*N);
	omega2 = (double *)malloc(sizeof(double)*N);
	
	kappa = (double complex *)malloc(sizeof(double complex)*N);
	kappa2 = (double *)malloc(sizeof(double)*N);
	scaleH = (double *)malloc(sizeof(double)*N);
	dlds = (double *)malloc(sizeof(double)*N);
	dldc2 = (double *)malloc(sizeof(double)*N);
	lsig = (double *)malloc(sizeof(double)*N);
	lc2 = (double *)malloc(sizeof(double)*N);
	d2lds = (double *)malloc(sizeof(double)*N);
	lom = (double *)malloc(sizeof(double)*N);
	dldom = (double *)malloc(sizeof(double)*N);
	d2dom = (double *)malloc(sizeof(double)*N);

	pres = (double *)malloc(sizeof(double)*N);
	lpres = (double *)malloc(sizeof(double)*N);
	temp = (double *)malloc(sizeof(double)*N);
	ltemp = (double *)malloc(sizeof(double)*N);
	dldpres = (double *)malloc(sizeof(double)*N);
	d2ldpres = (double *)malloc(sizeof(double)*N);
	dldtemp = (double *)malloc(sizeof(double)*N);
	d2ldtemp = (double *)malloc(sizeof(double)*N);


	nu = (double *)malloc(sizeof(double)*N);
	dldnu = (double *)malloc(sizeof(double)*N);
	lnu = (double *)malloc(sizeof(double)*N);


	return;

}
void free_globals(void) {

	free(D);
	free(D2);
	free(kernel);
	free(kernel0);
	free(weights);
	free(KD);
	free(DKD);
	free(H);
	free(HL);
	free(cwork);
	free(work);
	
	free(r);
	free(lr);
	free(c2);
	free(sigma);
	free(omega);
	free(kappa);
	free(omega2);
	free(omegap2);
	free(kappa2);
	free(scaleH);
	free(dlds);
	free(dldc2);
	free(lsig);
	free(lc2);
	free(dphi0dr);
	free(d2lds);
	free(lom);
	free(dldom);
	free(d2dom);
	
//#ifndef ISOTHERMAL
	free(pres);
	free(lpres);
	free(temp);
	free(d2ldpres);
	free(dldpres);
	free(ltemp);
	free(d2ldtemp);
	free(dldtemp);
//#endif
	
	free(nu); 
	free(dldnu);
	free(lnu);



	return;

}
