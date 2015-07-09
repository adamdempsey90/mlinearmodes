#include "eigen.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
#ifdef OPENMP
#include <omp.h>
#endif

typedef struct IntParams {
	int m;
	double r,rp,eps;
} IntParams;

double *kernel1, *kernel2;



double integrand1(double x, void *params);
double integrand2(double x, void *params);
double integrandom2(double x, void *params);
double integrandkap2(double x, void *params);



double integrand1(double x, void *params) {
	IntParams pars = * (IntParams *)params;
	double r1 = pars.r;
	double rp1 = pars.rp;
	double eps1 = pars.eps;
	int m1 = pars.m;


	return 2*cos(m1*x)/sqrt(r1*r1 + rp1*rp1 - 2*r1*rp1*cos(x) + eps1*r1*rp1);

}

double integrand2(double x, void *params) {
	IntParams pars;
	pars = *(IntParams *)params;
	double r1 = pars.r;
	double rp1 = pars.rp;
	double eps1 = pars.eps;
	int m1 = pars.m;

	double fac = (2*r1 + rp1*eps1 - 2*rp1*cos(x))*cos(m1*x);

	return -fac*pow(r1*r1+rp1*rp1-2*r1*rp1*cos(x) + eps1*r1*rp1,-1.5);
}

void integrateK(double r1, double rp1, double *ans1, double *ans2) {
	char errorstr[MAXSTRLEN];
	int nsubintervals = 3000;
	int status;
	double error, res1, res2;
	IntParams params1, params2;
	gsl_set_error_handler_off();
	gsl_integration_workspace *wksp1 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_integration_workspace *wksp2 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_function func1, func2;

	params1.m = mval;
	params1.eps = eps;
	params1.r = r1;
	params1.rp = rp1;

	params2.m = mval;
	params2.eps = eps;
	params2.r = r1;
	params2.rp = rp1;

	func1.function = &integrand1;
	func1.params = &params1;
	func2.function = &integrand2;
	func2.params = &params2;

	status = gsl_integration_qags(&func1,0,M_PI,0,tol, nsubintervals , wksp1,&res1,&error);
	if (status) {
		sprintf(errorstr,"kernel1, r=%lg, rp=%lg\n",r1,rp1);
		gsl_integration_error_handler(status,errorstr);
	}
	status =gsl_integration_qags(&func2,0,M_PI,0,tol, nsubintervals , wksp2,&res2,&error);
	if (status) {
		sprintf(errorstr,"kernel2, r=%lg, rp=%lg\n",r1,rp1);
		gsl_integration_error_handler(status,errorstr);
	}



	*ans1 = res1;
	*ans2 = res2;

	gsl_integration_workspace_free(wksp1);
	gsl_integration_workspace_free(wksp2);
	return;
}

void calc_kernels(void) {
	int i,j;
	double res1, res2, wghti, wghtj;

	printf("Calculating Kernels...\n");
// OPENMP
	for(i=0;i<N;i++) {
		wghti = -weights[i]*r[i]*r[i];
		for(j=i;j<N;j++) {
			wghtj = -weights[j]*r[j]*r[j];
			integrateK(r[i],r[j],&res1,&res2);
			kernel1[j + i*N] = res1*wghtj;
			kernel1[i + j*N] = res1*wghti; 		// Overwrite with same value for i=j
			kernel2[j + i*N] = res2*wghtj;
			kernel2[i + j*N] = res2*wghti;
		}
	}
/* Indirect Pontential for m=1 */

	if (mval==1) {
		printf("\tIncluding Indirect Potential for m=1 ...\n");
		for(i=0;i<N;i++) {
			for(j=0;j<N;j++) {
				kernel1[j+i*N] += M_PI*r[i]*weights[j];
				kernel2[j+i*N] += M_PI*weights[j];
			}
		}
	}

	return;
}


double integrandkap2(double x,void *params) {
	double r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,kval,ee,ek;
	double res, denom, wght;
	IntParams pars = * (IntParams *)params;
	r1 = pars.r;
	r2 = r1*r1;
	r3 = r1*r2;
	r4 = r1*r3;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp3 = rp1*rp2;
	rp4 = rp1*rp3;
	eps1 = pars.eps;
	eps2 = eps1*eps1;

	kval = sqrt(4*r1*rp1/(eps1*r1*rp1 + (r1+rp1)*(r1+rp1)));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	denom = r2*pow((r1-rp1)*(r1-rp1)+eps1*r1*rp1,2)*pow(eps1*r1*rp1+(r1+rp1)*(r1+rp1),1.5);
	res = 4*rp1*(rp1*(r2-rp2)*(r2-rp2)+r1*(-r4+-4*r2*rp2+rp4)*eps1-r4*rp1*eps2)*ee;
	res -= ((r1-rp1)*(r1-rp1)+eps1*r1*rp1)*(4*(r2-rp2)*(r2-rp2)+6*r1*rp1*(r2+rp2)*eps1+ 3*eps2*r2*rp2)*ek;

	wght = -rp2*sigma_func(rp1);

	return wght*res/denom;

}
double integrandom2(double x,void *params) {
	double r1,r2,rp1,rp2,eps1,kval,ee,ek;
	double res, denom,wght;
	IntParams pars = * (IntParams *)params;

	r1 = pars.r;
	r2 = r1*r1;

	rp1 = exp(x);
	rp2 = rp1*rp1;

	eps1 = pars.eps;



	kval = sqrt(4*r1*rp1/(eps1*r1*rp1 + (r1+rp1)*(r1+rp1)));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);


	denom = r2*((r1-rp1)*(r1-rp1)+r1*rp1*eps1)*sqrt((r1+rp1)*(r1+rp1)+eps1*r1*rp1);
	res = 2*(rp2-r2)*ee - 2*((r1-rp1)*(r1-rp1)+eps1*r1*rp1)*ek;
	wght = -rp2*sigma_func(rp1);

	// printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",r1,rp1,eps1,kval,ee,ek,denom,res,wght);

	return wght*res/denom;

}

void calc_omega_sg(double *omg2, double *kapg2) {
	int i,status;
	int nsubintervals = 2*N;
	double error, res1, res2, lowerbound, upperbound;
	IntParams params1, params2;
	char errorstr[MAXSTRLEN];
	gsl_set_error_handler_off();

	gsl_integration_workspace *wksp1 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_integration_workspace *wksp2 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_function func1, func2;

	lowerbound = log(r[0]);
	upperbound = log(r[N-1]);

#ifdef EXTENDINTEG
	lowerbound -= log(10.);
	upperbound += log(10.);
#endif


	func1.function = &integrandom2;
	func2.function = &integrandkap2;

	params1.m = 0;
	params1.eps = eps;
	params1.rp = 0;

	params2.m = 0;
	params2.eps = eps;
	params2.rp = 0;

	printf("Calculating omega and kappa corrections from self gravity...\n");
// OPENMP
	for(i=0;i<N;i++) {
		params1.r = r[i];
		params2.r = r[i];
		func1.params = &params1;
		func2.params = &params2;
		status=gsl_integration_qags(&func1,lowerbound,upperbound,0,tol, nsubintervals , wksp1,&res1,&error);
		if (status){
			sprintf(errorstr,"omega2, r=%lg, lower=%lg, upper=%lg,tol=%e\n",
										r[i],lowerbound,upperbound,tol);
			gsl_integration_error_handler(status,errorstr);
		}
		status=gsl_integration_qags(&func2,lowerbound,upperbound,0,tol, nsubintervals , wksp2,&res2,&error);
		if (status){
			sprintf(errorstr,"kappa2, r=%lg,lower=%lg, upper=%lg,tol=%e\n",
										r[i],lowerbound,upperbound,tol);
			gsl_integration_error_handler(status,errorstr);
		}

		if (res1 < 0) {
			fprintf(stderr, "omega_g^2 < 0!\n omp2=%lg, r=%lg",res1,r[i]);
		}
		if (res2 < 0) {
			fprintf(stderr, "kappa_g^2 < 0!\n omp2=%lg, r=%lg",res1,r[i]);
		}
		omg2[i] = res1;
		kapg2[i] = res2;

	}

	gsl_integration_workspace_free(wksp1);
	gsl_integration_workspace_free(wksp2);
	return;




}

void add_sg(double complex *mat, double complex *bcmat) {
	kernel1 = (double *)malloc(sizeof(double)*N*N);
	kernel2 = (double *)malloc(sizeof(double)*N*N);
	calc_kernels();
	int i,j,indx;

	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			indx = getindex2(i,j,N);
			mat[getindex4(i,j,0,2,NF,N)] += kernel2[indx]/(I*mval);
			mat[getindex4(i,j,1,2,NF,N)] += kernel1[indx]/r[i];
		}
	}


// output kernels here

	free(kernel1); free(kernel2);
	return;
}
