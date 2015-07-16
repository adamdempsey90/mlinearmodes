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




double integrand1(double x, void *params);
double integrand2(double x, void *params);
double integrandom2(double x, void *params);
double integrandkap2(double x, void *params);
void calc_kernels(double *kernel1, double *kernel2);
void integrateK(IntParams params,double *ans1, double *ans2);

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

void integrateK(IntParams params,double *ans1, double *ans2)  {
	char errorstr[MAXSTRLEN];
	int nsubintervals = 1000;
	int status;
	double error, res1, res2;
	double r1 = params.r;
	double rp1 = params.rp;
	gsl_set_error_handler_off();
	gsl_integration_workspace *wksp1 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_integration_workspace *wksp2 = gsl_integration_workspace_alloc( nsubintervals );
	gsl_function func1, func2;

	func1.function = &integrand1;
	func1.params = &params;
	func2.function = &integrand2;
	func2.params = &params;

	status = gsl_integration_qag(&func1,0,M_PI,0,tol, nsubintervals ,6, wksp1,&res1,&error);
	if (status) {
		sprintf(errorstr,"kernel1, r=%lg, rp=%lg\n",r1,rp1);
		gsl_integration_error_handler(status,errorstr);
	}
	status =gsl_integration_qag(&func2,0,M_PI,0,tol, nsubintervals ,6, wksp2,&res2,&error);
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

void calc_kernels(double *kernel1, double *kernel2) {
	int i,j;
	double res1, res2, wghti, wghtj;
	IntParams params;
	params.m = mval;
	params.eps = eps;
	printf("Calculating Kernels...\n");
// OPENMP
	for(i=0;i<N;i++) {
//		wghti = -weights[i]*r[i]*r[i];
		params.r = r[i];
		for(j=0;j<N;j++) {
			params.rp = r[j];
			wghtj = -weights[j]*r[j]*r[j];
			integrateK(params,&res1,&res2);
			kernel1[j + i*N] = res1*wghtj;
//			kernel1[i + j*N] = res1*wghti; 		// Overwrite with same value for i=j
			kernel2[j + i*N] = res2*wghtj;
//			kernel2[i + j*N] = res2*wghti;
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

	res = -(4*rp1*(-rp1*(r2 - rp2)*(r2-rp2) + r1*(r4 + 4*r2*rp2 - rp4)*eps1
					+ r4*rp1*eps2)*ee + ((r1 - rp1)*(r1-rp1)
					+ r1*rp1*eps1)*(4*(r2 - rp2)*(r2-rp2) + 6*r1*rp1*(r2 + rp2)*eps1
					 + 3*r2*rp2*eps2)*ek);
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
	res = 2*(r1-rp1)*(r1+rp1)*ee + 2*((r1-rp1)*(r1-rp1)+eps1*r1*rp1)*ek;

	wght = rp2*sigma_func(rp1);

	return wght*res/denom;

}

void calc_omega_sg(double *omg2, double *kapg2) {
	int i,status;
	int nsubintervals = 2*N;
	double error, res1, res2, lowerbound, upperbound;
	IntParams params;
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

	params.m = 0;
	params.eps = eps;
	params.rp = 0;

	printf("Calculating omega and kappa corrections from self gravity...\n");
// OPENMP
	for(i=0;i<N;i++) {
		params.r = r[i];
		func1.params = &params;
		func2.params = &params;
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
			fprintf(stderr, "omega_g^2 < 0!\n omp2=%lg, r=%lg\n",res1,r[i]);
		}
		if (res2 < 0) {
			fprintf(stderr, "kappa_g^2 < 0!\n omp2=%lg, r=%lg\n",res1,r[i]);
		}
		omg2[i] = res1;
		kapg2[i] = res2;

	}

	gsl_integration_workspace_free(wksp1);
	gsl_integration_workspace_free(wksp2);
	return;




}

void add_sg(double complex *mat, double complex *bcmat) {
	double *kernel1 = (double *)malloc(sizeof(double)*N*N);
	double *kernel2 = (double *)malloc(sizeof(double)*N*N);
	calc_kernels(kernel1, kernel2);
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
