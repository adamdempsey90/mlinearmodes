#include "eigen.h"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>


double integrand(double x, void *params);
void add_edge_sg(double complex *mat);

void read_kernel(void) {
	int i,j;
	FILE *f;

	f=fopen("kernel.dat","r");

	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			if (!fscanf(f, "%lg", &kernel[j+N*i]))
          		break;

		}
	}
	fclose(f);

	return;

}

double sg_integrand(double r1, double rp1, double eps1) {
	double res;
	double r2,r3,r4,r5,r6,r7,r8,rp2,rp3,rp4,rp5,rp6,rp7,rp8,eps2,eps4,eps6,eps8,kval,ek,ee;



		r2 = r1*r1;
		r3 = r2*r1;
		r4 = r2*r2;
		r5 = r4*r1;
		r6 = r5*r1;
		r7 = r6*r1;
		r8 = r7*r1;

		rp2 = rp1*rp1;
		rp3 = rp2*rp1;
		rp4 = rp2*rp2;
		rp5 = rp4*rp1;
		rp6 = rp5*rp1;
		rp7 = rp6*rp1;
		rp8 = rp7*rp1;


		eps2 = eps1*eps1;
		eps4 = eps2*eps2;
		eps6 = eps4*eps2;
		eps8 = eps6*eps2;


#ifdef CONSTSOFT
	kval = sqrt(4*r1*rp1/((r1+rp1)*(r1+rp1)+ eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	res = (2*(r2-rp2)*(r2-rp2)*(r4-r2*rp2+rp4)+(r2+rp2)*(7*r4-18*r2*rp2+7*rp4)*eps2
					+9*(r4+rp4)*eps4+5*(r2+rp2)*eps6+eps8)*-2*ee;
	res -= ( r_m_rp +eps2)*(2*(r2-rp2)*(r2-rp2)*(r2+rp2)+(5*r4+4*r2*rp2+5*rp4)*eps2
					+4*(r2+rp2)*eps4+eps6)*-2*ek;
  res /= (pow(r_m_rp+eps2,2)*pow(r_p_rp+eps2,1.5));

#else

		kval = sqrt(4*r1*rp1/((r1+rp1)*(r1+rp1)+ eps2*r1*rp1));
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

		res = -ee*(8*r8 + 8*rp8 + 32*r7*rp1*eps2 + 32*r1*rp7*eps2
					+2*r5*rp3*eps2 *(-24 + 17 *eps4) + 2*r3*rp5*eps2 *(-24 + 17*eps4)
					+r6*rp2*(-24 + 49*eps4) + r2*rp6*(-24 + 49*eps4)
					+r4*rp4 *(32 + 22*eps4 + 9*eps8));

		res += ((r1 - rp1)*(r1-rp1) + r1*rp1*eps2)*(8*(r2 - rp2)*(r2-rp2)*(r2 + rp2)
						+ 8*r1*rp1*(3*r4 + r2*rp2 + 3*rp4)*eps2
						+ 25*r2*rp2*(r2 + rp2)*eps4 + 9*r3*rp3*eps6)*ek;

		res /= (2*pow((r1-rp1)*(r1-rp1)+r1*rp1*eps2,2)*pow((r1+rp1)*(r1+rp1)+r1*rp1*eps2,1.5));
#endif


	return res;
}

void compute_kernels(void) {

#ifdef READKERNEL
	read_kernel();
#else

	int i,j, indx;
	double eps1;

#ifdef OPENMP
#pragma omp parallel private(indx,i,j,eps1) shared(kernel,N,r,scaleH)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;

#ifdef CONSTSOFT
		eps1 = eps;
#else
#ifdef SYMSOFT
		eps1 = eps * sqrt(scaleH[i]*scaleH[i] + scaleH[j]*scaleH[j]);
#else
#ifdef SYMSOFT2
//		eps1  = eps*sqrt(scaleH[i]*scaleH[j]);
		eps1 = eps;
#else
		eps1 = eps * scaleH[j];
#endif
#endif
#endif

		kernel[indx] = sg_integrand(r[i],r[j],eps1);

	}
#endif
	return;
}

void add_sg(double complex *mat, double complex *bcmat) {
	int i,j;
	int indx,mindx;
	double complex G;
	double lambda = 1;

	for(i=0;i<N;i++) {
		G = lambda*.5/ (omega[i]*r[i]*r[i]*r[i]);

#ifdef INDIRECT
		mat[ncols*i] -= G *3*M_PI*r[i]*r[i]*sigma[0];
		mat[N-1 + ncols*i] += G*3*M_PI*r[i]*r[i]*sigma[N-1];
#endif
		for(j=0;j<N;j++) {
			mindx = j + ncols*i;
			indx = j+N*i;
			mat[mindx] += G*kernel[indx]*weights[j]*r[j]*sigma[j];
		}


	}
	return;
}




void add_edge_sg(double complex *mat) {
	int i,j;
	double x1,x2,x4,kval,kern,rd,eps2,ee,ek, G;

	rd = r[N-1];

	j = N-1;

#ifdef OPENMP
#pragma omp parallel private(i,x1,x2,x4,kval, kern,eps2,ee,ek, G) shared(r,scaleH,eps,rd,mat,N)
#pragma omp for schedule(static)
#endif
	for(i=0;i<N;i++) {
		G = 1.0/(2*omega[i]*pow(r[i],3));
		x1 = r[i]/rd;
		x2 = x1*x1;
		x4 = x2*x2;

#ifdef CONSTSOFT
		eps2 = eps;
#else
		eps2 = eps * scaleH[N-1]/rd;
#endif

		eps2 *= eps2;

		kval = sqrt(4*x1/(eps2 + (1+x1)*(1+x1)));
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

		kern = -2*( ((1+eps2)*(1+eps2)+ 3*(eps2-1)*x2 + 2*x4)*ee - (eps2+(x1-1)*(x1-1))*(1+eps2+2*x2)*ek);
		kern /= ((eps2 + (x1-1)*(x1-1))*sqrt(eps2+(1+x1)*(1+x1)));

		mat[j+N*i] += G*rd*rd*sigma[N-1]*kern;
	}


	return;

}

double integrand(double x, void *params) {
	double r1,r2,r3,r4,r6,rp1,rp2,rp3,rp4,rp6,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee;
	double ans,norm;
//	int i=  *(int *)params;
	r1 = * (double *)params;
//	r1 = r[i];
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r2*r2;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp3 = rp1*rp2;
	rp4 = rp2*rp2;
	r6 = r4*r2;
	rp6 = rp4*rp2;

	eps1 = eps;

// #ifdef CONSTSOFT
// 	eps1 = eps;
// #else
// #ifdef SYMSOFT
// 	eps1 = eps * sqrt(scaleH_func(r1)*scaleH_func(r1) + scaleH_func(rp1)*scaleH_func(rp1));
// #else
// #ifdef SYMSOFT2
// //	eps1  = eps*sqrt(scaleH_func(r1)*scaleH_func(rp1));
// 	eps1 = eps;
// #else
// #ifdef LINEARSOFT
// 	eps1 = eps*rp1;
// #else
// 	eps1 = eps * scaleH_func(rp1);
// #endif
// #endif
// #endif
// #endif
	eps2 = eps1*eps1;
	eps4 = eps2*eps2;
	r_p_rp = r1 + rp1;
	r_p_rp *= r_p_rp;
	r_m_rp = r1 - rp1;
	r_m_rp *= r_m_rp;

#ifdef CONSTSOFT

	kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	norm = 2 * pow( (r1-rp1)*(r1-rp1) + eps2,-2)*pow( (r1+rp1)*(r1+rp1) + eps2,-1.5);

	ans = ee*(r6 + pow(rp2 + eps2,3) - r4*(rp2 + 5*eps2) - r2*(rp2 + eps2)*(rp2 + 5*eps2))
	-ek*(pow(r1 - rp1,2) + eps2)*(r4 -2*r2*rp2 + pow(rp2 + eps2,2));

#else
	kval = sqrt(4*r1*rp1/(r_p_rp + r1*rp1*eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	norm = pow( (r1-rp1)*(r1-rp1) + eps2*r1*rp1,-2)*pow( (r1+rp1)*(r1+rp1) + eps2*r1*rp1,-1.5);

	ans = 2*ee*(r6 + +rp6 - 8*r3*rp3*eps2  - ( r4*rp2 +r2*rp4)*(1+eps4))
	-ek*((r1-rp1)*(r1-rp1) + rp1*r1*eps2)*(2*(r2-rp2)*(r2-rp2) + 2*r1*rp1*(r2+rp2)*eps2+r2*rp2*eps4);

#endif



	ans *= norm;


	ans *= rp2*sigma_func(rp1);

	return ans;

}

void calc_omega_prec_grav(void) {


	int nsubintervals = 3000;
	double error, res, norm,a;
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( nsubintervals );
	gsl_function func;
	int i;

	if (analytic_potential()) {
		for(i=0;i<N;i++) {
			omega_prec[i] += omega_prec_grav_analytic(r[i]);
		}
#ifdef PLANETS
			for(i=0;i<NP;i++) {
				Planets[i].wp  += omega_prec_grav_analytic(Planets[i].position);
			}

#endif
	}
	else {

		func.function = &integrand;

		for(i=0;i<N;i++) {

			func.params = &r[i];
			gsl_integration_qags(&func, log(r[0]/10),log(r[N-1]*10),0,tol, nsubintervals , workspace,&res,&error);
			norm = 2*sqrt(r[i]);
			omega_prec[i] += res/norm;

		}
#ifdef PLANETS
	for(i=0;i<NP;i++) {
		a = Planets[i].position;
		func.params = &a;
		gsl_integration_qags(&func,log(r[0]/10),log(r[N-1]*10),0,tol, nsubintervals , workspace,&res,&error);
		norm = 2*omk_func(a)*a*a;
		Planets[i].wp += res/norm;
	}

#endif
	}

	gsl_integration_workspace_free(workspace);

	return;

}
