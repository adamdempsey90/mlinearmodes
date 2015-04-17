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
void compute_kernels(void) {

#ifdef READKERNEL
	read_kernel();
#else

	int i,j, indx;
	double r1,r2,r4,rp1,rp2,rp4,eps1,eps2,eps4,eps6,eps8,r_p_rp,r_m_rp,kval,ek,ee; 
	
#ifdef OPENMP
#pragma omp parallel private(indx,i,j,r1,r2,r4,rp1,rp2,rp4,eps1,eps2,eps4,eps6,eps8,r_p_rp,r_m_rp,kval,ek,ee) shared(kernel,N,r,scaleH,weights)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		

		r1 = r[i];
		r2 = r1*r1;
		r4 = r2*r2;
		rp1 = r[j];
		rp2 = rp1*rp1;
		rp4 = rp2*rp2;

#ifdef CONSTSOFT
		eps1 = eps;
#else
#ifdef SYMSOFT
		eps1 = eps * sqrt(scaleH[i]*scaleH[i] + scaleH[j]*scaleH[j]);
#else
		eps1 = eps *scaleH[j];
#endif
#endif

		eps2 = eps1*eps1;
		eps4 = eps2*eps2;
		eps6 = eps4*eps2;
		eps8 = eps6*eps2;
		
		r_p_rp = (r1+rp1)*(r1+rp1);
		r_m_rp = (r1-rp1)*(r1-rp1);
		
		
		kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
	
	
		kernel[indx] = (2*(r2-rp2)*(r2-rp2)*(r4-r2*rp2+rp4)+(r2+rp2)*(7*r4-18*r2*rp2+7*rp4)*eps2
						+9*(r4+rp4)*eps4+5*(r2+rp2)*eps6+eps8)*-2*ee;
		kernel[indx] -= ( r_m_rp +eps2)*(2*(r2-rp2)*(r2-rp2)*(r2+rp2)+(5*r4+4*r2*rp2+5*rp4)*eps2
						+4*(r2+rp2)*eps4+eps6)*-2*ek;
		kernel[indx] /= (pow(r_m_rp+eps2,2)*pow(r_p_rp+eps2,1.5));
			
		
	//	kernel[indx] *=  weights[j]*rp1*sigma[j];



		
	}	
#endif
	return;
}

void add_sg(double complex *mat, double complex *bcmat) {	
	int i,j;
	int indx;
	double complex G;
	
	
	for(i=0;i<N;i++) {
		G = .5/ ( omega[i]*r[i]*r[i]*r[i]);
		for(j=0;j<N;j++) {
			indx = j + N*i;
			mat[indx] += G*kernel[indx]*weights[j]*r[j]*sigma[j];		
		}
			
			
	}			
	return;
}




// void compute_kernels(void) {
// 
// #ifdef READKERNEL
// 	read_kernel();
// #else
// 
// 	int i,j, indx;
// 	double r1,r2,r4,rp1,rp2,rp4,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee; 
// 	
// #ifdef OPENMP
// #pragma omp parallel private(indx,i,j,r1,r2,r4,rp1,rp2,rp4,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee) shared(kernel,N,r,scaleH,weights)
// #pragma omp for schedule(static)
// #endif
// 	for(indx=0;indx<N*N;indx++) {
// 		i = indx/N;
// 		j = indx - i*N;
// 		
// 
// 		r1 = r[i];
// 		r2 = r1*r1;
// 		r4 = r2*r2;
// 		rp1 = r[j];
// 		rp2 = rp1*rp1;
// 		rp4 = rp2*rp2;
// 
// #ifdef CONSTSOFT
// 		eps1 = eps;
// #else
// #ifdef SYMSOFT
// 		eps1 = eps * sqrt(scaleH[i]*scaleH[i] + scaleH[j]*scaleH[j]);
// #else
// 		eps1 = eps *scaleH[j];
// #endif
// #endif
// 
// 		eps2 = eps1*eps1;
// 		eps4 = eps2*eps2;
// 		r_p_rp = r[i] + r[j];
// 		r_p_rp *= r_p_rp;
// 		r_m_rp = r[i] - r[j];
// 		r_m_rp *= r_m_rp;
// 		
// 		kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
// 		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
// 		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
// 			
// 		kernel[indx] = -2*(eps4 + 2*r4 - 3*r2*rp2 + rp4 + eps2*(3*r2 + 2*rp2))*ee;
// 		kernel[indx] += 2*(eps2 + r_m_rp)*(eps2 + 2*r2 + rp2)*ek;
// 		kernel[indx] /= ((eps2 + r_m_rp)*rp1*sqrt(eps2 +r_p_rp));
// 		
// 		kernel[indx] *=  weights[j]*rp2;
// 
// 
// 
// 		
// 	}	
// #endif
// 	return;
// }
// 
// void add_sg(double complex *mat, double complex *bcmat) {	
// 	int i,j;
// 	int l,indx;
// 	int il, lj;
// 	double complex res,G;
// 	
// 	
// 	for(i=0;i<N;i++) {
// 		G = .5/ ( omega[i]*r[i]*r[i]*r[i]);
// 		for(j=0;j<N;j++) {
// 			indx = j + N*i;
// 			res = 0;
// 			for(l=0;l<N;l++) {
// 				il = l + N*i;
// 				lj = j + N*l;
// 				res -= sigma[i]*(dlds[l]*Identity[il] + D[il])*kernel[lj];
// 			}
// 			mat[indx] += G*res;		
// 		}
// 			
// 			
// 	}			
// 
// #ifdef EDGEGRAVITY
// 		add_edge_sg(mat);
// #endif
// 	
// 	return;
// 	
// }	

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
	double r1,r2,r4,rp1,rp2,rp4,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee; 
	double ans;
	int i=  *(int *)params;
	
	r1 = r[i];
	r2 = r1*r1;
	r4 = r2*r2;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp4 = rp2*rp2;
#ifdef CONSTSOFT
	eps1 = eps;
#else
#ifdef SYMSOFT
	eps1 = eps * sqrt(scaleH_func(r1)*scaleH_func(r1) + scaleH_func(rp1)*scaleH_func(rp1));
#else
	eps1 = eps * scaleH_func(rp1);
#endif
#endif
	eps2 = eps1*eps1;
	eps4 = eps2*eps2;
	r_p_rp = r1 + rp1;
	r_p_rp *= r_p_rp;
	r_m_rp = r1 - rp1;
	r_m_rp *= r_m_rp;
	
	kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
		
	ans=2*(eps2 - 2 *eps1 *r1 - r2 + rp2)*(eps2 + 2 *eps1*r1 - r2 + rp2)*(eps2 + r2 + rp2)*ee; 
	ans -= 2*(eps2 + (r1 - rp1)*(r1-rp1))*(eps4 + r4 + 2*(eps1 - r1)*(eps1 + r1)*rp2 + rp4)*ek;
	ans /= (pow(eps2 + (r1 - rp1)*(r1-rp1),2) *pow(eps2 + (r1 + rp1)*(r1+rp1),1.5));
	
	ans *= -rp2*sigma_func(rp1);
	
	return ans;

}

void calc_omega_prec_grav(void) {


	int nsubintervals = 1000;
	double error, res;
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( nsubintervals );
	gsl_function func;
	int i;

	if (analytic_potential()) {
		for(i=0;i<N;i++) {
			omega_prec[i] += omega_prec_grav_analytic(r[i]);
		}
		
	}
	else { 

		func.function = &integrand;
	
		for(i=0;i<N;i++) {
	
			func.params = &i;
			gsl_integration_qags(&func, log(r[0]),log(r[N-1]),0,tol, nsubintervals , workspace,&res,&error);
		
			omega_prec[i] += res*(-.5/sqrt(r[i]));
		
		}
	}
	
	gsl_integration_workspace_free(workspace);
	
	
	return; 

}