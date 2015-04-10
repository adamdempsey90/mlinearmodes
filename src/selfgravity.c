#include "eigen.h"


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

	f=fopen("kernel0.dat","r");
	
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			if (!fscanf(f, "%lg", &kernel0[j+N*i])) 
          		break;
			
		}
	}
	fclose(f);
	
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
		eps2 = eps*scaleH[N-1]/rd;
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
	double r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee; 
	double ans;
	int i=  *(int *)params;
	
	r1 = r[i];
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r2*r2;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp3 = rp2*rp1;
	rp4 = rp2*rp2;
	eps1 = eps * scaleH_func(r1);
	eps2 = eps1*eps1;
	eps4 = eps2*eps2;
	eps6 = eps4 * eps2;
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
	double error;
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( nsubintervals );
	int i,j;

	
	
	gsl_function func;
	func.function = &integrand;
	
	for(i=0;i<N;i++) {
	
		func.params = &i;
		gsl_integration_qags(&func, log(r[0]),log(r[N-1]),0,tol, nsubintervals , workspace,&omega_prec[i],&error);
		
		omega_prec[i] *= (-.5/sqrt(r[i]));
		
	}

	gsl_integration_workspace_free(workspace);
	
	
	
	
	return; 

}