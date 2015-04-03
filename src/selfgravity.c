#include "eigen.h"


void compute_kernels(void) {
	int i,j, indx;
	double r2,r4,rp2,rp4,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee; 
	
#ifdef OPENMP
#pragma omp parallel private(indx,i,j,r2,r4,rp2,rp4,eps1,eps2,eps4,r_p_rp,r_m_rp,kval,ek,ee) shared(kernel,kernel0,N,r,scaleH,weights)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		
		
		
		
		r2 = r[i]*r[i];
		r4 = r2*r2;
		rp2 = r[j]*r[j];
		rp4 = rp2*rp2;
		eps1 = eps * scaleH[j];
		eps2 = eps1*eps1;
		eps4 = eps2*eps2;
		r_p_rp = r[i] + r[j];
		r_p_rp *= r_p_rp;
		r_m_rp = r[i] - r[j];
		r_m_rp *= r_m_rp;
		
		kval = 4*r[i]*r[j]/(r_p_rp + eps1);
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
			
		kernel[indx] = -(eps4 + 2*r4 - 3*r2*rp2 + rp4 + eps2*(3*r2 + 2*rp2))*ee;
		kernel[indx] += (eps2 + r_m_rp)*(eps2 + 2*r2 + rp2)*ek;
		kernel[indx] /= ((eps2 + r_m_rp)*r[j]*sqrt(eps2 +r_p_rp));
		kernel[indx] *= 2;
		
		kernel0[indx] = 2*((eps2 - r2+rp2)*ee - (eps2 + r_m_rp)*ek);
		kernel0[indx] /= (r[i]*(eps2 + r_m_rp)*sqrt(eps2 + r_p_rp));
		
		
		kernel[indx] *= weights[j]*r[j]*r[j];
		kernel0[indx] *= weights[j]*r[j]*r[j];
		
	}	

	return;
}




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
	