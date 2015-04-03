#include "eigen.h"


double *kernel,*kernel0, *kernel02;

void init_sg(void) {

	alloc_sg();

#ifdef READKERNEL
	read_kernel();
#else	
	compute_kernels();
	output_kernels();
#endif
		
	return;


}

void alloc_sg(void) {

	

	kernel = (double complex *)malloc(sizeof(double complex)*N*N);
	kernel0 = (double *)malloc(sizeof(double)*N*N);
	kernel02 = (double *)malloc(sizeof(double)*N*N);
	
	return;
}


void free_sg(void) {
	
	free(kernel);
	free(kernel0);
	free(kernel02);

	return;
}
void compute_kernels(void) {
	int i,j, indx;
	double r1,r2,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee; 
	
#ifdef OPENMP
#pragma omp parallel private(indx,i,j,r1,r2,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee) shared(kernel,kernel0,kernel02,N,r,scaleH,weights)
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
		rp3 = rp2*r[j];
		rp4 = rp2*rp2;
		eps1 = eps * scaleH[j];
		eps2 = eps1*eps1;
		eps4 = eps2*eps2;
		eps6 = eps4 * eps2;
		r_p_rp = r[i] + r[j];
		r_p_rp *= r_p_rp;
		r_m_rp = r[i] - r[j];
		r_m_rp *= r_m_rp;
		
		kval = 4*r1*rp1/(r_p_rp + eps2);
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
			
		kernel[indx] = -(eps4 + 2*r4 - 3*r2*rp2 + rp4 + eps2*(3*r2 + 2*rp2))*ee;
		kernel[indx] += (eps2 + r_m_rp)*(eps2 + 2*r2 + rp2)*ek;
		kernel[indx] /= ((eps2 + r_m_rp)*rp1*sqrt(eps2 +r_p_rp));
		kernel[indx] *= 2;
		
		kernel0[indx] = 2*((eps2 - r2+rp2)*ee - (eps2 + r_m_rp)*ek);
		kernel0[indx] /= (r1*(eps2 + r_m_rp)*sqrt(eps2 + r_p_rp));
		
		kernel02[indx] = -ee*(eps6 + eps4*(-2*r2+3*rp2)+pow(rp3-rp1*r2,2)+eps2*(-3*r4-4*r2*rp2+3*rp4));
		kernel02[indx] += ek*(eps2+r_m_rp)*(eps4+r_m_rp*r_p_rp + eps2*(r2+2*rp2));
		kernel02[indx] *= 4*r1;
		kernel02[indx] /= ( pow(eps2+r_m_rp,2) * pow(eps2+r_p_rp,1.5));
		
		kernel[indx] *= weights[j]*r[j]*r[j];
		kernel0[indx] *= weights[j]*r[j]*r[j];
		kernel02[indx] *= weights[j]*r[j]*r[j];
		
		
	}	

	return;
}


void output_kernel(void) {
	int i,j;
	FILE *f = fopen("kernel.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.20lg\t",creal(kernel[j+i*N]));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	
	f = fopen("kernel0.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.20lg\t",kernel0[j+i*N]);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	
	f = fopen("kernel02.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.20lg\t",kernel02[j+i*N]);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);

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
	
	f=fopen("kernel02.dat","r");
	
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			if (!fscanf(f, "%lg", &kernel02[j+N*i])) 
          		break;
			
		}
	}
	fclose(f);
	
	return;

}


void sg_omega_correction(void) {
	int i,j;
	
	double *omegag2 = (double *)malloc(sizeof(double)*N);
	double *kappag2 = (double *)malloc(sizeof(double)*N);
	
	for(i=0;i<N;i++) {

		
		omegag2[i] = 0;
		kappag2[i] = 0;
		for(j=0;j<N;j++) {
			omegag2[i] += kernel0[j+i*N]*sigma[j];
			kappag2[i] += kernel02[j+i*N]*sigma[j];
		}
		omegag2[i] /= r[i];
		kappag2[i] *= pow(r[i],-3);

		omega2[i] += omegag2[i];
		kappa2[i] += kappag2[i];
		
	}

	free(omegag2);
	free(kappag2);
	return;

}
	