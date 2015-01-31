#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>

#define NPHI 10000

#define OPENMP
double *r, *lr,  *weights, *kernel, *kernel0;
int N;
double eps = .1;
double h0 = .05;
double dlr,ri,ro;


void calc_weights(void);
double Kij(int i, int j);
double kernel_integrand(int m, double r, double rp, double rs, double phi);
double K0ij(int i, int j);

int main(int argc, char *argv[]) {

	int i,j,indx;

#ifdef OPENMP	
	omp_set_num_threads(8);
#endif
	N = atoi(argv[1]);
	ri = atof(argv[2]);
	ro = atof(argv[3]);
	
	
	dlr = (log(ro) - log(ri))/((float) N);
	
	printf("Using Nr=%d\nri=%lg\nro=%lg\ndlr=%lg\n",N,ri,ro,dlr);
	
	r = (double *)malloc(sizeof(double)*N);
	lr = (double *)malloc(sizeof(double)*N);
	weights = (double *)malloc(sizeof(double)*N);
	kernel = (double *)malloc(sizeof(double)*N*N);
	kernel0 = (double *)malloc(sizeof(double)*N*N);
	
	for(i=0;i<N;i++ ){
		lr[i] = log(ri) + i*dlr;
		r[i] = exp(lr[i]);
	}
	
	calc_weights();

#ifdef OPENMP
#pragma omp parallel private(indx,i,j) shared(kernel,N)
#pragma omp for schedule(static)
#endif
	for(indx = 0; indx < N*N; indx++) {
		i = indx/N;
		j = indx - N*i;

		kernel[indx] = Kij(i,j);
		kernel0[indx] = K0ij(i,j);
	}	
	
	
	FILE *f = fopen("kernel.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.12e\t",kernel[j+N*i]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
	
	f = fopen("kernel0.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.12e\t",kernel0[j+N*i]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
	 

	free(r); free(lr); free(weights); free(kernel); free(kernel0);
	return 1;

}

void calc_weights(void) {
/* O(N^(-4)) weights */
	int i;
	for(i=0;i<N;i++) {
		if (i==0 || i==N-1) {
			weights[i] = 3./8 * dlr;
		}
		else if (i==1 || i==N-2) {
			weights[i] = 7./6 * dlr;
		}
		else if (i==2 || i==N-3) {
			weights[i] = 23./24 * dlr;
		}
		else {
			weights[i] = dlr;
		}
	
	}
	
	return;

}

double Kij(int i, int j) {
	int indx;
	double dp = M_PI/NPHI;

	double result = 0;
	double rs2 = r[j]*r[j]*(1 + eps*eps*h0*h0);
	
/* phi_i = i*dp */
	result = (3./8)*(kernel_integrand(1,r[i],r[j],rs2,0) + kernel_integrand(1,r[i],r[j],rs2,dp*NPHI));
	result += (7./6)*(kernel_integrand(1,r[i],r[j],rs2,dp) + kernel_integrand(1,r[i],r[j],rs2,dp*(NPHI-1)));
	result += (23./24)*(kernel_integrand(1,r[i],r[j],rs2,dp*2) + kernel_integrand(1,r[i],r[j],rs2,dp*(NPHI-2)));
	
	for(indx=3;indx<NPHI-2;indx++) {
		result += kernel_integrand(1,r[i],r[j],rs2,dp*indx);
	}
	
	result *= dp*2;
	
	result = weights[j]*r[i]*r[i]*(r[j]*r[j]*result - M_PI*r[i]);
	
	return result;


}
double K0ij(int i, int j) {
	int indx;
	double dp = M_PI/NPHI;

	double result = 0;
	double rs2 = r[j]*r[j]*(1 + eps*eps*h0*h0);
	
	result = (3./8)*(kernel_integrand(0,r[i],r[j],rs2,0) + kernel_integrand(0,r[i],r[j],rs2,dp*NPHI));
	result += (7./6)*(kernel_integrand(0,r[i],r[j],rs2,dp) + kernel_integrand(0,r[i],r[j],rs2,dp*(NPHI-1)));
	result += (23./24)*(kernel_integrand(0,r[i],r[j],rs2,dp*2) + kernel_integrand(0,r[i],r[j],rs2,dp*(NPHI-2)));
	
	for(indx=3;indx<NPHI-2;indx++) {
		result += kernel_integrand(0,r[i],r[j],rs2,dp*indx);
	}
	
	result *= dp*2*r[j]*r[j]*weights[j];
	
	return result;


}


double kernel_integrand(int m, double r, double rp, double rs, double phi) {
	if (m==0) {
	
		return (rp*cos(phi) - r)*pow(r*r + rs - 2*r*rp*cos(phi),-1.5);
	}
	else if (m==1) {
		return cos(phi)/sqrt(r*r + rs - 2*r*rp*cos(phi));
	}
	else {
		return 0;
	}
}
