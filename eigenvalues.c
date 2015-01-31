#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define OPENMP
#ifdef OPENMP	
#include <omp.h>
#endif

#define NPHI 5000

int N;

double poly_n;
double Mdisk, eps, h0, dlr,rout; 

double *weights, *kernel,*kernel0, *H, *HL, *work;


double *c2, *sigma, *scaleH,  *r, *lr, *dlds, *dldc2, *lsig, *lc2, *d2lds;
double *omega,*kappa,*omega2,*kappa2;
double *dphi0dr;

double *D, *D2, *KD,*DKD;


void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA);
void alloc_globals(void);
void free_globals(void);
int init(double ri);
double Dij(int i, int j);
double D2ij(int i, int j);
void calc_weights(void) ;
void init_derivatives(void);
int calc_matrices(double complex *mat, double complex *bcmat);
void calc_coefficients(int i, double *A, double *B, double *C, double *G);
double K0ij(int i, int j);
double Kij(int i, int j);
double kernel_integrand(int m, double r, double rp, double rs, double phi);

void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat);
void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat);

void matmat(double  *A, double *B, double *C, 
					double alpha, double beta, int nA); 

void matvec(double  *A, double *B, double *C, 
					double alpha, double beta, int nA); 


void output(double complex *evals, double complex *evecs);

void output_globals(void);
void output_matrix(double complex *mat, double complex *bcmat);
void output_kernel(void);
					
int main(int argc, char *argv[]) {
	int i,j;
	double ri, ro;
	
	printf("Reading arguments...\n");
	if (argc < 9) {
		printf("\n\nToo Few Arguments!\n\n");
		return -1;
	}
	N = atoi(argv[1]);
	
	ri = atof(argv[2]);
	
	ro = atof(argv[3]);
	
	Mdisk = atof(argv[4]);
	eps = atof(argv[5]);
	h0 = atof(argv[6]);
	poly_n = atof(argv[7]);
	
	
	dlr = (log(ro) - log(ri))/((float) N);
	
	rout = 100;
	
	printf("Using the Parameters:\n \
		\tNr=%d\n \
		\tInner radius = %lg\n \
		\tOuter radius = %lg\n \
		\tLog spacing = %.3e\n \
		\tDisk Mass = %lg\n \
		\tPolytropic Index = %lg\n",
		N,ri,ro,dlr,Mdisk,poly_n);
	
#ifdef OPENMP
	omp_set_num_threads(atoi(argv[8]));
	printf("\t\tOpenMP threads = %d\n", atoi(argv[8]));
#endif	
	
	double complex *mat = (double complex *)malloc(sizeof(double complex)*N*N);	
	double complex *evals = (double complex *)malloc(sizeof(double complex)*N);	
	double complex *evecs = (double complex *)malloc(sizeof(double complex)*N*N);	
	double complex *bcmat = (double complex *)malloc(sizeof(double complex)*N*N);	  
  
  
  	printf("Allocating Arrays...\n");
  	alloc_globals();
  	
  	printf("Initializing Derivative Matrices...\n");
  	init_derivatives();
  	calc_weights();
  	printf("Initializing Variables...\n");
  	int nanflag = init(ri);
  	
  	if (nanflag == -1) {
  		printf("Aborting...\n");
  		free_globals();
		free(mat); free(evecs); free(evals); free(bcmat);
		return nanflag;
	}
  	
  	printf("Outputting Variables...\n");
  	output_globals();
  	output_kernel();
  	printf("Populating Matrix...\n");
  	nanflag = calc_matrices(mat,bcmat);
  	if (nanflag == -1) {
  		printf("Aborting...\n");
  		free_globals();
		free(mat); free(evecs); free(evals); free(bcmat);
		return nanflag;
	}
  	printf("Outputting Matrix...\n");
  	output_matrix(mat,bcmat);
  	printf("Solving For Eigenvalues and Eigenvectors...\n");
  	reigenvalues(mat,bcmat,evals,evecs,N);
  	
  	printf("Outputting Results...\n");
  	output(evals,evecs);

  	
  	printf("Freeing Arrays...\n");
  	
  	free_globals();
	free(mat); free(evecs); free(evals); free(bcmat);
	
  return 0;
}

void alloc_globals(void) {
	
	D = (double *)malloc(sizeof(double)*N*N);
	D2 = (double *)malloc(sizeof(double)*N*N);
	kernel = (double *)malloc(sizeof(double)*N*N);
	kernel0 = (double *)malloc(sizeof(double)*N*N);
	weights = (double *)malloc(sizeof(double)*N);
	KD = (double *)malloc(sizeof(double)*N*N);
	DKD = (double *)malloc(sizeof(double)*N*N);
	H = (double *)malloc(sizeof(double)*N*N);
	HL = (double *)malloc(sizeof(double)*N*N);
	work = (double *)malloc(sizeof(double)*N*N);
	
	r = (double *)malloc(sizeof(double)*N);
	lr = (double *)malloc(sizeof(double)*N);
	c2 = (double *)malloc(sizeof(double)*N);
	sigma = (double *)malloc(sizeof(double)*N);
	omega = (double *)malloc(sizeof(double)*N);
	omega2 = (double *)malloc(sizeof(double)*N);
	kappa = (double *)malloc(sizeof(double)*N);
	kappa2 = (double *)malloc(sizeof(double)*N);
	scaleH = (double *)malloc(sizeof(double)*N);
	dlds = (double *)malloc(sizeof(double)*N);
	dldc2 = (double *)malloc(sizeof(double)*N);
	lsig = (double *)malloc(sizeof(double)*N);
	lc2 = (double *)malloc(sizeof(double)*N);
	dphi0dr=(double *)malloc(sizeof(double)*N);
	d2lds = (double *)malloc(sizeof(double)*N);
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
	free(work);
	
	free(r);
	free(lr);
	free(c2);
	free(sigma);
	free(omega);
	free(kappa);
	free(omega2);
	free(kappa2);
	free(scaleH);
	free(dlds);
	free(dldc2);
	free(lsig);
	free(lc2);
	free(dphi0dr);
	free(d2lds);
	
	return;

}

int init(double ri) {

	int i,j,indx;
	double sigfac;
	double phifac = 0;
	
	
	
	for(i=0;i<N;i++) {
		lr[i] = log(ri) + i*dlr;
		r[i] = exp(lr[i]);
		
		scaleH[i] = h0*r[i];
		
		c2[i] = scaleH[i] * scaleH[i] / (r[i]*r[i]*r[i])*(1 - pow(r[i],-10))*(1-pow(r[i]/rout,10));
		sigma[i] = pow(c2[i],poly_n);
		
		omega[i] = pow(r[i],-1.5);
		omega2[i] = omega[i]*omega[i];
		lc2[i] = log(c2[i]);
		lsig[i] = log(sigma[i]);
		
		dlds[i] = 0;
		dldc2[i] = 0;
		d2lds[i] = 0;
		
		if (isnan(c2[i]) != 0) {
			printf("\n\n Detected NaN in c2 at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
		
		if (isnan(lc2[i]) != 0) {
			printf("\n\n Detected NaN in lc2 at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
		
		if (isnan(lsig[i]) != 0) {
			printf("\n\n Detected NaN in lsig at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
		

	}
	
	
	printf("Calculating Kernels\n");

#ifdef OPENMP
#pragma omp parallel private(indx,i,j) shared(kernel,kernel0,N)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		
		kernel[indx] = Kij(i,j);
		kernel0[indx] = K0ij(i,j);
	}
	printf("Finished part 1\n");
	
	sigfac=0;
	for(i=0;i<N;i++) {
		sigfac += weights[i]*sigma[i]*r[i]*r[i];	
		dphi0dr[i] = 0;
		for(j=0;j<N;j++) {
			dphi0dr[i] -= weights[j]*kernel0[j+N*i]*sigma[j];
		}
		if (isnan(dphi0dr[i]) != 0) {
			printf("\n\n Detected NaN in dphi0dr at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
	}
	
	sigfac = Mdisk/(2*M_PI*sigfac);
	
	for(i=0;i<N;i++) {
		sigma[i] *= sigfac;
		dphi0dr[i] *= sigfac;
		if (isnan(sigma[i]) != 0) {
			printf("\n\n Detected NaN in sigma at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
	}	
	
	
	printf("Finished part 2\n");
	matvec(D,lsig,dlds,1,0,N);
	matvec(D,lc2,dldc2,1,0,N);
	matvec(D2,lsig,d2lds,1,0,N);
	
	
	
	for(i=0;i<N;i++) {
		omega2[i] += (c2[i]*dlds[i] + dphi0dr[i])/(r[i]*r[i]);
		kappa2[i] = 4*omega2[i];
	}
	
	matvec(D,omega2,kappa2,1,1,N);
	
	for(i=0;i<N;i++) {
		kappa[i] = sqrt(kappa2[i]);
		omega[i] = sqrt(omega2[i]);
		if (isnan(kappa[i]) != 0) {
			printf("\n\n Detected NaN in kappa at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
		if (isnan(omega[i]) != 0) {
			printf("\n\n Detected NaN in omega at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}

	}
	
	
	
	return 0;
}

double Dij(int i, int j) {
	
	if (i==0) {
		if (j==0) {
			return -1.0;
		}
		else if (j==1) {
			return 1.0;
		}
		else {
			return 0;
		}
	}
	else if (i==(N-1)) {
		if (j==(N-2)) {
			return -1.0;
		}
		else if (j==(N-1)) {
			return 1.0;
		}
		else {
			return 0;
		}
	}
	else {
		if (i == (j+1)) {
			return -.5;
		}
		else if (i==(j-1)) {
			return .5;
		}
		else {
			return 0;
		}
	}
		
}

double D2ij(int i, int j) {
	
	if (i==0) {
		if (j==0) {
			return 1.0;
		}
		else if (j==1) {
			return -2.0;
		}
		else if (j==2) {
		
			return 1.0;
		}
		else {
			return 0;
		}
	}
	else if (i==(N-1)) {
		
		if (j==(N-3)) {
			return 1.0;
		
		}
		else if (j==(N-2)) {
			return -2.0;
		}
		else if (j==(N-1)) {
			return 1.0;
		}
		else {
			return 0;
		}
	}
	else {
		if (i == (j+1)) {
			return 1.0;
		}
		else if (i==(j-1)) {
			return 1.0;
		}
		else if (i==j) {
			
			return -2.0;
		
		}
		else {
			return 0;
		}
	}
		
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
void init_derivatives(void) {
	int i,j;
	double dfac = 1.0/dlr;
	double d2fac = dfac*dfac;
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			D[j + N*i] = Dij(i,j)*dfac;
			D2[j + N*i] = D2ij(i,j)*d2fac;
		
		}
	}
	return;

}
int calc_matrices(double complex *mat, double complex *bcmat) {
	int i,j,indx;
	
	double A, B, C, G;
/* Compute the matrix including all of its component matrices */

	
// #ifdef OPENMP
// #pragma omp parallel 
// #pragma omp for private(i,A,B,C,indx,j,
// #endif	
	for(i=0;i<N;i++ ) {
		calc_coefficients(i,&A,&B,&C,&G);

		for(j=0 ; j<N ; j++) {
			indx = j + N*i;
			
			if (i==j) {
				mat[indx] = A;
				bcmat[indx] = 1;
				KD[indx] = -sigma[j]* ( dlds[j] + D[indx]);
				HL[indx] = G*(2*r[i] + D[indx]);
			}
			else {
				mat[indx] = 0;
				bcmat[indx] = 0;
				KD[indx] = -D[indx]*sigma[j];
				HL[indx] = G*D[indx];
			}
			
			mat[indx] += B*D[indx] + C*D2[indx];
			
			H[indx] = 0;
			DKD[indx] = 0;
			
		
		}
	}
	printf("Done 1\n");
			
			
			
//	matmat(kernel,D,KD,1,0,N);
	printf("Done 2\n");
//	matmat(D,KD,DKD,1,0,N);
	matmat(kernel,KD,H,1,0,N);
	printf("Done 3\n");
	matmat(HL,H,DKD,1,0,N);
	
	
	for(i=0;i<N*N;i++) mat[i] += DKD[i];
	
	printf("Done 4\n");
	lagrangian_pressure_bc_inner(mat, bcmat);
	printf("Done 5\n");
	lagrangian_pressure_bc_outer(mat, bcmat);
	printf("Done 6\n");

	return 0;
}


void calc_coefficients(int i, double *A, double *B, double *C, double *G) {

	
	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) * ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega[i] - kappa[i];
	*B = (2 + dldc2[i] + dlds[i]) * (*C);
	*G = -1.0/(2*omega[i]*r[i]*r[i]*r[i]*r[i]);
	return; 

}


double K0ij(int i, int j) {
//	printf("start\n");
	int indx;
	double dp = M_PI/NPHI;

	double result = 0;
	double rs2 = r[j]*r[j] + eps*eps*scaleH[j]*scaleH[j];
	
/* phi_i = i*dp */
	result = (3./8)*(kernel_integrand(0,r[i],r[j],rs2,0) + kernel_integrand(0,r[i],r[j],rs2,dp*NPHI));
	result += (7./6)*(kernel_integrand(0,r[i],r[j],rs2,dp) + kernel_integrand(0,r[i],r[j],rs2,dp*(NPHI-1)));
	result += (23./24)*(kernel_integrand(0,r[i],r[j],rs2,dp*2) + kernel_integrand(0,r[i],r[j],rs2,dp*(NPHI-2)));
	
	for(indx=3;indx<NPHI-2;indx++) {
		result += kernel_integrand(0,r[i],r[j],rs2,dp*indx);
	}
	
	result *= dp*2*r[j]*r[j];
	
//	printf("%d\t%d\t%lg\n",i,j,result);
	return result;


}
double Kij(int i, int j) {
	int indx;
	double dp = M_PI/NPHI;

	double result = 0;
	double rs2 = r[j]*r[j] + eps*eps*scaleH[j]*scaleH[j];
	
/* phi_i = i*dp */
	result = (3./8)*(kernel_integrand(1,r[i],r[j],rs2,0) + kernel_integrand(1,r[i],r[j],rs2,dp*NPHI));
	result += (7./6)*(kernel_integrand(1,r[i],r[j],rs2,dp) + kernel_integrand(1,r[i],r[j],rs2,dp*(NPHI-1)));
	result += (23./24)*(kernel_integrand(1,r[i],r[j],rs2,dp*2) + kernel_integrand(1,r[i],r[j],rs2,dp*(NPHI-2)));
	
	for(indx=3;indx<NPHI-2;indx++) {
		result += kernel_integrand(1,r[i],r[j],rs2,dp*indx);
	}
	
	
	
	result *= dp*2;
//	result = weights[j]*r[i]*r[i]*(r[j]*r[j]*result - M_PI*r[i]);
	result = weights[j] * (r[j]*r[j]*result - M_PI*r[i]);	
	return result;


}

double kernel_integrand(int m, double r, double rp, double rs2, double phi) {
	if (m==0) {
	
		return (rp*cos(phi) - r)*pow(r*r + rs2 - 2*r*rp*cos(phi),-1.5);
	}
	else if (m==1) {
		return cos(phi)/sqrt(r*r + rs2 - 2*r*rp*cos(phi));
	}
	else {
		return 0;
	}
}

void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {
	int j,indx;
	for(j=0;j<N;j++) {
		indx = j;
		mat[indx] = D[indx];
		bcmat[indx] = 0;		
	}
	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {
	int j,indx;
	
	for(j=0;j<N;j++) {
		indx= j + N*(N-1);
		mat[indx] = D[indx];
		bcmat[indx] = 0;
	}
	return;
}


void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA) 
{
/* Computes the eigenvalues and right eigenvectors of the complex matrix A which is NxN.
	A . evecs = evals Q . evecs 
	This is essentially a wrapper for the ZGEEV LAPACK routine.
	
	INPUTS:
		The matrices M, Q, and evecs and the vector evals which are all overwritten
		
	OUTPUTS:
		The eigenvalues are stored in the evals array.
		The eigenvectors are stored in the ROWS of the evecs matrix	
*/
	int i,j;
	char JOBVL = 'N';
	char JOBVR = 'V';
	int INFO; 
	int LDA = nA;
	int LDVL = nA;
	int LDVR = nA;
	int LWORK = 2*nA;
	
	double *RWORK = (double *)malloc(sizeof(double)*2*nA);
	double complex *CWORK = (double complex *)malloc(sizeof(double complex)*2*nA);
	
	double complex *tA = (double complex *)malloc(sizeof(double complex)*nA*nA);
	double complex *tQ = (double complex *)malloc(sizeof(double complex)*nA*nA);
	

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) {
	  		tA[i+nA*j] = A[j+nA*i];
	  		tQ[i+nA*j] = Q[j+nA*i];
		}
	}


	zgeev_( &JOBVL, &JOBVR, &nA, tA, &LDA, evals, tQ, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );

	free(tA); free(tQ);
	free(RWORK); free(CWORK);
 	return;
}

void matmat(double  *A, double *B, double *C, 
					double alpha, double beta, int nA) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C.
	A,B, and C are all matrices.
	This is essenitally a wrapper for the ZGEMM BLAS routine 
*/
	int i,j;
	char TRANSA = 't';
	char TRANSB = 't';
	int m = nA;
	int n = nA;
	int k = nA;
	int LDA = nA;
	int LDB = nA;
	int LDC = nA;
		 
	
	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) work[i+N*j] = C[j + nA*i];
	}
	
	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+N*j];
	}
	
	dgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  work[i+N*j] = C[j + nA*i];
	}
	
	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = work[i+N*j];
	}
	return;

}


void matvec(double  *A, double*B, double *C, 
					double alpha, double beta, int nB) 
{
/* Performs \alpha * A.B + \beta * C and stores the output in C. 
	A is a matrix, B and C are vectors.
	This is essenitally a wrapper for the ZGEMV BLAS routine 
*/

	char TRANS = 't';
	int m = nB;
	int n = nB;
	int LDA = nB;
	int INCX = 1;
	int INCY = 1;
		 
	

	dgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}

void output_globals(void) {
	int i;
	
	FILE *f = fopen("globals.dat","w");
	
	
	fprintf(f,"# lr \t r \t omega \t c2 \t sigma \t H/r \t soft \t dldc2 \t dlds \t kappa\n");
	
	
	for(i=0;i<N;i++) {
		fprintf(f,"%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\n",
			lr[i],
			r[i],
			omega[i],
			c2[i],
			sigma[i],
			scaleH[i]/r[i],
			eps*scaleH[i],
			dldc2[i],
			dlds[i],
			kappa[i],
			d2lds[i]);
	}
			
			
		
			
	fclose(f);
	return;


}


void output_kernel(void) {
	int i,j;
	FILE *f = fopen("kernel.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%lg\t",kernel[j+i*N]);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	
	f = fopen("kernel0.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%lg\t",kernel0[j+i*N]);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);

	return;
}


void output_matrix(double complex *mat, double complex *bcmat) {
	int i,j;
	FILE *f = fopen("matrix.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%lg\t%lg\t",creal(mat[j+i*N]),cimag(mat[j+i*N]));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	
	f = fopen("bcmatrix.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%lg\t%lg\t",creal(bcmat[j+i*N]),cimag(bcmat[j+i*N]));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);

	return;
}

void output(double complex *evals, double complex *evecs) {
	int i,j;
	double complex evecij;
	FILE *f = fopen("eigen.dat","w");
	
	fprintf(f,"# lambda \t evecs along row\n");
	for(i=0;i<N;i++) {
		fprintf(f,"%.12lg\t%.12lg\t",creal(evals[i]),cimag(evals[i]));
		for(j=0;j<N;j++) {
			evecij = evecs[j+N*i];
// 			evecij /= sigma[i];
// 			evecij /= evecs[N*i];
// 			evecij *= .1;
// 			
			fprintf(f,"%.12lg\t%.12lg\t",creal(evecij),cimag(evecij));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	return;
}



