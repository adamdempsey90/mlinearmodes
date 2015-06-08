#include "eigen.h"

void normalize_evectors(double complex *evecs);

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
	int LDB = nA;
	int LDVL = nA;
	int LDVR = nA;
	int LWORK = 2*nA;

	double *RWORK = (double *)malloc(sizeof(double)*8*nA);
	double complex *CWORK = (double complex *)malloc(sizeof(double complex)*2*nA);

	double complex *tA = (double complex *)malloc(sizeof(double complex)*nA*nA);
	double complex *tQ = (double complex *)malloc(sizeof(double complex)*nA*nA);

	double complex *evals_alpha = (double complex *)malloc(sizeof(double complex)*nA);
	double complex *evals_beta = (double complex *)malloc(sizeof(double complex)*nA);

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++) {
	  		tA[i+nA*j] = A[j+nA*i];
	  		tQ[i+nA*j] = Q[j+nA*i];
		}
	}


#ifdef NOPRESSURE
	zgeev_( &JOBVL, &JOBVR, &nA, tA, &LDA, evals, tQ, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );
#else
	zggev_( &JOBVL, &JOBVR, &nA, tA, &LDA, tQ, &LDB, evals_alpha,evals_beta, NULL, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );

	for(i=0;i<nA;i++) {
		if (cabs(evals_beta[i]) != 0) {
			evals[i] = evals_alpha[i]/evals_beta[i];
		}
	}
#endif

#if defined(NORMALIZE_INT) || defined(NORMALIZE_MAX)
	printf("Normalzing eigenvectors...\n");
	normalize_evectors(evecs);
#endif


	free(tA); free(tQ);
	free(RWORK); free(CWORK);
	free(evals_alpha);
	free(evals_beta);
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

void cmatmat(double complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nA)
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
		for(j=0;j<nA;j++) cwork[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = cwork[i+N*j];
	}

	zgemm_(&TRANSA, &TRANSB, &m,&n,&k,&alpha,A,&LDA,B,&LDB,&beta,C,&LDC);


	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)  cwork[i+N*j] = C[j + nA*i];
	}

	for(i=0;i<nA;i++) {
		for(j=0;j<nA;j++)	C[i + nA*j] = cwork[i+N*j];
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

void cmatvec(double  complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nB)
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



	zgemv_(&TRANS, &m,&n,&alpha,A,&LDA,B,&INCX,&beta,C,&INCY);

	return;

}



void normalize_evectors(double complex *evecs) {
/* Normalize the eigenvectors */
/* Calculate the factor to normalize the disk eccentricity.
	 Each planet eccentricity will then be normalized by the same factor.
*/
	int i,j,indx;
	double norm;


#ifdef OPENMP
#pragma omp parallel private(i,j,norm,indx) shared(evecs,nrows,ncols)
#pragma omp for schedule(static)
#endif
	for(i=0;i<nrows;i++) {


			norm = 0;
#ifdef NORMALIZE_INT
			for(j=0;j<N;j++) {
				indx = j + ncols*i;

				norm += weights[j]*conj(evecs[indx])*evecs[indx];

			}
			norm = sqrt(norm);
#else
#ifdef NORMALIZE_MAX
			for(j=0;j<N;j++) {
				indx = j+ncols*i;
//				printf("%lg\t%lg\t%lg",norm,abs(evecs[indx]),fmax(norm,abs(evecs[indx])));
				norm = fmax(norm,abs(evecs[indx]));
			}
			if (norm == 0) norm = 1;


#endif
#endif
			for(j=0;j<ncols;j++) {
					indx = j + ncols*i;
					evecs[indx] /= norm;
			}
	}

	return;


}
