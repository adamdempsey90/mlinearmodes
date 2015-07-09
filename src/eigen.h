#ifndef HEADER_H
#define HEADER_H

#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


/* NF can be defined in defines.h based on EOS, either 3 or 4
* If it's not then define it here.
*/
#ifndef NF
#if defined(BAROTROPIC) || defined(ISOTHERMAL)
#define NF 3
#else
#define NF 4
#endif
#endif

#define MAXSTRLEN 300


int N, NP,nrows,ncols, nprocs;

double Mdisk, ri, ro,eps, h0, dlr, flare_index, sigma_index, temp_index,sigma0, tol;
int mval,mval2;

char outputname[MAXSTRLEN];


double *weights, *kernel, *work;
double complex *cwork;


double *r, *lr , *scaleH,*omega, *dldom, *d2ldom, *sigma, *dlds, *d2lds, *pres, *dldpres, *d2ldpres, *aspect_ratio;
double *temp , *dldtemp, *d2ldtemp, *c2, *dldc2, *d2ldc2, *kappa2;


double complex *coeffs_A , *coeffs_B, *coeffs_C;
double *D, *D2, *Identity;

double alpha_s, alpha_b;

double adi_gam, beta_cool;



#ifdef PLANETS
typedef struct Planet {

	double mass;
	double position;
	double hill;
	double wp;
	double complex *pot0;
	double complex *pot1;
	int index;
	int InteriorPlanet;
	int ExteriorPlanet;

} Planet;

Planet *Planets;

#endif

/* Global Prototypes */


/* error.c */
void gsl_integration_error_handler(int status, char *reason);

/* readinputs.c */
void read_arguments(int argc, char *argv[]);

/* alloc.c */
void alloc_globals(void);
void free_globals(void);

/* init.c */
int init(double ri,double ro);

/* derivatives.c */
void init_derivatives(void);

/* algo.c */
int calc_matrices(double complex *mat, double complex *bcmat) ;

/* boundary.c */
void set_bc(double complex *mat, double complex *bcmat);

/* profiles.c */
double sigma_func(double x);
double dlogsigma_func(double x);
double d2logsigma_func(double x);
double temp_func(double x);
double dlogtemp_func(double x);
double d2logtemp_func(double x);
double omk_func(double x);
double dlogomk_func(double x);
double d2logomk_func(double x);
double scaleH_func(double x);

/* output.c */
void output_globals(void);
void output_kernel(void);
void output_matrix(double complex *mat, double complex *bcmat);
void output(double complex *evals, double complex *evecs);
void output_derivatives(void);
void output_coefficients(void);

/* matrix.c */
int getindex4(int i, int j, int k, int l, int size_small, int size_large);
int getindex3(int i, int k, int l,int size);
int getindex2(int k, int l, int size);
void cmatvec(double  complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nB) ;
void matvec(double  *A, double*B, double *C,
					double alpha, double beta, int nB) ;
void cmatmat(double complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nA) ;
void matmat(double  *A, double *B, double *C,
					double alpha, double beta, int nA) ;
void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA) ;

/* selfgravity.c */
#ifdef SELFGRAVITY
void add_sg(double complex *mat, double complex *bcmat);
void calc_omega_sg(double *omg2, double *kapg2);
#endif

/* hdf5.c */
#ifdef HDF5_OUTPUT
void output_hdf5_file(double complex *mat,double complex *bcmat,double complex *evecs,double complex *evals);
#endif
#endif
