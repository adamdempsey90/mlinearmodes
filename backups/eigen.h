#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
#include <time.h>
#include <unistd.h>


#ifdef OPENMP	
#include <omp.h>
#endif

#define RDECAY 5


//#define INDIRECT


//#define FOURTHORDER


int N;

double poly_n;
double Mdisk, eps, h0, dlr,rout, flare_index, sigma_index, sigma0, tol; 

double *weights,*kernel0, *kernel02, *kernel,*work;
double complex *H, *HL, *KD,*DKD, *cwork;


double *c2, *sigma, *scaleH,  *r, *lr, *dlds, *dldc2, *lsig, *lc2, *d2lds, *lom, *dldom, *d2dom;
double *omega,*omega2,*kappa2, *omegap2, *kappap2, *omegag2, *kappag2, *omega_prec;
double complex *kappa;
double *dphi0dr;

double *D, *D2;

double *nu, *dldnu, *lnu;
double alpha_s, alpha_b;

double *pres, *lpres, *temp, *d2ldpres, *dldpres, *ltemp, *dldtemp, *d2ldtemp;
double adi_gam, beta_cool;




void print_time(double t);
void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA);
void alloc_globals(void);
void free_globals(void);
int init(double ri, double ro);
double Dij(int i, int j);
double D2ij(int i, int j);
void calc_weights(void) ;
void init_derivatives(void);
int calc_matrices(double complex *mat, double complex *bcmat);
void calc_coefficients(int i, double complex *A, double complex *B, double complex *C, double complex *G);
double K0ij(int i, int j);
double Kij(int i, int j);
double kernel_integrand(int m, double r, double rp, double rs, double phi);

void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat);
void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat);

void calc_sigmap(double complex *sigmap, double complex *evecs);

void matmat(double  *A, double *B, double *C, 
					double alpha, double beta, int nA); 

void matvec(double  *A, double *B, double *C, 
					double alpha, double beta, int nA); 


void output(double complex *evals, double complex *evecs, double complex *sigmap);

void output_globals(void);
void output_matrix(double complex *mat, double complex *bcmat);
void output_kernel(void);
void output_derivatives(void);
void cmatvec(double  complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta, int nB);
					
void cmatmat(double complex *A, double complex *B, double complex *C, 
					double complex alpha, double complex beta, int nA);					

double sigma_profile(double rval, double h, double rc,double beta);

void read_kernel(void);
void compute_kernels(void);
void add_edge_sg(double complex *mat);
#ifdef MLIN
double bump_function(double rval);
#endif

#ifdef TESTFUNCTION
double complex test_function(double rval);
void fill_mat(double complex *mat, double complex *bcmat);
#endif

double sigma_function(double rval);
double dlogsigma_dlogr(double rval);
double dlogsigma2_dlogr2(double rval);

void calc_epicyclic(void);
void output_omega_corrections(double *omegap2, double *omegag2, double *kappap2, double *kappag2) ;


double sigma_func(double x);
double temp_func(double x);
double omk_func(double x); 
double scaleH_func(double x); 
void calc_omega_prec_grav(void);
void calc_omega_prec_pres(void);

