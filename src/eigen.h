#ifndef HEADER_H
#define HEADER_H

#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>



#ifdef OPENMP
#include <omp.h>
#endif


int N, NP,nrows,ncols;

double Mdisk, eps, h0, dlr, flare_index, sigma_index, temp_index,sigma0, tol;

double *weights, *kernel, *work;
double complex *cwork;


double *r, *lr , *scaleH,*omega, *dldom, *d2ldom, *sigma, *dlds, *d2lds, *pres, *dldpres, *d2ldpres;
double *temp , *dldtemp, *d2ldtemp, *c2, *dldc2, *d2ldc2, *omega_prec;

double complex *coeffs_A , *coeffs_B, *coeffs_C;
double *D, *D2, *Identity;

double alpha_s, alpha_b;

double adi_gam, beta_cool;


#ifdef PLANETS
typedef struct Planet {

	double mass, position, hill,wp;
	double complex *pot0, *pot1;
	double complex e;
	int index,InteriorPlanet, ExteriorPlanet;

} Planet;

Planet *Planets;

#endif


void alloc_globals(void);
void free_globals(void);
int init(double ri,double ro);

void init_derivatives(void);


int calc_matrices(double complex *mat, double complex *bcmat) ;
void set_bc(double complex *mat, double complex *bcmat);

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


void output_globals(void);
void output_kernel(void);
void output_matrix(double complex *mat, double complex *bcmat);
void output(double complex *evals, double complex *evecs);
void output_derivatives(void);
void output_coefficients(void);

void cmatvec(double  complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nB) ;
void matvec(double  *A, double*B, double *C,
					double alpha, double beta, int nB) ;
void cmatmat(double complex *A, double complex *B, double complex *C,
					double complex alpha, double complex beta, int nA) ;
void matmat(double  *A, double *B, double *C,
					double alpha, double beta, int nA) ;
void reigenvalues(double complex *A, double complex *Q, double complex *evals, double complex *evecs, int nA) ;


#ifdef SELFGRAVITY
void compute_kernels(void);
void add_sg(double complex *mat, double complex *bcmat);
void calc_omega_prec_grav(void);
double sg_integrand(double r1, double rp1, double eps1);
#endif
int analytic_potential(void);
double omega_prec_grav_analytic(double x);
#endif

#ifdef PLANETS
void add_planets(double complex *mat, double complex *bcmat);
void init_planets(void);
void calc_planet_matrices(void);
void free_planets(void);
void output_planet_summary(void);
#endif
