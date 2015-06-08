#include "eigen.h"
#include <time.h>
#include <unistd.h>

void print_time(double t);

int main(int argc, char *argv[]) {
	int i;

	clock_t tic, toc;
	tic = clock();
	printf("Reading arguments...\n");

	for(i=1;i<argc;i++) {
		printf("%s ",argv[i]);
	}
	printf("\n");

	if ((argc < 15) && (argc != 3)) {
		printf("\n\nToo Few Arguments!\n\n");
		return -1;
	}

	read_arguments(argc,argv);






#if defined(ISOTHERMAL) || defined(BAROTROPIC)
	printf("\nUsing the Parameters:\n \
		\tNr=%d\n \
		\tInner radius = %lg\n \
		\tOuter radius = %lg\n \
		\tLog spacing = %.3e\n \
		\tDisk Mass = %lg\n \
		\tSigma Index = %lg\n \
		\tFlare Index = %lg\n \
		\tAlpha Viscosity = %lg\n \
		\tAdiabatic Index = 1\n \
		\tBeta Cooling = inf\n",
		N,ri,ro,dlr,Mdisk,sigma_index,flare_index,alpha_s);
#else
	printf("\nUsing the Parameters:\n \
		\tNr=%d\n \
		\tInner radius = %lg\n \
		\tOuter radius = %lg\n \
		\tLog spacing = %.3e\n \
		\tDisk Mass = %lg\n \
		\tSigma Index = %lg\n \
		\tFlare Index = %lg\n \
		\tAlpha Viscosity = %lg\n \
		\tAdiabatic Index = %lg\n \
		\tBeta Cooling = %lg\n",
		N,ri,ro,dlr,Mdisk,sigma_index,flare_index,alpha_s,adi_gam,beta_cool);
#endif

#ifdef OPENMP
	omp_set_num_threads(nprocs);
	printf("\t\t\tOpenMP threads = %d\n", nprocs);
#endif

	double complex *mat = (double complex *)malloc(sizeof(double complex)*nrows*ncols);
	double complex *evals = (double complex *)malloc(sizeof(double complex)*nrows);
	double complex *evecs = (double complex *)malloc(sizeof(double complex)*nrows*ncols);
	double complex *bcmat = (double complex *)malloc(sizeof(double complex)*nrows*ncols);




  	printf("Allocating Arrays...\n");

  	alloc_globals();




  	printf("Initializing Variables...\n");


  	init(ri,ro);


  	printf("Outputting Variables...\n");
  	output_globals();

  	printf("Outputting Derivative Matrices...\n");
  	output_derivatives();
#ifndef READKERNEL
	printf("Outputting Kernel...\n");
  	output_kernel();
#endif

#ifdef PLANETS
	printf("Calculating Matrices for %d Planets...\n",NP);
	calc_planet_matrices();
#endif

  	printf("Populating Matrix...\n");


   calc_matrices(mat,bcmat);


	printf("Outputting Coefficients...\n");
  	output_coefficients();
  	printf("Outputting Matrix...\n");
  	output_matrix(mat,bcmat);

  	printf("Solving For Eigenvalues and Eigenvectors...\n");
  	reigenvalues(mat,bcmat,evals,evecs,nrows);


  	printf("Outputting Results...\n");
  	output(evals,evecs);


  	printf("Freeing Arrays...\n");

  	free_globals();
#ifdef PLANETS
	free_planets();
#endif
	free(mat); free(evecs); free(evals); free(bcmat);

	toc = clock();
	print_time( (double)(toc - tic) / CLOCKS_PER_SEC );
  return 0;
}
void print_time(double t) {
	int hr, min;
	hr = (int)floor(t/(60.*60.));
	t -= hr*60*60;
	min = (int)floor(t/60);
	t -= min*60;


	if (hr==0) {
		if (min == 0) {
			printf("Total Runtime:\t%.3lgs\n",t);

		}
		else {
			printf("Total Runtime:\t%dm%.3lgs\n",min,t);
		}
	}
	else {
		printf("Total Runtime:\t%dh%dm%.3lgs\n",hr,min,t);
	}
	return;
}
