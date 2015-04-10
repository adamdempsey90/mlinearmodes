#include "eigen.h"

int main(int argc, char *argv[]) {
	int i,j;
	double ri, ro;
	clock_t tic, toc;
	tic = clock();
	printf("Reading arguments...\n");
	for(i=1;i<argc;i++) {
		printf("%s\t",argv[i]);
	}

#ifdef ISOTHERMAL	
	if (argc < 13) {
		printf("\n\nToo Few Arguments!\n\n");
		return -1;
	}
#else 
	if (argc < 15) {
		printf("\n\nToo Few Arguments!\n\n");
		return -1;
	}
#endif

#ifdef TESTFUNCTION
	printf("\n\n\n\n\n\n RUNNING IN TEST FUNCTION MODE \n\n\n\n\n\n");
#endif

	N = atoi(argv[1]);
	
	ri = atof(argv[2]);
	
	ro = atof(argv[3]);

#ifdef INPUTMASS	
	Mdisk = atof(argv[4]);
	sigma0 = 1;
#else
	sigma0 = atof(argv[4]);
	Mdisk = 1;
#endif

	eps = atof(argv[5]);
	h0 = atof(argv[6]);
	sigma_index = atof(argv[7]);
	flare_index = atof(argv[8]);

#ifdef VISCOSITY
	alpha_s = atof(argv[9]);
	alpha_b = atof(argv[10]);
#else
	alpha_s = 0;
	alpha_b = 0;
#endif

#if  defined(COOLING) || defined(ADIABATIC) 
	adi_gam = atof(argv[12]);
	beta_cool = atof(argv[13]);
#else
	adi_gam = 1;
#endif

#ifdef ADIABATIC 
	beta_cool = 0;
#endif

	tol = atof(argv[14]);
	dlr = (log(ro) - log(ri))/((float) N);
	
	rout = 100;



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
	omp_set_num_threads(atoi(argv[11]));
	printf("\t\t\tOpenMP threads = %d\n", atoi(argv[11]));
#endif	
	
	double complex *mat = (double complex *)malloc(sizeof(double complex)*N*N);	
	double complex *evals = (double complex *)malloc(sizeof(double complex)*N);	
	double complex *evecs = (double complex *)malloc(sizeof(double complex)*N*N);	
	double complex *bcmat = (double complex *)malloc(sizeof(double complex)*N*N);	
	
	double complex *sigmap = (double complex *)malloc(sizeof(double complex)*N*N);
	  
  
  
  	printf("Allocating Arrays...\n");
  	alloc_globals();
  	
  	printf("Initializing Derivative Matrices...\n");
  	init_derivatives();
  	output_derivatives();
  	calc_weights();
  	 	
  	printf("Initializing Variables...\n");
  	int nanflag = init(ri,ro);


  	if (nanflag == -1) {
  		printf("Aborting...\n");
  		free_globals();
		free(mat); free(evecs); free(evals); free(bcmat);
		return nanflag;
	}
  	
  	printf("Outputting Variables...\n");
  	output_globals();
#ifndef READKERNEL
  	output_kernel();
#endif
  	printf("Populating Matrix...\n");

#ifdef TESTFUNCTION 
	fill_mat(mat,bcmat);
#else
  	nanflag = calc_matrices(mat,bcmat);
  	if (nanflag == -1) {
  		printf("Aborting...\n");
  		free_globals();
		free(mat); free(evecs); free(evals); free(bcmat);
		return nanflag;
	}
#endif
  	printf("Outputting Matrix...\n");
  	output_matrix(mat,bcmat);
  	printf("Solving For Eigenvalues and Eigenvectors...\n");
  	reigenvalues(mat,bcmat,evals,evecs,N);

#ifndef TESTFUNCTION  	
  	calc_sigmap(sigmap, evecs);
#endif
  	printf("Outputting Results...\n");
  	output(evals,evecs,sigmap);

  	
  	printf("Freeing Arrays...\n");
  	
  	free_globals();
	free(mat); free(evecs); free(evals); free(bcmat); free(sigmap);
	
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