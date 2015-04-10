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
void alloc_globals(void) {
	
	D = (double *)malloc(sizeof(double)*N*N);
	D2 = (double *)malloc(sizeof(double)*N*N);
	kernel = (double *)malloc(sizeof(double)*N*N);
	kernel0 = (double *)malloc(sizeof(double)*N*N);
	kernel02 = (double *)malloc(sizeof(double)*N*N);
	weights = (double *)malloc(sizeof(double)*N);
	KD = (double complex*)malloc(sizeof(double complex)*N*N);
	DKD = (double complex*)malloc(sizeof(double complex)*N*N);
	H = (double complex *)malloc(sizeof(double complex)*N*N);
	HL = (double complex*)malloc(sizeof(double complex)*N*N);
	cwork = (double complex *)malloc(sizeof(double complex)*N*N);
	work = (double *)malloc(sizeof(double)*N*N);
	
	r = (double *)malloc(sizeof(double)*N);
	lr = (double *)malloc(sizeof(double)*N);
	c2 = (double *)malloc(sizeof(double)*N);
	sigma = (double *)malloc(sizeof(double)*N);
	omega = (double *)malloc(sizeof(double)*N);
	omega2 = (double *)malloc(sizeof(double)*N);
	omegap2 = (double *)malloc(sizeof(double)*N);
	omegag2 = (double *)malloc(sizeof(double)*N);
	kappa = (double complex *)malloc(sizeof(double complex)*N);
	kappa2 = (double *)malloc(sizeof(double)*N);
	kappap2 = (double *)malloc(sizeof(double)*N);
	kappag2 = (double *)malloc(sizeof(double)*N);
	scaleH = (double *)malloc(sizeof(double)*N);
	dlds = (double *)malloc(sizeof(double)*N);
	dldc2 = (double *)malloc(sizeof(double)*N);
	lsig = (double *)malloc(sizeof(double)*N);
	lc2 = (double *)malloc(sizeof(double)*N);
	dphi0dr=(double *)malloc(sizeof(double)*N);
	d2lds = (double *)malloc(sizeof(double)*N);
	lom = (double *)malloc(sizeof(double)*N);
	dldom = (double *)malloc(sizeof(double)*N);
	d2dom = (double *)malloc(sizeof(double)*N);

	pres = (double *)malloc(sizeof(double)*N);
	lpres = (double *)malloc(sizeof(double)*N);
	temp = (double *)malloc(sizeof(double)*N);
	ltemp = (double *)malloc(sizeof(double)*N);
	dldpres = (double *)malloc(sizeof(double)*N);
	d2ldpres = (double *)malloc(sizeof(double)*N);
	dldtemp = (double *)malloc(sizeof(double)*N);
	d2ldtemp = (double *)malloc(sizeof(double)*N);


	nu = (double *)malloc(sizeof(double)*N);
	dldnu = (double *)malloc(sizeof(double)*N);
	lnu = (double *)malloc(sizeof(double)*N);

	omega_prec =  (double *)malloc(sizeof(double)*N);

	return;

}
void free_globals(void) {

	free(D);
	free(D2);
	free(kernel);
	free(kernel0);
	free(kernel02);
	free(weights);
	free(KD);
	free(DKD);
	free(H);
	free(HL);
	free(cwork);
	free(work);
	
	free(r);
	free(lr);
	free(c2);
	free(sigma);
	free(omega);
	free(kappa);
	free(omega2);
	free(omegap2);
	free(omegag2);
	free(kappa2);
	free(kappap2);
	free(kappag2);
	free(scaleH);
	free(dlds);
	free(dldc2);
	free(lsig);
	free(lc2);
	free(dphi0dr);
	free(d2lds);
	free(lom);
	free(dldom);
	free(d2dom);
	
//#ifndef ISOTHERMAL
	free(pres);
	free(lpres);
	free(temp);
	free(d2ldpres);
	free(dldpres);
	free(ltemp);
	free(d2ldtemp);
	free(dldtemp);
//#endif
	
	free(nu); 
	free(dldnu);
	free(lnu);

	free(omega_prec);


	return;

}

int init(double ri,double ro) {

	int i,j,indx;
	double sigfac;
	double phifac = 0;
	
	
	for(i=0;i<N;i++) {
		lr[i] = log(ri) + i*dlr;
		r[i] = exp(lr[i]);
		
		scaleH[i] = h0*r[i] * pow(r[i],flare_index);


//		sigma[i] = sigma_profile(r[i],scaleH[i],30., -poly_n); //.5*(ri + exp( log(ri) + (N-1)*dlr)),-poly_n);
		
//		c2[i] = scaleH[i]*scaleH[i] / (r[i]*r[i]*r[i]);
		
		
		omega[i] = pow(r[i],-1.5);
		omega2[i] = omega[i]*omega[i];
		sigma[i] = sigma_function(r[i]);
#ifndef INPUTMASS
		sigma[i] *= sigma0;
#endif
		c2[i] = scaleH[i]*scaleH[i] * omega2[i];	
		temp[i] = c2[i]/adi_gam;
		
/* Set up special sigma and c^2 profiles */
#ifdef PAPALOIZOU		
		c2[i] = scaleH[i] * scaleH[i] / (r[i]*r[i]*r[i])*(1 - pow(r[i],-10))*(1-pow(r[i]/rout,10));
		sigma[i] = pow(c2[i],1.5);
		temp[i] = c2[i];
#else
#ifdef HEEMSKERK

		sigma[i] = pow(r[i]-ri + .05,2)*pow(ro - r[i] +.05,2.5);

#ifndef INPUTMASS
		sigma[i] *= sigma0;
#endif

		c2[i] = sigma[i]*h0*h0;
		temp[i] = c2[i];

#else
#ifdef MLIN

	
		c2[i] = h0 * pow(r[i],-.5*flare_index);
	
		scaleH[i] = c2[i]/omega[i];
		c2[i] *= c2[i];

	
//		sigfac = .05 * pow(2.0,2 - .5*flare_index - 1.5)/(2 * M_PI * bump_function(2.0));
//		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-2.0);
		sigfac = h0 /(2 * M_PI * bump_function(2.0));
		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-(1.5 + .5*flare_index));
	
		
		temp[i] = c2[i];

#endif // MLIN
#endif // HEEMSKERK
#endif // PAPALOIZOU


/*	Set pressure and temperature */

//		temp[i] = adi_gam * c2[i];
		pres[i] = sigma[i] * temp[i];

	
/*	Set viscosity if enabled, if not then make sure it's zero */
		
#ifdef VISCOSITY		
//		nu[i] = alpha_s * scaleH[i]*scaleH[i]*omega[i];
		nu[i] = alpha_s * temp[i]/omega[i];		
		if (alpha_s == 0) {
			lnu[i] = 0;
		}
		else { 
			lnu[i] = log(nu[i]);
		}
		
#else
		nu[i] = 0;
		lnu[i] = 0;
#endif



		lc2[i] = log(c2[i]);
		ltemp[i] = log(temp[i]);
		
		dlds[i] = 0;
		dldc2[i] = 0;
		d2lds[i] = 0;
		d2dom[i] = 0;
		dldom[i] = 0;
		dldnu[i] = 0;		
		dldpres[i] = 0;
		d2ldpres[i] = 0;
		dldtemp[i] = 0;
		d2ldtemp[i] = 0;
/* Check for bad values of c^2 and sigma */
		if (isnan(c2[i]) != 0) {
			printf("\n\n Detected NaN in c2 at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
		
		if (isnan(lc2[i]) != 0) {
			printf("\n\n Detected NaN in lc2 at i=%d, r=%.3lg, c2=%lg\n\n", i, r[i],c2[i]);
			return -1;
		}
		
		if (isnan(lsig[i]) != 0) {
			printf("\n\n Detected NaN in lsig at i=%d, r=%.3lg, sig=%lg\n\n", i, r[i],sigma[i]);
			return -1;
		}
		

	}

/* Normalize the density to give the prescribed disk mass, unless using Lin's profile which 
	is already normalized.
*/

#ifndef MLIN
	sigfac=0;
	for(i=0;i<N;i++) {
		sigfac += weights[i]*sigma[i]*r[i]*r[i];	

	}
	
	
#ifdef INPUTMASS
	sigfac = Mdisk/(2*M_PI*sigfac);	
#else
	Mdisk = 2*M_PI*sigfac;
	sigfac = 1;
	printf("Total Disk Mass = %lg\n", Mdisk);
#endif
	
	
#else
	sigfac = 1;
#endif	


	for(i=0;i<N;i++) {
		sigma[i] *= sigfac;
		lsig[i] = log(sigma[i]);
		pres[i] *= sigfac;
		lpres[i] = log(pres[i]);
#ifdef HEEMSKERK
		c2[i] *= sigfac;
		lc2[i] = log(c2[i]);
		temp[i] *= sigfac;
		ltemp[i] = log(temp[i]);
#endif
		if (isnan(sigma[i]) != 0) {
			printf("\n\n Detected NaN in sigma at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}
	}	
	

/* Initialize Kernels if we're not reading it from a file */
	
	printf("Calculating Kernels\n");

#ifdef SELFGRAVITY


#ifdef READKERNEL

/* Read the kernel from the kernel.dat and kernel0.dat files written by output_kerenels
	and found in the execution directory */

	read_kernel();

#else

// #ifdef OPENMP
// #pragma omp parallel private(indx,i,j) shared(kernel,kernel0,N)
// #pragma omp for schedule(static)
// #endif
// 	for(indx=0;indx<N*N;indx++) {
// 		i = indx/N;
// 		j = indx - i*N;
// 		
// 		kernel[indx] = Kij(i,j);
// 		kernel0[indx] = K0ij(i,j);
// 	}
	compute_kernels();
#endif	
// 	for(i=0;i<N;i++) {
// 		dphi0dr[i] = 0;
// 		for(j=0;j<N;j++) {
// //			dphi0dr[i] -= weights[j]*kernel0[j+N*i]*sigma[j];
// 			dphi0dr[i] -= kernel0[j+N*i]*sigma[j];
// 		}
// 		if (isnan(dphi0dr[i]) != 0) {
// 			printf("\n\n Detected NaN in dphi0dr at i=%d, r=%.3lg\n\n", i, r[i]);
// 			return -1;
// 		}
// 
// 	}
#endif
	

	printf("Calculating Derivatives\n");


#ifdef EXACTKAPPA
	for(i=0;i<N;i++) {
	
		dlds[i] = dlogsigma_dlogr(r[i]);
		d2lds[i] = dlogsigma2_dlogr2(r[i]);
	
		dldc2[i] = 2 * flare_index - 1;
		dldtemp[i] = dldc2[i];
		d2ldtemp[i] = 0;
#ifdef VISCOSITY
		dldnu[i] = 2 * flare_index + .5;
#endif
		
		dldpres[i] = dlds[i] + dldtemp[i];	
		d2ldpres[i] = d2lds[i];
		
	}
	

#else	
	matvec(D,lsig,dlds,1,0,N);
	matvec(D,lc2,dldc2,1,0,N);
	
	matvec(D2,lsig,d2lds,1,0,N);
#ifdef VISCOSITY
	matvec(D,lnu,dldnu,1,0,N);
#endif

	matvec(D,lpres,dldpres,1,0,N);
	matvec(D2,lpres,d2ldpres,1,0,N);
	matvec(D,ltemp,dldtemp,1,0,N);
	matvec(D2,ltemp,d2ldtemp,1,0,N);

#endif


/* Correct the background rotation for pressure and self gravity.
	Set epicyclic frequency 
*/

// 	for(i=0;i<N;i++) {
// 	
// 		kappa2[i] = pow(r[i],-3);
// 	
// #ifdef PRESSURECORRECTION
// 
// #ifdef BAROTROPIC
// 		omegap2[i] = c2[i]*dlds[i]/(r[i]*r[i]);
// #else
// 		omegap2[i] = dldpres[i] * temp[i]/(r[i]*r[i]);
// #endif
// 		omega2[i] += omegap2[i];
// 
// #endif
// 
// 
// 
// #if defined(SELFGRAVITY) && defined(GRAVITYCORRECTION)
// 		omega2[i] += dphi0dr[i]/(r[i]);
// #endif
// 
// #ifdef EXACTKAPPA
// #ifdef BAROTROPIC
// 		kappa2[i] += c2[i]*((2 + dldc2[i]) * dlds[i] + d2lds[i])/(r[i]*r[i]);
// #else
// 		kappa2[i] += temp[i]*((2 + dldtemp[i]) * dldpres[i] + d2ldpres[i])/(r[i]*r[i]);
// #endif
// #else
// 		kappa2[i] = 4*omega2[i];
// #endif
// 	}
// 
// #ifndef EXACTKAPPA	
// 	matvec(D,omega2,kappa2,1,1,N);
// #endif	
// 	for(i=0;i<N;i++) {
// 	
// #ifdef PRESSURECORRECTION
// 
// 		kappa[i] = csqrt(kappa2[i]);
// #else	
// 
// 		kappa[i] = omega[i];
// #endif
// 		omega[i] = sqrt(omega2[i]);
// 		lom[i] = log(omega[i]);
// 		if (isnan(kappa[i]) != 0) {
// 			printf("\n\n Detected NaN in kappa at i=%d, r=%.3lg, kap2=%.3lg\n\n", i, r[i],kappa2[i]);
// 			return -1;
// 		}
// 		if (isnan(omega[i]) != 0) {
// 			printf("\n\n Detected NaN in omega at i=%d, r=%.3lg, om2=%.3lg\n\n", i, r[i],omega2[i]);
// 			return -1;
// 		}
// 
// 	}
	
	
//	calc_epicyclic();

#ifndef INFINITEDISK
	calc_omega_prec_grav();
#endif
	calc_omega_prec_pres();
		
	matvec(D,lom,dldom,1,0,N);
	matvec(D2,lom,d2dom,1,0,N);
	
	
	
	return 0;
}


void calc_epicyclic(void) {
	int i,j;


	for(i=0;i<N;i++) {
	
#ifdef PRESSURECORRECTION	

#ifdef BAROTROPIC

		omegap2[i] = c2[i]*dlds[i]/(r[i]*r[i]);
		
		
#ifdef EXACTKAPPA
		kappap2[i] = c2[i]*((2 + dldc2[i]) * dlds[i] + d2lds[i])/(r[i]*r[i]);
#else
		kappap2[i] = 4*omegap2[i];
#endif


#else

		omegap2[i] = dldpres[i] * temp[i]/(r[i]*r[i]);
		
		
#ifdef EXACTKAPPA
		kappap2[i] = temp[i]*((2 + dldtemp[i]) * dldpres[i] + d2ldpres[i])/(r[i]*r[i]);
#else
		kappap2[i] = 4*omegap2[i];
#endif

#endif // BAROTROPIC

#else
		omegap2[i] = 0;
		kappap2[i] = 0;
		
#endif // PRESSURECORRECTION


#if defined(SELFGRAVITY) && defined(GRAVITYCORRECTION)
		
		omegag2[i] = 0;
		kappag2[i] = 0;
		dphi0dr[i] = 0;
		for(j=0;j<N;j++) {
			omegag2[i] += kernel0[j+i*N]*sigma[j];
			kappag2[i] += kernel02[j+i*N]*sigma[j];
			dphi0dr[i] += kernel[j+i*N]*sigma[j];
		}
//		omegag2[i] *= (-1.0/r[i]);
//		kappag2[i] *= -pow(r[i],-3);
#else
		omegag2[i] = 0;
		kappag2[i] = 0;		

#endif
		
		
	
	}
	
#ifndef EXACTKAPPA
	matvec(D,omegap2,kappap2,1,1,N);
#endif

	for(i=0;i<N;i++) {
	
		omega2[i] = pow(r[i],-3) + omegap2[i] + omegag2[i];
		kappa2[i] = pow(r[i],-3) + kappap2[i] + kappag2[i];
		
		kappa[i] = sqrt(kappa2[i]);
		omega[i] = sqrt(omega2[i]);
	}
	
	output_omega_corrections(omegap2,omegag2,kappap2,kappag2);
	

	return;

}


void output_omega_corrections(double *omegap2, double *omegag2, double *kappap2, double *kappag2) {
	FILE *f;
	int i;
	
	f = fopen("omegacorrecs.dat","w");
	
	for(i=0;i<N;i++) {
		fprintf(f,"%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\n",
				r[i],
				pow(r[i],-3),
				omegap2[i],
				omegag2[i],
				kappap2[i],
				kappag2[i]);
	}
	fclose(f);
	
	return;

}



double sigma_profile(double rval, double h, double rc, double beta) {
	double sigma_min, index, width,result, x;
	
	x = rval - rc;
	sigma_min = 1e-1;
	index = 2;
	width = 2*h;
	
//	printf("%lg\t%lg\t%lg\t%lg\n", x,sigma_min,index,width);
	result = (1 - (1-sigma_min)/(1 + pow(x/width,index)*exp(-(2.5*width/x)*(2.5*width/x))));
	
//	printf(" %lg\n", result);
	result *= pow(rval,beta);
	
	return result;
}


#ifdef SECONDORDER
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

#endif



#ifndef COMPTRAPZ
void calc_weights(void) {
/* O(N^(-4)) weights for numerical quadrature*/
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
#else
void calc_weights(void) {
/* Composite Trapezoid numerical quadrature*/
	int i;
	weights[0] = .5*dlr;
	weights[N-1] = .5*dlr;
	
	for(i=1;i<N-1;i++) {
		weights[i] = dlr;
	}

	return;

}


#endif

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
	
	double complex A, B, C, G;
/* Compute the matrix including all of its component matrices */

	
// #ifdef OPENMP
// #pragma omp parallel 
// #pragma omp for private(i,A,B,C,indx,j,
// #endif	

	for(i=0;i<N;i++ ) {
		calc_coefficients(i,&A,&B,&C,&G);

		for(j=0 ; j<N ; j++) {
			indx = j + N*i;
		
			
//			KD[indx] = G*sigma[j]*D[indx];


			mat[indx] = 0;
			bcmat[indx] = 0;
			
			
			if (i==j) {
				mat[indx] += A;
				bcmat[indx] += 1;

//				KD[indx] += G*sigma[j]*dlds[j];

			}
			
			mat[indx] += B*D[indx] + C*D2[indx];
	
		
		}
	}
			
			
			
//	matmat(kernel,D,KD,1,0,N);
//	matmat(D,KD,DKD,1,0,N);
	
//	cmatmat(HL,H,DKD,1,0,N);
	
	
#ifdef SELFGRAVITY	
	int l;
//	cmatmat(kernel,KD,DKD,1,0,N);
//	for(i=0;i<N*N;i++) mat[i] += DKD[i];
	for(i=0;i<N;i++) {
		G = 1.0/(2*omega[i]*r[i]*r[i]*r[i]);
		for(j=0;j<N;j++) {
			indx = j + N*i;
			KD[indx] = 0;
			for(l=0;l<N;l++) {
				KD[indx] += -sigma[j]*kernel[l + N*i]*D[l*N+j];
				if (l==j) {
					KD[indx] += -sigma[j]*kernel[l+N*i]*dlds[j];
				}	
			}
			mat[indx] += G * KD[indx];
		}
		
	}

#ifdef EDGEGRAVITY
	add_edge_sg(mat);
#endif

#endif

	

/* Set Boundary Conditions */

	lagrangian_pressure_bc_inner(mat, bcmat);
	lagrangian_pressure_bc_outer(mat, bcmat);

	return 0;
}


void calc_coefficients(int i, double complex *A, double complex *B, double complex *C, double complex *G) {



#ifdef BAROTROPIC
	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) * ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega_prec[i];
	*B = (*C) * (2 + dldc2[i] + dlds[i]) ;
	


#endif


#ifdef ISOTHERMAL

	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) *( 2 * dlds[i] + d2lds[i]) + omega_prec[i];
	*B = (*C) *( 2 + dlds[i]);
	


#endif




#ifdef ADIABATIC

	*C = temp[i]/(2*omega[i]*r[i]*r[i]);

	*B = (*C) * adi_gam*(2 + dldpres[i]);
	
	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega_prec[i];

	*C *= adi_gam;

#endif


#ifdef COOLING

//	double complex cool_fac = ( 1 + 1j*beta_cool)/(1 + beta_cool*beta_cool);
	double complex cool_fac = ( 1 + 1j*beta_cool*(adi_gam-1))/(1 + beta_cool*beta_cool*(adi_gam-1)*(adi_gam-1));
	
// 	*C = cool_fac *  temp[i]/(2*omega[i]*r[i]*r[i]);
// 
// 	*B = (*C) * ( adi_gam*(2 + dldpres[i]) + 1j*beta_cool*dldpres[i]);
// 	
// 	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];
// 
// 	*C *= adi_gam;

	*C = temp[i]/(2 * omega[i] * r[i] * r[i]);
	
	*A = omega_prec[i] + (*C) * ( 2 * dlds[i] + d2lds[i] 
				  + cool_fac * (d2ldtemp[i] + dldtemp[i]*(2 + dlds[i] + dldtemp[i])));
				  
	*B = (*C) * (2 + dlds[i] + cool_fac * ( adi_gam*dldtemp[i] + ( adi_gam - 1)*(dlds[i] + 2)));
	
	*C *= 1 + cool_fac * (adi_gam -1);
	
#endif







// #ifdef BAROTROPIC
// 	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
// 	
// 	*A = ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega[i]-kappa[i];
// 	*B = (2 + dldc2[i] + dlds[i]) ;
// 	
// 	*A *= (*C);
// 	*B *= (*C);
// #else
// 
// #ifdef ISOTHERMAL
// 
// 	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
// 	
// 	*A = 2 * dlds[i] + d2lds[i] + omega[i]-kappa[i];
// 	*B = 2 + dlds[i];
// 	
// 	*A *= (*C);
// 	*B *= (*C);
// 	
// #else // ADIABATIC OR COOLING
// 	double complex cool_fac = ( 1 + 1j*beta_cool)/(1 + beta_cool*beta_cool);
// 	
// 	*C = cool_fac *  temp[i]/(2*omega[i]*r[i]*r[i]);
// 
// 	*B = (*C) * ( adi_gam*(2 + dldpres[i]) + 1j*beta_cool*dldpres[i]);
// 	
// 	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];
// 
// 	*C *= adi_gam;
// 
// #endif
// #endif


#ifdef SELFGRAVITY
	*G = -1.0/(2*omega[i]*r[i]*r[i]*r[i]);
#else
	*G = 0;
#endif


#ifdef VISCOSITY
	double complex norm;
// 	double q = dldom[i];
// 	double qp = d2dom[i];
// 	double beta = dlds[i];
// 	double gam = dldnu[i];
// 	double betap = d2lds[i];
	
// 	norm = I/(r[i]*omega[i]);
// 	
// 	temp = -2.*(12.+9.*beta + 9.*gam + 7.*qp) + q*(-31.+4.*beta*(-2.+3*beta)-14.*gam-12.*betap);
// 
// 	*A += temp*(nu[i]*omega[i]*norm/(12*r[i]));
// 	
// 	temp = 31. + 14.*beta + 2*q*(11. + 6*beta) + 14.*gam;
// 	
// 	*B += -temp*norm*(nu[i]*omega[i]/(12*r[i]));
// 
// 	*C += -(nu[i]*omega[i]*norm/(6*r[i]))*(7 + 6*q);

	norm = temp[i]/(omega[i]*r[i]*r[i]);
/* Shear Viscosity */

// 	*A -= I * (alpha_s / 6.) * norm * (33.5 + dldtemp[i] - 2*dlds[i] + 18*d2lds[i]);
// 	*B -= I * (alpha_s / 3.) * norm * (-9.5 - 7*dldtemp[i] + 2*dlds[i]);
// 	*C -= I * (2 *alpha_s / 3.) * norm;

	*A -= I*alpha_s*(norm/4)*(9 + 8*dldtemp[i]  + dlds[i] + 6*d2lds[i]); 
	*B += I*(norm/6)*(3*alpha_b*(2 + 2*dldtemp[i] + dlds[i])+alpha_s*(32+24*dldtemp[i]+3*dlds[i]));
	*C += I*(norm/6)*(3*alpha_b - 8*alpha_s);

/* Bulk Viscosity */

// 	*B += I * alpha_b * norm * (2.5 + dldtemp[i] + dlds[i]);
// 	*C -= I * alpha_b * norm;

#endif


//	A = 0; B= 0; C=0; G=0;

	return; 

}


double K0ij(int i, int j) {
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
	result *= weights[j]*r[j]*r[j];

#ifdef INDIRECT
	result -= 3*M_PI*r[i]*weights[j];
#endif	
	return result;


}

void compute_kernels(void) {
	int i,j, indx;
	double r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee; 
	
#ifdef OPENMP
#pragma omp parallel private(indx,i,j,r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee) shared(kernel,kernel0,kernel02,N,r,scaleH,weights)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		

		r1 = r[i];
		r2 = r1*r1;
		r3 = r2*r1;
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
		
		kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
			
		kernel[indx] = -2*(eps4 + 2*r4 - 3*r2*rp2 + rp4 + eps2*(3*r2 + 2*rp2))*ee;
		kernel[indx] += 2*(eps2 + r_m_rp)*(eps2 + 2*r2 + rp2)*ek;
		kernel[indx] /= ((eps2 + r_m_rp)*rp1*sqrt(eps2 +r_p_rp));
		
// 		kernel0[indx] = 2*((eps2 - r2+rp2)*ee - (eps2 + r_m_rp)*ek);
// 		kernel0[indx] /= (r1*(eps2 + r_m_rp)*sqrt(eps2 + r_p_rp));
// 		
// 		kernel02[indx] = -ee*(eps6 + eps4*(-2*r2+3*rp2)+pow(rp3-rp1*r2,2)+eps2*(-3*r4-4*r2*rp2+3*rp4));
// 		kernel02[indx] += ek*(eps2+r_m_rp)*(eps4+r_m_rp*r_p_rp + eps2*(r2+2*rp2));
// 		kernel02[indx] *= -4*r1;
// 		kernel02[indx] /= ( pow(eps2+r_m_rp,2) * pow(eps2+r_p_rp,1.5));
		
		kernel[indx] *=  weights[j]*rp2;
// 		kernel0[indx] *=  -weights[j]*rp2/r1;
// 		kernel02[indx] *=  -weights[j]*rp2/r3;


		
	}	

	return;
}

void add_edge_sg(double complex *mat) {
	int i,j;
	double x1,x2,x4,kval,kern,rd,eps2,ee,ek, G;
	
	rd = r[N-1];
	
	j = N-1;
	
#ifdef OPENMP
#pragma omp parallel private(i,x1,x2,x4,kval, kern,eps2,ee,ek, G) shared(r,scaleH,eps,rd,mat,N)
#pragma omp for schedule(static)
#endif
	for(i=0;i<N;i++) {
		G = 1.0/(2*omega[i]*pow(r[i],3));
		x1 = r[i]/rd;
		x2 = x1*x1;
		x4 = x2*x2;
		eps2 = eps*scaleH[N-1]/rd;
		eps2 *= eps2;
		
		kval = sqrt(4*x1/(eps2 + (1+x1)*(1+x1)));
		ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
		ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
		
		kern = -2*( ((1+eps2)*(1+eps2)+ 3*(eps2-1)*x2 + 2*x4)*ee - (eps2+(x1-1)*(x1-1))*(1+eps2+2*x2)*ek);
		kern /= ((eps2 + (x1-1)*(x1-1))*sqrt(eps2+(1+x1)*(1+x1)));
		
		mat[j+N*i] += G*rd*rd*sigma[N-1]*kern;
	}
	
	
	return;

}


double kernel_integrand(int m, double r, double rp, double rs2, double phi) {
	double numerator, denominator;
	if (m==0) {
	
		return (rp*cos(phi) - r)*pow(r*r + rs2 - 2*r*rp*cos(phi),-1.5);
	}
	else if (m==1) {
//		return cos(phi)/sqrt(r*r + rs2 - 2*r*rp*cos(phi));	
		denominator = pow(r*r + rs2 - 2*r*rp*cos(phi), -1.5);
		numerator = cos(phi)*(2*rs2 + r*r - 3*r*rp*cos(phi));
		return numerator * denominator;
	}
	else {
		return 0;
	}
}

void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at inner boundary */

	int j,indx;
#ifdef COOLING
	double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
	eta_fac /= (1 - eta_fac);
#endif
	for(j=0;j<N;j++) {
		indx = j;
		mat[indx] = D[indx];
		bcmat[indx] = 0;		
		
		if (j==0) {
#ifdef COOLING
//			mat[indx] += I*beta_cool * dldpres[j] / adi_gam;
//			mat[indx] -= dldtemp[j]*(beta_cool*beta_cool - I*adi_gam*beta_cool)/(beta_cool*beta_cool + adi_gam*adi_gam);
			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
#endif	
#ifdef ISOTHERMAL
			mat[indx] -= dldtemp[j];
#endif

		}	
	}
	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {

/* Set up zero Lagrangian Pressure B.C at outer boundary */	
#ifdef COOLING
	double complex eta_fac = 1 - I*beta_cool*(adi_gam -1);
	eta_fac /= (1 - eta_fac);
#endif
	int j,indx;
	for(j=0;j<N;j++) {
		indx= j + N*(N-1);
	
		mat[indx] = D[indx];
		bcmat[indx] = 0;
		
		if (j==(N-1)) {
#ifdef COOLING
//			mat[indx] += I*beta_cool * dldpres[j] / adi_gam;
//			mat[indx] -= dldtemp[j] * (beta_cool*beta_cool - I*adi_gam*beta_cool)/(beta_cool*beta_cool + adi_gam*adi_gam);
			mat[indx] -= dldtemp[j]/(1 + eta_fac * adi_gam);
#endif	
#ifdef ISOTHERMAL
			mat[indx] -= dldtemp[j];
#endif

		}
		
	}
	return;
}

void calc_sigmap(double complex *sigmap, double complex *evecs) {
	int i,j,k,indx;
	
#ifdef OPENMP
#pragma omp parallel private(indx,i,j,k) shared(N,dlds,evecs,sigma,sigmap,D)
#pragma omp for schedule(static)
#endif
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		sigmap[indx] = -sigma[j]*evecs[indx]*dlds[j];
		
		for(k=0;k<N;k++) {
				sigmap[indx] += -sigma[j]*D[k+N*i]*evecs[k+N*i];
		}
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


//	zgeev_( &JOBVL, &JOBVR, &nA, tA, &LDA, evals, tQ, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );
	zggev_( &JOBVL, &JOBVR, &nA, tA, &LDA, tQ, &LDB, evals_alpha,evals_beta, NULL, &LDVL, evecs, &LDVR, CWORK, &LWORK, RWORK, &INFO );
	
	for(i=0;i<nA;i++) {
		if (cabs(evals_beta[i]) != 0) {
			evals[i] = evals_alpha[i]/evals_beta[i];
		}
	}
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

void output_globals(void) {
	int i;
	
	FILE *f = fopen("globals.dat","w");
	
	
	fprintf(f,"# lr \t r \t omega \t c2 \t sigma \t H/r \t soft \t dldc2 \t dlds \t kappa^2 \t d2lds \t dldom \t d2dom \t nu \t dldnu\n");
	
	
	for(i=0;i<N;i++) {
		fprintf(f,"%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t \
				   %.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t%.12lg\t \
				   %.12lg\t%.12lg\t%.12lg\t%.12lg\n",
			lr[i],
			r[i],
			omega[i],
			c2[i],
			sigma[i],
			scaleH[i]/r[i],
			pres[i],
			temp[i],
			eps*scaleH[i],
			omega_prec[i],
			dldc2[i],
			dlds[i],
			dldpres[i],
			kappa2[i],
			d2lds[i],
			d2ldpres[i],
			dldom[i],
			d2dom[i],
			nu[i],
			dldnu[i]);
	}
			
			
		
			
	fclose(f);
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

void output(double complex *evals, double complex *evecs, double complex *sigmap) {
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
	
	f = fopen("sigmap.dat","w");
	
	fprintf(f,"#sigmap along rows\n");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.12lg\t%.12lg\t",creal(sigmap[j+i*N]),cimag(sigmap[j+i*N]));
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
	
	return;
}

void output_derivatives(void) {
	FILE *f;
	
	int i,j;
	
	f=fopen("D.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.12lg\t",D[j+N*i]);
		}
		fprintf(f,"\n");
	}	
	fclose(f);
	
	f=fopen("D2.dat","w");
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			fprintf(f,"%.12lg\t",D2[j+N*i]);
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
	
	return;

}
	

#ifdef MLIN
double bump_function(double rval) {
	
	
	double delta1 = 5.0 * h0 * pow(1.0,1.5-.5*flare_index);
	double delta2 = 5.0 * h0 * pow(2.0,1.5-.5*flare_index);

	double fac1 = .5*(1-.1)*(1+tanh((rval-1.0)/delta1))+.1;
	double fac2 = .5*(1-.1)*(1-tanh((rval-2.0)/delta2))+.1;
	
	return (fac1*fac2);

}
#endif

	
	
#ifdef TESTFUNCTION

double complex test_function(double rval) {
	double complex result;
	
	result = 1 / pow(rval,2);
		
	return result;
}

void fill_mat(double complex *mat, double complex *bcmat) {
	int i,j,indx;

#ifdef OPENMP
#pragma omp parallel private(indx,i,j) shared(D2,r,mat)
#pragma omp for schedule(static)
#endif	
	for(indx=0;indx<N*N;indx++) {
		i = indx/N;
		j = indx - i*N;
		mat[indx] = D2[indx] * test_function(r[i]);
		if (i==j) {
			bcmat[indx] = 1;
		}
		else {
			bcmat[indx] = 0;
		}
	}
	
	for(j=0;j<N;j++) {
		indx= j;
	
		mat[indx] = D[indx];
		bcmat[indx] = 0;
		indx +=  N*(N-1);
		mat[indx] = D[indx];
		bcmat[indx] = 0;
	}
	
	return;
}

#endif


double sigma_function(double rval) {

	double ans = 0;

#ifdef EXPDECAY

	ans = pow(rval,sigma_index) * exp(-rval/RDECAY);

#else
#ifdef EXPDECAY2

	ans = pow(rval,sigma_index) * exp(-pow(rval/RDECAY,2));

#else
	
	ans = pow(rval,sigma_index);

#endif
#endif
	return ans;

}

double dlogsigma_dlogr(double rval) {
	double ans = 0;
#ifdef EXPDECAY

	ans = sigma_index - rval/RDECAY;
#else
#ifdef EXPDECAY2

	ans = sigma_index - 2*pow(rval/RDECAY,2);
#else

	ans = sigma_index;	

#endif
#endif
	
	return ans ;
	
}

double dlogsigma2_dlogr2(double rval) {
	double ans = 0;
#ifdef EXPDECAY
	ans = -rval/RDECAY;
#else
#ifdef EXPDECAY2

	ans = -4*pow(rval/RDECAY,2);
#else

	ans = 0;	
#endif
#endif
	return ans;
	
}


double integrand(double x, void *params) {
	double r1,r2,r3,r4,rp1,rp2,rp3,rp4,eps1,eps2,eps4,eps6,r_p_rp,r_m_rp,kval,ek,ee; 
	double ans;
	int i=  *(int *)params;
	
	r1 = r[i];
	r2 = r1*r1;
	r3 = r2*r1;
	r4 = r2*r2;
	rp1 = exp(x);
	rp2 = rp1*rp1;
	rp3 = rp2*rp1;
	rp4 = rp2*rp2;
	eps1 = eps * scaleH_func(r1);
	eps2 = eps1*eps1;
	eps4 = eps2*eps2;
	eps6 = eps4 * eps2;
	r_p_rp = r1 + rp1;
	r_p_rp *= r_p_rp;
	r_m_rp = r1 - rp1;
	r_m_rp *= r_m_rp;
	
	kval = sqrt(4*r1*rp1/(r_p_rp + eps2));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);
		
	ans=2*(eps2 - 2 *eps1 *r1 - r2 + rp2)*(eps2 + 2 *eps1*r1 - r2 + rp2)*(eps2 + r2 + rp2)*ee; 
	ans -= 2*(eps2 + (r1 - rp1)*(r1-rp1))*(eps4 + r4 + 2*(eps1 - r1)*(eps1 + r1)*rp2 + rp4)*ek;
	ans /= (pow(eps2 + (r1 - rp1)*(r1-rp1),2) *pow(eps2 + (r1 + rp1)*(r1+rp1),1.5));
	
	ans *= -rp2*sigma_func(rp1);
	
	return ans;

}

double sigma_func(double x) {

	return sigma0 * pow(x, sigma_index);

}

double temp_func(double x) {
	return h0*h0*pow(x,2*flare_index-1);
}

double omk_func(double x) {
	return pow(x,-1.5);
}
double scaleH_func(double x) {

	return h0*x*pow(x,flare_index);

}

void calc_omega_prec_grav(void) {
	int nsubintervals = 1000;
	double error;
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc( nsubintervals );
	int i,j;

	
	
	gsl_function func;
	func.function = &integrand;
	
	for(i=0;i<N;i++) {
	
		func.params = &i;
		gsl_integration_qags(&func, log(r[0]),log(r[N-1]),0,tol, nsubintervals , workspace,&omega_prec[i],&error);
		
		omega_prec[i] *= (-.5/sqrt(r[i]));
		
	}

	gsl_integration_workspace_free(workspace);
	
	
	
	
	return; 

}

void calc_omega_prec_pres(void) {
	int i;
	double norm, fac;
	double delta = 2*flare_index - 1;
	double mu = sigma_index;
	
	for(i=0;i<N;i++) {
		norm = -.5/sqrt(r[i]);
#ifdef ISOTHERMAL
	
		fac = (mu+delta)*(delta+1)*temp[i];

#endif

#ifdef COOLING
		
		fac =  (mu+delta)*(delta+1)*temp[i];

#endif

#ifdef BAROTROPIC
	
		fac = mu*(delta + 1)*c2[i];

#endif		
	
#ifdef INFINITEDISK
		fac -= sigma[i]*(1+mu)*(2+mu)*r[i]*27.5;
#endif	
		omega_prec[i] += fac*norm;
	
	}


	return;
}