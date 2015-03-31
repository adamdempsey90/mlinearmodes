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
	if (argc < 11) {
		printf("\n\nToo Few Arguments!\n\n");
		return -1;
	}
#else 
	if (argc < 13) {
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
	
	Mdisk = atof(argv[4]);
	eps = atof(argv[5]);
	h0 = atof(argv[6]);
	sigma_index = atof(argv[7]);
	flare_index = atof(argv[8]);

#ifdef VISCOSITY
	alpha_s = atof(argv[9]);
#else
	alpha_s = 0;
#endif

#if  defined(COOLING) || defined(ADIABATIC) 
	adi_gam = atof(argv[11]);
	beta_cool = atof(argv[12]);
#endif

#ifdef ADIABATIC 
	beta_cool = 0;
#endif

	
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
	omp_set_num_threads(atoi(argv[10]));
	printf("\t\t\tOpenMP threads = %d\n", atoi(argv[10]));
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
	kernel = (double complex *)malloc(sizeof(double complex)*N*N);
	kernel0 = (double *)malloc(sizeof(double)*N*N);
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
	kappa = (double complex *)malloc(sizeof(double complex)*N);
	kappa2 = (double *)malloc(sizeof(double)*N);
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
	free(kappa2);
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


/* Set up special sigma and c^2 profiles */
#ifdef PAPALOIZOU		
		c2[i] = scaleH[i] * scaleH[i] / (r[i]*r[i]*r[i])*(1 - pow(r[i],-10))*(1-pow(r[i]/rout,10));
		sigma[i] = pow(c2[i],1.5);
#else
#ifdef HEEMSKERK

		sigma[i] = pow(r[i]-ri + .05,2)*pow(ro - r[i] +.05,2.5);
		c2[i] = sigma[i]*h0*h0;


#else
#ifdef MLIN

	
		c2[i] = .05 * pow(r[i],-.5*flare_index);
	
		scaleH[i] = c2[i]/omega[i];
		c2[i] *= c2[i];

	
//		sigfac = .05 * pow(2.0,2 - .5*flare_index - 1.5)/(2 * M_PI * bump_function(2.0));
//		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-2.0);
		sigfac = .05 /(2 * M_PI * bump_function(2.0));
		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-(1.5 + .5*flare_index));
	
		
		

#else
		sigma[i] = pow(r[i],sigma_index);
		c2[i] = scaleH[i]*scaleH[i] * omega2[i];	
#endif // MLIN
#endif // HEEMSKERK
#endif // PAPALOIZOU


/*	Set pressure and temperature */

		temp[i] = c2[i];
		pres[i] = sigma[i] * temp[i];

	
/*	Set viscosity if enabled, if not then make sure it's zero */
		
#ifdef VISCOSITY		
		nu[i] = alpha_s * scaleH[i]*scaleH[i]*omega[i];
		
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
	
	printf("sig tot = %lg\n", 2*M_PI*sigfac);
	sigfac = Mdisk/(2*M_PI*sigfac);
	
	printf("sig tot = %lg\n", sigfac);
#else
	sigfac = 1;
#endif	
	printf("sigfac = %lg\n",sigfac);
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
#endif	
	for(i=0;i<N;i++) {
		dphi0dr[i] = 0;
		for(j=0;j<N;j++) {
			dphi0dr[i] -= weights[j]*kernel0[j+N*i]*sigma[j];
		}
		if (isnan(dphi0dr[i]) != 0) {
			printf("\n\n Detected NaN in dphi0dr at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}

	}
#endif
	

	printf("Calculating Derivatives\n");
	
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


/* Correct the background rotation for pressure and self gravity.
	Set epicyclic frequency 
*/

	for(i=0;i<N;i++) {

	
#ifdef PRESSURECORRECTION


#ifdef BAROTROPIC
#ifdef EXACTKAPPA
		omegap2[i] = c2[i]*sigma_index/(r[i]*r[i]);
#else
		omegap2[i] = c2[i]*dlds[i]/(r[i]*r[i]);
#endif
#else
#ifdef EXACTKAPPA
		omegap2[i] = (sigma_index + 2*flare_index - 1)*temp[i]/(r[i]*r[i]);
#else
		omegap2[i] = dldpres[i] * temp[i]/(r[i]*r[i]);
#endif
#endif
		omega2[i] += omegap2[i];
#endif



#if defined(SELFGRAVITY) && defined(GRAVITYCORRECTION)
		omega2[i] += dphi0dr[i]/(r[i]*r[i]);
#endif

#ifdef EXACTKAPPA
		kappa2[i] = omega2[i] + (2*flare_index +1)*omegap2[i];
#else
		kappa2[i] = 4*omega2[i];
#endif
	}

#ifndef EXACTKAPPA	
	matvec(D,omega2,kappa2,1,1,N);
#endif	
	for(i=0;i<N;i++) {
	
#ifdef PRESSURECORRECTION

		kappa[i] = csqrt(kappa2[i]);
#else	

		kappa[i] = omega[i];
		kappa2[i] = kappa[i]*kappa[i];
#endif
		omega[i] = sqrt(omega2[i]);
		lom[i] = log(omega[i]);
		if (isnan(kappa[i]) != 0) {
			printf("\n\n Detected NaN in kappa at i=%d, r=%.3lg, kap2=%.3lg\n\n", i, r[i],kappa2[i]);
			return -1;
		}
		if (isnan(omega[i]) != 0) {
			printf("\n\n Detected NaN in omega at i=%d, r=%.3lg, om2=%.3lg\n\n", i, r[i],omega2[i]);
			return -1;
		}

	}
	
	matvec(D,lom,dldom,1,0,N);
	matvec(D2,lom,d2dom,1,0,N);
	
	
	
	return 0;
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
		
			
			KD[indx] = G*sigma[j]*D[indx];
//			HL[indx] = G*D[indx];
			HL[indx] = 0;
			mat[indx] = 0;
			bcmat[indx] = 0;
			
			
			if (i==j) {
				mat[indx] += A;
				bcmat[indx] += 1;
//				KD[indx] = -sigma[j]* ( dlds[j] + D[indx]);
//				HL[indx] = G*(2. + D[indx]);
				KD[indx] += G*sigma[j]*dlds[j];
//				HL[indx] += 2*G;
			}
// 			else {
// 				mat[indx] = 0;
// 				bcmat[indx] = 0;
// 				KD[indx] = -D[indx]*sigma[j];
// 				HL[indx] = G*D[indx];
// 			}
			
			mat[indx] += B*D[indx] + C*D2[indx];
			
			H[indx] = 0;
			DKD[indx] = 0;
			
		
		}
	}
			
			
			
//	matmat(kernel,D,KD,1,0,N);
//	matmat(D,KD,DKD,1,0,N);
	
//	cmatmat(HL,H,DKD,1,0,N);
	

#ifdef SELFGRAVITY	
	cmatmat(kernel,KD,DKD,1,0,N);
	for(i=0;i<N*N;i++) mat[i] += DKD[i];
#endif


/* Set Boundary Conditions */

	lagrangian_pressure_bc_inner(mat, bcmat);
	lagrangian_pressure_bc_outer(mat, bcmat);

	return 0;
}


void calc_coefficients(int i, double complex *A, double complex *B, double complex *C, double complex *G) {



#ifdef BAROTROPIC
	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) * ( dlds[i]*(2 + dldc2[i]) + d2lds[i]) + omega[i]-kappa[i];
	*B = (*C) * (2 + dldc2[i] + dlds[i]) ;
	


#endif


#ifdef ISOTHERMAL

	*C = c2[i]/(2*omega[i]*r[i]*r[i]);
	
	*A = (*C) *( 2 * dlds[i] + d2lds[i]) + omega[i]-kappa[i];
	*B = (*C) *( 2 + dlds[i]);
	


#endif




#ifdef ADIABATIC

	*C = temp[i]/(2*omega[i]*r[i]*r[i]);

	*B = (*C) * adi_gam*(2 + dldpres[i]);
	
	*A = (*C) * ( (2 + dldpres[i])*dldpres[i] + d2ldpres[i]) + omega[i] - kappa[i];

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
	
	*A = omega[i] - kappa[i] + (*C) * ( 2 * dlds[i] + d2lds[i] 
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
	*G = -1.0/(2*omega[i]*r[i]*r[i]);
#else
	*G = 0;
#endif


#ifdef VISCOSITY
	double complex temp, norm;
	double q = dldom[i];
	double qp = d2dom[i];
	double beta = dlds[i];
	double gam = dldnu[i];
	double betap = d2lds[i];
	
	norm = I/(r[i]*omega[i]);
	
	temp = -2.*(12.+9.*beta + 9.*gam + 7.*qp) + q*(-31.+4.*beta*(-2.+3*beta)-14.*gam-12.*betap);

	*A += temp*(nu[i]*omega[i]*norm/(12*r[i]));
	
	temp = 31. + 14.*beta + 2*q*(11. + 6*beta) + 14.*gam;
	
	*B += -temp*norm*(nu[i]*omega[i]/(12*r[i]));

	*C += -(nu[i]*omega[i]*norm/(6*r[i]))*(7 + 6*q);
#endif

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
	result *= weights[j]*r[j]*r[j];

#ifdef INDIRECT
	result -= 3*M_PI*r[i]*weights[j];
#endif	
	return result;


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
			mat[indx] -= dldc2[j];
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
			mat[indx] -= dldc2[j];
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
				   %.12lg\t%.12lg\t%.12lg\n",
			lr[i],
			r[i],
			omega[i],
			c2[i],
			sigma[i],
			scaleH[i]/r[i],
			pres[i],
			temp[i],
			eps*scaleH[i],
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
	
	
	double delta1 = 5.0 * .05 * pow(1.0,1.5-.5*flare_index);
	double delta2 = 5.0 * .05 * pow(2.0,1.5-.5*flare_index);

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



