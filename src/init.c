#include "eigen.h"

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

	
		c2[i] = .05 * pow(r[i],-.5*flare_index);
	
		scaleH[i] = c2[i]/omega[i];
		c2[i] *= c2[i];

	
//		sigfac = .05 * pow(2.0,2 - .5*flare_index - 1.5)/(2 * M_PI * bump_function(2.0));
//		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-2.0);
		sigfac = .05 /(2 * M_PI * bump_function(2.0));
		sigma[i] = sigfac * bump_function(r[i]) * pow(r[i],-(1.5 + .5*flare_index));
	
		
		temp[i] = c2[i];

#endif // MLIN
#endif // HEEMSKERK
#endif // PAPALOIZOU


/*	Set pressure and temperature */

		temp[i] = adi_gam * c2[i];
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
	for(i=0;i<N;i++) {
		dphi0dr[i] = 0;
		for(j=0;j<N;j++) {
//			dphi0dr[i] -= weights[j]*kernel0[j+N*i]*sigma[j];
			dphi0dr[i] -= kernel0[j+N*i]*sigma[j];
		}
		if (isnan(dphi0dr[i]) != 0) {
			printf("\n\n Detected NaN in dphi0dr at i=%d, r=%.3lg\n\n", i, r[i]);
			return -1;
		}

	}
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

	for(i=0;i<N;i++) {
	
		kappa2[i] = pow(r[i],-3);
	
#ifdef PRESSURECORRECTION

#ifdef BAROTROPIC
		omegap2[i] = c2[i]*dlds[i]/(r[i]*r[i]);
#else
		omegap2[i] = dldpres[i] * temp[i]/(r[i]*r[i]);
#endif
		omega2[i] += omegap2[i];

#endif



#if defined(SELFGRAVITY) && defined(GRAVITYCORRECTION)
		omega2[i] += dphi0dr[i]/(r[i]);
#endif

#ifdef EXACTKAPPA
#ifdef BAROTROPIC
		kappa2[i] += c2[i]*((2 + dldc2[i]) * dlds[i] + d2lds[i])/(r[i]*r[i]);
#else
		kappa2[i] += temp[i]*((2 + dldtemp[i]) * dldpres[i] + d2ldpres[i])/(r[i]*r[i]);
#endif
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
