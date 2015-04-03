#include "eigen.h"

int init(double ri,double ro) {

	int i,indx;
	double sigfac;
	
	
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
		pres[i] = temp[i] * sigma[i];
	}
		
/* Set up special sigma and c^2 profiles */

#ifdef PAPALOIZOU
	set_papaloizou_profile();
#endif


#ifdef HEEEMSKERK
	set_papaloizou_profile();
#endif

#ifdef MLIN
	set_papaloizou_profile();
#endif

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
	}
	
		

	for(i=0;i<N;i++) {
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
		

	}

	
	



	printf("Calculating Derivatives\n");


#ifdef EXACTKAPPA
	for(i=0;i<N;i++) {
	
		dlds[i] = dlogsigma_dlogr(r[i]);
		d2lds[i] = dlogsigma2_dlogr2(r[i]);
	
		dldc2[i] = 2 * flare_index - 1;
		dldtemp[i] = dldc2[i];
		d2ldtemp[i] = 0;

		
		dldpres[i] = dlds[i] + dldtemp[i];	
		d2ldpres[i] = d2lds[i];
		
	}
	

#else	
	matvec(D,lsig,dlds,1,0,N);
	matvec(D,lc2,dldc2,1,0,N);
	
	matvec(D2,lsig,d2lds,1,0,N);


	matvec(D,lpres,dldpres,1,0,N);
	matvec(D2,lpres,d2ldpres,1,0,N);
	matvec(D,ltemp,dldtemp,1,0,N);
	matvec(D2,ltemp,d2ldtemp,1,0,N);

#endif


/* Correct the background rotation for pressure and self gravity.
	Set epicyclic frequency 
*/
	
	
	calc_epicyclic();
	
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


void calc_epicyclic(void) {
	int i;

	for(i=0;i<N;i++) {
		omega2[i] = pow(r[i],-3);
		kappa2[i] = omega2[i];
	}


#ifdef PRESSURECORRECTION	
	pressure_omega_correction();
#endif


#if defined(SELFGRAVITY) && defined(GRAVITYCORRECTION)
	sg_omega_correction();

#endif

	for(i=0;i<N;i++) {		
		kappa[i] = sqrt(kappa2[i]);
		omega[i] = sqrt(omega2[i]);
		lom[i] = log(omega[i]);
	}
	
	return;

}