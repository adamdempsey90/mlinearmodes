#include "eigen.h"


void init_globals(double ri, double ro);

void init_weights(void);
double calc_total_disk_mass(void);
void calc_omega_corrections(void);
void calc_omega_pres(double *omp2, double *kapp2);



int init(double ri,double ro) {
	init_derivatives();


	init_weights();

	init_globals(ri,ro);


	calc_omega_corrections();


	return 0;
}

void calc_omega_corrections(void) {
	int i;
	double *omp2 = (double *)malloc(sizeof(double)*N);
	double *omg2 = (double *)malloc(sizeof(double)*N);
	double *kapp2 = (double *)malloc(sizeof(double)*N);
	double *kapg2 = (double *)malloc(sizeof(double)*N);



	calc_omega_pres(omp2,kapp2);
#ifdef SELFGRAVITY
	calc_omega_sg(omg2,kapg2);
#else
	for(i=0;i<N;i++) {kapg2[i] =0; omg2[i] = 0;}
#endif
// #ifdef OPENMP
// #pragma omp for private(i)
// #endif
	for(i=0;i<N;i++) {
		kappa2[i] = omega[i]*omega[i] + kapp2[i] + kapg2[i];
		omega[i] = sqrt(omega[i]*omega[i] + omp2[i]+omg2[i]);

	}

	free(omp2); free(omg2); free(kapp2); free(kapg2);
	return;
}

void calc_omega_pres(double *omp2, double *kapp2) {
	int i;

// OPENMP
	for(i=0;i<N;i++) {
#ifdef BAROTROPIC
		omp2[i] = dlds[i]*c2[i]/(r[i]*r[i]);
		kapp2[i] = c2[i]*(dlds[i]*(2+dldc2[i]) + d2lds[i])/(r[i]*r[i]);
#endif
#if defined(ISOTHERMAL) || defined(COOLING)
		omp2[i]= temp[i]*(dldtemp[i] + dlds[i])/(r[i]*r[i]);
		kapp2[i] = temp[i]*((2+dldtemp[i])*(dldtemp[i]+dlds[i])+d2ldtemp[i]+d2lds[i])/(r[i]*r[i]);
#endif
	}

	return;

}
void init_globals(double ri, double ro) {
	int i;
// OPENMP
	for(i=0;i<N;i++) {


		lr[i] = log(ri) + i*dlr;
		r[i] = exp(lr[i]);

		omega[i] = omk_func(r[i]);
		dldom[i] = dlogomk_func(r[i]);
		d2ldom[i] = d2logomk_func(r[i]);

		scaleH[i] = scaleH_func(r[i]);
		aspect_ratio[i] = scaleH[i]/r[i];

		sigma[i] = sigma_func(r[i]);
		dlds[i] = dlogsigma_func(r[i]);
		d2lds[i] = d2logsigma_func(r[i]);

		temp[i] = temp_func(r[i]);
		dldtemp[i] = dlogtemp_func(r[i]);
		d2ldtemp[i] = d2logtemp_func(r[i]);



		c2[i] = adi_gam * temp[i];
		dldc2[i] = adi_gam * dldtemp[i];
		d2ldc2[i] = adi_gam * d2ldtemp[i];

		pres[i] = sigma[i] * temp[i];
		dldpres[i] = dlds[i] + dldtemp[i];
		d2ldpres[i] = d2lds[i] + d2ldtemp[i];

	}

#ifdef INPUTMASS
	sigma0 = Mdisk/calc_total_disk_mass();
	for(i=0;i<N;i++) {
		sigma[i] *= sigma0;
	}
#else
	Mdisk = calc_total_disk_mass();
#endif

	}



double calc_total_disk_mass(void) {
	int i;
	double res;
	res=0;
	for(i=0;i<N;i++) {
		res += weights[i]*sigma[i]*r[i]*r[i];
	}
	return res*2*M_PI;
}

void init_weights(void) {

	int i;

#ifdef CONSTWEIGHTS
	for(i=0;i<N;i++) weights[i] = dlr;
#endif
#ifdef COMPSIMPS
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


#endif

#ifdef COMPTRAPZ

	weights[0] = .5*dlr;
	weights[N-1] = .5*dlr;

	for(i=1;i<N-1;i++) {
		weights[i] = dlr;
	}

#endif



	return;

}
