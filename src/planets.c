#include "eigen.h"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>

#ifdef PLANETS



double complex *matpp, *matdp, *matpd;


void read_planets_file(void);
void calc_planet_potentials(void);
double planetpot0(double x, double a, double eps);
double planetpot1(double x, double a, double eps);
double regularize_radius(double x, int *indx);
void output_planet_summary(void);
void output_planet_potentials(void);


void init_planets(void) {
	int i;
	Planets = (Planet *)malloc(sizeof(Planet)*NP);

	for(i=0;i<NP;i++) {
		Planets[i].pot0 = (double complex *)malloc(sizeof(double complex)*(N+NP));
		Planets[i].pot1 = (double complex *)malloc(sizeof(double complex)*(N+NP));
	}

	matpp = (double complex *)malloc(sizeof(double complex)*NP*NP);
	matdp = (double complex *)malloc(sizeof(double complex)*N*NP);
	matpd = (double complex *)malloc(sizeof(double complex)*N*NP);


	read_planets_file();
	return;


}

void free_planets(void) {
	int i;

	for(i=0;i<NP;i++) {
		free(Planets[i].pot0);
		free(Planets[i].pot1);
	}
	free(Planets);

	free(matpp);
	free(matdp);
	free(matpd);

	return;
}



void read_planets_file(void) {
	double mp, a;
	int i,indx;
	FILE *f;
	f=fopen("planets.in","r");
	if (f==NULL) printf("Can't find Planets File!"); return;
	for(i=0;i<NP;i++) {
		if (!fscanf(f, "%lg %lg\n", &a,&mp))
          		break;
		Planets[i].mass = mp;
		Planets[i].position = regularize_radius(a,&indx);
		Planets[i].index = indx;

		Planets[i].hill = Planets[i].position*pow(Planets[i].mass/3.,1./3);
		Planets[i].InteriorPlanet = 0;
		Planets[i].ExteriorPlanet = 0;
		Planets[i].wp = 0;

		if (Planets[i].position < r[0]) {
			Planets[i].InteriorPlanet = 1;
		}

		if (Planets[i].position > r[N-1]) {
				Planets[i].ExteriorPlanet = 1;
		}


	}
	fclose(f);


	return;
}


void output_planet_summary(void) {
	int i;
	FILE *f;
	f = fopen("planet_summary.dat","w");

	printf("\n\nPlanet Summary \n");
	printf("PlanetNumber\tRadius\tMass\tHillRadius\twp\n");
	fprintf(f,"#PlanetNumber\tRadius\tMass\tHillRadius\twp\n");
	for(i=0;i<NP;i++) {
		printf("%d\t%lg\t%.2e\t%lg\t%lg\n",i,Planets[i].position,
							Planets[i].mass,Planets[i].hill,Planets[i].wp);
		fprintf(f,"%d\t%.12lg\t%.12e\t%.12lg\t%.12lg\n",i,Planets[i].position,
												Planets[i].mass,Planets[i].hill,Planets[i].wp);
	}

	printf("\n\n");
	fclose(f);
	return;


}

void output_planet_potentials(void) {
	int i,j;
	FILE *f;

	f = fopen("planet_pots.dat","w");

	for(j=0;j<(N+NP);j++) {
		for(i=0;i<NP;i++) {
			fprintf(f,"%.12lg\t%.12lg\t%.12lg\t%.12lg\t",
						creal(Planets[i].pot0[j]),cimag(Planets[i].pot0[j]),
						creal(Planets[i].pot1[j]),cimag(Planets[i].pot1[j]));
		}
		fprintf(f,"\n");
	}
	fclose(f);
	return;

}


void add_planets(double complex *mat, double complex *bcmat) {
	int i,j,indx;

/* Add Planet-Planet Matrix */
	for(i=0;i<NP;i++) {
		for(j=0;j<NP;j++) {
			indx = j+N+ nrows*(i+N);
			mat[indx] = matpp[j + NP*i];
			bcmat[indx] = 0;
			if (i==j) {
					bcmat[indx] = 1;
			}
		}
	}

/* Add Planet-Disk Matrix */
	for(i=0;i<N;i++) {
		for(j=0;j<NP;j++) {
			indx = j+N + nrows*i;
			mat[indx] = matpd[j+NP*i];
			bcmat[indx] = 0;
		}
	}
/* Add Disk-Planet Matrix */
		for(i=0;i<NP;i++) {
			for(j=0;j<N;j++) {
				indx = j + nrows*(i+N);
				mat[indx] = matdp[j+N*i];
				bcmat[indx] = 0;
			}
		}


	return;



}
void calc_planet_matrices(void) {
	int i,j,indx;
	double a, norm;


	printf("\tCalculating Planet Potentials...\n");
	calc_planet_potentials();
	printf("Filling Planet Matrices...\n");
/* Planet-Planet Matrix */

	for(i=0;i<NP;i++) {
		a = Planets[i].position;
		norm = 2*omk_func(a)*a*a*a;
		for(j=0;j<NP;j++) {
			indx = j + NP*i;
			matpp[indx] = (Planets[j].pot1[i+N])/norm;
			if (i==j) {
				matpp[indx] += Planets[i].wp;
			}

		}
	}
/* Planet-Disk Matrix */

	for(i=0;i<N;i++) {
		norm = 2*omk_func(r[i])*r[i]*r[i]*r[i];
		for(j=0;j<NP;j++) {
			indx = j+NP*i;
			matpd[indx] = (Planets[j].pot1[i])/norm;
		}
	}

/* Disk-Planet Matrix */

	for(i=0;i<NP;i++) {
		a = Planets[i].position;
		norm = 2*omk_func(a)*a*a*a;
		for(j=0;j<N;j++) {
			indx = j+N*i;
			norm = 2*omk_func(a);
			matdp[indx] =  sg_integrand(a, r[j],eps)*weights[j]*r[j]*sigma[j]/norm;

		}
	}

	return;
}

void calc_planet_potentials(void)  {
	int i,j;
	double a;

	for(i=0;i<NP;i++) {

/* Potential of Planet i in the disk at r[j] */
		for(j=0;j<N;j++) {
			Planets[i].pot0[j] = Planets[i].mass* planetpot0(r[j],Planets[i].position,
																Planets[i].hill);
			Planets[i].pot1[j] = Planets[i].mass*planetpot1(r[j],Planets[i].position,
																	Planets[i].hill);
		}

/* Potential of Planet i at Planet j */
		for(j=0;j<NP;j++) {
			if (i!=j) {
				Planets[i].pot0[j+N] = Planets[i].mass
										*planetpot0(Planets[j].position,Planets[i].position,
													Planets[i].hill);
				Planets[i].pot1[j+N] = Planets[i].mass
										*planetpot1(Planets[j].position,Planets[i].position,
													Planets[i].hill);
			}
			else {
				Planets[i].pot0[j+N] = 0;
				Planets[i].pot1[j+N] = 0;
			}
		}

	}

/* Calculate the precession frequency */
	for(i=0;i<NP;i++) {
		a = Planets[i].position;
		for(j=0;j<NP;j++) {
			Planets[i].wp += Planets[j].pot0[i+N]/(2 * omk_func(a)*a*a);
		}
	}

	output_planet_potentials();
	return;
}


double planetpot0(double x, double a, double eps) {
	double res;
	double kval, ee, ek;

	double eps2,a2,x2,x4;

	eps2 = eps*eps;
	a2 = a*a;
	x2 = x*x;
	x4 = x2*x2;


	kval = sqrt(4*a*x/(eps2+(a+x)*(a+x)));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	res = ee*(a2 + eps2+x2)*(a2+eps2+2*eps*x-x2)*(a2+eps2-2*eps*x-x2);
	res -= ek*(eps2+(a-x)*(a-x))*((a2+eps2)*(a2+eps2)-2*a2*x2+x4);
	res /= (M_PI*pow(eps2 + (a-x)*(a-x),2)*pow(eps2 + (a+x)*(a+x),1.5));

	return res;
}

double planetpot1(double x, double a, double eps) {
	double res;
	double kval,ee, ek;

	double eps2,eps4,a2,a4,x2,x4,x6,x8;

	eps2 = eps*eps;
	eps4 = eps2*eps2;
	a2 = a*a;
	a4 = a2*a2;
	x2 = x*x;
	x4 = x2*x2;
	x6 = x4*x2;
	x8 = x6*x2;




	kval = sqrt(4*a*x/(eps2+(a+x)*(a+x)));
	ek =  gsl_sf_ellint_Kcomp(kval, GSL_PREC_DOUBLE);
	ee = gsl_sf_ellint_Ecomp(kval,GSL_PREC_DOUBLE);

	res = -ee*(2*pow((a2+eps2),4) - (6*a2-7*eps2)*(a2+eps2)*(a2+eps2)*x2
		  +2*(4*a4-9*a2*eps2+5*eps4)*x4 +(-6*a2+7*eps2)*x6 + 2*x8);
	res += ek*(eps2 +(a-x)*(a-x))*(a2+eps2+x2)*(2*(a2+eps2)*(a2+eps2)
			+(-4*a2+eps2)*x2+2*x4);

	res /= (a*M_PI*pow(eps2+(a-x)*(a-x),2)*pow(eps2 + (a+x)*(a+x),1.5));

	return res;
}



double regularize_radius(double x, int *indx) {
	int i;

	if (x <= r[0]) {
		*indx = -2;
		return x;
	}

	if (x >= r[N-1]) {
		*indx=-1;
		return x;
	}

	for(i=1;i<(N-1);i++) {

		if (x <= r[i]) {
			if ( (x - r[i-1]) <= (r[i] - x)) {
				*indx = i-1;
				return r[i-1];
			}
			else {
				*indx = i;
				return r[i];
			}
		}
	}
	printf("\n\n\n\n Couldn't put the planet on the grid! \n\n\n\n");
	return -1;
}
#endif
