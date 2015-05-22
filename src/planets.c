#include "eigen.h"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>



typedef struct Planet {

	double mass, position, hill;
	double *pot0, *pot1;
	int InteriorPlanet, ExteriorPlanet

} Planet;

Planet *Planets;

void init_planets(void) {
	int i;
	Planets = (Planet *)malloc(sizeof(Planet)*NP);
	
	for(i=0;i<NP;i++) {
		Planets[i]->pot0 = (double *)malloc(sizeof(double)*N);
		Planet[i]->pot1 = (double *)malloc(sizeof(double)*N);
	}
	
	read_planets_file();
	init_planet_potentials();
	
	return;


}

void free_planets(void) {
	int ;
	
	for(i=0;i<NP;i++) {
		free(Planets[i]->pot0);
		free(Planets[i]->pot1);
	}
	free(Planets);

	return;
}



void read_planets_file(void) {
	double mp, a;
	
	FILE *f;
	f=fopen("planets.in","r");

	for(i=0;i<NP;i++) {
		if (!fscanf(f, "%lg %lg\n", &mp,&a)
          		break;
		Planets[i]->mass = mp;
		Planets[i]->position = regularize_radius(a);
		Planets[i]->hill = Planets[i]->position*pow(Planets[i]->mass/3.,1./3);
		Planets[i]->InteriorPlanet = 0;
		Planets[i]->ExteriorPlanet = 0;
		
		if (Planets[i]->position < r[0]) {
			Planets[i]->InteriorPlanet = 1;
		}
		
		if (Planets[i]->position > r[N-1]) {
				Planets[i]->ExteriorPlanet = 1;
		}
		

	}
	fclose(f);


	return;
}

void init_planet_potentials(void)  {
	int i,j;
	
	
	for(i=0;i<NP;i++) {
		for(j=0;j<N;j++) {
			Planets[i]->pot0[j] = Planets[i]->mass*planetpot0(r[j],Planets[i]->position,Planets[i]->hill);
			Planets[i]->pot1[j] = Planets[i]->mass*planetpot1(r[j],Planets[i]->position,Planets[i]->hill);



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
	res /= (M_Pi*pow(eps2 + (a-x)*(a-x),2)*pow(eps2 + (a+x)*(a+x),1.5));

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
		indx = NULL;
		return x;
	}
	
	if (x >= r[N-1]) {
		indx=NULL;
		return x;
	}
	
	for(i=1;i<(N-1);i++) {
		
		if (x <= r[i]) {
			ld = x - r[i-1];
			rd = r[i] - x;
			if ( (x - r[i-1]) <= (r[i] - x)) {
				indx = i-1;
				return r[i-1];
			}
			else {
				indx = i;
				return r[i];
			}
		}
	}
	printf("\n\n\n\n Couldn't put the planet on the grid! \n\n\n\n");
	return -1;
}

