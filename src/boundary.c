#include "eigen.h"



void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {
	int ind;
/* Set up zero Lagrangian Pressure B.C at inner boundary */


/* Zero first NF rows */
	for(ind=0; ind < NF*N*NF;ind++) {
		mat[ind] = 0;
		bcmat[ind] = 0;
	}
/* Add b.c condition */

	int i = 0;
	int j = 0;
#ifdef BAROTROPIC
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	mat[getindex4(i,j,0,0,NF,N)] = sigma[i]*dlds[i]/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];

#endif

#ifdef ISOTHERMAL
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	mat[getindex4(i,j,0,0,NF,N)] = sigma[i]*temp[i]*(dldtemp[i] + dlds[i])/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
#endif

#ifdef COOLING
	bcmat[getindex4(i,j,3,3,NF,N)] = 1;
	mat[getindex4(i,j,0,0,NF,N)] = pres[i]*dldpres[i]/(I*mval*r[i]);
	mat[getindex4(i,j,3,3,NF,N)] = omega[i];
#endif

	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {
	int ind;
/* Set up zero Lagrangian Pressure B.C at outer boundary */

	int i = N-1;
	int j = N-1;
/* Zero last NF rows */
	for(ind=getindex4(i,0,0,0,NF,N) ; ind < getindex4(i,j,NF-1,NF-1,NF,N);ind++) {
		mat[ind] = 0;
		bcmat[ind] = 0;
	}
/* Add b.c condition */


#ifdef BAROTROPIC
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*dlds[i]/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];

#endif

#ifdef ISOTHERMAL
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*temp[i]*(dldtemp[i] + dlds[i])/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
#endif

#ifdef COOLING
	bcmat[getindex4(i,j,3,3,NF,N)] = 1;
	mat[getindex4(i,j,3,0,NF,N)] = pres[i]*dldpres[i]/(I*mval*r[i]);
	mat[getindex4(i,j,3,3,NF,N)] = omega[i];
#endif

	return;
}

void set_bc(double complex *mat, double complex *bcmat) {


#ifndef NOBC
#ifndef NOPRESSURE

#ifdef ZEROBCIN
	zero_e_bc_inner(mat,bcmat);
#endif

#ifdef ZEROBCOUT
	zero_e_bc_outer(mat,bcmat);
#endif

#ifdef PRESBCIN
	lagrangian_pressure_bc_inner(mat, bcmat);
#endif

#ifdef PRESBCOUT
	lagrangian_pressure_bc_outer(mat, bcmat);
#endif

#ifdef GRADBCIN
	user_gradient_bc_outer(mat,bcmat,.5);
#endif


#endif
#endif
	return;
}
