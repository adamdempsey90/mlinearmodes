#include "eigen.h"



void lagrangian_pressure_bc_inner(double complex *mat, double complex *bcmat) {
	int ind,k,l;
/* Set up zero Lagrangian Pressure B.C at inner boundary */
	int i = 0;
	int j = 0;

	for(k=0;k<NF;k++) {
		for(j=0;j<N;j++) {
			for(l=0;l<NF;l++) {
				ind = getindex4(i,j,k,l,NF,N);
				mat[ind] = 0;
				if (j==0 && k==l) {
					bcmat[ind] = 1;
				}
				else {
					bcmat[ind] = 0;
				}
			}
		}
	}


	// for(ind=getindex4(i,0,0,0,NF,N) ; ind < getindex4(i,N-1,NF-1,NF-1,NF,N);ind++) {
	// 	mat[ind] = 0;
	// 	bcmat[ind] = 0;
	// }

/* Zero first NF rows */
	// for(ind=0; ind < NF*N*NF;ind++) {
	// 	mat[ind] = 0;
	// 	bcmat[ind] = 0;
	// }
/* Add b.c condition */

	j=0;
#ifdef BAROTROPIC
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	bcmat[getindex4(i,j,0,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*dlds[i]/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
//	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,2,0,NF,N)];
//	mat[getindex4(i,j,0,2,NF,N)] = mat[getindex4(i,j,2,2,NF,N)];

#endif

#ifdef ISOTHERMAL
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	bcmat[getindex4(i,j,0,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*(dldtemp[i] + dlds[i])/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,2,0,NF,N)];
	mat[getindex4(i,j,0,2,NF,N)] = mat[getindex4(i,j,2,2,NF,N)];
#endif

#ifdef COOLING
	bcmat[getindex4(i,j,3,3,NF,N)] = 1;
	bcmat[getindex4(i,j,0,3,NF,N)] = 1;
	mat[getindex4(i,j,3,0,NF,N)] = pres[i]*dldpres[i]/(I*mval*r[i]);
	mat[getindex4(i,j,3,3,NF,N)] = omega[i];
	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,3,0,NF,N)];
	mat[getindex4(i,j,0,3,NF,N)] = mat[getindex4(i,j,3,3,NF,N)];
#endif

	return;
}

void lagrangian_pressure_bc_outer(double complex *mat, double complex *bcmat) {
	int ind,k,l;
/* Set up zero Lagrangian Pressure B.C at outer boundary */

	int i = N-1;
	int j = N-1;

	for(k=0;k<NF;k++) {
		for(j=0;j<N;j++) {
			for(l=0;l<NF;l++) {
				ind = getindex4(i,j,k,l,NF,N);
				mat[ind] = 0;
				if (j==N-1 && k==l) {
					bcmat[ind] = 1;
				}
				else {
					bcmat[ind] = 0;
				}
			}
		}
	}
/* Zero last NF rows */
	// for(ind=getindex4(i,0,0,0,NF,N) ; ind < getindex4(i,N-1,NF-1,NF-1,NF,N);ind++) {
	// 	mat[ind] = 0;
	// 	bcmat[ind] = 0;
	// }
/* Add b.c condition */

	j=N-1;
#ifdef BAROTROPIC
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	bcmat[getindex4(i,j,0,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*dlds[i]/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,2,0,NF,N)];
	mat[getindex4(i,j,0,2,NF,N)] = mat[getindex4(i,j,2,2,NF,N)];

#endif

#ifdef ISOTHERMAL
	bcmat[getindex4(i,j,2,2,NF,N)] = 1;
	bcmat[getindex4(i,j,0,2,NF,N)] = 1;
	mat[getindex4(i,j,2,0,NF,N)] = sigma[i]*(dldtemp[i] + dlds[i])/(I*mval*r[i]);
	mat[getindex4(i,j,2,2,NF,N)] = omega[i];
	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,2,0,NF,N)];
	mat[getindex4(i,j,0,2,NF,N)] = mat[getindex4(i,j,2,2,NF,N)];
#endif

#ifdef COOLING
	bcmat[getindex4(i,j,3,3,NF,N)] = 1;
	bcmat[getindex4(i,j,0,3,NF,N)] = 1;
	mat[getindex4(i,j,3,0,NF,N)] = pres[i]*dldpres[i]/(I*mval*r[i]);
	mat[getindex4(i,j,3,3,NF,N)] = omega[i];
	mat[getindex4(i,j,0,0,NF,N)] = mat[getindex4(i,j,3,0,NF,N)];
	mat[getindex4(i,j,0,3,NF,N)] = mat[getindex4(i,j,3,3,NF,N)];
#endif
	return;
}

void zero_bc_outer(double complex *mat,double complex *bcmat) {
	int ind,k,l;
/* Set up zero u,v,s,p  B.C at outer boundary
 * Can be thought of as a no slip wall b.c */

	int i = N-1;
	int j = N-1;
/* Zero last NF rows */
for(j=0;j<N;j++) {
	for(k=0;k<NF;k++) {
		for(l=0;l<NF;l++) {
			ind = getindex4(i,j,k,l,NF,N);
			mat[ind] = 0;
			bcmat[ind] = 0;
		}
	}
}
/* Add b.c condition */
	for(k=0;k<NF;k++) {
		mat[getindex4(i,j,k,k,NF,N)] = 1;
	}
	return;
}

void zero_bc_inner(double complex *mat,double complex *bcmat) {
	int ind,k;
/* Set up zero u,v,s,p  B.C at inner boundary
 * Can be thought of as a no slip wall b.c */

	int i = 0;
	int j = 0;
/* Zero last NF rows */
	for(ind=getindex4(i,0,0,0,NF,N) ; ind < getindex4(i,N-1,NF-1,NF-1,NF,N);ind++) {
		mat[ind] = 0;
		bcmat[ind] = 0;
	}
/* Add b.c condition */
	for(k=0;k<NF;k++) {
		mat[getindex4(i,j,k,k,NF,N)] = 1;
	}
	return;
}
void reflecting_bc(double complex *mat, double complex *bcmat) {
	int j,k,l;
	for(j=0;j<N;j++) {
		for(k=0;k<NF;k++) {
			for(l=0;l<NF;l++) {
				mat[getindex4(0,j,k,l,NF,N)] = 0;
				mat[getindex4(N-1,j,k,l,NF,N)] = 0;
				bcmat[getindex4(0,j,k,l,NF,N)] = 0;
				bcmat[getindex4(N-1,j,k,l,NF,N)] = 0;
			}
		}
	}

	mat[getindex4(0,1,0,0,NF,N)] = 1;
	mat[getindex4(N-1,N-2,0,0,NF,N)] = 1;
	mat[getindex4(0,0,0,0,NF,N)] = 1;
	mat[getindex4(N-1,N-1,0,0,NF,N)] = 1;
	for(k=1;k<NF;k++) {
		mat[getindex4(0,1,k,k,NF,N)] = -1;
		mat[getindex4(N-1,N-2,k,k,NF,N)] = -1;
		mat[getindex4(0,0,k,k,NF,N)] = -1;
		mat[getindex4(N-1,N-1,k,k,NF,N)] = -1;
	}
	return;
}

void set_bc(double complex *mat, double complex *bcmat) {


#ifndef NOBC
#ifndef NOPRESSURE

#ifdef ZEROBCIN
	zero_bc_inner(mat,bcmat);
#endif

#ifdef ZEROBCOUT
	zero_bc_outer(mat,bcmat);
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

#ifdef REFLECTBC
	reflecting_bc(mat,bcmat);
#endif

#endif
#endif
	return;
}
