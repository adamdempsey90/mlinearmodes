#include "eigen.h"




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

