#include "eigen.h"



double Dij(int i, int j) {
// 	if (i == (j+1)) {
// 		return -.5;
// 	}
// 	else if (i==(j-1)) {
// 		return .5;
// 	}
// 	else {
// 		return 0;
// 	}
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

// 	if (i == (j+1)) {
// 		return 1.0;
// 	}
// 	else if (i==(j-1)) {
// 		return 1.0;
// 	}
// 	else if (i==j) {
//
// 		return -2.0;
//
// 	}
// 	else {
// 		return 0;
// 	}

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




void init_derivatives(void) {
	int i,j,k,l,mindx;
	double dfac = 1.0/dlr;
	double d2fac = dfac*dfac;
	double d1, d2;
	int idn,id3;
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			d1 = Dij(i,j)*dfac;
			d2 = D2ij(i,j)*d2fac;
			if (i==j) idn = 1;
			else idn = 0;
			for(k=0;k<NF;k++) {
				for(l=0;l<NF;l++) {
					if (k==l) id3 = 1;
					else id3 = 0;
					mindx = getindex4(i,j,k,l,NF,N);
					D[mindx] = d1*id3;
					D2[mindx] = d2*id3;
					Identity[mindx] = idn*id3;


				}
			}
		}
	}
	return;

}
