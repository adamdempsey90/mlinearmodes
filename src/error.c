#include "eigen.h"
#include <gsl/gsl_errno.h>

void bad_exit(void) {
  printf("Aborting..\n");
  free_globals();
#ifdef PLANETS
  free_planets();
#endif
  exit (-1);
}


void gsl_integration_error_handler(int status,char *reason) {
  fprintf(stderr, reason);
	if (status == GSL_SUCCESS){
		fprintf(stderr,"Success\n");
		bad_exit();
	}
	else if (status == GSL_FAILURE){
		fprintf(stderr,"Failure\n");
		bad_exit();
	}
	else if (status == GSL_CONTINUE){
		fprintf(stderr,"iteration has not converged\n");
		bad_exit();
	}
	else if (status == GSL_EDOM){
		fprintf(stderr,"input domain error, e.g sqrt(-1)\n");
		bad_exit();
	}
	else if (status == GSL_ERANGE){
		fprintf(stderr,"output range error, e.g. exp(1e100)\n");
		bad_exit();
	}
	else if (status == GSL_EFAULT){
		fprintf(stderr,"invalid pointer\n");
		bad_exit();
	}
	else if (status == GSL_EINVAL){
		fprintf(stderr,"invalid argument supplied by user\n");
		bad_exit();
	}
	else if (status == GSL_EFAILED){
		fprintf(stderr,"generic failure\n");
		bad_exit();
	}
	else if (status == GSL_EFACTOR){
		fprintf(stderr,"factorization failed\n");
		bad_exit();
	}
	else if (status == GSL_ESANITY){
		fprintf(stderr,"sanity check failed - shouldn't happen\n");
		bad_exit();
	}
	else if (status == GSL_ENOMEM){
		fprintf(stderr,"malloc failed\n");
		bad_exit();
	}
	else if (status == GSL_EBADFUNC){
		fprintf(stderr,"problem with user-supplied function\n");
		bad_exit();
	}
	else if (status == GSL_ERUNAWAY){
		fprintf(stderr,"iterative process is out of control\n");
		bad_exit();
	}
	else if (status == GSL_EMAXITER){
		fprintf(stderr,"exceeded max number of iterations\n");
		bad_exit();
	}
	else if (status == GSL_EZERODIV){
		fprintf(stderr,"tried to divide by zero\n");
		bad_exit();
	}
	else if (status == GSL_EBADTOL){
		fprintf(stderr,"user specified an invalid tolerance\n");
		bad_exit();
	}
	else if (status == GSL_ETOL){
		fprintf(stderr,"failed to reach the specified tolerance\n");
    fprintf(stderr, "Continuing..,\n");
//		bad_exit();
	}
	else if (status == GSL_EUNDRFLW){
		fprintf(stderr,"underflow\n");
		bad_exit();
	}
	else if (status == GSL_EOVRFLW){
		fprintf(stderr,"overflow\n");
		bad_exit();
	}
	else if (status == GSL_ELOSS){
		fprintf(stderr,"loss of accuracy\n");
		bad_exit();
	}
	else if (status == GSL_EROUND){
		fprintf(stderr,"failed because of roundoff error\n");
    fprintf(stderr,"Continuing anyway...\n");
//		bad_exit();
	}
	else if (status == GSL_EBADLEN){
		fprintf(stderr,"matrix, vector lengths are not conformant\n");
		bad_exit();
	}
	else if (status == GSL_ENOTSQR){
		fprintf(stderr,"matrix not square\n");
		bad_exit();
	}
	else if (status == GSL_ESING){
		fprintf(stderr,"apparent singularity detected\n");
		bad_exit();
	}
	else if (status == GSL_EDIVERGE){
		fprintf(stderr,"integral or series is divergent\n");
		bad_exit();
	}
	else if (status == GSL_EUNSUP){
		fprintf(stderr,"requested feature is not supported by the hardware\n");
		bad_exit();
	}
	else if (status == GSL_EUNIMPL){
		fprintf(stderr,"requested feature not (yet) implemented\n");
		bad_exit();
	}
	else if (status == GSL_ECACHE){
		fprintf(stderr,"cache table limit exceeded\n");
		bad_exit();
	}
	else if (status == GSL_ENOPROG){
		fprintf(stderr,"iteration is not making progress towards solution\n");
		bad_exit();
	}
	else if (status == GSL_ENOPROGJ){
		fprintf(stderr,"jacobian evaluations are not improving the solution\n");
		bad_exit();
	}
	else {
		fprintf(stderr,"GSL Error: Could not find error\n");
		bad_exit();
	}
	return;
}
