

void read_arguments(int argc, char *argv[]) {


  	N = atoi(argv[1]);
  	nrows = N; ncols = N;

  	ri = atof(argv[2]);

  	ro = atof(argv[3]);

  #ifdef INPUTMASS
  	Mdisk = atof(argv[4]);
  	sigma0 = 1;
  #else
  	sigma0 = atof(argv[4]);
  	Mdisk = 1;
  #endif

  	eps = atof(argv[5]);
  	h0 = atof(argv[6]);
  	sigma_index = atof(argv[7]);
  	flare_index = atof(argv[8]);
  	temp_index = 2*flare_index - 1;

  	alpha_s = atof(argv[9]);
  	alpha_b = atof(argv[10]);




  #if  defined(COOLING) || defined(ADIABATIC)
  	adi_gam = atof(argv[12]);
  	beta_cool = atof(argv[13]);
  #else
  	adi_gam = 1;
  #endif

  #ifdef ADIABATIC
  	beta_cool = 0;
  #endif

  	tol = atof(argv[14]);

  	dlr = (log(ro) - log(ri))/((float) N);

  #ifdef PLANETS
  	NP = atoi(argv[15]);
  	nrows += NP;
  	ncols += NP;
  #else
  	NP = 0;
  #endif

#ifdef OPENMP
  nprocs = atoi(argv[11]);
#else
  nprocs = 0;
#endif

    return;
}


void read_input_file(void) {
  char garbage[100];
	char *gchar;
  int read_res;
  FILE *f;

  f = fopen("params.in","r");

  if (f==NULL) printf("\n\nERROR Can't Find Input File!\n\n");

	gchar=fgets(garbage,sizeof(garbage),f);	// Input Parameters
	read_res=fscanf(f,"Nr = %d \n",&N);
  read_res=fscanf(f,"ri = %lg \n",&ri);
  read_res=fscanf(f,"ro = %lg \n",&ro);
  read_res=fscanf(f,"Mdisk = %lg \n",&Mdisk);
  read_res=fscanf(f,"rs = %lg \n",&eps);
  read_res=fscanf(f,"h0 = %lg \n",&h);
  read_res=fscanf(f,"sig_ind = %lg \n",&sigma_index);
  read_res=fscanf(f,"flare_ind = %lg \n",&flare_index);
  read_res=fscanf(f,"alpha_s = %lg \n",&alpha_s);
  read_res=fscanf(f,"alpha_b = %lg \n",&alpha_b);
  read_res=fscanf(f,"np = %d \n",&np);
  read_res=fscanf(f,"gam = %lg \n",&adi_gam);
  read_res=fscanf(f,"beta = %lg \n",&beta_cool);
  read_res=fscanf(f,"tol = %lg \n",&tol);
  read_res=fscanf(f,"Nplanets = %d \n",&NP);


  fclose(f);

  set_secondary_inputs();


  return;




}


void set_secondary_inputs(void) {

  nrows = N; ncols = N;
#ifdef PLANETS
  nrows += NP;
  ncols += NP;
#else
  NP = 0;
#endif


#ifdef INPUTMASS
  sigma0 = 1;
#else
  sigma0 = Mdisk;
  Mdisk = 1;
#endif

  temp_index = 2*flare_index - 1;

#if !defined(COOLING) && !defined(ADIABATIC)
  adi_gam = 1;
#endif

  return;
}
