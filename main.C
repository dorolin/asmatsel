// simulate the effects of divergent selection and assortative mating
// on patterns of genome differentiation and linkage disequilibrium

// compile: try
// g++ -Wall -O2 -o asmatsel main.C func.C -lm -lgsl -lgslcblas

// run:
// ./asmatsel -i ./map2c.txt -l 700 -x 1 -z 1 -b 10 -u 0 -o ./test2.out


// >>>> TODO:
//      - instead exiting when no suitable haplotypes for assortative
//        mating locus, continue simulation without it (also add
//        possibility to introduce assortative mating without having
//        selected loci at all)

// ----------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <getopt.h>

#include "asmatsel.H"

using namespace std;

gsl_rng * r;  // global state variable for random number generator

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {


  int ch = 0, x, X;
  int i, k, t, lastT;
  double a0, a1;
  gsl_vector * alleles;

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0;
  int ext_seed = -1;

  string infile = "undefined";
  string outfile = "simul.out";
  ofstream file;


  // fixed parameters
  int burn = 20000; // # generations in burn-in phase
  int subphase = 5000; // # generations between phases
  int ngen = 50000; // total # generations past burn-in
  int samp1 = 100; // sampling interval during phases
  int samp2 = 1000; // sampling interval past phases

  // default values
  struct parameters param;
  param.uL = 0;
  param.Zg = 1;
  param.bins = 10;
  param.N = 500;
  param.sel = 0.5;
  param.ass = 0.6;
  param.mig = 0.1;
  param.mu = 0.0001;
  param.apos = -1;
  param.nsel = 0;
  param.phases = 0;
  X = 1;

  // additional parameters controlling assortative mating evolution
  param.mutone = 0; // mutate just one allele (0, no; 1, yes)
  param.mutfrac = 1.0; // prob that A allele is on haplotype with derived alleles under divergent selection if mutone = 0



  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }

  while ((ch = getopt(argc, argv, "i:o:l:u:N:s:a:m:n:x:z:b:d:")) != -1){
    switch(ch){
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'l':
      param.iL = atoi(optarg);
      break;
    case 'u':
      param.uL = atoi(optarg);
      break;
    case 'N':
      param.N = atoi(optarg);
      break;
    case 's':
      param.sel = atof(optarg);
      break;
    case 'a':
      param.ass = atof(optarg);
      break;
    case 'm':
      param.mig = atof(optarg);
      break;
    case 'n':
      param.mu = atof(optarg);
      break;
    case 'x':
      X = atoi(optarg);
      break;
    case 'z':
      param.Zg = atoi(optarg);
      break;
    case 'b':
      param.bins = atoi(optarg);
      break;
    case 'd':
      ext_seed = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }



  if (infile.compare("undefined") == 0){
    cerr << "Error: infile must be specified\n" << endl;
    exit(1);
  }
  if(param.iL == 0){
    cerr << "Error: number of loci must be specified\n" << endl;
    exit(1);
  }
  if(param.bins < 1){
    cerr << "Error: need at least one locus per bin" << endl;
    exit(1);
  }
  if(param.Zg != 0 && param.bins < 2){
    cerr << "Error: need > 1 locus per bin to compute Zg" << endl;
    exit(1);
  }
  if(param.Zg != 0 && param.uL > 0){
    cerr << "Warning: number of unlinked loci was set to zero to calculate Zg" << endl;
    param.uL = 0;
  }


  // sum # loci in map and # unlinked loci
  param.L = param.iL + param.uL;

  // overwrite old outfile if any
  file.open(outfile.c_str());
  file.close();


  // set up gsl random number generation
  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);
  (ext_seed > -1) ? srand(ext_seed) : srand(time(NULL));
  rng_seed = rand();
  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
                               was seeded with result of time(NULL) */



  // allocate memory
  demecont grid[NDEMES];
  if(grid == NULL){
    cerr << "Error in allocating memory for grid (" << sizeof(struct demecont) << ", " <<  NDEMES << ")\n" << endl;
  }

  for (k=0; k<NDEMES; k++){
    grid[k].adults = gsl_matrix_calloc(param.L,param.N * 2);
    grid[k].chpnt = gsl_matrix_calloc(param.L,param.N * 2);
    grid[k].pool1 = gsl_matrix_calloc(param.L,param.N * 2);
    grid[k].pool2 = gsl_matrix_calloc(param.L,param.N * 2);
  }

  param.d = gsl_vector_calloc(param.iL);
  alleles = gsl_vector_alloc(param.N * 2);


  // read infile
  cerr << "\nInput file: " << infile << endl;
  getdata(infile, &param);


  // print settings
  cerr << "Deme size: " << param.N << endl;
  cerr << "Migration rate: " << param.mig << endl;
  cerr << "Mutation rate: " << param.mu << endl;
  if(param.nsel > 0) cerr << "Selection coefficient: " << param.sel << endl;
  if(param.apos >= 0) cerr << "Assortment strength: " << param.ass << endl;
  cerr << "Number of loci in map: " << param.iL << endl;
  cerr << "Number of unlinked loci: " << param.uL << endl;
  cerr << "Number of loci per bin: " << param.bins << endl;
  cerr << "Map length: " << param.mL << endl;
  if(param.apos >= 0) cerr << "Assortative mating locus at position: " << gsl_vector_get(param.d,param.apos) << endl;
  if(param.nsel == 1){
    cerr << "Locus under divergent selection at position: ";
    cerr << gsl_vector_get(param.d,gsl_vector_int_get(param.spos,0)) << endl;
  }
  else if(param.nsel > 1){
    cerr << "Loci under divergent selection at positions: ";
    for(i=0; i<param.nsel; i++){
      cerr << gsl_vector_get(param.d,gsl_vector_int_get(param.spos,i)) << " ";
    }
    cerr << endl;
  }
  // cerr << "\nPositions of loci in map: " << endl;
  // for(i=0; i<param.iL; i++){
  //   cerr << gsl_vector_get(param.d,i) << " ";
  // }
  // cerr << endl;


  // simulate
  x = 0;
  while(x < X){
    cerr << "\nReplicate " << x+1 << endl;
    cleargrid(grid);
    // burn-in
    for(t=0; t<burn; t++){
      simulate(grid,&param,-1);
    }
    // subphases
    for(i=0; i<param.phases; i++){
      cerr << "Phase " << i+1 << endl;
      changephase(grid,&param,i+1);
      copygrid(grid, 0); // first checkpoint
      lastT = burn+subphase*i;
      for(t=burn+subphase*i; t<burn+subphase*(1+i) && t<burn+ngen; t++){
    	if (t % samp1 == 0){ // writing step
	  //if(t > lastT && t <= lastT + samp1){ // checkpoint
	  a0 = 0;
	  a1 = 0;
	  if(param.nsel > 0 && i<param.nsel){ // selection phase
	    // verify that mutation still polymorphic
	    for (k=0; k<NDEMES; k++){
	      gsl_matrix_get_row(alleles,grid[k].adults,
				 gsl_vector_int_get(param.spos,i));
	      a0 += gsl_vector_max(alleles);
	      a1 -= gsl_vector_min(alleles);
	    }
	  }
	  else if(param.apos >= 0 && i==param.nsel){ // assort mat phase
	    // verify that mutation still polymorphic
	    for (k=0; k<NDEMES; k++){
	      gsl_matrix_get_row(alleles,grid[k].adults,param.apos);
	      a0 += gsl_vector_max(alleles);
	      a1 -= gsl_vector_min(alleles);
	    }
	  }
	  if(a0==0 || a1<0){ // go back to last checkpoint
	    copygrid(grid, 1);
	    t = lastT;
	  }
	  else{
	    writeOut(grid,&param,x,t,i,outfile); // write outfile
	    copygrid(grid, 0); // new checkpoint
	    lastT = t;
	  }
	  //}
	  //else{
	  //writeOut(grid,&param,x,t,i,outfile);
	  //}
    	}
    	simulate(grid,&param,i+1);
      }
    }
    // final phase
    cerr << "Final phase" << endl;
    for(t=burn+subphase*i; t<burn+ngen; t++){
      if (t % samp2 == 0){ // writing step
    	writeOut(grid,&param,x,t,i,outfile);
      }
      simulate(grid,&param,i);
    }
    writeOut(grid,&param,x,t,i,outfile);
    x++;
  }


  cerr << "\nOutput written to " << outfile << endl;

  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " <<
    (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;

  return 0;
}

