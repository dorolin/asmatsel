#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>

using namespace std;

// macro definitions
#define VERSION "1.00 -- Mar 02 2017"
#define NDEMES 2 // number of demes (don't change)
#define MAXTRIES 100000 // max number of tries for some while loops

extern gsl_rng * r;


// data structures

struct demecont{
  gsl_matrix * adults; // genotype data all inds (L x (N+N))
  gsl_matrix * chpnt; // genotype data all inds checkpoint (L x (N+N))
  gsl_matrix * pool1; // genotype data mating pool 1 (L x (N+N))
  gsl_matrix * pool2; // genotype data mating pool 1 (L x (N+N))
  int N1; // individuals in pool1
  int N2; // individuals in pool2
};


struct parameters{
  int N; // total carrying capacity
  int iL; // number of loci in infile
  int uL; // number of unlinked loci
  int Zg; // calculate Zg
  int bins; // binsize
  int L; // total number of loci
  gsl_vector * d; // genetic map distances (iL)
  double mL; // map length
  int apos; // position of assortative mating locus
  gsl_vector_int * spos; // positions of loci under divergent selection
  int nsel; // number of loci under divergent selection
  int phases; // number of phases
  double sel; // selection coefficient
  double ass; // assortative mating strength
  double mig; // migration rate
  double mu; // mutation rate
  int mutone;
  double mutfrac;
};





// function declarations
void usage(char * name);
void getdata(string file, parameters * param);
void cleargrid(demecont (&grid)[NDEMES]);
void simulate(demecont (&grid)[NDEMES], parameters * param, int sub);
void copygrid(demecont (&grid)[NDEMES], int dir);
void changephase(demecont (&grid)[NDEMES], parameters * param, int sub);
void makegamete(gsl_vector * gamete, gsl_matrix * deme, int ind,
		parameters * param, gsl_matrix * tmpind);
void writeOut(demecont (&grid)[NDEMES], parameters * param,
	      int rep, int gen, int sub, string file);
void calcAf(gsl_vector * pvec, gsl_matrix * deme, parameters * param);
void calcFst(gsl_vector * fstvec, gsl_vector * p0vec, gsl_vector * p1vec,
	     parameters * param);
void calcZg(gsl_vector * zgvec, gsl_vector * p0vec, gsl_vector * p1vec,
	     parameters * param);
