#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sort_vector.h>
#include <float.h>
#include <math.h>

#include "asmatsel.H"

using namespace std;


// -----------------
// usage
void usage(char * name){
  fprintf(stdout,"\n%s version %s\n\n", name, VERSION);
  fprintf(stdout, "Usage: balselabc [options]\n");
  fprintf(stdout, "-i     Infile\n");
  fprintf(stdout, "-o     Outfile\n");
  fprintf(stdout, "-l     Number of loci in infile\n");
  fprintf(stdout, "-u     Number of unlinked loci [0]\n");
  fprintf(stdout, "-N     Deme size [500]\n");
  fprintf(stdout, "-s     Selection coefficient [0.5]\n");
  fprintf(stdout, "-a     Strength assortative mating [0.6]\n");
  fprintf(stdout, "-m     Migration rate [0.1]\n");
  fprintf(stdout, "-n     Mutation rate [0.0001]\n");
  fprintf(stdout, "-x     Number of replicates [1]\n");
  fprintf(stdout, "-z     Calculate Zg [1] (forces -u 0 if -z 1)\n");
  fprintf(stdout, "-b     Number of loci per bin [10]\n");
  fprintf(stdout, "-d     Optionally specify external seed [int]\n");

  exit(1);
}


// -----------------
// getdata
//   Input file struct: ' ' delimited
//   No header
//   In rows: loci
//   1st column: map position (in cM)
//   2nd column: tags (0, neutral; 1, assortative mating locus;
//                     2, loci under divergent selection)
void getdata(string file, parameters * param){
  int l, tag;
  double pos;
  string line, element;
  ifstream infile;
  istringstream stream;

  gsl_vector_int * tmppos;
  gsl_vector_int * tmpord;
  tmppos = gsl_vector_int_calloc(param->iL);
  tmpord = gsl_vector_int_calloc(param->iL);

  // read data
  infile.open(file.c_str());
  if (!infile){
    cerr << "Cannot open file " << file << endl;
    exit(1);
  }

  // read and store map positions
  for(l=0; l<param->iL; l++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    // map position
    stream >> element;
    pos = atof(element.c_str());
    gsl_vector_set(param->d, l, pos);
    // read tags
    stream >> element;
    tag = atoi(element.c_str());
    if(tag == 0){
      ;
    }
    else if(tag == 1){
      if(param->apos >= 0){
	cerr << "Error: only one assortative mating locus allowed" << endl;
	exit(1);
      }
      else{
	param->apos = l;
	param->phases++;
      }
    }
    else if(tag >= 2){
      gsl_vector_int_set(tmppos,param->nsel,l); // store positions temporarily
      gsl_vector_int_set(tmpord,param->nsel,tag - 2); // store order
      param->nsel++;
      param->phases++;
    }
    else {
      cerr << "Error: unknown tag " << tag << " in " << file << endl;
      exit(1);
    }

  }
  infile.close();

  // copy positions of divergently selected loci and sort
  param->spos = gsl_vector_int_calloc(param->nsel);
  for (l=0; l<param->nsel; l++){
    gsl_vector_int_set(param->spos, gsl_vector_int_get(tmpord, l), gsl_vector_int_get(tmppos, l));
  }


  // get map length
  param->mL = gsl_vector_get(param->d, param->iL-1);


  gsl_vector_int_free(tmppos);
  gsl_vector_int_free(tmpord);
}


// -----------------
// cleargrid
void cleargrid(demecont (&grid)[NDEMES]){

  int k;

  for (k=0; k<NDEMES; k++){
    gsl_matrix_set_zero(grid[k].adults);
    gsl_matrix_set_zero(grid[k].chpnt);
    gsl_matrix_set_zero(grid[k].pool1);
    gsl_matrix_set_zero(grid[k].pool2);
    grid[k].N1 = 0;
    grid[k].N2 = 0;
  }

}



// -----------------
// simulate
void simulate(demecont (&grid)[NDEMES], parameters * param, int sub){

  int k, i, l, m, n = 0, s, pheno;
  int ag = 0, nmig, cnt = 0, trynext = 0, surv;
  uint join1 = 1;
  double nmut, a;

  gsl_vector * haplotype;
  gsl_vector * matgam;
  gsl_vector * patgam;
  gsl_matrix * tmpind;
  gsl_vector_int * migrants0;
  gsl_vector_int * migrants1;
  haplotype = gsl_vector_alloc(param->L);
  matgam = gsl_vector_alloc(param->L);
  patgam = gsl_vector_alloc(param->L);
  tmpind = gsl_matrix_alloc(param->L, 2);

  for (k=0; k<NDEMES; k++){

    // --------------
    // (1) assortment
    // --------------
    // Assign individuals to two mating pools -> this will increase
    //   inbreeding a bit even if ass = 0
    // cerr << "assortment" << endl;
    grid[k].N1 = 0;
    grid[k].N2 = 0;
    gsl_matrix_set_zero(grid[k].pool1);
    gsl_matrix_set_zero(grid[k].pool2);

    if(sub == param->phases && param->apos >= 0){ // assortative mating
      // fill two mating pools
      for(i=0; i<param->N; i++){
	// get assortment genotype
	ag = gsl_matrix_get(grid[k].adults,param->apos,2*i);
	ag += gsl_matrix_get(grid[k].adults,param->apos,2*i+1);
	if(ag == 0){
	  // prob to join pool 1 = (1+d)/2
	  join1 = gsl_ran_binomial(r,(1+param->ass)/2.0,1);
	}
	else if(ag == 1){
	  // prob to join pool 1 = 1/2
	  join1 = gsl_ran_binomial(r,1/2.0,1);
	}
	else{
	  // prob to join pool 1 = (1-d)/2
	  join1 = gsl_ran_binomial(r,(1-param->ass)/2.0,1);
	}
	// copy individual accordingly
	if(join1 == 1){
	  gsl_matrix_get_col(haplotype,grid[k].adults,2*i);
	  gsl_matrix_set_col(grid[k].pool1,grid[k].N1*2,haplotype);
	  gsl_matrix_get_col(haplotype,grid[k].adults,2*i+1);
	  gsl_matrix_set_col(grid[k].pool1,grid[k].N1*2+1,haplotype);
	  grid[k].N1++;
	}
	else{
	  gsl_matrix_get_col(haplotype,grid[k].adults,2*i);
	  gsl_matrix_set_col(grid[k].pool2,grid[k].N2*2,haplotype);
	  gsl_matrix_get_col(haplotype,grid[k].adults,2*i+1);
	  gsl_matrix_set_col(grid[k].pool2,grid[k].N2*2+1,haplotype);
	  grid[k].N2++;
	}
      }
    }
    else{ // random mating
      // fill just one mating pool
      gsl_matrix_memcpy(grid[k].pool1,grid[k].adults);
      grid[k].N1 = param->N;
    }

    if(grid[k].N1 + grid[k].N2 != param->N){
      cerr << "Error: something went wrong during assortment" << endl;
      exit(1);
    }


    // -------------------------
    // (2) mutation (all phases)
    // -------------------------
    // cerr << "mutation" << endl;

    // total number of mutation events
    nmut = gsl_ran_binomial(r, param->mu, param->L * param->N * 2);

    while(n < (int) nmut){
      trynext = 0;
      l = gsl_rng_uniform_int(r, param->L); // locus
      // // no mutations allowed at assortative mating and selected loci
      // //   prior to entering the respective phases
      // if(param->nsel > 0 && sub<param->nsel){ // prior selection
      // 	for(s=param->nsel; s>sub && s>0; s--){
      // 	  if(l == gsl_vector_int_get(param->spos,s-1)){
      // 	    trynext = 1;
      // 	  }
      // 	}
      // }
      // if(param->apos >= 0 && sub<param->phases &&
      // 	 l == param->apos){ // prior assort mat
      // 	trynext = 1;
      // }
      // no mutations allowed at assortative mating and selected loci
      //   during any phase
      for(s=0; s<param->nsel; s++){ // selection
	if(l == gsl_vector_int_get(param->spos,s)){
	  trynext = 1;
	}
      }
      if(param->apos >= 0 && l == param->apos){ // assort mat
	trynext = 1;
      }
      if(trynext == 1) continue;
      i = gsl_rng_uniform_int(r, param->N); // individual
      m = gsl_ran_binomial(r,1/2.0,1); // allele
      if(i<grid[k].N1){ // ind in pool1
	a = gsl_matrix_get(grid[k].pool1, l, i*2+m); // old allele
	(a == 0.0) ? a = 1.0 : a = 0.0; // switch
	gsl_matrix_set(grid[k].pool1, l, i*2+m, a);
      }
      else{ // ind in pool2
	i -= grid[k].N1;
	a = gsl_matrix_get(grid[k].pool2, l, i*2+m); // old allele
	(a == 0.0) ? a = 1.0 : a = 0.0; // switch
	gsl_matrix_set(grid[k].pool2, l, i*2+m, a);
      }
      n++;
    }



    // ------------------------------------------------
    // (3) recombination (all phases) and (4) selection
    // ------------------------------------------------
    // cerr << "recombination" << endl;
    gsl_matrix_set_zero(grid[k].adults);

    i = 0; // fill adult generation with same # inds as previous generation
    cnt = 0;
    while(i<param->N && cnt < MAXTRIES){
      surv = 1;
      // choose mom
      m = gsl_rng_uniform_int(r, param->N);
      if(m<grid[k].N1){ // mom in pool1
	// make mom's gamete
	makegamete(matgam, grid[k].pool1, m, param, tmpind);
	// choose dad (selfing allowed) and make gamete
	m = gsl_rng_uniform_int(r, grid[k].N1);
	makegamete(patgam, grid[k].pool1, m, param, tmpind);
      }
      else{ // mom in pool2
	m -= grid[k].N1;
	// make mom's gamete
	makegamete(matgam, grid[k].pool2, m, param, tmpind);
	// choose dad (selfing allowed) and make gamete
	m = gsl_rng_uniform_int(r, grid[k].N2);
	makegamete(patgam, grid[k].pool2, m, param, tmpind);
      }
      // apply selection on zygote
      if(sub > 0 && param->nsel > 0){ // there is selection
	pheno = 0;
	for(s=0; s<sub && s<param->nsel; s++){
	  pheno += gsl_vector_get(matgam,gsl_vector_int_get(param->spos,s));
	  pheno += gsl_vector_get(patgam,gsl_vector_int_get(param->spos,s));
	}
	// deme 1, 1 alleles selected against; deme 2, 0 alleles
	if(k == 1){
	  pheno = s * 2 - pheno;
	}
	a = pow((1.0 - param->sel), (double) pheno);
	if(gsl_rng_uniform(r) >= a){ // inviable zygote
	  surv = 0;
	}
      }
      if(surv == 1){
	// copy to adults
	gsl_matrix_set_col(grid[k].adults, i*2, matgam);
	gsl_matrix_set_col(grid[k].adults, i*2+1, patgam);
	i++;
	cnt--;
      }
      else{
	cnt++;
      }
    }
    if(cnt == MAXTRIES){
      cerr << "Error: not enough survivors" << endl;
      exit(1);
    }
  } // close deme loop

  // --------------------------
  // (5) migration (all phases)
  // --------------------------
  // cerr << "migration" << endl;

  nmig = (int) round(param->N * param->mig);

  if(nmig > 0){
    // set random indices
    migrants0 = gsl_vector_int_calloc(param->N);
    migrants1 = gsl_vector_int_calloc(param->N);
    for(n=0; n<param->N; n++){
      gsl_vector_int_set(migrants0,n,n);
    }
    gsl_vector_int_memcpy(migrants1,migrants0);
    gsl_ran_shuffle(r, migrants0->data, param->N, sizeof(int));
    gsl_ran_shuffle(r, migrants1->data, param->N, sizeof(int));

    // swap individuals
    for(n=0; n<nmig; n++){
      for(a=0; a<2; a++){
	gsl_matrix_get_col(matgam, grid[0].adults,
			   gsl_vector_int_get(migrants0,n) * 2 + a);
	gsl_matrix_get_col(patgam, grid[1].adults,
			   gsl_vector_int_get(migrants1,n) * 2 + a);
	gsl_matrix_set_col(grid[0].adults,
			   gsl_vector_int_get(migrants0,n) * 2 + a, patgam);
	gsl_matrix_set_col(grid[1].adults,
			   gsl_vector_int_get(migrants1,n) * 2 + a, matgam);
      }
    }

    gsl_vector_int_free(migrants0);
    gsl_vector_int_free(migrants1);
  }

  gsl_vector_free(haplotype);
  gsl_vector_free(matgam);
  gsl_vector_free(patgam);
  gsl_matrix_free(tmpind);
}



// -----------------
// makegamete
void makegamete(gsl_vector * gamete, gsl_matrix * deme, int ind,
		parameters * param, gsl_matrix * tmpind){

  int ncross = 0, id, n, pos, l;
  gsl_vector * crosspos;

  // copy data for ind to tmpind
  gsl_matrix_get_col(gamete, deme, ind * 2);
  gsl_matrix_set_col(tmpind, 0, gamete);
  gsl_matrix_get_col(gamete, deme, ind * 2 + 1);
  gsl_matrix_set_col(tmpind, 1, gamete);

  // number of crossovers
  ncross = gsl_ran_poisson(r, param->mL);

  // choose gamete to start with
  id = gsl_ran_binomial(r, 0.5, 1);

  if (ncross == 0){
    gsl_matrix_get_col(gamete, tmpind, id);
    // unlinked loci are overwritten below
  }
  else { // there are crossovers
    // positions of crossovers
    crosspos = gsl_vector_alloc(ncross);
    for(n=0; n<ncross; n++){
      gsl_vector_set(crosspos, n, gsl_ran_flat(r, 0, param->mL));
    }
    gsl_sort_vector(crosspos);

    n = 0;
    pos = 0;

    // copy genetic data
    while(n < ncross &&
	  gsl_vector_get(crosspos, n) < gsl_vector_get(param->d, 0)){
      n++; // ignore crossovers that are before first locus
    }
    for(l=0; l<param->iL; l++){
      if (n == ncross){ // this was last breakpoint
	while(pos<param->iL){
	  // copy all loci after last breakpoint
	  gsl_vector_set(gamete, pos, gsl_matrix_get(tmpind, pos, id));
	  pos++;
	}
	break;
      }
      else if (gsl_vector_get(crosspos, n) <= gsl_vector_get(param->d, l)){
	// copy all loci before crosspos to gamete
	while(pos<l){
	  // copy all loci before breakpoint
	  gsl_vector_set(gamete, pos, gsl_matrix_get(tmpind, pos, id));
	  pos++;
	}
	n++;
	id = (int) pow( (double) (id - 1), (double) id); // switch id
      }
      else{
	; // go to next locus
      }
    }
    gsl_vector_free(crosspos);
  }

  // unlinked loci
  for(l=param->iL; l<param->L; l++){
    id = gsl_ran_binomial(r, 0.5, 1);
    gsl_vector_set(gamete, l, gsl_matrix_get(tmpind, l, id));
  }

  // cerr << "made gamete" << endl;
}



// -----------------
// copygrid
void copygrid(demecont (&grid)[NDEMES], int dir){

  int k;

  if(dir == 0){ // copy from adults to chpnt
    for (k=0; k<NDEMES; k++){
      gsl_matrix_memcpy(grid[k].chpnt, grid[k].adults);
    }
  }
  else if(dir == 1){ // copy from chpnt to adults
    for (k=0; k<NDEMES; k++){
      gsl_matrix_memcpy(grid[k].adults, grid[k].chpnt);
    }
  }

}



// -----------------
// changephase
void changephase(demecont (&grid)[NDEMES], parameters * param, int sub){

  int i, m, k, mutate = 0, n, a, cnt = 0;
  gsl_vector * alleles;

  alleles = gsl_vector_alloc(param->N*2);

  if(param->nsel > 0 && sub<=param->nsel){ // selection phase
    // verify that there are no mutations present already
    for (k=0; k<NDEMES; k++){
      gsl_matrix_get_row(alleles,grid[k].adults,
			 gsl_vector_int_get(param->spos,sub-1));
      if(gsl_vector_max(alleles)>0){
	cerr << "Error: mutation already present\n" << endl;
	exit(1);
      }
    }

    // mutate allele 0 to 1 in deme 2
    i = gsl_rng_uniform_int(r, param->N); // individual
    m = gsl_ran_binomial(r,1/2.0,1); // allele
    gsl_matrix_set(grid[1].adults,
		   gsl_vector_int_get(param->spos,sub-1),i*2+m,1);
    cerr << "mutated selected locus " << sub << endl;
  }
  if(param->apos >= 0 && sub==param->phases){ // assort mat phase
    // verify that there are no mutations present already
    for (k=0; k<NDEMES; k++){
      gsl_matrix_get_row(alleles,grid[k].adults,param->apos);
      if(gsl_vector_max(alleles)>0){
	cerr << "Error: mutation already present\n" << endl;
	exit(1);
      }
    }
    if(param->mutone == 0){
      // find haplotypes carrying alt selected alleles, if any
      mutate = 0;
      for (k=0; k<NDEMES; k++){ // deme
	for (i=0; i<param->N; i++){ // individual
	  for (m=0; m<2; m++){ // allele
	    a = 0;
	    for(n=0; n<param->nsel; n++){ // test haplotype
	      a += gsl_matrix_get(grid[k].adults,
				  gsl_vector_int_get(param->spos,n), i*2+m);
	    }
	    if(param->nsel == a){
	      if(gsl_ran_binomial(r, param->mutfrac, 1) == 1){
		gsl_matrix_set(grid[k].adults,param->apos,i*2+m,1);
		mutate++;
	      }
	    }
	    else if(a == 0){
	      if(gsl_ran_binomial(r, param->mutfrac, 1) == 0){
		gsl_matrix_set(grid[k].adults,param->apos,i*2+m,1);
		mutate++;
	      }
	    }
	    else{
	      if(gsl_ran_binomial(r, a/(1.0 * param->nsel), 1) == 1){
		gsl_matrix_set(grid[k].adults,param->apos,i*2+m,1);
		mutate++;
	      }
	    }
	  }
	}
      }
      if(mutate == 0){
	cerr << "Error: couldn't find suitable haplotype to mutate assortative mating locus" << endl;
	exit(1);
      }
      else if(mutate == 1){
	cerr << "mutated one allele at assortative mating locus" << endl;
      }
      else{
	cerr << "mutated " << mutate << " alleles at assortative mating locus" << endl;
      }
    }
    else if(param->mutone == 1){
      while(mutate == 0 && cnt < MAXTRIES){
	cnt++;
	// find haplotype carrying alt selected alleles, if any
	i = gsl_rng_uniform_int(r, param->N); // individual
	m = gsl_ran_binomial(r,1/2.0,1); // allele
	a = 0;
	for(n=0; n<param->nsel; n++){
	  a += gsl_matrix_get(grid[1].adults,
			      gsl_vector_int_get(param->spos,n), i*2+m);
	}
	if(param->nsel - a > 0) continue;
	mutate = 1;
      }
      if(mutate == 0){
	cerr << "Error: couldn't find suitable haplotype to mutate assortative mating locus" << endl;
	exit(1);
      }
      // mutate allele 0 to 1 in deme 2
      gsl_matrix_set(grid[1].adults,param->apos,i*2+m,1);
      cerr << "mutated one allele at assortative mating locus" << endl;
    }
  }

  gsl_vector_free(alleles);

}



// -----------------
// writeOut
void writeOut(demecont (&grid)[NDEMES], parameters * param,
	      int rep, int gen, int sub, string file){

  ofstream outfile;
  int l;

  gsl_vector * p0;
  gsl_vector * p1;
  gsl_vector * fst;
  gsl_vector * Zg;

  p0 = gsl_vector_alloc(param->L);
  p1 = gsl_vector_alloc(param->L);
  fst = gsl_vector_alloc(param->L);
  Zg = gsl_vector_alloc(param->L);

  calcAf(p0, grid[0].adults, param);
  calcAf(p1, grid[1].adults, param);
  calcFst(fst, p0, p1, param);
  calcZg(Zg, p0, p1, param);

  outfile.open(file.c_str(), ios::out | ios::app );
  if (outfile.is_open()){
    // print allele frequency deme 1
    outfile << rep+1 << " " << gen << " " << sub+1;
    for(l=0; l<param->L; l++){
      outfile << " " << gsl_vector_get(p0, l);
    }
    // print allele frequency deme 2
    outfile << "\n" << rep+1 << " " << gen << " " << sub+1;
    for(l=0; l<param->L; l++){
      outfile << " " << gsl_vector_get(p1, l);
    }
    // print fst
    outfile << "\n" << rep+1 << " " << gen << " " << sub+1;
    for(l=0; l<param->L; l++){
      outfile << " " << gsl_vector_get(fst, l);
    }
    // print Zg
    outfile << "\n" << rep+1 << " " << gen << " " << sub+1;
    for(l=0; l<param->L; l++){
      outfile << " " << gsl_vector_get(Zg, l);
    }
    outfile << "\n";
    outfile.close();
  }
  else {
    cerr << "Cannot open file " << file << endl;
    exit(1);
  }

  gsl_vector_free(p0);
  gsl_vector_free(p1);
  gsl_vector_free(fst);
  gsl_vector_free(Zg);

}


// -----------------
// calcAf
void calcAf(gsl_vector * pvec, gsl_matrix * deme, parameters * param){

  int l, i;
  double p;

  gsl_vector * alleles;
  alleles = gsl_vector_alloc(param->N * 2);

  for(l=0; l<param->L; l++){
    gsl_matrix_get_row(alleles, deme, l);
    p = 0.0;
    for(i=0; i<param->N * 2; i++){
      p += gsl_vector_get(alleles, i);
    }
    p /= (double) (param->N * 2);
    gsl_vector_set(pvec, l, p);
  }

  gsl_vector_free(alleles);

}



// -----------------
// calcFst
// Hudson's Fst, formula as in Soria-Carrasco et al. 2014
void calcFst(gsl_vector * fstvec, gsl_vector * p0vec, gsl_vector * p1vec,
	     parameters * param){

  int l, cnt, n = 0;
  double fst, hw = 0.0, hb = 0.0, p0, p1;

  for(l=0, cnt=0; l<param->L; l++, cnt++){
    p0 = gsl_vector_get(p0vec, l);
    p1 = gsl_vector_get(p1vec, l);
    if (p0+p1 > 0 && p0+p1 < 2){ // ignore invariant loci
      hw += p0 * (1-p0) + p1 * (1-p1);
      hb += p0 * (1-p1) + p1 * (1-p0);
      n++;
    }
    if(cnt == param->bins - 1){ // calc average over bin
      hw /= (double) n;
      hb /= (double) n;
      fst = 1 - hw/hb;
      while(cnt >= 0){ // write
	gsl_vector_set(fstvec, l - cnt, fst);
	cnt--;
      }
      hw = 0.0;
      hb = 0.0;
      n = 0;
    }
  }
  while(cnt > 0){ // write nan for remaining loci, if any
    gsl_vector_set(fstvec, l - cnt, NAN);
    cnt--;
  }
}


// -----------------
// calcZg
// calculate Zg (Storz & Kelly 2008)
void calcZg(gsl_vector * zgvec, gsl_vector * p0vec, gsl_vector * p1vec,
	     parameters * param){

  int l, m, j, cnt, var = 0, n, nz = 0;
  double Zg, num0, num1, denom, pj, pm, tmp;

  gsl_vector * p0;
  gsl_vector * p1;

  n = (param->bins * (param->bins - 1)) / 2;

  p0 = gsl_vector_calloc(param->bins);
  p1 = gsl_vector_calloc(param->bins);

  for(l=0, cnt=0; l<param->L; l++, cnt++){
    pj = gsl_vector_get(p0vec, l);
    pm = gsl_vector_get(p1vec, l);
    if(pj+pm > 0 && pj+pm < 2){ // ignore invariant loci
      gsl_vector_set(p0, var, pj);
      gsl_vector_set(p1, var, pm);
      var++;
    }
    if(cnt == param->bins - 1){ // calc Zg for bin
      for(j=0; j<var - 1; j++){
	for(m=j+1; m<var; m++){
	  // mean allele frequencies both demes for locus pair
	  pj = (gsl_vector_get(p0, j) + gsl_vector_get(p1, j)) / 2.0;
	  pm = (gsl_vector_get(p0, m) + gsl_vector_get(p1, m)) / 2.0;
	  num0 = (1/2.0)*(gsl_vector_get(p0, j)*gsl_vector_get(p0, m) - pj*pm);
	  num1 = (1/2.0)*(gsl_vector_get(p1, j)*gsl_vector_get(p1, m) - pj*pm);
	  // (1/2.0) is for equal deme sizes
	  denom = pj * (1-pj) * pm * (1-pm);
	  tmp = pow(num0+num1, 2.0) / denom;
	  if(isfinite(tmp) && tmp <= 1){
	    // can be zero for fixed alleles or > 1 or Inf with some rounding issues -> get valid values to calculate their mean
	    Zg += tmp;
	    nz++;
	  }
	}
      }
      Zg /= (double) nz;
      while(cnt >= 0){ // write
	gsl_vector_set(zgvec, l - cnt, Zg);
	cnt--;
      }
      gsl_vector_set_zero(p0);
      gsl_vector_set_zero(p1);
      Zg = 0.0;
      nz = 0;
      var = 0;
    }
  }
  while(cnt > 0){ // write nan for remaining loci, if any
    gsl_vector_set(zgvec, l - cnt, NAN);
    cnt--;
  }

  gsl_vector_free(p0);
  gsl_vector_free(p1);

}




