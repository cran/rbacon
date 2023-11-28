/*

 BACON:

Implements the semiparametric autoregressive model for
radiocarbon chronologies, using the twalk.  See paper for mathematical details and
the files:

- bacon.h: This is the implementation, with the model definitions etc.

- cal.h: Reads and manages calibration curves and determinations
 g++ -I/usr/share/R/include -DNDEBUG -W -I/usr/include -I../inst/include  -I"/usr/lib/R/site-library/Rcpp/include"   -fpic  -O3 -pipe  -g  -c bacon.cpp -o bacon.o
- input.h, input.c: reads the input files and stores all data

- ranfun.h, twalk.h, Matrix.h: for some gsl interfaces for random number generation, the C++ twalk implementation and a simple Matrix class.

 */
#include <Rcpp.h>
//#include <stdio.h>
#include <math.h>
//#include <unistd.h>
//#include <string.h>
#include <string>

#include "bacon.h"
#include "input.h"
#include "cal.h"
#include "ranfun.h"
#include "Matrix.h"
#include "twalk.h"
//using namespace Rcpp;

#define MAXNUMOFCURVES 100
#define MAXNUMOFDETS  1000

#define BUFFSIZE 4000

//the "every" thinning subsampling parameter is EVERY_MULT*All.Dim()
//when we were saving only accepted iterations EVERY_MULT was 5
#define EVERY_MULT 25
//The burn in is BURN_IN_MULT*All.Dim(), this was 200
#define BURN_IN_MULT 3000

//Every how many iterations we expect and acceptance: inverse of the acceptance rate
#define ACCEP_EV 20

// [[Rcpp::export]]
int bacon( std::string inputfile1, std::string outputfile1, int ssize, std::string dircc) {// Command line: bacon inputfile outputfile
//
  //Rprintf("inputfile is: %s,\n outputfile: %s,\n ssize %d,\n dircc %s", inputfile1, outputfile1, ssize, dircc);
  char *inputfile = new char[inputfile1.length() + 1];
//  printf("%s\n", inputfile1.c_str());
  strcpy(inputfile, inputfile1.c_str());
  //JEV avoid warning int i=0; i<inputfile1.length() ; i++
//  for (unsigned int i=0; i<inputfile1.length() ; i++){   // MB remove this Oct 2020 as it just repeats printing the name/location of the bacon file?
//     printf("'%c' ", inputfile1.c_str()[i]);
//     Rprintf("'%c' ", inputfile1.c_str()[i]);
//    Rprintf("%c", inputfile1.c_str()[i]);
//  }

  char *outputfile = new char[outputfile1.length() + 1];
  //  printf("%s\n", inputfile1.c_str());
  strcpy(outputfile, outputfile1.c_str());

  // do stuff
  // delete [] inputfile;

//  char inputfile[]=inputfile1.c_str() ; //"EDENFULL_31.bacon";
//  char outputfile[]="test.out";
  //printf("%s %s %d \n",inputfile,outputfile,ssize);

  //	if (argc < 4) {
  // 	printf("Usage: bacon inputfile outputfile ssize\n");
  //
  // 	exit(0);
  // }

  // char  ax[BUFFSIZE];
  //Rprintf("antes!\n");

  //Program file
  //sprintf( ax, "Cores/%s/%s.bacon", argv[1], argv[2]);
  // sprintf( ax, "%s", argv[1]);
  //Read everything from the program file
  Input All( inputfile, MAXNUMOFCURVES, MAXNUMOFDETS, dircc);
  //Rprintf("despues!\n");

  //File to save the twalk output
  // sprintf( ax, "%s", argv[2]);

  // int ssize;
  // sscanf( argv[3], " %d", &ssize);

  //ssize is the final sample size needed
  //ssize = it/(ACCEP_EV * All.Dim() * EVERY_MULT) - BURN_IN_MULT
  //Then we let

  //this was the previous definition of number of iterations, when we were saving accepted iterations only
  //int it = ACCEP_EV * All.Dim() * EVERY_MULT * (ssize + BURN_IN_MULT);
  int it = EVERY_MULT * All.Dim() * (BURN_IN_MULT + ssize);

  //int every=  -1*EVERY_MULT*All.Dim(); // only accepted iterations
  int every =  EVERY_MULT*All.Dim();

  
  //Run the twalk
  All.RunTwalk( outputfile, it, every);
  All.PrintNumWarnings();

  All.outputFiles(outputfile1); // this was not present in rbacon's bacon.cpp! March 2021

  /*
  char  ax2[BUFFSIZE];
  //File to save the thinned twalk output
  sprintf( ax2, "%s", argv[2]);


  FILE *F, *G;
  if ((F = fopen( ax, "r")) == NULL) {
  printf("Could not open %s for reading.\n", ax);

  //exit(0);
  }
  if ((G = fopen( ax2, "w+")) == NULL) {
  printf("Could not open %s for writing.\n", ax2);

  //exit(0);
  }


  char ln[CHARBUFFER];
  int burnin= BURN_IN_MULT*All.Dim(), j=0;
  every = 1;
  //subsample the twalk output:
  printf("\nSubsampling %s, burnin= %d, every= %d and store results in %s ...\n", ax, burnin, every, ax2);
  while (!feof(F)) {
  fgets( ln, CHARBUFFER, F);
  if (j == burnin) break;
  j++;
  }

  j=0;
  int ss=0;
  while (!feof(F)) {

  fgets( ln, CHARBUFFER, F);
  if ((j % every) == 0) {

  fputs( ln, G);
  ss++;
  }
  j++;
  }
  printf("out.all size= %d\n", j);

  fclose(F);
  fclose(G);

  printf("Removing %s\n", ax);
  remove(ax);

  printf("Final sample size %d\n", ss);
  */

 Rprintf("bacon: burn in (initial iterations which will be removed): %d\n", All.Dim() * EVERY_MULT * BURN_IN_MULT);

if (Un01() < 0.5) {
	Rprintf(FAREWELL);} else
		if(Un01() < 0.5) {
			Rprintf("Ats us nai!\n");} else
				if(Un01() < 0.2) {
					Rprintf("... sizzle spatter sizzle...\n");} else
						if(Un01() < 0.2) {
							Rprintf("... adding maple...\n");} else
								if(Un01() < 0.5) {
						 			Rprintf("Looking good, turning off the fire\n\n");} else
										{Rprintf("Remember, never pour grease down the drain!\n");};

//  printf("bacon: suggested burn in= %d\n", All.Dim() * EVERY_MULT * BURN_IN_MULT);
 // printf(FAREWELL);

  return All.Dim() * EVERY_MULT * BURN_IN_MULT;

}
