/*
 *  hist1.c
 *
*
 */

#include <Rcpp.h>
//#include <stdio.h>
#include <math.h>
//#include <unistd.h>
//#include <string.h>
#include <string>

#include "ranfun.h"
#include "Matrix.h"
//    c0    thr1   alpha1      c1  m1   thr2   alpha2      c2  m2   thr3   alpha3      c3  m3    LogPost

class Hist {

protected:

	Matrix *OutB;
	SubMatrix Out;
	SubMatrix A;

	int it, K, cols, extrapol_warning;

	double d; //depth
	double th0, th1, Dth;

	//Based depth and increment between depths
	double c0, Dc;
	double c(int i) { return c0 + i*Dc; }

	double Model(double d, double thr, double alpha, double dr)
	{
		return thr + alpha*(d-dr);
	}

	double dr1(int i, int sec)  { return Out( i, 1 + sec*4 - (sec == 0 ? 1 : 2));}
	double thr(int i, int sec)   { return Out( i, 1 + sec*4); }
	double alpha(int i, int sec) { return Out( i, 1 + sec*4 + 1); }
	double dr(int i, int sec)     { return Out( i, 1 + sec*4 + 2); }

public:

	//read all data, # of simulations # of sections
	Hist( char *fnam, int itt, int KK, double cc0, double DDc) {

		it = itt;
		K = KK;
		c0 = cc0;
		Dc = DDc;
//		Rprintf("Dc = %f  DDc= %f\n", Dc, DDc); // tmp MB

		extrapol_warning = 0;

		//    th0  x's w   U
		cols = 1 + K + 1 + 1;

		OutB = new Matrix( it, cols);

		Out.Set( OutB, it, cols);
		Out.filescan(fnam);

		//A.Set( OutB, 20, cols);
		//A.print();
	}

	~Hist() {
		A.~SubMatrix();
		Out.~SubMatrix();
		delete OutB;
	}


	double Model(int i) {

						//c0
        if (fcmp( d, c0) == -1) { //d < c0
            Rprintf("hist: ERROR: d = %6.4f < c0= %6.4f!!\n", d, c0);
            Rcpp::stop("hist: ERROR: d = %6.4f < c0= %6.4f!!\n", d, c0);
            //exit(0);
		}

		double S=Out( i, 0); //th0
		if (fcmp(d, c(1)) == -1) // tmp MB to try and correct th0 bug
			return S + Out(i, 1)*(d-c(0)); // tmp MB to try and correct th0 bug

		for (int k=1; k<K; k++) {
			S += Out( i, k)*Dc;
			if (fcmp( d, c(k+1)) == -1)
				return S + Out( i, k+1)*(d-c(k));
		}

		if (extrapol_warning <= 0) {
            Rprintf("hist: WARNING: extrapolation, depth d = %f above cK = %f\n", d, c(K));
		}
		return S + Out( i, K)*(d-c(K));

	}


	// calculate the counts at depth d, n=number of divisions, buffer
	void GetHist( double dd, double n, int *hi) {

		th0=1e300, th1=-1e300;
		double th;

		d = dd;

		int j;
		for (j=0; j<n; j++)
			hi[j] = 0;

		for (int i=0; i<it; i++) {

			th = Model(i);

			th0 = fmin( th, th0);
			th1 = fmax( th, th1);
		}

		// Use the above and below integers
		th0 = floor(th0);
		th1 = ceil(th1);

		Dth = (th1-th0)/(double) n;

		for (int i=0; i<it; i++) {

			th = Model(i);

			j = (int) floor((th-th0)/Dth);

			if ((j < 0) || (n <= j))
                Rprintf("i= %3d, d= %6.4f, th= %6.4f, j= %d\n", i, d, th, j);

			hi[j]++;
		}

	}

	double GetTh0() { return th0; }
	double GetTh1() { return th1; }
	double GetDth() { return Dth; }
	int Get_extrapol_warnings() { return extrapol_warning; }
};


// [[Rcpp::export]]
void hist2(std::string MCMCsamplesfname1,int  samplesize, int c0 , double Dc, int K, int n, std::string fout1, int n_depths, std::string dfin1) {

  char *MCMCsamplesfname = new char[MCMCsamplesfname1.length() + 1];
  strcpy(MCMCsamplesfname, MCMCsamplesfname1.c_str());
  char *fout = new char[fout1.length() + 1];
  strcpy(fout, fout1.c_str());
  char *dfin = new char[dfin1.length() + 1];
  strcpy(dfin, dfin1.c_str());



	FILE *F;
    //JEV warning
    /*if (strcmp( fout, "-") == 0)
       //  F = stdout;
    else*/
		if ((F = fopen( fout, "w+")) == NULL) {
            Rprintf("Could not open file %s for writing.\n", fout); // JEV WARNING CHECK
            Rcpp::stop("Could not open file %s for writing.\n", fout);

            //exit(0);
		}

 	FILE *fr;
 	if ((fr = fopen(dfin, "r")) == NULL) {
//        Rprintf("hist: ERROR, depths file %s not found.", dfin);
        Rcpp::stop("hist: ERROR, depths file %s not found.", dfin);
// 		exit(0);
 	}

   	double *depth = new double[n_depths];
    int rtn = 0; //JEV warning
	for(int i=0; i<n_depths; i++) { // extract the depths
  	  rtn = fscanf(fr, " %lf", &depth[i]);
     if(rtn){}// JEV warning
   	}
	fclose(fr);

	Hist Hi( MCMCsamplesfname, samplesize, K, c0, Dc);

	int *hi = new int[n];

	fprintf( F, "### file: %s, it= %d\n\nhists <- NULL;\n\n", MCMCsamplesfname, samplesize);
	for (int i=0;  i<n_depths; i++) {
		Hi.GetHist( depth[i], n, hi);

		fprintf( F, "hists <- append( hists, pairlist(");
		fprintf( F, "list( d=%f, th0=%6.0f, th1=%6.0f, n=%d, Dh=%f, ss=%d, counts=c( ",
			depth[i], Hi.GetTh0(), Hi.GetTh1(), n, Hi.GetDth(), samplesize);

		for (int i=0; i<n-1; i++)
			fprintf( F, " %d,", hi[i]);
		fprintf( F, " %d))))\n", hi[n-1]);
	}

	if (0 < Hi.Get_extrapol_warnings())
        Rprintf("hist: WARNING: %d extrapolation warnings.\n", Hi.Get_extrapol_warnings());

	delete[] hi;
	delete[] depth;
	fclose(F);
}

