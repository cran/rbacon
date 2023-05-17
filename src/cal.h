
#ifndef CAL_H
#define CAL_H

//#include <stdio.h>
#include <math.h>
//#include <unistd.h>
#include <string>

#include "ranfun.h"
#include "Matrix.h"

//#define IntCal13FNAM "Curves/3Col_intcal13.14C"
#define IntCal13FNAM "3Col_intcal13.14C"
#define IntCal13ROWS 5141
#define IntCal13COLS 3

//#define Marine13FNAM "Curves/3Col_marine13.14C"
#define Marine13FNAM "3Col_marine13.14C"
#define Marine13ROWS 4801
#define Marine13COLS 3

//#define SHCal13FNAM "Curves/3Col_shcal13.14C"
#define SHCal13FNAM "3Col_shcal13.14C"
#define SHCal13ROWS 5141
#define SHCal13COLS 3

//#define IntCal20FNAM "Curves/3Col_intcal20.14C"
#define IntCal20FNAM "3Col_intcal20.14C"
#define IntCal20ROWS 9501
#define IntCal20COLS 3

//#define Marine20FNAM "Curves/3Col_marine20.14C"
#define Marine20FNAM "3Col_marine20.14C"
#define Marine20ROWS 5501
#define Marine20COLS 3

//#define SHCal20FNAM "Curves/3Col_shcal20.14C"
#define SHCal20FNAM "3Col_shcal20.14C"
#define SHCal20ROWS 9501
#define SHCal20COLS 3

#define GENCCMAXLINLEN 255
#define GENCCCOLS 3

#define POSTBOMBFNAMS	"None", \
                        "postbomb_NH1.14C", \
                        "postbomb_NH2.14C", \
                        "postbomb_NH3.14C", \
                        "postbomb_SH1-2.14C", \
                        "postbomb_SH3.14C"

/**** Template class to hold a calibration curve and perform all necessary calibrations. ****/
class Cal {

protected:

	int k;
	double mu, sig;

public:

    Cal(int kk) { k = kk;  }

	double GetSig() { return sig; }
	double GetMu() { return mu; }

    virtual const char *Name() = 0;

	virtual double cal(double theta) = 0;
	virtual double U( double y, double vr, double theta) = 0;
	virtual double Ut( double y, double vr, double theta, double a, double b) = 0;

	virtual double MinCal() = 0;
	virtual double MaxCal() = 0;

};


//Constant "calibration curve", no curve, basically
class ConstCal : public Cal {

public:

    ConstCal() : Cal(0) {
        Rprintf("Constant calibration curve.\n");
	}

	double cal(double theta) {
		mu = theta;
		sig = 0.0;

		return theta;
	}

     const char* Name() { return "Constant c. curve"; }

	double U( double y, double vr, double theta) { return 0.5*sqr(y - theta)/vr; }
	double Ut( double y, double vr, double theta, double a, double b) { return (a + 0.5)*log( b + 0.5*sqr(y - theta)/vr); }

	double MinCal() { return -1.0e300; }
	double MaxCal() { return 1.0e300; }
};




//Generic Cal. curve.  Three columns, no header, first col. ascending order BP
class GenericCal : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	int numrows, min, max, mid;
	char name[1024];
	double mincal, maxcal, const2;

public:

    GenericCal(const char *fnam, std::string ccdir ) : Cal(0) {

        std::string filename= std::string(ccdir)+std::string(fnam);

		//Count the number of lines in the file
		FILE *F;
        if ((F = fopen( filename.c_str(), "r")) == NULL) {
            REprintf("Cal: ERROR: Could not find generic cal. curve, file not found: %s\n", filename.c_str());
            Rcpp::stop("Cal: ERROR: Could not find generic cal. curve, file not found: %s\n", filename.c_str());
            //exit(0);
		}

		numrows=0;
		char ln[GENCCMAXLINLEN];
		while (!feof(F)) {
            char* result=fgets(ln, GENCCMAXLINLEN, F); //JEV warning
            if(result==NULL){}//JEV warning

			numrows++;
		}
		numrows--;
		fclose(F);

		//Read the file as usual:
		CCB = new Matrix( numrows, GENCCCOLS);

		CC.Set( CCB, CCB->nRow(), CCB->nCol());

        Rprintf("GenericCal: Reading from file: %s, %d rows, 3 cols.\n", filename.c_str(), numrows);

        if (CC.filescan((char*)filename.c_str()) == 0) {
          REprintf("Cal: ERROR: Could not find generic cal. curve, file not found: %s\n", filename.c_str());
          Rcpp::stop("Cal: ERROR: Could not find generic cal. curve, file not found: %s\n", filename.c_str());
        //	exit(0);
		}

		//Set mincal and maxcal
		mincal = CC(0,0);
		maxcal = CC(numrows-1,0);

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library

		// sprintf( name, "Generic cal. curve %s", filename.c_str());
		snprintf( name, sizeof(name), "Generic cal. curve %s", filename.c_str());  // MB Dec 2022
	}

	~GenericCal() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double cal(double theta) {


		//Find k, the correct knot.
		//We need k such that:  CC(k,0) <= theta < CC(k+1,0)


		if (fcmp(theta, mincal) == -1)
			//fprintf( stderr, "WARNING: Calibration attempted beyond Generic cal. curve %s limits, theta= %f\n", fnam, theta);
			k = 0; //extrapolation
		else
			if (fcmp(theta, maxcal) == -1) {
				//Binary search:

				min = 0;
				max = numrows-1;
				mid = (min + max)/2; //Integer division
				            //CC( mid, 0) <= theta < CC( mid+1, 0)
				while (!( (fcmp(CC( mid, 0), theta) <= 0) && (fcmp(theta, CC( mid+1, 0)) == -1) )) {

					if (fcmp(theta, CC( mid, 0)) == 1) //theta > CC( mid, 0)
						min = mid + 1;
					else
						max = mid - 1;
					mid = (min + max)/2;
				}
				k = mid;
			}
			else
			//fprintf( stderr, "WARNING: Calibration attempted beyond Generic cal. curve %s limits, theta= %f\n", fnam, theta);
				k = numrows-2; //extrapolation

		//printf(" %d: %6.3f %6.3f %6.3f\n", k, CC( k, 0), theta, CC( k+1, 0));

		mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/(CC(k+1,0)-CC(k,0));
		sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/(CC(k+1,0)-CC(k,0));

		return mu;
	}

    const char *Name() { return name; } //JEV warning


	double U( double y, double vr, double theta)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}


	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}


	double MinCal() { return mincal; }
	double MaxCal() { return maxcal; }
};


/* draft curve for IntCal20 */
class IntCal20 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	int Bomb;
	Cal *bombcc;
	char name[255];
	double mincal, const2;

public:

	IntCal20(int bomb, std::string ccdir) : Cal(IntCal20ROWS) {

	CCB = new Matrix( IntCal20ROWS, IntCal20COLS);
	CC.Set( CCB, CCB->nRow(), CCB->nCol());

	std::string filename= ccdir+IntCal20FNAM;

	Rprintf("IntCal20: Reading from file: %s\n", filename.c_str());


	if (CC.filescan((char*)filename.c_str()) == 0) {
		REprintf("Cal: ERROR: Could not find IntCal20 cal. curve, file not found: %s\n", filename.c_str());
		Rcpp::stop("Cal: ERROR: Could not find IntCal20 cal. curve, file not found: %s\n", filename.c_str());
		//  exit(0);
	}

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library

		const char *postbombfnam[] = { POSTBOMBFNAMS };

		/** Read bomb **/
		Bomb = bomb;
		if (Bomb == 0) {
			mincal = 0.0; // no bomb; 17 Dec 2018 changed -5.0 to 0.0
			//sprintf( name, "IntCal20");
			snprintf( name, sizeof(name), "IntCal20");
		}
		else
			if (Bomb < 6) { // curve number, not cal BP yr. 26 March 2019: Was Bomb < 5 but now there are 5 postbomb curves

			bombcc = new GenericCal(postbombfnam[Bomb], ccdir);
			mincal = bombcc->MinCal();
			//sprintf( name, "IntCal20+%s", postbombfnam[Bomb]);
			snprintf( name, sizeof(name), "IntCal20+%s", postbombfnam[Bomb]);
			}
			else {
				REprintf("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
				Rcpp::stop("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
//				exit(0);
			}
	}

	~IntCal20() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return mincal; }
	double MaxCal() { return 55000.0; } // was 50000.0

	virtual const char *Name() { return name; } //JEV warning


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
		if (fcmp(theta, -0.0) == -1) // value < 0 cal BP, postbomb
		{
			if (Bomb == 0) {
				//fprintf( stderr, "WARNING: Calibration attempted beyond IntCal20 cal. curve limits, theta= %f\n",theta);
				k = 0; //Extrapolation
				mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
				sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
			}
			else {
				bombcc->cal(theta);
				mu = bombcc->GetMu();
				sig = bombcc->GetSig();
			}

		}
		else { // until 5,000 cal BP, line 4999 (check!), values every year
			if (fcmp(theta, 5000.0) != 1) // value < 5,000 cal BP
				{
					k = 0 + (int) floor(theta/1.0); // 0 at start was 1
					mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/1.0;
					sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/1.0;
				}
				else // from 5,000 until 15,000 cal BP, line 7000 (check!), values every 5 years
					if (fcmp(theta, 15000.0) != 1)
						{
							k = 4999 + (int) floor((theta-5000.0)/5.0);
							mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
							sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
						}
						else // from 15,000 until until 25,000 cal BP, line 8000 (check!), values every 10 years
							if (fcmp(theta, 25000.0) != 1)
								{
									k = 7000 + (int) floor((theta-15000.0)/10.0);
									mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
									sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
								}
								else // from 25,000 until 50,000 cal BP, line 9250, values every 20 years
									if (fcmp(theta, 50000.0) != 1)
										{
											k = 8000 + (int) floor((theta-25000.0)/20.0);
											mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
											sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
										}
										else // value > 50,000 cal BP, linear "estimated" cc
											{
                                                //The  cc is a straight line from the value at
                                                //50000 to:
                                                double last_th  = 100000.0;
                                                double last_mu  = 95840.0;
                                                double last_sig = 10000.0;
												k = 9250; //Row at 50000 cal BP
                                                mu = CC(k,1) + (theta-CC(k,0))*(last_mu-CC(k,1))/(last_th-CC(k,0));
                                                sig = CC(k,2) + (theta-CC(k,0))*(last_sig-CC(k,2))/(last_th-CC(k,0));
											}
		}
		return mu;
	}

	double U( double y, double vr, double theta)
	{
		cal(theta);

		double tau = 1.0/(vr + sqr(sig));

		return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
		cal(theta);

		double tau = 1.0/(vr + sqr(sig));

		return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}
};







/////////////////////////////////////
// Marine20
////////////////////////////////////

class Marine20 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	double const2;

public:

    Marine20(std::string ccdir) : Cal(Marine20ROWS) {

		CCB = new Matrix( Marine20ROWS, Marine20COLS);

		CC.Set( CCB, CCB->nRow(), CCB->nCol());

        std::string filename= std::string(ccdir)+std::string(Marine20FNAM);

        Rprintf("Marine20: Reading from file: %s\n", filename.c_str());

        if (CC.filescan((char*)filename.c_str()) == 0) {
           REprintf("Cal: ERROR: Could not find Marine20 cal. curve, file not found: %s\n", filename.c_str());
           Rcpp::stop("Cal: ERROR: Could not find Marine20 cal. curve, file not found: %s\n", filename.c_str());

//			exit(0);
		}

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library


	}

	~Marine20() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return 0.0; }
	double MaxCal() { return 55000.0; } // was 50000.0

 const char* Name() { return "Marine20"; }


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
        if (fcmp(theta, 0.0) == -1)
        {
                //fprintf( stderr, "WARNING: Calibration attempted beyond marine20 cal. curve limits, theta= %f\n",theta);
                k = 0;
                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5;
                //sig <- CC(k,2); before JUDY
                 sig = CC(k,2);
        }
        else   // cal BP jumps are 10 yr throughout the curve
        {      
            if (fcmp(theta, 55000.0) != 1)
               {
                  k = 0 + (int) floor((theta)/10.0); 
                  mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                  sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
              }
                  else
                     {
                         k = Marine20ROWS - 2;
                         mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/100.0; // why 100?
                         sig =CC(k,2);
                     }
        }       
                return mu;
        }



	double U( double y, double vr, double theta)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}


	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}

};




//////////////
// SHCal20
//////////////

class SHCal20 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	int Bomb;
	Cal *bombcc;
	char name[255];
	double mincal, const2;

public:

	SHCal20(int bomb, std::string ccdir) : Cal(SHCal20ROWS) {

	CCB = new Matrix( SHCal20ROWS, SHCal20COLS);
	CC.Set( CCB, CCB->nRow(), CCB->nCol());

	std::string filename= ccdir+SHCal20FNAM;

	Rprintf("SHCal20: Reading from file: %s\n", filename.c_str());


	if (CC.filescan((char*)filename.c_str()) == 0) {
		REprintf("Cal: ERROR: Could not find SHCal20 cal. curve, file not found: %s\n", filename.c_str());
		Rcpp::stop("Cal: ERROR: Could not find SHCal20 cal. curve, file not found: %s\n", filename.c_str());
		//  exit(0);
	}

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library

		const char *postbombfnam[] = { POSTBOMBFNAMS };

		/** Read bomb **/
		Bomb = bomb;
		if (Bomb == 0) {
			mincal = 0.0; // no bomb; 17 Dec 2018 changed -5.0 to 0.0
			//sprintf( name, "SHCal20");
			snprintf( name, sizeof(name), "SHCal20");
		}
		else
			if (Bomb < 6) { // curve number, not cal BP yr. 26 March 2019: Was Bomb < 5 but now there are 5 postbomb curves

			bombcc = new GenericCal(postbombfnam[Bomb], ccdir);
			mincal = bombcc->MinCal();
			//sprintf( name, "SHCal20+%s", postbombfnam[Bomb]);
			snprintf( name, sizeof(name), "SHCal20+%s", postbombfnam[Bomb]);
			}
			else {
				REprintf("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
				Rcpp::stop("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
//				exit(0);
			}
	}

	~SHCal20() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return mincal; }
	double MaxCal() { return 55000.0; } // was 50000.0

	virtual const char *Name() { return name; } //JEV warning


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
		if (fcmp(theta, -0.0) == -1) // value < 0 cal BP, postbomb
		{
			if (Bomb == 0) {
				//fprintf( stderr, "WARNING: Calibration attempted beyond IntCal13 cal. curve limits, theta= %f\n",theta);
				k = 0; //Extrapolation
				mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
				sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
			}
			else {
				bombcc->cal(theta);
				mu = bombcc->GetMu();
				sig = bombcc->GetSig();
			}

		}
		else { // until 5,000 cal BP, line 4999 (check!), values every year
			if (fcmp(theta, 5000.0) != 1) // value < 5,000 cal BP
				{
					k = 0 + (int) floor(theta/1.0); // 0 at start was 1
					mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/1.0;
					sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/1.0;
				}
				else // from 5,000 until 15,000 cal BP, line 7000 (check!), values every 5 years
					if (fcmp(theta, 15000.0) != 1)
						{
							k = 4999 + (int) floor((theta-5000.0)/5.0);
							mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
							sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
						}
						else // from 15,000 until until 25,000 cal BP, line 8000 (check!), values every 10 years
							if (fcmp(theta, 25000.0) != 1)
								{
									k = 7000 + (int) floor((theta-15000.0)/10.0);
									mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
									sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
								}
								else // from 25,000 until 50,000 cal BP, line 9250, values every 20 years
									if (fcmp(theta, 50000.0) != 1)
										{
											k = 8000 + (int) floor((theta-25000.0)/20.0);
											mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
											sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
										}
									//	else // value > 55,000 cal BP, extrapolate
									//		{
									//			k = SHCal20ROWS - 2;
									//			mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/100.0;
									//			sig = CC(k,2);
									//		}
										else // value > 50,000 cal BP, linear "estimated" cc
											{
                                                //The  cc is a straight line from the value at
                                                //50000 to:
                                                double last_th  = 100000.0;
                                                double last_mu  = 95840.0;
                                                double last_sig = 10000.0;
												k = 9250; //Row at 50000 cal BP
                                                mu = CC(k,1) + (theta-CC(k,0))*(last_mu-CC(k,1))/(last_th-CC(k,0));
                                                sig = CC(k,2) + (theta-CC(k,0))*(last_sig-CC(k,2))/(last_th-CC(k,0));
											}

		}
		return mu;
	}

	double U( double y, double vr, double theta)
	{
		cal(theta);

		double tau = 1.0/(vr + sqr(sig));

		return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
		cal(theta);

		double tau = 1.0/(vr + sqr(sig));

		return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}
};



/* IntCal13 */
class IntCal13 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	int Bomb;
	Cal *bombcc;
	char name[255];
	double mincal, const2;

public:

    IntCal13(int bomb, std::string ccdir) : Cal(IntCal13ROWS) {

		CCB = new Matrix( IntCal13ROWS, IntCal13COLS);
		CC.Set( CCB, CCB->nRow(), CCB->nCol());

        std::string filename= ccdir+IntCal13FNAM;

        Rprintf("IntCal13: Reading from file: %s\n", filename.c_str());


        if (CC.filescan((char*)filename.c_str()) == 0) {
            REprintf("Cal: ERROR: Could not find IntCal13 cal. curve, file not found: %s\n", filename.c_str());
            Rcpp::stop("Cal: ERROR: Could not find IntCal13 cal. curve, file not found: %s\n", filename.c_str());
      //  exit(0);
    }

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library

        const char *postbombfnam[] = { POSTBOMBFNAMS };

		/** Read bomb **/
		Bomb = bomb;
		if (Bomb == 0) {
			mincal = 0.0; // no bomb; 17 Dec 2018 changed -5.0 to 0.0
			//sprintf( name, "IntCal13");
			snprintf( name, sizeof(name), "IntCal13");
		}
		else
			if (Bomb < 6) { // curve number, not cal BP yr. 26 March 2019: Was Bomb < 5 but now there are 5 postbomb curves

            bombcc = new GenericCal(postbombfnam[Bomb], ccdir);
			mincal = bombcc->MinCal();
			//sprintf( name, "IntCal13+%s", postbombfnam[Bomb]);
			snprintf( name, sizeof(name), "IntCal13+%s", postbombfnam[Bomb]);
			}
			else {
                REprintf("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
                Rcpp::stop("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
//				exit(0);
			}


	}

	~IntCal13() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return mincal; }
	double MaxCal() { return 50000.0; }

    virtual  const char *Name() { return name; } //JEV warning


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
        if (fcmp(theta, -0.0) == -1) // 17 Dec 2018 changed -5.0 to 0.0
        {
			if (Bomb == 0) {
				//fprintf( stderr, "WARNING: Calibration attempted beyond IntCal13 cal. curve limits, theta= %f\n",theta);
				k = 0; //Extrapolation
				mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
				sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
			}
			else {
				bombcc->cal(theta);
				mu = bombcc->GetMu();
				sig = bombcc->GetSig();
			}

        }
        else {
            if (fcmp(theta, 13900.0) != 1)
                {
                        k = 0 + (int) floor(theta/5.0); // 0 was 1
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
                        if (fcmp(theta, 25000.0) != 1)
                        {
                                k = 2780 + (int) floor((theta-13900.0)/10.0); // line 2781 (in c-speak 2780) is 13900 calBP - 10yr steps after this
                                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                                sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
                        }
                        else
                                if (fcmp(theta, 50000.0) != 1) // line 3891 (3890 in c-speak) is 25000 cal BP - 20yr steps after this
                                {
                                        k = 3890 + (int) floor((theta-25000.0)/20.0);
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
                                }
					else
						{
                        //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal13 cal. curve limits, theta= %f\n",theta);
							k = IntCal13ROWS - 2;
							mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/100.0;
							sig = CC(k,2);
						}
        }
            return mu;
	}

	double U( double y, double vr, double theta)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}
};





class Marine13 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	double const2;

public:

    Marine13(std::string ccdir) : Cal(Marine13ROWS) {

		CCB = new Matrix( Marine13ROWS, Marine13COLS);

		CC.Set( CCB, CCB->nRow(), CCB->nCol());

        std::string filename= std::string(ccdir)+std::string(Marine13FNAM);

        Rprintf("Marine13: Reading from file: %s\n", filename.c_str());

        if (CC.filescan((char*)filename.c_str()) == 0) {
           REprintf("Cal: ERROR: Could not find Marine13 cal. curve, file not found: %s\n", filename.c_str());
           Rcpp::stop("Cal: ERROR: Could not find Marine13 cal. curve, file not found: %s\n", filename.c_str());

//			exit(0);
		}

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library


	}

	~Marine13() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return 0.0; }
	double MaxCal() { return 50000.0; }


	//char *Name() { return "IntCal09Marine"; }
 const char* Name() { return "Marine13"; }


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
        if (fcmp(theta, 0.0) == -1)
        {
                //fprintf( stderr, "WARNING: Calibration attempted beyond marine13 cal. curve limits, theta= %f\n",theta);
                k = 0;
                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5;
                //sig <- CC(k,2); before JUDY
                 sig = CC(k,2);
        }
        else
        {       // MB added '{' on 30 Jan 2019 to avoid indentation warning while compiling on Fedora. This '{' was placed correctly for IntCal13 (line 318)
                if (fcmp(theta, 10500.0) != 1)
                {       //************** NB: 0 is the number of rows before cal year 0
                        k = 0 + (int) floor(theta/5.0);
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
                        if (fcmp(theta, 25000.0) != 1)
                        {
                                k = 2100 + (int) floor((theta-10500.0)/10.0); // steps of 10yr after 10500 (line 2100)
                                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                                sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
                        }
                        else
                                if (fcmp(theta, 50000.0) != 1)
                                {
                                        k = 3550 + (int) floor((theta-25000.0)/20.0); // steps of 20yr after 25000 (line 3550)
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
                                }
                                else
                                {
                        //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal13 cal. curve limits, theta= %f\n",theta);
                                        k = Marine13ROWS - 2;
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/100.0;
                                        sig =CC(k,2);
                                }
        }       // MB added '}' on 30 Jan 2019 to avoid indentation warning while compiling on Fedora
                return mu;
        }



	double U( double y, double vr, double theta)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}


	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}

};




class SHCal13 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	SubMatrix A;
	int Bomb;
	Cal *bombcc;
	char name[255];
	double mincal, const2;

public:

    SHCal13(int bomb, std::string ccdir) : Cal(SHCal13ROWS) {

		CCB = new Matrix( SHCal13ROWS, SHCal13COLS);
		CC.Set( CCB, CCB->nRow(), CCB->nCol());
        std::string filename= std::string(ccdir)+SHCal13FNAM;

        Rprintf("SHCal13: Reading from file: %s\n", filename.c_str());

            if (CC.filescan((char*)filename.c_str()) == 0) {
                REprintf("Cal: ERROR: Could not find SHCal13 cal. curve, file not found: %s\n", filename.c_str());
                Rcpp::stop("Cal: ERROR: Could not find SHCal13 cal. curve, file not found: %s\n", filename.c_str());
//			exit(0);
		}

		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library

        const char * postbombfnam[] = { POSTBOMBFNAMS };

		/** Read bomb **/
		Bomb = bomb;
		if (Bomb == 0) {
			mincal = -0.0; // no bomb, was -5.0 changed 17 Dec 2018
			//sprintf( name, "SHCal13");
			snprintf( name, sizeof(name), "SHCal13");
		}
		else
			if (Bomb < 6) { /// was Bomb < 5 but now there are 5 postbomb curves (March 2019)

            bombcc = new GenericCal(postbombfnam[Bomb],(char*)ccdir.c_str());
			mincal = bombcc->MinCal();
			//sprintf( name, "SHCal13+%s", postbombfnam[Bomb]);
			snprintf( name, sizeof(name), "SHCal13+%s", postbombfnam[Bomb]);
			}
			else {
                REprintf("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
                Rcpp::stop("Bacon: ERROR: Post bomb curve: 0 None, 1 NH1, 2 NH2, 3 NH3, 4 SH1-2, 5 SH3\n");
            //	exit(0);
			}


	}

	~SHCal13() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}


	double MinCal() { return mincal; }
	double MaxCal() { return 50000.0; }

    virtual const char *Name() { return name; } //JEV warning


//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001);
	double cal(double theta)
	{
        if (fcmp(theta, 0.0) == -1)
        {
			if (Bomb == 0) {
				//fprintf( stderr, "WARNING: Calibration attempted beyond SHCal13 cal. curve limits, theta= %f\n",theta);
				k = 0; //Extrapolation
				mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
				sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
			}
			else {
				bombcc->cal(theta);
				mu = bombcc->GetMu();
				sig = bombcc->GetSig();
			}

        }
        else {
                if (fcmp(theta, 13900.0) != 1)
                {				//************** NB: In the official SHCal13 year 0 is node 1 (node 0 is -5)
                        k = 0 + (int) floor(theta/5.0); // 0 was 1
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
                        if (fcmp(theta, 25000.0) != 1) // 10yr steps from 13900 on - line 2781 (2780 in c speak)
                        {
                                k = 2780 + (int) floor((theta-13900.0)/10.0);
                                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                                sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
                        }
                        else
                                if (fcmp(theta, 50000.0) != 1) // 20yr steps from 25000 on - line 3891 (3890 in c speak)
                                {
                                        k = 3890 + (int) floor((theta-25000.0)/20.0);
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig = CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
                                }
					else
						{
                        //fprintf( stderr, "WARNING: Calibration attempted beyond SHCal13 cal. curve limits, theta= %f\n",theta);
							k = SHCal13ROWS - 2;
							mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/100.0;
							sig = CC(k,2);
						}
        }
		return mu;
	}

	double U( double y, double vr, double theta)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)
	{
        cal(theta);

        double tau = 1.0/(vr + sqr(sig));

        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}
};


class Plum : public Cal {

	double alPhi, mPhi, alS, mS; //a priori pars form phi and PS
	double Al, theta0; //Detection limit.
	int radon;
	int nS; //number of supported observations
	Matrix *SB;
	SubMatrix S; //To hold supported data

public:
	// Parameters for the apriori for \Phi and P^S .  radon=0 No radon, radon =1 Radon case (A) and radon =2 Radon case (B).
	// and fnam the file name for the supportted data
    Plum( double alPhi_, double mPhi_, double alS_, double mS_, double Al_, double theta0_, int radon_, const char *fnam, std::string ccdir ) : Cal(0) {
      Rprintf("Calibration 'curve' used to handle 210Pb data (Plum).\n");
      //printf("Calibration 'curve' used to handle 210Pb data (Plum).\n");
		  alPhi = alPhi_;
  		mPhi = mPhi_;
  		alS = alS_;
  		mS = mS_;
  		Al = Al_;
  		radon = radon_;
      theta0 = theta0_;

      //printf("alPhi %lf mPhi %lf alS %lf ms %lf Al %lf radon %d theta0 %lf\n", alPhi, mPhi, alS, mS, Al, radon, theta0);

  		// Read the supported data file, two colouns only no header, y^S_i and s_i
      //std::string filename = std::string(ccdir)+std::string(fnam);
      std::string filename = std::string(fnam);


  		//Count the number of lines in the file
  		FILE *F;
      if ((F = fopen( filename.c_str(), "r")) == NULL) {
        REprintf("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        //printf("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        Rcpp::stop("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        //exit(0);
      }

      Rprintf("Supported data file %s\n", filename.c_str());

  		int numrows=0;
  		char ln[GENCCMAXLINLEN];
  		while (!feof(F)) {
        char* result=fgets(ln, GENCCMAXLINLEN, F); //JEV warning
        if(result==NULL){}//JEV warning
        numrows++;
		  }
		  numrows--;
		  fclose(F);

  		//Read the file as usual:
  		SB = new Matrix( numrows, 2);

  		S.Set( SB, SB->nRow(), SB->nCol());

      Rprintf("Plum: Reading supported data from file: %s, %d rows, 2 cols.\n", filename.c_str(), numrows);
      //printf("Plum: Reading supported data from file: %s, %d rows, 2 cols.\n", filename.c_str(), numrows);

      if (S.filescan((char*)filename.c_str()) == 0) {
        REprintf("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        //printf("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        Rcpp::stop("Plum: ERROR: Could not find supported data, file not found: %s\n", filename.c_str());
        //exit(0);
		  }



		  nS = S.nRow();
		  for (int j=0; j < nS; j++) {
        //printf("S(%d, %d) = %lf\n", j, 1, S(j,1));
        S( j, 1) = S( j, 1)*S( j, 1); //store the variances
		  }



	}

	double cal(double theta) {
		mu = theta; //The identity
		sig = 0.0;
		return theta;
	}

  const char* Name() { return "Plum"; }

	//Here is simply the identity
	double U( double y, double vr, double theta) { return 0.5*sqr(y - theta)/vr; }
	double Ut( double y, double vr, double theta, double a, double b) { return (a + 0.5)*log( b + 0.5*sqr(y - theta)/vr); }

	double MinCal() { return -1.0e300; }
	double MaxCal() { return 1.0e300; }

	int NS() { return nS; }
	double yS(int j) { return S( j, 0); } // the supported data
	double s2(int j) { return S( j, 1); }
	int MorePars() {
		if (radon < 2)
			return 2; // No Radon or Radon case (A), Phi and PS are the additional parameters
		else
			return 1 + nS; //Radon case (B), Phi and n_S additional parameters
	}
	double GetAl() { return Al; }
  double GetTheta0() { return theta0; }
  double GetAlPhi() { return alPhi; }
  double GetMPhi() { return mPhi; }
  double GetAlS() { return alS; }
  double GetMS() { return mS; }
};


/******* Class Det *******/


//Definition for class Det, to hold a generic "determination", could be non-radiocarbon
class Det {
protected:
	char *nm;	// lab number
	double y;		// mean
	double std;	// std dev

	double x;  //depth or anything else

	double deltaR;	// delta-R (reservoir correction)
	double deltaSTD;	// delta-R std dev. (reservoir correction)

    int is_210Pb; //flag to verify if this is a 210Pb datum is_201Pb=1 if it is a 210Pb datum
	double rho, delta; //rho (density) for 210Pb, delta (thickness)

	double a, b; //prior parameters for the t distribution

	Cal *cc;	// calibration curve to use

	double med;	// mean and variance and std dev. after reservoir correction
	double vr;
	double corrstd;

public:
	Det(char *enm, double ey, double estd, double xx, double edeltaR, double edeltaSTD, double ea, double eb, Cal *ecc) {

		//Read members
		nm = strdup(enm);
		y = ey;
		std = estd;
		x = xx;
		deltaR = edeltaR;
		deltaSTD = edeltaSTD;

		a = ea;
		b = eb;

		cc = ecc;

        if (strcmp( cc->Name(), "Plum") == 0) { //It is a 210Pb datum
            is_210Pb = 1;
            delta = deltaR;
            rho = deltaSTD;
            std = std*rho;
            med = y*rho;
            vr = sqr(std);
            corrstd = std;
        } else {
            is_210Pb = 0;
            rho = 0.0;
            delta = 0.0;
            med = y - deltaR;
            vr = sqr(std) + sqr(deltaSTD);
            corrstd = sqrt(vr);
        }

	}

	virtual void ShortOut() {
        Rprintf("%s: %6.1f+-%-6.1f d=%-g ResCorr=%6.1f+-%-6.1f a=%-g b=%-g cc=%s\n",
				nm, y, std, x, deltaR, deltaSTD, a, b, cc->Name());
	}

	double ChangeCorrMean(double y) { return (med = y); }
	double ChangeSd(double stdnew) {
		std = stdnew;
		vr = sqr(std) + sqr(deltaSTD);
		return std;
	}

	const char *labnm() { return nm; }
	double mean() { return y; }
	double sd() { return std; }
	double corr_mean() { return med; }
	double corr_vr() { return vr; }
	double res_mean() { return deltaR; }
	double res_std() { return deltaSTD; }
	double d() { return x; }

    int Is210Pb() { return is_210Pb; }
	double Delta210Pb() { return delta; }
	double Rho210Pb() { return rho; }
	const Cal *GetCC() { return cc; }

	//exp(-U) will be the likelihood for this determination.
	virtual double U(double theta) { return cc->U( med, vr, theta); }
	virtual double Ut(double theta) { return cc->Ut( med, vr, theta, a, b); }
};

class DetCensor : public Det {
    
public:
	DetCensor(char *enm, double ey, double estd, double xx, double edeltaR, double edeltaSTD, double ea, double eb, Cal *ecc) : Det( enm, ey, estd, xx, edeltaR, edeltaSTD, ea, eb, ecc) { }

	void ShortOut() {
        Rprintf("%s (censored later): %6.1f+-%-6.1f d=%-g ResCorr=%6.1f+-%-6.1f a=%-g b=%-g cc=%s\n",
				nm, y, std, x, deltaR, deltaSTD, a, b, cc->Name());
	}

	double U(double theta) {
        cc->cal(theta);

        double sigma = sqrt(vr + sqr(cc->GetSig()));
        
        return -log(1.0-NorF((y - cc->GetMu())/sigma));
    }
	double Ut(double theta) { 
        /* Not yet implemented, the stydent-t cdf:
        return -log(1.0-StTF((y - cc->GetMu())/cc->GetSig())); }
        we will use a Gaussian as above */
        return U(theta);
    }
};


class DetCensorE : public Det {
    
public:
	DetCensorE(char *enm, double ey, double estd, double xx, double edeltaR, double edeltaSTD, double ea, double eb, Cal *ecc) : Det( enm, ey, estd, xx, edeltaR, edeltaSTD, ea, eb, ecc) { }

	void ShortOut() {
        Rprintf("%s (censored earlier): %6.1f+-%-6.1f d=%-g ResCorr=%6.1f+-%-6.1f a=%-g b=%-g cc=%s\n",
				nm, y, std, x, deltaR, deltaSTD, a, b, cc->Name());
	}

	double U(double theta) {
        cc->cal(theta);

        double sigma = sqrt(vr + sqr(cc->GetSig()));
        
        return -log(NorF((y - cc->GetMu())/sigma));
    }
	double Ut(double theta) { 
        /* Not yet implemented, the stydent-t cdf:
        return -log(1.0-StTF((y - cc->GetMu())/cc->GetSig())); }
        we will use a Gaussian as above */
        return U(theta);
    }
};



//To hold a series of determinations
class Dets {

protected:

	Det **det;  //array with the determinations

	int m; //current number of determinations
	int max_m; //Maximum number of determinations

public:
	//The constructor only opens an array of pointers to Det structures
	Dets(int emax_m) {

		max_m = emax_m;

		m = 0;

		det = new Det * [max_m];
	}

	void AddDet(Det *de) {

		if (m == max_m) {

            REprintf("ERROR: Maximum number of determinations exceeded\n\n");
            Rcpp::stop("ERROR: Maximum number of determinations exceeded\n\n");

//			exit(0);
		}

		m++;
		det[m-1] = de;

        Rprintf("Added det: ");
		ShortOut(m-1);
	}

	int Size() { return m; }

	void ShortOut(int j) { det[j]->ShortOut(); }

	const char *labnm(int j) { return det[j]->labnm(); }
	double y(int j)  {  return det[j]->corr_mean(); }
	double sd(int j) {  return det[j]->sd(); }
	double vr(int j) {  return det[j]->corr_vr(); }
	double d(int j)  {  return det[j]->d(); }


	double SetY(int j, double y) {  return det[j]->ChangeCorrMean(y); }
	double SetSd(int j, double y) {  return det[j]->ChangeSd(y); }

  int Is210Pb(int j)  {  return det[j]->Is210Pb(); }
	double Rho210Pb(int j) { return det[j]->Rho210Pb(); }
	double Delta210Pb(int j) { return det[j]->Delta210Pb(); }
	const Cal *GetCC(int j) { return det[j]->GetCC();}

	double U(int j, double theta)  { return det[j]->U(theta); }
	double Ut(int j, double theta) { return det[j]->Ut(theta); }
};

#endif
