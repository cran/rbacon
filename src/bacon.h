/*
*********

Andres christen, May 2009.

Bacon

 BaconFix to be used with large K, number of sections

 BaconMov with moving borders, to be used with small K (<10). Not used any more (comment MB April 2019)
 NEW VERSION OF BACON: With Plum and USING THE ALPHA'S AS PARAMETERS
 instead of the x's ... that is: x_j = w x_{j+1} + (1-w) alpha_j

 Use the alpha_j's as parameter and then transform to x_j ... this is the way Nico and Maroc does it.

 The array X is is now used to communicate with the twalk.  This is then transalated to x (and thetas) in
 the SetThetas function, that now uses the object's variable x.  Everything remains the same afterwards.

 **/



#ifndef BACON_H
#define BACON_H

//This is my traditional farewell, may be changed to something more "serious"
#define FAREWELL "Eso es to...eso es to...eso es to...eso es toooodo amigos!\n"

//#include <stdio.h>
#include <math.h>
//#include <unistd.h>
#include <string.h>

#include "cal.h"
#include "ranfun.h"
#include "Matrix.h"
#include "twalk.h"              // twalk simulator class



#define CHARBUFFER 8000

#define LA_CONST 0.03114


class Bacon: public obj_fcn {

	public:
		Bacon(int dim) : obj_fcn(dim) { /*nothing, template class*/}

        void show_descrip() const { Rprintf("Bacon:\n"); }

		virtual double *Getx0() = 0;
		virtual double *Getxp0() = 0;

		virtual double Getc0() = 0;
		virtual double GetcK() = 0;
		virtual void ShowDescrip() = 0;
		virtual void PrintNumWarnings() = 0;

};


//Fixed number of sections
class BaconFix: public Bacon {
       protected:

			//Object that holds all determinations
			Dets *dets;

			int m, K; //m number of dets, K number of sections
			int H; //number of hiatuses
			double *h; //location of the hiatuses

			int useT; //=1 to use the t model, =0 to use the noermal model

			double w, w0, wp0;

			double *x, *X0, *Xp0, *theta;

			double MinYr, MaxYr;
			double MaxYrTheta0Plum;

			//Based depth and increment between depths
			double c0, Dc;
			virtual double c(int i) { return c0 + i*Dc; }

			double U, Uprior, Uli;
            void AccPars(int prime) { /*fprintf( F, "%f  %f  %f\n", Uprior, Uli, U);*/
            prime=0;}

			double *alpha, *beta; //prior pars for the acc gamma prior in each inter hiatus section
			double prioracU(int i, const double al) { return (1.0-alpha[i])*log(al) + beta[i]*al; }

			double priorPhiU(double *x) {
				double scale_fi = plumobj->GetMPhi()/plumobj->GetAlPhi();
				double shapefi = plumobj->GetAlPhi();
				//prior = prior -  (  (shapefi-1.)*log(param[0])-(param[0]/scale_fi) ) NOTE: prior for fi
				return  ( 1.0 - shapefi)*log( x[K+2] ) + ( x[K+2] / scale_fi ) ;
				//printf("Valor de Fi %lf\n",x[K+2] );
			}

			double priorPSU(double *x){
				double priorU = 0.0;
				double shapeAs = plumobj->GetAlS();
				double scale_As = plumobj->GetMS()/plumobj->GetAlS();
				double *PS = x + ( K+3 );
				for (int i = 0; i < GetnPs() ; i++) {
					priorU += (1.0 - shapeAs)*log( PS[i] ) + (PS[i]/scale_As );
					//printf("Valor de Ps %lf\n",PS[i]);
				}
				//priorU=1000*priorU;
				//for k in range(Ran):
				//	prior= prior -  (  (shapeAS-1.)*log(param[1+k])-(param[1+k]/scale_As) )NOTE: prior for supp

				return priorU;
			}

			double a, b; //a priori pars for the w beta prior
			double ds;
			double rsc, logrsc, logw;
			// ds=1.0, rsc=ds/Dc and logrsc=log(ds/Dc) set in the creator, lines 177 and 178 after reading Dc
			double priorwU(const double w) {
				logw=log(w);
			  return rsc*(1.0-a)*logw + (1.0-b)*log(1.0-exp(rsc*logw)) + (1.0-rsc)*logw - logrsc;

				//rsc = ds/Dc = 1/[(cm-c0)/K]
				//( ((1./by)-1.)*log(w)- log(by)+  ((1./by)*(shape1_m-1.))*log(w) +
				//(shape2_m - 1.)*log(1.-w^(1./by) ) )# prior for w	#

				//the last term (1.0-rsc)*logw - logrsc was missing,
				//see f(w), p.461, of the paper, jac: changed 22OCT2018
			}
			//here ds = 1.0 (in your depth units), it could be changed to a parameter


			double *ha, *hb; //a priori pars for the uniform prior on hiatus jumps in each inter hiatus.
//H Change			double priorHU(int i, const double x) { return (1.0-ha[i])*log(x) + hb[i]*Dc*x; }
			double priorHU(int i, const double x) { return 1.0; } //Uniform

			int WarnBeyondLimits;
			//Sets the thetas and verifies correct limits
			int SetThetas(double *X) {
				double S=X[0]; //th0
				double w = X[K+1];
				int rt = 1;

				for (int k=get_dim()-1; k>K; k--)
					x[k] = X[k]; //Copy all (Plum) parameters and w = X[K]
				x[0] = X[0]; //and th0
				theta[0] = x[0];

				x[K]  = X[K]; //= alpha[K]
				if (H == 0) {  //with no hiatus
				//if( true ){
					for (int k=K-1; k>0; k--) {
						x[k]  = w*x[k+1] + (1.0-w)*X[k]; //Create the x's
					}
				} else {

					//we go backwards until we find the hiatus
					int l=0;
					for (int k=K-1; k>0; k--) {
						if ((fcmp( c(k-1), h[l]) == -1) && (fcmp( h[l], c(k)) != 1)) { //forgets
							x[k] = X[k];

							l++; //jump to next hiatus, but max one hiatus in each section.
						}
						else //continue with  memory
							x[k] = w*x[k+1] + (1.0-w)*X[k];
					}

				}

				//Create the thetas, this is the same as in the old version, once x has been created:
				if ( (fcmp( theta[0], MinYr) == -1) || (fcmp( theta[0], MaxYrTheta0Plum ) == 1) ) {// [MinYr, MaxYrTheta0Plum]
					WarnBeyondLimits++;

					//beyond established limits
					rt = 0;
				}
				for (int k=1; k<K; k++) {
					S += x[k]*(c(k)-c(k-1)); //For fixed c's, Dc = c(k)-c(k-1)
					theta[k] = S;
				}
				//Last theta
				theta[K] = theta[K-1] + x[K]*(c(K)-c(K-1));
				if (fcmp( theta[K], MaxYr) == 1)
					WarnBeyondLimits++;
					//beyond established limits



				return rt;
			}

			int plumUsed, nPS, last210Pb; //=1 if PLum needs to be used, number of PS's parameters and number of the last 210Pb datum
			double /**PS,*/ phi;
			Plum *plumobj; //Pointer to Plum class

       public:
	         BaconFix( Dets *detsdets, int KK, int HH, double **hiatus_pars, double aa, double bb,
				double MMinYr, double MMaxYr, double th0, double thp0, double cc0, double cm, int uuseT,
				unsigned long int seed, int more_pars=0)
				: Bacon((KK+1) + 1 + more_pars) {


				//Use the student t model, 1, or the traditional normal model, 0
				useT = uuseT;

				dets = detsdets;

				m = dets->Size();
				//Minimum and maximum years
				MinYr = MMinYr;
				MaxYr = MMaxYr;

				WarnBeyondLimits = 0;

				K = KK; // Number of sections
				H = HH; // Number of hiatuses
				//   reg.   w

				a = aa;
				b = bb;

				//hiatus_pars is a pointer to 5 arrays of doubles of size H+1
				//containing:
				h     = hiatus_pars[0];		// hiatus depth(s)

				alpha = hiatus_pars[1];		// acc.shape
				beta  = hiatus_pars[2];		// acc.shape/acc.mean

				ha    = hiatus_pars[3]; 	// WAS hiatus.shape, now a dummy
				hb    = hiatus_pars[4]; 	// WAS hiatus.mean, now hiatus.max


				//Open memory for the two points in the parameter space
				//IN THE NEW BACON:
				//X[0] is theta[0], then alpha[1] ... alpha[K-], alpha[K]=x[K] and X[K+1]=w
				X0  = new double[get_dim()];
				Xp0 = new double[get_dim()];

				//AND FOR THE NEW BACON, x is now not passed to the twalk, but is a local variable:
				x = new double[get_dim()];

				//In the initial point, x will be translated from X0

				//These will hold the cal. years at each node
				//The translation from x's to thetas's is done in method insupport
				theta = new double[K+1];

				//Set the sections, locations for the c's
				c0 = cc0;
				Dc = (cm-c0)/(double) K;
				ds=1.0;
				rsc=ds/Dc;
				logrsc=log(ds/Dc);


				//Verify the ordering in the h's // disabled MB 13 May 2019 - JAC OK
				//The h's must be an array of size H+1!!! although there are only H hiatuses
				//		for (int k=0; k<H; k++) {
				//			if (fcmp( h[k], ((k == 0) ? c(K) :  h[k-1]) - Dc) != -1) { //we need only one per section
        		//                //REprintf("Bacon: ERROR: The hiatuses are not in descending order and/or less than %f\n", c(K));

				//				exit(0); //h[k] not in the correct order ... we need to have h[H-1] < ... < h[0] < c(K)
        		//				Rcpp::stop("Bacon: ERROR: The hiatuses are not in descending order and/or less than %f\n", //c(K)); // commented MB 11 May 2019
				//			}
				//		}

				if (H > 0)
					if (fcmp( h[H-1], c0) == -1) {
                        REprintf("Bacon: ERROR: The last hiatus location is not greater than %f\n", c0);

                    //	exit(0); //we need to have c0 < h[H-1]
                        Rcpp::stop("Bacon: ERROR: The last hiatus location is not greater than %f\n", c0);

					}
				h[H] = c0 - 2*Dc; //fix h[H] lower than the low limit for depths


				//Initial values for x0
				X0[0]  = th0;
				x[0] = X0[0];
				Xp0[0] = thp0;

				//Rprintf("Plum: Seed %d\n", seed);

				Seed(seed); //Set the Seed for random number generation

				//and for w, from its prior
				X0[K+1]  = BetaSim( a, b);
				x[K+1] = X0[K+1];
				Xp0[K+1] = BetaSim( a, b);
				w0 = X0[K+1];
				wp0 = Xp0[K+1]; //short names for the initial values

				//******************* NB ************************
				//the prior is scale=1/beta[0] ... however, to avoid models growing out of bounds
				//we prefer higher accumulation rates: scale=mult/beta[0]
				double mult=1.0;

				//initial values for the acc. rates
				X0[K]  = GammaSim( alpha[H], 1.0/beta[H]);
				Xp0[K] = GammaSim( alpha[H], 1.0/beta[H]);
				//x[K] = X0[K];
				if (H == 0) {  //with no hiatus
					for (int k=K-1; k>0; k--) {
						X0[k] = GammaSim( alpha[0], mult/beta[0]); //alpha[k]
						//x[k]  = w0*x[k+1] + (1.0-w0)*X0[k];
						Xp0[k] = GammaSim( alpha[0], mult/beta[0]);

					}
				} else {//initial values for the acc. rates, with hiatus

					//we go backwards until we find the hiatus
					int l=0;
					for (int k=K-1; k>0; k--) {
						if ((fcmp( c(k-1), h[l]) == -1) && (fcmp( h[l], c(k)) != 1)) { //if c_{k-1} < h_l & h_l !> c_k, forgets
							X0[k]  = GammaSim( ha[l], 1.0/(hb[l]*Dc) );
							//x0[k]  = GammaSim( alpha[l], mult/(beta[l]) ); // MB May 2019
							//x[k] = X0[k];
							l++; //jump to next hiatus, but max one hiatus in each section.
						} else { //continue with the memory
							X0[k]  = GammaSim( alpha[l], mult/beta[l]);
							//x[k] = w0*x[k+1] + (1.0-w0)*X0[k];
						}

					}

					l = 0; //do it again
					for (int k=K-1; k>0; k--) {
						if ((fcmp( c(k-1), h[l]) == -1) && (fcmp( h[l], c(k)) != 1)) { //forgets
							Xp0[k]  = GammaSim( ha[l], 1.0/(hb[l]*Dc) );
							//xp0[k]  = GammaSim( alpha[l], mult/(beta[l]) ); // MB Apr 2019
							l++; //jump to next hiatus, but max one hiatus in each section.
						} else{ //continue with the memory
							Xp0[k]  = GammaSim( alpha[l], mult/beta[l]);
						}
					}

				}

				if (more_pars != 0) { //Plum needs to be used!!!!
 				 plumUsed = 1;
 				 for (int j=0; j<m; j++)
 				 	 if (dets->Is210Pb(j) == 1) {
 						 plumobj = (Plum*) dets->GetCC(j); //Gets hold of the pointer to the Plum object in use.
 						 last210Pb = j;
 					 }
 				 nPS = more_pars - 1;

 				 /***initial values for phi **/

 				 double limitPhi=0.0;

				 for (int k=K; k>0; k--) {
					 X0[k]  = X0[k]*0.3;
					 Xp0[k] = Xp0[k]*0.3;
				 }

				 //limitPhi is the smallest value of phi
 				 SetThetas(Xp0); //Creates x from Xp0
 				 limitPhi = LA_CONST*plumobj->GetAl() * exp( LA_CONST * (G( dets->d(last210Pb), x) - x[0])  );
 				 Xp0[K+2] = 1.2*limitPhi + (0.8*limitPhi)*Un01(); //phi Random values between [1.2*limitPhi,2*limitPhi]

 				 SetThetas(X0); //Creates back x from X0
         limitPhi = LA_CONST*plumobj->GetAl() * exp( LA_CONST * (G( dets->d(last210Pb), x) - x[0]));
 				 X0[K+2]  = 1.2*limitPhi + (0.8*limitPhi)*Un01(); //phi Random values between [1.2*limitPhi,2*limitPhi]

				 //Initial values for PS
 				 for (int j = 0; j < nPS; j++) { //Random value
 					 X0[K+3+j]  =  Un01()*15.0;        //[0,15]
 					 Xp0[K+3+j] =  Un01()*15.0 + 15.0; //[15,30]
 				 }

 				 MaxYrTheta0Plum = MinYr + 0.04;

 			 } else {
 				 plumUsed = 0;
 				 plumobj = NULL;
 				 nPS = 0;
 				 phi = -1.0;
 				 MaxYrTheta0Plum = MaxYr;
 			 }
			}


			//Return the value of the PS parameter, in the 210Pb Plum dating form the vector of pars x.
			//In case 0 and 1 (nPS = 1), it is fixed to x[K+3], otherwise is x[K+3 + j]
			//It is called many times so we better use it inline
			inline double GetPS(int j, double *x) {
				if (nPS == 1)
					return x[K+3];
				else
					return x[K+3 + j];
			}

			double GetnPs(){
				return nPS;
			}

			//x[0] is theta[0], then x[1] ... x[K], x[K+1]=w
			//phi is x[K+2] in plum
			//x[K+3] ... x[K+nPS] is PS for support data in plum
	    int insupport(double *X) {


				//NOTE: Check the support for PS and phi
				if (plumUsed == 1) {

					for (int j=0; j<nPS; j++) {
						//printf("PS[%d]=%lf\n", j, x[K+3 + j]);
					 	//PS[j] = x[K+3 + j];
						if   (fcmp( X[K+3 + j], 0.0) != 1){  //PS out of support
							//Rprintf("Plum 1: PS out of support\n");
							return 0;
						}
					}

					phi = X[K+2];
					if   (fcmp( phi, 0.0) != 1){  //phi out of support
						//Rprintf("Plum 2: phi out of support\n");
						return 0;
					}

				} //endif of plumUsed



				w = X[K+1];
				if   ((fcmp( w, 0.0) != 1) || (fcmp( w, 1.0) != -1)){  //w out of support, should be <0, 1>
					//Rprintf("Bacon: w out of support, should be <0, 1> %.2lf\n", w);
					return 0;
				}



				if (fcmp( X[K], 0.0) != 1){ //acc. rate alpha_{K} <= 0, out of support
					//Rprintf("Bacon: acc. rate alpha_{K} <= 0, out of support\n");
					return 0;
				}

				//Set the thetas, return if chronology exeeds general limit
				int rt = SetThetas(X);
				if( rt == 0 ){
						//Rprintf("Theta out of support\n");
						return 0;
				}

				//Check forst that all the alphas >=0
				for (int k=1; k<K; k++) {

					if (fcmp( x[k], 0.0) != 1) { //alpha_k <= 0
						//Rprintf("Bacon: alpha_k <= 0 (%.2lf <= 0.0)::::%d\n", X[k], k);
						return 0;
					}
				}





				//x has been created from X, all the rest is the same

				if (H > 0) {
				//if( false ){
					//Additional checks if there are hiatuses

					//we go backwards until we find the hiatus
					int l=0;

					for (int k=K-1; k>0; k--) {
						//printf("B: %d  %f  %f\n", k, x[k], (x[k]-w*x[k+1])/(1.0-w));
						if ((fcmp( c(k-1), h[l]) == -1) && (fcmp( h[l], c(k)) != 1)) { //forgets
//H Change
							if ((fcmp( x[k], 0.0) != 1) || (fcmp( hb[l], x[k]) != 1)){ //we require 0.0 < x[k] < hb[l]
								//Rprintf("we require 0.0 < x[k] < hb[l], %.2lf < %.2lf < %.2lf\n", 0.0, x[k], hb[l]);
								return 0;
							}
							l++; //jump to next hiatus, but max one hiatus in each section.
						} else if (fcmp( (x[k]-w*x[k+1])/(1.0-w), 0.0) != 1) { //e_k <= 0
							//Rprintf("e_k <= 0 %.2lf <= 0\n", (x[k]-w*x[k+1])/(1.0-w));
							return 0; // do not accept proposal where x <= 0 ???
						}
					}
				}



				if (plumUsed == 1) { //Check the chronology limit
					phi = x[K+2];

					double plumchronolim = (1.0/LA_CONST)*log( phi / (plumobj->GetAl()*LA_CONST) );



					//printf("PLUM %lf %lf\n",  G( dets->d(last210Pb), x), plumchronolim);
					//fflush(stdout);

					// G(dets)-theta0 >= plumchronolim
					if (fcmp( G( dets->d(last210Pb), x)-x[0], plumchronolim) != -1){
						//Rprintf("Plum 4: The chronology at the last 210Pb datum is beyond the plum chonology limit\n");
						//printf("PLUM %lf %lf\n",  G( dets->d(last210Pb), x)-theta[0], plumchronolim);
						return 0; //The chronology at the last 210Pb datum is beyond the plum chonology limit
					}

					//printf("BACON %lf %lf\n",  G( dets->d(last210Pb), x)-theta[0], plumchronolim);

				}

				return rt;
			 }


			//we assume the correct thetas are in theta
			//G is only called right after insupport
			 virtual double G( const double d, const double *x) {

				int i = (int) floor((d-c0)/Dc);
				return theta[i] + x[i+1]*(d-c(i));

			 }

			//This is a polymorphic version of the above age-depth model
			//where only the theta at the nearest node is returned.
			//we assume the correct thetas are in theta!!!!
			virtual double G(const double d) {

				int i = (int) floor((d-c0)/Dc);
				return theta[i] + (d-c(i))*(theta[i+1]-theta[i])/Dc;

			}

			virtual double G_Plum( const double d, const double *x, const double delta, const double AS, const double phi) {

	 			double th1 = G( d-delta, x)-theta[0], th = G( d, x)-theta[0];
	 			//AS = rho_i * PSi
	 			return AS + (phi/LA_CONST)*(exp(-LA_CONST*th1) - exp(-LA_CONST*th));

	 		}



	         ~BaconFix(){
				delete x;
				delete X0;
				delete Xp0;
				delete theta;
				delete dets;
			 }


			double Getc0() { return c(0); }
			double GetcK() { return c(K); }

			virtual void ShowDescrip() {
                Rprintf("BaconFixed: Bacon jumps model with fixed c's.\n");
                Rprintf("            K= %d, H= %d, dim= %d, Seed= %ld, Dc=%f, c(0)= %f, c(K)= %f\n",
					K, H, get_dim(), GetSeed(), Dc, c(0), c(K));
			}

			void PrintNumWarnings() {
				if (WarnBeyondLimits != 0) {
                    Rprintf("bacon: %d WarnBeyondLimits warnings:\n", WarnBeyondLimits);
                    Rprintf("bacon: WARNING: calibration attempted beyond MinYr= %f or MaxYr= %f\n", MinYr, MaxYr);
				}
			}

			 double *Getx0() { return X0; }
			 double *Getxp0() { return Xp0; }

			virtual double eval(double *X, int prime) {
				//x has been created from X in insupport function with the SetThetas function
				//The rest is the same!! X is ignored here.

        prime=0; //avoid warning JEV

				Uprior = 0.0;
				Uli = 0.0;

				//printf(" Init=%f", U);

				//assuming insupport is called right before eval


				if (useT) { //uses t model

					for (int j=0; j<(m-1); j++) {

						if (dets->Is210Pb(j) == 1)
							Uli += dets->Ut( j, G_Plum( dets->d(j), x, dets->Delta210Pb(j), dets->Rho210Pb(j)*GetPS(j, x), phi));
						else
							Uli += dets->Ut( j, G( dets->d(j), x)); //likelihood

					//printf("%d  %f  %f  %f\n", j, dets->d(j), G( dets->d(j), x), Uli);
					}
					//Uli += dets->Ut( m-1, G( dets->d(m-1), x));
				} else { //uses standard normal model

					for (int j=0; j<(m-1); j++) {

						if (dets->Is210Pb(j) == 1)
							Uli += dets->U( j, G_Plum( dets->d(j), x, dets->Delta210Pb(j), dets->Rho210Pb(j)*GetPS(j, x), phi));
						else
							Uli += dets->U( j, G( dets->d(j), x)); //likelihood

					//printf("%d  %f  %f  %f\n", j, dets->d(j), G( dets->d(j), x), Uli);
					}
					//Uli += dets->U( m-1, G( dets->d(m-1), x));
				}


				if ( plumUsed == 1) {
					//Uli += additional term in likelihood involving PS, eq (6)
					for (int j=0; j < plumobj->NS(); j++) {
						Uli += sqr(plumobj->yS(j) - GetPS(j, x))/(2*plumobj->s2(j));
					}


					Uprior += priorPhiU(x);

					Uprior += priorPSU(x);

				}








				Uprior += priorwU(w); //prior for w
				//printf(" priorw=%f", Uprior);


				//Set the prior for all accumulation rates
				Uprior += prioracU( 0, x[K]); //prior for alpha_K
				if (H == 0) {
					//printf("A: %d  %f  %d\n", 0, x[0], K);
					for (int k=1; k<K; k++) {
						Uprior += prioracU( 0, (x[k]-w*x[k+1])/(1.0-w)); //prior for e_k
						//printf("%f %f\n", (x[k]-w*x[k+1])/(1.0-w), U);
					}
				}
				else {

					//we go backwards until we find the hiatus
					int l=0;

					for (int k=K-1; k>0; k--) {
						if ((fcmp( c(k-1), h[l]) == -1) && (fcmp( h[l], c(k)) != 1)) { //forgets
							Uprior += priorHU( l, x[k]); //prior for the hiatus jump in hiatus l
							l++; //jump to next hiatus, but max one hiatus in each section.
						}
						else
							Uprior += prioracU( l, (x[k]-w*x[k+1])/(1.0-w)); //prior for e_k in section l
					}

				}

				//printf(" Uprior=%f, U=%f\n", Uprior, Uprior+Uli);
				/*if (!prime) {
					for (int i=0; i<n; i++)
                    //	printf("%f  ", x[i]);
                //	printf("%f\n", U);
				}*/

				//Rprintf("Memory %.8lf\n", x[K+1]);

				U = Uprior + Uli;


				return U;
			}
};


#endif
