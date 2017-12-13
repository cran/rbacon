
/*
   This code is derived from the C++ version of t-Walk "cpptwalk-beta-1.0"
   Feb 2008 version, kindly supplied by J. Andres Christen.

   Tony Begg <tony.begg@dataventures.com>
   * 
   * Andres Christen 08DEC2016
*/


#include <time.h>
#include <math.h>
//#include <stdlib.h>

#ifndef RANFUN_H
#define RANFUN_H



inline double sqr( double x) { return (x*x);}

int fcmp (double x1, double x2, double epsilon = 0.00000000001); 

void Seed(unsigned long int s);

unsigned long int GetSeed();

double Un01();  /*Un01() */

double Unab(double a, double b);  /*U(a,b]*/

double NorSim(double m, double sigma); /*Normal*/

double GammaSim(double a, double b); /*Gamma: p(x) dx = K x^{a-1} e^{-x/b} dx*/

double BetaSim(double a, double b); /*Beta: p(x) dx = K x^{a-1} (1-x)^{b-1} dx */

#endif







 
