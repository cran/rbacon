/*
   This code is derived from the C++ version of t-Walk "cpptwalk-beta-1.0"
   Feb 2008 version, kindly supplied by J. Andres Christen.

   Tony Begg <tony.begg@dataventures.com>
   * 
   * Andres Christen 08DEC2016
*/





//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Rcpp.h>

#define New(objects, have, TYPE)                                            \
   (objects) = (TYPE *)malloc((have) * sizeof(*(objects)));                 \
   if ((objects) == NULL)                                                   \
   {                                                                        \
     /*fprintf(stderr, "New(): Out of Heap Memory (%d)\n", have); JEV WARNING CHECK   */          \
      Rcpp::stop("memory problem, stopping");                                                             \
   }

/*
   fcmp
   Copyright (c) 1998-2000 Theodore C. Belding
   University of Michigan Center for the Study of Complex Systems
   <mailto:Ted.Belding@umich.edu>
   <http://www-personal.umich.edu/~streak/>
  
   This file is part of the fcmp distribution. fcmp is free software;
   you can redistribute and modify it under the terms of the GNU Library
   General Public License (LGPL), version 2 or later.
*/
int fcmp(double x1, double x2, double epsilon = 0.00000000001)
{
   int exponent;
   double delta;
   double difference;
   
   /*
      Get exponent(max(fabs(x1), fabs(x2))) and store it in exponent.
   */
   frexp(fabs(x1) > fabs(x2) ? x1 : x2, &exponent);
   /*
      Do the comparison.
      delta = epsilon * pow(2, exponent)
   */
   delta = ldexp(epsilon, exponent); 
   difference = x1 - x2;
   if (difference > delta)
   {
      return 1; /* x1 > x2 */
   }
   else if (difference < -delta) 
   {
      return -1;  /* x1 < x2 */
   }
   else /* -delta <= difference <= delta */
   {
      return 0;  /* x1 == x2 */
   } 
}

/*
   Random is based on rng/taus.c from the GNU Scientific Library.
   Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
   
   NormalDev is based on randist/gauss.c from the GNU Scientific Library.
   Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006 James Theiler, Brian Gough
   Copyright (C) 2006 Charles Karney

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
struct Random
{
   unsigned long int s, s1, s2, s3;
};

typedef struct Random Random;

static unsigned long Random32(Random *rng)
{
#define MASK 0xffffffffUL
#define TAUSWORTHE(s,a,b,c,d) (((s & c) << d) & MASK) ^ ((((s << a) & MASK)^s) >> b)

   rng->s1 = TAUSWORTHE (rng->s1, 13, 19, 4294967294UL, 12);
   rng->s2 = TAUSWORTHE (rng->s2, 2, 25, 4294967288UL, 4);
   rng->s3 = TAUSWORTHE (rng->s3, 3, 11, 4294967280UL, 17);
   return (rng->s1 ^ rng->s2 ^ rng->s3);
}




/*
   Here is our random number generator object.  It is used for
   all our deviate needs.
*/
static Random RNG;


void RandomSeed(Random *rng, unsigned long int s)
{
   if (s == 0) s = (unsigned long int) time(NULL);      /* default seed is use a seed taken from the calendar */
   rng->s = s;
#define LCG(n) ((69069 * n) & 0xffffffffUL)
   rng->s1 = LCG(s);
   if (rng->s1 < 2) rng->s1 += 2UL;
   rng->s2 = LCG(rng->s1);
   if (rng->s2 < 8) rng->s2 += 8UL;
   rng->s3 = LCG(rng->s2);
   if (rng->s3 < 16) rng->s3 += 16UL;
   /* "warm it up" */
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   Random32(rng);
   return;
}

void Seed(unsigned long int s) {
	RandomSeed( &RNG, s);
}

unsigned long int GetSeed() {
	return RNG.s;
}
	



/*
   Uniform real random number between 0 and 1.
*/

inline double Random32_Un01(Random *rng)
{
   return (double) Random32(rng) / 4294967296.0;
}



double Un01()
{
   return Random32_Un01(&RNG);
}


/*
   Uniform real random number between values a and b where b > a.
*/
inline double Unab( double a, double b)
{
   return (b - a) * Un01() + a;
}

double NormalDev(Random *rng, double mean, double sigma)
{
   double u, v, x, y, Q;
   const double s = 0.449871;  /* Constants from Leva */
   const double t = -0.386595;
   const double a = 0.19600;
   const double b = 0.25472;
   const double r1 = 0.27597;
   const double r2 = 0.27846;

   do /* This loop is executed 1.369 times on average  */
   {
      /*
         Generate a point P = (u, v) uniform in a rectangle enclosing
         the K+M region v^2 <= - 4 u^2 log(u).
         u in (0, 1] to avoid singularity at u = 0.
      */
      u = 1.0 - Random32_Un01(rng);
      /*
         v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
         is rejected in the last part of the while clause.  The
         resulting normal deviate is strictly symmetric about 0
         (provided that v is symmetric once v = -0.5 is excluded).
      */
      v = Random32_Un01(rng) - 0.5;
      /*
         Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
         much (for efficiency).
      */
      v *= 1.7156;
      /*
         Compute Leva's quadratic form Q.
      */
      x = u - s;
      y = fabs(v) - t;
      Q = x * x + y * (a * y - b * x);
      /*
         Accept P if Q < r1 (Leva)
         Reject P if Q > r2 (Leva)
         Accept if v^2 <= -4 u^2 log(u) (K+M)
         This final test is executed 0.012 times on average.
      */
   }
   while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));
   return mean + (sigma * (v / u)); /* Return slope plus mean */
}

/* TAKEN FROM: GNU GSL  randist/gamma.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * *************
 * New version based on Marsaglia and Tsang, "A Simple Method for
 * generating gamma variables", ACM Transactions on Mathematical
 * Software, Vol 26, No 3 (2000), p363-372.
 *
 * Implemented by J.D.Lamb@btinternet.com, minor modifications for GSL
 * by Brian Gough
 */


double GammaDev(Random *rng, const double a, const double b)
{
  /* assume a > 0 */

  if (a < 1)
    {
      double u = Random32_Un01(rng); /*gsl_rng_uniform_pos (r);*/
      return GammaDev( rng, 1.0 + a, b) * pow(u, 1.0 / a);
    }

  {
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    while (1)
      {
        do
          {
            x = NormalDev(rng, 0.0, 1.0); /*gsl_ran_gaussian_ziggurat (r, 1.0);*/
            v = 1.0 + c * x;
          }
        while (v <= 0);

        v = v * v * v;
        u = Random32_Un01(rng); /*gsl_rng_uniform_pos (r);*/

        if (u < 1 - 0.0331 * x * x * x * x) 
          break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
    
    return b * d * v;
  }
}

/*Simple implementation of Beta random gen. using two gammas*/
double BetaDev( Random *rng, const double a, const double b) {

	double g1=GammaDev( rng, a, 1.0);
	double g2=GammaDev( rng, b, 1.0);

	return g1/(g1+g2);
}



double NorSim( double m, double sigma) { /*Normal*/
	return NormalDev( &RNG, m, sigma);
}

double GammaSim( double a, double b) { /*Gamma: p(x) dx = K x^{a-1} e^{-x/b} dx*/
	return GammaDev( &RNG, a, b);
}

double BetaSim( double a, double b) { /*Beta: p(x) dx = K x^{a-1} (1-x)^{b-1} dx */
	return BetaDev( &RNG, a, b);
}











#define min(x,y) ((x) < (y) ? (x) : (y))





/*
   Some useful functions for vectors of doubles.
*/
/*
   Allocate storage for a vector of length n.
*/  
double *Vector(int n)
{
   double *v;

   New(v, n, double);
   return v;
}

/*
   Free storage for a vector.
*/  
void FreeVector(double *v)
{
   free((void *)v);
}

/*
   Copy vector src to vector dest.
*/
void VectorCopy(double *src, double *dest, int n)
{
   int i;

   for (i = 0; i < n; i++) dest[i] = src[i];
}

/*
   Subtract vector v2 from v1, elementwise, and store result in diff.
*/
void VectorDiff(double *v1, double *v2, int n, double *diff)
{
   int i;

   for (i = 0; i < n; i++) diff[i] = v1[i] - v2[i];
}

/*
   Find the maximum absolute value of a vector (masked by phi) and
   return its index.
*/
void VectorIndexMax(double *v, int n, int *ix, int *phi)
{
   int i, max_i;

   max_i = 0;
   for (i = 0; i < n; i++)
   {
      max_i = ((fcmp((double)phi[max_i] * fabs(v[max_i]),
                     (double)phi[i] * fabs(v[i])) == -1) ? i : max_i);
   }
   *ix = max_i;
}

/*
   Compare vector v1 with v2, elementwise, and return 1 if all
   elements agree (are equal), else return 0.
*/
int VectorCmp(double *v1, double *v2, int n)
{
   int i = 0;

   while ((fcmp(v1[i], v2[i]) == 0) && (i < n)) i++;
   if (i == n) return 1;
   return 0;
}

/*
   Print out a vector on the file pointer fp in a particular format.
   Note that large n will lead to very long lines.
*/ 
void VectorPrint(FILE *fp, double *v, int n)
{
   int i;

   fprintf(fp, "\n");
   for (i = 0; i < n; i++)
   {
      fprintf(fp, "\t%13.6g", v[i]);
   }
}


