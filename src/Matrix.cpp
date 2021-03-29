/*
 *  Matrix.c
 *  
 *
 *  Created by and copyright of J. Andres Christen on 17/04/2007.
 *
 */

//#include <stdio.h>
#include <math.h>

#include "Matrix.h"


gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2) {
	gsl_matrix *ot = new gsl_matrix;
	ot->size1 = n1;
	ot->size2 = n2;
	ot->tda = n2;
	ot->data = new double[n1*n2];
	ot->owner = 1;
	return ot;
}

void gsl_matrix_free (gsl_matrix * m) {
	delete [] m->data;
	delete m; 
}


void gsl_matrix_set_all (gsl_matrix * m, double x) {
    for (unsigned int i=0; i<m->size1; i++)
        for (unsigned int j=0; j<m->size2; j++)
			m->data[ (m->tda)*i + j] = x;
}

int gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src) {
	if ((dest->size1 != src->size1) || (dest->size2 != src->size2)) {
            REprintf("ERROR: copy only allowed for same size matrices.");
            Rcpp::stop("ERROR: copy only allowed for same size matrices.");
        //	exit(1);
		}
    for (unsigned int i=0; i<dest->size1; i++)
        { for (unsigned int j=0; j<dest->size2; j++)
            { dest->data[ (dest->tda)*i + j] = src->data[ (src->tda)*i + j]; }}
    return 1;
}
