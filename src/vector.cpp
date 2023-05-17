/*
 *  vector.cpp
 *
 *
 *  Created by jac on 29/04/2009.
 *
 */

#include "twalk.h"


//---------------------------------------------------------------------------
// Creates a vector
double *vector(int ncols) {
        double *cabeza;
        cabeza = (double *)malloc(ncols*sizeof(double));

        return cabeza;
}


//---------------------------------------------------------------------------
// Destroys a vector
void free_vector(double *cabeza) { /*free((void *) cabeza);*/ }
//---------------------------------------------------------------------------


// void fver_vector(FILE *fptr, double *u, int m) { // this is how this function was written in bacon's vector.cpp
//         fprintf(fptr, "\n");
//         for(int i=0; i<(m); i++) // rplum says m-1!
//             fprintf(fptr, "\t%13.6g", u[i]);
// 		//      fprintf(fptr, "%13.6g", u[m-1]); // this from rplum. Not sure if it makes any difference
// }

void fver_vector(FILE *fptr, double *u, int m) { // these lines are from plum's vector.cpp
         fprintf(fptr, "\n");
         for(int i=0; i<m-1; i++) // plum says m-1, bacon (m). Plum runs go bad if set to (m)
         //for(int i=0; i<m; i++)             
         	fprintf(fptr, "%13.6g\t", u[i]);
         fprintf(fptr, "%13.6g", u[m-1]); // plum has this line, but bacon not
}


/*Returns 1 if v = u, component to component*/
int vector_cmp(double *v, double *u, int n) {
        int i = 0;
        while ((fcmp( v[i], u[i]) == 0) && (i<n))
                i++;

        if (i==n)
                return 1;
        else
                return 0;
}

//copia vector v1 en v2
void cp_vector(double *v1, double *v2, int n) {
        for(int k=0; k<n ; k++)
               v2[k] = v1[k];

}

//---------------------------------------------------------------------------
void resta_vector(double *v1, double *v2, int n, double *res) {
        for(int k=0; k<n; k++) {
                res[k] = v1[k]-v2[k];
        }
}

/*encuentra el  maximo  ai en valor abs de una vector y devuelve su indice*/
void indice_max_vector(double *v, int n, int &indice, int *phi) {

        indice =0;

        for(int i=0; i<n; i++)   ///       A                  <         B
                indice = ( (fcmp( phi[indice]*fabs(v[indice]), phi[i]*fabs(v[i])) == -1) ? i : indice);
//	return indice;

}


int sum(int *v, int n) {
  int summ = 0;
  for(int i=0; i<n; i++)
    summ+=v[i];

  return summ;
}
