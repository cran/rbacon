/* Original simple matrix class taken from http://www.algarcia.org/nummeth/Cpp/Matrix.h */

#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>  // Defines the assert function.

#include <string.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <Rcpp.h>

#define BUFF_SIZE 4096


typedef struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  //gsl_block * block;
  int owner;
} gsl_matrix;


gsl_matrix * gsl_matrix_alloc (size_t n1, size_t n2);
void gsl_matrix_free (gsl_matrix * m);
void gsl_matrix_set_all (gsl_matrix * m, double x);
int gsl_matrix_memcpy (gsl_matrix * dest, const gsl_matrix * src);

//*********************************************************************
class Matrix {
protected:

	// Matrix data.
	gsl_matrix *ma;
	char *header;

public:

// Default Constructor. Creates a 1 by 1 matrix; sets value to zero.
Matrix () {
  ma = gsl_matrix_alloc( 1, 1);// Allocate memory
  set(0.0);                // Set value of data_[0] to 0.0

  header = NULL;
}

// Regular Constructor. Creates an nR by nC matrix; sets values to zero.
// If number of columns is not specified, it is set to 1.
Matrix(int nR, int nC = 1, char *head=NULL) {
  assert(nR > 0 && nC > 0);    // Check that nC and nR both > 0.
  ma = gsl_matrix_alloc( nR, nC);// Allocate memory
  set(0.0);                    // Set values of data_[] to 0.0

  if (head != NULL)
	header = strdup(head);
  else
	header = NULL;
}


//copy constructor
Matrix(Matrix &mat)  {
Rprintf("Matrix::Matrix 2 ...\n");
  ma = gsl_matrix_alloc( mat.nRow(), mat.nCol());// Allocate memory
  copy(mat);

  if (mat.Header() != NULL)
	header = strdup(mat.Header());
  else
	header = NULL;

}


// Destructor. Called when a Matrix object goes out of scope or deleted.
~Matrix() {
	//printf("Matrix::~Matrix\n");
	if (ma != NULL)
		  gsl_matrix_free(ma);   // Release allocated memory
	if (header != NULL) {
		free(header);
	}
}



// Assignment operator function.
// Overloads the equal sign operator to work with
// Matrix objects.
Matrix& operator=(const Matrix& mat) {
  if( this == &mat ) return *this;  // If two sides equal, do nothing.
  this->copy(mat);                  // Copy right hand side to l.h.s.
  return *this;
}



// Set function. Sets all elements of a matrix to a given value.
void set(double value) {
   gsl_matrix_set_all( ma, value);
}


//info
// Simple "get" functions. Return number of rows or columns.
int nRow() const { return ma->size1; }
int nCol() const { return ma->size2; }


// Parenthesis operator function.
// Allows access to values of Matrix via (i,j) pair.
// Example: a(1,1) = 2*b(2,3);
// If column is unspecified, take as 1.
double& operator() (int i, int j = 0) {
  assert(i >= 0 && (size_t)i < ma->size1);          // Bounds checking for rows
  assert(j >= 0 && (size_t)j < ma->size2);          // Bounds checking for columns
  return ma->data[ (ma->tda)*i + j];  // Access appropriate value
}

// Parenthesis operator function (const version).
const double& operator() (int i, int j = 0) const{
  assert(i >= 0 && (size_t)i < ma->size1);          // Bounds checking for rows
  assert(j >= 0 && (size_t)j < ma->size2);          // Bounds checking for columns
  return ma->data[(ma->tda)*i + j];  // Access appropriate value
}

// Allows access to values of Matrix via ele(i,j) pair.
// If column is unspecified, take as 1.
double ele(int i, int j = 0) {
  assert(i >= 0 && (size_t)i < ma->size1);          // Bounds checking for rows
  assert(j >= 0 && (size_t)j < ma->size2);          // Bounds checking for columns
  return ma->data[ (ma->tda)*i + j];  // Access appropriate value
}



int filescan(char *fnam, int file_header=0) {
	FILE *F;

	if ((F = fopen( fnam, "r")) == NULL)
	{
        Rprintf( "File %s not found\n", fnam);

		return 0;
	}
	else
	{
		if (file_header == 1) {
			header = (char *) malloc((size_t) BUFF_SIZE);

			header = fgets( header, (size_t) BUFF_SIZE, F);
			header[strlen(header)-1] = '\0'; //remove the \n
		}

		//gsl_matrix_fscanf( F, ma);
		int k=0;
		double tmp;
		//printf("\n");
		while (fscanf( F, " %lf", &tmp) == 1) {
            if ((size_t) k >= ma->size1*ma->size2) {
                REprintf("ERROR: Reading matrix/table from file larger than previously opened.\n");
				return 0;
			}
			ma->data[k] = tmp;
			k++;
			//printf("%13.6g ", tmp);
			//if ((k % ma->size2) == 0)
			//	printf("\n");
		};
		//printf("\nk=%d, size1*size2=%d\n", k, ma->size1*ma->size2);
        if ((size_t) k < ma->size1*ma->size2) {
                Rprintf("WARNING: Read matrix/table from file smaller than previously opened.\n");
		}

		fclose(F);
		return 1;
	}
}

/*********************** Views *********************/
//Get the gsl matrix
gsl_matrix *Ma() { return ma; }


// Copy function.
// Copies values from one Matrix object to another.
void copy(const Matrix& mat) {
    gsl_matrix_memcpy( ma, mat.ma);
}


const char *Header() { return header; }


}; // Class Matrix








class SubMatrix : public Matrix {

private:
	Matrix *Parent;

public:

	SubMatrix() {
	//delay the opening of this submatrix
		ma = NULL;
		Parent = NULL;
		header = NULL;
	}

	~SubMatrix() {
		ma = NULL; //we avoid the base class call to free with this
		if (header != NULL) {
			free(header);
			header=NULL;
		}

	}


	const char *SetHeader(char *head) {
		if (header != NULL)
			free(header);
		header = NULL;
		if (head != NULL)
			header = strdup(head);
		return header;
	}

/***************** WARNING **************************/
/** ONLY 'SUBMATRICES' USED ARE EQUAL TO THE PARENT MATRIX **/
/** ONLY THOSE ARE USED HERE, JUST TO DELAY THE DEFINITION OF THE MATRIX **/
	void Set( Matrix *mat, size_t n1, size_t n2) {
		//view = gsl_matrix_submatrix ( mat->Ma(), (size_t) 0, (size_t) 0, n1, n2);
        if ((n1 != (size_t)mat->nRow()) || (n2 != (size_t)mat->nCol())) {
            REprintf("ERROR: resizing of submatrix not allowed.");
            Rcpp::stop("ERROR: resizing of submatrix not allowed.");
            //exit(1);
		}
		ma = mat->Ma(); //&view.matrix;
		Parent = mat;
	}



}; // class SubMatrix




#endif
