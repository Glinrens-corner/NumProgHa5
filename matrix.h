
#ifndef MATRIX_H
#define MATRIX_H

#include "basic.h"

typedef struct
{
  int rows;
  int cols;
  int ld;

  double *a;
  bool issub;
} matrix;

typedef matrix *pmatrix;

typedef struct
{
  int rows;

  double *x;
} vector;

typedef vector *pvector;

typedef struct
{
  pvector a; //main diagonal
  pvector b; //secondary diagonal(s)
  int dim;
} symtridiag;

typedef symtridiag *psymtridiag;


/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pmatrix
new_matrix(int rows, int cols);

pmatrix
new_zero_matrix(int rows, int cols);

pmatrix
new_identity_matrix(int rows);

pmatrix
new_sub_matrix(pmatrix src, int rows, int roff,
               int cols, int coff);

pmatrix
clone_matrix(pmatrix a);

void copy_matrix(pmatrix src, pmatrix dest);

void del_matrix(pmatrix a);

pvector
new_vector(int rows);

pvector
new_zero_vector(int rows);

void copy_vector(pvector src, pvector dest);

void del_vector(pvector x);

psymtridiag
new_symtridiag(int rows);

void del_symtridiag(psymtridiag a);


/* ------------------------------------------------------------
 * Example matrix
 * ------------------------------------------------------------ */

pmatrix
new_diaghilbert_matrix(int rows);

pmatrix
new_hilbert_matrix(int rows);

psymtridiag
new_wave1d_tridiag(int rows);


/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void clear_matrix(pmatrix x);

void clear_vector(pvector x);

void print_matrix(pmatrix a);

void print_symtridiag(psymtridiag a);

void print_vector(pvector a);

void random_vector(pvector x);

double
normmax_vector(pvector x);

double
normmax_diff_vector(pvector x, pvector y);

double
norm2_diff_vector(const pvector x, const pvector y);

//Perform update y <- y + ax for symmetric tridiagonal matrix a
void eval_symtridiag(psymtridiag a, pmatrix x, pmatrix y);

#endif
