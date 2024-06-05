
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"
#include "miniblas.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pmatrix
new_matrix(int rows, int cols)
{

  pmatrix a;

  a = (pmatrix)malloc(sizeof(matrix));
  a->a = (double *)malloc(sizeof(double) * rows * cols);
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;
  a->issub = false;

  return a;
}

pmatrix
new_zero_matrix(int rows, int cols)
{

  pmatrix a;

  a = calloc(1, sizeof(matrix));
  a->a = calloc(rows * cols, sizeof(double));
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;
  a->issub = false;

  return a;
}

pmatrix
new_identity_matrix(int rows)
{

  pmatrix a;
  int i;

  a = new_zero_matrix(rows, rows);

  for (i = 0; i < rows; i++)
  {
    a->a[i + i * rows] = 1.0;
  }

  return a;
}

pmatrix
clone_matrix(pmatrix a)
{
  pmatrix b;
  int rows, cols, lda, ldb, i, j;

  rows = a->rows;
  cols = a->cols;
  lda = a->ld;

  b = new_zero_matrix(rows, cols);
  ldb = b->ld;

  for (j = 0; j < cols; j++)
  {
    for (i = 0; i < rows; i++)
    {
      b->a[i + j * ldb] = a->a[i + j * lda];
    }
  }

  return b;
}

void copy_matrix(pmatrix src, pmatrix dest)
{
  int rows = src->rows;
  int cols = src->cols;
  int ldsrc = src->ld;
  int lddst = dest->ld;

  int i, j;

  assert(dest->rows == rows);
  assert(dest->cols == cols);

  for (j = 0; j < cols; j++)
  {
    for (i = 0; i < rows; i++)
    {
      dest->a[i + j * lddst] = src->a[i + j * ldsrc];
    }
  }
}

pmatrix
new_sub_matrix(pmatrix src, int rows, int roff,
               int cols, int coff)
{
  int ld = src->ld;

  pmatrix a;

  assert(src != NULL);
  assert(roff + rows <= src->rows);
  assert(coff + cols <= src->cols);

  a = malloc(sizeof(matrix));

  a->a = src->a + roff + ld * coff;
  a->ld = src->ld;
  a->rows = rows;
  a->cols = cols;
  a->issub = true;

  return a;
}

void del_matrix(pmatrix a)
{
  if (!a->issub)
  {
    free(a->a);
  }
  free(a);
}

pvector
new_vector(int rows)
{

  pvector x;

  x = (pvector)malloc(sizeof(vector));
  x->x = (double *)malloc(sizeof(double) * rows);
  x->rows = rows;

  return x;
}

pvector
new_zero_vector(int rows)
{

  pvector x;
  int i;

  x = new_vector(rows);

  for (i = 0; i < rows; i++)
  {
    x->x[i] = 0.0;
  }

  return x;
}

void copy_vector(pvector src, pvector dest)
{
  int n = src->rows;

  int i;

  assert(dest->rows == n);

  for (i = 0; i < n; i++)
  {
    dest->x[i] = src->x[i];
  }
}

void del_vector(pvector x)
{
  free(x->x);
  free(x);
}

psymtridiag
new_symtridiag(int rows)
{
  psymtridiag a;

  assert(rows > 0);

  a = (psymtridiag)malloc(sizeof(symtridiag));
  a->dim = rows;
  a->a = new_vector(rows);
  a->b = new_vector(rows-1);

  return a;
}

void del_symtridiag(psymtridiag a)
{
  assert(a != NULL);

  del_vector(a->a);
  del_vector(a->b);
  free (a);
}


/* ------------------------------------------------------------
 * Example matrix
 * ------------------------------------------------------------ */

pmatrix
new_diaghilbert_matrix(int rows)
{

  pmatrix a;
  double *aa;
  double sum;
  int lda;
  int i, j;

  a = new_matrix(rows, rows);
  aa = a->a;
  lda = a->ld;

  for (j = 0; j < rows; j++)
  {
    sum = 1.0;
    for (i = 0; i < j; i++)
    {
      aa[i + j * lda] = 1.0 / (1.0 + i + j);
      sum += fabs(aa[i + j * lda]);
    }
    for (i = j + 1; i < rows; i++)
    {
      aa[i + j * lda] = 1.0 / (1.0 + i + j);
      sum += fabs(aa[i + j * lda]);
    }
    aa[j + j * lda] = sum;
  }

  return a;
}

pmatrix
new_hilbert_matrix(int rows)
{

  pmatrix a;
  double *aa;
  int lda;
  int i, j;

  a = new_matrix(rows, rows);
  aa = a->a;
  lda = a->ld;

  for (j = 0; j < rows; j++)
  {
    for (i = 0; i < rows; i++)
    {
      aa[i + j * lda] = 1.0 / (1.0 + i + j);
    }
  }

  return a;
}

psymtridiag
new_wave1d_tridiag(int rows)
{
  int i;
  real hinv;

  psymtridiag a;

  hinv = 1.0 * (rows + 1);
  hinv *= hinv;

  a = new_symtridiag(rows);

  for (i = 0; i < rows - 1; i++)
  {
    a->b->x[i] = -1.0 * hinv;
    a->a->x[i] = 2.0 * hinv;
  }
  a->a->x[rows - 1] = 2.0 * hinv;

  return a;
}

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void clear_matrix(pmatrix a)
{

  double *aa = a->a;
  int rows = a->rows;
  int cols = a->cols;
  int lda = a->ld;
  int i, j;

  for (j = 0; j < cols; j++)
    for (i = 0; i < rows; i++)
      aa[i + j * lda] = 0.0;
}

void clear_vector(pvector x)
{

  double *xx = x->x;
  int rows = x->rows;
  int i;

  for (i = 0; i < rows; i++)
    xx[i] = 0.0;
}

void print_matrix(pmatrix a)
{

  int i, j;
  int lda = a->ld;

  printf("Matrix (%d, %d)\n", a->rows, a->cols);

  for (i = 0; i < a->rows; i++)
  {
    printf("( %e ", a->a[i]);
    for (j = 1; j < a->cols; j++)
    {
      printf(", %e ", a->a[i + j * lda]);
    }
    printf(")\n");
  }
}

void print_symtridiag(psymtridiag a)
{
  int i, j;
  int rows = a->dim;

  printf("Symmetric tridiagonal matrix (%d, %d)\n", rows, rows);
  for (i = 0; i < rows; i++)
  {
    printf("( ");
    for (j = 0; j < i - 1; j++)
    {
      printf("         ");
    }

    if (i > 0)
      printf("%+8.2f ", a->b->x[i - 1]);

    printf("%+8.2f ", a->a->x[i]);

    if (i < rows - 1)
      printf("%+8.2f ", a->b->x[i]);

    for (j = i + 2; j < rows; j++)
    {
      printf("         ");
    }
    printf(")\n");
  }
}


void print_vector(pvector x)
{

  int i;
  int dim = x->rows;

  printf("Vector (%d)\n", dim);

  printf("( %e ", x->x[0]);
  for (i = 1; i < dim; i++)
  {
    printf(", %e ", x->x[i]);
  }
  printf(")\n");
}

void random_vector(pvector x)
{
  double *xx = x->x;
  int rows = x->rows;
  int i;

  for (i = 0; i < rows; i++)
  {
    xx[i] = 2.0 * rand() / RAND_MAX - 1.0;
  }
}

double
normmax_vector(const pvector x)
{

  double norm = 0.0;
  double *xx = x->x;
  int i;

  for (i = 0; i < x->rows; i++)
  {
    norm = (norm < fabs(xx[i]) ? fabs(xx[i]) : norm);
  }

  return norm;
}

double
normmax_diff_vector(const pvector x, const pvector y)
{

  double error = 0.0;
  int i;

  assert(x->rows == y->rows);

  for (i = 0; i < x->rows; i++)
  {
    error = ((fabs(x->x[i] - y->x[i])) > error ? (fabs(x->x[i] - y->x[i])) : error);
  }

  return error;
}

double
norm2_diff_vector(const pvector x, const pvector y)
{

  double *diff;
  double norm;
  int i;

  assert(x->rows == y->rows);
  diff = (double *)malloc(sizeof(double) * x->rows);

  for (i = 0; i < x->rows; i++)
  {
    diff[i] = x->x[i] - y->x[i];
  }

  norm = nrm2(x->rows, diff, 1);

  return norm;
}

void eval_symtridiag(psymtridiag a, pmatrix x, pmatrix y)
{
    assert(a != NULL);
    assert(x != NULL);

    int rows = a->dim;
    int m = x->cols;
    real *xx = x->a;
    int xld = x->ld;
    real *yx = y->a;
    int yld = y->ld;
    real *ax = a->a->x;
    real *bx = a->b->x;

    assert(x->rows == rows);
    assert(y->rows == rows);
    assert(y->cols == m);

    int k, j;

    for (j = 0; j < m; j++)
    {
        for (k = 0; k < 1; k++)
        {
            yx[k + j * yld] += ax[k] * xx[k + j * xld] + bx[k] * xx[(k + 1) + j * xld];
        }
        for (k = 1; k < rows - 1; k++)
        {
            yx[k + j * yld] += bx[k - 1] * xx[(k - 1) + j * xld] + ax[k] * xx[k + j * xld] + bx[k] * xx[(k + 1) + j * xld];
        }
        for (k = rows - 1; k < rows; k++)
        {
            yx[k + j * yld] += bx[k - 1] * xx[(k - 1) + j * xld] + ax[k] * xx[k + j * xld];
        }
    }
}

