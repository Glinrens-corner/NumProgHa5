
#ifndef FACTORIZATIONS_H
#define FACTORIZATIONS_H

#include <assert.h>

/* Additional libraries */
#include "basic.h"    /* basic types and time measurement */
#include "miniblas.h" /* our simple BLAS version */
#include "matrix.h"   /* matrix functions */

/* ---------------------------
 *  RDR^T factorization
 * --------------------------- */

/**
 * @brief Compute RDR^T decompositon of a given Matrix A.
 *        After completion the upper triangulatr part of the matrix A
 *        is being overwritten by the diagonal matrix D and the unit upper
 *        triangular matrix R.
 *
 * @param Aa Pointer to the first element of A
 * @param ldA Leading dimension of A
 * @param rows Number of rows of A (also number of columns)
 */
HEADER_PREFIX
void decomp_rdrt(real *Aa, int ldA, int rows);

/**
 * @brief Evaluate an upper triangular matrix R with a vector x. The result
 *        y = \alpha R x overwrites the input vector x.
 *
 * @param unit If true a unit upper triangular matrix is assumed, otherwise
 *             diagonal entries will be used.
 * @param alpha Scaling factor
 * @param R Upper triangular matrix R
 * @param x Initially the input vector x, after completion the result vector y
 */
HEADER_PREFIX
void eval_upper(bool unit, real alpha, pmatrix R, pvector x);

/**
 * @brief Evaluate an upper triangular matrix R with a vector x. The result
 *        y = \alpha R^T x overwrites the input vector x.
 *
 * @param unit If true a unit upper triangular matrix is assumed, otherwise
 *             diagonal entries will be used.
 * @param alpha Scaling factor
 * @param R Upper triangular matrix R
 * @param x Initially the input vector x, after completion the result vector y
 */
HEADER_PREFIX
void eval_uppertrans(bool unit, real alpha, pmatrix R, pvector x);

/**
 * @brief Evaluate diagonal matrix D with a vector x. The result
 *        y = \alpha D x overwrites the input vector x.
 *
 * @param alpha Scaling factor
 * @param R Diagonal matrix D
 * @param x Initially the input vector x, after completion the result vector y
 */
HEADER_PREFIX
void eval_diag(real alpha, pmatrix R, pvector x);

/**
 * @brief Evaluate the decomposition (R D R^T) with a vector x. The result
 *        y = \alpha (R D R^T) x overwrites the input vector x.
 *
 * @param alpha Scaling factor
 * @param R Upper triangular matrix R and diagonal matrix D
 * @param x Initially the input vector x, after completion the result vector y
 */
HEADER_PREFIX
void eval_rdrt(real alpha, pmatrix R, pvector x);

/**
 * @brief Solve the equation R x = b with R being an upper triangular matrix
 *
 * @param unit If 'true' the value 1 is assumed to be on every diagonal entry of
 *             the matrix. Otherwise the diagonal elements from R will be used
 *             instead.
 * @param R Upper triangular matrix
 * @param x On entry the right-hand-side b is stored in this vector.
 *          After completion of the routine the solution x to the above system
 *          is being stored.
 */
HEADER_PREFIX
void solve_upper(bool unit, pmatrix R, pvector x);

/**
 * @brief Solve the equation D x = b with D being a diagonal matrix
 *
 * @param D Diagonal matrix
 * @param x On entry the right-hand-side b is stored in this vector.
 *          After completion of the routine the solution x to the above system
 *          is being stored.
 */
HEADER_PREFIX
void solve_diag(pmatrix R, pvector x);

/**
 * @brief Solve the equation R^T x = b with R being an upper triangular matrix
 *
 * @param unit If 'true' the value 1 is assumed to be on every diagonal entry of
 *             the matrix. Otherwise the diagonal elements from R will be used
 *             instead.
 * @param R Upper triangular matrix
 * @param x On entry the right-hand-side b is stored in this vector.
 *          After completion of the routine the solution x to the above system
 *          is being stored.
 */
HEADER_PREFIX
void solve_uppertrans(bool unit, pmatrix R, pvector x);

/**
 * @brief Solve the equation R D R^T x = b with R being an upper triangular matrix,
 *        D being a diagonal matrix.
 *
 * @param R Upper triangular matrix 'R' is stored in the upper triagonal part,
 *          The matrix 'D' is stored on its diagonal.
 * @param x On entry the right-hand-side b is stored in this vector.
 *          After completion of the routine the solution x to the above system
 *          is being stored.
 */
HEADER_PREFIX
void solve_rdrt(pmatrix R, pvector x);

/**
 * @brief Solve linear system with unit upper triangular matrix R of the form:
 *        R X^T = B^T ,
 *        with R having dimensions rows x rows,
 *          X, B having dimensions cols x rows.
 *        Initially the right-hand-side Matrix B is being stored in X,
 *        after completion the solution to the system is stored in X.
 *
 * @param rows Number of rows of R and columns of X
 * @param cols Number of rows of X
 * @param Ra Pointer to the first element of R
 * @param ldR Leading dimension of R
 * @param Xa Pointer to the first element of X and B respectively
 * @param ldX Leading dimension of X and B respectively
 */
HEADER_PREFIX
void block_rsolve(int rows, int cols, real *Ra, int ldR, real *Xa, int ldX);

/**
 * @brief Solve linear system with diagonal matrix D of the form:
 *        D X = B,
 *        with R having dimensions rows x rows,
 *          X, B having dimensions cols x rows.
 *        Initially the right-hand-side Matrix B is being stored in X,
 *        after completion the solution to the system is stored in X.
 *
 * @param rows Number of rows of D and columns of X
 * @param cols Number of rows of X
 * @param Ra Pointer to the first element of D
 * @param ldR Leading dimension of D
 * @param Xa Pointer to the first element of X and B respectively
 * @param ldX Leading dimension of X and B respectively
 */
HEADER_PREFIX
void block_diagsolve(int rows, int cols, real *Ra, int ldR, real *Xa, int ldX);

/**
 * @brief Compute matrix-matrix products of the form
 *           X <-- X + \alpha A D B^T,
 *        where X is of dimension rows x cols,
 *              A is of dimension rows x k,
 *              B is of dimension cols x k and
 *              D is of dimension k x k and is a diagonal matrix.
 *
 * @param rows Number of rows of matrix A and X
 * @param cols Number of columns of matrix X and rows of B
 * @param k Intermediate dimension of matrices A,D,B
 * @param alpha Scaling factor
 * @param Aa Pointer to first element of matrix A
 * @param ldA Leading dimension of A
 * @param Da Pointer to first element of matrix D
 * @param ldD Leading dimension of D
 * @param Ba Pointer to first element of matrix B
 * @param ldB Leading dimension of B
 * @param Xa Pointer to first element of matrix X
 * @param ldX Leading dimension of X
 * @return HEADER_PREFIX
 */
HEADER_PREFIX
void addmul_rdrt(int rows, int cols, int k, real alpha, real *Aa, int ldA, real *Da, int ldD, real *Ba, int ldB, real *Xa, int ldX);

/* ------------------------------------------------------------
 * RDR^T decomposition for tridiagonal matrices
 * ------------------------------------------------------------ */

/**
 * @brief Compute RDR^T decompositon of a given symmetric  tridiagonal matrix A.
 *        After completion the upper triangulatr part of the matrix A
 *        is being overwritten by the diagonal matrix D and the unit upper
 *        triangular matrix R. (compare decomp_rdrt)
 *
 * @param a Coefficients and dimensions of symmetric tridiagonal matrix A
 */
HEADER_PREFIX
void decomp_rdrt_symtridiag(psymtridiag a);

/**
 * @brief Evaluate the decomposition (R D R^T) of a symmetric tridiagonal matrix with a matrix X.
 *        The result Y = (R D R^T) X overwrites the input matrix X.
 *
 * @param a (R D R^T) decomposed symmetric tridiagonal matrix.
 * @param x Initially the input matrix X,  after completion the result matrix Y
 */
void eval_rdrt_symtridiag(psymtridiag a, pmatrix x);

/**
 * @brief Solve the equation R D R^T X = B with R being an upper triangular matrix,
 *        D being a diagonal matrix and X and B are general matrices.
 *
 * @param a (R D R^T) decomposed symmetric tridiagonal matrix.
 * @param x Initially the input matrix B,  after completion the result matrix X
 */
void solve_rdrt_symtridiag(psymtridiag a, pmatrix x);

/* ------------------------------------------------------------
 * QR factorization
 * ------------------------------------------------------------ */

// Compute QR decomposition of a matrix A
// tau stores scaling factors
// work is auxiliary storage of dimension lwork

/**
 * @brief Compute a QR-decomposition of a given matrix A. Scaling coefficients
 *        are being stored within the vector 'tau'.
 *        work is intermediate storage to apply the householder reflections.
 *        its size is determined by 'lsize'
 * 
 * @param A The matrix which is being factorized: after completion the upper
 *          triangular matrxi will be stored within the upper triangular part
 *          of A. The lower triangular part will be used to store the Householder
 *          vectors except their first component, which is being assumed to be 1.
 * @param tau Vector containing als the scaling coeffiecents for the
 *            Householder reflections
 * @param lwork Size of 'work'
 * @param work Intermediate storage that is needed to apply the Householder 
 *             reflections to remainder of the matrix.
 */
void qrdecomp(pmatrix A, real *tau, int lwork, real *work);

/**
 * @brief Solves Q y = x or Q^T y = x, respectively and overwrite input x 
 *        with the result y. The matrix Q is stored in terms of Householder
 *        reflections within the matrix A. 
 * 
 * @param trans if true, then solve for Q^T y = x, otherwise solve for Q y = x.
 * @param A Matrix containing the Householder reflectiosns of Q
 * @param tau Vector containing als the scaling coeffiecents for the
 *            Householder reflections
 * @param x Initially the right-hand-side of the linear system 'x', after
 *          completion the solution vector 'y'.
 */
void qsolve(bool trans, pmatrix A, const real *tau, pvector x);

/**
 * @brief Compute orthogonal matrix Q from QR decomposition stored in A.
 * 
 * @param A Matrix containing the Householder reflections that are needed to 
 *          construct Q
 * @param tau Vector containing als the scaling coeffiecents for the
 *            Householder reflections
 * @param Q Matrix the will be overwritten with the Matrix Q of the
 *          QR-factorization stored in A.
 * @param lwork Size of 'work'
 * @param work Intermediate storage that is needed to apply the Householder 
 *             reflections.
 */
void qexpand(pmatrix A, const real *tau, pmatrix Q, int lwork, real *work);

#endif /* FACTORIZATIONS_H */
