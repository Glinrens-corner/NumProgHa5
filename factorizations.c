
#include "factorizations.h"

/* ---------------------------
 *  RDR^T factorization
 * --------------------------- */

void decomp_rdrt(real *Aa, int ldA, int rows)
{
  int k;

  for (k = rows; k-- > 0;)
  {
    scal(k, 1.0 / Aa[k + k * ldA], Aa + k * ldA, 1);
    syr(false, k, -Aa[k + k * ldA], Aa + k * ldA, 1, Aa, ldA);
  }
}

void eval_upper(bool unit, real alpha, pmatrix R, pvector x)
{
  int i, j;
  real sum;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  if (unit)
  {
    for (i = 0; i < n; i++)
    {
      sum = Xx[i] + dot(n - i - 1, Ra + i + (i + 1) * ldR, ldR, Xx + (i + 1), 1);
      Xx[i] = alpha * sum;
    }
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      sum = 0.0;
      for (j = i; j < n; j++)
      {
        sum += Ra[i + j * ldR] * Xx[j];
      }
      Xx[i] = alpha * sum;
    }
  }
}

void eval_uppertrans(bool unit, real alpha, pmatrix R, pvector x)
{
  int i;
  real sum;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  if (unit)
  {
    for (i = n; i-- > 0;)
    {
      sum = Xx[i] + dot(i, Ra + i * ldR, 1, Xx, 1);
      Xx[i] = alpha * sum;
    }
  }
  else
  {
    for (i = n; i-- > 0;)
    {
      sum = dot(i + 1, Ra + i * ldR, 1, Xx, 1);
      Xx[i] = alpha * sum;
    }
  }
}

void eval_diag(real alpha, pmatrix R, pvector x)
{
  int i;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  assert(n == R->cols);
  assert(n = x->rows);

  for (i = 0; i < n; i++)
  {
    Xx[i] *= alpha * Ra[i + i * ldR];
  }
}

void eval_rdrt(real alpha, pmatrix R, pvector x)
{
  eval_uppertrans(true, 1.0, R, x);
  eval_diag(alpha, R, x);
  eval_upper(true, 1.0, R, x);
}

void solve_upper(bool unit, pmatrix R, pvector x)
{
  int i;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  if (unit)
  {
    for (i = n; i-- > 0;)
    {
      axpy(i, -Xx[i], Ra + i * ldR, 1, Xx, 1);
    }
  }
  else
  {
    for (i = n; i-- > 0;)
    {
      Xx[i] /= Ra[i + i * ldR];
      axpy(i, -Xx[i], Ra + i * ldR, 1, Xx, 1);
    }
  }
}

void solve_diag(pmatrix R, pvector x)
{
  int i;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  assert(n == R->cols);
  assert(n = x->rows);

  for (i = 0; i < n; i++)
  {
    Xx[i] /= Ra[i + i * ldR];
  }
}

void solve_uppertrans(bool unit, pmatrix R, pvector x)
{
  int i;

  int n = R->rows;
  real *Ra = R->a;
  int ldR = R->ld;
  real *Xx = x->x;

  if (unit)
  {
    for (i = 0; i < n; i++)
    {
      axpy(n - i - 1, -Xx[i], Ra + i + (i + 1) * ldR, ldR, Xx + i + 1, 1);
    }
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      Xx[i] /= Ra[i + i * ldR];
      axpy(n - i - 1, -Xx[i], Ra + i + (i + 1) * ldR, ldR, Xx + i + 1, 1);
    }
  }
}

void solve_rdrt(pmatrix R, pvector x)
{
  solve_upper(true, R, x);
  solve_diag(R, x);
  solve_uppertrans(true, R, x);
}

void block_rsolve(int rows, int cols, real *Ra, int ldR, real *Xa, int ldX)
{
  int j, k;

  for (k = cols; k-- > 0;)
  {
    for (j = 0; j < k; j++)
    {
      axpy(rows, -Ra[j + k * ldR], Xa + k * ldX, 1, Xa + j * ldX, 1);
    }
  }
}

void block_diagsolve(int rows, int cols, real *Ra, int ldR, real *Xa, int ldX)
{
  int i;

  for (i = 0; i < cols; i++)
  {
    scal(rows, 1.0 / Ra[i + i * ldR], Xa + i * ldX, 1);
  }
}

void addmul_rdrt(int rows, int cols, int k, real alpha, real *Aa, int ldA, real *Da, int ldD, real *Ba, int ldB, real *Xa, int ldX)
{
  int i;

  for (i = 0; i < k; i++)
  {
    ger(rows, cols, alpha * Da[i + i * ldD], Aa + i * ldA, 1, Ba + i * ldB, 1, Xa, ldX);
  }
}

/* ------------------------------------------------------------
 * Householder elementary reflections
 * ------------------------------------------------------------ */

static void refl_generate(int n, real *a, int inca, real *tau)
{
  real na, v1;

  na = nrm2(n, a, inca);

  if (na == 0.0)
    *tau = 0.0;
  else
  {
    *tau = 1.0 + fabs(a[0]) / na;

    if (a[0] < 0.0)
    {
      v1 = a[0] - na;
      a[0] = na;
    }
    else
    {
      v1 = a[0] + na;
      a[0] = -na;
    }

    scal(n - 1, 1.0 / v1, a + inca, inca);
  }
}

static void refl_apply(int rows, int cols, real tau,
                       const real *v, int incv,
                       real *work, real *A, int ldA)
{
  int j;

  /* Compute work = A^* v */
  for (j = 0; j < cols; j++)
    work[j] = A[j * ldA];

  gemv(true, rows - 1, cols, 1.0, A + 1, ldA,
       v + incv, incv, work, 1);

  /* Compute A -= tau v work^* */
  axpy(cols, -tau, work, 1, A, ldA);

  ger(rows - 1, cols, -tau,
      v + incv, incv, work, 1, A + 1, ldA);
}

/* ------------------------------------------------------------
 * QR decomposition
 * ------------------------------------------------------------ */

void qrdecomp(pmatrix A, real *tau, int lwork, real *work)
{
  int k, rows, cols, ldA;
  real *Aa = A->a;
  rows = A->rows;
  cols = A->cols;
  ldA = A->ld;

  assert(lwork >= cols - 1);

  for (k = 0; k < rows && k < cols; k++)
  {
    refl_generate(rows - k, Aa + k + k * ldA, 1, tau + k);

    refl_apply(rows - k, cols - k - 1, tau[k],
               Aa + k + k * ldA, 1, work, Aa + k + (k + 1) * ldA, ldA);
  }
}

void qsolve(bool trans, pmatrix A, const real *tau,
            pvector x)
{
  real work;
  int k, rows, cols, ldA;
  real *Aa = A->a;
  rows = A->rows;
  cols = A->cols;
  ldA = A->ld;
  real *xx = x->x;

  if (trans)
  {
    for (k = (rows < cols ? rows : cols); k-- > 0;)
    {
      refl_apply(rows - k, 1, tau[k],
                 Aa + k + k * ldA, 1, &work, xx + k, rows);
    }
  }
  else
  {
    for (k = 0; k < rows && k < cols; k++)
    {
      refl_apply(rows - k, 1, tau[k],
                 Aa + k + k * ldA, 1, &work, xx + k, rows);
    }
  }
}

void qexpand(pmatrix A, const real *tau, pmatrix Q, int lwork, real *work)
{
  int i, j, k;
  int rows, cols, ldA, ldQ;
  real *Aa = A->a;
  real *Qa = Q->a;
  rows = A->rows;
  cols = A->cols;
  ldA = A->ld;
  ldQ = Q->ld;

  assert(rows >= cols);
  assert(rows = Q->rows);
  assert(cols = Q->cols);

  for (j = 0; j < cols; j++)
  {
    for (i = 0; i < j; i++)
    {
      Qa[i + j * ldQ] = 0.0;
    }
    Qa[j + j * ldQ] = 1.0;
    for (i = j + 1; i < rows; i++)
    {
      Qa[i + j * ldQ] = 0.0;
    }
  }

  for (k = cols; k-- > 0;)
  {
    refl_apply(rows - k, cols - k, tau[k],
               Aa + k + k * ldA, 1, work, Qa + k + k * ldQ, ldQ);
  }
}

/* ------------------------------------------------------------
 * RDR^T decomposition for tridiagonal matrices
 * ------------------------------------------------------------ */

void decomp_rdrt_symtridiag(psymtridiag a)
{

  int n = a->dim;
  real *ax = a->a->x;
  real *bx = a->b->x;

  int k;

  for (k = n; k-- > 1;)
  {
    bx[k - 1] /= ax[k];
    ax[k - 1] -= bx[k - 1] * ax[k] * bx[k - 1];
  }
}

void eval_rdrt_symtridiag(psymtridiag a, pmatrix x)
{
  int n = a->dim;
  real *ax = a->a->x;
  real *bx = a->b->x;
  real *xa = x->a;
  int cols = x->cols;
  int ldX = x->ld;

  int i, k;

  assert(x->rows == n);

  // eval R^T
  for (k = n; k-- > 1;)
  {
    for (i = 0; i < cols; i++)
    {
      xa[k + i * ldX] += bx[k - 1] * xa[(k - 1) + i * ldX];
    }
  }

  // eval D
  for (k = 0; k < n; k++)
  {
    for (i = 0; i < cols; i++)
    {
      xa[k + i * ldX] *= ax[k];
    }
  }

  // eval R
  for (k = 0; k < n - 1; k++)
  {
    for (i = 0; i < cols; i++)
    {
      xa[k + i * ldX] += bx[k] * xa[(k + 1) + i * ldX];
    }
  }
}

void solve_rdrt_symtridiag(psymtridiag a, pmatrix x)
{

  int rows = a->dim;
  real *ax = a->a->x;
  real *bx = a->b->x;
  real *xx = x->a;
  int xld = x->ld;
  int m = x->cols;

  assert(x->rows == rows);

  int k, j;

  // solve R
  for (k = rows - 1; k-- > 0;)
  {
    for (j = 0; j < m; j++)
    {
      xx[k + j * xld] -= bx[k] * xx[(k + 1) + j * xld];
    }
  }

  // solve D
  for (k = 0; k < rows; k++)
  {
    for (j = 0; j < m; j++)
    {
      xx[k + j * xld] /= ax[k];
    }
  }

  // solve R^T
  for (k = 1; k < rows; k++)
  {
    for (j = 0; j < m; j++)
    {
      xx[k + j * xld] -= bx[k - 1] * xx[(k - 1) + j * xld];
    }
  }
}
