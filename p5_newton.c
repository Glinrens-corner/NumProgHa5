#include "factorizations.h"
#include "matrix.h"
#include "basic.h"
#include "miniblas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief function pointer type for evaluation of the function 'f'
 *
 * @param data untyped pointer to data structure containing
 *             necessary information of evaluating f.
 * @param z Position at which we want to evaluate f.
 * @param y Result of the evaluation: y = f(z).
 */
typedef void (*eval)(void *data, pvector z, pvector y);

/**
 * @brief function pointer type for evaluation of the Jacobian of the function 'f'
 *
 * @param data untyped pointer to data structure containing
 *             necessary information of evaluating Df.
 * @param z Position at which we want to evaluate Df.
 * @param y Result of the evaluation: Df = Df(z).
 */
typedef void (*eval_Df)(void *data, pvector z, pmatrix Df);

struct _grav
{
  real alpha2; // Squared angular velocity
  real gamma;  // Gravitational constant
  int n;       // Number of point masses
  real **x;    // Locations of point masses
  real *m;     // Point masses
};
typedef struct _grav grav;
typedef grav *pgrav;



/*
 *  copy one vector from another one.
 */
static void vec_copy(int n, real* vec, const real* source){
  for (int ientry =0; ientry < n ; ientry++){
    vec[ientry]=source[ientry];
  }
}

/*
 *  substract one vector from another one
 * 
 *  vec is the vector to be updated
 */
static void vec_substract(int n, real* vec, const real* subtrahend){
  for (int ientry =0; ientry < n ; ientry++){
    vec[ientry] -= subtrahend[ientry];
  }
}

/*
 *  add  one vector to another one
 * 
 *  vec is the vector to be updated
 */
static void vec_add(int n, real* vec, const real* summand){
  for (int ientry =0; ientry < n ; ientry++){
    vec[ientry] += summand[ientry];
  }
}

void eval_f_grav(void *data, pvector z, pvector fz)
{
  pgrav ourData = (pgrav) data;
  
  real *masses = ourData->m;
  real **locations = ourData->x;
  real * res = fz->x;
  vec_copy(3, res, z->x);
  scal(3, ourData->alpha2, res, 1); // now alpha^2*z is on fz
  
  for (int ibody=0 ; ibody <ourData->n; ibody++ ){
    real difference[3]; 
    vec_copy(3, difference, locations[ibody]);
    vec_substract(3, difference, z->x);
    real norm = nrm2(3,difference, 1);
    real scaling_factor = ourData->gamma*masses[ibody]/(norm*norm*norm);
    scal(3,scaling_factor, difference, 1);
    vec_add(3,res,difference);
  }
}


void eval_Df_grav(void *data, pvector z, pmatrix Dfz)
{
  pgrav ourData = (pgrav) data;
  
  real *masses = ourData->m;
  real **locations = ourData->x;
  real * invec = z->x;
  //                         ( alpha^2    0       0    )
  // (d_zi (alpha^2 zj)) z = (   0     alpha^2    0    ) z = alpha^2 z
  //                         (   0        0    alpha^2 )

  clear_matrix(Dfz);
  for(int idim=0; idim < 3 ; idim ++) {
    Dfz->a[idim*3+idim] += ourData->alpha2;
  }
  
  for (int ibody=0 ; ibody <ourData->n; ibody++ ){
    real difference[3]; 
    vec_copy(3, difference, locations[ibody]);
    vec_substract(3, difference, invec);
  
    real norm = nrm2(3,difference, 1);
   
    for(int idim=0; idim < 3 ; idim ++) {
      Dfz->a[idim*3+idim] -= ourData->gamma*masses[ibody]/(norm*norm*norm);
    }


    real factor = ourData->gamma*masses[ibody]*3.0/(norm*norm*norm*norm*norm);
 
    for (int icol=0; icol <3 ; icol++) {
      
      for (int irow=0; irow <3 ; irow++) {
	Dfz->a[icol*3+irow] += factor * difference[icol]* difference[irow];
      }
    } 
  }
}

unsigned int
iterate_newton(real eps, int maxiter,
               pvector fz, void *data_f, eval eval_f,
               pmatrix Dfz, void *data_Df, eval_Df eval_Df,
               pvector z)
{

  
  real stepwidth = 1.0;
  unsigned int niterations = 0;
  real error = eps+1.0;
  while(error > eps && niterations < 40 ){
    eval_f(data_f, z, fz);
    eval_Df(data_Df,z,Dfz);
    decomp_rdrt(Dfz->a, Dfz->ld,Dfz->rows);
    solve_rdrt(Dfz, fz);
    error = nrm2(3, fz->x, 1)/nrm2(3, z->x, 1);
    axpy(3,-1.0*stepwidth , fz->x,1,z->x,1);
    niterations ++;
  }
  return niterations;
}

int main(void)
{
  pmatrix Dfz;
  pvector z, fz;
  real *m, **x, eps, alpha2, gamma, dist;
  unsigned int i, n, dim, maxiter;
  eval eval_f;
  eval_Df eval_Df;
  pgrav g;
  pstopwatch sw;
  real t;
  real error;

  dim = 3;       // working in 3D space
  n = 2;         // Using only two gravitational objects: sun and earth.
  maxiter = 100; // max number of iterations for Newtons method.
  eps = 1.0e-12; // relative accuracy for Newtons method.

  sw = new_stopwatch();

  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  Setting up point masses:\n");
  printf("--------------------------------------------------------------------------------\n" NORMAL);

  start_stopwatch(sw);

  z = new_vector(dim);             // position vector
  fz = new_vector(dim);            // force vector
  Dfz = new_matrix(dim, dim);      // Jacobian matrix
  g = (pgrav)malloc(sizeof(grav)); // structure for storing all necessary information

  m = (real *)malloc(n * sizeof(real));    // masses
  x = (real **)malloc(n * sizeof(real *)); // positions
  x[0] = (real *)malloc(n * dim * sizeof(real));
  for (i = 1; i < n; i++)
  {
    x[i] = x[0] + i * dim;
  }

  
  dist = 1.496e11;     // distance sun-earth
  gamma = 6.67430e-11; // gravitational constant
  m[0] = 1.989e30;     // mass sun
  m[1] = 5.9723e24;    // mass earth
  alpha2 = gamma * (m[0] + m[1]) / (dist * dist * dist);

  x[0][0] = -dist * m[1] / (m[0] + m[1]); // position sun
  x[0][1] = 0.0;
  x[0][2] = 0.0;

  x[1][0] = dist * m[0] / (m[0] + m[1]); // position earth
  x[1][1] = 0.0;
  x[1][2] = 0.0;

  // clear auxiliary objects
  clear_vector(z);
  clear_vector(fz);
  clear_matrix(Dfz);

  // fill data structure
  eval_f = &eval_f_grav;
  eval_Df = &eval_Df_grav;
  g->n = n;
  g->alpha2 = alpha2;
  g->gamma = gamma;
  g->x = x;
  g->m = m;

  t = stop_stopwatch(sw);
  printf("  %.3f ms\n", t * 1.0e3);

  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  Performing Newton iteration for Lagrangian points:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);
  
  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  L_1:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  z->x[0] = 1.0e11;
  z->x[1] = 0.0;
  z->x[2] = 0.0;
  printf("Initial point: (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  start_stopwatch(sw);
  i = iterate_newton(eps, maxiter, fz, g, eval_f, Dfz, g, eval_Df, z);
  t = stop_stopwatch(sw);
  printf("  %u iterations needed.\n", i);
  printf("  %.3f ms\n", t * 1.0e3);
  printf("L_1 = (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  error = (148108114680.8629150391 - z->x[0]) * (148108114680.8629150391 - z->x[0]);
  error += (0.0 - z->x[1]) * (0.0 - z->x[1]);
  error += (0.0 - z->x[2]) * (0.0 - z->x[2]);
  error = sqrt(error) / (148108114680.8629150391 * 148108114680.8629150391 + 0.0 * 0.0 + 0.0 * 0.0);
  printf("rel. error: %s%.5e%s\n", error < 1.0e-15 ? BGREEN : error < 1.0e-12 ? BYELLOW
                                                                              : BRED,
         error, NORMAL);
  printf("\n");
  
  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  L_2:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  z->x[0] = 1.5e11;
  z->x[1] = 0.0;
  z->x[2] = 0.0;
  printf("Initial point: (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  start_stopwatch(sw);
  i = iterate_newton(eps, maxiter, fz, g, eval_f, Dfz, g, eval_Df, z);
  t = stop_stopwatch(sw);
  printf("  %u iterations needed.\n", i);
  printf("  %.3f ms\n", t * 1.0e3);
  printf("L_2 = (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  error = (151100965998.0512084961 - z->x[0]) * (151100965998.0512084961 - z->x[0]);
  error += (0.0 - z->x[1]) * (0.0 - z->x[1]);
  error += (0.0 - z->x[2]) * (0.0 - z->x[2]);
  error = sqrt(error) / (151100965998.0512084961 * 151100965998.0512084961 + 0.0 * 0.0 + 0.0 * 0.0);
  printf("rel. error: %s%.5e%s\n", error < 1.0e-15 ? BGREEN : error < 1.0e-12 ? BYELLOW
                                                                              : BRED,
         error, NORMAL);
  printf("\n");
  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  L_3:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  z->x[0] = -3.0e11;
  z->x[1] = 0.0;
  z->x[2] = 0.0;
  printf("Initial point: (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  start_stopwatch(sw);
  i = iterate_newton(eps, maxiter, fz, g, eval_f, Dfz, g, eval_Df, z);
  t = stop_stopwatch(sw);
  printf("  %u iterations needed.\n", i);
  printf("  %.3f ms\n", t * 1.0e3);
  printf("L_3 = (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  error = (-149600187164.4413146973 - z->x[0]) * (-149600187164.4413146973 - z->x[0]);
  error += (0.0 - z->x[1]) * (0.0 - z->x[1]);
  error += (0.0 - z->x[2]) * (0.0 - z->x[2]);
  error = sqrt(error) / (-149600187164.4413146973 * -149600187164.4413146973 + 0.0 * 0.0 + 0.0 * 0.0);
  printf("rel. error: %s%.5e%s\n", error < 1.0e-15 ? BGREEN : error < 1.0e-12 ? BYELLOW
                                                                              : BRED,
         error, NORMAL);
  printf("\n");
  
  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  L_4:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  z->x[0] = 1.0e11;
  z->x[1] = 1.0e11;
  z->x[2] = 0.0;
  printf("Initial point: (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  start_stopwatch(sw);
  i = iterate_newton(eps, maxiter, fz, g, eval_f, Dfz, g, eval_Df, z);
  t = stop_stopwatch(sw);
  printf("  %u iterations needed.\n", i);
  printf("  %.3f ms\n", t * 1.0e3);
  printf("L_4 = (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);

  error = (74799550930.9104614258 - z->x[0]) * (74799550930.9104614258 - z->x[0]);
  error += (129557400332.1854400635 - z->x[1]) * (129557400332.1854400635 - z->x[1]);
  error += (0.0 - z->x[2]) * (0.0 - z->x[2]);
  error = sqrt(error) / (74799550930.9104614258 * 74799550930.9104614258 + 129557400332.1854400635 * 129557400332.1854400635 + 0.0 * 0.0);
  printf("rel. error: %s%.5e%s\n", error < 1.0e-15 ? BGREEN : error < 1.0e-12 ? BYELLOW
                                                                              : BRED,
         error, NORMAL);
  printf("\n");

  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  L_5:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  z->x[0] = 1.0e11;
  z->x[1] = -1.0e11;
  z->x[2] = 0.0;
  printf("Initial point: (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  start_stopwatch(sw);
  i = iterate_newton(eps, maxiter, fz, g, eval_f, Dfz, g, eval_Df, z);
  t = stop_stopwatch(sw);
  printf("  %u iterations needed.\n", i);
  printf("  %.3f ms\n", t * 1.0e3);
  printf("L_5 = (%.3f, %.3f, %.3f)\n", z->x[0], z->x[1], z->x[2]);
  error = (74799550930.9104614258 - z->x[0]) * (74799550930.9104614258 - z->x[0]);
  error += (-129557400332.1854400635 - z->x[1]) * (-129557400332.1854400635 - z->x[1]);
  error += (0.0 - z->x[2]) * (0.0 - z->x[2]);
  error = sqrt(error) / (74799550930.9104614258 * 74799550930.9104614258 + -129557400332.1854400635 * -129557400332.1854400635 + 0.0 * 0.0);
  printf("rel. error: %s%.5e%s\n", error < 1.0e-15 ? BGREEN : error < 1.0e-12 ? BYELLOW
                                                                              : BRED,
         error, NORMAL);
  printf("\n");

  printf(BCYAN "--------------------------------------------------------------------------------\n");
  printf("  Cleaning up:\n");
  printf("--------------------------------------------------------------------------------\n\n" NORMAL);

  free(g);
  free(x[0]);
  free(x);
  free(m);
  del_matrix(Dfz);
  del_vector(z);
  del_vector(fz);
  del_stopwatch(sw);

  return EXIT_SUCCESS;
}
