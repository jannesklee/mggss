/*****************************************************************************
 *                                                                           *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grid.h"                                 // grid struct            //
#include "routines.h"                             // gauss_seidel routine   //

void initialize(double *, double *, unsigned int, unsigned int);

int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid          G;                                // grid                   //
  double        *u, *v;                           // solution u and rhs v   //
  unsigned int  n = 199;                          // res. has to be odd     //
  double        eps0 = 1e-10;                      // error-limit            //
  double        eps  = 0.0;                       // error                  //
  unsigned int  k = 1;                            // 1: const. dens. circle //
                                                  // 2: exponential distr.  //

  //---------------------------- allocations -------------------------------//
  u = (double *) calloc((n+2)*(n+2), sizeof(double));
  v = (double *) calloc((n+2)*(n+2), sizeof(double));

  //--------------------------- initialization -----------------------------//
  initialize(u, v, n, k);

  // create grid with initial values for poisson equation
  G = Grid_Create();
  Grid_Set(G, u, v, n);

  //----------------------------- main loop --------------------------------//
  do {
    MG_Method(G, &eps);
  } while(eps > eps0);
  Output(G);

  //--------------------------- deallocations ------------------------------//
  free(u);
  free(v);
  Grid_Destroy(&G);
  return 0;
}

void initialize(double *u, double *v, unsigned int n, unsigned int k) {
  double        r;                                // radius for init.       //
  unsigned int  i,j;

  if (k== 1) { // enhanced density on circle
    // initial values & boundaries
    for (i = 0; i < n+2; i++) {
      for (j = 0; j < n+2; j++) {
        r = sqrt((1./((n-1.)*(n-1.))*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2))));
        if (r < 0.1) {
          v[i+j*(n+2)] = 100.;
        } else {
          v[i+j*(n+2)] = 0.0;
        }
        u[i+j*(n+2)] = 0.0;
      }
    }
  }
  else if (k == 2) { // gaussian distribution
    for (i = 0; i < n+2; i++) {
      for (j = 0; j < n+2; j++) {
        r = sqrt((1./((n-1.)*(n-1.))*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2))));
        v[i+j*(n+2)] = exp(-1000*r*r);
      }
    }
  }
  else {
    // do nothing
  }
}
