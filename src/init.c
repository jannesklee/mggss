/*****************************************************************************
 *                                                                           *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grid.h"                                 // grid struct            //
#include "routines.h"                             // gauss_seidel routine   //

void Output(Grid);


int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid          G;                                // grid                   //
  double        *u, *v;                           // solution u and rhs v   //
  unsigned int  n = 100;                          // mesh resolution        //
  double        eps0 = 1e-10;                     // error                  //
  double        eps=0.0;                          // error                  //
  unsigned int  i, j;                             // loop-indices           //

  //---------------------------- allocations -------------------------------//
  u = (double *) malloc((n+2)*(n+2)*sizeof(double));
  v = (double *) malloc((n+2)*(n+2)*sizeof(double));

  //--------------------------- initialization -----------------------------//
  // initial values
  for (i = 0; i < n+2; i++) {
    for (j = 0; j < n+2; j++) {
      v[i+j*(n+2)] = 10.*exp(-((i-50.)*(i-50.)+(j-50.)*(j-50.))
          *pow((1./(n+1.)/0.1),2.));
      u[i+j*(n+2)] = 0.0;
    }
  }

  // create grid with initial values for poisson equation
  G = Grid_Create();
  Grid_Set(G, u, v, n);

  //----------------------------- main loop --------------------------------//
  i = 0;
  do {
    MG_Method(G, &eps);
    i++;
  } while(eps > eps0);

  printf("%d \n \n \n \n", i);

  Output(G);

  //--------------------------- deallocations ------------------------------//
  free(u);
  free(v);
  Grid_Destroy(&G);
  return 0;
}
