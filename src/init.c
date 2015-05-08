/*****************************************************************************
 *                                                                           *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grid.h"                                 // grid struct            //
#include "routines.h"                             // gauss_seidel routine   //


int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid          G;                                // grid                   //
  double        *u, *v;                           // solution u and rhs v   //
  unsigned int  n = 200;                          // mesh resolution        //
  double        eps0 = 1e-11;                     // error                  //
  double        eps  = 0.0;                       // error                  //
  unsigned int  i, j;                             // loop-indices           //
  double        r;                                // radiuas for init.      //

  //---------------------------- allocations -------------------------------//
  u = (double *) malloc((n+2)*(n+2)*sizeof(double));
  v = (double *) malloc((n+2)*(n+2)*sizeof(double));

  //--------------------------- initialization -----------------------------//
  // initial values
  for (i = 0; i < n+2; i++) {
    for (j = 0; j < n+2; j++) {
      r = sqrt((1./((n-1.)*(n-1.))*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2))));
      if (r < 0.1) {
        v[i+j*(n+2)] = 1000.;
      } else {
        v[i+j*(n+2)] = 0.0;
      }
      u[i+j*(n+2)] = 0.0;
    }
  }

  // create grid with initial values for poisson equation
  G = Grid_Create();                              // fine grid              //
  Grid_Set(G, u, v, n);

  //----------------------------- main loop --------------------------------//
  i = 0;
  do {
    MG_Method(G, &eps);
    i++;
  } while(eps > eps0);
  Output(G);
  //--------------------------- deallocations ------------------------------//

  free(u);
  free(v);
  Grid_Destroy(&G);
  return 0;
}
