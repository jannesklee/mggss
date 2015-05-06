/*****************************************************************************
 *                                                                           *
 ****************************************************************************/

#include <stdlib.h>                               // allocate, ...          //
#include <math.h>                                 // math library           //
#include "grid.h"                                 // linear systems struct  //
#include "routines.h"                             // gauss_seidel routine   //

void Output(Grid);


int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid          G;                                // grid                   //
  double        * u, * v;                         // solution u and rhs v   //
  unsigned int  n = 99;                           // mesh resolution        //
  //double        eps0 = 1e-10;                     // error                  //
  double        eps;                              // stepwidth              //
  unsigned int  i, j;                             // loop-indices           //


  //---------------------------- allocations -------------------------------//
  u = (double *) malloc((n+2)*(n+2)*sizeof(double));
  v = (double *) malloc((n+2)*(n+2)*sizeof(double));

  // TODO: Initialisierung in extra funktion
  //--------------------------- initialization -----------------------------//
  // initial values
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      u[i+j*(n+2)] = 4.0;
      v[i+j*(n+2)] = exp(-pow((1./(n+1.)*(i-50.)*(j-50.)),2.));
    }
  }
  // boundaries
  for (i = 0; i < n+2; i++) {
    u[i+(n+1)*(n+2)]  = 0.0;
    u[(n+1)+i*(n+2)]  = 0.0;
    u[i]              = 0.0;
    u[i*(n+2)]        = 0.0;
    v[i+(n+1)*(n+2)]  = 0.0;
    v[(n+1)+i*(n+2)]  = 0.0;
    v[i]              = 0.0;
    v[i*(n+2)]        = 0.0;
  }

  // create grid with initial values for poisson equation
  G = Grid_Create();
  Grid_Set(G, u, v, n);


  Output(G);

  eps = 0.0;
  Gauss_Seidel(G, &eps);

  Output(G);

  Restriction(G);

  Output(G);

  eps = 0.0;
  Gauss_Seidel(G, &eps);

  Output(G);

  Prolongation(G);

  Output(G);


  Grid_Destroy(&G);
  return 0;
}
