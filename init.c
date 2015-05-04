/*****************************************************************************
 *                                                                           *
 ****************************************************************************/

#include <stdio.h>                                // standard input-output  //
#include <stdlib.h>                               //                        //
#include "grid.h"                                 // linear systems struct  //


int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid G;                                         // system of linear eqs.  //
  double * u, * v;                                // solution u and rhs v   //
  unsigned int n = 9;                             // mesh resolution        //
  double eps = 1e-3;                              // error                  //
  double h;                                       // stepwidth              //
  unsigned int i, j;                              // loop-indices           //

  u = (double *) malloc((n+2)*(n+2)*sizeof(double));
  v = (double *) malloc((n+2)*(n+2)*sizeof(double));

  //--------------------------- initialization -----------------------------//
  // initial values
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      u[i+j*(n+2)] = 4.0;
    }
  }
  // boundaries
  for (i = 0; i < n+2; i++) {
    u[i+(n+1)*(n+2)]  = 0.0;
    u[(n+1)+i*(n+2)]  = 0.0;
    u[i]              = 0.0;
    u[i*(n+2)]        = 0.0;
  }

  // create grid with initial values for poisson equation
  G = Grid_Create();
  Grid_Set(G, u, v, n);
  //------------------------------------------------------------------------//

  v = Grid_Get_u(G);

  for (i = 0; i < n+2; i++) {
    for (j = 0; j < n+2; j++) {
      printf("%f\n", v[i+j*(n+2)]);
    }
    printf("\n");
  }

//  while (norm(u) < eps) {
//    Gauss_Seidel(G);
//    u = Grid_Get_u(G);
//  }


  Grid_Destroy(&G);
  return 0;
}

//double norm(u)
//{
//  //TODO: Find a useful norm
//  return u[2];
//}
