#include <stdlib.h>
#include "math.h"                                 // fmax, fabs             //
#include "grid.h"
#include "routines.h"


//---------------------------- gauss-seidel --------------------------------//
void Gauss_Seidel (Grid G, double *eps) {
  // declarations
  double *u, *v;                                  // sol, rhs               //
  double h;                                       // stepwidth              //
  unsigned int n;                                 // number of grid points  //
  double tmp;                                     // temporary variable     //
  unsigned int i, j;                              // indices                //

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);


  // main-part of gauss-seidel
  *eps = 0.0;
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      tmp = u[i+j*(n+2)];
      u[i+j*(n+2)] = h*h/4.*(v[i+j*(n+2)] + u[(i+1)+j*(n+2)]
          + u[(i-1)+j*(n+2)] + u[i+(j+1)*(n+2)] + u[i+(j-1)*(n+2)]);
      (*eps) = fmax((*eps),fabs(u[i+j*(n+2)]-tmp));
    }
  }

  // write solution in grid
  Grid_Set_u(G,u);
}

//---------------------------- restriction ---------------------------------//
void Restriction (Grid G) {
  double *r, *u2, *r2;
//  double *v2;
  unsigned int n, n2, i, j;

  // get old grid size
  n = Grid_Get_n(G);

  // allocate variables for new (coarser) grid
  n2 = n/2;
  u2 = (double *) malloc((n2+2)*(n2+2)*sizeof(double));
//  v2 = (double *) malloc((n2+2)*(n2+2)*sizeof(double));
  r2 = (double *) malloc((n2+2)*(n2+2)*sizeof(double));

  // safe residuel of old grid
  r = AddEval(1.0,G);                             // get residual           //

  // loop over even, even combinations
  for (j = 2; j < n+1; j=j+2) {
    for (i = 2; i < n+1; i=i+2) {
      r2[i/2+j/2*(n2+2)] = 0.25 * (r[i+j*(n+2)] + 0.5*(r[(i+1)+j*(n+2)] +
          r[(i-1)+j*(n+2)] + r[i+(j+1)*(n+2)] + r[i+(j-1)*(n+2)]) +
          0.25*(r[(i+1)+(j+1)*(n+2)] + r[(i-1)+(j+1)*(n+2)] +
          r[(i+1)+(j-1)*(n+2)] + r[(i-1)+(j-1)*(n+2)]));
    }
  }

  //TODO: Zusammenhang residuum und Lösung. Hier umrechnen oder später?
  // overwrite grid
  Grid_Set(G,u2,r2,n2);
}

//void Prolongation(Grid G) {
//  double *r, *u2, *v2, *r2;
//  unsigned int n, n2, i, j;
//
//  u = Grid_Get_u(G);
//  v = Grid_Get_v(G);
//  n = Grid_Get_n(G);
//
//
//}

double * AddEval(double alpha, Grid G) {
  //--------------------------- declarations -------------------------------//
  double *u, *v;                                  // sol, rhs               //
  double h;                                       // stepwidth              //
  unsigned int n,i,j;                             // number of grid points  //

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  alpha = alpha/(h*h);

  for (j = 1; j < n+1; j++){
    for (i = 1; i < n+1; i++){
      v[i+j*(n+2)] = v[i+j*(n+2)] + alpha* (4.*u[i+j*(n+2)] -
                      u[(i+1)+j*(n+2)] - u[(i-1)+j*(n+2)] -
                      u[i+(j+1)*(n+2)] - u[i+(j-1)*(n+2)]);
    }
  }

  return v;
}
