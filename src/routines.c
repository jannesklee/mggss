
#include "math.h"                                 // fmax, fabs             //
#include "grid.h"

void Gauss_Seidel (Grid G) {
  double *u, *v;                                  // sol, rhs               //
  double h;                                       // stepwidth              //
  unsigned int n;                                 // number of grid points  //
  double tmp;                                     // temporary variable     //
  double eps;

  unsigned int i, j;                              // indices                //

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  // TODO: error here oder outside?
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      tmp = u[i+j*(n+2)];
      u[i+j*(n+2)] = h*h/4.*(v[i+j*(n+2)] + u[(i+1)+j*(n+2)]
          + u[(i-1)+j*(n+2)] + u[i+(j+1)*(n+2)] + u[i+(j-1)*(n+2)]);
      eps = fmax(eps,fabs(u[i+j*(n+2)]-tmp));
    }
  }

  Grid_Set_u(G,u);
}

void Restriction (Grid G) {
  // TODO: write this in jdetail
  double *u, *v;
//  double h;
  unsigned int n;

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
//  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  n = n/2;

  //restriction part

  Grid_Set(G,u,v,n);
}

double * AddEval(double alpha, Grid G) {
  double *u, *v;                                  // sol, rhs               //
  double h;                                       // stepwidth              //
  unsigned int n,i,j;                             // number of grid points  //

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  for (j = 1; j < n+1; j++){
    for (i = 1; i < n+1; i++){
      v[i+j*(n+2)] = v[i+j*(n+2)] +
        alpha/(h*h)*(4*u[i+j*(n+2)] - u[(i+1)+j*(n+2)] - u[(i-1)+j*(n+2)] -
                                      u[i+(j+1)*(n+2)] - u[i+(j-1)*(n+2)]);
    }
  }

  return v;
}
