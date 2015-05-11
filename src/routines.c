#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "grid.h"
#include "routines.h"

//------------------------- multigrid-method -------------------------------//
void MG_Method(Grid G, double *eps) {
  double *u, *u_save, *v, *v_save;
  unsigned int n;
  unsigned int i, j;

  n = Grid_Get_n(G);
  u = Grid_Get_u(G);

  u_save = (double *) calloc((n+2)*(n+2), sizeof(double));
  v_save = (double *) calloc((n+2)*(n+2), sizeof(double));

  //----------------------------- finer grid -------------------------------//
  // pre-smoothing
  Gauss_Seidel(G);

  // save u and v in temporary variables
  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      u_save[i+j*(n+2)] = u[i+j*(n+2)];
      v_save[i+j*(n+2)] = v[i+j*(n+2)];
    }
  }
  // calculate residual r and write it in v
  AddEval(1.0, G);
  //------------------------------------------------------------------------//

  //--------------------------- coarser grid -------------------------------//
  // restricts only v, because it is assumed to use the residue!
  Restriction(G);

  // set u (now the error) to zero
  n = Grid_Get_n(G);
  u = Grid_Get_u(G);
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      u[i+j*(n+2)] = 0.0;
    }
  }
  Grid_Set_u(G, u);

  // recalculate on coarser grid (TODO: here eventuelly not an iterative scheme?)
  Gauss_Seidel(G);

  //------------------------------------------------------------------------//

  //--------------------------- finer grid ---------------------------------//
  // prolongate array
  Prolongation(G);

  // finally calculate solution
  u = Grid_Get_u(G);
  n = Grid_Get_n(G);
  v = Grid_Get_v(G);
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      u[i+j*(n+2)] = u_save[i+j*(n+2)] + u[i+j*(n+2)];
      v[i+j*(n+2)] = v_save[i+j*(n+2)];
    }
  }

  // post-smoothing
  Grid_Set_v(G, v);
  Grid_Set_u(G, u);
  Gauss_Seidel(G);

  // calculate error and set back to zero before doing so
  *eps = 0.0;
  for (j = 1; j < n+1; j++) {
    for (i = 1; i < n+1; i++) {
      (*eps) = fmax((*eps),fabs(u[i+j*(n+2)]-u_save[i+j*(n+2)]));
    }
  }

  free(u_save);
  free(v_save);
}

//---------------------------- gauss-seidel --------------------------------//
void Gauss_Seidel (Grid G) {
  double *u, *v;                                  // solution, rhs          //
  double h;                                       // stepwidth              //
  unsigned int n;                                 // number of grid points  //
  unsigned int i, j;

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  // main-part of gauss-seidel
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      u[i+j*(n+2)] = 0.25*(h*h*v[i+j*(n+2)] + u[(i+1)+j*(n+2)]
          + u[(i-1)+j*(n+2)] + u[i+(j+1)*(n+2)] + u[i+(j-1)*(n+2)]);
    }
  }

  Grid_Set_u(G,u);
}

//---------------------------- restriction ---------------------------------//
void Restriction (Grid G) {
  double *u, *u_c;                                // values at              //
  double *v, *v_c;                                // current & coarser grid //
  unsigned int n, n_c;                            // number of grid-points  //
  unsigned int i, j;

  // get and allocate values on current and coarser grid
  n = Grid_Get_n(G);
  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  n_c = (n-1)/2;
  u_c = (double *) calloc((n_c+2)*(n_c+2), sizeof(double));
  v_c = (double *) calloc((n_c+2)*(n_c+2), sizeof(double));

  // loop over even, even combination in order to define coarse grid
  for (j = 2; j < n+1; j = j+2) {
    for (i = 2; i < n+1; i = i+2) {
      v_c[i/2+j/2*(n_c+2)] = 0.25*(v[i+j*(n+2)] + 0.5*(v[(i+1)+j*(n+2)] +
          v[(i-1)+j*(n+2)] + v[i+(j+1)*(n+2)] + v[i+(j-1)*(n+2)]) +
          0.25*(v[(i+1)+(j+1)*(n+2)] + v[(i-1)+(j+1)*(n+2)] +
          v[(i+1)+(j-1)*(n+2)] + v[(i-1)+(j-1)*(n+2)]));
    }
  }

  // copy values in first square of v array
  for (j = 1; j < n_c+1; j++) {
    for (i = 1; i < n_c+1; i++) {
      v[i+j*(n_c+2)] = v_c[i+j*(n_c+2)];
    }
  }

  // overwrite grid
  Grid_Set(G, u, v, n_c);
  free(u_c);
  free(v_c);
}

//--------------------------- prolongation ---------------------------------//
void Prolongation(Grid G) {
  double *v, *v_f;                                // values at              //
  double *u, *u_f;                                // current & finer grid   //
  unsigned int n, n_f;                            // number of grid points  //
  unsigned int i, j;

  // get and allocate values on current and finer grid
  n = Grid_Get_n(G);
  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  n_f = n*2+1;
  u_f = (double *) calloc((n_f+2)*(n_f+2), sizeof(double));
  v_f = (double *) calloc((n_f+2)*(n_f+2), sizeof(double));

  // pad fine temporary array with zeros at non-defined places
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      if ((j%2 == 0) && (i%2 == 0)) {
        u_f[i+j*(n_f+2)] = u[i/2+j/2*(n+2)];
      }else {
        u_f[i+j*(n_f+2)] = 0.0;
      }
    }
  }

  // set everything to zero to avoid side-effects
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      v[i+j*(n+2)] = 0.0;
      u[i+j*(n+2)] = 0.0;
    }
  }

  Grid_Set(G,u_f,v_f,n_f);

  // loop over all elements in finer array
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      v[i+j*(n_f+2)] = v_f[i+j*(n_f+2)] + 0.5*(v_f[(i+1)+j*(n_f+2)]
          + v_f[(i-1)+j*(n_f+2)] + v_f[i+(j+1)*(n_f+2)]
          + v_f[i+(j-1)*(n_f+2)]) + 0.25*(v_f[(i+1)+(j+1)*(n_f+2)]
          + v_f[(i+1)+(j-1)*(n_f+2)] + v_f[(i-1)+(j+1)*(n_f+2)]
          + v_f[(i-1)+(j-1)*(n_f+2)]);
      u[i+j*(n_f+2)] = u_f[i+j*(n_f+2)] + 0.5*(u_f[(i+1)+j*(n_f+2)]
          + u_f[(i-1)+j*(n_f+2)] + u_f[i+(j+1)*(n_f+2)]
          + u_f[i+(j-1)*(n_f+2)]) + 0.25*(u_f[(i+1)+(j+1)*(n_f+2)]
          + u_f[(i+1)+(j-1)*(n_f+2)] + u_f[(i-1)+(j+1)*(n_f+2)]
          + u_f[(i-1)+(j-1)*(n_f+2)]);
    }
  }

  // overwrite grid
  Grid_Set(G,u,v,n_f);

  // deallocation
  free(v_f);
  free(u_f);
}

//----------------------------- addeval ------------------------------------//
void AddEval(double alpha, Grid G) {
  double *u, *v;
  double h;
  unsigned int n;
  unsigned int i,j;

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  alpha = alpha/(h*h);

  for (j = 1; j < n+1; j++){
    for (i = 1; i < n+1; i++){
      v[i+j*(n+2)] = v[i+j*(n+2)] - alpha*(4.*u[i+j*(n+2)] -
                       u[(i+1)+j*(n+2)] - u[(i-1)+j*(n+2)] -
                       u[i+(j+1)*(n+2)] - u[i+(j-1)*(n+2)]);
    }
  }

  Grid_Set_v(G,v);
}

//-------------------------------- output ----------------------------------//
void Output(Grid G){
  double *u, *v;
  unsigned int i,j,n;

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  n = Grid_Get_n(G);

  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      printf("%.18f, %.18f, %.18f, %.18f\n",
          1./(n-1)*(i-1.), 1./(n-1)*(j-1.), u[i+j*(n+2)], v[i+j*(n+2)]);
    }
    printf("\n");
  }
  printf("\n");
}
