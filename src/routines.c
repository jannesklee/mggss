#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "grid.h"
#include "routines.h"

//------------------------- multigrid-method -------------------------------//
void MG_Method(Grid G, double *eps) {
  double *u, *v, *u_c, *v_c;
  double *u_save, *v_save;
  unsigned int n;
  unsigned int i, j;
  Grid H;

  n = Grid_Get_n(G);
  u = Grid_Get_u(G);

  u_c = (double *) calloc((n/2+2)*(n/2+2), sizeof(double));
  v_c = (double *) calloc((n/2+2)*(n/2+2), sizeof(double));
  u_save = (double *) calloc((n+2)*(n+2), sizeof(double));
  v_save = (double *) calloc((n+2)*(n+2), sizeof(double));

  // create coarse grid
  H = Grid_Create();
  Grid_Set(H,u_c,v_c,n/2);

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
  AddEval(-1.0, G);
  //------------------------------------------------------------------------//
  *eps = MaxNorm(Grid_Get_v(G),Grid_Get_n(G));

  //--------------------------- coarser grid -------------------------------//
  // restricts only v, because it is assumed to use the residue!
  Restriction(G, H);

  // recalculate on coarser grid
  Gauss_Seidel(H);
  //------------------------------------------------------------------------//

  // here a recursive call of Multigrid has to be made

  //--------------------------- finer grid ---------------------------------//
  // prolongate array from coarser grid H to finer grid G
  Prolongation(H, G);

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

//  // calculate error and set back to zero before doing so
//  *eps = 0.0;
//  for (j = 1; j < n+1; j++) {
//    for (i = 1; i < n+1; i++) {
//      (*eps) = fmax((*eps),fabs(u[i+j*(n+2)]-u_save[i+j*(n+2)]));
//    }
//  }

  free(u_save);
  free(v_save);
  free(u_c);
  free(v_c);
  Grid_Destroy(&H);
}

//---------------------------- gauss-seidel --------------------------------//
void Gauss_Seidel (Grid G) {
  double *u, *v;                                  // solution, rhs          //
  double h, h2;                                   // stepwidth              //
  unsigned int n;                                 // number of grid points  //
  unsigned int i, j;

  u = Grid_Get_u(G);
  v = Grid_Get_v(G);
  h = Grid_Get_h(G);
  n = Grid_Get_n(G);

  h2 = h*h;
  // main-part of gauss-seidel
  for (i = 1; i < n+1; i++) {
    for (j = 1; j < n+1; j++) {
      u[i+j*(n+2)] = 0.25*(h2*v[i+j*(n+2)] + u[(i+1)+j*(n+2)]
          + u[(i-1)+j*(n+2)] + u[i+(j+1)*(n+2)] + u[i+(j-1)*(n+2)]);
    }
  }

  Grid_Set_u(G,u);
}

//---------------------------- restriction ---------------------------------//
void Restriction (Grid G, Grid H) {               // finer and coarser grid //
  double *u_c;                                    // values at              //
  double *v, *v_c;                                // current & coarser grid //
  unsigned int n, n_c;                            // number of grid-points  //
  unsigned int i, j;

  // get and allocate values on current and coarser grid
  n = Grid_Get_n(G);
  v = Grid_Get_v(G);
  u_c = Grid_Get_u(H);
  v_c = Grid_Get_v(H);
  n_c = Grid_Get_n(H);

  // loop over even, even combination in order to define coarse grid
  for (j = 2; j < n+1; j = j+2) {
    for (i = 2; i < n+1; i = i+2) {
      v_c[i/2+j/2*(n_c+2)] = 0.25*(v[i+j*(n+2)] + 0.5*(v[(i+1)+j*(n+2)] +
          v[(i-1)+j*(n+2)] + v[i+(j+1)*(n+2)] + v[i+(j-1)*(n+2)]) +
          0.25*(v[(i+1)+(j+1)*(n+2)] + v[(i-1)+(j+1)*(n+2)] +
          v[(i+1)+(j-1)*(n+2)] + v[(i-1)+(j-1)*(n+2)]));
    }
  }

//  // copy values in first square of v array
//  for (j = 1; j < n_c+1; j++) {
//    for (i = 1; i < n_c+1; i++) {
//      v[i+j*(n_c+2)] = v_c[i+j*(n_c+2)];
//    }
//  }

  // set u manually to zero
  for (j = 0; j < n_c+2; j++) {
    for (i = 0; i < n_c+2; i++) {
      u_c[i+j*(n_c+2)] = 0.0;
    }
  }

  // write values in coarse grid
  Grid_Set(H, u_c, v_c, n_c);
}

//--------------------------- prolongation ---------------------------------//
void Prolongation(Grid H, Grid G) {
  double *v_f;                              // values at              //
  double *u_c, *u_f, *u_tmp;                              // current & finer grid   //
  unsigned int n_c, n_f;                          // number of grid points  //
  unsigned int i, j;

  // get and allocate values on current and finer grid
  n_c = Grid_Get_n(H);
  u_c = Grid_Get_u(H);
  n_f = Grid_Get_n(G);
  u_f = Grid_Get_u(G);
  v_f = Grid_Get_v(G);

  u_tmp = (double *) calloc((n_f+2)*(n_f+2), sizeof(double));

  // pad fine temporary array with zeros at non-defined places
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      if ((j%2 == 0) && (i%2 == 0)) {
        u_f[i+j*(n_f+2)] = u_c[i/2+j/2*(n_c+2)];
      }else {
        u_f[i+j*(n_f+2)] = 0.0;
      }
    }
  }

//  // set everything to zero to avoid side-effects
//  for (j = 0; j < n_c+2; j++) {
//    for (i = 0; i < n_c+2; i++) {
//      u_c[i+j*(n+2)] = 0.0;
//    }
//  }

//  Grid_Set(G,u_f,v_f,n_f);

  // loop over all elements in finer array
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      u_tmp[i+j*(n_f+2)] = u_f[i+j*(n_f+2)] + 0.5*(u_f[(i+1)+j*(n_f+2)]
          + u_f[(i-1)+j*(n_f+2)] + u_f[i+(j+1)*(n_f+2)]
          + u_f[i+(j-1)*(n_f+2)]) + 0.25*(u_f[(i+1)+(j+1)*(n_f+2)]
          + u_f[(i+1)+(j-1)*(n_f+2)] + u_f[(i-1)+(j+1)*(n_f+2)]
          + u_f[(i-1)+(j-1)*(n_f+2)]);
    }
  }
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      u_f[i+j*(n_f+2)] = u_tmp[i+j*(n_f+2)];
    }
  }

  // write prolongated values in grid
  Grid_Set(G,u_f,v_f,n_f);
  free(u_tmp);
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
      v[i+j*(n+2)] = v[i+j*(n+2)] + alpha*(4.*u[i+j*(n+2)] -
                       u[(i+1)+j*(n+2)] - u[(i-1)+j*(n+2)] -
                       u[i+(j+1)*(n+2)] - u[i+(j-1)*(n+2)]);
    }
  }

  Grid_Set_v(G,v);
}

//----------------------------- maximum norm -------------------------------//
double MaxNorm(double * u, unsigned int n){
  double umax;
  unsigned int i, j;

  umax = 0.0;
  for (j = 1; j < n+1; j++) {
    for (i = 1; i < n+1; i++) {
      umax = fmax(umax,fabs(u[i+j*(n+2)]));
    }
  }

  return umax;
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
