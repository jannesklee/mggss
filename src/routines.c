#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "grid.h"
#include "routines.h"

//------------------------- multigrid-method -------------------------------//
void MG_Method(Grid G, double *eps) {
  double *u, *u_tmp;
  unsigned int n;
  unsigned int i, j;

  n = Grid_Get_n(G);
  u = Grid_Get_u(G);

  u_tmp = (double *) malloc((n+2)*(n+2)*sizeof(double));

  // pre-smoothing
  Gauss_Seidel(G, eps); // gauss-seidel with residual


  // save u in u_tmp
  u = Grid_Get_u(G);
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      u_tmp[i+j*(n+2)] = u[i+j*(n+2)];
    }
  }

  // calculate residual r and write it in v
  AddEval(1.0, G);

  Restriction(G);

  // set u (now the error) to zero
  n = Grid_Get_n(G);
  u = Grid_Get_u(G);
  memset(u, 0, n*sizeof(u[0]));
  Grid_Set_u(G, u);

  // recalculate on coarser grid (TODO: here eventuelly exact solver?)
  Gauss_Seidel(G, eps);

  // prolongate array
  Prolongation(G);

  // finally calculate solution
  u = Grid_Get_u(G);
  for (j = 0; j < n+2; j++) {
    for (i = 0; i < n+2; i++) {
      u_tmp[i+j*(n+2)] = u_tmp[i+j*(n+2)] + u[i+j*(n+2)];
    }
  }

  Grid_Set_u(G,u_tmp);

  // post-smoothing
  Gauss_Seidel(G, eps);
}

//---------------------------- gauss-seidel --------------------------------//
void Gauss_Seidel (Grid G, double *eps) {
  double *u, *v;                                  // solution, rhs          //
  double h;                                       // stepwidth              //
  unsigned int n;                                 // number of grid points  //
  double tmp;
  unsigned int i, j;

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

  Grid_Set_u(G,u);
}

//---------------------------- restriction ---------------------------------//
void Restriction (Grid G) {
  double *r, *r_c;                                // actual & coarser grid  //
  double *u_c;                                    // coarser grid           //
  unsigned int n, n_c;                            // number of grid-points  //
  unsigned int i, j;

  // get grid sizes of actual and calculate coarser grid
  n = Grid_Get_n(G);
  n_c = n/2;

  // get the residual (equals v)
  r = Grid_Get_v(G);

  // allocate residual for coarser grid
  u_c = (double *) malloc((n_c+2)*(n_c+2)*sizeof(double));
  r_c = (double *) malloc((n_c+2)*(n_c+2)*sizeof(double));

  // loop over even, even combinations in order to define coarse grid
  for (j = 2; j < n+1; j = j+2) {
    for (i = 2; i < n+1; i = i+2) {
      r_c[i/2+j/2*(n_c+2)] = 0.25*(r[i+j*(n+2)] + 0.5*(r[(i+1)+j*(n+2)] +
          r[(i-1)+j*(n+2)] + r[i+(j+1)*(n+2)] + r[i+(j-1)*(n+2)]) +
          0.25*(r[(i+1)+(j+1)*(n+2)] + r[(i-1)+(j+1)*(n+2)] +
          r[(i+1)+(j-1)*(n+2)] + r[(i-1)+(j-1)*(n+2)]));
    }
  }

  // overwrite grid
  Grid_Set(G, u_c, r_c, n_c);
}

//--------------------------- prolongation ---------------------------------//
void Prolongation(Grid G) {
  double *r, *r_f;                                // residues both grids    //
  double *tmp_f, *u_f;                            // finer grid             //
  unsigned int n, n_f;                            // number of points       //
  unsigned int i, j;

  // actual (coarser) grid
  r = Grid_Get_v(G);
  n = Grid_Get_n(G);

  // finer grid
  n_f = n*2;
  r_f = (double *) malloc((n_f+2)*(n_f+2)*sizeof(double));
  tmp_f = (double *) malloc((n_f+2)*(n_f+2)*sizeof(double));
  u_f = (double *) malloc((n_f+2)*(n_f+2)*sizeof(double));

  // pad fine array with zeros at non-defined places
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      if ((j%2 == 0) && (i%2 == 0)) {
        tmp_f[i+j*(n_f+2)] = r[i/2+j/2*(n+2)];
      }else {
        tmp_f[i+j*(n_f+2)] = 0;
      }
    }
  }

  Grid_Set(G,u_f,tmp_f,n_f);

  // loop over all elements in finer array
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      r_f[i+j*(n_f+2)] = tmp_f[i+j*(n_f+2)] + 0.5*(tmp_f[(i+1)+j*(n_f+2)]
          + tmp_f[(i-1)+j*(n_f+2)] + tmp_f[i+(j+1)*(n_f+2)]
          + tmp_f[i+(j-1)*(n_f+2)]) + 0.25*(tmp_f[(i+1)+(j+1)*(n_f+2)]
          + tmp_f[(i+1)+(j-1)*(n_f+2)] + tmp_f[(i-1)+(j+1)*(n_f+2)]
          + tmp_f[(i-1)+(j-1)*(n_f+2)]);
    }
  }

  // overwrite grid
  Grid_Set(G,u_f,r_f,n_f);
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
