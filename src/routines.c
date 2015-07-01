/*****************************************************************************
*                                                                            *
* MGGSS - Multigrid-Gauss-Seidel-Solver                                      *
* Copyright (C) <2015>                                                       *
* Jannes Klee <jklee@astrophysik.uni-kiel.de>                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*                                                                            *
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "math.h"
#include "grid.h"
#include "routines.h"

//------------------------- multigrid-method -------------------------------//
void MG_Method(Grid G, double *eps, unsigned int k, unsigned int gamma) {
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

  // restriction - writes the residue in v and sets u to zero
  Restriction(G, H);
  *eps = MaxNorm(Grid_Get_v(G),Grid_Get_n(G));

  // smooth coarser error
  Gauss_Seidel(H);

  // recursive call of the multigrid method
  if (k < gamma) {
    MG_Method(H, eps, k + 1, gamma);
  }

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

  // calculate residual r and write it in v
  AddEval(-1.0, G);

  // loop over even, even combination in order to define coarse grid
  for (j = 2; j < n+1; j = j+2) {
    for (i = 2; i < n+1; i = i+2) {
      v_c[i/2+j/2*(n_c+2)] = 0.25*(v[i+j*(n+2)] + 0.5*(v[(i+1)+j*(n+2)] +
          v[(i-1)+j*(n+2)] + v[i+(j+1)*(n+2)] + v[i+(j-1)*(n+2)]) +
            0.5*(v[(i+1)+(j+1)*(n+2)] + v[(i-1)+(j-1)*(n+2)]));
    }
  }

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
  double *v_f;                                    // values at              //
  double *u_c, *u_f;                      // current & finer grid   //
  unsigned int n_c, n_f;                          // number of grid points  //
  unsigned int i, j;

  // get and allocate values on current and finer grid
  n_c = Grid_Get_n(H);
  u_c = Grid_Get_u(H);
  n_f = Grid_Get_n(G);
  u_f = Grid_Get_u(G);
  v_f = Grid_Get_v(G);

  // loop over all elements in finer array
  for (j = 1; j < n_f+1; j++) {
    for (i = 1; i < n_f+1; i++) {
      if ((j%2 == 0) && (i%2 == 0)) {
        u_f[i+j*(n_f+2)] = u_c[i/2+j/2*(n_c+2)];
      }
      else if ((j%2 == 0) && (i%2 != 0)) {
        u_f[i+j*(n_f+2)] = 0.5*(u_c[(i+1)/2+j/2*(n_c+2)] +
                                u_c[(i-1)/2+j/2*(n_c+2)]);
      }
      else if ((j%2 != 0) && (i%2 == 0)) {
        u_f[i+j*(n_f+2)] = 0.5*(u_c[i/2+(j+1)/2*(n_c+2)] +
                                u_c[i/2+(j-1)/2*(n_c+2)]);
      }
      else if ((j%2 != 0) && (i%2 != 0)) {
        u_f[i+j*(n_f+2)] = 0.5*(u_c[(i+1)/2+(j+1)/2*(n_c+2)] +
                                u_c[(i-1)/2+(j-1)/2*(n_c+2)]);
      }
    }
  }

  // write prolongated values in grid
  Grid_Set(G,u_f,v_f,n_f);
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
