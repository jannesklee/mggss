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
#include <math.h>
#include "grid.h"                                 // grid struct            //
#include "routines.h"                             // gauss_seidel routine   //

void initialize(double *, double *, unsigned int, unsigned int);

int main (void) {
  //--------------------------- declarations -------------------------------//
  Grid          G;                                // grid                   //
  double        *u, *v;                           // solution u and rhs v   //
  unsigned int  n = 199;                          // res. has to be odd     //
  double        eps0 = 1e-10;                     // error-limit            //
  double        eps  = 0.0;                       // error                  //
  double        eps_old = 0.0;
  unsigned int  k = 1;                            // 1: const. dens. circle //
                                                  // 2: exponential distr.  //
  unsigned int  gamma = 6;                        // numb. of mesh refinem. //

  //---------------------------- allocations -------------------------------//
  u = (double *) calloc((n+2)*(n+2), sizeof(double));
  v = (double *) calloc((n+2)*(n+2), sizeof(double));

  //--------------------------- initialization -----------------------------//
  initialize(u, v, n, k);

  // create grid with initial values for poisson equation
  G = Grid_Create();
  Grid_Set(G, u, v, n);

  //----------------------------- main loop --------------------------------//
  do {
    eps_old = eps;
    MG_Method(G, &eps, 0, gamma);
    printf("Error: %.2e \t Error_Fractional: %.2e \n", eps, eps/eps_old);
  } while(eps > eps0);
  Output(G);

  //--------------------------- deallocations ------------------------------//
  free(u);
  free(v);
  Grid_Destroy(&G);
  return 0;
}

void initialize(double *u, double *v, unsigned int n, unsigned int k) {
  double        r;                                // radius for init.       //
  unsigned int  i,j;

  if (k== 1) { // enhanced density on circle
    // initial values & boundaries
    for (i = 0; i < n+2; i++) {
      for (j = 0; j < n+2; j++) {
        r = sqrt((1./((n-1.)*(n-1.))*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2))));
        if (r < 0.1) {
          v[i+j*(n+2)] = 100.;
        } else {
          v[i+j*(n+2)] = 0.0;
        }
        u[i+j*(n+2)] = 0.0;
      }
    }
  }
  else if (k == 2) { // gaussian distribution
    for (i = 0; i < n+2; i++) {
      for (j = 0; j < n+2; j++) {
        r = sqrt((1./((n-1.)*(n-1.))*((i-n/2)*(i-n/2)+(j-n/2)*(j-n/2))));
        v[i+j*(n+2)] = exp(-1000*r*r);
      }
    }
  }
  else {
    // do nothing
  }
}
