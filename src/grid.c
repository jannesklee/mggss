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

#include <stdio.h>                                // printf, scanf, NULL    //
#include <stdlib.h>                               // malloc, calloc, free   //
#include "grid.h"

//------------------------------ structs -----------------------------------//
struct Grid {
  double *u;                                      // solution pointer       //
  double *v;                                      // right hand side point. //
  unsigned int  n;                                // number of grid points  //
  double  h;                                      // stepwidth              //
};
//--------------------------------------------------------------------------//

//--------------------------- getters, setters -----------------------------//
void *  Grid_Get_u (Grid G) {return G->u;}        // get solution           //
void *  Grid_Get_v (Grid G) {return G->v;}        // get right-hand-side    //
unsigned int Grid_Get_n (Grid G) {return G->n;}   // get number grid points //
double  Grid_Get_h (Grid G) {return G->h;}        // get stepwidth          //

void    Grid_Set   (Grid G, double * u, double * v, unsigned int n) {
  G->u = u;                                       // set solution           //
  G->v = v;                                       // set right-hand-side    //
  G->n = n;                                       // set number grid points //
  G->h = 1./(n+1.);                               // set stepwidth          //
}
void    Grid_Set_u (Grid G, double * u) {G->u = u;}
void    Grid_Set_v (Grid G, double * v) {G->v = v;}
void    Grid_Set_n (Grid G, unsigned int n) {
  G->n = n;                                       // set n and implicitely  //
  G->h = 1./(n+1.);                               // also h                 //
}

//--------------------------------------------------------------------------//

//------------------------ memory (de-)allocation --------------------------//
Grid Grid_Create (void) {
  return calloc(1,sizeof(struct Grid));
}

void Grid_Destroy (Grid * G) {
  free(*G);
  *G = NULL;
}
//--------------------------------------------------------------------------//
