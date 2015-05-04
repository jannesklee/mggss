#include <stdio.h>                                // printf, scanf, NULL    //
#include <stdlib.h>                               // malloc, calloc, free   //
#include "grid.h"

//------------------------------ structs -----------------------------------//
struct Grid {
  double *u;                                      // solution pointer       //
  double *v;                                      // right hand side point. //
  double  h;                                      // stepwidth              //
};
//--------------------------------------------------------------------------//

//--------------------------- getters, setters -----------------------------//
double  Grid_Get_h (Grid G) {return G->h;}        // get stepwidth          //
void *  Grid_Get_u (Grid G) {return G->u;}        // get solution           //
void *  Grid_Get_v (Grid G) {return G->v;}        // get right-hand-side    //

void    Grid_Set   (Grid G, double * u, double * v, unsigned int N) {
  G->u = u;                                       // set solution           //
  G->v = v;                                       // set right-hand-side    //
  G->h = 1./(N+1.);                               // set stepwidth          //
}
void    Grid_Set_u (Grid G, double * u) {G->u = u;}
void    Grid_Set_v (Grid G, double * v) {G->v = v;}

//--------------------------------------------------------------------------//

//-------------------------- memory allocation -----------------------------//
Grid Grid_Create (void) {
  return calloc(1,sizeof(struct Grid));
}

void Grid_Destroy (Grid * G) {
  *G = NULL;
}
//--------------------------------------------------------------------------//
