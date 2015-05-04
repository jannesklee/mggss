#ifndef Header_Grid

#define Header_Grid
//------------------------------ structs -----------------------------------//
typedef struct Grid * Grid;
//--------------------------------------------------------------------------//

//--------------------------- getters, setters -----------------------------//
double  Grid_Get_h (Grid);
void *  Grid_Get_u (Grid);
void *  Grid_Get_v (Grid);

void Grid_Set (Grid, double *, double *, unsigned int);
//--------------------------------------------------------------------------//

//-------------------------- memory allocation -----------------------------//
Grid Grid_Create (void);
void Grid_Destroy (Grid *);
//--------------------------------------------------------------------------//

#endif
