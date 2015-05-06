#ifndef Header_Grid

#define Header_Grid
//------------------------------ structs -----------------------------------//
typedef struct Grid * Grid;
//--------------------------------------------------------------------------//

//--------------------------- getters, setters -----------------------------//
double        Grid_Get_h (Grid);
void *        Grid_Get_u (Grid);
void *        Grid_Get_v (Grid);
unsigned int  Grid_Get_n (Grid);

void          Grid_Set (Grid, double *, double *, unsigned int);
void          Grid_Set_u (Grid, double *);
void          Grid_Set_v (Grid, double *);
//--------------------------------------------------------------------------//

//-------------------------- memory allocation -----------------------------//
Grid Grid_Create (void);
void Grid_Destroy (Grid *);
//--------------------------------------------------------------------------//

#endif
