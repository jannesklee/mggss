#ifndef Header_Routines

#define Header_Routines

//----------------------------- functions ----------------------------------//
void Gauss_Seidel(Grid, double *);

void Restriction(Grid);

void Prolongation(Grid);

double * AddEval(double, Grid);
//--------------------------------------------------------------------------//
#endif
