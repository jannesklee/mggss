#ifndef Header_Routines

#define Header_Routines

//----------------------------- functions ----------------------------------//
void MG_Method(Grid, double *);
void Gauss_Seidel(Grid, double *);

void Restriction(Grid);

void Prolongation(Grid);

void AddEval(double, Grid);

void Output(Grid);
//--------------------------------------------------------------------------//
#endif
