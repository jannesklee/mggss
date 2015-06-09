#ifndef Header_Routines

#define Header_Routines

//----------------------------- functions ----------------------------------//
void MG_Method(Grid, double *);
void Gauss_Seidel(Grid);

void Restriction(Grid, Grid);

void Prolongation(Grid, Grid);

void AddEval(double, Grid);

double MaxNorm(double *, unsigned int);

void Output(Grid);
//--------------------------------------------------------------------------//
#endif
