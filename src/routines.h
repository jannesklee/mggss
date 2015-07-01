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

#ifndef Header_Routines

#define Header_Routines

//----------------------------- functions ----------------------------------//
void MG_Method(Grid, double *, unsigned int, unsigned int);

void Gauss_Seidel(Grid);

void Restriction(Grid, Grid);

void Prolongation(Grid, Grid);

void AddEval(double, Grid);

double MaxNorm(double *, unsigned int);

void Output(Grid);
//--------------------------------------------------------------------------//
#endif
