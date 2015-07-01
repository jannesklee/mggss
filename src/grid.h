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
void          Grid_Set_n (Grid, unsigned int);
//--------------------------------------------------------------------------//

//-------------------------- memory allocation -----------------------------//
Grid Grid_Create (void);
void Grid_Destroy (Grid *);
//--------------------------------------------------------------------------//

#endif
