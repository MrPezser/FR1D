//
// Created by Tsail on 7/26/2023.
//


#ifndef FR1D_SPATIALDISCRETIZATION_H
#define FR1D_SPATIALDISCRETIZATION_H

#include "indexing.h"
#include "EulerFlux.h"
#include <cmath>

//Swap this out to change equations. If you want to compare different PDEs, you've come to the wrong place buddy.
#define GAM (1.4)
#define MU (1e-3)
#define A (1.0)
#define FACEFLUX(a, b, c) AdvectionFaceFlux(a, b, c) //LeerFlux(GAM, a, b, c)
#define FLUX(a,b)  AdvectionFlux(a, b)  //EulerFlux(GAM, a, b)
//
void CalcDudt(int nx, int ndegr, int nvar, double dx, double* u, const double* Dmatrix, const double* Dradau, double* dudt );

#endif //FR1D_SPATIALDISCRETIZATION_H
