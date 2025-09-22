//
// Created by Tsail on 6/24/2023.
//

#ifndef FR1D_EULERFLUX_H
#define FR1D_EULERFLUX_H

#include "indexing.h"

void EulerFlux(double gam, const double *u, double* flux);
void RoeFDS(double gam, const double* uL, const double *uR, double* roeFlux);
void LeerFlux(double gam, const double* uL, const double* uR, double *flux);
void LDFSS(double gam, const double* uL, const double* uR, double *flux);
void getPrimativesPN(double gam, const double *unkel, double *rho, double *v, double *p, double *c, double *M);


#endif //FR1D_EULERFLUX_H
