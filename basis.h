//
// Created by Tsail on 6/16/2023.
//

#ifndef FR1D_BASIS_H
#define FR1D_BASIS_H

void GenerateLobattoPoints(int numPoints, double* points);
void GenerateLagrangeDMatrix( int ndegr, const double* loPoints, double* Dmatrix);
void GenerateRadauDerivatives(int ndegr, const double *x, double* Dradau);

#endif //FR1D_BASIS_H
