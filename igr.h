//
// Header file for igr.cpp
//
#ifndef IGR_H
#define IGR_H


void CalculateIGRSigma(int nelem, int ndegr, int nvar, double dx, double alpha, const double* Dmatrix, const double *Dradau, double* u, double* du, double* igr_sigma);

#endif //IGR_H
