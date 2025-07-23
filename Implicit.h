//
// Sail 7/22/25
//

#ifndef IMPLICIT_H
#define IMPLICIT_H

void SemiImplicitSolve(const int nx, const int ndegr,      const int nvar,  const double dx,
    double* u, const double* Dmatrix,const double dt, double* dudt);

#endif // IMPLICIT_H
