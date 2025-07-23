#include "Implicit.h"
#include "indexing.h"
#include "SpatialDiscretization.h"
#include "LUtools.h"
#include <cstdio>

void BuildJacobian(
    const double dt, const int ndegr, const int nvar, const double* u,
    const double dx, const double* Dmatrix, double** D
                    ) {
    //build semi implicit jacobian (I/dt - 0.5 dF^D)
    // ONLY TREATING DISCONTINUOUS FLUX IMPLICITLY
    // ONLY FOR SCALAR ADVECTION EQUATION

    double dti = 1.0/dt;
    size_t ld = ndegr*nvar; ///length d

    for (int i=0; i<ld; i++){
        for (int j=0; j<ld; j++){
            D[i][j] = 0.0;
        }
        D[i][i] = dti;
    }

    // Calculate discontinuous flux jacobian
    // - 0.5  d f_i / d u_j
    // conversion from ref space to dimensional included in here

    for (int i=0; i<ndegr; i++){
        for (int j=0; j<ndegr; j++){
             D[i][j] += (1.0 / dx) * A * Dmatrix[iup(i, j, ndegr)];
         }
    }
}

void SemiImplicitSolve(
    const int nx, const int ndegr,      const int nvar,  const double dx,
    double* u,    const double* Dmatrix,const double dt, double* dudt
                       ) {
    // Routine to calculate the element-local jacobian and solve the
    // semi implicit system on each element
   
   if (nvar != 1) {printf("Wrong options selected, many such cases :(\n"); exit(0);}  

    auto D = (double**)malloc((nvar*ndegr) * sizeof(double*));
    for (int id=0; id<nvar*ndegr; id++) {
        D[id] = (double*)malloc((nvar*ndegr) * sizeof(double));
    }

    double delta_u[ndegr*nvar];

    for (int ielem=0; ielem<nx; ielem++){
        // Build the jacobian matrix.. interestingly is constant for advection eqn
        BuildJacobian(dt, ndegr, nvar, u, dx, Dmatrix, D);
        
        int N = ndegr*nvar;
        int P[N+1]{}; //permutation vector for pivoting

         //get the rhs block needed
        double *b = &(dudt[iu3(ielem, 0, 0, ndegr)]);

        //Evaluate the jacobian / Implicit matrix
        LUPDecompose(D, N, 1.0e-12, P);
        LUPSolve(D, P, b, N, delta_u);

        for (int idegr=0; idegr<ndegr; idegr++) {
            u[iu3(ielem,idegr,0,ndegr)] += delta_u[idegr];
        }
    }

}

