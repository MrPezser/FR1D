#include "Implicit.h"
#include "indexing.h"
#include "SpatialDiscretization.h"
#include "LUtools.h"
#include "EulerFlux.h"

#include <cstdio>
#include <cstdlib>

void BuildJacobian(
    const double dt, const int ndegr, const int nvar, const double* u,
    const double dx, const double* Dmatrix, double** D
                    ) {
    //build semi implicit jacobian (I/dt - 0.5 dF^D)
    // ONLY TREATING DISCONTINUOUS FLUX IMPLICITLY
    // ONLY FOR SCALAR ADVECTION EQUATION
    double gam = 1.4;
    double dti = 1.0/dt;
    size_t ld = ndegr*nvar; ///length d

    for (int i=0; i<ld; i++){
        for (int j=0; j<ld; j++){
            D[i][j] = 0.0;
        }
        D[i][i] = dti;
    }

    if (nvar ==3) {
        //
        for (int idegr=0; idegr<ndegr; idegr++){
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            //for each pair of degree dependance, calculate the flux jacobian and multiply by lagrange slope        
            //
            // Calculate u at x_j
            const double *uj = &(u[iup(jdegr,0,nvar)]); 
            double rho, v, p, c, M;
            getPrimativesPN(gam, uj, &rho, &v, &p, &c, &M);
            //
            // calculate dfj/duj
            double df1du1, df2du1, df3du1, df1du2, df2du2, df3du2, df1du3, df2du3, df3du3;
            double flux[nvar], fluxp1[nvar], fluxp2[nvar], fluxp3[nvar];
            // Numerical difference to get flux derivative
            EulerFlux(gam, uj, flux);
            //
            double utest[nvar];
            utest[0] = uj[0]*(1.0+1.0e-8);
            utest[1] = uj[1];
            utest[2] = uj[2];
            EulerFlux(gam, utest, fluxp1);
            //
            utest[0] = uj[0];
            utest[1] = uj[1]*(1.0+1.0e-8) + 1.0e-16;
            utest[2] = uj[2];
            EulerFlux(gam, utest, fluxp2);
            //
            utest[0] = uj[0];
            utest[1] = uj[1];
            utest[2] = uj[2]*(1.0+1.0e-8);
            EulerFlux(gam, utest, fluxp3);
            //
            df1du1 = (fluxp1[0]-flux[0])/(uj[0]*1.0e-8);
            df1du2 = (fluxp2[0]-flux[0])/(uj[1]*1.0e-8 + 1.0e-16);
            df1du3 = (fluxp3[0]-flux[0])/(uj[2]*1.0e-8);
            //
            df2du1 = (fluxp1[1]-flux[1])/(uj[0]*1.0e-8);
            df2du2 = (fluxp2[1]-flux[1])/(uj[1]*1.0e-8 + 1.0e-16);
            df2du3 = (fluxp3[1]-flux[1])/(uj[2]*1.0e-8);
            //
            df3du1 = (fluxp1[2]-flux[2])/(uj[0]*1.0e-8);
            df3du2 = (fluxp2[2]-flux[2])/(uj[1]*1.0e-8 + 1.0e-16);
            df3du3 = (fluxp3[2]-flux[2])/(uj[2]*1.0e-8);
            /*
            df1du1 = 0.0; // v
            df1du2 = 1.0;
            df1du3 = 0.0;
            //
            df2du1 = p/rho; // + v*v
            df2du2 = v; //(3.0 - gam) * v;
            df2du3 = gam-1.0;
            //
            df3du1 = v * (p) / rho; // + v*uj[2]/rho
            df3du2 = (uj[2] + p) / rho;//(gam*uj[2]/rho) + ((gam-1.0)/2.0)*v*v;
            df3du3 = gam*v;
            */
            //
            // Use FR to transform to dfi/duj
            int ii = iup(idegr,0,nvar);
            int jj = iup(jdegr,0,nvar);
            double lagrange_slope = Dmatrix[(idegr,jdegr,ndegr)];
            double coeff = ALPHA * lagrange_slope / dx;
            //
            D[ii+0][jj+0] += coeff * df1du1;
            D[ii+0][jj+1] += coeff * df1du2;
            D[ii+0][jj+2] += coeff * df1du3;
            //
            D[ii+1][jj+0] += coeff * df2du1;
            D[ii+1][jj+1] += coeff * df2du2;
            D[ii+1][jj+2] += coeff * df2du3;
            //
            D[ii+2][jj+0] += coeff * df3du1;
            D[ii+2][jj+1] += coeff * df3du2;
            D[ii+2][jj+2] += coeff * df3du3;
        }
        }
    } else {
        // Calculate discontinuous flux jacobian
        // - 0.5  d f_i / d u_j
        // conversion from ref space to dimensional included in here
 
        for (int i=0; i<ndegr; i++){
            for (int j=0; j<ndegr; j++){
                 D[i][j] += ALPHA * (1.0 / dx) * A * Dmatrix[iup(i, j, ndegr)];
            }
        }
    }
}

void SemiImplicitSolve(
    const int nx, const int ndegr,      const int nvar,  const double dx,
    double* u,    const double* Dmatrix,const double dt, double* dudt
                       ) {
    // Routine to calculate the element-local jacobian and solve the
    // semi implicit system on each element
   
    auto D = (double**)malloc((nvar*ndegr) * sizeof(double*));
    for (int id=0; id<nvar*ndegr; id++) {
        D[id] = (double*)malloc((nvar*ndegr) * sizeof(double));
    }

    double delta_u[ndegr*nvar];

    for (int ielem=0; ielem<nx; ielem++){
        // Build the jacobian matrix.. interestingly is constant for advection eqn
        double *ui = &(u[iu3(ielem, 0, 0, ndegr)]);
        BuildJacobian(dt, ndegr, nvar, ui, dx, Dmatrix, D);
        
        int N = ndegr*nvar;
        int P[N+1]{}; //permutation vector for pivoting

         //get the rhs block needed
        double *b = &(dudt[iu3(ielem, 0, 0, ndegr)]);

        //Evaluate the jacobian / Implicit matrix
        LUPDecompose(D, N, 1.0e-12, P);
        LUPSolve(D, P, b, N, delta_u);

        for (int idegr=0; idegr<ndegr; idegr++) {
            for (int kvar=0; kvar<nvar; kvar++) {
                u[iu3(ielem,idegr,kvar,ndegr)] += delta_u[iup(idegr, kvar, nvar)];
            }
        }
    }

}

