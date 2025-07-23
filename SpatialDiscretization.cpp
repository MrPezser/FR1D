//
// Created by Tsail on 7/26/2023.
//
#include "SpatialDiscretization.h"

#include <cmath>
#include <stdexcept>

#include "indexing.h"
#include "EulerFlux.h"

double AdvectionFaceFlux(const double a, const double uL, const double uR) {
    /// Finds the upwind flux according to the 1D advection equation
    ///a = wave speed | uL & uR are the left and right state variables
    if (a > 0){
        return a * uL;
    } else {
        return a * uR;
    }
}

double AdvectionFlux(const double a, const double u) {
    /// Finds the flux according to the 1D advection equation
    ///a = wave speed | u is state variable

    return a * u;
}


void FluxFaceCorrection(int nx, int ndegr, int nvar, double* u, const double* Dradau, double* fcorr_xi, double* dudt){
    double fL[nvar], fR[nvar], common_flux[nvar];
    double gam = 1.4;

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //extrapolation boundary condition
            //iep1 = iface;
            //periodic bc
            iep1 = 0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }

        //Calculate the common flux at the face
        if (nvar == 3){
            RoeFDS(gam,&u[iu3(ielem, ndegr-1, 0, ndegr)], &u[iu3(iep1, 0, 0, ndegr)], common_flux);
        }else{
            common_flux[0] = AdvectionFaceFlux(A, u[iu3(ielem, ndegr-1, 0, ndegr)], u[iu3(iep1, 0, 0, ndegr)]);
        }

        //Calculate the element's local fluxes at the face
        if (nvar == 3){
            EulerFlux(gam,&u[iu3(ielem , ndegr-1 , 0, ndegr)], &fL[0]);
            EulerFlux(gam,&u[iu3(iep1  , 0       , 0, ndegr)], &fR[0]);
        } else {
            fL[0] = AdvectionFlux(A, u[iu3(ielem, ndegr-1, 0, ndegr)]);
            fR[0] = AdvectionFlux(A, u[iu3(iep1 , 0      , 0, ndegr)]);
        }

        for (int inode=0; inode<ndegr; inode++){   //loop over solution points
            for (int kvar = 0; kvar<nvar; kvar++) {   //loop over variables/equations (used when Euler Equation fluxes enabled)
                //Element left of face flux slope
                fcorr_xi[iu3(ielem,ndegr-1-inode, kvar, ndegr)] -= (common_flux[kvar] - fL[kvar]) * (Dradau[inode]);
                //Element right of face flux slope
                fcorr_xi[iu3(iep1, inode        , kvar, ndegr)] += (common_flux[kvar] - fR[kvar]) * (Dradau[inode]);
            }
        }
    }
}

void DiscontinuousFlux(const int nx, const int np, const int ndegr, const int nvar, const double* u, const double* Dmatrix, double* fcorr_xi) {
    double flux_node[np];
    double fxi[nvar];
    double gam = 1.4;

    for (int inode=0; inode<ndegr; inode++){
        for (int kvar=0; kvar<nvar;kvar++){
            flux_node[iup(inode, kvar, nvar)] = 0.0;
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++) {

        for (int inode = 0; inode < ndegr; inode++) {

            if (nvar == 3) {
                double flux[nvar];
                EulerFlux(gam,&u[iu3(ielem, inode, 0, ndegr)], &flux[0]);
                //
                flux_node[iup(inode, 0, nvar)] = flux[0];
                flux_node[iup(inode, 1, nvar)] = flux[1];
                flux_node[iup(inode, 2, nvar)] = flux[2];
            } else {
                flux_node[inode] = AdvectionFlux(A, u[iu3(ielem, inode , 0, ndegr)]);
            }
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int inode = 0; inode < ndegr; inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                fxi[kvar] = 0.0;
            }

            //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
            for (int jnode = 0; jnode < ndegr; jnode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {
                    fxi[kvar] += flux_node[iup(jnode, kvar, nvar)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                fcorr_xi[iu3(ielem, inode, kvar, ndegr)] += fxi[kvar];
            }
        }
    }

}

void CalcDudt(const int nx, const int ndegr, const int nvar, const double dx, double* u, const double* Dmatrix, const double* Dradau, double* dudt ){
    ///Calculates the solution update given a function to find flux
    int nu = nx * ndegr * nvar;
    int np = ndegr * nvar;
    double fcorr_xi[nu], ucorr_xi[nu], ucorr_x[nu], ucx_xi[nu];


    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
        ucorr_xi[i] = 0.0;
        ucx_xi[i] = 0.0;
    }

    //Flux Reconstruction (P_ndegr)
    FluxFaceCorrection(nx, ndegr, nvar, u, Dradau, fcorr_xi, dudt);
    DiscontinuousFlux(nx, np, ndegr, nvar, u, Dmatrix, fcorr_xi);

    for (int ielem=0; ielem<nx; ielem++) {
        for (int inode=0; inode<ndegr ;inode++) {
            for (int kvar = 0; kvar<nvar; kvar++) {
                //convert from the reference element to the real element and negate to put make it dudt (= -f_x)
                dudt[iu3(ielem, inode, kvar, ndegr)] += -(2/ dx) * (fcorr_xi[iu3(ielem,  inode, kvar, ndegr)] - MU*ucx_xi[iu3(ielem,  inode, kvar, ndegr)]);

               /* //fix end points for shock-tube like problem
               if (ielem == 0 || ielem == nx-1){
                   dudt[iu3(ielem, inode, kvar, ndegr)] = 0;
               }*/

                if (__isnan(dudt[iu3(ielem, inode, kvar, ndegr)])){
                    throw std::overflow_error("dudt NAN\n");
                }
            }
        }
    }

}
