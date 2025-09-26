//
// Created by Tsail on 7/26/2023.
//
#include "SpatialDiscretization.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

#include "indexing.h"
#include "EulerFlux.h"
#include "igr.h"

#define RHOLIMIT (1e-2)
#define ELIMIT (1e-2)

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

void StateDerivatives(const double dx, const int nx, const int np, const int ndegr, const int nvar, const double* u, const double* Dmatrix, double* du_dx) {
    double uxi[nvar];
    double gam = 1.4;

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++) {

        //compute discontinuous solution slope (in reference element of width 2)
        for (int inode = 0; inode < ndegr; inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                uxi[kvar] = 0.0;
            }
            //Use the Lagrange polynomial derivative matrix to perform the MatVec
            for (int jnode = 0; jnode < ndegr; jnode++) {
                for (int kvar = 0; kvar < nvar; kvar++) {
                    uxi[kvar] += u[iu3(ielem, jnode , kvar, ndegr)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                du_dx[iu3(ielem, inode, kvar, ndegr)] += 2.0 * uxi[kvar] / dx;
            }
        }
    }

}

void FluxFaceCorrection(int nx, int ndegr, int nvar, double* u, const double* Dradau, double* fcorr_xi, double* igr_sigma){
    double fL[nvar], fR[nvar], common_flux[nvar];
    double gam = 1.4;

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //extrapolation boundary condition
            iep1 = iface;
            //periodic bc
            //iep1 = 0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }
        //
        double sigmaface[2];
        sigmaface[0] = igr_sigma[iup(ielem,ndegr-1,ndegr)];
        sigmaface[1] = igr_sigma[iup(iep1 ,0      ,ndegr)];


        //Calculate the common flux at the face
        if (nvar == 3){
            //LDFSS(gam, &u[iu3(ielem, ndegr-1, 0, ndegr)], &u[iu3(iep1, 0, 0, ndegr)], common_flux);
            LeerFlux(gam, &u[iu3(ielem, ndegr-1, 0, ndegr)], &u[iu3(iep1, 0, 0, ndegr)], common_flux, sigmaface);
        }else{
            common_flux[0] = AdvectionFaceFlux(A, u[iu3(ielem, ndegr-1, 0, ndegr)], u[iu3(iep1, 0, 0, ndegr)]);
        }

        //Calculate the element's local fluxes at the face
        if (nvar == 3){
            EulerFlux(gam,&u[iu3(ielem , ndegr-1 , 0, ndegr)], &fL[0], sigmaface[0]);
            EulerFlux(gam,&u[iu3(iep1  , 0       , 0, ndegr)], &fR[0], sigmaface[1]);
            //
            //LDFSS(gam, &u[iu3(ielem , ndegr-1 , 0, ndegr)], &u[iu3(ielem , ndegr-1 , 0, ndegr)], &fL[0]);
            //LDFSS(gam, &u[iu3(iep1  , 0       , 0, ndegr)], &u[iu3(iep1  , 0       , 0, ndegr)], &fR[0]);
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

void DiscontinuousFlux(const int nx, const int np, const int ndegr, const int nvar, const double* u, const double* Dmatrix, double* fcorr_xi, double* igr_sigma) {
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
                double sig = igr_sigma[iup(ielem,inode,ndegr)];
                EulerFlux(gam,&u[iu3(ielem, inode, 0, ndegr)], &flux[0], sig);
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

void CalcDudt(const int nx, const int ndegr, const int nvar, const double dx, double* u, const double* Dmatrix, const double* Dradau, double* igr_sigma, double* dudt ){
    ///Calculates the solution update given a function to find flux
    int nu = nx * ndegr * nvar;
    int np = ndegr * nvar;
    double fcorr_xi[nu], ucorr_xi[nu], ucorr_x[nu], u_x[nu];
    
    double alpha = ALPH * dx*dx;

    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
        ucorr_xi[i] = 0.0;
        u_x[i] = 0.0;
    }

    //State Reconstruction
    StateDerivatives(dx, nx, np, ndegr, nvar, u, Dmatrix, u_x);
    CalculateIGRSigma(nx, ndegr, nvar, dx, alpha, Dmatrix, u, u_x, igr_sigma);

    //Flux Reconstruction (P_ndegr) 
    FluxFaceCorrection(nx, ndegr, nvar, u, Dradau, fcorr_xi, igr_sigma);
    DiscontinuousFlux(nx, np, ndegr, nvar, u, Dmatrix, fcorr_xi, igr_sigma);

    for (int ielem=0; ielem<nx; ielem++) {
        for (int inode=0; inode<ndegr ;inode++) {
            for (int kvar = 0; kvar<nvar; kvar++) {
                //convert from the reference element to the real element and negate to put make it dudt (= -f_x)
                dudt[iu3(ielem, inode, kvar, ndegr)] += -(2.0/ dx) * (fcorr_xi[iu3(ielem,  inode, kvar, ndegr)]);

               /* //fix end points for shock-tube like problem
               if (ielem == 0 || ielem == nx-1){
                   dudt[iu3(ielem, inode, kvar, ndegr)] = 0;
               }*/

                if (__isnan(dudt[iu3(ielem, inode, kvar, ndegr)])){
                    printf("ielem:%d \tinode:%d \tkvar:%d\n", ielem, inode, kvar);
                    throw std::overflow_error("dudt NAN\n");
                }
            }
        }
    }

}

bool IsBad(double* u) {
    if (u[0] < RHOLIMIT) {return true;}

    double p = (GAM - 1) * (u[2] - (0.5 * u[1] * u[1] / u[0]));
    if (u[2] < ELIMIT) {return true;}

    return false;
}

void StupidSmoothSolution(int nelem, int ndegr, int nvar, double* u) {
    // dumb averaging based smoother, doesn't work well at all
    // Smoothing strength
    int nsmooth = 100000.0; //guessing
    double ksmooth = 2.0 / nsmooth; // guessing
    //
    int npoin = nelem*ndegr;
    double unew[npoin*NVAR];
    // Apply diffusion kernel to all elements, ignoring first and last
    for (int it=0; it<nsmooth; it++) {
        for (int ip=0; ip<npoin; ip++){
            //
            double *ui = &u[ip*NVAR];
            double *uim = &u[(ip-1)*NVAR];
            double *uip = &u[(ip+1)*NVAR];
            //
            for (int k=0; k<NVAR; k++) {
                if (ip==0 || ip==npoin-1) {
                    unew[ip+k] = ui[k];
                } else {
                    unew[ip+k] = ((1.0-2.0*ksmooth)*ui[k] + ksmooth*(uim[k] + uip[k]));
                }
            }
        }
    memcpy(u,unew,npoin*NVAR*sizeof(double));
    }

}

void LimitSolution(int nelem, int ndegr, int nvar, double* u) {
    //return;
    // Limit by scaling
    // Setp 1: Find cell average (NOTE TO ADD WEIGHTING TO SOLPTS)
    // Step 2: Find the ratio: RR = (rmin-rlim)/(rbar-rmin)
    // Step 3: Set r = rbar + RR
    // Recalculate
    // Step 4:If Still bad, repeat with energy
    bool badcell;

    for (int ielem=0; ielem<nelem; ielem++) {
        double *ui,rhomin,emin;
        ui = &(u[iu3(ielem,0,0,ndegr)]);
        rhomin = FP_NAN;
        
        // Step -1 - find bad cells
        badcell = false;
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            badcell = IsBad(ui+jdegr*nvar) or badcell;
            //
            if (badcell && rhomin==FP_NAN){
                rhomin = (ui+jdegr*nvar)[0];
            } else if (badcell){
                rhomin = fmin(rhomin,(ui+jdegr*nvar)[0]);
            }
            //
        }
        if (not badcell) {continue;}
        // Step 1: Find cell average rho
        double rhoave = 0.0;
        // Currently doing the stupid thing of equal weighting
        // This assumption starts out okay but gets very bad for higher orders
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            rhoave += ui[jdegr*nvar] / (ndegr);
        }
        //rhoave = fmax(rhoave, 1.5*RHOLIMIT);
        ASSERT(rhoave >= RHOLIMIT, "average density below limit")
        // Step 2 - calculate RR
        double rr = (rhoave-RHOLIMIT) / (rhoave-rhomin);
        //
        // Step 3 - SQUASH ALL RESISTANCE
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            (ui+jdegr*nvar)[0] = fmax(RHOLIMIT, rhoave  + ((ui+jdegr*nvar)[0] - rhoave) * rr);
        }
        //
        // ==================================================  Recheck if cell bad
        badcell = false;
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            badcell = IsBad(ui+jdegr*nvar) or badcell;
            //
            if (badcell && rhomin==FP_NAN){
                emin = (ui+jdegr*nvar)[2];
            } else if (badcell){
                emin = fmin(rhomin,(ui+jdegr*nvar)[2]);
            }
            //
        }
        if (not badcell) {continue;}
        //
        double eave = 0.0;
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            eave += ui[jdegr*nvar + 2] / (ndegr);
        }
        //eave = fmax(eave, 1.5*ELIMIT);
        ASSERT(eave >= ELIMIT, "average energy below limit")
        double er = (eave-ELIMIT) / (eave-emin);
        //
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            (ui+jdegr*nvar)[2] = fmax(ELIMIT, eave  + ((ui+jdegr*nvar)[2] - eave) * er);
        }
        //
        badcell = false;
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            badcell = IsBad(ui+jdegr*nvar) or badcell;
            if (badcell) {printf("u:%e,  %e,  %e\n",(ui+jdegr*nvar)[0],(ui+jdegr*nvar)[1],(ui+jdegr*nvar)[2]);}
            //
        }
        ASSERT(not badcell, "the cell is not fixed dummy");        
    }

}
