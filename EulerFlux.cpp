//
// Created by Tsail on 6/24/2023.
//

#include <cmath>
#include <stdexcept>

#include "EulerFlux.h"

using namespace std;

void getPrimativesPN(const double gam, const double *unkel, double *rho, double *v, double *p, double *c, double *M) {
    // Gets the primative variables at a specific node
    rho[0]       = unkel[0];
    double rhov  = unkel[1];
    double rhoe  = unkel[2];

    //Density Limiter
    //rho[0] = fmax(1e-8, rho[0]);

    //break down into primatives
    v[0] = rhov / rho[0];
    double v2 = v[0]*v[0];

    p[0] = (gam - 1) * (rhoe - (0.5 * rhov * v[0]));

    //Pressure Limit
    p[0] = fmax(1e-8, p[0]);

    c[0] = sqrt(gam * p[0] / fmax(1e-8, rho[0]));
    M[0] = v[0] / c[0];

    if (isnan(M[0])) {throw overflow_error("Nonpositive Density\n");}
}

void EulerFlux(double gam, const double *u, double* flux){
    flux[0] = 0.0;
    flux[1] = 0.0;
    flux[2] = 0.0;

    double rho, v, p, c, M;
    getPrimativesPN(gam, u, &rho, &v, &p, &c, &M);

    flux[0] = u[1];
    flux[1] = (u[1]*v) + p;
    flux[2] = v * (u[2] + p);

}

void RoeFDS(double gam, const double* uL, const double *uR, double* roeFlux){

    double rhoL, vL, pL, cL, ML;
    getPrimativesPN(gam, uL, &rhoL, &vL, &pL, &cL, &ML);
    double hL = (uL[2] + pL) / uL[0];

    double rhoR, vR, pR, cR, MR;
    getPrimativesPN(gam, uR, &rhoR, &vR, &pR, &cR, &MR);
    double hR = (uR[2] + pR) / uR[0];

    //find Roe-averaged state vars
    double R, rho, v, h0, c, v2;
    R = sqrt(rhoR / rhoL);
    rho = R * rhoL;
    v  = (R*vR + vL) / (R+1);
    h0 = (R*hR + hL) / (R+1);
    v2 = v*v;
    c  = sqrt((gam-1.0)*(h0 - 0.5*v2));

    double deltaU = vR - vL;
    double deltaP = pR - pL;

    //find the approx flux jacobian eigenvalues with the Harten correction
    double lam[3], phi[3];
    lam[0] = v;
    lam[1] = v + c;
    lam[2] = v - c;

    double eps = 1e-6;//deltaU;

    for (int i=0; i<3; i++){
        if (fabs(lam[i]) < eps){
            phi[i] = 0.5*(lam[i]*lam[i]/eps + eps);//(lam[i]*lam[i] + eps*eps) / (2*eps);
        } else {
            phi[i] = fabs(lam[i]);
        }
    }

    //transformed conservative vars
    double dw1, dw2, dw3;
    dw1 = (rhoR-rhoL) - (deltaP / (c*c));
    dw2 = deltaU + (deltaP / (rho*c));
    dw3 = deltaU - (deltaP / (rho*c));

    //eigenvectors of approx jac matrix
    double e1[3], e2[3], e3[3];

    e1[0] = 1;
    e1[1] = v;
    e1[2] = 0.5*v2;

    double coeff = 0.5*rho/c;
    e2[0] = coeff;
    e2[1] = coeff*(v + c);
    e2[2] = coeff*(h0 + (v*c));

    coeff = -coeff;
    e3[0] = coeff;
    e3[1] = coeff * (v - c);
    e3[2] = coeff * (h0 - (v*c));

    //Calculate left and right flux states
    /*auto*fluxLeft = (double*)malloc(3*sizeof(double));
    auto*fluxRight = (double*)malloc(3*sizeof(double));
    EulerFlux(gam, uL, fluxLeft);
    EulerFlux(gam, uR, fluxRight);*/
    double fluxLeft[3], fluxRight[3];
    fluxLeft[0]  = uL[1];
    fluxRight[0] = uR[1];

    fluxLeft[1]  = (((gam - 1) / gam) * rhoL * hL)  +  ((0.5 * (gam + 1) / gam) * uL[1] * vL);
    fluxRight[1] = (((gam - 1) / gam) * rhoR * hR)  +  ((0.5 * (gam + 1) / gam) * uR[1] * vR);

    fluxLeft[2]  = uL[1]*hL;
    fluxRight[2] = uR[1]*hR;

    //Compile the final Roe flux


    roeFlux[0] = 0.5 * (fluxRight[0] + fluxLeft[0] - (e1[0]*phi[0]*dw1) - (e2[0]*phi[1]*dw2) - (e3[0]*phi[2]*dw3));
    roeFlux[1] = 0.5 * (fluxRight[1] + fluxLeft[1] - (e1[1]*phi[0]*dw1) - (e2[1]*phi[1]*dw2) - (e3[1]*phi[2]*dw3));
    roeFlux[2] = 0.5 * (fluxRight[2] + fluxLeft[2] - (e1[2]*phi[0]*dw1) - (e2[2]*phi[1]*dw2) - (e3[2]*phi[2]*dw3));

}

double F1pm(const int isPlus, const double M, const double rho, const double c){

    if (isPlus==1){
        if (M<=-1.0) {
            //F1+
            return 0.0;
        }
        if (M>=1.0) {
            return rho * M * c;
        }
        return 0.25*rho*c*(M+1)*(M+1);
    }

    if (isPlus==0) {
        if (M <= -1.0) {
            //F1-
            return rho * M * c;// *-1
        }
        if (M >= 1.0) {
            return 0.0;
        }
        return -0.25 * rho * c * (M - 1) * (M - 1);
    }

    return NAN;
}

void LeerFlux(double gam, const double* uL, const double* uR, double *flux){
    double F1L, F1R, fPlus[3], fMins[3];

    flux[0] = 0.0;
    flux[1] = 0.0;
    flux[2] = 0.0;

    double rhoL, vL, pL, cL, ML;
    getPrimativesPN(gam, uL, &rhoL, &vL, &pL, &cL, &ML);

    double rhoR, vR, pR, cR, MR;
    getPrimativesPN(gam, uR, &rhoR, &vR, &pR, &cR, &MR);


    if (ML >= 1.0){
        EulerFlux(gam, uL, flux);
        return;
    }

    if (MR <= -1.0){
        EulerFlux(gam, uR, flux);
        return;
    }

    F1L = F1pm(1, ML, rhoL, cL);
    F1R = F1pm(0, MR, rhoR, cR);

    if(isnan(F1L+F1R)){
        throw overflow_error("getting NAN mass flux!");
    }

    fPlus[0] = F1L;
    double A = ((gam - 1) * vL) + (2.0 * cL); //vL
    fPlus[1] = F1L * A / gam;
    fPlus[2] = F1L * A*A * 0.5 / (gam*gam - 1.0);


    fMins[0] = F1R;
    A = -((gam - 1) * vL) - (2.0 * cR); //vR
    fMins[1] = F1R * A / gam;
    fMins[2] = F1R * A*A * 0.5 / (gam*gam - 1.0);

    flux[0] = fPlus[0] + fMins[0];
    flux[1] = fPlus[1] + fMins[1];
    flux[2] = fPlus[2] + fMins[2];
}
