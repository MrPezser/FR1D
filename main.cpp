#include <iostream>
#include <cmath>

#include "indexing.h"
#include "basis.h"

//Swap this out to change equations. If you want to compare different PDEs, you've come to the wrong place buddy.
#define FACEFLUX(a, b, c) AdvectionFaceFlux(a, b, c)
#define FLUX(a,b) AdvectionFlux(a,b)

using namespace std;

double Initialize(double x){
    ///This defines the intial state of the solution space

    /*//Gaussian bump and step combo
    if (x < 0.6) {
        double beta = 0.01;
        return 1 + exp(-(x-0.3)*(x-0.3) / beta);
    } else {
        if (x < 0.8) {
            return 2.0;
        } else {
            return 1.0;
        }
    }
     */

    //return sin(2.0*M_PI*x);

    return exp(-40*(x-0.5)*(x-0.5));
}

double AdvectionFaceFlux(const double a, const double uL, const double uR) {
    /// Finds the upwind flux according to the 1D advection equation
    ///a = wave speed | uL & uR are the left and right state variables
    if (a > 0){
        return a * uL;
    } else {
        return a * uR;
    }
}

double AdvectionFlux(const double a, double u) {
    /// Finds the flux according to the 1D advection equation
    ///a = wave speed | u is state variable

    return a * u;
}

void CalcDudt(const int nx, const int ndegr, const double a, const double dx, const double* u, const double* Dmatrix, const double* Dradau, double* dudt ){
    ///Calculates the solution update given a function to find flux
    double uL, uR, flux;
    int nu = ndegr * nx;

    //Initialize dudt
    for (int i=0;i<nu; i++){
        dudt[i] = 0;
    }

    /*
    //Finite Volume
    //loop through faces (convention is left,-, face of element i)
    for (int iface=0; iface<nx; iface++){
        if (iface == 0) {
            //Periodic boundary condition
            ifm1 = nx-1;
        } else {
            //Interior Cell
            ifm1 = iface-1;
        }

        //Get the left and right states
        uL = u[ifm1];
        uR = u[iface];

        //Calculate the flux at the face
        flux = FACEFLUX(a, uL, uR);

        //Add the flux contribution to the RHS
        dudt[ifm1]    -= (flux / dx);
        dudt[iface]   += (flux / dx);
    }*/

    //Flux Reconstruction (P1)
    //find&store the common upwind fluxes at each face
    auto common_flux = (double*)malloc(nu*sizeof(double));

    for (int iface=0; iface<nx; iface++){
        int ielem, iem1;
        ielem = iface;
        if (iface == 0) {
            //Periodic boundary condition
            iem1 = nx-1;
        } else {
            //Interior Cell
            iem1 = iface-1;
        }
        //Get the left and right states
        uL = u[iup(iem1 , ndegr-1 , ndegr)];
        uR = u[iup(ielem, 0       , ndegr)];

        //Calculate the flux at the face
        common_flux[iface] = FACEFLUX(a, uL, uR);
    }

    //find corrected flux and calc dudt
    for (int ielem=0; ielem<nx; ielem++){
        ///using linear 'hat' basis functions for P1
        ///reference slope = 1/dx
        ///Solution points = cell end points
        double fL, fR, fxi[ndegr], fcorr_xi[ndegr], flux_node[ndegr];
        int iface, ifp1;

        iface = ielem; // left face
        if (ielem == nx-1) {ifp1 = 0;} else {ifp1 = iface + 1;}

        for (int inode=0; inode<ndegr; inode++){
            flux_node[inode] = FLUX(a, u[iup(ielem, inode ,ndegr)]);
        }

        fL = flux_node[0];
        fR = flux_node[ndegr-1];

        //flux slope (in reference element of width 2) -----   for p1 this is constant
        for (int inode=0; inode<ndegr ;inode++){
            //Use the Lagrange polynomial derivative to perform the MatVec
            fxi[inode] = 0.0;
            for (int jnode=0; jnode<ndegr; jnode++) {
                fxi[inode] += flux_node[jnode] * Dmatrix[iup(inode, jnode, ndegr)];
            }
        }

        for (int inode=0; inode<ndegr; inode++){
            fcorr_xi[inode] = fxi[inode] + (common_flux[iface]-fL)*Dradau[inode] + (common_flux[ifp1]-fR)*Dradau[ndegr-inode-1];
        }

        for (int j=0; j<ndegr ;j++) {
            dudt[iup(ielem, j, ndegr)] = -(2 / dx) * fcorr_xi[j];
        }
    }

    free(common_flux);
}

int main() {
    ///hardcoded inputs
    int     nx = 10;           //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied domain from x=0 to x=1

    int  ndegr = 4;             //Degrees of freedom per element
    int nu = nx * ndegr;

    double cfl = 0.01 / (ndegr*ndegr);          //CFL Number
    double a = 1.0;             //Wave Speed

    double tmax = 50.0;
    double dt = (cfl * dx) / a;
    int niter = ceil(tmax/dt);  //Guess number of iterations required to get to the given tmax


    //Find the solution points in reference space
    auto* xi = (double*)malloc(ndegr*sizeof(double));
    GenerateLobattoPoints(ndegr, xi);

    //Find the derivative matrix of the associated lagrange polynomials
    //Derivatie of L_j at position x_i
    auto* Dmatrix = (double*)malloc(ndegr*ndegr*sizeof(double));
    GenerateLagrangeDMatrix(ndegr, xi, Dmatrix);

    //Find the derivatives of the Radau polynomial of appropriate degree
    auto* Dradau = (double*)malloc((1+ndegr)*sizeof(double));
    GenerateRadauDerivatives(ndegr, xi, Dradau);

    //for (int i=0;i<ndegr; i++){
    //    printf("xi[%d] = %f \t Dradau[%d] = %f\n",i, xi[i], i, Dradau[i]);
    //}
    //return 0;

    //Allocate Arrays
    auto* x = (double*)malloc(nx*sizeof(double));
    auto* u = (double*)malloc(nu*sizeof(double));
    auto* u0 = (double*)malloc(nu*sizeof(double));
    auto* dudt = (double*)malloc(nu*sizeof(double));

    int itype = 0; //nodal basis


    //Generate Grid (currently uniform) & initialize solution
    for (int i=0; i<nx; i++){
        //defining x position of cell centers
        x[i] = (i+0.5) * dx;

        for (int j=0; j<ndegr; j++) {
            //Working Solution
            u[iup(i, j, ndegr)]  = Initialize(x[i] + xi[j]*(0.5 * dx));
            //Initial Solution
            u0[iup(i, j, ndegr)] = Initialize(x[i] + xi[j]*(0.5 * dx));
        }
    }


    // Begin Time Marching
    for (int iter=0; iter<niter; iter++){

        //Explicit Euler
        CalcDudt(nx, ndegr, a, dx, u, Dmatrix, Dradau, dudt);

        for (int i=0; i<nu; i++){
            u[i] += dt * dudt[i];
            if (isnan(u[i])){
                throw overflow_error("Kaboom!\n");
            }
        }

        if (iter % 50000 == 0){printf("%7.2f%% Complete\n", 100.0*(double)iter/(double)niter);}
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Final Solution
    FILE* fout = fopen("waveout.tec", "w");
    fprintf(fout, "x\tu\tu0\n");

    for (int i=0;i<nx;i++) {
        for (int j=0; j<ndegr; j++) {
            fprintf(fout, "%f\t%f\t%f\n", x[i] + xi[j]*(0.5 * dx), u[iup(i,j,ndegr)], u0[iup(i,j,ndegr)]);
        }
    }

    fclose(fout);

    free(xi);
    free(u);
    free(u0);
    free(dudt);
    free(Dmatrix);
    free(Dradau);
}
