#include <iostream>
#include <cmath>

#include "indexing.h"

//Swap this out to change equations. If you want to compare different PDEs, you've come to the wrong place buddy.
#define FACEFLUX(a, b, c) AdvectionFaceFlux(a, b, c)
#define FLUX(a,b) AdvectionFlux(a,b)

using namespace std;

double Initialize(double x){
    ///This defines the intial state of the solution space

    //Gaussian bump and step combo
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

void CalcDudt(const int nx, const int ndegr, const double a, const double dx, const double dt, const double* u, double* dudt ){
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
        uL = u[iup(iem1,  1, 2)];
        uR = u[iup(ielem, 0, 2)];

        //Calculate the flux at the face
        common_flux[iface] = FACEFLUX(a, uL, uR);
    }

    //find corrected flux and calc dudt
    for (int ielem=0; ielem<nx; ielem++){
        ///using linear 'hat' basis functions for P1
        ///reference slope = 1/dx
        ///Solution points = cell end points
        double fL, fR, fxi, fcorr_xi[2];
        int iface, ifp1;
        iface = ielem; // left face
        if (ielem == nx-1) {ifp1 = 0;} else {ifp1 = iface + 1;}

        fL = FLUX(a, u[iup(ielem,0,2)]);
        fR = FLUX(a, u[iup(ielem,1,2)]);

        //flux slope (in reference element of width 2) -----   for p1 this is constant
        fxi = 0.5 * (fR - fL);

        /*
        //make corrections to the flux gradient using the correction function and common flux
        //using g_corr = Radau to simulate DG, P1 requires us to use 2nd order radau
        // !!! LEFT boundary uses RIGHT radau poly !!!! and likewise
        //Rr2 = (1/2) * (1.5xi^2 - xi- 0.5);
        //gr_xi =  1.5xi - 0.5; = -2.0 (@xi = -1)   = 1.0   (@xi = 1)     LEFT BOUNDARY
        //gl_xi = -1.5xi - 0.5; =  1.0 (@xi = -1)   = -2.0  (@xi = 1)     RIGHT BOUNDARY
        //The discrete flux fxn at the left and right boundaries is on the solution points
         */
        fcorr_xi[0] = (fxi + (common_flux[iface] - fL)*(-2.0) + (common_flux[ifp1] - fR)*( 1.0));   //xi = -1
        fcorr_xi[1] = (fxi + (common_flux[iface] - fL)*( 1.0) + (common_flux[ifp1] - fR)*(-2.0));   //xi =  1

        dudt[iup(ielem, 0, 2)] =  -(2/dx) * fcorr_xi[0];
        dudt[iup(ielem, 1, 2)] =  -(2/dx) * fcorr_xi[1];
    }

    free(common_flux);
}

int main() {
    ///hardcoded inputs
    int     nx = 100;            //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied domain from x=0 to x=1

    double cfl = 0.03;          //CFL Number
    double a = 1.0;             //Wave Speed

    double tmax = 1.0;
    double dt = (cfl * dx) / a;
    int niter = ceil(tmax/dt);


    int  ndegr = 2;             //Degrees of freedom per element
    int nu = nx * ndegr;

    //Allocate Arrays
    auto* x = (double*)malloc(nx*sizeof(double));
    auto* u = (double*)malloc(nu*sizeof(double));
    auto* u0 = (double*)malloc(nu*sizeof(double));
    auto* dudt = (double*)malloc(nu*sizeof(double));

    //Generate Grid (currently uniform) & initialize solution
    for (int i=0; i<nx; i++){
        //defining x position of cell centers
        x[i] = (i+0.5) * dx;

        u[iup(i,0,ndegr)] = Initialize(x[i]-(0.5*dx));
        u[iup(i,1,ndegr)] = Initialize(x[i]+(0.5*dx));
        u0[iup(i,0,ndegr)] = Initialize(x[i]-(0.5*dx));
        u0[iup(i,1,ndegr)] = Initialize(x[i]+(0.5*dx));
    }

    // Begin Time Marching
    for (int iter=0; iter<niter; iter++){

        //Explicit Euler
        CalcDudt(nx, ndegr, a, dx, dt, u, dudt);
        for (int i=0; i<nu; i++){
            u[i] += dt * dudt[i];
        }
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Final Solution
    FILE* fout = fopen("waveout.tec", "w");
    fprintf(fout, "x\tu\tu0\n");

    for (int i=0;i<nx;i++) {
        fprintf(fout, "%f\t%f\t%f\n", x[i]-(0.5*dx), u[iup(i,0,ndegr)], u0[iup(i,0,ndegr)]);
        fprintf(fout, "%f\t%f\t%f\n", x[i]+(0.5*dx), u[iup(i,1,ndegr)], u0[iup(i,1,ndegr)]);
    }

    fclose(fout);
}
