#include <iostream>
#include <cmath>

//Swap this out to change equations. If you want to compare different PDEs, you've come to the wrong place buddy.
#define FLUX(a, b, c) AdvectionFlux(a, b, c)

using namespace std;

double Initialize(double x){
    ///This defines the intial state of the solution space

    //Gaussian bump and step combo
    if (x < 0.8) {
        double beta = 0.01;
        return exp(-(x-0.4)*(x-0.4) / beta);
    } else {
        return 1.0;
    }
}

double AdvectionFlux(const double a, const double uL, const double uR) {
    /// Finds the upwind flux according to the 1D advection equation
    ///a = wave speed | uL & uR are the left and right state variables
    if (a > 0){
        return a * uL;
    } else {
        return a * uR;
    }
}

void CalcDudt(const int nx, const double a, const double dx, const double* u, double* dudt ){
    ///Calculates the solution update given a function to find flux
    double uL, uR, flux;
    int ifm1;
    //Initialize dudt
    for (int i=0;i<nx; i++){
        dudt[i] = 0;
    }


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
        flux = FLUX(a, uL, uR);

        //Add the flux contribution to the RHS
        dudt[ifm1]    -= (flux / dx);
        dudt[iface]   += (flux / dx);
    }
}

int main() {
    ///hardcoded inputs
    int     nx = 100;            //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied range from x=0 to x=1

    double cfl = 0.75;          //CFL Number
    double a = 1.0;             //Wave Speed

    double tmax = 1.0;
    double dt = cfl * dx / a;
    int niter = round(tmax/dt);


    int  ndegr = 1;             //Degrees of freedom per element
    int nu = nx * ndegr;

    //Allocate Arrays
    auto* x = (double*)malloc(nx*sizeof(double));
    auto* u = (double*)malloc(nu*sizeof(double));
    auto* dudt = (double*)malloc(nu*sizeof(double));

    //Generate Grid (currently uniform) & initialize solution
    for (int i=0; i<nx; i++){
        //defining x position of cell centers
        x[i] = (i+0.5) * dx;

        u[i] = Initialize(x[i]);
    }

    // Begin Time Marching
    for (int iter=0; iter<niter; iter++){

        //Explicit Euler
        CalcDudt(nx, a, dx, u, dudt);
        for (int i=0; i<nu; i++){
            u[i] += dt * dudt[i];
        }
    }

    //Printout Final Solution
    FILE* fout = fopen("waveout.dat", "w");
    fprintf(fout, "x\tu\n");

    for (int i=0;i<nx;i++) {
        fprintf(fout, "%f\t%f\n", x[i]-(0.5*dx), u[i]);
        fprintf(fout, "%f\t%f\n", x[i]+(0.5*dx), u[i]);
    }

    fclose(fout);
}
