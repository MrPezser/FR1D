#include <iostream>
#include <cmath>

#include "indexing.h"
#include "basis.h"
#include "EulerFlux.h"

//Swap this out to change equations. If you want to compare different PDEs, you've come to the wrong place buddy.
#define GAM 1.4
#define FACEFLUX(a, b, c) RoeFDS(GAM, a, b, c) //AdvectionFaceFlux(a, b, c)
#define FLUX(a,b)  EulerFlux(GAM, a, b) //AdvectionFlux(a, b)

using namespace std;

void veccopy(double* a, const double* b, size_t n){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

double Initialize(double x){
    ///This defines the intial state of the solution space

    /*
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
    }*/

    //return 2.0 + sin(2.0*M_PI*x);

    return 1.0 + exp(-40*(x-0.5)*(x-0.5));
}

void InitializeEuler(double x, double* u){
    double rho = 1.0;
    double v = 1.0;
    double p = 1.0;

    //rho = Initialize(x);

    if (x < 0.5){
        //rho = 1.0;
        v = 0.0;
        //p = 1.0;
    } else {
        rho = 0.125;
        v = 0.0;
        p = 0.1;
    }

    u[0] = rho;                             //rho
    u[1] = rho * v;                         //rho V
    u[2] = 0.5*rho*v*v + (p/(GAM-1.0));       //rho e
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

void CalcDudt(const int nx, const int ndegr, const int nvar, const double a, const double dx, const double* u, const double* Dmatrix, const double* Dradau, double* dudt ){
    ///Calculates the solution update given a function to find flux
    int nu = nx * ndegr * nvar;
    int np = ndegr * nvar;
    double fL[nvar], fR[nvar], common_flux[nvar], fcorr_xi[nu], fxi[nvar];
    double flux_node[np];

    for (int inode=0; inode<ndegr; inode++){
        for (int kvar=0; kvar<nvar;kvar++){
            flux_node[iup(inode, kvar, nvar)] = 0.0;
        }
    }

    //Initialize dudt & corrected flux slope
    for (int i=0;i<nu; i++){
        dudt[i] = 0.0;
        fcorr_xi[i] = 0.0;
    }

    //Flux Reconstruction (Pn)

    //Compute face contributions to the corrected flux gradient
    for (int iface=0; iface<nx; iface++){
        int ielem, iep1;
        ielem = iface;
        if (iface == nx-1) {
            //Periodic boundary condition
            iep1 = 0;
        } else {
            //Interior Cell
            iep1 = iface+1;
        }


        //Calculate the common flux at the face
        FACEFLUX(&u[iu3(ielem, ndegr-1, 0, ndegr)], &u[iu3(iep1, 0, 0, ndegr)], common_flux);
        //common_flux[0] = FACEFLUX(a, u[iu3(iem1, ndegr-1, 0, ndegr)], u[iu3(ielem, 0, 0, ndegr)]);

        //Calculate the element's local fluxes at the face
        FLUX(&u[iu3(ielem , ndegr-1 , 0, ndegr)], &fL[0]);
        FLUX(&u[iu3(iep1, 0       , 0, ndegr)], &fR[0]);
        //fL[0] = FLUX(a, u[iu3(iem1, ndegr-1, 0, ndegr)]);
        //fR[0] = FLUX(a, u[iu3(ielem, 0, 0, ndegr)]);


        for (int inode=0; inode<ndegr; inode++){
            for (int kvar = 0; kvar<nvar; kvar++) {
                fcorr_xi[iu3(ielem,ndegr-1-inode, kvar, ndegr)] -= (common_flux[kvar] - fL[kvar]) * Dradau[inode];   //ndegr-1-
                fcorr_xi[iu3(iep1, inode        , kvar, ndegr)] += (common_flux[kvar] - fR[kvar]) * Dradau[inode];
            }
        }
    }

    //cell-internal components and calculating dudt
    for (int ielem=0; ielem<nx; ielem++){

        for (int inode=0; inode<ndegr; inode++){
            double flux[3];
            FLUX(&u[iu3(ielem, inode , 0, ndegr)], &flux[0]);
            //flux_node[inode][0] = FLUX(a, u[iu3(ielem, inode , 0, ndegr)]);

            flux_node[iup(inode, 0, nvar)] = flux[0];
            flux_node[iup(inode, 1, nvar)] = flux[1];
            flux_node[iup(inode, 2, nvar)] = flux[2];
        }

        //compute discontinuous flux slope (in reference element of width 2)
        for (int inode=0; inode<ndegr ;inode++) {
            for (int kvar = 0; kvar < nvar; kvar++) {
                fxi[kvar] = 0.0;
            }

            //Use the Lagrange polynomial derivative matrix to perform the MatVec -> interior contribution to flux slope
            for (int jnode = 0; jnode < ndegr; jnode++) {
                //if (jnode == inode) continue;
                for (int kvar = 0; kvar < nvar; kvar++) {
                    fxi[kvar] += flux_node[iup(jnode, kvar, nvar)] * Dmatrix[iup(inode, jnode, ndegr)];
                }
            }

            // add the interior contribution to the flux slope
            for (int kvar = 0; kvar < nvar; kvar++) {
                fcorr_xi[iu3(ielem, inode, kvar, ndegr)] += fxi[kvar];
            }
        }

        for (int inode=0; inode<ndegr ;inode++) {
            for (int kvar = 0; kvar<nvar; kvar++) {

                //convert from the reference element to the real element and negate to put make it dudt (= -f_x)
                dudt[iu3(ielem, inode, kvar, ndegr)] = -fcorr_xi[iu3(ielem,  inode, kvar, ndegr)] * (2/ dx);

                 //fix end points for shock-tube like problem
                if (ielem == 0 || ielem == nx-1){
                    dudt[iu3(ielem, inode, kvar, ndegr)] = 0;
                }

                if (isnan(dudt[iu3(ielem, inode, kvar, ndegr)])){
                    throw overflow_error("dudt NAN\n");
                }
            }

        }
    }

}

int main() {
    ///hardcoded inputs
    int     nx = 50;           //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied domain from x=0 to x=1

    int ndegr = 2;             //Degrees of freedom per element
    int nvar = 3;               //Number of variables
    int nu = nx * ndegr * nvar;

    double cfl = 0.01 / (ndegr*ndegr);          //CFL Number
    double a = 1.0;             //Wave Speed

    double tmax = 0.1;
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

    //Allocate Arrays
    auto* x = (double*)malloc(nx*sizeof(double));
    auto* u = (double*)malloc(nu*sizeof(double));
    auto* u0 = (double*)malloc(nu*sizeof(double));
    auto* dudt = (double*)malloc(nu*sizeof(double));


    //Generate Grid (currently uniform) & initialize solution
    for (int i=0; i<nx; i++){
        //defining x position of cell centers
        x[i] = (i+0.5) * dx;
        for (int j=0; j<ndegr; j++) {
            InitializeEuler(x[i] + xi[j]*(0.5 * dx), &u[iu3(i, j, 0, ndegr)]);
            //u[iu3(i,j,0, ndegr)] = Initialize(x[i] + xi[j]*(0.5 * dx));
        }
    }

    veccopy(u0,u,nu);

    // Begin Time Marching (3 stage TVD RK)
    auto* u_tmp = (double*)malloc(nu*sizeof(double));

    for (int iter=0; iter<niter; iter++){
        veccopy(u_tmp, u, nu);
        //1st stage
        CalcDudt(nx, ndegr, nvar, a, dx, u, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            //u_tmp[i] += dt * dudt[i];
            u[i] += dt * dudt[i];

            if (isnan(u[i])){
                throw overflow_error("Kaboom!\n");
            }
        }
        /*
        //2nd stage
        CalcDudt(nx, ndegr, nvar, a, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u_tmp[i] = 0.75*u[i] + 0.25*( u_tmp[i] + dt*dudt[i]);
        }

        //3rd stage
        CalcDudt(nx, ndegr, nvar, a, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u[i] = (1.0/3.0)*u[i] + (2.0/3.0)*(u_tmp[i] + dt*dudt[i]);
        }*/

        if (iter % 100 == 0){printf("iter:%10d\t%7.2f%% Complete\n",iter, 100.0*(double)iter/(double)niter);}
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Final Solution
    FILE* fout = fopen("waveout.tec", "w");
    fprintf(fout, "x\tu\tu0\n");

    for (int i=0;i<nx;i++) {
        for (int j=0; j<ndegr; j++) {
            double xj = x[i] + xi[j]*(0.5 * dx);
            fprintf(fout, "%f\t%f\t%f\t%f\t", xj , u[iu3(i,j,0,ndegr)],  u[iu3(i,j,1,ndegr)],  u[iu3(i,j,2,ndegr)]);

            double rho, v, p, c, M;
            getPrimativesPN(GAM , &u[iu3(i,j,0,ndegr)], &rho, &v, &p, &c, &M);

            fprintf(fout, "%f\t%f\t%f\n", p, c, M);
            //fprintf(fout, "%f\t%f\t", xj , u[iu3(i,j,0,ndegr)]);
            //fprintf(fout, "%f\n", u0[iu3(i,j,0,ndegr)]);
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
