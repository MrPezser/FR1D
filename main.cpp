#include <iostream>
#include <cmath>

#include "indexing.h"
#include "basis.h"
#include "EulerFlux.h"
#include "SpatialDiscretization.h"


using namespace std;

void veccopy(double* a, const double* b, size_t n){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

double Initialize(double x){
    ///This defines the intial state of the solution space


    //Gaussian bump and step combo

    if (x < 0.6) {
        double beta = 0.005;
        return 1.0;// + exp(-(x-0.3)*(x-0.3) / beta);
    } else {
        if (x < 0.8) {
            return 2.0;
        } else {
            return 1.0;
        }
    }


    //return 2.0 + sin(2.0*M_PI*x);

    //return 1.0 + exp(-40*(x-0.5)*(x-0.5));

    if (x<0.5){
        return 1.0;
    } else {
        return 0.0;
    }
}

void InitializeEuler(double x, double* u){
    double rho = 1.0;
    double v = 1.0;
    double p = 1.0;

    //rho = Initialize(x) + 10.0;

    if (x < 0.5){
        rho = 1.0;
        v = 0.0;
        p = 1.0;
    } else {
        rho = 1.0 / 211.0 ;//0.125;
        v = 0.0;
        p = rho;
    }


    u[0] = rho;                             //rho
    u[1] = rho * v;                         //rho V
    u[2] = 0.5*rho*v*v + (p/(GAM-1.0));       //rho e
}

int main() {
    ///hardcoded inputs
    int     nx = 2;           //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied domain from x=0 to x=1

    int ndegr = 20;             //Degrees of freedom per element
    int nvar = NVAR;              //Number of variables
    int nu = nx * ndegr * nvar;

    double cfl = 0.1/(ndegr*ndegr);          //CFL Number

    double tmax = 10.4;
    double dt = (cfl * dx); // /a;
    int niter = ceil(tmax/dt);  //Guess number of iterations required to get to the given tmax //10*3*80

    //Find the solution points in reference space
    auto* xi = (double*)malloc(ndegr*sizeof(double));
    GenerateLobattoPoints(ndegr, xi);

    //Find the derivative matrix of the associated lagrange polynomials
    //Derivatie of L_j at position x_i
    auto* Dmatrix = (double*)malloc(ndegr*ndegr*sizeof(double));
    GenerateLagrangeDMatrix(ndegr, xi, Dmatrix);

    //Find the derivatives of the Radau polynomial of appropriate degree
    auto* Dradau = (double*)malloc(2*(1+ndegr)*sizeof(double));
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
            //InitializeEuler(x[i], &u[iu3(i, j, 0, ndegr)]);    //+ xi[j]*(0.5 * dx)

            if (x[i] > 0.6) {
                u[iu3(i, j, 0, ndegr)] = Initialize(x[i] + (0.0 * dx)); //enforce a sharp initial discon
            } else {
                u[iu3(i, j, 0, ndegr)] = Initialize(x[i] + xi[j]*(0.5 * dx));
            }
            //Initialize(x[i] + xi[j]*(0.5 * dx));  //will allow a slop initial cond.
        }
    }

    veccopy(u0,u,nu);

    // Begin Time Marching (3 stage TVD RK)
    auto* u_tmp = (double*)malloc(nu*sizeof(double));

    for (int iter=0; iter<niter; iter++){
        veccopy(u_tmp, u, nu);
        //1st stage
        CalcDudt(nx, ndegr, nvar, dx, u, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            //u_tmp[i] += dt * dudt[i];
            u[i] += dt * dudt[i];

            if (isnan(u[i])){
                throw overflow_error("Kaboom!\n");
            }
        }



        //2nd stage
        CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u_tmp[i] = 0.75*u[i] + 0.25*( u_tmp[i] + dt*dudt[i]);
        }

        //3rd stage
        CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
        for (int i=0; i<nu; i++){
            u[i] = (1.0/3.0)*u[i] + (2.0/3.0)*(u_tmp[i] + dt*dudt[i]);
        }

        if (iter % 100 == 0){printf("iter:%10d\t%7.2f%% Complete\n",iter, 100.0*(double)iter/(double)niter);}
    }

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Final Solution
    FILE* fout = fopen("waveout.tec", "w");
    fprintf(fout, "x\tu\tu0\n");

    for (int i=0;i<nx;i++) {
        for (int j=0; j<max(ndegr,2); j++) {
            double xj;

            if (nvar == 1) {
                if (ndegr == 1) {
                    double xii = -1.0 + 2.0 * j;
                    xj = x[i] + xii * (0.5 * dx);
                    fprintf(fout, "%f\t%f\n", xj, u[iu3(i, 0, 0, ndegr)]);
                } else {
                    xj = x[i] + xi[j] * (0.5 * dx);
                    fprintf(fout, "%f\t%f\n", xj, u[iu3(i, j, 0, ndegr)]);
                }
            } else {
                xj = x[i] + xi[j] * (0.5 * dx);
                fprintf(fout, "%f\t%f\t%f\t%f\t", xj, u[iu3(i, j, 0, ndegr)], u[iu3(i, j, 1, ndegr)],
                        u[iu3(i, j, 2, ndegr)]);

                double rho, v, p, c, M;
                getPrimativesPN(GAM, &u[iu3(i, j, 0, ndegr)], &rho, &v, &p, &c, &M);

                fprintf(fout, "%f\t%f\t%f\n", p, c, M);
            }
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
