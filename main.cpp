#include <iostream>
#include <cmath>

#include "indexing.h"
#include "basis.h"
#include "EulerFlux.h"
#include "SpatialDiscretization.h"
#include "initial.h"
#include "Implicit.h"

#define EEULER   (1)
#define TVDRK3   (3)
#define SEMIIMP  (10)
#define SITVDRK3 (30)

using namespace std;

void veccopy(double* a, const double* b, size_t n){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

double dt_from_CFL(double dt0, int nelem, int ndegr, double* u){
  double gam = 1.4;
  double dt = 9999999999.0;

  for (int ielem=0; ielem<nelem; ielem++){
    for (int j=0; j<ndegr; j++){
      //
      double *uij = &(u[iu3(ielem,j,0,ndegr)]);
      //
      double rho, v, p, c, M;
      getPrimativesPN(gam, uij, &rho, &v, &p, &c, &M);
      //
      dt = fmin(dt,dt0/c);
      //printf("ix,jdeg,dt,dt0,c : %d,%d,%f,%f,%f\n",ielem,j,dt,dt0,c);
      //printf("uij:rho,v,c,M: %f,%f,%f : %f,%f,%f,%f\n",uij[0], uij[1], uij[2],rho,p,c,M);
    }
  }

  return dt;

}

int main() {
    ///hardcoded inputs
    int     nx = 750;           //Number of elements, nx+1 points
    double  dx = 1.0 / nx;      //Implied domain from x=0 to x=1

    int ndegr = 3;             //Degrees of freedom per element
    int nvar = NVAR;              //Number of variables
    int nu = nx * ndegr * nvar;

    double cfl = 0.1/(ndegr*ndegr);          //CFL Number

    double tmax = 0.2;
    double dt = (cfl * dx); // /a;
    double dt0 = dt;
    int niter = ceil(tmax/dt);  //Guess number of iterations required to get to the given tmax //10*3*80
    int mxiter = 1e6;

    double gam = 1.4;//warning, not global
    
    int timestepping = TVDRK3;
    //timestepping = SEMIIMP;
    //timestepping = SITVDRK3;
    //timestepping = EEULER;

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
            if (nvar ==3) {
                InitializeEuler(x[i]+ xi[j]*(0.5 * dx), &u[iu3(i, j, 0, ndegr)]);
            } else {

                //if (x[i] > 0.6) {
                //    u[iu3(i, j, 0, ndegr)] = Initialize(x[i] + (0.0 * dx)); //enforce a sharp initial discon
                //} else {
                    u[iu3(i, j, 0, ndegr)] = Initialize(x[i] + xi[j]*(0.5 * dx));
                //}
            //Initialize(x[i] + xi[j]*(0.5 * dx));  //will allow a slop initial cond.
            }
        }
    }

    veccopy(u0,u,nu);

    // Begin Time Marching 
    auto* u_tmp = (double*)malloc(nu*sizeof(double));
    
    double soltime = 0.0;
    
    //SmoothSolution(nx, ndegr, nvar, u);

    for (int iter=0; iter<mxiter; iter++){
        //break;
        if (nvar ==3) {
            dt = dt_from_CFL(dt0, nx, ndegr, u);
        }
        soltime += dt;

        if (timestepping == EEULER) {
            //1st stage
            CalcDudt(nx, ndegr, nvar, dx, u, Dmatrix, Dradau, dudt);
            for (int i=0; i<nu; i++){
                //u_tmp[i] += dt * dudt[i];
                u[i] += dt * dudt[i];
       
                if (isnan(u[i])){
                    throw overflow_error("Kaboom!\n");
                }
            }
            if (nvar ==3) {
            LimitSolution(nx, ndegr, nvar, u);
            }
        }

        if (timestepping == TVDRK3) {
            veccopy(u_tmp, u, nu);
            //1st stage
            CalcDudt(nx, ndegr, nvar, dx, u, Dmatrix, Dradau, dudt);
            for (int i=0; i<nu; i++){
                u_tmp[i] += dt * dudt[i];
        
                if (isnan(u[i])){
                    throw overflow_error("Kaboom!\n");
                }
            }if (nvar ==3) {
            LimitSolution(nx, ndegr, nvar, u_tmp);
            }
            //2nd stage
            CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
            for (int i=0; i<nu; i++){
                u_tmp[i] = 0.75*u[i] + 0.25*( u_tmp[i] + dt*dudt[i]);
            }
            if (nvar ==3) {
            LimitSolution(nx, ndegr, nvar, u_tmp);
            }
            //3rd stage
            CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
            for (int i=0; i<nu; i++){
                u[i] = (1.0/3.0)*u[i] + (2.0/3.0)*(u_tmp[i] + dt*dudt[i]);
            }
            if (nvar ==3) {
            LimitSolution(nx, ndegr, nvar, u);
            }
        }

        if (timestepping == SEMIIMP) {
            // Do regular explicit euler, but with implicit discontinuous flux
            CalcDudt(nx, ndegr, nvar, dx, u, Dmatrix, Dradau, dudt);
            //
            // Implicit Solve
            SemiImplicitSolve(nx, ndegr, nvar, dx, u, Dmatrix, dt, dudt);
        }

        if (timestepping == SITVDRK3) {
            veccopy(u_tmp, u, nu);
            //1st stage
            CalcDudt(nx, ndegr, nvar, dx, u, Dmatrix, Dradau, dudt);
            // Implicit Solve
            SemiImplicitSolve(nx, ndegr, nvar, dx, u_tmp, Dmatrix, dt, dudt);
            LimitSolution(nx, ndegr, nvar, u_tmp);

            //2nd stage
            CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
            SemiImplicitSolve(nx, ndegr, nvar, dx, u_tmp, Dmatrix, dt, dudt);
            for (int i=0; i<nu; i++){
                u_tmp[i] = 0.75*u[i] + 0.25*(u_tmp[i]);
            }
            LimitSolution(nx, ndegr, nvar, u_tmp);
         
            //3rd stage
            CalcDudt(nx, ndegr, nvar, dx, u_tmp, Dmatrix, Dradau, dudt);
            SemiImplicitSolve(nx, ndegr, nvar, dx, u_tmp, Dmatrix, dt, dudt);
            for (int i=0; i<nu; i++){
                u[i] = (1.0/3.0)*u[i] + (2.0/3.0)*(u_tmp[i]);
            }
            LimitSolution(nx, ndegr, nvar, u);
        }

        if (timestepping == -999) { //abandoned for now
            // butcher tableau for the implicit method (kennedy DIRK methods for ode a review, 2016, p73)
            // 0  |  0         0    0
            // 2g |  g         g    0
            // c3 |  c3-a32-g  a32  g
            //       1-b2-b3   b2   b3
            // What to do for explicit portion?? - need to adjust coeffs to maintain consistency or something like that
            //
            double g,c3,a32,b1,b2,b3,t2,t3;
            g = (3.0 + sqrt(3.0)) / 6.0;
            c3 = 1.0; // 1.0=stiffly accurate, 
            // other options: internal stability:(3+sqrt(3))/3 - sqrt((2+sqrt(3))/3),  min error:(15+7sqrt(3))/(12(2+sqrt(3)))
            a32 = c3*(c3-2*g)/(4*g);
            b2 = (-2.0+3*c3)/(12*(c3-2*g)*g);
            b3 = (1.0-3*g)/(3*c3*(c3-2*g));
            b1 = 1.0 - b2 - b3;
            //
            // Algorithm Sketch fuly implicit ODE
            // F1 = F[t0        , u0                                    ]
            // F2 = F[t0 + 2g*dt, u0+dt*g*(F1 + F2)                     ]
            // F3 = F[t0 + c3*dt, u0+dt*((c3-a32-g)*F1 + a32*F2 + g*F3) ]
            // uout = u + dt*(b1*F1 + b2*F2 + b3*F3);
            //
            // Algorithm Sketch PDE with linear expansion of implicit terms
            // FN = delta_u/(delta_t*??) <- cN or bN?
            // F1 = F[t0        , u0                                    ]
            // F2 = F[t0 + 2g*dt, u0+dt*g*(2*F1 + dF/du|1*delta_u|2)    ]
            // F3 = F[t0 + c3*dt, u0+dt*((c3-a32-g)*F1 + a32*dFdu|1*delta_u|3+2 + g*dFdu|2*delta_u|3) ]
            // uout = u + dt*(b1*F1 + b2*F2 + b3*F3);
            //
            
        }
        
        if (nvar ==3) {
          if (iter % 1 == 0){printf("iter:%10d\tdt:%7.2e\t%7.2f%% Complete\n",iter, dt, 100.0*soltime/tmax);}
        }else{
          if (iter % 1 == 0){printf("iter:%10d\t%7.2f%% Complete\n",iter, 100.0*soltime/tmax);}
        }
        if (soltime >= tmax) {break;}
    //}

    //Printout Solution
    FILE* fout = fopen("waveout.tec", "w");
    if (nvar == 3) {
      fprintf(fout, "x\trho\trhou\trhoe\tp\tc\tM\n");
    } else {
      fprintf(fout, "x\tu\tu0\n");
    }

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
                getPrimativesPN(gam, &u[iu3(i, j, 0, ndegr)], &rho, &v, &p, &c, &M);

                fprintf(fout, "%f\t%f\t%f\n", p, c, M);
            }
        }
    }

    fclose(fout);

    }//

    printf("iter=%d\tdt=%f\n", niter, dt);

    //Printout Solution
    FILE* fout = fopen("waveout.tec", "w");
    if (nvar == 3) {
      fprintf(fout, "x\trho\trhou\trhoe\tp\tc\tM\n");
    } else {
      fprintf(fout, "x\tu\tu0\n");
    }

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
                getPrimativesPN(gam, &u[iu3(i, j, 0, ndegr)], &rho, &v, &p, &c, &M);

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
