#include "initial.h"
#include "SpatialDiscretization.h"
#include <cmath>

double Initialize(double x){
    ///This defines the intial state of the solution space


    //Gaussian bump and step combo

    if (x < 0.6) {
        double beta = 0.005;
        return 1.0 + exp(-(x-0.3)*(x-0.3) / beta);
    } else {
        if (x < 0.8) {
            return 2.0;
        } else {
            return 1.0;
        }
    }
    return 0.0;


    //return 2.0 + sin(2.0*M_PI*x);

    //return 1.0 + exp(-40*(x-0.5)*(x-0.5));

    //if (x<0.5){
    //    return 1.0;
    //} else {
    //    return 0.0;
    //}
}

void InitializeEuler(double x, double* u){
    double rho = 1.0;
    double v = 1.0;
    double p = 1.0;

    rho = Initialize(x);

    
    if (x < 0.5){
        rho = 1.0;
        v = 0.01;
        p = 1.0;
    } else {
        rho = 0.125;
        v = 0.01;
        p = 0.1;
    }
    


    u[0] = rho;                             //rho
    u[1] = rho * v;                         //rho V
    u[2] = 0.5*rho*v*v + (p/(GAM-1.0));       //rho e
}
