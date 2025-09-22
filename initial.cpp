#include "initial.h"
#include "SpatialDiscretization.h"
#include <cmath>

double sigmoid(double x0, double h, double w, double x){
    // returns the value of a sigmoid function at x 
    // centered around x0, with a height of h and width of w
    double xn = 1.5 * 12.0 * (x - x0) / w;
    double sign = 1.0 / (1.0 + exp(-xn));
    sign *= h;
    return sign;
}

double Initialize(double x){
    ///This defines the intial state of the solution space


    //Gaussian bump and step combo

    if (x < 0.6) { // 0.0 to 0.6 = gaussian bump
        double beta = 0.005;
        return 1.0 + exp(-(x-0.3)*(x-0.3) / beta);
    }
    if (x < 0.65) { // 0.6 to 0.65 = 1st sigmoid
        return 1.0 + sigmoid(0.625, 1, 0.05, x);
    } 
    if (x < 0.8) {
        return 2.0;
    }
    if (x < 0.85) {
        return 1.0 + sigmoid(-0.825, 1, 0.05, -x);
    }
    
    return 1.0;


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
    
    //if (x >= 0.6) {rho = 1.0;}

    double prog = sigmoid(0.5, 1.0, 0.1, x);
    rho = (1.0-prog)*1.0 + prog*0.125;
    p   = (1.0-prog)*1.0 + prog*0.1;
    v   = (1.0-prog)*0.0 + prog*0.0;
    /*
    if (x < 0.5){
        rho = 1.0;
        v = 0.01;
        p = 1.0;
    } else {
        rho = 0.125;
        v = 0.01;
        p = 0.1;
    }
    */


    u[0] = rho;                             //rho
    u[1] = rho * v;                         //rho V
    u[2] = 0.5*rho*v*v + (p/(GAM-1.0));       //rho e
}
