//
// Created by Tsail on 6/16/2023.
//

#include <cmath>
#include <cstdio>
#include <cmath>
#include <random>
#include "basis.h"
#include "indexing.h"
#include "SailorMath.h"

using namespace std;

double EvalLobatto(int numPoints, double x){
    double p0, p1, p2, pkm2;
    p0 = 1.0;
    p1 = x;
    //This finds the value of the jth order lobatto polynomial at the given points (i.e. the value we want to =0)
    //method is that Lo_k(x) = P_legendre_k(x) - P_l_k-2(x)

    if (numPoints == 1){return 1.0;}

    for (int j = 2; j <= numPoints; ++j) {
        p2 = ((2.0 * j - 1.0) * x * p1 - (j - 1.0) * p0) / j;
        if (j == numPoints){
            pkm2 = p0;
        }

        p0 = p1;
        p1 = p2;
    }

    return p2-pkm2;
}

void GenerateLobattoPoints(int numPoints, double* points){
    int isign;
    double x;
    double dxi = 1e-8;
    double pointMinSep = 1e-5;
    double x0 = -1.0+pointMinSep;
    double tolerance = 1e-10;
    int maxIterations = 100;
    int outmaxIter = ceil(2.0 / pointMinSep);

    // Number of points should be at least 2
    if (numPoints < 2) {
        points[0] = 0.0;
        return;
    }

    // Compute Lobatto quadrature points
    points[0] = -1.0;
    points[numPoints - 1] = 1.0;

    if (EvalLobatto(numPoints, x0) > 0){
        isign = 1;
    } else {
        isign = -1;
    }

    for (int i = 1; i < numPoints-1; ++i) {


        for (int outIter = 0; outIter < outmaxIter; outIter++) {
            x = x0;
            double p, xitest, p2, dp;
            double p0 = EvalLobatto(numPoints, x);

            if ((p0>=0 && isign==-1)  ||  (p0<=0 && isign==1)) {
                // Newton Iterate to refine the root estimate of the lobatto fxn
                for (int iter = 0; iter < maxIterations; ++iter) {
                    xitest = (x + dxi);
                    p = EvalLobatto(numPoints, x);
                    p2 = EvalLobatto(numPoints, xitest);

                    dp = (p2 - p) / (dxi);
                    x -= p / dp;

                    // Check convergence
                    if (std::abs(p) < tolerance) {
                        break;
                    }
                    if (iter == maxIterations - 1) {
                        printf("Error finding Lobatto Point\n");
                    }
                }
                if (p0>0){
                    isign = 1;
                } else {
                    isign = -1;
                }
                break;
            } else {
                x0 += pointMinSep;
                continue;
            }

            /*
            //check if it is a duplicate point
            int iDupe = 0;

            for (int i2=0; i2<i; i2++){
                if ( fabs(x - points[i2]) < pointMinSep || x >= (1.0-pointMinSep) ){
                    iDupe = 1;
                }
            }

            if (iDupe == 1) {
                x0 += dxi;

                if (outIter == outmaxIter-1){
                    printf("Unlucky\n");
                }
                //we'll find it eventually
                continue;
            }
            break;
             */
        }

        // Save Lobatto quadrature point
        points[i] = x;
    }
}


double LagrangeDerivative(int i, int j, int ndegr, const double *x) {
    ///Modified from Thomas: https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
 //computes the derivative of L_j at a given x_i
    double z = x[i];
    double y = 0.0;
    for (int l=0; l<ndegr;l++) {
        if (l != j) {
            double k = 1 / (x[j] - x[l]);
            for (int m=0; m<ndegr; m++) {
                if((m != j) && (m != l)) {
                    k = k * (z - x[m]) / (x[j] - x[m]);
                }
            }
            y = y + k;
        }
    }
    return y;
}

void GenerateLagrangeDMatrix( int ndegr, const double* loPoints, double* Dmatrix){
//computes the derivative of L_j at a given x_i
    for (int i=0; i<ndegr; i++){
        for (int j=0; j<ndegr; j++){
            //cry
            Dmatrix[iup(i,j, ndegr)] = LagrangeDerivative(i,j,ndegr, loPoints);
        }
    }
}


void LegendreCoefficients(int ndegr, double* coeff){
    ///Find the coefficients of the Legendre Polynomial of degree ndegr

    for (int idegr=0; idegr<=ndegr; idegr++){
        double temp = 0.5*(double)(ndegr+idegr-1);
        coeff[idegr] = pow(2.0, (double)ndegr) * BinomialCoefficient((double)ndegr, idegr) * BinomialCoefficient(temp, ndegr);
    }

}

void GenerateRadauDerivatives(int ndegr, const double *x, double* Dradau){
    ///Finds the derivatives of the Right Radau polynomial of degree ndegr at points x
    //Radau_R = ((-1^k)/2) * (P_k - P_k-1), left Radau is just mirrored around x=0


    //Monotobasis functions
    for (int ipoin=0; ipoin<ndegr; ipoin++) {
        int gorder = ndegr;

        //Dradau[ipoin] = - 0.5 * M_PI_2 * cos(x[ipoin]*M_PI_2);
        //Dradau[ipoin+ndegr] = 0.5 * (- sin(x[ipoin]*M_PI_2) + 1.0);

        //Dradau[ipoin] = -0.5 * gorder * pow(0.5*(1-x[ipoin]), gorder-1); //derivatives
        //Dradau[ipoin+ndegr] = pow(0.5*(1-x[ipoin]), gorder);//function values

        //Dradau[ipoin] = 0.0;
        //Dradau[ipoin+ndegr] = 0.0;
    }
    //Dradau[0] = -2.0 / (x[1]-x[0]);
    //Dradau[0+ndegr] = 1.0;
    //return;

    if (ndegr==2){
        Dradau[0] = -20/8.0;//2.124 - 1.0/2.0;
        Dradau[1] =  12.0/8.0;//1.0/2.124 - 1.0/2.0;
        //Dradau[2] = 1.0;
        //Dradau[3] = 0.0;
    }
    return;

    //Get coefficients of the p_k and p_k-1 legendre polynomials
    auto* coeffpk   = (double*)malloc((ndegr+1)*sizeof(double));
    LegendreCoefficients(ndegr, coeffpk);

    auto* coeffpkm1 = (double*)malloc((ndegr+1)*sizeof(double));
    coeffpkm1[ndegr] = 0; // initialize the highest degree which will not be used
    LegendreCoefficients(ndegr-1, coeffpkm1);


    //Construct the Radau Polynomial
    double coeffRR[ndegr+1];
    for (int idegr=0; idegr<=ndegr; idegr++){

        coeffRR[idegr] = 0.5*(coeffpk[idegr] - coeffpkm1[idegr]);
        if (ndegr % 2) {
            coeffRR[idegr] *= -1.0;
        }
    }

    //Compute the analytical derivative using the power-rule
    double coeffDrad[ndegr];
    for (int idegr=0; idegr < ndegr; idegr++) {
        //also store the value of the radau polynomial
        Dradau[idegr+ndegr] = 0.0;
        for (int ipoin=0; ipoin<ndegr+1; ipoin++){
            Dradau[idegr+ndegr] += coeffRR[ipoin] * pow(x[idegr], (double)ipoin);  //degr and point are swapped on this
        }

        coeffDrad[idegr] = (idegr+1)*coeffRR[idegr+1];
    }

    //evaluate the derivative of the polynomial at the solution points
    for (int ipoin=0; ipoin<ndegr; ipoin++){
        Dradau[ipoin] = 0.0;

        for (int idegr=0; idegr<ndegr; idegr++){
            Dradau[ipoin] += coeffDrad[idegr] * pow(x[ipoin], (double)idegr);
        }


    }

    free(coeffpk);
    free(coeffpkm1);
}
