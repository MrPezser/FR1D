//
// Information Geometric Regularization -> form of shock-capturing
//  Cao, Schafer arXive 2023/2024
//
// Trying this guy out - Sailor 09/20/2025
#include <cstdio>
//
#include "igr.h"
#include "indexing.h"
//


void IGRRHS(int nelem, int ndegr, int nvar, double alpha, double* u, double* du, double* igr_rhs) {
    // Calculates the right hand side / residual of the IGR function
    // a*(tr^2([Du]) + tr([Du]^2))
    //
    for (int ielem=0; ielem<nelem; ielem++){
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            //
            int ind = iu3(ielem, jdegr, 0, ndegr);
            // Chain rule d(rho u) / dx = rho dudx + u drhodx
            //            dudx = (drhoudx-udrhodx)/rho
            double dveldx = (du[ind+1] - u[1]*du[0]/u[0]) / u[ind];
            //
            igr_rhs[iup(ielem, jdegr, ndegr)] = alpha * 2.0 * dveldx * dveldx;
        }
    }
}
//
//
//
//void InitializeJacOps(){
    // maybe there's a way to put stuff in here idk.
//}
//
//
//
void UpdateJacobian(int nelem, int ndegr, int nvar, double dx, double alpha, double *sigma, double* u,const double* Dmatrix, \
                    double *jac){
    //
    // Calculates the IGR Jacobian
    //   Currently just doing local-to-element communication
    //
    //   sigma = current value of the IGR entropic pressure
    //
    int npoin = nelem*ndegr;
    int njac = 3*ndegr;
    double tmp[npoin*njac]; // sparce reresentation of block diagonal matrix
    //
    for (int i=0;i<npoin*njac;i++){
        jac[i] = 0.0;
        tmp[i] = 0.0;
    }
    // Loop through elements
    for (int ielem=0; ielem<nelem; ielem++){
        // Solving local problem on element
        //
        // Entries for finding gradient of sigma for each node on each element
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            //Use the Lagrange polynomial derivative matrix to perform the MatVec
            for (int knode = 0; knode < ndegr; knode++) {
                //
                // Scaling factor
                double k = 2.0 * alpha / (dx * u[iu3(ielem,knode,0,ndegr)]);
                //
                int ii = iup(ielem,jdegr,ndegr);
                //int ij = iup(ielem,knode,ndegr);
                tmp[iup(ii,knode+ndegr,njac)] += k * Dmatrix[iup(jdegr, knode, ndegr)];
            }
        }
        //
        // Entries for dependance on adjacent nodes
        // --Approach, Interface sig = 0.5 * (s_left + s_right)
        //          ds += g * (s_comm - s_int)
        //          coeff =  g * (0.5*(s_I + s_E) - s_I)
        //          coeff =  g * (0.5*s_E - 0.5*s_I)
        //          _I =   i_0   or   i_n-1
        //          _E = i-1_n-1 or i+1_0
        /*
        for (int idegr=0; idegr<ndegr; idegr++){
            // idegr - affected DoF
            int ii = iup(ielem,idegr,ndegr);
            //
            // Left Interface
            double g = Dradau[idegr];
            if (ielem>0) {
            // s_E = s(ielem-1,ndegr-1)
            tmp[iup(ii,ndegr-1,njac)] =  0.5 * g
            // s_I = s(ielem,0)
            tmp[iup(ii,ndegr,njac)]   = -0.5 * g
            }
            //
            // Right Interface
            g = -Dradau[ndegr-1-idegr];
            if (ielem<nelem-1) {
            // s_E = s(ielem+1,0)
            tmp[iup(ii,2*ndegr,njac)] =  0.5 * g
            // s_I = s(ielem,ndegr-1)
            tmp[iup(ii,ndegr+ndegr-1,njac)]   = -0.5 * g
            }
        }
         */
        //
        // Now find the divergence of our scaled gradient
        // I.E. need to multiply the div matrix by the matrix we currently have
        for (int row=0; row<ndegr; row++){ //   i - ish
            double sxi = 0.0;
            //Use the Lagrange polynomial derivative matrix to perform MatMat
            for (int col = 0; col < ndegr; col++) { //    j - ish
                //
                int ii = iup(ielem,row,ndegr);
                //int jj = iup(ielem,col,ndegr);
                //
                for (int k=0; k<ndegr; k++){
                    //
                    int kk = iup(ielem,k,ndegr);
                    //
                    double bkj = tmp[iup(kk,col+ndegr,njac)];
                    double aik = 2.0 * Dmatrix[iup(row, k, ndegr)] / dx;
                    jac[iup(ii,col+ndegr,njac)] -= aik * bkj;
                    // Really hope that's right
                }
            }// column loop
        }// row loop
        //
        // Finally, add the diagonal component, I/rho
        for (int jdegr=0; jdegr<ndegr; jdegr++){
            int ii = iup(ielem,jdegr,ndegr);
            jac[iup(ii,jdegr+ndegr,njac)] += 1.0 / u[iu3(ielem,jdegr,0,ndegr)];
        }
        //
    }// Element loop
}
//
//
//
void GaussSeidelIteration(int nelem, int ndegr, int ngs, double *rhs, double* jac, double* sigma){
    // Perform ngs number of gauss sidel iterations
    // recall the sparce representation of a block diagonal matrix [npoin] X [ndegr]
    int njac = 3*ndegr;
    //
    for (int gsiter=0; gsiter<ngs; gsiter++){
        // Do the GS solve
        for (int ielem=0; ielem<nelem; ielem++){
            for (int idegr=0; idegr<ndegr; idegr++){
                int ii = iup(ielem,idegr,ndegr);
                //
                double tmp = rhs[ii];
                // Contribution from current element
                for (int jdegr=0; jdegr<ndegr; jdegr++){
                    //
                    if (jdegr==idegr) {continue;} // i =/= j
                    int ij = iup(ielem,jdegr,ndegr);
                    tmp -= jac[iup(ii,jdegr+ndegr,njac)] * sigma[ij];
                    //
                }
                /*
                // Contribution from i-1 element
                if (ielem > 0) {
                for (int jdegr=0; jdegr<ndegr; jdegr++){
                    //
                    int ij = iup(ielem-1,jdegr,ndegr);
                    tmp -= jac[iup(ii,jdegr,njac)] * sigma[ij];
                    //
                }
                }
                // Contribution from i+1 element
                if (ielem < nelem-1) {
                for (int jdegr=0; jdegr<ndegr; jdegr++){
                    //
                    int ij = iup(ielem+1,jdegr,ndegr);
                    tmp -= jac[iup(ii,jdegr+2*ndegr,njac)] * sigma[ij];
                    //
                }
                }
                */
                //
                sigma[ii] = tmp / jac[iup(ii,idegr+ndegr,njac)];
            }// degr loop
        }// elem loop

    }// Gauss Seidel Iter Loop
}
//
//
//
void CalculateIGRSigma(int nelem, int ndegr, int nvar, double dx, double alpha,\
                       const double *Dmatrix, const double *Dradau, double* u, \
                       double* du, double* igr_sigma){
    // Calulate the entropic pressure needed to regularize the solution
    // 
    int num_fp_iter = 3;
    int npoin = nelem*ndegr;
    double igr_rhs[npoin];
    //
    // Step 1, calculate the RHS vector
    //printf("Calculating IGR RHS\n");
    IGRRHS(nelem, ndegr, nvar, alpha, u, du, igr_rhs);
    // 
    // Steb 2 Update the Jacobian Matrix
    //printf("Updating IGR Jacobian\n");
    double jac[npoin*ndegr*3];
    UpdateJacobian(nelem, ndegr, nvar, dx, alpha, igr_sigma, u, Dmatrix, jac);
    //
    // Step 3 Perform Fixed-Point Iterations to Approximately Solve Equation
    //printf("Performing IGR Iterations\n");
    GaussSeidelIteration(nelem, ndegr, num_fp_iter, igr_rhs, jac, igr_sigma);
    //
    return;
    FILE* fdebug = fopen("debug.log", "w");
    fprintf(fdebug,"Debug Log rhs:\n");
    for (int i=0; i<npoin; i++){
        fprintf(fdebug,"%e\n",igr_rhs[i]);
    }
    fprintf(fdebug,"Debug Log jac:\n");
    for (int i=0; i<npoin; i++){
        for (int j=0; j<ndegr*3; j++){
            fprintf(fdebug,"%12.4e",jac[iup(i,j,ndegr)]);
        }
        fprintf(fdebug,"\n");
    }
    fprintf(fdebug,"Debug Log sigma:\n");
    for (int i=0; i<npoin; i++){
        fprintf(fdebug,"%e\n",igr_sigma[i]);
    }
    fclose(fdebug);
    

}







