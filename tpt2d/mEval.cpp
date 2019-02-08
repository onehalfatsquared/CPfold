/*********************************************************************
 * mEval.cpp
 *
 * This mex file computes the  potential for a system of N particles  
 * in 3D with a Morse interaction potential with given range and depth. 
 * Includes functionality for up to 3 particle types and variable
 * interaction strength between particles.   
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <math.h>
#include <algorithm>

/*
 * Function:  euDist
 * --------------------
 * computes euclidean distance b/w particles i and j. Also gives unit vector
 *
 *  particles: N by 3 array of positions. Z: stores unit vector. 
 *
 */

double euDist(double* particles, int i, int j, int N, double* Z){
    Z[0] = particles[i]-particles[j];
    Z[1] = particles[N+i]-particles[N+j];

    double R = sqrt(Z[0]*Z[0]+Z[1]*Z[1]);
    Z[0] /= R; Z[1] /= R;
    return R;
}

/*
 * Function:  morseEval
 * --------------------
 * computes and returns value of Morse potential for distance r, range rho
 * and well-depth E. 
 *
 */

double morseP(double r, double rho, double E){
    double Y = exp(-(rho) * ((r) - 1));
    return -2 * (rho) * (E) * ( Y*Y -Y);
}

// Main function 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *particles_m, *rho_m, *E_m, *N_m, *P_m, *u_m;
    double *particles, *rho, *E, *P, *u;
    const double rep = 0;
    int N;
    int i,j;

//associate inputs
    particles_m = mxDuplicateArray(prhs[0]);
    rho_m = mxDuplicateArray(prhs[1]);
    E_m = mxDuplicateArray(prhs[2]);
    N_m = mxDuplicateArray(prhs[3]);
    P_m = mxDuplicateArray(prhs[4]);

//associate pointers
    particles = mxGetPr(particles_m);
    rho = mxGetPr(rho_m);
    E = mxGetPr(E_m);
    N = round(*(mxGetPr(N_m)));
    P = mxGetPr(P_m);

//associate outputs
    u_m = plhs[0] = mxCreateDoubleScalar(1);
    u = (mxGetPr(u_m));

//compute the gradient
    double S = 0; 
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            if(j!=i){
                double* Z = new double[2];
                double r = euDist(particles, i, j, N, Z);
                int minimum = std::min(i,j);
                int maximum = std::max(i,j);
                if(P[N*maximum+minimum]==1){
                    double Eeff = E[N*maximum+minimum];
                    S += morseP(r, *rho, Eeff); 
                }
                else if(r<1){
                    S -= rep/r; 
                }
              delete []Z;
            }
        }
    }
    *u = S; 
}

