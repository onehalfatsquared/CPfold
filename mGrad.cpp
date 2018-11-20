/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 *
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <math.h>
#include <algorithm>


double euDist(double* particles, int i, int j, int N, double* Z){
    Z[0] = particles[i]-particles[j];
    Z[1] = particles[N+i]-particles[N+j];
    Z[2] = particles[2*N+i]-particles[2*N+j];

    double R = sqrt(Z[0]*Z[0]+Z[1]*Z[1]+Z[2]*Z[2]);
    Z[0] /= R; Z[1] /= R; Z[2] /= R;
    return R;
}

double morseEval(double r, double rho, double E){
    double Y = exp(-(rho) * ((r) - 1));
    return -2 * (rho) * (E) * ( Y*Y -Y);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *particles_m, *rho_m, *E_m, *N_m, *P_m, *g_m;
    double *particles, *rho, *E, *P, *g;
    const double rep = 3;
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
    g_m = plhs[0] = mxCreateDoubleMatrix(3*N,1,mxREAL);
    g = mxGetPr(g_m);

//do something
    for(i=0;i<N;i++){
        double* S = new double[3]; S[0] = S[1] = S[2] = 0;
        for(j=0;j<N;j++){
            if(j!=i){
                double* Z = new double[3];
                double r = euDist(particles, i, j, N, Z);
                printf("%.4f",r);
                int minimum = std::min(i,j);
                int maximum = std::max(i,j);
                if(P[N*maximum+minimum]==1){
                    for(int k=0; k<3;k++)
                        S[k] = S[k] + morseEval(r, *rho, *E) * Z[k];
                }
                else if(r<1){
                    for(int k=0; k<3;k++)
                        S[k] = S[k] - rep * Z[k] / (r*r);
                }
              delete []Z;
            }
        }
      for(int k=0; k<3;k++)
            g[3*i+k] = S[k];
      delete []S;
    }
}

