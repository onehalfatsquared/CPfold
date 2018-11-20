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
#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *r_m, *rho_m, *E_m, *U_m;
    double *r, *rho, *E, *U;
    double Y;

//associate inputs
    r_m = mxDuplicateArray(prhs[0]);
    rho_m = mxDuplicateArray(prhs[1]);
    E_m = mxDuplicateArray(prhs[2]);


//associate outputs
    U_m = plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

//associate pointers
    r = mxGetPr(r_m);
    rho = mxGetPr(rho_m);
    E = mxGetPr(E_m);
    U = mxGetPr(U_m);

//do something
    Y = exp(-(*rho) * ((*r) - 1));
    *U = -2 * (*rho) * (*E) * ( Y*Y -Y);
    return;
}
