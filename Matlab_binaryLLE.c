/* $Revision: 1.8.6.2 $ */
/*=========================================================
 * Example for illustrating how to use pass complex data 
 * from MATLAB to C and back again
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2006 The MathWorks, Inc.
 *=======================================================*/
#include "mex.h"

/*#include "PR_EoS.h" */

/* The gateway routine. */
/* to compute binary liquid liquid equilibrium
   input: 
          1) P: pressure
          2) T: temperature
          3) Pc vector: critical pressure
          4) Tc vector: critical temperature
          5) w vector: acentric factor
          6) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details         
   output:
          0) xa vector: molar concentration
          1) xb vector: molar concentration
          2) n_miscible integer: 1:single phase 0:two-phase
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *w, *xa, *xb;
    mwSize  rows, cols, n;
    
    int *ntk;
    double *tk;

    double *n_miscible;
    int i, nn_miscible;

    int bKijSet = 1;
    double  *kij_tmp;
    
    /* check for the proper number of arguments */
    if(nrhs != 6)
      mexErrMsgTxt("6 inputs required.");
    if(nlhs > 3)
      mexErrMsgTxt("Too many output arguments.");
    /*Check that 1st - 2nd inputs are scalars*/
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1  
      || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)
      mexErrMsgTxt("1st - 3rd inputs must be scalars.");

    rows=mxGetM(prhs[2]);
    cols=mxGetN(prhs[2]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("3rd input must be row or column vectors.");

    if(  mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols  
      || mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols  
      || mxGetM(prhs[5]) != rows || mxGetN(prhs[5]) != cols)
      mexErrMsgTxt("4th to 6th inputs must have the same shape.");

    /* get pointers to the inputs */
    P  = mxGetPr(prhs[0]); 
    T  = mxGetPr(prhs[1]);
    Pc = mxGetPr(prhs[2]);
    Tc = mxGetPr(prhs[3]);
    w  = mxGetPr(prhs[4]);
    tk = mxGetPr(prhs[5]);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
    }

    kij_tmp = (double*)malloc(sizeof(double)*n*n);

    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    xa = mxGetPr(plhs[0]);
    xb = mxGetPr(plhs[1]);
    n_miscible = mxGetPr(plhs[2]);

    /* call the C subroutine */
    calculate_kij_(&bKijSet,T,&n,Pc,Tc,w,ntk,kij_tmp);
    binarylle_(P,T,&n,Pc,Tc,w,ntk,xa,xb,&nn_miscible);

    *n_miscible = nn_miscible;

    free(kij_tmp);
    free(ntk);
    return;
}
