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
/* to compute binary interaction parameters
   input: 
          1) T: temperature
          2) Pc vector: critical pressure
          3) Tc vector: critical temperature
          4) w vector: acentric factor
          5) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details         
   output:
          0) kij matrix: binary interaction parameters
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *T, *Pc, *Tc, *w, *kij;
    mwSize  rows, cols, n;
    
    int *ntk;
    double *tk;

    int i;

    int bKijSet = 1;
    
    /* check for the proper number of arguments */
    if(nrhs != 5)
      mexErrMsgTxt("5 inputs required.");
    if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
    /*Check that 1st input is scalar*/
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1)
      mexErrMsgTxt("1st input must be scalars.");

    rows=mxGetM(prhs[1]);
    cols=mxGetN(prhs[1]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("2nd input must be row or column vectors.");

    if(  mxGetM(prhs[2]) != rows || mxGetN(prhs[2]) != cols  
      || mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols  
      || mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols)
      mexErrMsgTxt("3rd to 5th inputs must have the same shape.");

    /* get pointers to the inputs */
    T  = mxGetPr(prhs[0]);
    Pc = mxGetPr(prhs[1]);
    Tc = mxGetPr(prhs[2]);
    w  = mxGetPr(prhs[3]);
    tk = mxGetPr(prhs[4]);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
    }

    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    kij = mxGetPr(plhs[0]);

    /* call the C subroutine */
    calculate_kij_(&bKijSet,T,&n,Pc,Tc,w,ntk,kij);

    free(ntk);
    return;
}
