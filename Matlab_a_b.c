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
/* to compute a and b for all species using PR EoS
   input: 
          1) T: temperature
          2) Pc vector: critical pressure
          3) Tc vector: critical temperature
          4) w vector: acentric factor
   output:
          1) a: molecular interaction parameter
          2) b: excluded volume parameter    
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *T, *Pc, *Tc, *w, *coef_ab, *a, *b;
    mwSize  rows, cols, n;
       
    int i;
    
    /* check for the proper number of arguments */
    if(nrhs != 4)
      mexErrMsgTxt("4 inputs required.");
    if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1)     
      mexErrMsgTxt("1st - 2nd inputs must be scalars.");

    rows=mxGetM(prhs[1]);
    cols=mxGetN(prhs[1]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("2nd input must be row or column vectors.");

    if(  mxGetM(prhs[2]) != rows || mxGetN(prhs[2]) != cols  
      || mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols)  
      mexErrMsgTxt("3rd to 4th inputs must have the same shape.");

    /* get pointers to the inputs */
    T  = mxGetPr(prhs[0]);
    Pc = mxGetPr(prhs[1]);
    Tc = mxGetPr(prhs[2]);
    w  = mxGetPr(prhs[3]);

    coef_ab = (double*)malloc(sizeof(double)*n);

    for (i=0; i<n; i++) {
      coef_ab[i]=-1;
    }    

    /* create a new array and set the output pointer to it */
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    a = mxGetPr(plhs[0]);
    b = mxGetPr(plhs[1]);    

    /* call the C subroutine */
    calculate_a_b_(T,&n,Pc,Tc,w,coef_ab,a,b);

    free(coef_ab);
    return;
}
