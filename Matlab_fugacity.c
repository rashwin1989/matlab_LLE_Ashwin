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
/* to compute fugacity coefficient
   input: 
          1) P: pressure
          2) T: temperature
          3) Pc vector: critical pressure
          4) Tc vector: critical temperature
          5) w vector: acentric factor
          6) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details
          7) x vector: molar fractions
   output:
          0) ln(phi)
          1) phi
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *w, *x, *lnphi, *phi;
    mwSize  rows, cols, n;
    
    int *ntk;
    double *tk;

    int i;

    int bKijSet = 1;
    double  *kij_tmp;
    
    /* check for the proper number of arguments */
    if(nrhs != 7)
      mexErrMsgTxt("7 inputs required.");
    if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    /* check form of inputs */
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1  
      || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1  )
      mexErrMsgTxt("1st - 2nd inputs must be scalars.");

    rows=mxGetM(prhs[2]);
    cols=mxGetN(prhs[2]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("3rd input must be row or column vectors.");

    if(  mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols  
      || mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols  
      || mxGetM(prhs[5]) != rows || mxGetN(prhs[5]) != cols  
      || mxGetM(prhs[6]) != rows || mxGetN(prhs[6]) != cols )
      mexErrMsgTxt("4th to 7th inputs must have the same shape.");

    /* get pointers to the inputs */
    P  = mxGetPr(prhs[0]); 
    T  = mxGetPr(prhs[1]);
    Pc = mxGetPr(prhs[2]);
    Tc = mxGetPr(prhs[3]);
    w  = mxGetPr(prhs[4]);
    tk = mxGetPr(prhs[5]);
    x  = mxGetPr(prhs[6]);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
    }

    kij_tmp = (double*)malloc(sizeof(double)*n*n);

    /* create a new array and set the output pointer to it */
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    lnphi = mxGetPr(plhs[0]);
    phi   = mxGetPr(plhs[1]);
    
    /* call the C subroutine */
    calculate_kij_(&bKijSet,T,&n,Pc,Tc,w,ntk,kij_tmp);
    fugacities_only_(P,T,&n,Pc,Tc,w,ntk,x,lnphi,phi);

    free(ntk);
    free(kij_tmp);
    return;
}
