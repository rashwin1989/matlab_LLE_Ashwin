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
/* to compute matrix of binary diffusion co-efficients
   input: 
          1) P: pressure
          2) T: temperature
          3) Pc vector: critical pressure
          4) Tc vector: critical temperature
          5) Vc vector: critical volume (cm^3/mol)    
          6) w vector: acentric factor
          7) tk vector: type of k, binary interaction parameters
                        check subroutine calculate_kij in PR_EoS.f90 for details
          8) MW vector: molecular weight (g/mol)          
          9) x vector: molar fractions
   output:
          0) Dij (length n*n stored in form idx = i + j*n)
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *Vc, *w, *MW, *x, *Dij;
    mwSize  rows, cols, n;

    int *ntk;
    double *tk, *coef_ab;

    int i;
    
    /* check for the proper number of arguments */
    if(nrhs != 9)
      mexErrMsgTxt("7 inputs required.");
    if(nlhs > 1)
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
      || mxGetM(prhs[6]) != rows || mxGetN(prhs[6]) != cols 
      || mxGetM(prhs[7]) != rows || mxGetN(prhs[6]) != cols 
      || mxGetM(prhs[8]) != rows || mxGetN(prhs[6]) != cols )
      mexErrMsgTxt("4th to 9th inputs must have the same shape.");

    /* get pointers to the inputs */
    P  = mxGetPr(prhs[0]); 
    T  = mxGetPr(prhs[1]);
    Pc = mxGetPr(prhs[2]);
    Tc = mxGetPr(prhs[3]);
    Vc = mxGetPr(prhs[4]);
    w  = mxGetPr(prhs[5]);
    tk = mxGetPr(prhs[6]);
    MW = mxGetPr(prhs[7]);
    x  = mxGetPr(prhs[8]);

    coef_ab = (double*)malloc(sizeof(double)*n);

    ntk= (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) {
      ntk[i]=tk[i];
      coef_ab[i]=-1;
    }

    /* create a new array and set the output pointer to it */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    Dij = mxGetPr(plhs[0]);
    
    /* call the C subroutine */
    new_tlsm_diffusion_krishna_model_(P,T,&n,Pc,Tc,Vc,w,ntk,coef_ab,MW,x,Dij);

    free(coef_ab);
    free(ntk);
    return;
}
