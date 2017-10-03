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
/* to compute viscosity and thermal conductivity
   input: 
          1) P: pressure
          2) T: temperature
          3) Pc vector: critical pressure
          4) Tc vector: critical temperature
          5) Vc vector: critical volume (cm^3/mol)
          6) w vector: acentric factor          
          7) MW vector: MW (g/mol)
          8) k vector: association factor
          9) dm vector: dipole moment (Debye)
          10) x vector: molar fractions          
          11) V: molar volume (m^3/mol)
          12) Tb vector: boiling point (degC)
          13) SG vector: specific gravity
          14) H8 vector: ideal gas H, Cp type and index
   output:
          1) mu: viscosity (Pa-s)
          2) lambda: thermal conductivity (W/m-K)
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double  *P, *T, *Pc, *Tc, *Vc, *w, *MW, *k, *dm, *x, *V, *Tb, *SG, *H8, *mu, *lambda;
    mwSize  rows, cols, n;
        
    double Cv_i, H_i, Cv_m;
    double R_gas = 8.3144621;
    int i;
    
    /* check for the proper number of arguments */
    if(nrhs != 14)
      mexErrMsgTxt("14 inputs required.");
    if(nlhs > 2)
      mexErrMsgTxt("Too many output arguments.");
    /*Check that 1st two inputs are scalars*/
    if(  mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 1  
      || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1  
      || mxGetM(prhs[10]) != 1 || mxGetN(prhs[10]) != 1  )
      mexErrMsgTxt("1st - 2nd and 11th inputs must be scalars.");

    rows=mxGetM(prhs[2]);
    cols=mxGetN(prhs[2]);
    n = rows>cols?rows:cols;
    if (n!=rows*cols || n==1) 
      mexErrMsgTxt("3rd input must be row or column vectors.");

    if(  mxGetM(prhs[3]) != rows || mxGetN(prhs[3]) != cols  
      || mxGetM(prhs[4]) != rows || mxGetN(prhs[4]) != cols  
      || mxGetM(prhs[5]) != rows || mxGetN(prhs[5]) != cols  
      || mxGetM(prhs[6]) != rows || mxGetN(prhs[6]) != cols 
      || mxGetM(prhs[7]) != rows || mxGetN(prhs[7]) != cols
      || mxGetM(prhs[8]) != rows || mxGetN(prhs[8]) != cols
      || mxGetM(prhs[9]) != rows || mxGetN(prhs[9]) != cols
      || mxGetM(prhs[11]) != rows || mxGetN(prhs[11]) != cols
      || mxGetM(prhs[12]) != rows || mxGetN(prhs[12]) != cols
      || mxGetM(prhs[13]) != rows || mxGetN(prhs[13]) != cols)
      mexErrMsgTxt("4th to 10th and 12th to 14th inputs must have the same shape.");

    /* get pointers to the inputs */
    P  = mxGetPr(prhs[0]); 
    T  = mxGetPr(prhs[1]);
    Pc = mxGetPr(prhs[2]);
    Tc = mxGetPr(prhs[3]);
    Vc = mxGetPr(prhs[4]);
    w  = mxGetPr(prhs[5]);
    MW = mxGetPr(prhs[6]);
    k  = mxGetPr(prhs[7]);
    dm  = mxGetPr(prhs[8]);
    x  = mxGetPr(prhs[9]);
    V  = mxGetPr(prhs[10]);    
    Tb  = mxGetPr(prhs[11]);    
    SG  = mxGetPr(prhs[12]);    
    H8  = mxGetPr(prhs[13]);    

    /* create a new array and set the output pointer to it */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mu = mxGetPr(plhs[0]);
    lambda = mxGetPr(plhs[1]);

    /* call the C subroutine */
    Cv_m = 0.0;
    for(i=0;i<n;i++)
    {
        ig_cp_h_(&Tb[i],&Tc[i],&SG[i],&H8[i],&MW[i],T,&H_i,&Cv_i);
        Cv_i -= R_gas;
        Cv_m += x[i]*Cv_i;
    }
    vis_n_cond_(P,T,&n,Pc,Tc,Vc,w,MW,k,dm,x,&Cv_m,V,lambda,mu);
    
    return;
}
