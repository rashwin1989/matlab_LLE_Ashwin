#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MACROS.H"
#include "utility_func.h"

extern int step;
extern FILE *f_output;

void normalized_fractions(int n, double *x){
  normalized_fractions_dst(n, x, x);
}
void normalized_fractions_dst(int n, double *x, double *x_dst)
{
  int i;
  double sum_x = 0.0;
  double xi;

  for (i=0; i<n; i++) {
    xi = x[i];

    if (xi < 0 ) {
      MY_PRINTF("Error: fractions x_%d %9.2le became negative!\n", i, xi);
      //exit(-1);
      x[i] = 0.0;
      xi = x[i];
    }

    sum_x += xi;
  }

  if (sum_x <=0 ) {
    MY_PRINTF("Error: fractions became negative or zero!\n");
    exit(-1);
  }

  sum_x = 1.0/sum_x;
  for (i=0; i<n; i++) x_dst[i] = x[i]*sum_x;
}

// because in normalized_fractions(), x has been checked for discrepancy
// x2y() & y2x() will not check fractions' positive property
double x2y(int n, double *MW, double *x, double *y){
  int i;
  double sum_MWx = 0.0, r1_sum;
  for (i=0; i<n; i++) sum_MWx += MW[i]*x[i];

  r1_sum = 1.0/sum_MWx;

  for (i=0; i<n; i++) y[i] = MW[i]*x[i] * r1_sum;

  return sum_MWx;
}
void y2x(int n, double *MW, double *x, double *y){
  int i;
  double sum_y_MW = 0.0;
  for (i=0; i<n; i++) sum_y_MW += y[i]/MW[i];

  sum_y_MW = 1.0/sum_y_MW;

  for (i=0; i<n; i++) x[i] = y[i]/MW[i] * sum_y_MW;
}

double calc_rho_from_V(int n, double *x, double *MW, double V){
  int i;
  double sum = 0.0;

  for (i=0;i<n;i++) {         //Mm = sum(M2.*x);
    sum += MW[i]*x[i];        //rho = 1e-3*Mm/Vm_real;
  }

  return 1e-3*sum/V;
}
double calc_ave_MW(int n, double *x, double *MW){
  int i;
  double sum = 0.0;
  for (i=0;i<n;i++) {         //Mm = sum(M2.*x);
    sum += MW[i]*x[i];        
  }
  return sum;
}

double calc_CvIG_from_CpIG(int n, double *x, double *CpIG){
  int i;
  double sum = 0.0;

  for (i=0;i<n;i++) {
    sum += CpIG[i]*x[i];        //CvIG(i)=sum(CpIG.*x)-8.314;
  }

  return sum - 8.3144621;
}

void   calc_rhoY(int n, double rho, double *y, double *rhoY){
  int i;
  for (i=0;i<n;i++) {
    rhoY[i] = rho*y[i];
  }
}
