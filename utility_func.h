#ifndef     _UTILITY_FUNC_H_
#define     _UTILITY_FUNC_H_

void normalized_fractions(int n, double *x);
void normalized_fractions_dst(int n, double *x, double *x_dst);

double x2y(int n, double *MW, double *x, double *y);
void   y2x(int n, double *MW, double *x, double *y);

double calc_rho_from_V(int n, double *x, double *MW, double V);
double calc_CvIG_from_CpIG(int n, double *x, double *CpIG);
void   calc_rhoY(int n, double rho, double *y, double *rhoY);
double calc_ave_MW(int n, double *x, double *MW);

#endif
