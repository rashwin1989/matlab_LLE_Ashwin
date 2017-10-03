#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MACROS.H"
#include "PR_EoS.h"
#include "vis_n_therm_new.h"
#include "utility_func.h"

int main(int argc, char *argv[])
{  
    #include "DATA.H"
    int i, j, step;    
    step = 0;
  
    FILE *fo = fopen("LLE.dat", "w");
    FILE *fy1 = fopen("y1.dat", "w");
    FILE *fy2 = fopen("y2.dat", "w");
  
    #include "READ.H"    

    // normalize input mass fractions and calculate mole fractions  
    normalized_fractions(n, Y01);
    normalized_fractions(n, Y02);
    y2x(n, MW, x01, Y01);
    y2x(n, MW, x02, Y02);
    MW01 = calc_ave_MW(n, x01, MW);
    MW02 = calc_ave_MW(n, x02, MW);
    c01 = m012/MW01/(m012/MW01 + 1.0/MW02);

    for(j=0;j<n;j++) XC0[j] = c01*x01[j] + (1.0 - c01)*x02[j];

    #include "PRINT_INPUT.H"

    for (T=T0; T<=T1;)
    {    
        printf("\n--------------------------------------------------------------------\n");
        printf("T = %lf\n\n", T);

        bKijSet = 1;
        calculate_kij_(&bKijSet, &T, &n, Pc, Tc, w, tk, kij_tmp);

        n_miscible = -1;
        species_lle4_(&P, &T, &n, Pc, Tc, w, tk, x01, x02, &c01, &c1, x1, x2, &n_miscible); 

        x2y(n, MW, x1, Y1);
        x2y(n, MW, x2, Y2);
        MW1 = calc_ave_MW(n, x1, MW);
        MW2 = calc_ave_MW(n, x2, MW);
        m12 = MW1*c1/(MW1*c1 + MW2*(1.0 - c1));
        for(j=0;j<n;j++) XC[j] = c1*x1[j] + (1.0 - c1)*x2[j];

        printf("T %lf\n", T);
        printf("m1 %lf\n", m12);
        printf("c1 %lf\n", c1);
        printf("1-phase? %d\n", n_miscible);
        printf("Y1    ");
        for(j=0;j<n;j++) printf("%lf ", Y1[j]);
        printf("\nY2    ");
        for(j=0;j<n;j++) printf("%lf ", Y2[j]);
        printf("\nx1    ");
        for(j=0;j<n;j++) printf("%lf ", x1[j]);
        printf("\nx2    ");
        for(j=0;j<n;j++) printf("%lf ", x2[j]);
        printf("\nXC0    ");
        for(j=0;j<n;j++) printf("%lf ", XC0[j]);
        printf("\nXC    ");
        for(j=0;j<n;j++) printf("%lf ", XC[j]);
        printf("\n\n--------------------------------------------------------------------\n");

        fprintf(fo, "%lf ", T);
        fprintf(fo, "%lf ", m12);
        fprintf(fo, "%d\n", n_miscible);            

        for (j=0; j<n;  j++)
        {
            fprintf(fy1, "%19.12le ", Y1[j]);
            fprintf(fy2, "%19.12le ", Y2[j]);
        }
        fprintf(fy1, "\n");
        fprintf(fy2, "\n");
    
        T+=1;
        step+=1;
    }

    fclose(fo);  
    fclose(fy1);  
    fclose(fy2);  
    #include "CLEAN.H"
    return 0;
}
