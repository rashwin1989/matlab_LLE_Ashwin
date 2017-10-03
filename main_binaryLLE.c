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
  
    FILE *fo = fopen("binaryLLE.dat", "w");
  
    #include "READ.H"        

    #include "PRINT_INPUT.H"

    for (T=T0; T<=T1;)
    {    
        printf("\n--------------------------------------------------------------------\n");
        printf("T = %lf\n\n", T);

        bKijSet = 1;
        calculate_kij_(&bKijSet, &T, &n, Pc, Tc, w, tk, kij_tmp);

        n_miscible = -1;
        binarylle_(&P, &T, &n, Pc, Tc, w, tk, x1, x2, &n_miscible); 

        x2y(n, MW, x1, Y1);
        x2y(n, MW, x2, Y2);        

        printf("T %lf\n", T);        
        printf("1-phase? %d\n", n_miscible);        
        printf("x1    ");
        for(j=0;j<n;j++) printf("%lf ", x1[j]);
        printf("\nx2    ");
        for(j=0;j<n;j++) printf("%lf ", x2[j]);        
        printf("\n\n--------------------------------------------------------------------\n");

        fprintf(fo, "%lf ", T);
        fprintf(fo, "%d ", n_miscible);  
        fprintf(fo, "%19.12le ", x1[1]);
        fprintf(fo, "%19.12le\n", x2[1]);
    
        T+=1;
        step+=1;
    }

    fclose(fo);  
    #include "CLEAN.H"
    return 0;
}
