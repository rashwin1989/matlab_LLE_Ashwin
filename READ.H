printf("Begin reading input files\n"); 

// I. READ INPUT FILE
int n_species_db, iread;
char c[50];
FILE *f;

int n_pseudospecies_db, n_purespecies_db;
f = fopen("species.dat+", "r");
iread = fscanf(f, "%d %d", &n_pseudospecies_db, &n_purespecies_db);
n_species_db = n_pseudospecies_db + n_purespecies_db;
fclose(f);

f = fopen("input.dat+", "r");

iread=fscanf(f, "%s %d", c, &n);  // num of species

// II. ALLOCATE MEMORY 
// 1. Inputs
_NEW_(idx_species, int, n);
_NEW_(idx_groups, int, n);    
_NEW_(Y01, double, n);
_NEW_(Y02, double, n);

iread=fscanf(f, "%s", c);         // index of species
for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_species[i]);

iread=fscanf(f, "%s", c);         // index of PPR78 groups
for (i=0; i<n; i++) iread=fscanf(f, "%d", &idx_groups[i]);

iread=fscanf(f, "%s", c);         // left hand side phase mass fractions
for (i=0; i<n; i++) iread=fscanf(f, "%lf", &Y01[i]);
iread=fscanf(f, "%s", c);         // right hand side phase mass fractions
for (i=0; i<n; i++) iread=fscanf(f, "%lf", &Y02[i]);

iread=fscanf(f, "%s %lf", c, &m012); // initial oil water mass ratio
  
iread=fscanf(f, "%s %lf",c, &P );  // uniform pressure (Pa)
iread=fscanf(f, "%s %lf",c, &T0);  // initial oil temperature
iread=fscanf(f, "%s %lf",c, &T1);  // initial water temperature  

fclose(f);

// II. ALLOCATE MEMORY
// 2. species
_NEW_(tk, int, n);
_NEW_(Pc, double, n);
_NEW_(Tc, double, n);
_NEW_(Vc, double, n);
_NEW_(w,  double, n);
_NEW_(MW, double, n);
_NEW_(dm, double, n);
_NEW_(k,  double, n);
_NEW_(H8, double, n);
_NEW_(Tb, double, n);
_NEW_(SG, double, n);
// 3. Working variables
_NEW_(x01, double, n);
_NEW_(x02, double, n);
_NEW_(x1, double, n);
_NEW_(x2, double, n);
_NEW_(Y1, double, n);
_NEW_(Y2, double, n);
_NEW_(kij_tmp, double, n*n);
_NEW_(XC0, double, n);
_NEW_(XC, double, n);
  
// III. READ SPECIES

f = fopen("petro.dat+","r");

for (j=0;j<n;j++) {
    for (i=0;i<n_species_db;i++) {
        double tmp;
        double xPNA[3];

        //[x' y' Tb' Tc' Pc' w' MW' Vc' SG' H_8' xP' xN' xA'];
        iread=fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &tmp, &tmp, &Tb[j], &Tc[j], &Pc[j], &w[j], &MW[j], &Vc[j], 
        &SG[j], &H8[j], &xPNA[0], &xPNA[1], &xPNA[2]);

        Pc[j] *= 1e5;  // bar --> Pa

        if (idx_groups[j]>=0) {
            // the first condition
            // we use PPR78 to calculate BIPs
            tk[j] = -(n_pseudospecies_db*idx_groups[j]+i+1);
        }
        else {
            // in the following
            // Guang's BIPs are used
            // NOTE that there are only 3 species
            // water, toluene, n-decane
            tk[j] = -idx_groups[j];
        }

        //x1[j] *= xPNA[ idx_groups[j] ];

        dm[j] = 0.0; 
        k [j] = 0.0;

        if (i==idx_species[j]) {
            // SPECIAL DATA POINTS BUILT IN THE CODE
            // 1. water
            if (fabs(Tc[j]-647)<0.5 && fabs(Pc[j]-220.6e5)<0.5e5) {
                Vc[j] = 57.1; // cm^3/mol
                dm[j] = 1.85; // Debye
                k [j] = 0.076;

            }
        
            if (H8[j]<0)
                // 2. toluene
                if (fabs(Tc[j]-591.8)<0.5 && fabs(Pc[j]-41.1e5)<0.5e5) {
                    Vc[j] = 316.0;
                    dm[j] = 0.36; // Debye

                    // 3. n-decane
                }else if (fabs(Tc[j]-617.7)<0.5 && fabs(Pc[j]-21.1e5)<0.5e5) {
                    Vc[j] = 600.0;

                    // 4. n-C30
                }else if (fabs(Tc[j]-844.0)<0.5 && fabs(Pc[j]- 8.0e5)<0.5e5) {
                    Vc[j] = 1805.0;

                    // 5. n-C50
                }else if (fabs(Tc[j]-1073.6)<0.5&& fabs(Pc[j]-3.51e5)<0.5e5) {
                    Vc[j] = 2999;

                    // 6. Benzene-C10
                }else if (fabs(Tc[j]-753.0)<0.5 && fabs(Pc[j]-17.7e5)<0.5e5) {
                    Vc[j] = 813;

                    // 7. Benzene-C30
                }else if (fabs(Tc[j]-965.6)<0.5 && fabs(Pc[j]- 5.3e5)<0.5e5) {
                    Vc[j] = 1943.5;

                    // 8. Naphthalene-C10
                }else if (fabs(Tc[j]-859.0)<0.5 && fabs(Pc[j]-15.8e5)<0.5e5) {
                    Vc[j] = 1070;

                    // 9. Naphthalene-C12
                }else if (fabs(Tc[j]-854.6)<0.5 && fabs(Pc[j]-13.0e5)<0.5e5) {
                    Vc[j] = 1081.5;

                    // 10. Benzene
                }else if (fabs(Tc[j]-562.0)<0.5 && fabs(Pc[j]-49.0e5)<0.5e5) {
                    Vc[j] = 256.0;

                    // 11. o-xylene
                }else if (fabs(Tc[j]-630.3)<0.5 && fabs(Pc[j]-37.3e5)<0.5e5) {
                    Vc[j] = 370.0;

                    // 12. p-xylene
                }else if (fabs(Tc[j]-616.2)<0.5 && fabs(Pc[j]-35.1e5)<0.5e5) {
                    Vc[j] = 378.0;

                    // 13. 1,3,5-trimethylbenzene
                }else if (fabs(Tc[j]-637.3)<0.5 && fabs(Pc[j]-31.3e5)<0.5e5) {
                    Vc[j] = 433.0;

                    // 14. naphthalene
                }else if (fabs(Tc[j]-748.4)<0.5 && fabs(Pc[j]-40.5e5)<0.5e5) {
                    Vc[j] = 407.0;

                    // 15. 1-methylnaphthalene
                }else if (fabs(Tc[j]-772.0)<0.5 && fabs(Pc[j]-36.0e5)<0.5e5) {
                    Vc[j] = 465;

                    // 16. anthracene
                }else if (fabs(Tc[j]-873.0)<0.5 && fabs(Pc[j]-29.0e5)<0.5e5) {
                    Vc[j] = 554;

                    // 17. 1,2-diphenylethane
                }else if (fabs(Tc[j]-780.0)<0.5 && fabs(Pc[j]-26.5e5)<0.5e5) {
                    Vc[j] = 616;

                    // 18. pyrene
                }else if (fabs(Tc[j]-936.0)<0.5 && fabs(Pc[j]-26.1e5)<0.5e5) {
                    Vc[j] = 660;

                    // 19. n-C16
                }else if (fabs(Tc[j]-723.0)<0.5 && fabs(Pc[j]-14.0e5)<0.5e5) {
                    Vc[j] = 969.2;

                    // 20. trans-decalin
                }else if (fabs(Tc[j]-687.0)<0.5 && fabs(Pc[j]-32.0e5)<0.5e5) {
                    Vc[j] = 480.0;

                    // 21. butylbenzene
                }else if (fabs(Tc[j]-660.5)<0.5 && fabs(Pc[j]-28.9e5)<0.5e5) {
                    Vc[j] = 497.0;

                    // 22. hexylbenzene
                }else if (fabs(Tc[j]-698.0)<0.5 && fabs(Pc[j]-23.8e5)<0.5e5) {
                    Vc[j] = 593.0;
                }

            // WHEN EQUALS, WE FOUND THE DATA 
            // SO QUIT THE READING FOR SPECIES I
            break;
        }
    }
    rewind(f);
}
fclose(f);

printf("Done reading input files\n\n");

