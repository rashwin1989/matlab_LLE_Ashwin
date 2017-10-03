#include "data.h"

//-------------------------------------------------------------------------------
// Data for species parameters
//-------------------------------------------------------------------------------
extern
int     n,             // num of species
       *idx_species,   // idx of species according to database petro.dat+
       *idx_groups,    // idx of PPR78 groups according to database groups.dat+
       *tk;            // type of BIP
extern
double *Pc, *Tc, *Vc,  // critical properties (Pa, K, cm^3/mol)
       *w,             // acentric factor
       *MW,            // molecular weight (g/mol)
       *dm,            // dipole moment (Debye)
       *k,             // association parameter
       *H8,            // molar enthalpy @ T=.8Tc (Btu/lb), which can be 0
       *Tb,            // boiling point (C)
       *SG;            // specific gravity
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// Data for thermodynamic variables
//-------------------------------------------------------------------------------
extern
double P, T, T0, T1,   // pressure, temperature
       *x1,            // mole fractions: oil phase
       *x2,            // mole fractions: water phase
       *x01,           // initial mole fractions: oil phase
       *x02,           // initial mole fractions: water phase
       *Y1,            // mass fractions: oil phase
       *Y2,            // mass fractions: water phase
       *Y01,           // initial mass fractions: oil phase
       *Y02,           // initial mass fractions: water phase
       c1,             // total oil phase mole fraction
       c01,            // total initial oil phase mole fraction
       m12,            // oil water mass ratio 
       m012,           // initial oil water mass ratio
       *kij_tmp,       // BIP temp variable
       MW1,
       MW2,
       MW01,
       MW02;
extern
int    bKijSet,        // to set kij's state
       n_miscible;     // >0: miscible; <=0: immiscible
//-------------------------------------------------------------------------------

