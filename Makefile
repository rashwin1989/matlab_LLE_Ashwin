# COMPILER
CC = gcc                 
F77= gfortran            
FL = $(CC)#Linker        
NOFOR_MAIN =
CFLAGS_ADDON = -g 
FFLAGS_ADDON = -g
LIB_ADDON =

#  Intel C compiler
#CC = icc #-openmp 
#F77= ifort #-openmp 
#FL = $(F77)#Linker 
#NOFOR_MAIN = -nofor_main
#CFLAGS_ADDON = -static-intel
#FFLAGS_ADDON = -static-intel
#LIB_ADDON = -Wl,-Bstatic 

#  mixed
#CC = gcc
#F77= ifort
#FL = $(F77)#Linker 
#NOFOR_MAIN = -nofor_main
#CFLAGS_ADDON = -static
#FFLAGS_ADDON = -static-intel
#LIB_ADDON = -Wl,-Bstatic 

default: LLE all
all: LLE binaryLLE

r: c all #cp

#-DCUTOFF_MAX=1e6

CFLAGS = -O0 $(CFLAGS_ADDON)
FFLAGS = -O0 $(FFLAGS_ADDON) -fPIC 

LIB= -lm -lgfortran $(LIB_ADDON) \
     -L./gsl \
     -lgsl -lgslcblas

CPPFLAGS = -I./ \
           -I./gsl \
	   -I/usr/include/x86_64-linux-gnu/ 

SRC = vis_n_therm_new.o PR_EoS.o utility_func.o   

HEADER = CLEAN.H MACROS.H PRINT_INPUT.H READ.H DATA.H


LLE: $(SRC) main_LLE.o 
	$(FL) $(NOFOR_MAIN) $(CFLAGS) -o $@ $(SRC) main_LLE.o $(LIB)

binaryLLE: $(SRC) main_binaryLLE.o 
	$(FL) $(NOFOR_MAIN) $(CFLAGS) -o $@ $(SRC) main_binaryLLE.o $(LIB)

run: 
	./LLE

%.o: %.f
	$(F77) -c $(FFLAGS) -o $@ $<
%.o: %.f90
	$(F77) -c $(FFLAGS) -o $@ $<
%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

#-Wno-unused-result 

main_LLE.o: $(HEADER)
main_binaryLLE.o: $(HEADER)
utility_func.o: utility_func.c utility_func.h

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: 
	- $(RM) *.o
	- $(RM) LLE 
	- $(RM) binaryLLE

c: purge
cm: c default

