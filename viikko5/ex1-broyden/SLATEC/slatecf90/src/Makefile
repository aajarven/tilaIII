
###############################################################################
#                                                                             #
# Creates the SLATEC library from separate Fortran source files obtained from #
# http://people.scs.fsu.edu/~burkardt/f_src/slatec/slatec.html                #
# and split by f90split.                                                      #
#                                                                             #
# A. Kuronen, 2009                                                            #
#                                                                             #
###############################################################################

#########################
F90=gfortran
FFLAGS=-O3 #-ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time
AR=ar
#########################

#########################
#F90=ifort
#FFLAGS=-O2
#AR=xiar
#########################

LIB=libslatec.a
OBJ=$(patsubst %.f90,%.o,$(wildcard *.f90))

lib: $(OBJ)
	rm -f $(LIB)
	$(AR) qc $(LIB) $(OBJ)


# gfortran + optimization + DLAMC1 -> problems
#DLAMC1.o: DLAMC1.f90
#	$(F90) -O0 -c $<

%.o: %.f90
	$(F90) $(FFLAGS) -c $<


clean:
	rm -f *.o

cleanall:
	rm -f *.o $(LIB)
