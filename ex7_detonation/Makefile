#F90=ifort
#MPIF90=mpif90
#F90FLAGS=-r8 -O2
#MKLLINKS=-I$(MKLROOT)/include/intel64/lp64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#a.out : main.f90 weno5.o
#	$(MPIF90) -o a $(F90FLAGS) main.f90 weno5.o		
#weno5.o : weno5.f90
#	$(F90) -c $(F90FLAGS) weno5.f90
#clean:
#	rm *.o


# Compiler settings
F90      = gfortran
MPIF90   = mpif90
F90FLAGS = -O3 -fdefault-real-8  -ffree-line-length-none

# Targets
a.out : main.f90 weno5.o
	$(MPIF90) $(F90FLAGS) -o a main.f90 weno5.o

weno5.o : weno5.f90
	$(F90) $(F90FLAGS) -c weno5.f90

clean:
	rm -f *.o a.out
