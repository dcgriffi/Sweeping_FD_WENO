#a.out : main.f90 weno5.o
#	mpif90 -o a -r8 -fast main.f90 weno5.o		
#weno5.o : weno5.f90
#	pgf90 -c -r8 -fast weno5.f90

# Compiler settings
F90      = gfortran
MPIF90   = mpif90
F90FLAGS = -O3 -fdefault-real-8 -g -ffree-line-length-none

# Targets
a.out : main.f90 weno5.o
	$(MPIF90) $(F90FLAGS) -o a main.f90 weno5.o

weno5.o : weno5.f90
	$(F90) $(F90FLAGS) -c weno5.f90

clean:
	rm -f *.o a.out
