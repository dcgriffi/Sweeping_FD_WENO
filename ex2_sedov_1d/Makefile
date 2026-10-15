F90=gfortran
MPIF90=mpif90
F90FLAGS= -O2
a.out : main.f90 weno5.o
	$(F90) -o a $(F90FLAGS) main.f90 weno5.o		
weno5.o : weno5.f90
	$(F90) -c $(F90FLAGS) weno5.f90
clean:
	rm *.o

