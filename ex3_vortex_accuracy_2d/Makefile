# a.out : main.f90 weno5.o
# 	ifort -o a -r8 -O3 -vec_report0 main.f90 weno5.o
# 		
# weno5.o : weno5.f90
# 	ifort -c -r8 -O3 -vec_report0 weno5.f90
#a.out : main.f90 weno5.o
#	pgf90 -o a -r8 -fast main.f90 weno5.o		
#weno5.o : weno5.f90
#	pgf90 -c -r8 -fast weno5.f90


F90=gfortran
MPIF90=mpif90
F90FLAGS= -O3 -fdefault-real-8
a.out : main.f90 weno5.o
	$(F90) -o a $(F90FLAGS) main.f90 weno5.o		
weno5.o : weno5.f90
	$(F90) -c $(F90FLAGS) weno5.f90
clean:
	rm *.o

