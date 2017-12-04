CXX = mpic++
F90 = mpif90
F77 = mpif77
CC  = mpicc

CFLAGS = -g -Wall -I../include
FFLAGS = -g
LINK   = -L../src -lBuffer -lmpi_cxx -lstdc++

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
%.o:%.f90
	$(F90) $(FFLAGS) -c $<

%.o:%.f
	$(F77) $(FFLAGS) -c $<
