OBJECT	= gpe
OBJS	= constants.o derivs.o error.o gpe.o ic.o io.o parameters.o \
	  solve.o variables.o
FC	= sunf95
FFLAGS	= -fast #-fsimple=0
#FC	= gfortran
#FFLAGS	= -O3 -march=core2 -mfpmath=sse
LDFLAGS	= -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal \
	  -ldl -lnsl -lutil -lm -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
INCLUDE	= -I/usr/lib/openmpi/include
#-----------------------------------------------------------------------
%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $*.f90

all:    $(OBJECT)

clean :
	rm -f $(OBJECT) *.o *.mod
#-----------------------------------------------------------------------
$(OBJECT): $(OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) $(OBJS) $(LDFLAGS) -o $@

derivs.o: parameters.o
error.o: parameters.o
gpe.o: derivs.o error.o ic.o io.o parameters.o solve.o variables.o
ic.o: constants.o error.o parameters.o
io.o: error.o ic.o parameters.o variables.o
solve.o: derivs.o error.o ic.o parameters.o variables.o
variables.o: derivs.o ic.o parameters.o
