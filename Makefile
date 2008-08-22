OBJECT	= gpe
OBJS	= constants.o derivs.o gpe.o ic.o io.o parameters.o \
	  solve.o variables.o
FC	= sunf95
FFLAGS	= -fast -xtarget=pentium4 -fsimple=0 -xarch=generic
LDFLAGS	= -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal \
	  -ldl -lnsl -lutil -lm -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
INCLUDE	= /usr/lib/openmpi/include
#-----------------------------------------------------------------------
%.o : %.f90
	$(FC) $(FFLAGS) -I$(INCLUDE) -c $*.f90

all:    $(OBJECT)

clean :
	rm -f $(OBJECT) *.o *.mod
#-----------------------------------------------------------------------
$(OBJECT): $(OBJS)
	$(FC) $(FFLAGS) -I$(INCLUDE) $(OBJS) $(LDFLAGS) -o $@

derivs.o: parameters.o
gpe.o: derivs.o ic.o io.o parameters.o solve.o variables.o
ic.o: constants.o parameters.o variables.o
io.o: ic.o parameters.o variables.o
solve.o: derivs.o ic.o parameters.o variables.o
variables.o: derivs.o parameters.o
