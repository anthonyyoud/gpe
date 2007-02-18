# $Id: makefile,v 1.7 2007-02-18 18:30:30 najy2 Exp $
#----------------------------------------------------------------------------

OUTDIR		= ./
OBJECTS		= parameters.o constants.o derivs.o variables.o ic.o io.o \
                  solve.o gpe.o
FFLAGS	        = -O2 -w95 -tpp7 -xW -unroll -vec_report0
#FFLAGS	        = -O2
#FFLAGS	        = -O0 -w95
#FFLAGS	        = -g -O0 -d0 -CA -CB -CS -CU -CV
#FFLAGS	        = -g -O0 -check
LINKFLAGS	=
#LINKFLAGS	= -i_dynamic
COMPILER	= mpif90
#COMPILER	= gfortran
LIBS            =
COMPFLAGS       = $(FFLAGS)
LIBS            = -L$(FFTWHOME)/lib -lsrfftw_mpi -lsfftw_mpi \
                                    -lsrfftw -lsfftw

#-----------------------------------------------------------------------
all:	$(OBJECTS)
	$(COMPILER) $(COMPFLAGS) $(LINKFLAGS)-o \
        $(OUTDIR)gpe $(OBJECTS) $(LIBS)

#-----------------------------------------------------------------------
parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) -c parameters.f90

#-----------------------------------------------------------------------
constants.o : constants.f90
	$(COMPILER) $(COMPFLAGS) -c constants.f90

#-----------------------------------------------------------------------
ic.o : ic.f90
	$(COMPILER) $(COMPFLAGS) -c ic.f90

#-----------------------------------------------------------------------
io.o : io.f90
	$(COMPILER) $(COMPFLAGS) -c io.f90

#-----------------------------------------------------------------------
variables.o : variables.f90
	$(COMPILER) $(COMPFLAGS) -c variables.f90

#-----------------------------------------------------------------------
derivs.o : derivs.f90
	$(COMPILER) $(COMPFLAGS) -c derivs.f90

#-----------------------------------------------------------------------
gpe.o : gpe.f90
	$(COMPILER) $(COMPFLAGS) -c gpe.f90
#-----------------------------------------------------------------------
solve.o : solve.f90
	$(COMPILER) $(COMPFLAGS) -c solve.f90
#-----------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.out *.pc *.pcl *.il
