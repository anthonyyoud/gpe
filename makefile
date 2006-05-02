OUTDIR		= ./
OBJECTS		= parameters.o derivs.o ic.o variables.o io.o \
                  solve.o gpe.o
FFLAGS	        = -O2 -w95 -tpp7 -xW -unroll -vec_report0
#FFLAGS	        = -O2
#FFLAGS	        = -O0 -w95
#FFLAGS	        = -pg -d0 -CA -CB -CS -CU -CV
LINKFLAGS	=
#LINKFLAGS	= -i_dynamic
COMPILER	= ifort
#COMPILER	= gfortran
LDBLAS          = 
LDSCALA         = 
LDBLACS         = 

LIBS            = $(LDSCALA) $(LDBLACS) $(LDBLAS)
COMPFLAGS       = $(FFLAGS)
#-----------------------------------------------------------------------
all:	$(OBJECTS)
	$(COMPILER) $(COMPFLAGS) $(LINKFLAGS) -o \
        $(OUTDIR)gpe \
        $(OBJECTS) $(LIBS)

#-----------------------------------------------------------------------
parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) -c parameters.f90

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
