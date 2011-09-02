#Copyright 2011 Anthony Youd/Newcastle University
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# vim: set noexpandtab
#----------------------------------------------------------------------------

OBJECT	= gpe
OBJS	= constants.o derivs.o error.o gpe.o ic.o io.o mpi.o parameters.o \
        solve.o variables.o
#FC	= sunf95
#FFLAGS	= -fast -xmodel=medium #-fsimple=0
FC	= gfortran
FFLAGS	= -O3 -march=pentium4 -mfpmath=sse #-mcmodel=medium
#FFLAGS	= -O3 -march=core2 -mfpmath=sse #-mcmodel=medium
ifeq ($(precision), single)
  LDFFTW = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
else
  LDFFTW = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
endif
LDFLAGS	= -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal \
          -ldl -lnsl -lutil -lm $(LDFFTW)
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
gpe.o: derivs.o error.o ic.o io.o mpi.o parameters.o solve.o variables.o
ic.o: constants.o error.o parameters.o ic.in
io.o: error.o ic.o parameters.o variables.o
mpi.o: parameters.o
parameters.o: parameters.in
solve.o: derivs.o error.o ic.o parameters.o variables.o
variables.o: derivs.o ic.o parameters.o
