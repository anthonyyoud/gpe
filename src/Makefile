# vim: set noexpandtab
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
#
#----------------------------------------------------------------------------

include Makefile.inc

OBJECT	= gpe
OBJS	= derivs.o error.o gpe.o ic.o io.o mpi.o parameters.o solve.o variables.o
ifeq ($(precision), single)
  LDFFTW = -L$(FFTW_HOME)/lib -lfftw3f -lfftw3
  LD2DECOMP = -L$(DECOMP_HOME_SINGLE)/lib -l2decomp_fft
  DECOMP_INC = -I$(DECOMP_HOME_SINGLE)/include
else
  LDFFTW = -L$(FFTW_HOME)/lib -lfftw3
  LD2DECOMP = -L$(DECOMP_HOME_DOUBLE)/lib -l2decomp_fft
  DECOMP_INC = -I$(DECOMP_HOME_DOUBLE)/include
endif
LDFLAGS	= $(LD2DECOMP) $(LDFFTW) -lm
INCLUDE	= $(DECOMP_INC)
#-----------------------------------------------------------------------
%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $*.f90

all:    $(OBJECT)

clean :
	rm -f $(OBJECT) *.o *.mod
#-----------------------------------------------------------------------
$(OBJECT): $(OBJS)
	$(FC) $(OBJS) $(LDFLAGS) -o $@

derivs.o: parameters.o
error.o: parameters.o
gpe.o: derivs.o error.o ic.o io.o mpi.o parameters.o solve.o variables.o
ic.o: error.o parameters.o ic.in
io.o: error.o ic.o parameters.o variables.o
mpi.o: parameters.o
parameters.o: parameters.in
solve.o: derivs.o error.o ic.o parameters.o variables.o
variables.o: derivs.o ic.o parameters.o
