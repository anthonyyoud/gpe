!Copyright 2011 Anthony Youd/Newcastle University
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

module constants
  ! Constants for use in FFTW routines.  Original comments below:
  ! This file contains PARAMETER statements for various constants
  ! that can be passed to FFTW routines.  You should include
  ! this file in any FORTRAN program that calls the fftw_f77
  ! routines (either directly or with an #include statement
  ! if you use the C preprocessor).
  implicit none
  save
  
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1

  integer, parameter :: FFTW_REAL_TO_COMPLEX = -1 
  integer, parameter :: FFTW_COMPLEX_TO_REAL = 1

  integer, parameter :: FFTW_ESTIMATE = 0
  integer, parameter :: FFTW_MEASURE = 1

  integer, parameter :: FFTW_OUT_OF_PLACE = 0
  integer, parameter :: FFTW_IN_PLACE = 8
  integer, parameter :: FFTW_USE_WISDOM = 16

  integer, parameter :: FFTW_THREADSAFE = 128

  ! Constants for the MPI wrappers:
  integer, parameter :: FFTW_TRANSPOSED_ORDER = 1
  integer, parameter :: FFTW_NORMAL_ORDER = 0
  integer, parameter :: FFTW_SCRAMBLED_INPUT = 8192
  integer, parameter :: FFTW_SCRAMBLED_OUTPUT = 16384

end module constants
