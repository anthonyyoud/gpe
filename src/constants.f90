!
! Copyright (c) Anthony J. Youd/Newcastle University 2011
!
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
