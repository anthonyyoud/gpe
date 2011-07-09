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

module mpi
  ! Routines to do with MPI.
  implicit none

  private
  public :: mpi_constants_precision
  
  contains

  subroutine mpi_constants_precision()
    ! Define new MPI constants depending on the selected precision.
    use parameters
    implicit none

    if (pr > 4) then
      gpe_mpi_real     = MPI_DOUBLE_PRECISION
      gpe_mpi_2real    = MPI_2DOUBLE_PRECISION
      gpe_mpi_complex  = MPI_DOUBLE_COMPLEX
      gpe_mpi_2complex = MPI_2DOUBLE_COMPLEX
    else
      gpe_mpi_real     = MPI_REAL
      gpe_mpi_2real    = MPI_2REAL
      gpe_mpi_complex  = MPI_COMPLEX
      gpe_mpi_2complex = MPI_2COMPLEX
    end if

    return
  end subroutine mpi_constants_precision

end module mpi
