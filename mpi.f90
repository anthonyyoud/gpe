! $Id: mpi.f90 590 2010-05-18 15:04:46Z najy2 $
!----------------------------------------------------------------------------

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
