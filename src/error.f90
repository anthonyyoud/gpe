! $Id$
!----------------------------------------------------------------------------

module error
  ! Routines for error handling
  implicit none

  private
  public :: emergency_stop
  
  contains

  subroutine emergency_stop(err_string)
    use parameters
    implicit none

    character(*), intent(in) :: err_string

    if (myrank == 0) then
      ! Print error to screen, remove RUNNING and write error to ERROR.
      print*, err_string
      open (99, file = 'RUNNING')
      close (99, status = 'delete')
      open (96, file = 'ERROR')
      write (96, *) err_string
      close (96)
    end if

    ! Stop the MPI process grid.
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_FINALIZE(ierr)

    ! Emergency stop.
    stop
  end subroutine emergency_stop

end module error
