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
