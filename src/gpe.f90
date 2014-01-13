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

program gpe
  ! Code to solve the Gross-Pitaevskii equation in 3D.  Parallelised using MPI.
  ! See the README file for a description and usage instructions
  use decomp_2d
  use derivs
  use error
  use ic
  use io
  use mpi
  use parameters
  use solve
  use variables
  implicit none

  integer :: n=0, m=0, ps=0
  real (pr) :: norm=0.0_pr, prev_norm=1.0_pr, relnorm=0.0_pr
  type (var) :: psi, test
  logical :: run_exist, state_exist

  ! Initialise the MPI process grid
  call MPI_INIT(ierr)
  call decomp_2d_init(nx, ny, nz, nyprocs, nzprocs)
  js = xstart(2)-1; je = xend(2)-1
  ks = xstart(3)-1; ke = xend(3)-1
  xs1 = zstart(1) - 1; xe1 = zend(1) - 1
  ys1 = zstart(2) - 1; ye1 = zend(2) - 1
  zs1 = zstart(3) - 1; ze1 = zend(3) - 1

  myrank = nrank
  !call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  !call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  ! Get unit numbers so that files can be opened on each process
  call get_unit_no()
  ! Read the run-time parameters from run.in.
  call read_run_params()
  ! Assign the replacement MPI constants which allow for real/double precision.
  call mpi_constants_precision()
  ! Setup the lookup table for neighbouring processes
  call setup_itable()
  ! Calculate the start and end array indices on each process
  !call para_range(0, nz1, nzprocs, myrankz, ks, ke)
  !call para_range(0, ny1, nyprocs, myranky, js, je)
  ! Calculate the array dimensions on each process
  call array_len()
  ! Get the neighbouring process rank
  call neighbours()
  ! Get directory names for each process
  call get_dirs()

  ! Allocate array dimensions on each process
  allocate(psi%old(0:nx1,js-2:je+2,ks-2:ke+2))
  if (diagnostic) then
    allocate(psi%old2(0:nx1,js-2:je+2,ks-2:ke+2))
  end if
  allocate(psi%new(0:nx1,js:je,ks:ke))
  allocate(test%new(0:nx1,js:je,ks:ke))
  ! (work send and receive arrays)
  allocate(works1(0:nx1,2,kks:kke))
  allocate(works2(0:nx1,2,kks:kke))
  allocate(workr1(0:nx1,2,kks:kke))
  allocate(workr2(0:nx1,2,kks:kke))
  allocate(ave(0:nx1,js:je,ks:ke))

  ! Check which time stepping scheme we're using
  if (.not. pp_filtered_surface) then
    call print_runtime_info()
  end if
  
  ! Set the time step
  if (real_time) then
    dt = cmplx(tau, 0.0_pr, pr)
  else
    dt = cmplx(0.0_pr, -tau, pr)
  end if
  
  ! Open runtime files
  call open_files()

  ! Get the space mesh
  call get_grid()

  ! Get the cutoff wavenumber and amplitude
  call get_kc_amp()

  ! Post-process filtered isosurfaces
  if (pp_filtered_surface) then
    if (myrank == 0) then
      print*, 'Saving filtered isosurfaces'
    end if
    call pp_save_filter()
    if (myrank == 0) then
      open (99, file = 'RUNNING')  !delete 'RUNNING' to finish run
      close (99, status = 'delete')
    end if
  end if

  if (.not. pp_filtered_surface) then
  ! Get the initial conditions
  call ics(psi%new)
  if (save_3d) then
    call idl_surface(psi%new)
  end if
  if (eqn_to_solve == 4 .and. .not. real_time) then
    call get_norm(psi%new, prev_norm)
    if (renorm) call renormalise(psi%new, prev_norm)
    call save_norm(t, prev_norm, relnorm)
  end if
  call condensed_particles(t, psi%new)
  
  ! If this is not a restart...
  if (.not. restart) then
    inquire(file=end_state_file, exist=state_exist)
    ! Exit if not doing restart but end_state.dat exists
    if (state_exist) then
      call emergency_stop('ERROR: restart=.false. but ' &
        //end_state_file//' exists.')
    end if
  end if
  
  ps = int(t/p_save)
  n = int(t/save_rate2)
  m = int(t/save_rate3)
  ave = 0.0_pr
  
  psi%old = 0.0_pr
 
  ! Begin real time loop
  do while (t+im_t <= end_time)
    ! Check to see whether the RUNNING file exists
    if (myrank == 0) then
      inquire(file='RUNNING', exist=run_exist)
      ! If it doesn't then end the run
      if (.not. run_exist) then
        print*
        print*, 'Stop requested.'
        ! Flag to send to all processes
        end_proc = 1
      end if
    end if

    ! Send the end_proc flag to all processes
    call MPI_BCAST(end_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    ! If the end_proc flag has been set, then end the run
    if (end_proc == 1) then
      call end_state(psi%new, 1)
      if (myrank == 0) then
        open (99, file = 'RUNNING')  !delete 'RUNNING' to finish run
        close (99, status = 'delete')
        exit
      end if
      exit
    end if

    ! Update the variable
    if (diagnostic) then
      psi%old2 = psi%old
    end if
    psi%old(:,js:je,ks:ke) = psi%new

    ! Send and receive the variable so that all processes have the boundary
    ! data from neighbouring processes.  Since the data needing to be sent in
    ! the y-direction is not contiguous in memory, it must first be packed into
    ! a contiguous array, then sent, and finally unpacked.
    call send_recv_z(psi%old)
    call pack_y(psi%old)
    call send_recv_y()
    call unpack_y(psi%old)

    ! Save time-series data
    if (modulo(t+im_t, abs(dt)*save_rate) < abs(dt)) then
      call save_energy(t, psi%old)
      !call save_mass(t, psi%new)
      !call save_momentum(t, psi%old)
      call save_time(t, psi%new)
      if (save_ll) then
        call save_linelength(psi%old, 0)
      end if
    end if
    
    if (t+im_t >= save_rate3*m) then
      call condensed_particles(t, psi%new)
      if (save_pdf) then
        call save_velocity_pdf(psi%old)
      end if
      if (save_vcf) then
        call save_vel_corr(psi%old)
      end if
      m = m+1
    end if
    
    if (t+im_t >= save_rate2*n) then
      if (save_zeros) then
        ! Find the zeros of the wave function
        call get_zeros(psi%old)
        call get_extra_zeros(psi%old)
        call get_re_im_zeros(psi%old)
      end if
      if (save_contour) then
        ! Save 2D contour data
        call save_surface(psi%new)
      end if
      if (save_3d) then
        ! Save 3D isosurface data for use in IDL
        call idl_surface(psi%new)
      end if
      if (save_average) then
        ! Save time-averaged data
        call average(psi%new)
      end if
      n = n+1
    end if

    ! Periodically save the state
    if (t+im_t >= p_save*ps) then
      ! 0 is a flag which means the run has not yet ended
      call end_state(psi%new, 0)
      ps = ps+1
    end if

    ! Call the solver subroutine to solve the equation
    call solver(psi%old, psi%new)

    ! Check to see whether the file SAVE exists in the run directory, and save
    ! 3D isosurfaces if it does.
    call save_run(psi%new)

    if (t+im_t >= save_rate2*n) then
      if (diagnostic) then
        call diag(psi%old2, psi%old, psi%new)
      end if
    end if
    
    ! Calculate the norm and renormalise if necessary.
    if (eqn_to_solve == 4) then
      call get_norm(psi%new, norm)
      relnorm = abs((norm-prev_norm)/prev_norm)
      if (modulo(t+im_t, abs(dt)*save_rate) < abs(dt)) then
        call save_norm(t, norm, relnorm)
      end if
    end if

    ! Calculate the relative difference of successive norms.  If the difference
    ! is small enough switch to real time
    if (eqn_to_solve == 4 .and. .not. real_time) then
      if (relnorm < 1e-12_pr) then
        if (stop_imag) exit  ! End run if stop after imaginary time requested.
        real_time = .true.
        im_t = 0.0_pr
        ps = int(t/p_save)
        n = int(t/save_rate2)
        m = int(t/save_rate3)
        if (save_3d) then
          call idl_surface(psi%new)
        end if
        dt = cmplx(abs(aimag(dt)), abs(real(dt, pr)), pr)
        if (myrank == 0) then
          print*, 'Switching to real time'
          print*, 't =', im_t
          print*, 'p =', p
        end if
      end if
      if (renorm) call renormalise(psi%new, norm)
      if (imprint_vl) psi%new = imprint_vortex_line(psi%new)
    end if
    
    ! Update the norm
    prev_norm = norm

    ! Update the time index
    p = p+1
  end do
  
  ! Time loop finished so cleanly end the run
  call end_state(psi%new, 1)
  if (save_3d) then
    call idl_surface(psi%new)
  end if
  if (save_contour) then
    call save_surface(psi%new)
  end if
  end if

  ! Close runtime files
  call close_files()
  
  ! Deallocate arrays
  deallocate(psi%old)
  if (diagnostic) then
    deallocate(psi%old2)
  end if
  deallocate(psi%new)
  deallocate(works1)
  deallocate(works2)
  deallocate(workr1)
  deallocate(workr2)

  if (myrank == 0) then
    print*, 'DONE!'
  end if

  ! Stop the MPI process grid
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)

end program gpe
