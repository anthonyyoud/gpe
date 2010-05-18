! $Id$
!----------------------------------------------------------------------------

program gpe
  ! Code to solve the Gross-Pitaevskii equation in 3D.  Parallelised using MPI.
  ! See the README file for a description and usage instructions
  use derivs
  use error
  use ic
  use io
  use parameters
  use solve
  use variables
  implicit none

  integer    :: n=0, m=0, ps=0
  real       :: norm=0.0, prev_norm=0.0
  type (var) :: psi, test
  logical    :: run_exist, state_exist

  ! Initialise the MPI process grid
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  ! Setup the lookup table for neighbouring processes
  call setup_itable()
  ! Calculate the start and end array indices on each process
  call para_range(0, nz1, nzprocs, myrankz, ksta, kend)
  call para_range(0, ny1, nyprocs, myranky, jsta, jend)
  ! Calculate the array dimensions on each process
  call array_len()
  ! Get the neighbouring process rank
  call neighbours()
  ! Get unit numbers so that files can be opened on each process
  call get_unit_no()
  ! Get directory names for each process
  call get_dirs()

  ! Allocate array dimensions on each process
  allocate(psi%old(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
  if (diagnostic) then
    allocate(psi%old2(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
  end if
  allocate(psi%new(0:nx1,jsta:jend,ksta:kend))
  allocate(test%new(0:nx1,jsta:jend,ksta:kend))
  ! (work send and receive arrays)
  allocate(works1(0:nx1,2,kksta:kkend))
  allocate(works2(0:nx1,2,kksta:kkend))
  allocate(workr1(0:nx1,2,kksta:kkend))
  allocate(workr2(0:nx1,2,kksta:kkend))
  allocate(ave(0:nx1,jsta:jend,ksta:kend))

  ! Check which time stepping scheme we're using
  if (.not. pp_filtered_surface) then
    if (myrank == 0) then
      select case (scheme)
        case ('euler')
          print*, 'Explicit Euler time stepping'
        case ('rk4')
          print*, 'Explicit fourth order Runge-Kutta time stepping'
        case ('rk_adaptive')
          print*, 'Explicit fifth order &
                  &Runge-Kutta-Fehlberg adaptive time stepping'
        case default
          call emergency_stop('ERROR: Unrecognised time stepping scheme.')
      end select
    end if

    if (myrank == 0) then
      select case (eqn_to_solve)
        case (1)
          print*, 'Solving CASE 1'
        case (2)
          print*, 'Solving CASE 2'
        case (3)
          print*, 'Solving CASE 3'
      end select
    end if
  end if
  
  ! Set the time step
  if (real_time) then
    dt = cmplx(tau, 0.0)
  else
    dt = cmplx(0.0, -tau)
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
  if (eqn_to_solve == 4 .and. .not. real_time) then
    call get_norm(psi%new, prev_norm)
    call renormalise(psi%new, prev_norm)
    call save_norm(t, prev_norm)
  end if
  if (save_3d) then
    call idl_surface(psi%new)
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
  ave = 0.0
  
  psi%old = 0.0
 
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
    psi%old(:,jsta:jend,ksta:kend) = psi%new

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

    if (t+im_t >= save_rate2*n) then
      if (diagnostic) then
        call diag(psi%old2, psi%old, psi%new)
      end if
    end if
    
    ! Calculate the norm
    if (eqn_to_solve == 4) then
      call get_norm(psi%new, norm)
      if (modulo(t+im_t, abs(dt)*save_rate) < abs(dt)) then
        call save_norm(t, norm)
      end if
      if (.not. real_time) then
        call renormalise(psi%new, norm)
      end if
    end if

    ! Make sure all processes know what the norm is
    !call MPI_BCAST(norm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Calculate the relative difference of successive norms.  If the difference
    ! is small enough switch to real time
    if (eqn_to_solve == 4 .and. .not. real_time) then
      if (abs(norm-prev_norm)/abs(prev_norm) < 1e-8) then
        real_time = .true.
        im_t = 0.0
        ps = int(t/p_save)
        n = int(t/save_rate2)
        m = int(t/save_rate3)
        if (save_3d) then
          call idl_surface(psi%new)
        end if
        dt = cmplx(abs(aimag(dt)), abs(real(dt)))
        if (myrank == 0) then
          print*, 'Switching to real time'
          print*, 't =', im_t
          print*, 'p =', p
        end if
      end if
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
