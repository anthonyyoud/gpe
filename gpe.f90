program gpe
  ! Code to solve the Gross-Pitaevskii equation in 3D.  Parallelised using MPI.
  ! See the README file for a description and usage instructions
  use parameters
  use ic
  use io
  use derivs
  use solve
  use variables
  implicit none

  integer    :: p_start=0, n=0
  real       :: norm=0.0, prev_norm=0.0
  type (var) :: psi
  logical    :: run_exist, state_exist

  ! Initialise the MPI process grid
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  ! Find out on which host each process is running
  call system('hostname')

  ! Setup the lookup table for neighbouring processes
  call setup_itable()
  ! Calculate the start and end array indices on each process
  call para_range(0, nz1, nzprocs, myrankz, ksta, kend)
  call para_range(0, ny1, nyprocs, myranky, jsta, jend)
  ! Calculate the array dimensions on each process
  call array_len(jlen, klen)
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
  ! (work send and receive arrays)
  allocate(works1(0:nx1,2,kksta:kkend))
  allocate(works2(0:nx1,2,kksta:kkend))
  allocate(workr1(0:nx1,2,kksta:kkend))
  allocate(workr2(0:nx1,2,kksta:kkend))

  ! Check which time stepping scheme we're using
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
        STOP 'ERROR: Unrecognised time stepping scheme'
    end select
  end if
  
  ! Set the time step
  dt = time_step
  if (real_time) then
    dt = cmplx(abs(aimag(time_step)), 0.0)
  end if
  
  ! Open runtime files
  call open_files()

  ! Get the space mesh
  call get_grid()

  ! Get the initial conditions
  call ics(psi%new, p_start)
  
  ! Call the FFT routine to transform the initial condition
  !call fft(psi%new)

  ! If this is not a restart...
  if (.not. restart) then
    inquire(file=proc_dir//'end_state.dat', exist=state_exist)
    ! Exit if not doing restart but end_state.dat exists
    if (state_exist) stop 'ERROR: restart=.false. but end_state.dat exists.'
  end if
  
  p = p_start
  
  ! Calculate the norm of the initial condition
  call get_norm(psi%new, prev_norm)

  n = int(t/save_rate2)
 
  ! Begin real time loop
  do while (t <= end_time)

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
      call end_state(psi%new, p, 1)
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
    !if (mod(p, save_rate) == 0) then
      call save_energy(t, psi%old)
      !call save_mass(t, psi%new)
      call save_momentum(t, psi%old)
      call save_time(t, psi%new)
      call save_linelength(t, psi%old)
      !call condensed_particles(t, psi%new)
    end if
    
    !if (modulo(t+im_t, abs(dt)*save_rate2) < abs(dt)) then
    !if (mod(p, save_rate2) == 0) then
    if (t+im_t >= save_rate2*n) then
      if (save_zeros) then
        ! Find the zeros of the wave function
        call get_zeros(psi%old, p)
        call get_extra_zeros(psi%old, p)
        call get_re_im_zeros(psi%old, p)
      end if
      if (save_contour) then
        ! Save 2D contour data
        call save_surface(p, psi%new)
      end if
      if (save_3d) then
        ! Save 3D isosurface data for use in IDL
        call idl_surface(p, psi%new)
      end if
      n = n+1
    end if

    ! Periodically save the state
    !if (mod(p, save_rate2) == 0) then
    if (t+im_t >= 10.0*n) then
      ! 0 is a flag which means the run has not yet ended
      call end_state(psi%new, p, 0)
    end if

    ! Call the solver subroutine to solve the equation
    call solver(psi%old, psi%new)
    
    if (t+im_t >= save_rate2*n) then
      if (diagnostic) then
        call diag(psi%old2, psi%old, psi%new, p)
      end if
    end if
    
    ! Calculate the norm
    call get_norm(psi%new, norm)

    ! Make sure all process know what the norm is
    call MPI_BCAST(norm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Calculate the relative difference of successive norms.  If the difference
    ! is small enough switch to real time
    if (.not. real_time) then
      if (abs(norm-prev_norm)/abs(prev_norm) < 1e-8) then
        real_time = .true.
        switched = .true.
        im_t = 0.0
        call idl_surface(p, psi%new)
        call save_surface(p, psi%new)
        if (myrank == 0) then
          print*, 'Switching to real time'
        end if
      end if
    end if
    
    ! Flag to denote whether this is the first time step after switching to
    ! real time
    if (switched) then
      dt = cmplx(abs(aimag(dt)), abs(real(dt)))
      switched = .false.
    end if
    
    ! Update the norm
    prev_norm = norm

    ! Update the time index
    p = p+1

  end do
  
  ! Time loop finished so cleanly end the run
  call end_state(psi%new, p, 1)
  call idl_surface(p, psi%new)
  call save_surface(p, psi%new)

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
