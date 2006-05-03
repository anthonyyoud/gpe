program gpe
  ! Code to solve the Gross-Pitaevskii equation in 3D
  use parameters
  use ic
  use io
  use derivs
  use solve
  use variables
  implicit none

  integer :: p=0, p_start=0, j, k
  real :: norm=0.0, prev_norm=0.0
  type (var) :: psi
  logical :: run_exist, state_exist
  complex, dimension(0:nx1,0:ny1,0:nz1) :: tmp, tmp_var

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call system('hostname')

  call setup_itable()
  call para_range(0, nz1, nzprocs, myrankz, ksta, kend)
  call para_range(0, ny1, nyprocs, myranky, jsta, jend)
  call array_len(jlen, klen)
  call neighbours()
  call get_unit_no()
  call get_dirs()

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  allocate(psi%old(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
  allocate(psi%new(0:nx1,jsta:jend,ksta:kend))
  allocate(works1(0:nx1,2,kksta:kkend))
  allocate(works2(0:nx1,2,kksta:kkend))
  allocate(workr1(0:nx1,2,kksta:kkend))
  allocate(workr2(0:nx1,2,kksta:kkend))

  if (myrank == 0) then
    do k=-1,nzprocs
      print*, (itable(j,k), j=-1,nyprocs)
    end do
  end if
  print*, jsta, jend, ksta, kend, kksta, kkend, jlen, klen

  print*

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
  if (myrank == 0) then
    call open_files()
  end if

  ! Get the space mesh
  call get_grid()

  ! Get the initial conditions
  call ics(psi%new, p_start)
  !call get_bcs(psi%new)

  ! Get zeros of the real and imaginary part
  !call get_zeros(psi%new, p)
  !call get_re_im_zeros(psi%new, p)
  !call get_phase_zeros(psi%new, p)
  !call save_surface(p, psi%new)
  !call idl_surface(p, psi%new)
  !stop

  ! Get the uniform velocity U
  !call get_U(psi%new, U)

  ! If this is not a restart...
  if (.not. restart) then
    inquire(file=proc_dir//'end_state.dat', exist=state_exist)
    ! Exit if not doing restart but end_state.dat exists
    if (state_exist) stop 'ERROR: restart=.false. but end_state.dat exists.'
  end if
  
  p = p_start
  
  ! Calculate the norm of the initial condition
  call get_norm(psi%new, prev_norm)
  if (myrank == 0) then
    print*, "NORM:", prev_norm
  end if
  !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !stop
 
  ! Begin real time loop
  do while (t <= end_time)

    ! Check to see whether the RUNNING file exists
    if (myrank == 0) then
      inquire(file='RUNNING', exist=run_exist)
      ! If it doesn't then end the run
      if (.not. run_exist) then
        print*
        print*, 'Stop requested.'
        end_proc = 1
      end if
    end if

    call MPI_BCAST(end_proc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
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
    psi%old(:,jsta:jend,ksta:kend) = psi%new
    
    call send_recv_z(psi%old)
    call pack_y(psi%old)

    call send_recv_y()
    call unpack_y(psi%old)

    tmp_var = 0.0
    do k=ksta,kend
      do j=jsta,jend
        tmp_var(:,j,k) = psi%new(:,j,k)
      end do
    end do

    call MPI_REDUCE(tmp_var, tmp, nx*ny*nz, &
                    MPI_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                    
    ! Save time-series data
    if (modulo(t+im_t, abs(dt)*save_rate) < abs(dt)) then
    !if (mod(p, save_rate) == 0) then
      call save_energy(t, psi%old)
      call save_momentum(t, psi%old)
      if (myrank == 0) then
        !call save_time(t, psi%new)
        call save_time(t, tmp)
        !call save_energy(t, psi%old)
        !call save_momentum(t, psi%old)
        !call save_linelength(t, psi%new)
        !call save_linelength(t, tmp)
      end if
    end if
    
    !if (modulo(t+im_t, abs(dt)*save_rate2) < abs(dt)) then
    if (mod(p, save_rate2) == 0) then
      if (myrank == 0) then
        ! Find the zeros of the wave function
        !call get_zeros(psi%new, p)
        !call get_extra_zeros(psi%new, p)
        !call get_re_im_zeros(psi%new, p)
        call get_zeros(tmp, p)
        call get_extra_zeros(tmp, p)
        call get_re_im_zeros(tmp, p)
        !call get_phase_zeros(psi%new, p)
        if (save_contour) then
          ! Save 2D contour data
          !call save_surface(p, psi%new)
          call save_surface(p, tmp)
        end if
        if (idl_contour) then
          ! Save 3D isosurface data for use in IDL
          !call idl_surface(p, psi%new)
          call idl_surface(p, tmp)
        end if
      end if
    end if

    ! Periodically save the statea
    if (mod(p, save_rate2) == 0) then
      call end_state(psi%new, p, 0)
    end if

    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !call send_recv_z(psi%old)
    !call pack_y(psi%old)

    !call send_recv_y()
    !call unpack_y(psi%old)

    ! Call the solver subroutine to solve the equation
    call solver(psi%old, psi%new)
    !call get_bcs(psi%new)
    
    ! Calculate the norm
    call get_norm(psi%new, norm)
    !call get_norm(tmp, norm)

    call MPI_BCAST(norm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Calculate the relative difference of successive norms.  If the difference
    ! is small enough switch to real time
    if (.not. real_time) then
      if (abs(norm-prev_norm)/abs(prev_norm) < 1e-8) then
        real_time = .true.
        switched = .true.
        im_t = 0.0
        if (myrank == 0) then
          !call save_surface(p, psi%new)
          !call idl_surface(p, psi%new)
          call save_surface(p, tmp)
          call idl_surface(p, tmp)
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
    
    !print*, dt, t, im_t, norm, abs(norm-prev_norm)/abs(prev_norm)
    prev_norm = norm

    p = p+1

  end do
  
  ! Time loop finished so cleanly end the run
  call end_state(psi%new, p, 1)
  if (myrank == 0) then
    !call save_surface(p, psi%new)
    !call idl_surface(p, psi%new)
    !call end_state(tmp, p, 1)
    call save_surface(p, tmp)
    call idl_surface(p, tmp)
  end if

  ! Close runtime files
  if (myrank == 0) then
    call close_files()
  end if

  call MPI_FINALIZE(ierr)

end program gpe
