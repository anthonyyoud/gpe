program gpe
  ! Code to solve the Gross-Pitaevskii equation in 3D
  use parameters
  use ic
  use io
  use derivs
  use solve
  use variables
  implicit none

  integer :: p=0, p_start=0
  real :: norm=0.0, prev_norm=0.0
  type (var) :: psi
  logical :: run_exist, state_exist

  print*

  ! Check which time stepping scheme we're using
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
    inquire(file='end_state.dat', exist=state_exist)
    ! Exit if not doing restart but end_state.dat exists
    if (state_exist) stop 'ERROR: restart=.false. but end_state.dat exists.'
  end if
  
  p = p_start
  
  ! Calculate the norm of the initial condition
  call get_norm(psi%new, prev_norm)
 
  ! Begin real time loop
  do while (t <= end_time)

    ! Check to see whether the RUNNING file exists
    inquire(file='RUNNING', exist=run_exist)
    ! If it doesn't then end the run
    if (.not. run_exist) then
      print*
      print*, 'Stop requested.'
      call end_state(psi%new, p, 1)
      exit
    end if

    ! Update the variable
    psi%old = psi%new

    ! Save time-series data
    if (modulo(t+im_t, abs(dt)*save_rate) < abs(dt)) then
    !if (mod(p, save_rate) == 0) then
      call save_time(t, psi%new)
      call save_energy(t, psi%new)
      call save_momentum(t, psi%new)
      call save_linelength(t, psi%new)
    end if
    
    !if (modulo(t+im_t, abs(dt)*save_rate2) < abs(dt)) then
    if (mod(p, save_rate2) == 0) then
      ! Find the zeros of the wave function
      call get_zeros(psi%new, p)
      call get_extra_zeros(psi%new, p)
      call get_re_im_zeros(psi%new, p)
      !call get_phase_zeros(psi%new, p)
      if (save_contour) then
        ! Save 2D contour data
        call save_surface(p, psi%new)
      end if
      if (idl_contour) then
        ! Save 3D isosurface data for use in IDL
        call idl_surface(p, psi%new)
      end if
    end if
    
    ! Periodically save the state
    if (mod(p, save_rate2) == 0) then
      call end_state(psi%new, p, 0)
    end if

    ! Call the solver subroutine to solve the equation
    call solver(psi%old, psi%new)
    !call get_bcs(psi%new)
    
    ! Calculate the norm
    call get_norm(psi%new, norm)

    ! Calculate the relative difference of successive norms.  If the difference
    ! is small enough switch to real time
    if (.not. real_time) then
      if (abs(norm-prev_norm)/abs(prev_norm) < 1e-8) then
        real_time = .true.
        switched = .true.
        im_t = 0.0
        call save_surface(p, psi%new)
        call idl_surface(p, psi%new)
        print*, 'Switching to real time'
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
  call save_surface(p, psi%new)
  call idl_surface(p, psi%new)

  ! Close runtime files
  call close_files()

end program gpe
