module solve
  ! Routines to actually solve the equation
  use parameters
  implicit none

  private
  public :: solver

  contains

  subroutine solver(var_in, var_out)
    ! General solver routine which calls a specific scheme
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: var_in
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: var_out

    select case (scheme)
      case ('euler')
        ! Explicit Euler time stepping
        call euler(var_in, var_out)
      case ('rk4')
        ! Explicit fourth order Runge-Kutta time stepping
        call rk4(var_in, var_out)
      case ('rk_adaptive')
        ! Explicit fifth order Runge-Kutta-Fehlberg adaptive time stepping
        call rkqs(var_in, var_out)
      case default
        ! Error if scheme not recognised
        STOP 'ERROR: Unrecognised time stepping scheme'
    end select

    return
  end subroutine solver

  subroutine euler(old, new)
    ! Explicit Euler time stepping
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: old
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: new
    complex, dimension(0:nx,0:ny,0:nz) :: rhs

    call get_rhs(old, rhs)

    new = old + dt*rhs

    t = t+real(dt)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine euler
  
  subroutine rk4(old, new)
    ! Explicit fourth order Runge-Kutta time stepping
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: old
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: new
    complex, dimension(4,0:nx,0:ny,0:nz) :: k

    call get_rhs(old,k(1,:,:,:))
    k(1,:,:,:) = dt*k(1,:,:,:)
    call get_rhs(old+0.5*k(1,:,:,:), k(2,:,:,:))
    k(2,:,:,:) = dt*k(2,:,:,:)
    call get_rhs(old+0.5*k(2,:,:,:), k(3,:,:,:))
    k(3,:,:,:) = dt*k(3,:,:,:)
    call get_rhs(old+k(3,:,:,:), k(4,:,:,:))
    k(4,:,:,:) = dt*k(4,:,:,:)

    new = old + k(1,:,:,:)/6.0 + k(2,:,:,:)/3.0 + &
                k(3,:,:,:)/3.0 + k(4,:,:,:)/6.0

    t = t+real(dt)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine rk4

  subroutine rkqs(in_var, out_var)
    ! Adaptive Runge-Kutta-Fehlberg time stepping
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: out_var
    complex, dimension(0:nx,0:ny,0:nz) :: tmp_var, err_var, psi_scal
    real :: errmax, tnew
    complex :: dt_temp, dt_next, dt_did

    ! General error condition
    errcon = (5.0/safety)**(1.0/dt_increase)

    ! Get the RHS for use as a scaling variable
    call get_rhs(in_var, psi_scal)
    ! Specific scaling condition
    psi_scal = abs(in_var) + dt*abs(psi_scal) + 1e-30

    do
      ! Do a Runge-Kutta step
      call rkck(in_var, tmp_var, err_var)
      ! Get the maximum error over the whole field
      errmax = maxval(abs(err_var/psi_scal))/eps
      ! Step succeeded so exit
      if (errmax <= 1.0) exit
      ! Step didn't succeed so decrease the time step
      dt_temp = safety*dt*(errmax**dt_decrease)
      ! Don't decrease the time step by more than a factor of ten
      if (.not. real_time) then
        dt = cmplx(0.0,-sign(max(abs(dt_temp), 0.1*abs(dt)), abs(dt)))
      else
        dt = cmplx(sign(max(abs(dt_temp), 0.1*abs(dt)), abs(dt)), 0.0)
      end if
      ! New time
      tnew = t+real(dt)
      if (real_time) then
        if (tnew == t) then
          ! Guard against infinitesimal time steps
          print*, 'WARNING: Timestep underflow in rkqs()'
        end if
      end if
    end do
    
    if (errmax > errcon) then
      ! Increase the time step
      dt_next = safety*dt*(errmax**dt_increase)
    else
      ! But not by more than a factor of 5
      dt_next = 5.0*dt
    end if

    ! Time step that was actually performed
    dt_did = dt
    ! Update the time
    t = t+real(dt)
    ! Time step to try next time
    dt = dt_next
    ! Update imaginary time
    im_t = im_t+abs(aimag(dt))
    
    ! Update the variable
    out_var = tmp_var

    !if (abs(dt) <= 1e-5) then
      write (13, '(3e17.9)') t, im_t, abs(dt)
    !end if

    return
  end subroutine rkqs

  subroutine rkck(old, new, err)
    ! Explicit fifth order Runge--Kutta--Fehlberg time stepping
    use parameters
    implicit none

    real, parameter :: b21      = 0.2
    real, parameter :: b31      = 0.075
    real, parameter :: b32      = 0.225
    real, parameter :: b41      = 0.3
    real, parameter :: b42      = -0.9
    real, parameter :: b43      = 1.2
    real, parameter :: b51      = -11.0 / 54.0
    real, parameter :: b52      = 2.5
    real, parameter :: b53      = -70.0 / 27.0
    real, parameter :: b54      = 35.0 / 27.0
    real, parameter :: b61      = 1631.0 / 55296.0
    real, parameter :: b62      = 175.0 / 512.0
    real, parameter :: b63      = 575.0 / 13824.0
    real, parameter :: b64      = 44275.0 / 110592.0
    real, parameter :: b65      = 253.0 / 4096.0
    real, parameter :: c1       = 37.0 / 378.0
    real, parameter :: c2       = 0.0
    real, parameter :: c3       = 250.0 / 621.0
    real, parameter :: c4       = 125.0 / 594.0
    real, parameter :: c5       = 0.0
    real, parameter :: c6       = 512.0 / 1771.0
    real, parameter :: dc1      = c1 - 2825.0 / 27648.0
    real, parameter :: dc2      = c2 - 0.0
    real, parameter :: dc3      = c3 - 18575.0 / 48384.0
    real, parameter :: dc4      = c4 - 13525.0 / 55296.0
    real, parameter :: dc5      = c5 - 277.0 / 14336.0
    real, parameter :: dc6      = c6 - 0.25
    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: old
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: new, err
    complex, dimension(6,0:nx,0:ny,0:nz) :: k

    call get_rhs(old,k(1,:,:,:))
    k(1,:,:,:) = dt*k(1,:,:,:)
    call get_rhs(old+b21*k(1,:,:,:), k(2,:,:,:))
    k(2,:,:,:) = dt*k(2,:,:,:)
    call get_rhs(old+b31*k(1,:,:,:)+&
                     b32*k(2,:,:,:), k(3,:,:,:))
    k(3,:,:,:) = dt*k(3,:,:,:)
    call get_rhs(old+b41*k(1,:,:,:)+&
                     b42*k(2,:,:,:)+&
                     b43*k(3,:,:,:), k(4,:,:,:))
    k(4,:,:,:) = dt*k(4,:,:,:)
    call get_rhs(old+b51*k(1,:,:,:)+&
                     b52*k(2,:,:,:)+&
                     b53*k(3,:,:,:)+&
                     b54*k(4,:,:,:), k(5,:,:,:))
    k(5,:,:,:) = dt*k(5,:,:,:)
    call get_rhs(old+b61*k(1,:,:,:)+&
                     b62*k(2,:,:,:)+&
                     b63*k(3,:,:,:)+&
                     b64*k(4,:,:,:)+&
                     b65*k(5,:,:,:), k(6,:,:,:))
    k(6,:,:,:) = dt*k(6,:,:,:)

    new = old + c1*k(1,:,:,:) + c2*k(2,:,:,:) + &
                c3*k(3,:,:,:) + c4*k(4,:,:,:) + &
                c5*k(5,:,:,:) + c6*k(6,:,:,:)

    err = dc1*k(1,:,:,:) + dc2*k(2,:,:,:) + &
          dc3*k(3,:,:,:) + dc4*k(4,:,:,:) + &
          dc5*k(5,:,:,:) + dc6*k(6,:,:,:)

    return
  end subroutine rkck
  
  subroutine get_rhs(in_var, rhs)
    ! Get the right-hand-side of the equation
    use parameters
    use variables
    use derivs
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: in_var
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: rhs
    complex, dimension(0:nx,0:ny,0:nz) :: dpsidx, dpsidz
    real :: U=0.0 !0.18
    
    call deriv_x(in_var, dpsidx)
    !call deriv_z(in_var, dz)
    
    rhs = 0.5*eye * ( laplacian(in_var) + (1.0-abs(in_var)**2)*in_var ) + &
          U*dpsidx
    !rhs = 0.5 * ( laplacian(in_var) + (1.0-abs(in_var)**2)*in_var ) - &
    !      eye*U*dpsidx

    !write (90, '(2e17.9)'), rhs(nx/2,ny/2,nz/2)

    return
  end subroutine get_rhs
    
end module solve
