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

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: var_in
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: var_out

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

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: old
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: new
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: rhs

    call get_rhs(old, rhs)

    new = old(:,jsta:jend,ksta:kend) + dt*rhs

    t = t+real(dt)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine euler
  
  subroutine rk4(old, new)
    ! Explicit fourth order Runge-Kutta time stepping
    use parameters
    use variables
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: old
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: new
    complex, dimension(4,0:nx1,jsta-2:jend+2,ksta-2:kend+2) :: kk
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2) :: tmp
    integer :: j, k

    call get_rhs(old,kk(1,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(1,:,j,k) = dt*kk(1,:,j,k)
        tmp(:,j,k) = old(:,j,k)+0.5*kk(1,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(2,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(2,:,j,k) = dt*kk(2,:,j,k)
        tmp(:,j,k) = old(:,j,k)+0.5*kk(2,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(3,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(3,:,j,k) = dt*kk(3,:,j,k)
        tmp(:,j,k) = old(:,j,k)+kk(3,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(4,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(4,:,j,k) = dt*kk(4,:,j,k)
      end do
    end do

    do k=ksta,kend
      do j=jsta,jend
        new(:,j,k) = old(:,j,k) + kk(1,:,j,k)/6.0 + kk(2,:,j,k)/3.0 + &
                                  kk(3,:,j,k)/3.0 + kk(4,:,j,k)/6.0
      end do
    end do

    t = t+real(dt)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine rk4

  subroutine rkqs(in_var, out_var)
    ! Adaptive Runge-Kutta-Fehlberg time stepping. From Numerical Recipes
    ! section 16.2, page 708
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: tmp_var, err_var, psi_scal
    real :: tnew
    real, dimension(2) :: errmaxs, errmaxr
    complex :: dt_temp, dt_next, dt_did
    integer :: j, k

    ! General error condition
    errcon = (5.0/safety)**(1.0/dt_increase)

    ! Get the RHS for use as a scaling variable
    call get_rhs(in_var, psi_scal)
    ! Specific scaling condition
    psi_scal = abs(in_var(:,jsta:jend,ksta:kend)) + dt*abs(psi_scal) + 1e-30

    do
      ! Do a Runge-Kutta step
      call rkck(in_var, tmp_var, err_var)
      ! Get the maximum error over the whole field.
      errmaxs(1) = maxval(abs(err_var/psi_scal))/eps
      errmaxs(2) = 0.0

      ! errmaxs/r must be 2-element arrays since the MPI function MPI_MAXLOC
      ! calculates the array maximum as well as its location
      call MPI_ALLREDUCE(errmaxs, errmaxr, 1, MPI_2REAL, MPI_MAXLOC, &
                         MPI_COMM_WORLD, ierr)
                      
      ! Step succeeded so exit
      if (errmaxr(1) <= 1.0) exit
      ! Step didn't succeed so decrease the time step
      dt_temp = safety*dt*(errmaxr(1)**dt_decrease)
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
          if (myrank == 0) then
            ! Guard against infinitesimal time steps
            print*, 'WARNING: Timestep underflow in rkqs()'
          end if
        end if
      end if
    end do
    
    if (errmaxr(1) > errcon) then
      ! Increase the time step
      dt_next = safety*dt*(errmaxr(1)**dt_increase)
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

    if (myrank == 0) then
      !if (abs(dt) <= 1e-5) then
        write (11, '(3e17.9)') t, im_t, abs(dt)
      !end if
    end if

    return
  end subroutine rkqs

  subroutine rkck(old, new, err)
    ! Explicit fifth order Runge-Kutta-Fehlberg time stepping.  From Numerical
    ! Recipes, section 16.2, page 708
    use parameters
    use variables
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
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: old
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: new, err
    complex, dimension(6,0:nx1,jsta-2:jend+2,ksta-2:kend+2) :: kk
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2) :: tmp
    integer :: j, k

    call get_rhs(old,kk(1,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(1,:,j,k) = dt*kk(1,:,j,k)
        tmp(:,j,k) = old(:,j,k)+b21*kk(1,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(2,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(2,:,j,k) = dt*kk(2,:,j,k)
        tmp(:,j,k) = old(:,j,k)+b31*kk(1,:,j,k)+&
                                b32*kk(2,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(3,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(3,:,j,k) = dt*kk(3,:,j,k)
        tmp(:,j,k) = old(:,j,k)+b41*kk(1,:,j,k)+&
                                b42*kk(2,:,j,k)+&
                                b43*kk(3,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(4,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(4,:,j,k) = dt*kk(4,:,j,k)
        tmp(:,j,k) = old(:,j,k)+b51*kk(1,:,j,k)+&
                                b52*kk(2,:,j,k)+&
                                b53*kk(3,:,j,k)+&
                                b54*kk(4,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(5,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(5,:,j,k) = dt*kk(5,:,j,k)
        tmp(:,j,k) = old(:,j,k)+b61*kk(1,:,j,k)+&
                                b62*kk(2,:,j,k)+&
                                b63*kk(3,:,j,k)+&
                                b64*kk(4,:,j,k)+&
                                b65*kk(5,:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk(6,:,jsta:jend,ksta:kend))
    do k=ksta,kend
      do j=jsta,jend
        kk(6,:,j,k) = dt*kk(6,:,j,k)
      end do
    end do

    do k=ksta,kend
      do j=jsta,jend
        new(:,j,k) = old(:,j,k) + c1*kk(1,:,j,k) + c2*kk(2,:,j,k) + &
                                  c3*kk(3,:,j,k) + c4*kk(4,:,j,k) + &
                                  c5*kk(5,:,j,k) + c6*kk(6,:,j,k)

        err(:,j,k) = dc1*kk(1,:,j,k) + dc2*kk(2,:,j,k) + &
                     dc3*kk(3,:,j,k) + dc4*kk(4,:,j,k) + &
                     dc5*kk(5,:,j,k) + dc6*kk(6,:,j,k)
      end do
    end do

    return
  end subroutine rkck
  
  subroutine get_rhs(in_var, rhs)
    ! Get the right-hand-side of the equation
    use parameters
    use variables
    use derivs
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: rhs
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: dpsidx, dpsidz
    integer :: j, k
    real :: U=0.18
    
    call deriv_x(in_var, dpsidx)
    !call deriv_z(in_var, dz)
    
    rhs = 0.5*eye * ( laplacian(in_var) + &
                    (1.0-abs(in_var(:,jsta:jend,ksta:kend))**2)*&
                             in_var(:,jsta:jend,ksta:kend) ) + &
                     U*dpsidx
    !rhs = 0.5 * ( laplacian(in_var) + (1.0-abs(in_var)**2)*in_var ) - &
    !      eye*U*dpsidx

    return
  end subroutine get_rhs
    
end module solve
