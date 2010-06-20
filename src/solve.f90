! $Id$
!----------------------------------------------------------------------------

module solve
  ! Routines to actually solve the equation
  use parameters
  implicit none

  private
  public :: solver

  contains

  subroutine solver(var_in, var_out)
    ! General solver routine which calls a specific scheme
    use error
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: var_in
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: var_out

    select case (scheme)
      case ('euler')
        ! Explicit Euler time stepping
        call euler(var_in, var_out)
      case ('rk2')
        ! Explicit second order Runge-Kutta time stepping
        call rk2(var_in, var_out)
      case ('rk4')
        ! Explicit fourth order Runge-Kutta time stepping
        call rk4(var_in, var_out)
      case ('rk45')
        ! Explicit fifth order Runge-Kutta-Fehlberg adaptive time stepping
        call rkqs(var_in, var_out)
      case default
        ! Error if scheme not recognised
        call emergency_stop('ERROR: Unrecognised time stepping scheme.')
    end select

    return
  end subroutine solver

! ***************************************************************************  

  subroutine euler(old, new)
    ! Explicit Euler time stepping
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: old
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: new
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: rhs

    call get_rhs(old, rhs)

    new = old(:,js:je,ks:ke) + dt*rhs

    t = t+real(dt, pr)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine euler

! ***************************************************************************  

  subroutine rk2(old, new)
    ! Explicit second order Runge-Kutta time stepping
    use parameters
    use variables
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: old
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: new
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk1
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk2
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: tmp
    integer :: j, k

    call get_rhs(old,kk1(:,js:je,ks:ke))
    tmp(:,js:je,ks:ke) = old(:,js:je,ks:ke) + 0.5_pr*dt*kk1(:,js:je,ks:ke)
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk2(:,js:je,ks:ke))

    new(:,js:je,ks:ke) = old(:,js:je,ks:ke) + dt*kk2(:,js:je,ks:ke)

    t = t+real(dt, pr)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine rk2
  
! ***************************************************************************  

  subroutine rk4(old, new)
    ! Explicit fourth order Runge-Kutta time stepping
    use parameters
    use variables
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: old
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: new
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk1
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk2
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk3
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk4
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: tmp
    integer :: j, k

    call get_rhs(old,kk1(:,js:je,ks:ke))
    tmp(:,js:je,ks:ke) = old(:,js:je,ks:ke) + 0.5_pr*dt*kk1(:,js:je,ks:ke)
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk2(:,js:je,ks:ke))
    tmp(:,js:je,ks:ke) = old(:,js:je,ks:ke) + 0.5_pr*dt*kk2(:,js:je,ks:ke)
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk3(:,js:je,ks:ke))
    tmp(:,js:je,ks:ke) = old(:,js:je,ks:ke) + dt*kk3(:,js:je,ks:ke)
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk4(:,js:je,ks:ke))

    new(:,js:je,ks:ke) = old(:,js:je,ks:ke) + &
      dt*(kk1(:,js:je,ks:ke)/6.0_pr + &
          kk2(:,js:je,ks:ke)/3.0_pr + &
          kk3(:,js:je,ks:ke)/3.0_pr + &
          kk4(:,js:je,ks:ke)/6.0_pr)

    t = t+real(dt, pr)
    im_t = im_t+abs(aimag(dt))

    return
  end subroutine rk4

! ***************************************************************************  

  subroutine rkqs(in_var, out_var)
    ! Adaptive Runge-Kutta-Fehlberg time stepping. From Numerical Recipes
    ! section 16.2, page 708
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: out_var
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: tmp_var, err_var, psi_scal
    real (pr) :: tnew
    real (pr), dimension(2) :: errmaxs, errmaxr
    complex (pr) :: dt_temp, dt_next, dt_did
    integer :: j, k

    ! General error condition
    errcon = (5.0_pr/safety)**(1.0_pr/dt_increase)

    ! Get the RHS for use as a scaling variable
    call get_rhs(in_var, psi_scal)
    ! Specific scaling condition
    psi_scal = abs(in_var(:,js:je,ks:ke)) + dt*abs(psi_scal) + 1e-30_pr

    do
      ! Do a Runge-Kutta step
      call rkck(in_var, tmp_var, err_var)
      ! Get the maximum error over the whole field.
      errmaxs(1) = maxval(abs(err_var/psi_scal))/eps
      errmaxs(2) = 0.0_pr

      ! errmaxs/r must be 2-element arrays since the MPI function MPI_MAXLOC
      ! calculates the array maximum as well as its location
      call MPI_ALLREDUCE(errmaxs, errmaxr, 1, gpe_mpi_2real, MPI_MAXLOC, &
        MPI_COMM_WORLD, ierr)
                      
      ! Step succeeded so exit
      if (errmaxr(1) <= 1.0_pr) exit
      ! Step didn't succeed so decrease the time step
      dt_temp = safety*dt*(errmaxr(1)**dt_decrease)
      ! Don't decrease the time step by more than a factor of ten
      if (.not. real_time) then
        dt = cmplx(0.0_pr,-sign(max(abs(dt_temp), 0.1_pr*abs(dt)), abs(dt)), pr)
      else
        dt = cmplx(sign(max(abs(dt_temp), 0.1_pr*abs(dt)), abs(dt)), 0.0_pr, pr)
      end if
      ! New time
      tnew = t+real(dt, pr)
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
      dt_next = 5.0_pr*dt
    end if

    ! Time step that was actually performed
    dt_did = dt
    ! Update the time
    t = t+real(dt, pr)
    ! Time step to try next time
    dt = dt_next
    ! Update imaginary time
    im_t = im_t+abs(aimag(dt))
    
    ! Update the variable
    out_var = tmp_var

    if (myrank == 0) then
      !if (abs(dt) <= 1e-2_pr) then
      if (mod(p, save_rate) == 0) then
        open (11, status='unknown', position='append', file='timestep.dat')
        write (11, '(3e17.9,i10)') t, im_t, abs(dt), p
        close (11)
      end if
    end if

    return
  end subroutine rkqs

! ***************************************************************************  

  subroutine rkck(old, new, err)
    ! Explicit fifth order Runge-Kutta-Fehlberg time stepping.  From Numerical
    ! Recipes, section 16.2, page 708
    use parameters
    use variables
    implicit none

    real (pr), parameter :: b21 = 0.2_pr
    real (pr), parameter :: b31 = 0.075_pr
    real (pr), parameter :: b32 = 0.225_pr
    real (pr), parameter :: b41 = 0.3_pr
    real (pr), parameter :: b42 = -0.9_pr
    real (pr), parameter :: b43 = 1.2_pr
    real (pr), parameter :: b51 = -11.0_pr / 54.0_pr
    real (pr), parameter :: b52 = 2.5_pr
    real (pr), parameter :: b53 = -70.0_pr / 27.0_pr
    real (pr), parameter :: b54 = 35.0_pr / 27.0_pr
    real (pr), parameter :: b61 = 1631.0_pr / 55296.0_pr
    real (pr), parameter :: b62 = 175.0_pr / 512.0_pr
    real (pr), parameter :: b63 = 575.0_pr / 13824.0_pr
    real (pr), parameter :: b64 = 44275.0_pr / 110592.0_pr
    real (pr), parameter :: b65 = 253.0_pr / 4096.0_pr
    real (pr), parameter :: c1  = 37.0_pr / 378.0_pr
    real (pr), parameter :: c2  = 0.0_pr
    real (pr), parameter :: c3  = 250.0_pr / 621.0_pr
    real (pr), parameter :: c4  = 125.0_pr / 594.0_pr
    real (pr), parameter :: c5  = 0.0_pr
    real (pr), parameter :: c6  = 512.0_pr / 1771.0_pr
    real (pr), parameter :: dc1 = c1 - 2825.0_pr / 27648.0_pr
    real (pr), parameter :: dc2 = c2 - 0.0_pr
    real (pr), parameter :: dc3 = c3 - 18575.0_pr / 48384.0_pr
    real (pr), parameter :: dc4 = c4 - 13525.0_pr / 55296.0_pr
    real (pr), parameter :: dc5 = c5 - 277.0_pr / 14336.0_pr
    real (pr), parameter :: dc6 = c6 - 0.25_pr
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: old
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: new, err
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk1
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk2
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk3
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk4
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk5
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: kk6
    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: tmp
    integer :: j, k

    call get_rhs(old,kk1(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk1(:,j,k) = dt*kk1(:,j,k)
        tmp(:,j,k) = old(:,j,k)+b21*kk1(:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk2(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk2(:,j,k) = dt*kk2(:,j,k)
        tmp(:,j,k) = old(:,j,k)+b31*kk1(:,j,k)+&
                                b32*kk2(:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk3(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk3(:,j,k) = dt*kk3(:,j,k)
        tmp(:,j,k) = old(:,j,k)+b41*kk1(:,j,k)+&
                                b42*kk2(:,j,k)+&
                                b43*kk3(:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk4(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk4(:,j,k) = dt*kk4(:,j,k)
        tmp(:,j,k) = old(:,j,k)+b51*kk1(:,j,k)+&
                                b52*kk2(:,j,k)+&
                                b53*kk3(:,j,k)+&
                                b54*kk4(:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk5(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk5(:,j,k) = dt*kk5(:,j,k)
        tmp(:,j,k) = old(:,j,k)+b61*kk1(:,j,k)+&
                                b62*kk2(:,j,k)+&
                                b63*kk3(:,j,k)+&
                                b64*kk4(:,j,k)+&
                                b65*kk5(:,j,k)
      end do
    end do
    call send_recv_z(tmp)
    call pack_y(tmp)
    call send_recv_y()
    call unpack_y(tmp)
    call get_rhs(tmp, kk6(:,js:je,ks:ke))
    do k=ks,ke
      do j=js,je
        kk6(:,j,k) = dt*kk6(:,j,k)
      end do
    end do

    do k=ks,ke
      do j=js,je
        new(:,j,k) = old(:,j,k) + c1*kk1(:,j,k) + c2*kk2(:,j,k) + &
                                  c3*kk3(:,j,k) + c4*kk4(:,j,k) + &
                                  c5*kk5(:,j,k) + c6*kk6(:,j,k)

        err(:,j,k) = dc1*kk1(:,j,k) + dc2*kk2(:,j,k) + &
                     dc3*kk3(:,j,k) + dc4*kk4(:,j,k) + &
                     dc5*kk5(:,j,k) + dc6*kk6(:,j,k)
      end do
    end do

    return
  end subroutine rkck
  
! ***************************************************************************  

  subroutine get_rhs(in_var, rhs)
    ! Get the right-hand-side of the equation
    use parameters
    use variables
    use derivs
    use ic, only : x, y, z, wall, sphere, sphere2, Vtrap
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: rhs
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: dpsidx, dpsidz
    real (pr), dimension(0:nx1,js:je,ks:ke) :: diss
    real (pr) :: a, b, c
    integer :: i, j, k
    
    a = x(nx-3)
    b = y(ny-3)
    c = z(nz-3)
    
    do k=ks,ke
      do j=js,je
        do i=0,nx1
          diss(i,j,k) = diss_amp*((1.0_pr+tanh(x(i)-a)*tanh(x(i)+a)) + &
                                  (1.0_pr+tanh(y(j)-b)*tanh(y(j)+b)) + &
                                  (1.0_pr+tanh(z(k)-c)*tanh(z(k)+c)))
        end do
      end do
    end do
      
    call deriv_x(in_var, dpsidx)
    !call deriv_z(in_var, dz)
    
    select case (eqn_to_solve)
      case (1) !-2i*dpsi/dt + 2iU*dpsi/dx = del^2(psi) + (1-|psi|^2)psi
        rhs = 0.5_pr*(eye+diss) * ( laplacian(in_var) + &
          (1.0_pr-abs(in_var(:,js:je,ks:ke))**2) * &
          in_var(:,js:je,ks:ke) ) + Urhs*dpsidx

      case (2) !i*dpsi/dt = -del^2(psi) + |psi|^2*psi
        rhs = eye * ( laplacian(in_var) - &
          (abs(in_var(:,js:je,ks:ke))**2)*&
          in_var(:,js:je,ks:ke) ) + Urhs*dpsidx

      case (3) !case(1)*sphere()
        rhs = (0.5_pr*(eye+diss) * ( laplacian(in_var) + &
          (1.0_pr-abs(in_var(:,js:je,ks:ke))**2)*&
          in_var(:,js:je,ks:ke) ) + Urhs*dpsidx) * sphere()

      case (4)
        !i*dpsi/dt = -0.5_pr*del^2(psi) + Vtrap*psi + g(|psi|^2)psi - mu*psi
        rhs = -eye * ( -0.5_pr*laplacian(in_var) + &
          Vtrap()*in_var(:,js:je,ks:ke) + &
          g*(abs(in_var(:,js:je,ks:ke)**2) * &
          in_var(:,js:je,ks:ke)) - &
          mu*in_var(:,js:je,ks:ke) )
    end select

    return
  end subroutine get_rhs
    
! ***************************************************************************  

end module solve
