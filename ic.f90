module ic
  ! Routines to do with setting up the initial condition
  use parameters
  implicit none

  private
  public :: get_grid, ics, get_bcs
  real, dimension(0:nx1), public :: x
  real, dimension(0:ny1), public :: y
  real, dimension(0:nz1), public :: z

  contains

  subroutine get_grid()
    ! Define the space mesh on a box from xl to xr
    use parameters
    implicit none

    integer :: i, j, k

    do i=0,nx1
      x(i) = xl + real(i)*dx
    end do

    do j=0,ny1
      y(j) = yl + real(j)*dy
    end do
    
    do k=0,nz1
      z(k) = zl + real(k)*dz
    end do

    return
  end subroutine get_grid

  subroutine ics(out_var, p)
    ! Set up the initial conditions
    use parameters
    implicit none

    integer, intent(out) :: p
    complex, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: out_var
    logical :: state_exist
    real, dimension(0:nx1,0:ny1,0:nz1) :: real_rand, imag_rand

    ! If this is a restarted run...
    if (restart) then
      real_time = .true.
      inquire(file='end_state.dat', exist=state_exist)
      ! Exit if doing restart but end_state.dat does not exist
      if (.not. state_exist) stop 'ERROR: restart=.true.&
                                   &but end_state.dat does not exist.'
      print*, 'Getting restart conditions'
      ! Get saved data since this is a restart
      call state_restart(out_var, p)
    else
      ! Not a restart so define an initial condition
      !out_var = cmplx(fermi(),0.0)
      !out_var = vortex_pair()
      !out_var = vortex_ring(vr1%x0, vr1%r0)
      out_var = vortex_line(vl1) * &
                vortex_ring(vr1%x0, vr1%r0) * &
                vortex_ring(vr2%x0, vr2%r0)
      !out_var = pade_pulse_ring('pulse', vr1%x0, vr1%r0) * &
      !          pade_pulse_ring('pulse', vr2%x0, vr2%r0)
      !out_var = vortex_line(vl1) * &
      !          pade_pulse_ring('pulse', vr%x0, vr%r0)
      !out_var = pade_pulse_ring('ring', vr%x0, vr%r0)
      !out_var = pade_pulse_ring('pulse', vr%x0, vr%r0)
      !out_var = vortex_line(vl1) * &
      !          vortex_line(vl3)
      !          vortex_line(vl3) * &
      !          vortex_line(vl4)
    end if
  
    return
  end subroutine ics
  
  subroutine state_restart(out_var, p)
    ! Get restart data
    use parameters
    implicit none

    integer, intent(out) :: p
    complex, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: out_var
    integer :: nx_prev, ny_prev, nz_prev, alloc_err
    real :: dt_prev
    complex, dimension(:,:,:), allocatable :: var_prev

    open (50, file = 'end_state.dat', form='unformatted')

    ! Read in the grid sizes and time step from the previous run
    read (50) nx_prev
    read (50) ny_prev
    read (50) nz_prev
    read (50) p
    read (50) t
    read (50) dt_prev
    
    ! If the previous grid size is not the same as the new one then interpolate
    ! onto new grid
    if ((nx_prev /= nx) .or. (ny_prev /= ny) .or. (nz_prev /= nz)) then
      print*, 'Interpolating onto new grid...'
      allocate(var_prev(0:nx_prev-1,0:ny_prev-1,0:nz_prev-1), stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating allocation error' 
      read (50) var_prev
      call inter(var_prev, nx_prev, ny_prev, nz_prev, out_var)
      if (allocated(var_prev)) deallocate(var_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
    else
      ! Otherwise just read in the data
      read (50) out_var
    end if

    close (50)

    select case (scheme)
      case ('rk_adaptive')
        ! Get the previous time step for adaptive Runge-Kutta
        dt = dt_prev
      case default
        ! If restart time step is not the same as that of the previous run then
        ! adjust the time step index p
        if (dt_prev /= dt) then
          p = p * dt_prev / dt
        end if
    end select

    return
  end subroutine state_restart
  
  subroutine inter(prev_var, nxp, nyp, nzp, out_var)
    ! Setup for interpolating onto new grid
    use parameters
    implicit none

    integer, intent(in) :: nxp, nyp, nzp
    complex, dimension(0:nxp-1,0:nyp-1,0:nzp-1), intent(in) :: prev_var
    complex, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: out_var
    integer:: i, j, k
    real :: dx_prev, dy_prev, dz_prev
    real, dimension(0:nxp-1) :: x_prev
    real, dimension(0:nyp-1) :: y_prev
    real, dimension(0:nzp-1) :: z_prev

    ! Previous grid
    dx_prev = (xr-xl) / nxp
    dy_prev = (yr-yl) / nyp
    dz_prev = (zr-zl) / nzp

    ! Previous coordinates
    do i=0,nxp-1
      x_prev(i) = xl + real(i) * dx_prev
    end do

    do j=0,nyp-1
      y_prev(j) = yl + real(j) * dy_prev
    end do
    
    do k=0,nzp-1
      z_prev(k) = zl + real(k) * dz_prev
    end do

    ! Call the interpolation routine
    call inter_var(prev_var, x_prev, y_prev, z_prev, nxp, nyp, nzp, out_var)

    return
  end subroutine inter

  subroutine inter_var(in_var, x_prev, y_prev, z_prev, nxp, nyp, nzp, out_var)
    ! Bilinearly interpolate onto a new grid
    ! See Numerical Recipes in Fortran77 Chap. 3.6 p.116
    use parameters
    implicit none
    
    integer, intent(in) :: nxp, nyp, nzp
    complex, dimension(0:nxp-1,0:nyp-1,0:nzp-1), intent(in)  :: in_var
    real, dimension(0:nxp-1), intent(in) :: x_prev
    real, dimension(0:nyp-1), intent(in) :: y_prev
    real, dimension(0:nzp-1), intent(in) :: z_prev
    complex, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: out_var
    integer :: i, j, k, i2, j2, k2, i2plus, j2plus, k2plus
    real :: c1, c2, c3, one_c1, one_c2, one_c3

    do k=0,nz1   !new 'z' index
      k2 = int(nzp*k/nz)   !old 'z' index
      do j=0,ny1   !new 'y' index
        j2 = int(nyp*j/ny)   !old 'y' index
        do i=0,nx1   !new 'x' index
          i2 = int(nxp*i/nx)   !old 'x' index
          if (i == nx1) then
            i2plus = 0
          else 
            i2plus = i2+1
          end if
          if (j == ny1) then
            j2plus = 0
          else 
            j2plus = j2+1
          end if
          if (k == nz1) then
            k2plus = 0
          else 
            k2plus = k2+1
          end if
          ! Interpolating constants
          c1 = (x(i)-x_prev(i2)) / (x_prev(i2plus)-x_prev(i2))
          c2 = (y(j)-y_prev(j2)) / (y_prev(j2plus)-y_prev(j2))
          c3 = (z(k)-z_prev(k2)) / (z_prev(k2plus)-z_prev(k2))
          one_c1 = 1.0-c1
          one_c2 = 1.0-c2
          one_c3 = 1.0-c3
          ! Do the interpolation
          out_var(i,j,k) = one_c1*one_c2*one_c3*in_var(i2,j2,k2) + &
                           c1*one_c2*one_c3*in_var(i2+1,j2,k2) + &
                           c1*c2*one_c3*in_var(i2+1,j2+1,k2) + &
                           one_c1*c2*one_c3*in_var(i2,j2+1,k2) + &
                           one_c1*one_c2*c3*in_var(i2,j2,k2+1) + &
                           c1*one_c2*c3*in_var(i2+1,j2,k2+1) + &
                           c1*c2*c3*in_var(i2+1,j2+1,k2+1) + &
                           one_c1*c2*c3*in_var(i2,j2+1,k2+1)
        end do
      end do
    end do

    return
  end subroutine inter_var

  function fermi()
    ! Thomas-Fermi initial condition
    use parameters
    implicit none

    real, parameter :: c=450.0
    real, dimension(0:nx1,0:ny1,0:nz1) :: fermi, r
    real :: mu
    integer :: i, j, k

    mu = sqrt(c/pi)

    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          r(i,j,k) = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
          if (r(i,j,k)<=sqrt(2.0*mu)) then
            fermi(i,j,k) = sqrt((2.0*mu - r(i,j,k)**2)/(2.0*c))
          else
            fermi(i,j,k) = 0.0
          end if
        end do
      end do
    end do
    
    return
  end function fermi

  function ei(theta)
    ! exp(i*theta) where theta is the argument arctan(y/x)
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: ei
    real, dimension(0:nx1,0:ny1,0:nz1), intent(in) :: theta

    ei = cos(theta) + eye*sin(theta)

    return
  end function ei
  
  function amp(r)
    ! Amplitude of a vortex line
    use parameters
    implicit none

    real, dimension(0:nx1,0:ny1,0:nz1) :: amp
    real, dimension(0:nx1,0:ny1,0:nz1), intent(in) :: r
    real, parameter :: c1 = -0.7
    real, parameter :: c2 = 1.15

    amp = 1.0 - exp(c1*r**c2)

    return
  end function amp
  
  function vortex_line(vl)
    ! Vortex line initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: vortex_line
    type (line_param), intent(in) :: vl
    real, dimension(0:nx1,0:ny1,0:nz1) :: r, theta

    call get_r(vl%x0, vl%y0, vl%amp, vl%ll, r)
    call get_theta(vl%x0, vl%y0, vl%amp, vl%ll, vl%sgn, theta)

    vortex_line = amp(r)*ei(theta)

    return
  end function vortex_line
  
  function vortex_ring(x0, r0)
    ! Vortex ring initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: vortex_ring
    real, intent(in) :: x0, r0
    real, dimension(0:ny1,0:nz1) :: s
    real, dimension(0:nx1,0:ny1,0:nz1) :: rr1, rr2, d1, d2
    integer :: i, j, k

    call get_s(s)
    
    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          d1(i,j,k) = sqrt( (x(i)-x0)**2 + (s(j,k)+r0)**2 )
          d2(i,j,k) = sqrt( (x(i)-x0)**2 + (s(j,k)-r0)**2 )
        end do
      end do
    end do
    
    call get_rr(d1,rr1)
    call get_rr(d2,rr2)
    
    rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
                 (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
    rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
                 (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )

    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          vortex_ring(i,j,k) = rr1(i,j,k)*((x(i)-x0)+eye*(s(j,k)+r0)) * &
                               rr2(i,j,k)*((x(i)-x0)-eye*(s(j,k)-r0))
        end do
      end do
    end do

    return
  end function vortex_ring
  
  function pade_pulse_ring(pulse_or_ring, x0, r0)
    ! Pade approximation vortex ring or pulse initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: pade_pulse_ring
    real, intent(in) :: x0, r0
    character(*), intent(in) :: pulse_or_ring
    real, dimension(0:nx1) :: x2
    real, dimension(0:ny1) :: y2
    real, dimension(0:ny1,0:nz1) :: s, s2
    real, dimension(0:nx1,0:ny1,0:nz1) :: uu, vv, denom
    real, dimension(9) :: a
    real, parameter :: pow = 7.0/4.0
    real :: one_2u, U, m
    integer :: i, j, k

    call get_s(s)
    call get_consts(pulse_or_ring, a, U, m)
    
    one_2u = 1.0-2.0*U**2
    x2 = (x-x0)**2
    y2 = (y)**2
    s2 = (s-r0)**2
    
    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          denom(i,j,k) = (1.0 + a(8)*x2(i) + a(7)*s2(j,k) + &
                          a(9)*(x2(i)+one_2u*s2(j,k))**2)**pow
          
          uu(i,j,k) = 1.0 + ( a(1) + a(3)*x2(i) + a(2)*s2(j,k) + &
                              m*(a(9)**pow)*U*(2.0*x2(i) - &
                              one_2u*s2(j,k)) ) * (x2(i)+one_2u*s2(j,k)) / &
                              denom(i,j,k)

          vv(i,j,k) = (x(i)-x0) * ( a(4) + a(6)*x2(i) + a(5)*s2(j,k) - &
                                    m*(a(9)**pow)*(x2(i) + &
                                    one_2u*s2(j,k))**2 ) / denom(i,j,k)
        end do
      end do
    end do
    
    pade_pulse_ring = uu + eye*vv

    contains

    subroutine get_consts(flag, consts, vel, mom)
      ! Get constants required to set up the vortex ring or pulse
      implicit none

      character(*), intent(in) :: flag
      real, dimension(9), intent(out) :: consts
      real, intent(out) :: vel, mom

      select case (flag)
        case ('ring')
          consts = (/-1.1,0.0170524,0.0289452,&
                    -0.953,-0.0049767,-0.0594346,&
                    0.04,0.21,0.00612/)
          vel = 0.6
          mom = 8.97
        case ('pulse')
          consts = (/-0.79792,0.00388059,0.00882276,&
                     -0.7981,-0.012783,-0.0574092,&
                     0.0399,0.199,0.0058/)
          vel = 0.63
          mom = 8.37
      end select

      return
    end subroutine get_consts

  end function pade_pulse_ring
  
  function vortex_pair()
    ! Vortex pair initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: vortex_pair
    real, dimension(0:nx1) :: x2, x4
    real, dimension(0:ny1) :: y2, y4
    real, dimension(0:nx1,0:ny1) :: uu, vv, denom
    
    real, parameter :: a0 = -1.14026
    real, parameter :: a1 = -0.150112
    real, parameter :: a2 = -0.0294564
    real, parameter :: b0 = -0.830953
    real, parameter :: b1 = -0.108296
    real, parameter :: b2 = -0.073641
    real, parameter :: c0 = 0.35022
    real, parameter :: c1 = 0.03032
    real, parameter :: c2 = 0.15905
    real, parameter :: c3 = 0.04123
    real, parameter :: c4 = 0.01402

    integer :: i, j, k

    x2 = x**2
    x4 = x**4
    y2 = y**2
    y4 = y**4
    
    do j=0,ny1
      do i=0,nx1
        denom(i,j) = 1.0 + c0*x2(i) + c1*x4(i) + c2*y2(j) + &
                           c3*x2(i)*y2(j) + c4*y(j)**4
        
        uu(i,j) = 1.0 + (a0 + a1*x2(i) + a2*y2(j)) / denom(i,j)

        vv(i,j) = x(i) * (b0 + b1*x2(i) + b2*y2(j)) / denom(i,j)
      end do
    end do
    
    do k =0,nz1
      vortex_pair(:,:,k) = uu(:,:) + eye*vv(:,:)
    end do

    return
  end function vortex_pair

  subroutine get_r(x0, y0, a, ll, r)
    ! Get the cylindrical-polar radius r**2=x**2+y**2
    use parameters
    implicit none

    real, intent(in) :: x0, y0, a, ll
    real, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: r
    integer :: i, j, k

    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          r(i,j,k) = sqrt((x(i)-x0)**2 + (y(j)-y0-a*cos(2.0*pi*z(k)/ll))**2)
        end do
      end do
    end do

    return
  end subroutine get_r
  
  subroutine get_s(s)
    ! Another radial variable
    use parameters
    implicit none

    real, dimension(0:ny1,0:nz1), intent(out) :: s
    integer :: j, k

    do k=0,nz1
      do j=0,ny1
        s(j,k) = sqrt(y(j)**2 + z(k)**2)
      end do
    end do

    return
  end subroutine get_s

  subroutine get_theta(x0, y0, a, ll, sgn, theta)
    ! Get the argument theta=arctan(y/x)
    use parameters
    implicit none

    real, intent(in) :: x0, y0, a, ll, sgn
    real, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: theta
    integer :: i, j, k

    do k=0,nz1
      do j=0,ny1
        do i=0,nx1
          theta(i,j,k) = sgn * atan2(y(j)-y0-a*cos(2.0*pi*z(k)/ll), x(i)-x0)
        end do
      end do
    end do

    return
  end subroutine get_theta

  subroutine get_rr(r,rr)
    ! R in psi=R(r)exp(i*theta)
    use parameters
    implicit none

    real, dimension(0:nx1,0:ny1,0:nz1), intent(in) :: r
    real, dimension(0:nx1,0:ny1,0:nz1), intent(out) :: rr
    integer :: i, j
    
    rr = sqrt( ((0.3437+0.0286*r**2)) / &
                (1.0+(0.3333*r**2)+(0.0286*r**4)) )

    return
  end subroutine get_rr
    
  subroutine get_bcs(in_var)
    ! Zero boundary conditions
    use parameters
    implicit none

    complex, dimension(0:nx1,0:ny1,0:nz1) :: in_var

    in_var(0:2,:,:) = 0.0
    in_var(:,0:2,:) = 0.0
    in_var(:,:,0:2) = 0.0
    in_var(nx-3:nx1,:,:) = 0.0
    in_var(:,ny-3:ny1,:) = 0.0
    in_var(:,:,nz-3:nz1) = 0.0

    return
  end subroutine get_bcs

end module ic
