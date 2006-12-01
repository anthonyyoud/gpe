! $Id: ic.f90,v 1.43 2006-12-01 12:52:14 n8049290 Exp $
!----------------------------------------------------------------------------

module ic
  ! Routines to do with setting up the initial condition
  use parameters
  implicit none

  private
  public :: get_grid, ics, fft, get_kc_amp, wall, sphere

  real, dimension(0:nx1), public :: x
  real, dimension(0:ny1), public :: y
  real, dimension(0:nz1), public :: z

  contains

  subroutine get_grid()
    ! Define the space mesh on a box from xl to xr.  x, y, and z are defined in
    ! full on all processes because they are not 2D arrays.  If they are
    ! distributed and any calculations are done with them (such as a sum) then
    ! there are contributions from the processes over which x, y, z are NOT
    ! distributed leading to erroneous results.
    use parameters
    implicit none

    integer :: i, j, k

    ! x-coordinate
    do i=0,nx1
      x(i) = xl + real(i)*dx
    end do

    ! y-coordinate
    do j=0,ny1
      y(j) = yl + real(j)*dy
    end do
    
    ! z-coordinate
    do k=0,nz1
      z(k) = zl + real(k)*dz
    end do

    return
  end subroutine get_grid

! ***************************************************************************  

  subroutine ics(out_var, p)
    ! Set up the initial conditions
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    integer,                                       intent(out) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend)              :: tmp_var
    logical :: state_exist

    ! If this is a restarted run...
    if (restart) then
      real_time = .true.
      inquire(file=proc_dir//'end_state.dat', exist=state_exist)
      ! Exit if doing restart but end_state.dat does not exist
      if (.not. state_exist) stop 'ERROR: restart=.true.&
                                   &but end_state.dat does not exist.'
      if (myrank == 0) then
        print*, 'Getting restart conditions'
      end if
      ! Get saved data since this is a restart
      call state_restart(tmp_var, p)
      out_var = tmp_var !*vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir)
      !out_var = tmp_var*vortex_line(vl1)
    else
      ! Not a restart so define an initial condition
      !out_var = cmplx(fermi(),0.0)
      !out_var = vortex_pair()
      !out_var = vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir) !* &
      !          vortex_ring(vr2%x0, vr2%y0, vr2%z0, vr2%r0, vr2%dir) * &
      !          vortex_ring2(vr3%x0, vr3%y0, vr3%z0, vr3%r0, vr3%dir) * &
      !          vortex_ring2(vr4%x0, vr4%y0, vr4%z0, vr4%r0, vr4%dir)
      !out_var = vortex_line(vl1) * &
      !          vortex_ring(vr1%x0, vr1%z0, vr1%r0) * &
      !          vortex_ring(vr2%x0, vr2%z0, vr2%r0)
      !out_var = pade_pulse_ring('pulse', vr1%x0, vr1%y0, vr1%z0, vr1%r0) * &
      !          pade_pulse_ring('pulse', vr2%x0, vr2%y0, vr2%z0, vr2%r0)
      !out_var = vortex_line(vl1) * &
      !          pade_pulse_ring('pulse', vr1%x0, vr1%y0, vr1%z0, vr1%r0)
      !out_var = pade_pulse_ring('ring', vr1%x0, vr1%y0, vr1%z0, vr1%r0)
      !out_var = pade_pulse_ring('pulse', vr1%x0, vr1%y0, vr1%z0, vr1%r0)
      !out_var = vortex_line(vl1) * &
      !          vortex_line(vl2) * &
      !          vortex_line(vl3) !* &
      !          vortex_line(vl4)
      !out_var = sphere() !* &
                !vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir)
      !out_var = wall() !* &
      !          vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir)
      call random_phase(tmp_var)
      out_var = tmp_var !* vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir)
      !out_var = vortex_ring(vr1%x0, vr1%y0, vr1%z0, vr1%r0, vr1%dir) * &
      !          vortex_ring(vr2%x0, vr2%y0, vr2%z0, vr2%r0, vr2%dir) * &
      !          vortex_ring(vr3%x0, vr3%y0, vr3%z0, vr3%r0, vr3%dir) !* &
      !          vortex_ring(vr4%x0, vr4%y0, vr4%z0, vr4%r0, vr4%dir) * &
      !          vortex_ring(vr5%x0, vr5%y0, vr5%z0, vr5%r0, vr5%dir) * &
      !          vortex_ring(vr6%x0, vr6%y0, vr6%z0, vr6%r0, vr6%dir) * &
      !          vortex_ring(vr7%x0, vr7%y0, vr7%z0, vr7%r0, vr7%dir) !* &
      !          vortex_line(vl1)
    end if
  
    return
  end subroutine ics
  
! ***************************************************************************  

  subroutine state_restart(out_var, p)
    ! Get restart data
    use parameters
    use variables
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    integer,                                       intent(out) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: out_var1, out_var2
    complex :: dt_prev
    integer :: nx_prev, ny_prev, nz_prev, undef_int
    real :: undef_real

    out_var2 = 1.0
    
    open (unit_no, file=proc_dir//'end_state.dat', form='unformatted')

    ! Read in the distributed data from the previous run
    read (unit_no) nx_prev
    read (unit_no) ny_prev
    read (unit_no) nz_prev
    read (unit_no) p
    read (unit_no) t
    read (unit_no) dt_prev
    read (unit_no) out_var1

    close (unit_no)
    
    if (saved_restart) then
      open (unit_no, file=proc_dir//'end_state_filtered.dat', &
                     form='unformatted')

      ! Read in the filtered distributed data from the previous run
      read (unit_no) undef_int
      read (unit_no) undef_int
      read (unit_no) undef_int
      read (unit_no) undef_int
      read (unit_no) undef_real
      read (unit_no) undef_real
      read (unit_no) out_var2

      close (unit_no)
    end if

    out_var = out_var1*out_var2

    if (real(dt_prev) == 0.0) then
      dt_prev = cmplx(aimag(dt_prev), real(dt_prev))
    end if

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

! ***************************************************************************  

  subroutine random_phase(out_var)
    ! Strongly nonequilibrium initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: a
    integer, dimension(1) :: seed
    real :: phi, rand
    integer :: i, j, k, ii2, jj2, kk2, irand
    
    if (myrank == 0) then
      write (97, *) 'Complex amplitude=', comp_amp
      write (97, *) 'Cutoff wavenumber=', kc2
    end if
    
    if (myrank == 0) then
      call random_seed()
      call random_number(rand)
      irand = int(rand*100000)
    end if

    call MPI_BCAST(irand, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    seed(1) = myrank+irand
    call random_seed(put=seed)

    do k=ksta,kend
      if (k <= nz1/2+1) then
        kk2 = k**2
      else
        kk2 = (nz1-k+1)**2
      end if
      do j=jsta,jend
        if (j <= ny1/2+1) then
          jj2 = j**2
        else
          jj2 = (ny1-j+1)**2
        end if
        do i=0,nx1
          if (i <= nx1/2+1) then
            ii2 = i**2
          else
            ii2 = (nx1-i+1)**2
          end if
          if ((ii2 + jj2 + kk2) <= kc2) then
            call random_number(phi)
            a(i,j,k) = comp_amp*exp(2.0*pi*eye*phi)
          else
            a(i,j,k) = 0.0
          end if
        end do
      end do
    end do

    call fft(a, out_var, 'forward', .false.)
  
    return
  end subroutine random_phase

! ***************************************************************************  

  subroutine get_kc_amp()
    ! Get the cutoff wavenumber and amplitude of the random phase IC
    use parameters
    implicit none
    
    real :: ev
    
    ev = ((0.5*nx/pi)**2)*enerv
    comp_amp = sqrt(((0.11095279993106999*nv**2.5)*(nx*ny*nz))/(ev**1.5))
    kc2 = (1.666666666666666*ev)/nv

    return
  end subroutine get_kc_amp

! ***************************************************************************  

  function fermi()
    ! Thomas-Fermi initial condition
    use parameters
    implicit none

    real, parameter                            :: c=450.0
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: fermi, r
    real                                       :: mu
    integer                                    :: i, j, k

    mu = sqrt(c/pi)

    do k=ksta,kend
      do j=jsta,jend
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

! ***************************************************************************  

  function ei(theta)
    ! exp(i*theta) where theta is the argument arctan(y/x)
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend)             :: ei
    real,    dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: theta
    integer                                                   :: j, k

    ei = cos(theta) + eye*sin(theta)

    return
  end function ei
  
! ***************************************************************************  

  function amp(r)
    ! Amplitude of a vortex line
    use parameters
    implicit none

    real, dimension(0:nx1,jsta:jend,ksta:kend)             :: amp
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: r
    real, parameter                                        :: c1 = -0.7
    real, parameter                                        :: c2 = 1.15
    integer                                                :: j, k

    amp = 1.0 - exp(c1*r**c2)

    return
  end function amp
  
! ***************************************************************************  

  function vortex_line(vl)
    ! Vortex line initial condition
    use parameters
    implicit none

    complex,           dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_line
    type (line_param), intent(in)                           :: vl
    real,              dimension(0:nx1,jsta:jend,ksta:kend) :: r, theta
    integer                                                 :: j, k

    call get_r(vl%x0, vl%y0, vl%amp, vl%ll, r)
    call get_theta(vl%x0, vl%y0, vl%amp, vl%ll, vl%sgn, theta)

    vortex_line = amp(r)*ei(theta)
    
    return
  end function vortex_line
  
! ***************************************************************************  

  function vortex_ring(x0, y0, z0, r0, dir)
    ! Vortex ring initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_ring
    real,    intent(in)                           :: x0, y0, z0, r0, dir
    real,    dimension(jsta:jend,ksta:kend)       :: s
    real,    dimension(0:nx1,jsta:jend,ksta:kend) :: rr1, rr2, d1, d2
    integer                                       :: i, j, k

    call get_s(s, y0, z0)
    
    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          d1(i,j,k) = sqrt( (scal*(x(i)-x0))**2 + (scal*(s(j,k)+r0))**2 )
          d2(i,j,k) = sqrt( (scal*(x(i)-x0))**2 + (scal*(s(j,k)-r0))**2 )
        end do
      end do
    end do
    
    call get_rr(d1,rr1)
    call get_rr(d2,rr2)
    
    rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
                        (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
    rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
                        (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )

    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          vortex_ring(i,j,k) = rr1(i,j,k)*((scal*(x(i)-x0))+dir*eye*(scal*(s(j,k)+r0))) * &
                               rr2(i,j,k)*((scal*(x(i)-x0))-dir*eye*(scal*(s(j,k)-r0)))
        end do
      end do
    end do

    return
  end function vortex_ring
  
! ***************************************************************************  

  function vortex_ring2(x0, y0, z0, r0, dir)
    ! Vortex ring initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_ring2
    real,    intent(in)                           :: x0, y0, z0, r0, dir
    real,    dimension(0:nx1,ksta:kend)           :: s
    real,    dimension(0:nx1,jsta:jend,ksta:kend) :: rr1, rr2, d1, d2
    integer                                       :: i, j, k

    do k=ksta,kend
      do i=0,nx1
        s(i,k) = sqrt((x(i)-x0)**2 + (z(k)-z0)**2)
      end do
    end do
    
    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          d1(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)+r0)**2 )
          d2(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)-r0)**2 )
        end do
      end do
    end do
    
    rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
                        (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
    rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
                        (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )

    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          vortex_ring2(i,j,k) = rr1(i,j,k)*((y(j)-y0)+dir*eye*(s(i,k)+r0)) * &
                                rr2(i,j,k)*((y(j)-y0)-dir*eye*(s(i,k)-r0))
        end do
      end do
    end do

    return
  end function vortex_ring2
  
! ***************************************************************************  

  function pade_pulse_ring(pulse_or_ring, x0, y0, z0, r0)
    ! Pade approximation vortex ring or pulse initial condition
    use parameters
    implicit none

    complex,      dimension(0:nx1,jsta:jend,ksta:kend) :: pade_pulse_ring
    real,         intent(in)                           :: x0, y0, z0, r0
    character(*), intent(in)                           :: pulse_or_ring
    real,         dimension(0:nx1,jsta:jend,ksta:kend) :: uu, vv, denom
    real,         dimension(jsta:jend,ksta:kend)       :: s, s2
    real,         dimension(0:nx1)                     :: x2
    real,         dimension(9)                         :: a
    real,         parameter                            :: pow = 7.0/4.0
    real                                               :: one_2u, U, m
    integer                                            :: i, j, k

    call get_s(s, y0, z0)
    call get_consts(pulse_or_ring, a, U, m)
    
    one_2u = 1.0-2.0*U**2
    x2 = (x-x0)**2
    s2 = (s-r0)**2
    
    do k=ksta,kend
      do j=jsta,jend
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

      character(*),               intent(in)  :: flag
      real,         dimension(9), intent(out) :: consts
      real,                       intent(out) :: vel, mom

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
  
! ***************************************************************************  

  function vortex_pair()
    ! Vortex pair initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_pair
    real,    dimension(0:nx1,jsta:jend)           :: uu, vv, denom
    real,    dimension(0:nx1)                     :: x2, x4
    real,    dimension(0:ny1)                     :: y2, y4
    real,    parameter                            :: a0 = -1.14026
    real,    parameter                            :: a1 = -0.150112
    real,    parameter                            :: a2 = -0.0294564
    real,    parameter                            :: b0 = -0.830953
    real,    parameter                            :: b1 = -0.108296
    real,    parameter                            :: b2 = -0.073641
    real,    parameter                            :: c0 = 0.35022
    real,    parameter                            :: c1 = 0.03032
    real,    parameter                            :: c2 = 0.15905
    real,    parameter                            :: c3 = 0.04123
    real,    parameter                            :: c4 = 0.01402
    integer                                       :: i, j, k

    x2 = x**2
    x4 = x**4
    y2 = y**2
    y4 = y**4
    
    do j=jsta,jend
      do i=0,nx1
        denom(i,j) = 1.0 + c0*x2(i) + c1*x4(i) + c2*y2(j) + &
                           c3*x2(i)*y2(j) + c4*y(j)**4
        
        uu(i,j) = 1.0 + (a0 + a1*x2(i) + a2*y2(j)) / denom(i,j)

        vv(i,j) = x(i) * (b0 + b1*x2(i) + b2*y2(j)) / denom(i,j)
      end do
    end do
    
    do k=ksta,kend
      do j=jsta,jend
        vortex_pair(:,j,k) = uu(:,j) + eye*vv(:,j)
      end do
    end do

    return
  end function vortex_pair

! ***************************************************************************  

  function sphere()
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: sphere
    real, parameter :: rad = 10.0
    real, parameter :: eps = 2.0
    integer :: i, j, k

    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          sphere(i,j,k) = 0.5*(1.0+tanh(sqrt(x(i)**2+y(j)**2+z(k)**2)-rad-eps))
          !sphere(i,j,k) = max(0.5*(1.0+&
          !                         tanh(sqrt(x(i)**2+y(j)**2+z(k)**2)-rad)-&
          !                         eps),0.0)
        end do
      end do
    end do

    return
  end function sphere
  
! ***************************************************************************  

  function wall()
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: wall
    real, parameter :: pos = -30.0
    integer :: i, j, k

    !do k=ksta,kend
    !  wall(:,:,k) = 0.5*(1.0+tanh(z(k)-pos))
    !end do
    
    do k=ksta,kend
      if (z(k) <= pos) then
        wall(:,:,k) = 0.0
      else
        wall(:,:,k) = 1.0
      end if
    end do

    return
  end function wall

! ***************************************************************************  

  subroutine get_r(x0, y0, a, ll, r)
    ! Get the cylindrical-polar radius r**2=x**2+y**2
    use parameters
    implicit none

    real,                                       intent(in)  :: x0, y0, a, ll
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: r
    integer                                                 :: i, j, k

    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          r(i,j,k) = sqrt(scal*(x(i)-x0)**2 + &
                         (scal*(y(j)-y0-a*cos(2.0*pi*z(k)/ll))**2))
        end do
      end do
    end do

    return
  end subroutine get_r
  
! ***************************************************************************  

  subroutine get_s(s, y0, z0)
    ! Another radial variable
    use parameters
    implicit none

    real,                                 intent(in)  :: y0, z0
    real, dimension(jsta:jend,ksta:kend), intent(out) :: s
    integer                                           :: j, k

    do k=ksta,kend
      do j=jsta,jend
        s(j,k) = sqrt((y(j)-y0)**2 + (z(k)-z0)**2)
      end do
    end do

    return
  end subroutine get_s

! ***************************************************************************  

  subroutine get_theta(x0, y0, a, ll, sgn, theta)
    ! Get the argument theta=arctan(y/x)
    use parameters
    implicit none

    real, intent(in) :: x0, y0, a, ll, sgn
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: theta
    integer                                                 :: i, j, k

    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          theta(i,j,k) = sgn * atan2(scal*(y(j)-y0-a*cos(2.0*pi*z(k)/ll)), &
                                     scal*(x(i)-x0))
        end do
      end do
    end do

    return
  end subroutine get_theta

! ***************************************************************************  

  subroutine get_rr(r,rr)
    ! R in psi=R(r)exp(i*theta)
    use parameters
    implicit none

    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in)  :: r
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: rr
    integer                                                 :: i, j, k
    
    rr = sqrt( ((0.3437+0.0286*r**2)) / &
                (1.0+(0.3333*r**2)+(0.0286*r**4)) )

    return
  end subroutine get_rr
    
! ***************************************************************************  

  subroutine fft(in_var, out_var, dir, particles)
    ! Calculate the FFT (or inverse) of a variable
    use parameters
    use constants
    use variables, only : unit_no
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    character(*), intent(in) :: dir
    logical, intent(in) :: particles
    complex, dimension(0:nx1,0:ny1,0:nz1) :: tmp_var, tmp
    complex, allocatable, dimension(:)    :: local_data, work
    integer :: i, j, k, plan, iplan
    integer :: loc_nz, loc_z_sta, loc_ny, loc_y_sta, tot_loc

    ! Set up a temporary array so that the data can be eventually redistributed
    ! as required by fftw
    tmp_var = 0.0
    do k=ksta,kend
      do j=jsta,jend
        tmp_var(:,j,k) = in_var(:,j,k)
      end do
    end do 
    
    ! Make sure all processes have the temporary variable
    call MPI_ALLREDUCE(tmp_var, tmp, nx*ny*nz, & 
                       MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
                       
    select case (dir)
      case ('forward')
        ! Create the plan for a forward transform
        call fftw3d_f77_mpi_create_plan(plan, MPI_COMM_WORLD, nx, ny, nz, &
                                        FFTW_FORWARD, FFTW_ESTIMATE)
      case ('backward')
        ! Create the plan for a backward transform
        call fftw3d_f77_mpi_create_plan(plan, MPI_COMM_WORLD, nx, ny, nz, &
                                        FFTW_BACKWARD, FFTW_ESTIMATE)
      case default
        stop 'ERROR: Not a valid direction for FFT!'
    end select

    ! Get the sizes of the arrays local to each process
    call fftwnd_f77_mpi_local_sizes(plan, loc_nz, loc_z_sta, &
                                          loc_ny, loc_y_sta, &
                                          tot_loc)

    ! Allocate the array sizes
    allocate(local_data(0:tot_loc-1))
    allocate(work(0:tot_loc-1))

    ! Distribute the data into a local array as required by fftw
    do k=0,loc_nz-1
      do j=0,ny1
        do i=0,nx1
          local_data((k*ny+j)*nx+i) = tmp(i,j,k+loc_z_sta)
        end do
      end do
    end do

    ! Call the parallel FFT routine to do the transform
    call fftwnd_f77_mpi(plan, 1, local_data, work, 1, FFTW_TRANSPOSED_ORDER)
    
    ! Set up another temporary array to (eventually) get the transformed data
    ! back into the original distribution scheme
    tmp = 0.0
    do k=0,loc_ny-1
      do j=0,nz1
        do i=0,nx1
          tmp(i,k+loc_y_sta,j) = local_data((k*ny+j)*nx+i)
        end do
      end do
    end do
      
    ! Make sure all processes have the temporary variable
    call MPI_ALLREDUCE(tmp, tmp_var, nx*ny*nz, &
                       MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Renormalise
    tmp_var = tmp_var / sqrt(real(nx*ny*nz))

    ! Distribute the transformed data
    do k=ksta,kend
      do j=jsta,jend
        out_var(:,j,k) = tmp_var(:,j,k)
      end do
    end do

    !if (dir == 'backward') then
    !  open (60, file=proc_dir//'fft0000000.dat', form='unformatted')
    !  write (60) nx, ny, nz
    !  write (60) nyprocs, nzprocs
    !  write (60) jsta, jend, ksta, kend
    !  write (60) abs(out_var)**2
    !  write (60) x
    !  write (60) y
    !  write (60) z
    !  close (60)
    !end if
        
    
    ! Deallocate the allocated arrays
    deallocate(local_data)
    deallocate(work)
    
    ! Destroy the transform plan
    call fftwnd_f77_mpi_destroy_plan(plan)

    return
  end subroutine fft
  
! ***************************************************************************  

end module ic
