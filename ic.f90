module ic
  ! Routines to do with setting up the initial condition
  use parameters
  implicit none

  private
  public :: get_grid, ics, fft
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

  subroutine ics(out_var, p)
    ! Set up the initial conditions
    use parameters
    implicit none

    integer, intent(out) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
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
      !out_var = pade_pulse_ring('ring', vr1%x0, vr1%r0)
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
    use variables
    implicit none

    integer, intent(out) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: out_var
    integer :: nx_prev, ny_prev, nz_prev
    real :: dt_prev

    open (unit_no, file=proc_dir//'end_state.dat', form='unformatted')

    ! Read in the distributed data from the previous run
    read (unit_no) nx_prev
    read (unit_no) ny_prev
    read (unit_no) nz_prev
    read (unit_no) p
    read (unit_no) t
    read (unit_no) dt_prev
    read (unit_no) out_var

    close (unit_no)

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
  
  function fermi()
    ! Thomas-Fermi initial condition
    use parameters
    implicit none

    real, parameter :: c=450.0
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: fermi, r
    real :: mu
    integer :: i, j, k

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

  function ei(theta)
    ! exp(i*theta) where theta is the argument arctan(y/x)
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: ei
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: theta
    integer :: j, k

    ei = cos(theta) + eye*sin(theta)

    return
  end function ei
  
  function amp(r)
    ! Amplitude of a vortex line
    use parameters
    implicit none

    real, dimension(0:nx1,jsta:jend,ksta:kend) :: amp
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: r
    real, parameter :: c1 = -0.7
    real, parameter :: c2 = 1.15
    integer :: j, k

    amp = 1.0 - exp(c1*r**c2)

    return
  end function amp
  
  function vortex_line(vl)
    ! Vortex line initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_line
    type (line_param), intent(in) :: vl
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: r, theta
    integer :: j, k

    call get_r(vl%x0, vl%y0, vl%amp, vl%ll, r)
    call get_theta(vl%x0, vl%y0, vl%amp, vl%ll, vl%sgn, theta)

    vortex_line = amp(r)*ei(theta)
    
    return
  end function vortex_line
  
  function vortex_ring(x0, r0)
    ! Vortex ring initial condition
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_ring
    real, intent(in) :: x0, r0
    real, dimension(jsta:jend,ksta:kend) :: s
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: rr1, rr2, d1, d2
    integer :: i, j, k

    call get_s(s)
    
    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          d1(i,j,k) = sqrt( (x(i)-x0)**2 + (s(j,k)+r0)**2 )
          d2(i,j,k) = sqrt( (x(i)-x0)**2 + (s(j,k)-r0)**2 )
        end do
      end do
    end do
    
    call get_rr(d1,rr1)
    call get_rr(d2,rr2)
    
    rr1(:,j,k) = sqrt( ((0.3437+0.0286*d1(:,j,k)**2)) / &
                        (1.0+(0.3333*d1(:,j,k)**2)+(0.0286*d1(:,j,k)**4)) )
    rr2(:,j,k) = sqrt( ((0.3437+0.0286*d2(:,j,k)**2)) / &
                        (1.0+(0.3333*d2(:,j,k)**2)+(0.0286*d2(:,j,k)**4)) )

    do k=ksta,kend
      do j=jsta,jend
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

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: pade_pulse_ring
    real, intent(in) :: x0, r0
    character(*), intent(in) :: pulse_or_ring
    real, dimension(0:nx1) :: x2
    real, dimension(0:ny1) :: y2
    real, dimension(jsta:jend,ksta:kend) :: s, s2
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: uu, vv, denom
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

    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_pair
    real, dimension(0:nx1) :: x2, x4
    real, dimension(0:ny1) :: y2, y4
    real, dimension(0:nx1,jsta:jend) :: uu, vv, denom
    
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

  subroutine get_r(x0, y0, a, ll, r)
    ! Get the cylindrical-polar radius r**2=x**2+y**2
    use parameters
    implicit none

    real, intent(in) :: x0, y0, a, ll
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: r
    integer :: i, j, k

    do k=ksta,kend
      do j=jsta,jend
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

    real, dimension(jsta:jend,ksta:kend), intent(out) :: s
    integer :: j, k

    do k=ksta,kend
      do j=jsta,jend
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
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: theta
    integer :: i, j, k

    do k=ksta,kend
      do j=jsta,jend
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

    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: r
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: rr
    integer :: i, j, k
    
    rr = sqrt( ((0.3437+0.0286*r**2)) / &
                (1.0+(0.3333*r**2)+(0.0286*r**4)) )

    return
  end subroutine get_rr
    
  subroutine fft(in_var)
    use parameters
    use constants
    use variables, only : unit_no
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(inout) :: in_var
    complex, dimension(0:nx1,0:ny1,0:nz1) :: tmp_var, tmp
    complex, allocatable, dimension(:) :: local_data, work
    integer :: i, j, k
    integer :: plan, iplan
    integer :: loc_nz, loc_z_sta, loc_ny, loc_y_sta, tot_loc

    tmp_var = 0.0
    do k=ksta,kend
      do j=jsta,jend
        tmp_var(:,j,k) = in_var(:,j,k)
      end do
    end do 
    
    call MPI_ALLREDUCE(tmp_var, tmp, nx*ny*nz, & 
                       MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
                       
    open (unit_no, file=proc_dir//'orig.dat')
    do j=0,ny1
      write (unit_no, '(3e17.9)') (y(j), z(k), abs(tmp(nx/2,j,k)), k=0,nz1)
      write (unit_no, *)
    end do
    close (unit_no)

    call fftw3d_f77_mpi_create_plan(plan, MPI_COMM_WORLD, nx, ny, nz, &
                                    FFTW_FORWARD, FFTW_ESTIMATE)

    call fftwnd_f77_mpi_local_sizes(plan, loc_nz, loc_z_sta, &
                                          loc_ny, loc_y_sta, &
                                          tot_loc)

    print*, loc_nz, loc_z_sta, loc_ny, loc_y_sta, tot_loc

    allocate(local_data(0:tot_loc-1))
    allocate(work(0:tot_loc-1))

    do k=0,loc_nz-1
      do j=0,ny1
        do i=0,nx1
          local_data((k*ny+j)*nx+i) = tmp(i,j,k+loc_z_sta)
          !print*, i, j, k, (k*ny+j)*nx+i
        end do
      end do
    end do

    call fftwnd_f77_mpi(plan, 1, local_data, work, 1, FFTW_TRANSPOSED_ORDER)
    !call fftwnd_f77_mpi(plan, 1, local_data, work, 1, FFTW_NORMAL_ORDER)
    
    !tmp = 0.0
    !do k=0,loc_ny-1
    !  do j=0,nz1
    !    do i=0,nx1
    !      tmp(i,k+loc_y_sta,j) = local_data((k*ny+j)*nx+i)
    !      !print*, i, j, k, (k*ny+j)*nx+i
    !    end do
    !  end do
    !end do
    
    !do k=0,nz1
    !  do j=0,ny1
    !    do i=0,nx1
    !      tmp(i,j,k) = local_data((k*ny+j)*nx+i)
    !      !print*, i, j, k, (k*ny+j)*nx+i
    !    end do
    !  end do
    !end do

    !call MPI_ALLREDUCE(tmp, tmp_var, nx*ny*nz, &
    !                   MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

    !do k=ksta,kend
    !  do j=jsta,jend
    !    out_var(:,j,k) = tmp_var(:,j,k)
    !  end do
    !end do

    !open (unit_no, file=proc_dir//'fft.dat')
    !do j=0,ny1
    !  write (unit_no, '(3e17.9)') (y(j), z(k), abs(tmp_var(nx/2,j,k)), k=0,nz1)
    !  write (unit_no, *)
    !end do
    !close (unit_no)
    
    !if (myrank == 0) then
    !  open (55, file='fft.dat')
    !  do j=0,ny1
    !    write (55, '(3e17.9)') (y(j), z(k), abs(tmp_var(nx/2,j,k)), k=0,nz1)
    !    write (55, *)
    !  end do
    !end if
    
    call fftwnd_f77_mpi_destroy_plan(plan)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    call fftw3d_f77_mpi_create_plan(iplan, MPI_COMM_WORLD, nx, nz, ny, &
                                    FFTW_BACKWARD, FFTW_ESTIMATE)

    !call fftwnd_f77_mpi_local_sizes(iplan, loc_ny, loc_y_sta, &
    !                                       loc_nz, loc_z_sta, &
    !                                       tot_loc)
    !                                      
    !print*, loc_nz, loc_z_sta, loc_ny, loc_y_sta, tot_loc

    !allocate(local_data(0:tot_loc-1))
    !allocate(work(0:tot_loc-1))
    !
    call fftwnd_f77_mpi(iplan, 1, local_data, work, 1, FFTW_TRANSPOSED_ORDER)
    
    tmp = 0.0
    do k=0,loc_nz-1
      do j=0,ny1
        do i=0,nx1
          tmp(i,j,k+loc_z_sta) = local_data((k*ny+j)*nx+i)
          !print*, i, j, k, (k*ny+j)*nx+i
        end do
      end do
    end do
    
    call MPI_ALLREDUCE(tmp, tmp_var, nx*ny*nz, &
                       MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

    tmp_var = tmp_var / (nx*ny*nz)
    in_var = 0.0
    do k=ksta,kend
      do j=jsta,jend
        in_var(:,j,k) = tmp_var(:,j,k)
      end do
    end do

    open (unit_no, file=proc_dir//'fft.dat')
    do j=0,ny1
      write (unit_no, '(3e17.9)') (y(j), z(k), abs(tmp_var(nx/2,j,k)), k=0,nz1)
      write (unit_no, *)
    end do
    close (unit_no)

    deallocate(local_data)
    deallocate(work)

    call fftwnd_f77_mpi_destroy_plan(iplan)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    return
  end subroutine fft
  
end module ic
