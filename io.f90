! $Id: io.f90,v 1.36 2006-12-01 12:52:14 n8049290 Exp $
!----------------------------------------------------------------------------

module io
  ! Routines for input/output
  use parameters
  implicit none

  private
  public :: open_files, close_files, save_time, save_energy, &
            save_surface, idl_surface, end_state, get_zeros, get_re_im_zeros, &
            get_extra_zeros, save_linelength, save_momentum, &
            get_dirs, diag, condensed_particles, average
  
  contains

  function itos(n)
    ! Convert an integer into a string of length 7
    implicit none

    integer, intent(in) :: n 
    integer             :: i, n_, d(7)
    character(7)        :: itos
    character           :: c(0:9) = (/'0','1','2','3','4','5',&
                                      '6','7','8','9'/)

    n_ = n
    do i = 7, 1, -1
      d(i) = mod(n_,10)
      n_ = n_ / 10
    end do

    itos = c(d(1))//c(d(2))//c(d(3))//c(d(4))//c(d(5))//c(d(6))//c(d(7))

    return
  end function itos

! ***************************************************************************  

  subroutine open_files()
    ! Open runtime files
    implicit none

    if (myrank == 0) then
      open (10, file='u_time.dat', status='unknown')
      open (11, file='timestep.dat', status='unknown')
      open (12, file='energy.dat', status='unknown')
      open (13, file='linelength.dat', status='unknown')
      open (14, file='momentum.dat', status='unknown')
      open (15, file='mass.dat', status='unknown')
      open (17, file='eta_time.dat', status='unknown')
      open (18, file='filtered_ll.dat', status='unknown')
      open (97, file='misc.dat', status='unknown')
      open (99, file='RUNNING')
      close (99)
    end if

    return
  end subroutine open_files

! ***************************************************************************  

  subroutine close_files()
    ! Close runtime files
    implicit none

    if (myrank == 0) then
      close (10)
      close (11)
      close (12)
      close (13)
      close (14)
      close (15)
      close (17)
      close (18)
      close (97)
    end if

    return
  end subroutine close_files

! ***************************************************************************  

  subroutine get_dirs()
    ! Get named directories where each process can write its own data
    use parameters
    implicit none

    ! Directory names are proc** with ** replaced by the rank of each process
    if (myrank < 10) then
      write (proc_dir(6:6), '(1i1)') myrank
      write (proc_dir(5:5), '(1i1)') 0
    else
      write (proc_dir(5:6), '(1i2)') myrank
    end if

    return
  end subroutine get_dirs

! ***************************************************************************  

  subroutine save_time(time, in_var)
    ! Save time-series data
    use parameters
    use variables, only : get_phase, get_density
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: phase, density
    complex, dimension(3) :: tmp, var
    real :: xpos, ypos, zpos
    integer :: i, j, k

    call get_phase(in_var, phase)
    call get_density(in_var, density)

    xpos = nx/2
    ypos = ny/2
    zpos = nz/2

    tmp = 0.0
    ! Find out on which process the data occurs and copy it into a temporary
    ! array
    do k=ksta,kend
      do j=jsta,jend
        do i=0,nx1
          if ((i==xpos) .and. (j==ypos) .and. (k==zpos)) then
            tmp(1) = in_var(xpos,ypos,zpos)
            tmp(2) = density(xpos,ypos,zpos)
            tmp(3) = phase(xpos,ypos,zpos)
          end if
        end do
      end do
    end do

    ! Make sure process 0 has the correct data to write
    call MPI_REDUCE(tmp, var, 3, MPI_COMPLEX, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)

    ! Write the data to file
    if (myrank == 0) then
      write (10, '(6e17.9)') time, im_t, real(var(1)), &
                             aimag(var(1)), real(var(2)), real(var(3))
    end if

    return
  end subroutine save_time

! ***************************************************************************  

  subroutine save_energy(time, in_var)
    ! Save the energy
    use parameters
    use variables, only : energy
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx1,jsta-2:jsta+2,ksta-2:kend+2), intent(in) :: in_var
    real :: E

    call energy(in_var, E)
    
    if (myrank == 0) then
      write (12, '(2e17.9)') time, E/real(nx*ny*nz)
    end if

    return
  end subroutine save_energy
  
! ***************************************************************************  

  subroutine condensed_particles(time, in_var)
    ! Calculate the mass, calculate the total energy and temperature, save
    ! spectra and save filtered isosurface
    use parameters
    use ic, only : fft, x, y, z
    use variables, only : mass, unit_no, send_recv_z, send_recv_y, &
                          pack_y, unpack_y
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: a, filtered
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2) :: a_tmp
    real :: M, n0, temp, temp2, tot, tmp, &
            rho0, E0, H, k2, k4, kx, ky, kz, kc, dk
    integer :: i, j, k, ii, jj, kk, ii2, jj2, kk2, V

    kc = pi !sqrt(kc2)
    dk = 2.0*kc/real(nx1)
    V = nx*ny*nz
    
    call fft(in_var, a, 'backward', .true.)

    a_tmp = 0.0
    a_tmp(:,jsta:jend,ksta:kend) = a
    call send_recv_z(a_tmp)
    call pack_y(a_tmp)
    call send_recv_y()
    call unpack_y(a_tmp)
    
    !call fft(a, in_var, 'forward', .true.)
    !open (unit_no, status='unknown', file=proc_dir//'fft'//itos(p)//'.dat', &
    !      form='unformatted')
    !
    !write (unit_no) nx, ny, nz
    !write (unit_no) nyprocs, nzprocs
    !write (unit_no) jsta, jend, ksta, kend
    !write (unit_no) abs(in_var)
    !write (unit_no) x
    !write (unit_no) y
    !write (unit_no) z

    !close (unit_no)

    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !stop
    
    ! Calculate the number of condensed particles
    do k=ksta,kend
      do j=jsta,jend
        if ((j==0) .and. (k==0)) then
          n0 = abs(a(0,0,0))**2
          exit
        end if
      end do
    end do

    ! Calculate the mass
    call mass(in_var, M)
    
    call MPI_BCAST(n0, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

    tmp = 0.0

    ! Density of condensed particles
    rho0 = n0/V
    
    ! Total energy <H>, and temperature T (?)
    tmp = 0.0
    do k=ksta,kend
      kz = -kc+real(k)*dk
      do j=jsta,jend
        ky = -kc+real(j)*dk
        do i=0,nx1
          kx = -kc+real(i)*dk
          k2 = -(1.0/12.0) * &  ! minus sign by comparison with spectral
                ( ((-2.0*cos(2.0*kx*dx) + 32.0*cos(kx*dx) - 30.0) / dx2) + &
                  ((-2.0*cos(2.0*ky*dy) + 32.0*cos(ky*dy) - 30.0) / dy2) + &
                  ((-2.0*cos(2.0*kz*dz) + 32.0*cos(kz*dz) - 30.0) / dz2) )
          if (k2==0.0) cycle
          k4 = k2**2
          tmp = tmp + (k2+rho0)/(k4+2.0*rho0*k2)
        end do
      end do
    end do
    
    call MPI_REDUCE(tmp, tot, 1, MPI_REAL, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
      temp = (M-n0)/tot
      E0 = (1.0/real(2*V))*(M**2+(M-n0)**2)
      H = E0 + temp*real(V-1)
      temp2 = ((M/(8.0*xr*yr*zr))-(n0/V))/tot
      write (15, '(8e17.9)') time, M/(8.0*xr*yr*zr), n0/V, M, n0, &
                             temp, temp2, H/V
    end if

    if (save_spectrum) then
      ! Save the spectrum
      call spectrum(a)
    end if
    
    if (save_filter) then
      ! Save a filtered isosurface
      call filtered_surface(a, 0)
    end if

    ! Save the linelength of a filtered isosurface
    call save_linelength(t, a_tmp, 1)
    
    return
  end subroutine condensed_particles
  
! ***************************************************************************  

  subroutine save_momentum(time, in_var)
    ! Save the momentum
    use parameters
    use variables, only : momentum
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    real, dimension(3) :: P

    call momentum(in_var, P)
    
    if (myrank == 0) then
      write (14, '(4e17.9)') time, P(1), P(2), P(3)
    end if

    return
  end subroutine save_momentum

! ***************************************************************************  

  subroutine save_surface(p, in_var)
    ! Save 2D surface data for use in gnuplot.  The data is saved separately on
    ! each process so a shell script must be used to plot it
    use parameters
    use variables
    use ic, only : x, y, z
    implicit none

    integer, intent(in) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: phase, density
    real :: zpos
    integer :: i, j, k

    ! Get the phase and the density
    call get_phase(in_var, phase)
    call get_density(in_var, density)
    
    zpos = nz/2

    ! Write each process's own data to file, but only if 'zpos' resides on that
    ! particular process
    do k=ksta,kend
      if (k==zpos) then
        open (unit_no, status='unknown', file=proc_dir//'u'//itos(p)//'.dat')
        do i=0,nx1
          write (unit_no, '(6e17.9)') (x(i), y(j), density(i,j,zpos), &
                                  phase(i,j,zpos), real(in_var(i,j,zpos)), &
                                  aimag(in_var(i,j,zpos)), j=jsta,jend)
          write (unit_no, *)
        end do
        close (unit_no)
        exit
      end if
    end do
    
    return
  end subroutine save_surface
  
! ***************************************************************************  

  subroutine idl_surface(p, in_var)
    ! Save 3D isosurface data for use in IDL.  As for the gnuplot plots, this
    ! data is saved separately for each process.  It must be read in through
    ! IDL
    use parameters
    use variables, only : unit_no
    use ic, only : x, y, z
    implicit none

    integer, intent(in) :: p
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: density
    integer :: i, j, k

    density = abs(in_var)**2

    call get_minmax(density, 'dens')

    open (unit_no, status='unknown', file=proc_dir//'dens'//itos(p)//'.dat', &
          form='unformatted')
    
    write (unit_no) nx, ny, nz
    write (unit_no) nyprocs, nzprocs
    write (unit_no) jsta, jend, ksta, kend
    write (unit_no) density
    write (unit_no) x
    write (unit_no) y
    write (unit_no) z

    close (unit_no)

    return
  end subroutine idl_surface

! ***************************************************************************  
  
  subroutine filtered_surface(a, flag)
    ! Save a filtered 3D isosurface.  High-frequency harmonics are filtered
    use parameters
    use ic, only : fft, x, y, z
    use variables, only : unit_no
    implicit none

    integer, intent(in) :: flag
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(inout) :: a
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: filtered
    integer :: i, j, k, ii2, jj2, kk2, k2

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
          k2 = ii2 + jj2 + kk2
          a(i,j,k) = a(i,j,k)*max(1.0-(real(k2)/kc2),0.0)
        end do
      end do
    end do
      
    call fft(a, filtered, 'forward', .true.)
    
    select case (flag)
      case (0)
        open (unit_no, status='unknown', &
                      file=proc_dir//'filtered'//itos(p)//'.dat', &
                      form='unformatted')
        write (unit_no) nx, ny, nz
        write (unit_no) nyprocs, nzprocs
        write (unit_no) jsta, jend, ksta, kend
        write (unit_no) abs(filtered)**2
        write (unit_no) x
        write (unit_no) y
        write (unit_no) z
      case (1)
        open (unit_no, status='unknown', &
        file=proc_dir//'end_state_filtered.dat', &
                      form='unformatted')
        write (unit_no) nx
        write (unit_no) ny
        write (unit_no) nz
        write (unit_no) p
        write (unit_no) t
        write (unit_no) dt
        write (unit_no) filtered
      case default
        STOP 'ERROR:  Invalid flag (filtered_surface)'
    end select
    
    close (unit_no)

    call get_minmax(abs(filtered)**2, 'filtered')

    return
  end subroutine filtered_surface

! ***************************************************************************  

  subroutine spectrum(a)
    ! Calculate and save the spectrum
    use parameters
    use ic, only : fft, x, y, z
    use variables, only : unit_no
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: a
    real :: log2k
    integer :: i, j, k, m, ii2, jj2, kk2, k2
    integer, parameter :: nshells = 7
    real, dimension(nshells) :: eta, tot_eta
    integer, dimension(nshells) :: nharm, tot_nharm

    eta = 0.0
    nharm = 0
    
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
          k2 = ii2 + jj2 + kk2
          if (sqrt(real(k2)) == 0.0) cycle
          log2k = log( 0.5*sqrt(real(k2))/pi ) / log(2.0)
          !log2k = log( sqrt(real(k2)) ) / log(2.0)
          !print*, i, j, k, sqrt(real(k2)), log2k
          do m=1,nshells
            if ((abs(log2k) < real(m)) .and. (abs(log2k) >= real(m-1))) then
              eta(m) = eta(m) + abs(a(i,j,k))**2
              nharm(m) = nharm(m) + 1
              exit
            end if
          end do
          !open (unit_no, position='append', &
          !               file=proc_dir//'spectrum'//itos(p)//'.dat')
          !write (unit_no, '(2e17.9)') sqrt(real(k2)), abs(a(i,j,k))**2
          !close (unit_no)
        end do
      end do
    end do

    call MPI_REDUCE(eta, tot_eta, nshells, MPI_REAL, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(nharm, tot_nharm, nshells, MPI_INTEGER, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr)

    if (myrank == 0) then
      do m=1,nshells
        if (tot_eta(m) == 0.0) cycle
        tot_eta(m) = tot_eta(m)/real(tot_nharm(m))
      end do
    
      !open (49, file='comb_spect.dat')
      open (unit_no, file=proc_dir//'spectrum'//itos(p)//'.dat')
      do m=1,nshells
        write (unit_no, '(2i9,e17.9)') m, tot_nharm(m), tot_eta(m)
      end do
      close (unit_no)

      write (17, '(4e17.9)') t, tot_eta(1), tot_eta(2), tot_eta(3)
    end if

    !open (unit_no, file=proc_dir//'spectrum'//itos(p)//'.dat')
    !
    !do m=1,nshells
    !  write (unit_no, '(2i9,e17.9)') m, nharm(m), eta(m)
    !end do
    !
    !close (unit_no)

    return
  end subroutine spectrum
  
! ***************************************************************************  

  subroutine get_minmax(in_var, var)
    ! Find the minimum and maximum values of a variable over time, and save the
    ! overall maximum to file
    use parameters
    implicit none

    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(2) :: maxs, maxr, mins, minr
    character(*), intent(in) :: var
    integer :: k

    ! maxs/r and mins/r are arrays of length 2 because the MPI functions find
    ! the max/min value as well as its location
    ! Find max/min on each process
    do k=ksta,kend
      if (k /= nz/2) cycle
      maxs(1) = maxval(in_var(:,:,nz/2))
      maxs(2) = 0.0
      mins(1) = minval(in_var(:,:,nz/2))
      mins(2) = 0.0
    end do
    
    ! Find max/min over whole array
    call MPI_ALLREDUCE(maxs, maxr, 1, MPI_2REAL, MPI_MAXLOC, &
                       MPI_COMM_WORLD, ierr)
                       
    call MPI_ALLREDUCE(mins, minr, 1, MPI_2REAL, MPI_MINLOC, &
                       MPI_COMM_WORLD, ierr)

    ! Update max/min if correct conditions are met
    select case (var)
      case ('dens')
        if (minr(1) < minvar(1)) then
          minvar(1) = minr(1)
        end if

        if (maxr(1) > maxvar(1)) then
          maxvar(1) = maxr(1)
        end if
        
        ! Save current max/min to file
        if (myrank == 0) then
          !print*, 'dens', minvar(1), maxvar(1)
          open (16, file='minmax_'//var//'.dat', form='unformatted')
          write (16) minvar(1)
          write (16) maxvar(1)
          close (16)
        end if
        
      case ('ave')
        if (minr(1) < minvar(2)) then
          minvar(2) = minr(1)
        end if

        if (maxr(1) > maxvar(2)) then
          maxvar(2) = maxr(1)
        end if
        
        ! Save current max/min to file
        if (myrank == 0) then
          !print*, 'ave', minvar(2), maxvar(2)
          open (16, file='minmax_'//var//'.dat', form='unformatted')
          write (16) minvar(2)
          write (16) maxvar(2)
          close (16)
        end if
        
      case ('filtered')
        if (minr(1) < minvar(3)) then
          minvar(3) = minr(1)
        end if

        if (maxr(1) > maxvar(3)) then
          maxvar(3) = maxr(1)
        end if
        
        ! Save current max/min to file
        if (myrank == 0) then
          !print*, 'filtered', minvar(3), maxvar(3)
          open (16, file='minmax_'//var//'.dat', form='unformatted')
          write (16) minvar(3)
          write (16) maxvar(3)
          close (16)
        end if
      case default
        STOP 'ERROR: Unrecognised variable (get_minmax)'
    end select
    
    return
  end subroutine get_minmax
    
! ***************************************************************************  

  subroutine end_state(in_var, p, flag)
    ! Save variables for use in a restarted run.  Each process saves its own
    ! bit
    use parameters
    use ic, only : fft
    use variables, only : unit_no
    implicit none

    integer, intent(in) :: p, flag
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: a
    integer :: j, k

    open (unit_no, file=proc_dir//'end_state.dat', form='unformatted')

    write (unit_no) nx
    write (unit_no) ny
    write (unit_no) nz
    write (unit_no) p
    write (unit_no) t
    write (unit_no) dt
    write (unit_no) in_var

    close (unit_no)

    ! Write the variables at the last save
    if (myrank == 0) then
      open (98, file = 'save.dat')
      write (98, *) 'Periodically saved state'
      write (98, *) 't=', t
      write (98, *) 'dt=', dt
      write (98, *) 'nx=', nx
      write (98, *) 'ny=', ny
      write (98, *) 'nz=', nz
      write (98, *) 'p=', p
      close (98)
    end if
    
    ! flag = 1 if the run has been ended
    if (flag == 1) then
      ! Save a final filtered isosurface
      call fft(in_var, a, 'backward', .true.)
      call filtered_surface(a, flag)
      if (myrank == 0) then
        ! Delete RUNNING file to cleanly terminate the run
        open (99, file = 'RUNNING')
        close (99, status = 'delete')
      end if
    end if
    
    !if (myrank == 0) then
    !  ! flag = 1 if the run has been ended
    !  if (flag == 1) then
    !    ! Delete RUNNING file to cleanly terminate the run
    !    open (99, file = 'RUNNING')
    !    close (99, status = 'delete')
    !  end if
    !end if

    return
  end subroutine end_state
  
! ***************************************************************************  

  subroutine get_zeros(in_var, p)
    ! Find all the zeros of the wavefunction by determining where the real and
    ! imaginary parts simultaneously go to zero.  This routine doesn't find
    ! them all though - get_extra_zeros below finds the rest
    
    ! All the zeros routines are horrible - I'm sure there is a more efficient
    ! way of calculating them
    use parameters
    use variables, only : re_im, unit_no
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real :: zero
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: denom
    integer :: i, j, k
    !integer, parameter :: z_start=nz/2, z_end=nz/2
    integer :: z_start, z_end

    ! Decide whether to find the zeros over the whole 3D box or just over a 2D
    ! plane
    z_start=ksta
    z_end=kend

    allocate(var%re(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
    allocate(var%im(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
    
    open (unit_no, status='unknown', file=proc_dir//'zeros'//itos(p)//'.dat')
    
    var%re = real(in_var)
    var%im = aimag(in_var)

    write (unit_no, *) "# i,j,k --> i+1,j,k"

    !do k=nz/2+0, nz/2+0
    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if (((var%re(i,j,k) == 0.0) .and. &
               (var%im(i,j,k) == 0.0)) .or. &
              ((var%re(i,j,k) == 0.0) .and. &
               (var%im(i,j,k)*var%im(i+1,j,k) < 0.0)) .or. &
              ((var%im(i,j,k) == 0.0) .and. &
               (var%re(i,j,k)*var%re(i+1,j,k) < 0.0)) .or. &
              ((var%re(i,j,k)*var%re(i+1,j,k) < 0.0) .and. &
               (var%im(i,j,k)*var%im(i+1,j,k) < 0.0)) ) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i+1,j,k)*x(i)/denom(i,j,k) + &
                    var%re(i,j,k)*x(i+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (unit_no, *) "# i,j,k --> i,j+1,k"

    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if (((var%re(i,j,k) == 0.0) .and. &
               (var%im(i,j,k) == 0.0)) .or. &
              ((var%re(i,j,k) == 0.0) .and. &
               (var%im(i,j,k)*var%im(i,j+1,k) < 0.0)) .or. &
              ((var%im(i,j,k) == 0.0) .and. &
               (var%re(i,j,k)*var%re(i,j+1,k) < 0.0)) .or. &
              ((var%re(i,j,k)*var%re(i,j+1,k) < 0.0) .and. &
               (var%im(i,j,k)*var%im(i,j+1,k) < 0.0)) ) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i,j+1,k)*y(j)/denom(i,j,k) + &
                    var%re(i,j,k)*y(j+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (unit_no, *) "# i,j,k --> i,j,k+1"

      do k=z_start,z_end
        if ((k==0) .or. (k==nz1)) cycle
        do j=jsta,jend
          if ((j==0) .or. (j==ny1)) cycle
          do i=1,nx1-1
            if (((var%re(i,j,k) == 0.0) .and. &
                 (var%im(i,j,k) == 0.0)) .or. &
                ((var%re(i,j,k) == 0.0) .and. &
                 (var%im(i,j,k)*var%im(i,j,k+1) < 0.0)) .or. &
                ((var%im(i,j,k) == 0.0) .and. &
                 (var%re(i,j,k)*var%re(i,j,k+1) < 0.0)) .or. &
                ((var%re(i,j,k)*var%re(i,j,k+1) < 0.0) .and. &
                 (var%im(i,j,k)*var%im(i,j,k+1) < 0.0)) ) then
              denom(i,j,k) = var%re(i,j,k)-var%re(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -var%re(i,j,k+1)*z(k)/denom(i,j,k) + &
                      var%re(i,j,k)*z(k+1)/denom(i,j,k)
              write (unit_no, '(3e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if

    close (unit_no)

    deallocate(var%re)
    deallocate(var%im)

    return
  end subroutine get_zeros

! ***************************************************************************  

  subroutine get_extra_zeros(in_var, p)
    ! Find the zeros that the get_zeros routine did not pick up
    use parameters
    use variables, only : re_im, unit_no
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real, dimension(4) :: zero
    real, dimension(2) :: m
    real, dimension(4,0:nx1,jsta:jend,ksta:kend) :: denom
    real :: xp, yp, zp
    integer :: i, j, k
    !integer, parameter :: z_start=nz/2, z_end=nz/2
    integer :: z_start, z_end

    z_start=ksta
    z_end=kend

    allocate(var%re(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
    allocate(var%im(0:nx1,jsta-2:jend+2,ksta-2:kend+2))

    ! Write these new zeros to the same file as for the get_zeros routine
    open (unit_no, status='old', position='append', &
                   file=proc_dir//'zeros'//itos(p)//'.dat')

    var%re = real(in_var)
    var%im = aimag(in_var)
    
    write (unit_no, *) "# i,j,k --> i+1,j,k --> i+1,j+1,k --> i,j+1,k"
    
    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if (((var%re(i,j,k)*var%re(i+1,j,k) < 0.0) .and. &
               (var%im(i,j,k)*var%im(i+1,j,k) >= 0.0)) .and. &
              ((var%im(i+1,j,k)*var%im(i+1,j+1,k) < 0.0) .and. &
               (var%re(i+1,j,k)*var%re(i+1,j+1,k) >= 0.0)) .and. &
              ((var%re(i+1,j+1,k)*var%re(i,j+1,k) < 0.0) .and. &
               (var%im(i+1,j+1,k)*var%im(i,j+1,k) >= 0.0)) .and. &
              ((var%im(i,j+1,k)*var%im(i,j,k) < 0.0) .and. &
               (var%re(i,j+1,k)*var%re(i,j,k) >= 0.0)) ) then
            denom(1,i,j,k) = var%re(i,j,k)-var%re(i+1,j,k)
            zero(1) = -var%re(i+1,j,k)*x(i)/denom(1,i,j,k) + &
                       var%re(i,j,k)*x(i+1)/denom(1,i,j,k)
            denom(2,i,j,k) = var%im(i+1,j,k)-var%im(i+1,j+1,k)
            zero(2) = -var%im(i+1,j+1,k)*y(j)/denom(2,i,j,k) + &
                       var%im(i+1,j,k)*y(j+1)/denom(2,i,j,k)
            denom(3,i,j,k) = var%re(i+1,j+1,k)-var%re(i,j+1,k)
            zero(3) = -var%re(i,j+1,k)*x(i+1)/denom(3,i,j,k) + &
                       var%re(i+1,j+1,k)*x(i)/denom(3,i,j,k)
            denom(4,i,j,k) = var%im(i,j+1,k)-var%im(i,j,k)
            zero(4) = -var%im(i,j,k)*y(j+1)/denom(4,i,j,k) + &
                       var%im(i,j+1,k)*y(j)/denom(4,i,j,k)
            m(1) = (y(j)-y(j+1))/(zero(1)-zero(3))
            m(2) = (zero(2)-zero(4))/(x(i+1)-x(i))
            xp = (zero(4)-x(i)*m(2)-y(j)+zero(1)*m(1))/(m(1)-m(2))
            yp = xp*m(1)+y(j)-zero(1)*m(1)
            
            write (unit_no, '(3e17.9)') xp, yp, z(k)
          else if (((var%im(i,j,k)*var%im(i+1,j,k) < 0.0) .and. &
                    (var%re(i,j,k)*var%re(i+1,j,k) >= 0.0)) .and. &
                   ((var%re(i+1,j,k)*var%re(i+1,j+1,k) < 0.0) .and. &
                    (var%im(i+1,j,k)*var%im(i+1,j+1,k) >= 0.0)) .and. &
                   ((var%im(i+1,j+1,k)*var%im(i,j+1,k) < 0.0) .and. &
                    (var%re(i+1,j+1,k)*var%re(i,j+1,k) >= 0.0)) .and. &
                   ((var%re(i,j+1,k)*var%re(i,j,k) < 0.0) .and. &
                    (var%im(i,j+1,k)*var%im(i,j,k) >= 0.0)) ) then
            denom(1,i,j,k) = var%im(i,j,k)-var%im(i+1,j,k)
            zero(1) = -var%im(i+1,j,k)*x(i)/denom(1,i,j,k) + &
                       var%im(i,j,k)*x(i+1)/denom(1,i,j,k)
            denom(2,i,j,k) = var%re(i+1,j,k)-var%re(i+1,j+1,k)
            zero(2) = -var%re(i+1,j+1,k)*y(j)/denom(2,i,j,k) + &
                       var%re(i+1,j,k)*y(j+1)/denom(2,i,j,k)
            denom(3,i,j,k) = var%im(i+1,j+1,k)-var%im(i,j+1,k)
            zero(3) = -var%im(i,j+1,k)*x(i+1)/denom(3,i,j,k) + &
                       var%im(i+1,j+1,k)*x(i)/denom(3,i,j,k)
            denom(4,i,j,k) = var%re(i,j+1,k)-var%re(i,j,k)
            zero(4) = -var%re(i,j,k)*y(j+1)/denom(4,i,j,k) + &
                       var%re(i,j+1,k)*y(j)/denom(4,i,j,k)
            m(1) = (y(j)-y(j+1))/(zero(1)-zero(3))
            m(2) = (zero(2)-zero(4))/(x(i+1)-x(i))
            xp = (zero(4)-x(i)*m(2)-y(j)+zero(1)*m(1))/(m(1)-m(2))
            yp = xp*m(1)+y(j)-zero(1)*m(1)

            write (unit_no, '(3e17.9)') xp, yp, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (unit_no, *) "# i,j,k --> i+1,j,k --> i+1,j,k+1 --> i,j,k+1"
      do k=z_start,z_end
        if ((k==0) .or. (k==nz1)) cycle
        do j=jsta,jend
          if ((j==0) .or. (j==ny1)) cycle
          do i=1,nx1-1
            if (((var%re(i,j,k)*var%re(i+1,j,k) < 0.0) .and. &
                 (var%im(i,j,k)*var%im(i+1,j,k) >= 0.0)) .and. &
                ((var%im(i+1,j,k)*var%im(i+1,j,k+1) < 0.0) .and. &
                 (var%re(i+1,j,k)*var%re(i+1,j,k+1) >= 0.0)) .and. &
                ((var%re(i+1,j,k+1)*var%re(i,j,k+1) < 0.0) .and. &
                 (var%im(i+1,j,k+1)*var%im(i,j,k+1) >= 0.0)) .and. &
                ((var%im(i,j,k+1)*var%im(i,j,k) < 0.0) .and. &
                 (var%re(i,j,k+1)*var%re(i,j,k) >= 0.0)) ) then
              denom(1,i,j,k) = var%re(i,j,k)-var%re(i+1,j,k)
              zero(1) = -var%re(i+1,j,k)*x(i)/denom(1,i,j,k) + &
                         var%re(i,j,k)*x(i+1)/denom(1,i,j,k)
              denom(2,i,j,k) = var%im(i+1,j,k)-var%im(i+1,j,k+1)
              zero(2) = -var%im(i+1,j,k+1)*z(k)/denom(2,i,j,k) + &
                         var%im(i+1,j,k)*z(k+1)/denom(2,i,j,k)
              denom(3,i,j,k) = var%re(i+1,j,k+1)-var%re(i,j,k+1)
              zero(3) = -var%re(i,j,k+1)*x(i+1)/denom(3,i,j,k) + &
                         var%re(i+1,j,k+1)*x(i)/denom(3,i,j,k)
              denom(4,i,j,k) = var%im(i,j,k+1)-var%im(i,j,k)
              zero(4) = -var%im(i,j,k)*z(k+1)/denom(4,i,j,k) + &
                         var%im(i,j,k+1)*z(k)/denom(4,i,j,k)
              m(1) = (z(k)-z(k+1))/(zero(1)-zero(3))
              m(2) = (zero(2)-zero(4))/(x(i+1)-x(i))
              xp = (zero(4)-x(i)*m(2)-z(k)+zero(1)*m(1))/(m(1)-m(2))
              zp = xp*m(1)+z(k)-zero(1)*m(1)
              
              write (unit_no, '(3e17.9)') xp, y(j), zp
            else if (((var%im(i,j,k)*var%im(i+1,j,k) < 0.0) .and. &
                      (var%re(i,j,k)*var%re(i+1,j,k) >= 0.0)) .and. &
                     ((var%re(i+1,j,k)*var%re(i+1,j,k+1) < 0.0) .and. &
                      (var%im(i+1,j,k)*var%im(i+1,j,k+1) >= 0.0)) .and. &
                     ((var%im(i+1,j,k+1)*var%im(i,j,k+1) < 0.0) .and. &
                      (var%re(i+1,j,k+1)*var%re(i,j,k+1) >= 0.0)) .and. &
                     ((var%re(i,j,k+1)*var%re(i,j,k) < 0.0) .and. &
                      (var%im(i,j,k+1)*var%im(i,j,k) >= 0.0)) ) then
              denom(1,i,j,k) = var%im(i,j,k)-var%im(i+1,j,k)
              zero(1) = -var%im(i+1,j,k)*x(i)/denom(1,i,j,k) + &
                         var%im(i,j,k)*x(i+1)/denom(1,i,j,k)
              denom(2,i,j,k) = var%re(i+1,j,k)-var%re(i+1,j,k+1)
              zero(2) = -var%re(i+1,j,k+1)*z(k)/denom(2,i,j,k) + &
                         var%re(i+1,j,k)*z(k+1)/denom(2,i,j,k)
              denom(3,i,j,k) = var%im(i+1,j,k+1)-var%im(i,j,k+1)
              zero(3) = -var%im(i,j,k+1)*x(i+1)/denom(3,i,j,k) + &
                         var%im(i+1,j,k+1)*x(i)/denom(3,i,j,k)
              denom(4,i,j,k) = var%re(i,j,k+1)-var%re(i,j,k)
              zero(4) = -var%re(i,j,k)*z(k+1)/denom(4,i,j,k) + &
                         var%re(i,j,k+1)*z(k)/denom(4,i,j,k)
              m(1) = (z(k)-z(k+1))/(zero(1)-zero(3))
              m(2) = (zero(2)-zero(4))/(x(i+1)-x(i))
              xp = (zero(4)-x(i)*m(2)-z(k)+zero(1)*m(1))/(m(1)-m(2))
              zp = xp*m(1)+z(k)-zero(1)*m(1)

              write (unit_no, '(3e17.9)') xp, y(j), zp
            end if
          end do
        end do
      end do
      
      write (unit_no, *) "# i,j,k --> i,j,k+1 --> i,j+1,k+1 --> i,j+1,k"
      do k=z_start,z_end
        if ((k==0) .or. (k==nz1)) cycle
        do j=jsta,jend
          if ((j==0) .or. (j==ny1)) cycle
          do i=1,nx1-1
            if (((var%re(i,j,k)*var%re(i,j,k+1) < 0.0) .and. &
                 (var%im(i,j,k)*var%im(i,j,k+1) >= 0.0)) .and. &
                ((var%im(i,j,k+1)*var%im(i,j+1,k+1) < 0.0) .and. &
                 (var%re(i,j,k+1)*var%re(i,j+1,k+1) >= 0.0)) .and. &
                ((var%re(i,j+1,k+1)*var%re(i,j+1,k) < 0.0) .and. &
                 (var%im(i,j+1,k+1)*var%im(i,j+1,k) >= 0.0)) .and. &
                ((var%im(i,j+1,k)*var%im(i,j,k) < 0.0) .and. &
                 (var%re(i,j+1,k)*var%re(i,j,k) >= 0.0)) ) then
              denom(1,i,j,k) = var%re(i,j,k)-var%re(i,j,k+1)
              zero(1) = -var%re(i,j,k+1)*z(k)/denom(1,i,j,k) + &
                         var%re(i,j,k)*z(k+1)/denom(1,i,j,k)
              denom(2,i,j,k) = var%im(i,j,k+1)-var%im(i,j+1,k+1)
              zero(2) = -var%im(i,j+1,k+1)*y(j)/denom(2,i,j,k) + &
                         var%im(i,j,k+1)*y(j+1)/denom(2,i,j,k)
              denom(3,i,j,k) = var%re(i,j+1,k+1)-var%re(i,j+1,k)
              zero(3) = -var%re(i,j+1,k)*z(k+1)/denom(3,i,j,k) + &
                         var%re(i,j+1,k+1)*z(k)/denom(3,i,j,k)
              denom(4,i,j,k) = var%im(i,j+1,k)-var%im(i,j,k)
              zero(4) = -var%im(i,j,k)*y(j+1)/denom(4,i,j,k) + &
                         var%im(i,j+1,k)*y(j)/denom(4,i,j,k)
              m(1) = (y(j)-y(j+1))/(zero(1)-zero(3))
              m(2) = (zero(2)-zero(4))/(z(k+1)-z(k))
              zp = (zero(4)-z(k)*m(2)-y(j)+zero(1)*m(1))/(m(1)-m(2))
              yp = zp*m(1)+y(j)-zero(1)*m(1)
              
              write (unit_no, '(3e17.9)') x(i), yp, zp
            else if (((var%im(i,j,k)*var%im(i,j,k+1) < 0.0) .and. &
                      (var%re(i,j,k)*var%re(i,j,k+1) >= 0.0)) .and. &
                     ((var%re(i,j,k+1)*var%re(i,j+1,k+1) < 0.0) .and. &
                      (var%im(i,j,k+1)*var%im(i,j+1,k+1) >= 0.0)) .and. &
                     ((var%im(i,j+1,k+1)*var%im(i,j+1,k) < 0.0) .and. &
                      (var%re(i,j+1,k+1)*var%re(i,j+1,k) >= 0.0)) .and. &
                     ((var%re(i,j+1,k)*var%re(i,j,k) < 0.0) .and. &
                      (var%im(i,j+1,k)*var%im(i,j,k) >= 0.0)) ) then
              denom(1,i,j,k) = var%im(i,j,k)-var%im(i,j,k+1)
              zero(1) = -var%im(i,j,k+1)*z(k)/denom(1,i,j,k) + &
                         var%im(i,j,k)*z(k+1)/denom(1,i,j,k)
              denom(2,i,j,k) = var%re(i,j,k+1)-var%re(i,j+1,k+1)
              zero(2) = -var%re(i,j+1,k+1)*y(j)/denom(2,i,j,k) + &
                         var%re(i,j,k+1)*y(j+1)/denom(2,i,j,k)
              denom(3,i,j,k) = var%im(i,j+1,k+1)-var%im(i,j+1,k)
              zero(3) = -var%im(i,j+1,k)*z(k+1)/denom(3,i,j,k) + &
                         var%im(i,j+1,k+1)*z(k)/denom(3,i,j,k)
              denom(4,i,j,k) = var%re(i,j+1,k)-var%re(i,j,k)
              zero(4) = -var%re(i,j,k)*y(j+1)/denom(4,i,j,k) + &
                         var%re(i,j+1,k)*y(j)/denom(4,i,j,k)
              m(1) = (y(j)-y(j+1))/(zero(1)-zero(3))
              m(2) = (zero(2)-zero(4))/(z(k+1)-z(k))
              zp = (zero(4)-z(k)*m(2)-y(j)+zero(1)*m(1))/(m(1)-m(2))
              yp = zp*m(1)+y(j)-zero(1)*m(1)
              
              write (unit_no, '(3e17.9)') x(i), yp, zp
            end if
          end do
        end do
      end do
    end if
  
    close (unit_no)

    return
  end subroutine get_extra_zeros

! ***************************************************************************  

  subroutine get_re_im_zeros(in_var, p)
    ! Find where the real and imaginary parts separately go to zero
    use parameters
    use variables, only : re_im, unit_no
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real :: zero
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: denom
    integer :: i, j, k
    !integer, parameter :: z_start=nz/2, z_end=nz/2
    integer :: z_start, z_end

    z_start=ksta
    z_end=kend

    allocate(var%re(0:nx1,jsta-2:jend+2,ksta-2:kend+2))
    allocate(var%im(0:nx1,jsta-2:jend+2,ksta-2:kend+2))

    open (unit_no, status='unknown', &
                   file=proc_dir//'re_zeros'//itos(p)//'.dat')

    var%re = real(in_var)
    var%im = aimag(in_var)

    write (unit_no, *) "# i,j,k --> i+1,j,k"

    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if ((var%re(i,j,k) == 0.0) .or. &
              (var%re(i,j,k)*var%re(i+1,j,k) < 0.0)) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i+1,j,k)*x(i)/denom(i,j,k) + &
                    var%re(i,j,k)*x(i+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (unit_no, *) "# i,j,k --> i,j+1,k"

    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if ((var%re(i,j,k) == 0.0) .or. &
              (var%re(i,j,k)*var%re(i,j+1,k) < 0.0)) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i,j+1,k)*y(j)/denom(i,j,k) + &
                    var%re(i,j,k)*y(j+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (unit_no, *) "# i,j,k --> i,j,k+1"

      do k=z_start,z_end
        if ((k==0) .or. (k==nz1)) cycle
        do j=jsta,jend
          if ((j==0) .or. (j==ny1)) cycle
          do i=1,nx1-1
            if ((var%re(i,j,k) == 0.0) .and. &
                (var%re(i,j,k)*var%re(i,j,k+1) < 0.0)) then
              denom(i,j,k) = var%re(i,j,k)-var%re(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -var%re(i,j,k+1)*z(k)/denom(i,j,k) + &
                      var%re(i,j,k)*z(k+1)/denom(i,j,k)
              write (unit_no, '(3e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if

    close (unit_no)

    ! Barrier here to make sure some processes don't try to open the new file
    ! without it having been previously closed
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    open (unit_no, status='unknown', &
                   file=proc_dir//'im_zeros'//itos(p)//'.dat')
                      
    write (unit_no, *) "# i,j,k --> i+1,j,k"
    
    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if ((var%im(i,j,k) == 0.0) .or. &
              (var%im(i,j,k)*var%im(i+1,j,k) < 0.0)) then
            denom(i,j,k) = var%im(i,j,k)-var%im(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%im(i+1,j,k)*x(i)/denom(i,j,k) + &
                    var%im(i,j,k)*x(i+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (unit_no, *) "# i,j,k --> i,j+1,k"

    do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        do i=1,nx1-1
          if ((var%im(i,j,k) == 0.0) .or. &
              (var%im(i,j,k)*var%im(i,j+1,k) < 0.0)) then
            denom(i,j,k) = var%im(i,j,k)-var%im(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%im(i,j+1,k)*y(j)/denom(i,j,k) + &
                    var%im(i,j,k)*y(j+1)/denom(i,j,k)
            write (unit_no, '(3e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (unit_no, *) "# i,j,k --> i,j,k+1"

      do k=z_start,z_end
      if ((k==0) .or. (k==nz1)) cycle
        do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
          do i=1,nx1-1
            if ((var%im(i,j,k) == 0.0) .and. &
                (var%im(i,j,k)*var%im(i,j,k+1) < 0.0)) then
              denom(i,j,k) = var%im(i,j,k)-var%im(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -var%im(i,j,k+1)*z(k)/denom(i,j,k) + &
                      var%im(i,j,k)*z(k+1)/denom(i,j,k)
              write (unit_no, '(3e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if

    close (unit_no)

    return
  end subroutine get_re_im_zeros
  
! ***************************************************************************  

  subroutine save_linelength(t, in_var, flag)
    ! Save the total vortex line length
    use parameters
    use variables, only : linelength
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    real, intent(in) :: t
    integer :: flag
    real :: tmp, length

    ! Get the line length on each individual process
    tmp = linelength(t, in_var)

    ! Sum the line length over all processes and send it to process 0
    call MPI_REDUCE(tmp, length, 1, MPI_REAL, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ierr) 

    if (myrank == 0) then
      if (flag == 0) then
        ! Write the unfiltered line length
        write (13, '(2e17.9)') t, length
      else if (flag == 1) then
        ! Write the filtered line length
        write (18, '(2e17.9)') t, length
      end if
    end if

    return
  end subroutine save_linelength

! ***************************************************************************  

  subroutine diag(old2, old, new, p)
    use parameters
    use variables
    use ic, only : x, y, z
    implicit none
    
    integer, intent(in) :: p
    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: old, &
                                                                         old2
    complex, dimension(0:nx1,jsta:jend,ksta:kend),         intent(in) :: new
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: lhs, rhs
    real :: zpos
    integer :: i, j, k
    
    lhs = 0.5*(new-old2(:,jsta:jend,ksta:kend))/dt
    
    rhs = 0.5*(eye+0.01) * ( laplacian(old) + &
                           (1.0-abs(old(:,jsta:jend,ksta:kend))**2)*&
                                    old(:,jsta:jend,ksta:kend) )

    zpos = nz/2

    do k=ksta,kend
      if (k==zpos) then
        open (unit_no, status='unknown', file=proc_dir//'diag'//itos(p)//'.dat')
        do i=0,nx1
          write (unit_no, '(3e17.9)') (x(i), y(j), &
                 abs(rhs(i,j,zpos)-lhs(i,j,zpos)), j=jsta,jend)
          write (unit_no, *)
        end do
        close (unit_no)
        exit
      end if
    end do

    return
  end subroutine diag
    
! ***************************************************************************  

  subroutine average(in_var)
    ! Save time-averaged data
    use parameters
    use ic, only : x, y, z
    use variables, only : unit_no
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    integer :: j, k

    ave = ave + abs(in_var)**2
    
    open (unit_no, status='unknown', file=proc_dir//'ave'//itos(p)//'.dat', &
          form='unformatted')

    write (unit_no) nx, ny, nz
    write (unit_no) nyprocs, nzprocs
    write (unit_no) jsta, jend, ksta, kend
    write (unit_no) ave / real(snapshots)
    write (unit_no) x
    write (unit_no) y
    write (unit_no) z

    close (unit_no)

    call get_minmax(ave / real(snapshots), 'ave')

    snapshots = snapshots+1

    return
  end subroutine average
  
end module io
