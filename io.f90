module io
  use parameters
  implicit none

  private
  public :: open_files, close_files, save_time, save_energy, save_surface, &
            idl_surface, end_state, get_zeros, get_re_im_zeros, &
            get_phase_zeros, get_extra_zeros, save_linelength, save_momentum

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

  subroutine open_files()
    ! Open runtime files
    implicit none

    open (10, file='u_time.dat', status='unknown')
    open (13, file='timestep.dat', status='unknown')
    open (14, file='energy.dat', status='unknown')
    open (20, file='linelength.dat', status='unknown')
    open (21, file='momentum.dat', status='unknown')
    !open (90, file='diag.dat', status='unknown')
    open (99, file='RUNNING')
    close (99)

    return
  end subroutine open_files

  subroutine close_files()
    ! Close runtime files
    implicit none

    close (10)
    close (13)
    close (14)
    close (20)
    close (21)
    !close (90)

    return
  end subroutine close_files

  subroutine save_time(time, in_var)
    ! Save time-series data
    use parameters
    use variables, only : get_phase, get_density
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, dimension(0:nx,0:ny,0:nz) :: phase, density
    real :: xpos, ypos, zpos

    call get_phase(in_var, phase)
    call get_density(in_var, density)

    xpos = nx/2
    ypos = ny/2
    zpos = nz/2

    write (10, '(6e17.9)') time, im_t, real(in_var(xpos,ypos,zpos)), &
                           aimag(in_var(xpos,ypos,zpos)), &
                           density(xpos,ypos,zpos), &
                           phase(xpos,ypos,zpos)

    return
  end subroutine save_time

  subroutine save_energy(time, in_var)
    ! Save the energy
    use parameters
    use variables, only : energy
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real :: E

    call energy(in_var, E)
    
    write (14, '(2e17.9)') time, E

    return
  end subroutine save_energy
  
  subroutine save_momentum(time, in_var)
    ! Save the momentum
    use parameters
    use variables, only : momentum
    implicit none

    real, intent(in) :: time
    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, dimension(3) :: P

    call momentum(in_var, P)
    
    write (21, '(4e17.9)') time, P(1), P(2), P(3)

    return
  end subroutine save_momentum

  subroutine save_surface(p, in_var)
    ! Save 2D surface data for use in gnuplot
    use parameters
    use variables
    use ic, only : x, y, z
    implicit none

    integer, intent(in) :: p
    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, dimension(0:nx,0:ny,0:nz) :: phase, density
    real :: zpos
    integer :: i, j, k

    open (11, status = 'unknown', file = 'u'//itos(p)//'.dat')
    
    ! Get the phase and the density
    call get_phase(in_var, phase)
    call get_density(in_var, density)
    
    zpos = nz/2 !22 !44 !22
    
    do i=0,nx
      write (11, '(6e17.9)') (x(i), y(j), density(i,j,zpos), &
                              phase(i,j,zpos), real(in_var(i,j,zpos)), &
                              aimag(in_var(i,j,zpos)), j=0,ny)
      write (11, *)
    end do
    
    !do j=0,ny
    !  write (11, '(6e17.9)') (x(i), y(j), density(i,j,zpos), &
    !                          phase(i,j,zpos), real(in_var(i,j,zpos)), &
    !                          aimag(in_var(i,j,zpos)), i=0,nx)
    !  write (11, *)
    !end do
    
    !do i=0,nx
    !  write (11, '(4e17.9)') (x(i), z(k), density(i,ny/2,k), &
    !                          phase(i,ny/2,k), k=0,nz)
    !  write (11, *)
    !end do
    
    !do j=0,ny
    !  write (11, '(4e17.9)') (y(j), z(k), density(nx/2,j,k), &
    !                          phase(nx/2,j,k), k=0,nz)
    !  write (11, *)
    !end do

    close (11)

    return
  end subroutine save_surface
  
  subroutine idl_surface(p, in_var)
    ! Save 3D isosurface data for use in IDL
    use parameters
    use ic, only : x, y, z
    implicit none

    integer, intent(in) :: p
    complex, intent(in) :: in_var(0:nx,0:ny,0:nz)
    integer :: i, j, k

    open (12, status = 'unknown', file = 'u_idl'//itos(p)//'.dat', &
          form = 'unformatted')
    
    write (12) nx, ny, nz
    write (12) abs(in_var)
    write (12) x
    write (12) y
    write (12) z

    close (12)

    return
  end subroutine idl_surface

  subroutine end_state(in_var, p, flag)
    ! Save variables for use in a restarted run
    use parameters
    implicit none

    integer, intent(in) :: p, flag
    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    integer :: j, k

    open (50, file = 'end_state.dat', form='unformatted')

    write(50) nx
    write(50) ny
    write(50) nz
    write(50) p
    write(50) t
    write(50) dt
    write(50) in_var

    close (50)

    ! Write the variables at the last save
    open (98, file = 'save.dat')
    write (98, *) 'Periodically saved state'
    write (98, *) 't=', t
    write (98, *) 'dt=', dt
    write (98, *) 'nx=', nx
    write (98, *) 'ny=', ny
    write (98, *) 'nz=', nz
    write (98, *) 'p=', p
    close (98)
    
    if (flag == 1) then
      ! Delete RUNNING file to cleanly terminate the run
      open (99, file = 'RUNNING')
      close (99, status = 'delete')
    end if

    return
  end subroutine end_state
  
  subroutine get_zeros(in_var, p)
    use parameters
    use variables, only : re_im
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real :: zero
    real, dimension(0:nx,0:ny,0:nz) :: denom
    integer :: i, j, k
    !integer, parameter :: z_start=nz/2, z_end=nz/2
    integer, parameter :: z_start=1, z_end=nz1

    open (15, status = 'unknown', file = 'zeros'//itos(p)//'.dat')

    var%re = real(in_var)
    var%im = aimag(in_var)

    write (15, *) "# i,j,k --> i+1,j,k"

    !do k=nz/2+0, nz/2+0
    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
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
            write (15, '(5e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (15, *) "# i,j,k --> i,j+1,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
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
            write (15, '(5e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (15, *) "# i,j,k --> i,j,k+1"

      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
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
              write (15, '(5e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if

    close (15)

    return
  end subroutine get_zeros

  subroutine get_extra_zeros(in_var, p)
    use parameters
    use variables, only : re_im
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real, dimension(4) :: zero
    real, dimension(2) :: m
    real, dimension(4,0:nx,0:ny,0:nz) :: denom
    real :: xp, yp, zp
    integer :: i, j, k
    !integer, parameter :: z_start=nz/2, z_end=nz/2
    integer, parameter :: z_start=1, z_end=nz1

    open (15, status = 'old', position='append', &
              file = 'zeros'//itos(p)//'.dat')

    var%re = real(in_var)
    var%im = aimag(in_var)
    
    write (15, *) "# i,j,k --> i+1,j,k --> i+1,j+1,k --> i,j+1,k"
    
    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
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
            
            write (15, '(5e17.9)') xp, yp, z(k)
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

            write (15, '(5e17.9)') xp, yp, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (15, *) "# i,j,k --> i+1,j,k --> i+1,j,k+1 --> i,j,k+1"
      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
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
              
              write (15, '(5e17.9)') xp, y(j), zp
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

              write (15, '(5e17.9)') xp, y(j), zp
            end if
          end do
        end do
      end do
      
      write (15, *) "# i,j,k --> i,j,k+1 --> i,j+1,k+1 --> i,j+1,k"
      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
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
              
              write (15, '(5e17.9)') x(i), yp, zp
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
              
              write (15, '(5e17.9)') x(i), yp, zp
            end if
          end do
        end do
      end do
    end if
  
    close (15)

    return
  end subroutine get_extra_zeros

  subroutine get_re_im_zeros(in_var, p)
    use parameters
    use variables, only : re_im
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    integer, intent(in) :: p
    type (re_im) :: var
    real :: zero
    real, dimension(0:nx,0:ny,0:nz) :: denom
    integer :: i, j, k
    integer, parameter :: z_start=nz/2, z_end=nz/2
    !integer, parameter :: z_start=1, z_end=nz1

    open (16, status = 'unknown', file = 're_zeros'//itos(p)//'.dat')
    open (17, status = 'unknown', file = 'im_zeros'//itos(p)//'.dat')

    var%re = real(in_var)
    var%im = aimag(in_var)

    write (16, *) "# i,j,k --> i+1,j,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((var%re(i,j,k) == 0.0) .or. &
              (var%re(i,j,k)*var%re(i+1,j,k) < 0.0)) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i+1,j,k)*x(i)/denom(i,j,k) + &
                    var%re(i,j,k)*x(i+1)/denom(i,j,k)
            write (16, '(5e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (16, *) "# i,j,k --> i,j+1,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((var%re(i,j,k) == 0.0) .or. &
              (var%re(i,j,k)*var%re(i,j+1,k) < 0.0)) then
            denom(i,j,k) = var%re(i,j,k)-var%re(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%re(i,j+1,k)*y(j)/denom(i,j,k) + &
                    var%re(i,j,k)*y(j+1)/denom(i,j,k)
            write (16, '(5e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (16, *) "# i,j,k --> i,j,k+1"

      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
            if ((var%re(i,j,k) == 0.0) .and. &
                (var%re(i,j,k)*var%re(i,j,k+1) < 0.0)) then
              denom(i,j,k) = var%re(i,j,k)-var%re(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -var%re(i,j,k+1)*z(k)/denom(i,j,k) + &
                      var%re(i,j,k)*z(k+1)/denom(i,j,k)
              write (16, '(5e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if
    
    write (17, *) "# i,j,k --> i+1,j,k"
    
    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((var%im(i,j,k) == 0.0) .or. &
              (var%im(i,j,k)*var%im(i+1,j,k) < 0.0)) then
            denom(i,j,k) = var%im(i,j,k)-var%im(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%im(i+1,j,k)*x(i)/denom(i,j,k) + &
                    var%im(i,j,k)*x(i+1)/denom(i,j,k)
            write (17, '(5e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (17, *) "# i,j,k --> i,j+1,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((var%im(i,j,k) == 0.0) .or. &
              (var%im(i,j,k)*var%im(i,j+1,k) < 0.0)) then
            denom(i,j,k) = var%im(i,j,k)-var%im(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -var%im(i,j+1,k)*y(j)/denom(i,j,k) + &
                    var%im(i,j,k)*y(j+1)/denom(i,j,k)
            write (17, '(5e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (17, *) "# i,j,k --> i,j,k+1"

      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
            if ((var%im(i,j,k) == 0.0) .and. &
                (var%im(i,j,k)*var%im(i,j,k+1) < 0.0)) then
              denom(i,j,k) = var%im(i,j,k)-var%im(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -var%im(i,j,k+1)*z(k)/denom(i,j,k) + &
                      var%im(i,j,k)*z(k+1)/denom(i,j,k)
              write (17, '(5e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if

    close (16)
    close (17)

    return
  end subroutine get_re_im_zeros
  
  subroutine get_phase_zeros(in_var, p)
    use parameters
    use variables, only : get_phase
    use ic, only : x, y, z
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    integer, intent(in) :: p
    real :: zero
    real, dimension(0:nx,0:ny,0:nz) :: denom, phase
    integer :: i, j, k
    integer, parameter :: z_start=nz/2, z_end=nz/2
    !integer, parameter :: z_start=1, z_end=nz1

    open (18, status = 'unknown', file = 'phase_zeros'//itos(p)//'.dat')

    call get_phase(in_var, phase)

    write (18, *) "# i,j,k --> i+1,j,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((phase(i,j,k) == 0.0) .or. &
              (phase(i,j,k)*phase(i+1,j,k) < 0.0)) then
            denom(i,j,k) = phase(i,j,k)-phase(i+1,j,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -phase(i+1,j,k)*x(i)/denom(i,j,k) + &
                    phase(i,j,k)*x(i+1)/denom(i,j,k)
            write (18, '(5e17.9)') zero, y(j), z(k)
          end if
        end do
      end do
    end do
    
    write (18, *) "# i,j,k --> i,j+1,k"

    do k=z_start, z_end
      do j=1,ny1
        do i=1,nx1
          if ((phase(i,j,k) == 0.0) .or. &
              (phase(i,j,k)*phase(i,j+1,k) < 0.0)) then
            denom(i,j,k) = phase(i,j,k)-phase(i,j+1,k)
            if (denom(i,j,k) == 0.0) cycle
            zero = -phase(i,j+1,k)*y(j)/denom(i,j,k) + &
                    phase(i,j,k)*y(j+1)/denom(i,j,k)
            write (18, '(5e17.9)') x(i), zero, z(k)
          end if
        end do
      end do
    end do
    
    if (z_start /= z_end) then
      write (18, *) "# i,j,k --> i,j,k+1"

      do k=z_start, z_end
        do j=1,ny1
          do i=1,nx1
            if ((phase(i,j,k) == 0.0) .and. &
                (phase(i,j,k)*phase(i,j,k+1) < 0.0)) then
              denom(i,j,k) = phase(i,j,k)-phase(i,j,k+1)
              if (denom(i,j,k) == 0.0) cycle
              zero = -phase(i,j,k+1)*z(k)/denom(i,j,k) + &
                      phase(i,j,k)*z(k+1)/denom(i,j,k)
              write (18, '(5e17.9)') x(i), y(j), zero
            end if
          end do
        end do
      end do
    end if
    
    close (18)

    return
  end subroutine get_phase_zeros

  subroutine save_linelength(t, in_var)
    use parameters
    use variables, only : linelength
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, intent(in) :: t

    write (20, '(2e17.9)') t, linelength(t, in_var)

    return
  end subroutine save_linelength

end module io
