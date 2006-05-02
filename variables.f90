module variables
  ! Routines to do with setting up variables and operations on them
  use parameters
  implicit none

  private
  public :: laplacian, get_density, get_phase, get_norm, get_U, &
            energy, momentum, linelength

  type, public :: var
    complex, dimension(0:nx,0:ny,0:nz) :: new
    complex, dimension(0:nx,0:ny,0:nz) :: old
  end type var

  type, public :: deriv
    complex, dimension(0:nx,0:ny,0:nz) :: x
    complex, dimension(0:nx,0:ny,0:nz) :: y
    complex, dimension(0:nx,0:ny,0:nz) :: z
    complex, dimension(0:nx,0:ny,0:nz) :: xx
    complex, dimension(0:nx,0:ny,0:nz) :: yy
    complex, dimension(0:nx,0:ny,0:nz) :: zz
  end type deriv

  type, public :: re_im
    real, dimension(0:nx,0:ny,0:nz) :: re
    real, dimension(0:nx,0:ny,0:nz) :: im
  end type re_im
  
  real, parameter, private   :: c1 = 3.0/8.0, &
                                c2 = 7.0/6.0, &
                                c3 = 23.0/24.0

  !complex, dimension(0:nx,0:ny,0:nz), save, public :: U

  contains

  function laplacian(in_var)
    ! Laplacian in cartesian coordinates
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    complex, dimension(0:nx,0:ny,0:nz) :: laplacian
    type (deriv) :: d

    call deriv_xx(in_var,d%xx)
    call deriv_yy(in_var,d%yy)
    call deriv_zz(in_var,d%zz)
    
    laplacian = d%xx + d%yy + d%zz

    return
  end function laplacian

  subroutine get_phase(in_var, phase)
    ! Phase
    use parameters
    implicit none

    complex, dimension(:,:,:), intent(in) :: in_var
    real, dimension(0:nx,0:ny,0:nz), intent(out) :: phase

    phase = atan2(aimag(in_var)+1.0e-6, real(in_var))

    return
  end subroutine get_phase

  subroutine get_density(in_var, density)
    ! Density
    use parameters
    implicit none

    complex, dimension(:,:,:), intent(in) :: in_var
    real, dimension(0:nx,0:ny,0:nz), intent(out) :: density

    density = abs(in_var)

    return
  end subroutine get_density

  subroutine get_norm(in_var, norm)
    ! Calculate the norm
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, intent(out) :: norm
    real :: int_z
    real, dimension(0:nx,0:ny,0:nz) :: int_var
    real, dimension(0:ny,0:nz) :: int_x
    real, dimension(0:nz) :: int_y
    
    int_var = abs(in_var)**2
    
    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    norm = int_z

    return
  end subroutine get_norm

  subroutine get_U(in_var, U)
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: U
    type (deriv) :: d

    call deriv_x(in_var, d%x)
    call deriv_xx(in_var, d%xx)
    call deriv_yy(in_var, d%yy)

    U = -0.5*eye * (d%xx + d%yy + (1.0-abs(in_var)**2)*in_var ) / d%x
    print*, U
    stop

    return
  end subroutine get_U

  subroutine energy(in_var, E)
    ! Calculate the energy
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, intent(out) :: E
    real :: int_z
    real, dimension(0:nx,0:ny,0:nz) :: int_var
    real, dimension(0:ny,0:nz) :: int_x
    real, dimension(0:nz) :: int_y
    complex, dimension(0:nx,0:ny,0:nz) :: dpsidx

    call deriv_x(in_var, dpsidx)
    
    int_var = abs(dpsidx)**2

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    E = 0.5 * int_z

    return
  end subroutine energy
  
  subroutine momentum(in_var, P)
    ! Calculate the momentum
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, dimension(3), intent(out) :: P
    real :: int_z
    real, dimension(0:nx,0:ny,0:nz) :: int_var
    real, dimension(0:ny,0:nz) :: int_x
    real, dimension(0:nz) :: int_y
    type (deriv) :: dpsi, dpsistar

    call deriv_x(in_var, dpsi%x)
    call deriv_y(in_var, dpsi%y)
    call deriv_z(in_var, dpsi%z)
    call deriv_x(conjg(in_var), dpsistar%x)
    call deriv_y(conjg(in_var), dpsistar%y)
    call deriv_z(conjg(in_var), dpsistar%z)
    
    int_var = (real(in_var)-1.0)*aimag(dpsistar%x) - &
               aimag(in_var)*real(dpsi%x)

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(1) = int_z
    
    int_var = (real(in_var)-1.0)*aimag(dpsistar%y) - &
               aimag(in_var)*real(dpsi%y)

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(2) = int_z
    
    int_var = (real(in_var)-1.0)*aimag(dpsistar%z) - &
               aimag(in_var)*real(dpsi%z)

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(3) = int_z

    return
  end subroutine momentum

  subroutine integrate_x(in_var, x_int)
    ! Integrate a (3D) variable in x
    use parameters
    implicit none

    real, dimension(0:nx,0:ny,0:nz), intent(in) :: in_var
    real, dimension(0:ny,0:nz), intent(out) :: x_int
    integer :: j, k

    do k =0,nz
      do j=0,ny
        x_int(j,k) = (c1*in_var(0,j,k) + &
                      c2*in_var(1,j,k) + &
                      c3*in_var(2,j,k) + &
                      sum(in_var(3:nx-3,j,k)) + &
                      c3*in_var(nx-2,j,k) + &
                      c2*in_var(nx-1,j,k) + &
                      c1*in_var(nx,j,k)) * dx
      end do
    end do

    return
  end subroutine integrate_x
  
  subroutine integrate_y(in_var, y_int)
    ! Integrate a (2D) variable in y
    use parameters
    implicit none

    real, dimension(0:ny,0:nz), intent(in) :: in_var
    real, dimension(0:nz), intent(out) :: y_int
    integer :: k

    do k =0,nz
      y_int(k) = (c1*in_var(0,k) + &
                  c2*in_var(1,k) + &
                  c3*in_var(2,k) + &
                  sum(in_var(3:ny-3,k)) + &
                  c3*in_var(ny-2,k) + &
                  c2*in_var(ny-1,k) + &
                  c1*in_var(ny,k)) * dy
    end do

    return
  end subroutine integrate_y
  
  subroutine integrate_z(in_var, z_int)
    ! Integrate a (1D) variable in z
    use parameters
    implicit none

    real, dimension(0:nz), intent(in) :: in_var
    real, intent(out) :: z_int

    z_int = (c1*in_var(0) + &
             c2*in_var(1) + &
             c3*in_var(2) + &
             sum(in_var(3:nz-3)) + &
             c3*in_var(nz-2) + &
             c2*in_var(nz-1) + &
             c1*in_var(nz)) * dz

    return
  end subroutine integrate_z

  function linelength(time, psi)
    use parameters
    implicit none
  
    real :: linelength
    complex, intent(in) :: psi(0:nx,0:ny,0:nz)
    real, intent(in) :: time
    real :: xx, yy, h, a, den
    integer :: i, j, k, n, m, p, lv, lu
    real :: x(6), y(6), z(6), u(5), v(5), xu(6), xv(6), yu(6), yv(6)
  
  ! MUST have dx=dy=dz
  
    linelength = 0.0
    !do i=1,nx1
    !  do j=1,ny1
    !    do k=1,nz1
    do k=1,nz1
      do j=1,ny1
        iloop: do i=1,nx1
          m=1
          nloop: do n=1,6
  ! 1-front,2-back, 3-right, 4-left, 5-top, 6-bottom
            select case (n)
              case (1)
                u(1) = real(psi(i,j,k))
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i+1,j,k))
                v(2) = aimag(psi(i+1,j,k))
                u(4) = real(psi(i,j+1,k))
                v(4) = aimag(psi(i,j+1,k))
                u(3) = real(psi(i+1,j+1,k))
                v(3) = aimag(psi(i+1,j+1,k))
              case (2) 
                u(1) = real(psi(i,j,k+1))
                v(1) = aimag(psi(i,j,k+1))
                u(2) = real(psi(i+1,j,k+1))
                v(2) = aimag(psi(i+1,j,k+1))
                u(4) = real(psi(i,j+1,k+1))
                v(4) = aimag(psi(i,j+1,k+1))
                u(3) = real(psi(i+1,j+1,k+1))
                v(3) = aimag(psi(i+1,j+1,k+1))
              case(3)
                u(1) = real(psi(i+1,j,k))
                v(1) = aimag(psi(i+1,j,k))
                u(2) = real(psi(i+1,j,k+1))
                v(2) = aimag(psi(i+1,j,k+1))
                u(4) = real(psi(i+1,j+1,k))
                v(4) = aimag(psi(i+1,j+1,k))
                u(3) = real(psi(i+1,j+1,k+1))
                v(3) = aimag(psi(i+1,j+1,k+1))
              case (4)
                u(1) = real(psi(i,j,k))
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i,j,k+1))
                v(2) = aimag(psi(i,j,k+1))
                u(4) = real(psi(i,j+1,k))
                v(4) = aimag(psi(i,j+1,k))
                u(3) = real(psi(i,j+1,k+1))
                v(3) = aimag(psi(i,j+1,k+1))
              case (5)
                u(1) = real(psi(i,j+1,k))
                v(1) = aimag(psi(i,j+1,k))
                u(2) = real(psi(i+1,j+1,k))
                v(2) = aimag(psi(i+1,j+1,k))
                u(4) = real(psi(i,j+1,k+1))
                v(4) = aimag(psi(i,j+1,k+1))
                u(3) = real(psi(i+1,j+1,k+1))
                v(3) = aimag(psi(i+1,j+1,k+1))
              case (6)
                u(1) = real(psi(i,j,k))
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i+1,j,k))
                v(2) = aimag(psi(i+1,j,k))
                u(4) = real(psi(i,j,k+1))
                v(4) = aimag(psi(i,j,k+1))
                u(3) = real(psi(i+1,j,k+1))
                v(3) = aimag(psi(i+1,j,k+1))
            end select
            u(5) = u(1)
            v(5) = v(1)
  
  ! find zero line in u and v
  ! deal with planes first
            if ((u(1)==0.0 .and. u(2)==0.0 .and. u(3)==0.0 .and. &
                 u(4)==0.0 .and. &
                (v(1)*v(2)<0.0 .or. v(2)*v(3)<0.0 .or. &
                 v(3)*v(4)<0.0 .or. v(4)*v(1)<0.0)) .or. &
                (v(1)==0.0 .and. v(2)==0.0 .and. v(3)==0.0 .and. &
                 v(4)==0.0 .and. &
                (u(1)*u(2)<0.0 .or. u(2)*u(3)<0.0 .or. &
                 u(3)*u(4)<0.0 .or. u(4)*u(1)<0.0))) then
              linelength = linelength+0.5
              print*, 'cycle iloop'
              a = 1.0
              cycle iloop
              !goto 200
            end if
  ! deal with edge
            do p=1,4
              if(u(p)==0.0 .and. u(p+1)==0.0 .and. &
                 v(p)==0.0 .and. v(p+1)==0.0) then
                linelength = linelength+0.25
                print*, 'cycle iloop'
                a = 1.0
                cycle iloop
                !goto 200
              end if
            end do
  
            lu=1
            do p=1,4
              if(u(p)==0.0 .and. u(p+1)==0.0 .and. &
                 v(p)==0.0 .and. v(p+1)==0.0) then
                m=m+1
                print*, 'exit nloop'
                exit nloop
                !goto 100
              else if (u(p)==0.0) then
                select case (p)
                  case (1)
                    xu(lu)=0.0
                    yu(lu)=0.0
                  case (2)
                    xu(lu)=1.0
                    yu(lu)=0.0
                  case (3)
                    xu(lu)=1.0
                    yu(lu)=1.0
                  case (4)
                    xu(lu)=0.0
                    yu(lu)=1.0
                end select
              else if (u(p)*u(p+1)<0.0) then
                select case (p)
                  case (1)
                    xu(lu) = abs(u(p)/(u(p)-u(p+1)))
                    yu(lu) = 0.0
                  case (2)
                    xu(lu) = 1.0
                    yu(lu) = abs(u(p)/(u(p)-u(p+1)))
                  case (3)
                    xu(lu) = abs(u(p+1)/(u(p)-u(p+1)))
                    yu(lu) = 1.0
                  case (4)
                    xu(lu) = 0.0
                    yu(lu) = abs(u(p+1)/(u(p)-u(p+1)))
                end select
  
                if (lu==1 .or. xu(lu)/=xu(lu-1) .or. yu(lu)/=yu(lu-1)) then
                  lu = lu+1
                end if
              end if
            end do
            lv = 1
            do p=1,4
              if (v(p)==0.0) then
                select case (p)
                  case (1)
                    xv(lv) = 0.0
                    yv(lv) = 0.0
                  case (2)
                    xv(lv) = 1.0
                    yv(lv) = 0.0
                  case (3)
                    xv(lv) = 1.0
                    yv(lv) = 1.0
                  case (4)
                    xv(lv) = 0.0
                    yv(lv) = 1.0
                end select
              else if (v(p)*v(p+1)<0.0) then
                select case (p)
                  case (1)
                    xv(lv) = abs(v(p)/(v(p)-v(p+1)))
                    yv(lv) = 0.0
                  case (2)
                    xv(lv) = 1.0
                    yv(lv) = abs(v(p)/(v(p)-v(p+1)))
                  case (3)
                    xv(lv) = abs(v(p+1)/(v(p)-v(p+1)))
                    yv(lv) = 1.0
                  case (4)
                    xv(lv) = 0.0
                    yv(lv) = abs(v(p+1)/(v(p)-v(p+1)))
                end select
  
                if (lv==1 .or. xv(lv)/=xv(lv-1) .or. yv(lv)/=yv(lv-1)) then
                  lv = lv+1
                end if
              end if
            end do ! in p
            if (lu>2 .and. lv>2) then
              den = xv(2)*(yu(1)-yu(2))+xv(1)*(yu(2)-yu(1))+ &
                   (xu(1)-xu(2))*(yv(1)-yv(2))
              !if (den==0.0)  then
              !  write (*,*) i, j, k, xu(1), yu(1), xu(2), yu(2), &
              !              xv(1), yv(1), xv(2), yv(2)
              !  print*, 'ZERO DENOM IN linelength'
              !end if
              xx = (xu(1)*(xv(2)*(yv(1)-yu(2))+xv(1)*(yu(2)-yv(2)))+ &
                    xu(2)*(xv(2)*(yu(1)-yv(1))+xv(1)*(yv(2)-yu(1))))/den
              yy = (xv(2)*(yu(1)-yu(2))*yv(1)+xu(1)*yu(2)*yv(1)- &
                    xv(1)*yu(1)*yv(2)-xu(1)*yu(2)*yv(2)+xv(1)*yu(2)*yv(2)+ &
                    xu(2)*yu(1)*(yv(2)-yv(1)))/den
  
              if (xx>=0.0 .and. xx<=1.0 .and. yy>=0.0 .and. yy<=1.0) then
  ! found zero inside the square
                select case (n)
                  case (1)
                    x(m) = xx
                    y(m) = yy
                    z(m) = 0.0
                  case (2)
                    x(m) = xx
                    y(m) = yy
                    z(m) = 1.0
                  case (3)
                    x(m) = 1.0
                    y(m) = yy
                    z(m) = xx
                  case (4)
                    x(m) = 0.0
                    y(m) = yy
                    z(m) = xx
                  case (5)
                    x(m) = xx
                    y(m) = 1.0
                    z(m) = yy
                  case (6)
                    x(m) = xx
                    y(m) = 0.0
                    z(m) = yy
                end select
  
                !if (k==41 .and. i==41 .and. j==41) then
                !  print*, 'together=', n, m, den, x(m), y(m), z(m)
                !end if
                m=m+1
              end if
            end if
          end do nloop
  100     if (m>1) then      ! found zero in at least two sides
            a = 1.0 !scale
            if (x(1)==x(2) .and. (x(1)==0.0 .or. x(1)==1.0) .or. &
                y(1)==y(2) .and. (y(1)==0.0 .or. y(1)==1.0) .or. &
                y(1)==y(2) .and. (y(1)==0.0 .or. y(1)==1.0)) then
              a = a*0.5
              if ((x(1)==x(2) .and. (x(1)==0.0 .or. x(1)==1.0) .and. &
                  (y(1)==y(2) .and. (y(1)==0.0 .or. y(1)==1.0) .or. &
                   y(1)==y(2) .and. (y(1)==0.0 .or. y(1)==1.0))) .or. &
                  (y(1)==y(2) .and. (y(1)==0.0 .or. y(1)==1.0) .and. &
                  (z(1)==z(2) .and. (z(1)==0.0 .or. z(1)==1.0)))) then
                a=a*0.5
              end if
            end if
            linelength = linelength+a*sqrt((x(1)-x(2))**2+(y(1)-y(2))**2+ &
                                           (z(1)-z(2))**2)
          end if
  200     a = 1.0
        end do iloop
      end do
    end do
    linelength = linelength*dx
  
    return
  end function linelength
  
end module variables
