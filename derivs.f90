module derivs
  ! Routines to calculate finite-difference derivatives
  use parameters
  implicit none

  private
  public :: deriv_x, deriv_y, deriv_z, deriv_xx, deriv_yy, deriv_zz

  contains

  subroutine deriv_x(f,fx)
    ! First x-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fx
    integer :: i
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do i=0,nx
          call second_index(i, nx, minus, plus)

          fx(i,:,:) = ( f(plus(1),:,:) - f(minus(1),:,:) ) / (2.0*dx)
        end do
      case (4)
        ! Fourth order centred difference
        do i=0,nx
          call fourth_index(i, nx, minus, plus)

          fx(i,:,:) = ( -f(plus(2),:,:) + &
                     8.0*f(plus(1),:,:) - &
                     8.0*f(minus(1),:,:) + &
                         f(minus(2),:,:) ) / (12.0*dx)
        end do
    end select

    return
  end subroutine deriv_x
  
  subroutine deriv_y(f,fy)
    ! First y-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fy
    integer :: j
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do j=0,ny
          call second_index(j, ny, minus, plus)

          fy(:,j,:) = ( f(:,plus(1),:) - f(:,minus(1),:) ) / (2.0*dy)
        end do
      case (4)
        ! Fourth order centred difference
        do j=0,ny
          call fourth_index(j, ny, minus, plus)

          fy(:,j,:) = ( -f(:,plus(2),:) + &
                     8.0*f(:,plus(1),:) - &
                     8.0*f(:,minus(1),:) + &
                         f(:,minus(2),:) ) / (12.0*dy)
        end do
    end select

    return
  end subroutine deriv_y
  
  subroutine deriv_z(f,fz)
    ! First y-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fz
    integer :: k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=0,nz
          call second_index(k, nz, minus, plus)

          fz(:,:,k) = ( f(:,:,plus(1)) - f(:,:,minus(1)) ) / (2.0*dz)
        end do
      case (4)
        ! Fourth order centred difference
        do k=0,nz
          call fourth_index(k, nz, minus, plus)

          fz(:,:,k) = ( -f(:,:,plus(2)) + &
                     8.0*f(:,:,plus(1)) - &
                     8.0*f(:,:,minus(1)) + &
                         f(:,:,minus(2)) ) / (12.0*dz)
        end do
    end select

    return
  end subroutine deriv_z

  subroutine deriv_xx(f,fxx)
    ! Second x-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fxx
    integer :: i
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do i=0,nx
          call second_index(i, nx, minus, plus)

          fxx(i,:,:) = ( f(plus(1),:,:) - &
                         2.0*f(i,:,:) + &
                         f(minus(1),:,:) ) / dx2
        end do
      case (4)
        ! Fourth order centred difference
        do i=0,nx
          call fourth_index(i, nx, minus, plus)

          fxx(i,:,:) = ( -f(plus(2),:,:) + &
                     16.0*f(plus(1),:,:) - &
                     30.0*f(i,:,:) + &
                     16.0*f(minus(1),:,:) - &
                          f(minus(2),:,:) ) / (12.0*dx2)
        end do
    end select

    return
  end subroutine deriv_xx

  subroutine deriv_yy(f,fyy)
    ! Second y-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fyy
    integer :: j
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do j=0,ny
          call second_index(j, ny, minus, plus)

          fyy(:,j,:) = ( f(:,plus(1),:) - &
                         2.0*f(:,j,:) + &
                         f(:,minus(1),:) ) / dy2
        end do
      case (4)
        ! Fourth order centred difference
        do j=0,ny
          call fourth_index(j, ny, minus, plus)
          
          fyy(:,j,:) = ( -f(:,plus(2),:) + &
                     16.0*f(:,plus(1),:) - &
                     30.0*f(:,j,:) + &
                     16.0*f(:,minus(1),:) - &
                          f(:,minus(2),:) ) / (12.0*dy2)
        end do
    end select

    return
  end subroutine deriv_yy
  
  subroutine deriv_zz(f,fzz)
    ! Second z-derivative
    use parameters
    implicit none

    complex, dimension(0:nx,0:ny,0:nz), intent(in)  :: f
    complex, dimension(0:nx,0:ny,0:nz), intent(out) :: fzz
    integer :: k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=0,nz
          call second_index(k, nz, minus, plus)

          fzz(:,:,k) = ( f(:,:,plus(1)) - &
                         2.0*f(:,:,k) + &
                         f(:,:,minus(1)) ) / dz2
        end do
      case (4)
        ! Fourth order centred difference
        do k=0,nz
          call fourth_index(k, nz, minus, plus)
          
          fzz(:,:,k) = ( -f(:,:,plus(2)) + &
                     16.0*f(:,:,plus(1)) - &
                     30.0*f(:,:,k) + &
                     16.0*f(:,:,minus(1)) - &
                          f(:,:,minus(2)) ) / (12.0*dz2)
        end do
    end select

    return
  end subroutine deriv_zz

  subroutine second_index(indx, n, minus, plus)
    ! Determine the indices at the boundaries depending on whether periodic or
    ! reflective boundaries are chosen in the case of second order differences
    implicit none

    integer, intent(in)  :: indx, n
    integer, dimension(2), intent(out) :: minus, plus
    
    select case (bcs)
      case (1)
        ! periodic BCs
        if (indx == 0) then
          minus(1) = n-1
          plus(1) = 1
        else if (indx == n) then
          minus(1) = n-1
          plus(1) = 1
        else
          minus(1) = indx-1
          plus(1) = indx+1
        end if
      case (2)
        ! reflective BCs
        if (indx == 0) then
          minus(1) = 1
          plus(1) = 1
        else if (indx == n) then
          minus(1) = n-1
          plus(1) = n-1
        else
          minus(1) = indx-1
          plus(1) = indx+1
        end if
    end select

    return
  end subroutine second_index
  
  subroutine fourth_index(indx, n, minus, plus)
    ! Determine the indices at the boundaries depending on whether periodic or
    ! reflective boundaries are chosen in the case of fourth order differences
    implicit none

    integer, intent(in)  :: indx, n
    integer, dimension(2), intent(out) :: minus, plus
    
    select case (bcs)
      case (1)
        ! periodic BCs
        if (indx == 0) then
          minus(1) = n-1
          minus(2) = n-2
          plus(1) = 1
          plus(2) = 2
        else if (indx == 1) then
          minus(1) = 0
          minus(2) = n-1
          plus(1) = 2
          plus(2) = 3
        else if (indx == n-1) then
          minus(1) = n-2
          minus(2) = n-3
          plus(1) = n
          plus(2) = 1
        else if (indx == n) then
          minus(1) = n-1
          minus(2) = n-2
          plus(1) = 1
          plus(2) = 2
        else
          minus(1) = indx-1
          minus(2) = indx-2
          plus(1) = indx+1
          plus(2) = indx+2
        end if
      case (2)
        ! reflective BCs
        if (indx == 0) then
          minus(1) = 1
          minus(2) = 2
          plus(1) = 1
          plus(2) = 2
        else if (indx == 1) then
          minus(1) = 0
          minus(2) = 1
          plus(1) = 2
          plus(2) = 3
        else if (indx == n-1) then
          minus(1) = n-2
          minus(2) = n-3
          plus(1) = n
          plus(2) = n-1
        else if (indx == n) then
          minus(1) = n-1
          minus(2) = n-2
          plus(1) = n-1
          plus(2) = n-2
        else
          minus(1) = indx-1
          minus(2) = indx-2
          plus(1) = indx+1
          plus(2) = indx+2
        end if
    end select

    return
  end subroutine fourth_index

end module derivs
