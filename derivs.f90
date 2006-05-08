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

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fx
    integer :: i, j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            do i=0,nx1
              call second_index(i, nx1, minus, plus, .true.)

              fx(i,j,k) = ( f(plus(1),j,k) - f(minus(1),j,k) ) / (2.0*dx)
            end do
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            do i=0,nx1
              call fourth_index(i, nx1, minus, plus, .true.)

              fx(i,j,k) = ( -f(plus(2),j,k) + &
                         8.0*f(plus(1),j,k) - &
                         8.0*f(minus(1),j,k) + &
                             f(minus(2),j,k) ) / (12.0*dx)
            end do
          end do
        end do
    end select

    return
  end subroutine deriv_x
  
  subroutine deriv_y(f,fy)
    ! First y-derivative
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fy
    integer :: j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call second_index(j, ny1, minus, plus, .false.)

            fy(:,j,k) = ( f(:,plus(1),k) - f(:,minus(1),k) ) / (2.0*dy)
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call fourth_index(j, ny1, minus, plus, .false.)

            fy(:,j,k) = ( -f(:,plus(2),k) + &
                       8.0*f(:,plus(1),k) - &
                       8.0*f(:,minus(1),k) + &
                           f(:,minus(2),k) ) / (12.0*dy)
          end do
        end do
    end select

    return
  end subroutine deriv_y
  
  subroutine deriv_z(f,fz)
    ! First z-derivative
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fz
    integer :: j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call second_index(k, nz1, minus, plus, .false.)

            fz(:,j,k) = ( f(:,j,plus(1)) - f(:,j,minus(1)) ) / (2.0*dz)
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call fourth_index(k, nz1, minus, plus, .false.)

            fz(:,j,k) = ( -f(:,j,plus(2)) + &
                       8.0*f(:,j,plus(1)) - &
                       8.0*f(:,j,minus(1)) + &
                           f(:,j,minus(2)) ) / (12.0*dz)
          end do
        end do
    end select

    return
  end subroutine deriv_z

  subroutine deriv_xx(f,fxx)
    ! Second x-derivative
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fxx
    integer :: i, j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            do i=0,nx1
              call second_index(i, nx1, minus, plus, .true.)

              fxx(i,j,k) = ( f(plus(1),j,k) - &
                         2.0*f(i,j,k) + &
                             f(minus(1),j,k) ) / dx2
            end do
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            do i=0,nx1
              call fourth_index(i, nx1, minus, plus, .true.)

              fxx(i,j,k) = ( -f(plus(2),j,k) + &
                         16.0*f(plus(1),j,k) - &
                         30.0*f(i,j,k) + &
                         16.0*f(minus(1),j,k) - &
                              f(minus(2),j,k) ) / (12.0*dx2)
            end do
          end do
        end do
    end select

    return
  end subroutine deriv_xx

  subroutine deriv_yy(f,fyy)
    ! Second y-derivative
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fyy
    integer :: j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call second_index(j, ny1, minus, plus, .false.)

            fyy(:,j,k) = ( f(:,plus(1),k) - &
                       2.0*f(:,j,k) + &
                           f(:,minus(1),k) ) / dy2
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call fourth_index(j, ny1, minus, plus, .false.)
          
            fyy(:,j,k) = ( -f(:,plus(2),k) + &
                       16.0*f(:,plus(1),k) - &
                       30.0*f(:,j,k) + &
                       16.0*f(:,minus(1),k) - &
                            f(:,minus(2),k) ) / (12.0*dy2)
          end do
        end do
    end select

    return
  end subroutine deriv_yy
  
  subroutine deriv_zz(f,fzz)
    ! Second z-derivative
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: f
    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: fzz
    integer :: j, k
    integer, dimension(2) :: minus=0, plus=0

    select case (order)
      case (2)
        ! Second order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call second_index(k, nz1, minus, plus, .false.)

            fzz(:,j,k) = ( f(:,j,plus(1)) - &
                       2.0*f(:,j,k) + &
                           f(:,j,minus(1)) ) / dz2
          end do
        end do
      case (4)
        ! Fourth order centred difference
        do k=ksta,kend
          do j=jsta,jend
            call fourth_index(k, nz1, minus, plus, .false.)
          
            fzz(:,j,k) = ( -f(:,j,plus(2)) + &
                       16.0*f(:,j,plus(1)) - &
                       30.0*f(:,j,k) + &
                       16.0*f(:,j,minus(1)) - &
                            f(:,j,minus(2)) ) / (12.0*dz2)
          end do
        end do
    end select

    return
  end subroutine deriv_zz

  subroutine second_index(indx, n, minus, plus, x_deriv)
    ! Determine the indices at the boundaries depending on whether periodic or
    ! reflective boundaries are chosen in the case of second order differences
    implicit none

    integer, intent(in) :: indx, n
    logical, intent(in) :: x_deriv
    integer, dimension(2), intent(out) :: minus, plus
    
    select case (bcs)
      case (1)
        ! periodic BCs
        minus(1) = indx-1
        plus(1) = indx+1
        if (x_deriv) then
          ! Only need to define BCs for the x-direction which is not
          ! calculated in parallel
          if (indx == 0) then
            minus(1) = n
            plus(1) = 1
          end if
          if (indx == n) then
            minus(1) = n-1
            plus(1) = 0
          end if
        end if
      case (2)
        ! reflective BCs
        minus(1) = indx-1
        plus(1) = indx+1
        if (indx == 0) then
          minus(1) = 1
          plus(1) = 1
        end if
        if (indx == n) then
          minus(1) = n-1
          plus(1) = n-1
        end if
    end select

    return
  end subroutine second_index
  
  subroutine fourth_index(indx, n, minus, plus, x_deriv)
    ! Determine the indices at the boundaries depending on whether periodic or
    ! reflective boundaries are chosen in the case of fourth order differences
    implicit none

    integer, intent(in) :: indx, n
    logical, intent(in) :: x_deriv
    integer, dimension(2), intent(out) :: minus, plus
    
    select case (bcs)
      case (1)
        ! periodic BCs
        minus(1) = indx-1
        minus(2) = indx-2
        plus(1) = indx+1
        plus(2) = indx+2
        if (x_deriv) then
          ! Only need to define BCs for the x-direction which is not
          ! calculated in parallel
          if (indx == 0) then
            minus(1) = n
            minus(2) = n-1
            plus(1) = 1
            plus(2) = 2
          end if
          if (indx == 1) then
            minus(1) = 0
            minus(2) = n
            plus(1) = 2
            plus(2) = 3
          end if
          if (indx == n-1) then
            minus(1) = n-2
            minus(2) = n-3
            plus(1) = n
            plus(2) = 0
          end if
          if (indx == n) then
            minus(1) = n-1
            minus(2) = n-2
            plus(1) = 0
            plus(2) = 1
          end if
        end if
      case (2)
        ! reflective BCs
        minus(1) = indx-1
        minus(2) = indx-2
        plus(1) = indx+1
        plus(2) = indx+2
        if (indx == 0) then
          minus(1) = 1
          minus(2) = 2
          plus(1) = 1
          plus(2) = 2
        end if
        if (indx == 1) then
          minus(1) = 0
          minus(2) = 1
          plus(1) = 2
          plus(2) = 3
        end if
        if (indx == n-1) then
          minus(1) = n-2
          minus(2) = n-3
          plus(1) = n
          plus(2) = n-1
        end if
        if (indx == n) then
          minus(1) = n-1
          minus(2) = n-2
          plus(1) = n-1
          plus(2) = n-2
        end if
    end select

    return
  end subroutine fourth_index

end module derivs
