module variables
  ! Routines to do with setting up variables and operations on them
  use parameters
  implicit none

  private
  public :: laplacian, get_density, get_phase, get_norm, &
            energy, momentum, linelength, setup_itable, para_range, &
            array_len, neighbours, send_recv_y, send_recv_z, pack_y, &
            unpack_y, get_unit_no

  type, public :: var
    complex, allocatable, dimension(:,:,:) :: new
    complex, allocatable, dimension(:,:,:) :: old
    complex, allocatable, dimension(:,:,:) :: old2
  end type var

  type, public :: deriv
    complex, allocatable, dimension(:,:,:) :: x
    complex, allocatable, dimension(:,:,:) :: y
    complex, allocatable, dimension(:,:,:) :: z
    complex, allocatable, dimension(:,:,:) :: xx
    complex, allocatable, dimension(:,:,:) :: yy
    complex, allocatable, dimension(:,:,:) :: zz
  end type deriv

  type, public :: re_im
    real, allocatable, dimension(:,:,:) :: re
    real, allocatable, dimension(:,:,:) :: im
  end type re_im

  integer, dimension(-1:nyprocs, -1:nzprocs), public :: itable
  integer, public :: unit_no
  
  ! Constants for numerical integration
  real, parameter, private :: c1 = 3.0/8.0, &
                                c2 = 7.0/6.0, &
                                c3 = 23.0/24.0

  contains

  function laplacian(in_var)
    ! Laplacian in cartesian coordinates
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: laplacian
    type (deriv) :: d

    allocate(d%xx(0:nx1,jsta:jend,ksta:kend))
    allocate(d%yy(0:nx1,jsta:jend,ksta:kend))
    allocate(d%zz(0:nx1,jsta:jend,ksta:kend))

    call deriv_xx(in_var,d%xx)
    call deriv_yy(in_var,d%yy)
    call deriv_zz(in_var,d%zz)
    
    laplacian = d%xx + d%yy + d%zz

    deallocate(d%xx)
    deallocate(d%yy)
    deallocate(d%zz)
    
    return
  end function laplacian

  subroutine get_unit_no()
    ! Get the unit number which each process can write to
    implicit none

    unit_no = myrank+20

    return
  end subroutine get_unit_no

  subroutine setup_itable()
    ! Set up the lookup table for neighbouring processes
    use parameters
    implicit none

    integer :: j, k, irank

    ! Initially set each process position to null
    itable = MPI_PROC_NULL

    irank = 0

    ! Fill the lookup table
    do k=0,nzprocs-1
      do j=0,nyprocs-1
        itable(j,k) = irank
        if (myrank == irank) then
          myranky = j
          myrankz = k
        end if
        irank = irank+1
      end do
    end do

    ! If we have periodic BCs then fill in the boundary processes
    if (bcs == 1) then
      itable(-1,:) = itable(nyprocs-1,:)
      itable(nyprocs,:) = itable(0,:)
      itable(:,-1) = itable(:,nzprocs-1)
      itable(:,nzprocs) = itable(:,0)
    end if

    return
  end subroutine setup_itable

  subroutine para_range(n1, n2, nprocs, irank, ista, iend)
    ! Determine the start and end indices of the arrays on each process
    implicit none

    integer, intent(in) :: n1, n2, nprocs, irank
    integer, intent(out) :: ista, iend
    integer, dimension(2) :: iwork

    iwork(1) = (n2-n1+1)/nprocs
    iwork(2) = mod(n2-n1+1, nprocs)
    ista = irank*iwork(1)+n1+min(irank, iwork(2))
    iend = ista+iwork(1)-1
    if (iwork(2) > irank) iend = iend+1

    return
  end subroutine para_range

  subroutine array_len(jlen, klen)
    ! Determine the length of each array dimension on each process.  This 
    ! allows for the possibility that there are references to (for example) 
    ! j+1, k+1 simultaneously
    use parameters
    implicit none

    integer, intent(out) :: jlen, klen

    jlen = jend-jsta+1
    kksta = max(0, ksta-1)
    kkend = min(nz1, kend+1)
    klen = kkend-kksta+1

    return
  end subroutine array_len

  subroutine neighbours()
    ! Determine neighbouring processes for each process
    use parameters
    implicit none

    znext = itable(myranky, myrankz+1)
    zprev = itable(myranky, myrankz-1)
    ynext = itable(myranky+1, myrankz)
    yprev = itable(myranky-1, myrankz)

    return
  end subroutine neighbours

  subroutine send_recv_z(in_var)
    ! Send and receive boundary elements of each process to and from
    ! neighbouring processes in the z-direction.  This data is contiguous in
    ! memory and so can be immediately be sent
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    integer, dimension(4) :: jsend, jrecv

    ! Send the data two away from the boundary (to allow for fourth order
    ! derivatives)
    call MPI_ISEND(in_var(0,jsta,kend-1), nx*jlen, MPI_COMPLEX, znext, 1, &
                   MPI_COMM_WORLD, jsend(1), ierr)
    call MPI_ISEND(in_var(0,jsta,ksta+1), nx*jlen, MPI_COMPLEX, zprev, 1, &
                   MPI_COMM_WORLD, jsend(2), ierr)
                   
    ! Receive the data two away from the boundary
    call MPI_IRECV(in_var(0,jsta,ksta-2), nx*jlen, MPI_COMPLEX, zprev, 1, &
                   MPI_COMM_WORLD, jrecv(1), ierr)
    call MPI_IRECV(in_var(0,jsta,kend+2), nx*jlen, MPI_COMPLEX, znext, 1, &
                   MPI_COMM_WORLD, jrecv(2), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(jsend(1), istatus, ierr)
    call MPI_WAIT(jsend(2), istatus, ierr)
    call MPI_WAIT(jrecv(1), istatus, ierr)
    call MPI_WAIT(jrecv(2), istatus, ierr)

    ! Send the data one away from the boundary
    call MPI_ISEND(in_var(0,jsta,kend), nx*jlen, MPI_COMPLEX, znext, 1, &
                   MPI_COMM_WORLD, jsend(3), ierr)
    call MPI_ISEND(in_var(0,jsta,ksta), nx*jlen, MPI_COMPLEX, zprev, 1, &
                   MPI_COMM_WORLD, jsend(4), ierr)
                   
    ! Receive the data one away from the boundary
    call MPI_IRECV(in_var(0,jsta,ksta-1), nx*jlen, MPI_COMPLEX, zprev, 1, &
                   MPI_COMM_WORLD, jrecv(3), ierr)
    call MPI_IRECV(in_var(0,jsta,kend+1), nx*jlen, MPI_COMPLEX, znext, 1, &
                   MPI_COMM_WORLD, jrecv(4), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(jsend(3), istatus, ierr)
    call MPI_WAIT(jsend(4), istatus, ierr)
    call MPI_WAIT(jrecv(3), istatus, ierr)
    call MPI_WAIT(jrecv(4), istatus, ierr)

    return
  end subroutine send_recv_z

  subroutine send_recv_y()
    ! Send and receive boundary elements of each process to and from
    ! neighbouring processes in the y-direction.  This data is NOT contiguous
    ! in memory and so must first be packed into a contiguous array (see pack_y
    ! and unpack_y below)
    use parameters
    implicit none

    integer, dimension(2) :: ksend, krecv

    ! Send the boundary data to neighbouring processes
    call MPI_ISEND(works1(0,1,kksta), nx*2*klen, MPI_COMPLEX, ynext, 1, &
                   MPI_COMM_WORLD, ksend(1), ierr)
    call MPI_ISEND(works2(0,1,kksta), nx*2*klen, MPI_COMPLEX, yprev, 1, &
                   MPI_COMM_WORLD, ksend(2), ierr)

    ! Receive the boundary data from neighbouring processes
    call MPI_IRECV(workr1(0,1,kksta), nx*2*klen, MPI_COMPLEX, yprev, 1, &
                   MPI_COMM_WORLD, krecv(1), ierr)
    call MPI_IRECV(workr2(0,1,kksta), nx*2*klen, MPI_COMPLEX, ynext, 1, &
                   MPI_COMM_WORLD, krecv(2), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(ksend(1), istatus, ierr)
    call MPI_WAIT(ksend(2), istatus, ierr)
    call MPI_WAIT(krecv(1), istatus, ierr)
    call MPI_WAIT(krecv(2), istatus, ierr)

    return
  end subroutine send_recv_y

  subroutine pack_y(in_var)
    ! Pack the non-contiguous boundary data for the y-direction into a
    ! contiguous array
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    integer :: k

    if (myranky /= nyprocs-1) then
      do k=kksta,kkend
        works1(:,:,k) = in_var(:,jend-1:jend,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! pack their boundary data
      if (myranky == nyprocs-1) then
        do k=kksta,kkend
          works1(:,:,k) = in_var(:,jend-1:jend,k)
        end do
      end if
    end if

    if (myranky /= 0) then
      do k=kksta,kkend
        works2(:,:,k) = in_var(:,jsta:jsta+1,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! pack their boundary data
      if (myranky == 0) then
        do k=kksta,kkend
          works2(:,:,k) = in_var(:,jsta:jsta+1,k)
        end do
      end if
    end if

    return
  end subroutine pack_y

  subroutine unpack_y(in_var)
    ! Unpack the contiguous boundary data back into the non-contiguous array
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(out) :: in_var
    integer :: i, k
    
    if (myranky /= 0) then
      do k=kksta,kkend
        in_var(:,jsta-2:jsta-1,k) = workr1(:,:,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! unpack their boundary data
      if (myranky == 0) then
        do k=kksta,kkend
          in_var(:,jsta-2:jsta-1,k) = workr1(:,:,k)
        end do
      end if
    end if

    if (myranky /= nyprocs-1) then
      do k=kksta,kkend
        in_var(:,jend+1:jend+2,k) = workr2(:,:,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! unpack their boundary data
      if (myranky == nyprocs-1) then
        do k=kksta,kkend
          in_var(:,jend+1:jend+2,k) = workr2(:,:,k)
        end do
      end if
    end if

    return
  end subroutine unpack_y

  subroutine get_phase(in_var, phase)
    ! Calculate the phase
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: phase
    integer :: j, k

    phase = atan2(aimag(in_var)+1.0e-6, real(in_var))

    return
  end subroutine get_phase

  subroutine get_density(in_var, density)
    ! Calculate the density
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(out) :: density
    integer :: j, k

    density = abs(in_var)

    return
  end subroutine get_density

  subroutine get_norm(in_var, norm)
    ! Calculate the norm
    use parameters
    implicit none

    complex, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, intent(out) :: norm
    real :: int_z
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: int_var
    real, dimension(jsta:jend,ksta:kend) :: int_x
    real, dimension(ksta:kend) :: int_y
    integer :: j, k
    
    int_var = abs(in_var)**2
    
    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    norm = int_z

    return
  end subroutine get_norm

  subroutine energy(in_var, E)
    ! Calculate the energy
    use parameters
    use derivs
    implicit none

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    real, intent(out) :: E
    real :: int_z
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: int_var
    real, dimension(jsta:jend,ksta:kend) :: int_x
    real, dimension(ksta:kend) :: int_y
    complex, dimension(0:nx1,jsta:jend,ksta:kend) :: dpsidx
    integer :: j, k

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

    complex, dimension(0:nx1,jsta-2:jend+2,ksta-2:kend+2), intent(in) :: in_var
    real, dimension(3), intent(out) :: P
    real :: int_z
    real, dimension(0:nx1,jsta:jend,ksta:kend) :: int_var
    real, dimension(jsta:jend,ksta:kend) :: int_x
    real, dimension(ksta:kend) :: int_y
    type (deriv) :: dpsi, dpsistar
    integer :: j, k

    allocate(dpsi%x(0:nx1,jsta:jend,ksta:kend))
    allocate(dpsi%y(0:nx1,jsta:jend,ksta:kend))
    allocate(dpsi%z(0:nx1,jsta:jend,ksta:kend))
    allocate(dpsistar%x(0:nx1,jsta:jend,ksta:kend))
    allocate(dpsistar%y(0:nx1,jsta:jend,ksta:kend))
    allocate(dpsistar%z(0:nx1,jsta:jend,ksta:kend))
    
    call deriv_x(in_var, dpsi%x)
    call deriv_y(in_var, dpsi%y)
    call deriv_z(in_var, dpsi%z)
    call deriv_x(conjg(in_var), dpsistar%x)
    call deriv_y(conjg(in_var), dpsistar%y)
    call deriv_z(conjg(in_var), dpsistar%z)
    
    do k=ksta,kend
      do j=jsta,jend
        int_var(:,j,k) = (real(in_var(:,j,k))-1.0)*aimag(dpsistar%x(:,j,k)) - &
                          aimag(in_var(:,j,k))*real(dpsi%x(:,j,k))
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(1) = int_z
    
    do k=ksta,kend
      do j=jsta,jend
        int_var(:,j,k) = (real(in_var(:,j,k))-1.0)*aimag(dpsistar%y(:,j,k)) - &
                          aimag(in_var(:,j,k))*real(dpsi%y(:,j,k))
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(2) = int_z
    
    do k=ksta,kend
      do j=jsta,jend
        int_var(:,j,k) = (real(in_var(:,j,k))-1.0)*aimag(dpsistar%z(:,j,k)) - &
                          aimag(in_var(:,j,k))*real(dpsi%z(:,j,k))
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    P(3) = int_z

    deallocate(dpsi%x)
    deallocate(dpsi%y)
    deallocate(dpsi%z)
    deallocate(dpsistar%x)
    deallocate(dpsistar%y)
    deallocate(dpsistar%z)

    return
  end subroutine momentum

  subroutine integrate_x(in_var, x_int)
    ! Integrate a (3D) variable in x.  The x-direction is not parallelised so
    ! this integration is straight forward
    use parameters
    implicit none

    real, dimension(0:nx1,jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(jsta:jend,ksta:kend), intent(out) :: x_int
    integer :: j, k

    do k=ksta,kend
      do j=jsta,jend
        x_int(j,k) = (c1*in_var(0,j,k) + &
                      c2*in_var(1,j,k) + &
                      c3*in_var(2,j,k) + &
                      sum(in_var(3:nx-4,j,k)) + &
                      c3*in_var(nx-3,j,k) + &
                      c2*in_var(nx-2,j,k) + &
                      c1*in_var(nx-1,j,k)) * dx
      end do
    end do

    return
  end subroutine integrate_x
  
  subroutine integrate_y(in_var, y_int)
    ! Integrate a (2D) variable in y
    use parameters
    implicit none

    real, dimension(jsta:jend,ksta:kend), intent(in) :: in_var
    real, dimension(ksta:kend), intent(out) :: y_int
    real, dimension(jsta:jend,ksta:kend) :: tmp_var
    real, dimension(0:nz1) :: tmp, total
    integer :: j, k, irank

    ! Create a temporary variable on which to perform operations
    tmp_var = in_var
    
    ! Update the elements which must be multiplied by the integrating constants
    do j=jsta,jend
      if (j==0) tmp_var(j,:) = c1*in_var(j,:)
      if (j==1) tmp_var(j,:) = c2*in_var(j,:)
      if (j==2) tmp_var(j,:) = c3*in_var(j,:)
      if (j==ny-3) tmp_var(j,:) = c3*in_var(j,:)
      if (j==ny-2) tmp_var(j,:) = c2*in_var(j,:)
      if (j==ny-1) tmp_var(j,:) = c1*in_var(j,:)
    end do
      
    ! Sum the variable on individual processes
    tmp = 0.0
    do k=ksta,kend
      tmp(k) = sum(tmp_var(:,k))
    end do

    ! Sum the variable over all processes and make sure each process has the
    ! result
    call MPI_ALLREDUCE(tmp, total, nz, MPI_REAL, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

    ! Calculate the final integrated result
    do k=ksta,kend
      y_int(k) = total(k) * dy
    end do

    return
  end subroutine integrate_y
  
  subroutine integrate_z(in_var, z_int)
    ! Integrate a (1D) variable in z
    use parameters
    implicit none

    real, dimension(ksta:kend), intent(in) :: in_var
    real, intent(out) :: z_int
    real, dimension(ksta:kend) :: tmp_var
    real :: tmp
    integer :: k

    ! Create a temporary variable
    tmp_var = in_var
    
    ! Update the elements which must be multiplied by the integrating constants
    do k=ksta,kend
      if (k==0) tmp_var(k) = c1*in_var(k)
      if (k==1) tmp_var(k) = c2*in_var(k)
      if (k==2) tmp_var(k) = c3*in_var(k)
      if (k==nz-3) tmp_var(k) = c3*in_var(k)
      if (k==nz-2) tmp_var(k) = c2*in_var(k)
      if (k==nz-1) tmp_var(k) = c1*in_var(k)
    end do
    
    ! Calculate the sum on each individual process but DO NOT include those
    ! processes over which the variable is not distributed
    tmp = sum(tmp_var)
    if (myranky /= 0) tmp = 0.0

    ! Calculate the sum over all processes and make sure each process has the
    ! result
    call MPI_ALLREDUCE(tmp, z_int, 1, MPI_REAL, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)

    ! Calculate the final result
    z_int = z_int*dz

    return
  end subroutine integrate_z

  function linelength(time, psi)
    ! Calculate the total line length of vortices over the whole box
    ! individually on each process.  Serial algorithm written by Natalia
    ! Berloff in Fortran77.  Updated by me for Fortran 90
    use parameters
    implicit none
  
    real :: linelength
    complex, intent(in) :: psi(0:nx1,jsta-2:jend+2,ksta-2:kend+2)
    real, intent(in) :: time
    real :: xx, yy, h, a, den
    integer :: i, j, k, n, m, p, lv, lu
    real :: x(6), y(6), z(6), u(5), v(5), xu(6), xv(6), yu(6), yv(6)
  
  ! MUST have dx=dy=dz
  
    linelength = 0.0
    do k=ksta,kend
      if ((k==0) .or. (k==nz1)) cycle
      do j=jsta,jend
        if ((j==0) .or. (j==ny1)) cycle
        iloop: do i=1,nx1-1
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
              !print*, 'cycle iloop'
              a = 1.0
              cycle iloop
              !goto 200
            end if
  ! deal with edge
            do p=1,4
              if(u(p)==0.0 .and. u(p+1)==0.0 .and. &
                 v(p)==0.0 .and. v(p+1)==0.0) then
                linelength = linelength+0.25
                !print*, 'cycle iloop'
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
