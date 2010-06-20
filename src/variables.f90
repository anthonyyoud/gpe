! $Id$
!----------------------------------------------------------------------------

module variables
  ! Routines to do with setting up variables and operations on them
  use parameters
  implicit none

  private
  public :: laplacian, get_density, get_phase, get_pdf_velocity, get_norm, &
    get_pdf, get_vcf, energy, mass, momentum, linelength, setup_itable, &
    para_range, array_len, neighbours, send_recv_y, send_recv_z, pack_y, &
    unpack_y, renormalise, imprint_vortex_line

  type, public :: var
    complex (pr), allocatable, dimension(:,:,:) :: new
    complex (pr), allocatable, dimension(:,:,:) :: old
    complex (pr), allocatable, dimension(:,:,:) :: old2
  end type var

  type, public :: deriv
    complex (pr), allocatable, dimension(:,:,:) :: x
    complex (pr), allocatable, dimension(:,:,:) :: y
    complex (pr), allocatable, dimension(:,:,:) :: z
    complex (pr), allocatable, dimension(:,:,:) :: xx
    complex (pr), allocatable, dimension(:,:,:) :: yy
    complex (pr), allocatable, dimension(:,:,:) :: zz
  end type deriv

  type, public :: re_im
    real (pr), allocatable, dimension(:,:,:) :: re
    real (pr), allocatable, dimension(:,:,:) :: im
  end type re_im

  integer, dimension(-1:nyprocs, -1:nzprocs), private :: itable
  
  ! Constants for numerical integration
  real (pr), parameter, private :: c1 = 3.0_pr/8.0_pr, &
                                   c2 = 7.0_pr/6.0_pr, &
                                   c3 = 23.0_pr/24.0_pr

  contains

  function laplacian(in_var)
    ! Laplacian in cartesian coordinates
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: laplacian, dxx, dyy, dzz

    call deriv_xx(in_var,dxx)
    call deriv_yy(in_var,dyy)
    call deriv_zz(in_var,dzz)
    
    laplacian = dxx + dyy + dzz

    return
  end function laplacian

! ***************************************************************************  

  function grad(in_var)
    ! Gradient in cartesian coordinates
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: derx, dery, derz
    complex (pr), dimension(3,0:nx1,js:je,ks:ke) :: grad

    call deriv_x(in_var,derx)
    call deriv_y(in_var,dery)
    call deriv_z(in_var,derz)
    
    grad(1,:,:,:) = derx
    grad(2,:,:,:) = dery
    grad(3,:,:,:) = derz

    return
  end function grad

! ***************************************************************************  

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

! ***************************************************************************  

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

! ***************************************************************************  

  subroutine array_len()
    ! Determine the length of each array dimension on each process.  This 
    ! allows for the possibility that there are references to (for example) 
    ! j+1, k+1 simultaneously
    use parameters
    implicit none

    jlen = je-js+1
    kks = max(0, ks-1)
    kke = min(nz1, ke+1)
    klen = kke-kks+1

    return
  end subroutine array_len

! ***************************************************************************  

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

! ***************************************************************************  

  subroutine send_recv_z(in_var)
    ! Send and receive boundary elements of each process to and from
    ! neighbouring processes in the z-direction.  This data is contiguous in
    ! memory and so can be immediately be sent
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    integer, dimension(4) :: jsend, jrecv

    ! Send the data two away from the boundary (to allow for fourth order
    ! derivatives)
    call MPI_ISEND(in_var(0,js,ke-1), nx*jlen, gpe_mpi_complex, &
      znext, 1, MPI_COMM_WORLD, jsend(1), ierr)
    call MPI_ISEND(in_var(0,js,ks+1), nx*jlen, gpe_mpi_complex, &
      zprev, 1, MPI_COMM_WORLD, jsend(2), ierr)
                 
    ! Receive the data two away from the boundary
    call MPI_IRECV(in_var(0,js,ks-2), nx*jlen, gpe_mpi_complex, &
      zprev, 1, MPI_COMM_WORLD, jrecv(1), ierr)
    call MPI_IRECV(in_var(0,js,ke+2), nx*jlen, gpe_mpi_complex, &
      znext, 1, MPI_COMM_WORLD, jrecv(2), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(jsend(1), istatus, ierr)
    call MPI_WAIT(jsend(2), istatus, ierr)
    call MPI_WAIT(jrecv(1), istatus, ierr)
    call MPI_WAIT(jrecv(2), istatus, ierr)

    ! Send the data one away from the boundary
    call MPI_ISEND(in_var(0,js,ke), nx*jlen, gpe_mpi_complex, &
      znext, 1, MPI_COMM_WORLD, jsend(3), ierr)
    call MPI_ISEND(in_var(0,js,ks), nx*jlen, gpe_mpi_complex, &
      zprev, 1, MPI_COMM_WORLD, jsend(4), ierr)
                   
    ! Receive the data one away from the boundary
    call MPI_IRECV(in_var(0,js,ks-1), nx*jlen, gpe_mpi_complex, &
      zprev, 1, MPI_COMM_WORLD, jrecv(3), ierr)
    call MPI_IRECV(in_var(0,js,ke+1), nx*jlen, gpe_mpi_complex, &
      znext, 1, MPI_COMM_WORLD, jrecv(4), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(jsend(3), istatus, ierr)
    call MPI_WAIT(jsend(4), istatus, ierr)
    call MPI_WAIT(jrecv(3), istatus, ierr)
    call MPI_WAIT(jrecv(4), istatus, ierr)

    return
  end subroutine send_recv_z

! ***************************************************************************  

  subroutine send_recv_y()
    ! Send and receive boundary elements of each process to and from
    ! neighbouring processes in the y-direction.  This data is NOT contiguous
    ! in memory and so must first be packed into a contiguous array (see pack_y
    ! and unpack_y below)
    use parameters
    implicit none

    integer, dimension(2) :: ksend, krecv

    ! Send the boundary data to neighbouring processes
    call MPI_ISEND(works1(0,1,kks), nx*2*klen, gpe_mpi_complex, &
      ynext, 1, MPI_COMM_WORLD, ksend(1), ierr)
    call MPI_ISEND(works2(0,1,kks), nx*2*klen, gpe_mpi_complex, &
      yprev, 1, MPI_COMM_WORLD, ksend(2), ierr)

    ! Receive the boundary data from neighbouring processes
    call MPI_IRECV(workr1(0,1,kks), nx*2*klen, gpe_mpi_complex, &
      yprev, 1, MPI_COMM_WORLD, krecv(1), ierr)
    call MPI_IRECV(workr2(0,1,kks), nx*2*klen, gpe_mpi_complex, &
      ynext, 1, MPI_COMM_WORLD, krecv(2), ierr)

    ! Wait until all send/receive operations have been completed
    call MPI_WAIT(ksend(1), istatus, ierr)
    call MPI_WAIT(ksend(2), istatus, ierr)
    call MPI_WAIT(krecv(1), istatus, ierr)
    call MPI_WAIT(krecv(2), istatus, ierr)

    return
  end subroutine send_recv_y

! ***************************************************************************  

  subroutine pack_y(in_var)
    ! Pack the non-contiguous boundary data for the y-direction into a
    ! contiguous array
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    integer :: k

    if (myranky /= nyprocs-1) then
      do k=kks,kke
        works1(:,:,k) = in_var(:,je-1:je,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! pack their boundary data
      if (myranky == nyprocs-1) then
        do k=kks,kke
          works1(:,:,k) = in_var(:,je-1:je,k)
        end do
      end if
    end if

    if (myranky /= 0) then
      do k=kks,kke
        works2(:,:,k) = in_var(:,js:js+1,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! pack their boundary data
      if (myranky == 0) then
        do k=kks,kke
          works2(:,:,k) = in_var(:,js:js+1,k)
        end do
      end if
    end if

    return
  end subroutine pack_y

! ***************************************************************************  

  subroutine unpack_y(in_var)
    ! Unpack the contiguous boundary data back into the non-contiguous array
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(out) :: in_var
    integer :: i, k
    
    if (myranky /= 0) then
      do k=kks,kke
        in_var(:,js-2:js-1,k) = workr1(:,:,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! unpack their boundary data
      if (myranky == 0) then
        do k=kks,kke
          in_var(:,js-2:js-1,k) = workr1(:,:,k)
        end do
      end if
    end if

    if (myranky /= nyprocs-1) then
      do k=kks,kke
        in_var(:,je+1:je+2,k) = workr2(:,:,k)
      end do
    end if

    if (bcs == 1) then
      ! If using periodic BCs then make sure that the boundary processes also
      ! unpack their boundary data
      if (myranky == nyprocs-1) then
        do k=kks,kke
          in_var(:,je+1:je+2,k) = workr2(:,:,k)
        end do
      end if
    end if

    return
  end subroutine unpack_y

! ***************************************************************************  

  subroutine get_phase(in_var, phase)
    ! Calculate the phase
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var
    real (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: phase
    integer :: j, k

    phase = atan2(aimag(in_var)+1.0e-6_pr, real(in_var, pr))

    return
  end subroutine get_phase

! ***************************************************************************  

  subroutine get_pdf_velocity(in_var, vx, vy, vz, vmean, vstdev)
    ! Calculate the velocity
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    real (pr), allocatable, dimension(:) :: vx, vy, vz
    real (pr), dimension(3), intent(out) :: vmean, vstdev
    complex (pr), dimension(0:nx1,js:je,ks:ke) :: tmp_vx, tmp_vy, tmp_vz
    real (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: phase
    real (pr), dimension(3,0:nx1,js:je,ks:ke) :: vel, gradient, gradstar
    integer, dimension(3) :: valid_vel
    integer :: i

    gradient = grad(in_var)
    gradstar = grad(conjg(in_var))

    do i=1,3
      where (abs(in_var(:,js:je,ks:ke))**2 > 0.01_pr)
        vel(i,:,:,:) = -0.5_pr*eye * ( &
        conjg(in_var(:,js:je,ks:ke))*gradient(i,:,:,:) - &
        in_var(:,js:je,ks:ke)*gradstar(i,:,:,:) ) / &
        abs(in_var(:,js:je,ks:ke))**2
      elsewhere
        vel(i,:,:,:) = 1e3_pr
      end where
      valid_vel(i) = count(vel(i,:,:,:) < 0.5e3_pr)
      vmean(i) = mean(pack(vel(i,:,:,:), vel(i,:,:,:) < 0.5e3_pr))
      vstdev(i) = stdev(pack(vel(i,:,:,:), vel(i,:,:,:) < 0.5e3_pr), vmean(i))
    end do

    allocate(vx(valid_vel(1)))
    allocate(vy(valid_vel(2)))
    allocate(vz(valid_vel(3)))

    vx = pack(vel(1,:,:,:), vel(1,:,:,:) < 0.5e3_pr)
    vy = pack(vel(2,:,:,:), vel(2,:,:,:) < 0.5e3_pr)
    vz = pack(vel(3,:,:,:), vel(3,:,:,:) < 0.5e3_pr)

    return
  end subroutine get_pdf_velocity

! ***************************************************************************  
  
  subroutine get_pdf(in_var, pdf, max_vel)
    ! Calculate a PDF
    use parameters
    implicit none

    real (pr), dimension(:), intent(in) :: in_var
    real (pr), dimension(-nbins/2+1:nbins/2), intent(out) :: pdf
    real (pr), intent(out) :: max_vel
    integer, dimension(-nbins/2+1:nbins/2) :: hist, total_hist
    real (pr), dimension(2) :: maxs, maxr, mins, minr
    integer :: i, tmp_total, total

    hist = 0

    ! Find max/min on each process
    maxs(1) = maxval(in_var)
    maxs(2) = 0.0_pr
    mins(1) = minval(in_var)
    mins(2) = 0.0_pr

    ! Find max/min over whole array
    call MPI_ALLREDUCE(maxs, maxr, 1, gpe_mpi_2real, MPI_MAXLOC, &
      MPI_COMM_WORLD, ierr)
                       
    call MPI_ALLREDUCE(mins, minr, 1, gpe_mpi_2real, MPI_MINLOC, &
      MPI_COMM_WORLD, ierr)

    ! Maximum permissible value is maximum of the absolute values
    max_vel = max(abs(maxr(1)), abs(minr(1)))

    ! Find the total number of values satisfying the conditions on each process
    do i=-nbins/2+1,nbins/2
      hist(i) = count( &
        in_var > 2.0_pr*real(i-1, pr)*max_vel/real(nbins, pr) .and. &
        in_var <= 2.0_pr*real(i, pr)*max_vel/real(nbins, pr))
    end do

    ! Sum up the counts across all processes
    call MPI_ALLREDUCE(hist, total_hist, nbins, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    ! Overall total count
    total = sum(total_hist)

    ! Actual PDF
    pdf = total_hist/real(total, pr)

    return
  end subroutine get_pdf

! ***************************************************************************  

  subroutine get_vcf(in_var, f)
    ! Calculate the velocity correlation function.
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    real (pr), dimension(0:nx1), intent(out) :: f
    real (pr), dimension(0:nx1) :: Qxx, total_Qxx
    real (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2) :: phase
    real (pr), dimension(3,0:nx1,js:je,ks:ke) :: vel, gradient, gradstar
    integer, dimension(0:nx1) :: num
    integer :: i, r

    gradient = grad(in_var)
    gradstar = grad(conjg(in_var))

    do i=1,3
      vel(i,:,:,:) = -0.5_pr*eye * ( &
      conjg(in_var(:,js:je,ks:ke))*gradient(i,:,:,:) - &
      in_var(:,js:je,ks:ke)*gradstar(i,:,:,:) ) / &
      abs(in_var(:,js:je,ks:ke))**2
    end do

    Qxx = 0.0_pr
    num = 0
    do r=0,nx1
      do i=0,nx1
        if (i+r <= nx1) then
          num(r) = num(r)+1
          Qxx(r) = Qxx(r) + sum(vel(1,i,:,:) * vel(1,i+r,:,:))
        end if
      end do
    end do

    call MPI_ALLREDUCE(Qxx, total_Qxx, nx, gpe_mpi_real, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    Qxx = total_Qxx/real(num*ny*nz, pr)
    f = Qxx/Qxx(0)

    return
  end subroutine get_vcf

! ***************************************************************************  

  subroutine get_density(in_var, density)
    ! Calculate the density
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var
    real (pr), dimension(0:nx1,js:je,ks:ke), intent(out) :: density
    integer :: j, k

    density = abs(in_var)**2

    return
  end subroutine get_density

! ***************************************************************************  

  subroutine get_norm(in_var, norm)
    ! Calculate the norm
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var
    real (pr), intent(out) :: norm
    real (pr) :: int_z
    real (pr), dimension(0:nx1,js:je,ks:ke) :: int_var
    real (pr), dimension(js:je,ks:ke) :: int_x
    real (pr), dimension(ks:ke) :: int_y
    integer :: j, k
    
    int_var = abs(in_var)**2
    
    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    norm = int_z

    return
  end subroutine get_norm

! ***************************************************************************  

  subroutine energy(in_var, E)
    ! Calculate the energy
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    real (pr), intent(inout) :: E
    real (pr), dimension(0:nx1,js:je,ks:ke) :: int_var
    real (pr) :: tmp

    int_var = real( -laplacian(in_var) * conjg(in_var(:,js:je,ks:ke)) + &
      0.5_pr*abs(in_var(:,js:je,ks:ke))**4, pr )

    tmp = 0.0_pr

    tmp = sum(int_var)

    call MPI_ALLREDUCE(tmp, E, 1, gpe_mpi_real, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    E = E*dx*dy*dz
    
    return
  end subroutine energy
  
! ***************************************************************************  

  subroutine mass(in_var, M)
    ! Calculate the mass
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var
    real (pr), intent(inout) :: M
    real (pr) :: tmp
    real (pr), dimension(0:nx1,js:je,ks:ke) :: int_var

    int_var = abs(in_var)**2

    tmp = 0.0_pr
    tmp = sum(int_var)
          
    call MPI_ALLREDUCE(tmp, M, 1, gpe_mpi_real, MPI_SUM, MPI_COMM_WORLD, ierr)

    M = M*dx*dy*dz
                       
    return
  end subroutine mass
  
! ***************************************************************************  

  subroutine momentum(in_var, mom)
    ! Calculate the momentum
    use parameters
    use derivs
    implicit none

    complex (pr), dimension(0:nx1,js-2:je+2,ks-2:ke+2), intent(in) :: in_var
    real (pr), dimension(3), intent(out) :: mom
    real (pr) :: int_z
    real (pr), dimension(0:nx1,js:je,ks:ke) :: int_var
    real (pr), dimension(js:je,ks:ke) :: int_x
    real (pr), dimension(ks:ke) :: int_y
    type (deriv) :: dpsi, dpsistar
    integer :: j, k

    allocate(dpsi%x(0:nx1,js:je,ks:ke))
    allocate(dpsi%y(0:nx1,js:je,ks:ke))
    allocate(dpsi%z(0:nx1,js:je,ks:ke))
    allocate(dpsistar%x(0:nx1,js:je,ks:ke))
    allocate(dpsistar%y(0:nx1,js:je,ks:ke))
    allocate(dpsistar%z(0:nx1,js:je,ks:ke))
    
    call deriv_x(in_var, dpsi%x)
    call deriv_y(in_var, dpsi%y)
    call deriv_z(in_var, dpsi%z)
    call deriv_x(conjg(in_var), dpsistar%x)
    call deriv_y(conjg(in_var), dpsistar%y)
    call deriv_z(conjg(in_var), dpsistar%z)
    
    do k=ks,ke
      do j=js,je
        int_var(:,j,k) = (real(in_var(:,j,k), pr)-1.0_pr) * &
          aimag(dpsistar%x(:,j,k)) - &
          aimag(in_var(:,j,k))*real(dpsi%x(:,j,k), pr)
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    mom(1) = int_z
    
    do k=ks,ke
      do j=js,je
        int_var(:,j,k) = (real(in_var(:,j,k), pr)-1.0_pr) * &
          aimag(dpsistar%y(:,j,k)) - &
          aimag(in_var(:,j,k))*real(dpsi%y(:,j,k), pr)
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    mom(2) = int_z
    
    do k=ks,ke
      do j=js,je
        int_var(:,j,k) = (real(in_var(:,j,k), pr)-1.0_pr) * &
          aimag(dpsistar%z(:,j,k)) - &
          aimag(in_var(:,j,k))*real(dpsi%z(:,j,k), pr)
      end do
    end do

    call integrate_x(int_var, int_x)
    call integrate_y(int_x, int_y)
    call integrate_z(int_y, int_z)
    
    mom(3) = int_z

    deallocate(dpsi%x)
    deallocate(dpsi%y)
    deallocate(dpsi%z)
    deallocate(dpsistar%x)
    deallocate(dpsistar%y)
    deallocate(dpsistar%z)

    return
  end subroutine momentum

! ***************************************************************************  

  subroutine renormalise(var, norm)
    ! Renormalise when running in imaginary time
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(inout) :: var
    real (pr), intent(in) :: norm

    var = sqrt(nn) * var / sqrt(norm)

    return
  end subroutine renormalise

! ***************************************************************************  

  subroutine integrate_x(in_var, x_int)
    ! Integrate a (3D) variable in x.  The x-direction is not parallelised so
    ! this integration is straight forward
    use parameters
    implicit none

    real (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var
    real (pr), dimension(js:je,ks:ke), intent(out) :: x_int
    integer :: j, k

    do k=ks,ke
      do j=js,je
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
  
! ***************************************************************************  

  subroutine integrate_y(in_var, y_int)
    ! Integrate a (2D) variable in y
    use parameters
    implicit none

    real (pr), dimension(js:je,ks:ke), intent(in) :: in_var
    real (pr), dimension(ks:ke), intent(out) :: y_int
    real (pr), dimension(js:je,ks:ke) :: tmp_var
    real (pr), dimension(0:nz1) :: tmp, total
    integer :: j, k, irank

    ! Create a temporary variable on which to perform operations
    tmp_var = in_var
    
    ! Update the elements which must be multiplied by the integrating constants
    do j=js,je
      if (j==0) tmp_var(j,:) = c1*in_var(j,:)
      if (j==1) tmp_var(j,:) = c2*in_var(j,:)
      if (j==2) tmp_var(j,:) = c3*in_var(j,:)
      if (j==ny-3) tmp_var(j,:) = c3*in_var(j,:)
      if (j==ny-2) tmp_var(j,:) = c2*in_var(j,:)
      if (j==ny-1) tmp_var(j,:) = c1*in_var(j,:)
    end do
      
    ! Sum the variable on individual processes
    tmp = 0.0_pr
    do k=ks,ke
      tmp(k) = sum(tmp_var(:,k))
    end do

    ! Sum the variable over all processes and make sure each process has the
    ! result
    call MPI_ALLREDUCE(tmp, total, nz, gpe_mpi_real, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    ! Calculate the final integrated result
    do k=ks,ke
      y_int(k) = total(k) * dy
    end do

    return
  end subroutine integrate_y
  
! ***************************************************************************  

  subroutine integrate_z(in_var, z_int)
    ! Integrate a (1D) variable in z
    use parameters
    implicit none

    real (pr), dimension(ks:ke), intent(in) :: in_var
    real (pr), intent(out) :: z_int
    real (pr), dimension(ks:ke) :: tmp_var
    real (pr) :: tmp
    integer :: k

    ! Create a temporary variable
    tmp_var = in_var
    
    ! Update the elements which must be multiplied by the integrating constants
    do k=ks,ke
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
    if (myranky /= 0) tmp = 0.0_pr

    ! Calculate the sum over all processes and make sure each process has the
    ! result
    call MPI_ALLREDUCE(tmp, z_int, 1, gpe_mpi_real, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    ! Calculate the final result
    z_int = z_int*dz

    return
  end subroutine integrate_z

! ***************************************************************************  

  function mean(in_var)
    ! Calculate the mean of a vector of values.  Use pack in the function call
    ! if the data are not in a 1D vector.
    use parameters
    implicit none

    real (pr), dimension(:), intent(in) :: in_var
    real (pr) :: mean, tmp_sum, total_sum
    integer :: tmp_size, total_size

    tmp_sum = sum(in_var)
    tmp_size = size(in_var)

    call MPI_ALLREDUCE(tmp_sum, total_sum, 1, gpe_mpi_real, MPI_SUM, &
      MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(tmp_size, total_size, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    mean = total_sum / real(total_size, pr)

    return
  end function mean

! ***************************************************************************  

  function stdev(in_var, mean)
    ! Calculate the standard deviation of a vector of values.  Use pack in the
    ! function call if the data are not in a 1D vector.  Mean must already have
    ! been calculated.
    use parameters
    implicit none

    real (pr), dimension(:), intent(in) :: in_var
    real (pr), intent(in) :: mean
    real (pr) :: stdev, tmp_sumsq, total_sumsq
    integer :: tmp_size, total_size

    tmp_sumsq = sum(in_var**2)
    tmp_size = size(in_var)

    call MPI_ALLREDUCE(tmp_sumsq, total_sumsq, 1, gpe_mpi_real, &
      MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(tmp_size, total_size, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, ierr)

    stdev = sqrt( (total_sumsq / real(total_size, pr)) - mean**2 )

    return
  end function stdev

! ***************************************************************************  

  function linelength(time, psi)
    ! Calculate the total line length of vortices over the whole box
    ! individually on each process.  Serial algorithm written by Natalia
    ! Berloff in Fortran77.  Updated by me for Fortran 90
    use parameters
    implicit none
  
    real (pr) :: linelength
    complex (pr), intent(in) :: psi(0:nx1,js-2:je+2,ks-2:ke+2)
    real (pr), intent(in) :: time
    real (pr) :: xx, yy, h, a, den
    integer :: i, j, k, n, m, q, lv, lu
    real (pr) :: x(6), y(6), z(6), u(5), v(5), xu(6), xv(6), yu(6), yv(6)
  
  ! MUST have dx=dy=dz
  
    linelength = 0.0_pr
    do k=ks,ke
      if ((k==0) .or. (k==nz1)) cycle
      do j=js,je
        if ((j==0) .or. (j==ny1)) cycle
        iloop: do i=1,nx1-1
          m=1
          nloop: do n=1,6
  ! 1-front,2-back, 3-right, 4-left, 5-top, 6-bottom
            select case (n)
              case (1)
                u(1) = real(psi(i,j,k), pr)
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i+1,j,k), pr)
                v(2) = aimag(psi(i+1,j,k))
                u(4) = real(psi(i,j+1,k), pr)
                v(4) = aimag(psi(i,j+1,k))
                u(3) = real(psi(i+1,j+1,k), pr)
                v(3) = aimag(psi(i+1,j+1,k))
              case (2) 
                u(1) = real(psi(i,j,k+1), pr)
                v(1) = aimag(psi(i,j,k+1))
                u(2) = real(psi(i+1,j,k+1), pr)
                v(2) = aimag(psi(i+1,j,k+1))
                u(4) = real(psi(i,j+1,k+1), pr)
                v(4) = aimag(psi(i,j+1,k+1))
                u(3) = real(psi(i+1,j+1,k+1), pr)
                v(3) = aimag(psi(i+1,j+1,k+1))
              case(3)
                u(1) = real(psi(i+1,j,k), pr)
                v(1) = aimag(psi(i+1,j,k))
                u(2) = real(psi(i+1,j,k+1), pr)
                v(2) = aimag(psi(i+1,j,k+1))
                u(4) = real(psi(i+1,j+1,k), pr)
                v(4) = aimag(psi(i+1,j+1,k))
                u(3) = real(psi(i+1,j+1,k+1), pr)
                v(3) = aimag(psi(i+1,j+1,k+1))
              case (4)
                u(1) = real(psi(i,j,k), pr)
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i,j,k+1), pr)
                v(2) = aimag(psi(i,j,k+1))
                u(4) = real(psi(i,j+1,k), pr)
                v(4) = aimag(psi(i,j+1,k))
                u(3) = real(psi(i,j+1,k+1), pr)
                v(3) = aimag(psi(i,j+1,k+1))
              case (5)
                u(1) = real(psi(i,j+1,k), pr)
                v(1) = aimag(psi(i,j+1,k))
                u(2) = real(psi(i+1,j+1,k), pr)
                v(2) = aimag(psi(i+1,j+1,k))
                u(4) = real(psi(i,j+1,k+1), pr)
                v(4) = aimag(psi(i,j+1,k+1))
                u(3) = real(psi(i+1,j+1,k+1), pr)
                v(3) = aimag(psi(i+1,j+1,k+1))
              case (6)
                u(1) = real(psi(i,j,k), pr)
                v(1) = aimag(psi(i,j,k))
                u(2) = real(psi(i+1,j,k), pr)
                v(2) = aimag(psi(i+1,j,k))
                u(4) = real(psi(i,j,k+1), pr)
                v(4) = aimag(psi(i,j,k+1))
                u(3) = real(psi(i+1,j,k+1), pr)
                v(3) = aimag(psi(i+1,j,k+1))
            end select
            u(5) = u(1)
            v(5) = v(1)
  
  ! find zero line in u and v
  ! deal with planes first
            if ((u(1)==0.0_pr .and. u(2)==0.0_pr .and. u(3)==0.0_pr .and. &
                 u(4)==0.0_pr .and. &
                (v(1)*v(2)<0.0_pr .or. v(2)*v(3)<0.0_pr .or. &
                 v(3)*v(4)<0.0_pr .or. v(4)*v(1)<0.0_pr)) .or. &
                (v(1)==0.0_pr .and. v(2)==0.0_pr .and. v(3)==0.0_pr .and. &
                 v(4)==0.0_pr .and. &
                (u(1)*u(2)<0.0_pr .or. u(2)*u(3)<0.0_pr .or. &
                 u(3)*u(4)<0.0_pr .or. u(4)*u(1)<0.0_pr))) then
              linelength = linelength+0.5_pr
              !print*, 'cycle iloop'
              a = 1.0_pr
              cycle iloop
            end if
  ! deal with edge
            do q=1,4
              if(u(q)==0.0_pr .and. u(q+1)==0.0_pr .and. &
                 v(q)==0.0_pr .and. v(q+1)==0.0_pr) then
                linelength = linelength+0.25_pr
                !print*, 'cycle iloop'
                a = 1.0_pr
                cycle iloop
              end if
            end do
  
            lu=1
            do q=1,4
              if(u(q)==0.0_pr .and. u(q+1)==0.0_pr .and. &
                 v(q)==0.0_pr .and. v(q+1)==0.0_pr) then
                m=m+1
                print*, 'exit nloop'
                exit nloop
              else if (u(q)==0.0_pr) then
                select case (q)
                  case (1)
                    xu(lu)=0.0_pr
                    yu(lu)=0.0_pr
                  case (2)
                    xu(lu)=1.0_pr
                    yu(lu)=0.0_pr
                  case (3)
                    xu(lu)=1.0_pr
                    yu(lu)=1.0_pr
                  case (4)
                    xu(lu)=0.0_pr
                    yu(lu)=1.0_pr
                end select
              else if (u(q)*u(q+1)<0.0_pr) then
                select case (q)
                  case (1)
                    xu(lu) = abs(u(q)/(u(q)-u(q+1)))
                    yu(lu) = 0.0_pr
                  case (2)
                    xu(lu) = 1.0_pr
                    yu(lu) = abs(u(q)/(u(q)-u(q+1)))
                  case (3)
                    xu(lu) = abs(u(q+1)/(u(q)-u(q+1)))
                    yu(lu) = 1.0_pr
                  case (4)
                    xu(lu) = 0.0_pr
                    yu(lu) = abs(u(q+1)/(u(q)-u(q+1)))
                end select
  
                if (lu==1 .or. xu(lu)/=xu(lu-1) .or. yu(lu)/=yu(lu-1)) then
                  lu = lu+1
                end if
              end if
            end do
            lv = 1
            do q=1,4
              if (v(q)==0.0_pr) then
                select case (q)
                  case (1)
                    xv(lv) = 0.0_pr
                    yv(lv) = 0.0_pr
                  case (2)
                    xv(lv) = 1.0_pr
                    yv(lv) = 0.0_pr
                  case (3)
                    xv(lv) = 1.0_pr
                    yv(lv) = 1.0_pr
                  case (4)
                    xv(lv) = 0.0_pr
                    yv(lv) = 1.0_pr
                end select
              else if (v(q)*v(q+1)<0.0_pr) then
                select case (q)
                  case (1)
                    xv(lv) = abs(v(q)/(v(q)-v(q+1)))
                    yv(lv) = 0.0_pr
                  case (2)
                    xv(lv) = 1.0_pr
                    yv(lv) = abs(v(q)/(v(q)-v(q+1)))
                  case (3)
                    xv(lv) = abs(v(q+1)/(v(q)-v(q+1)))
                    yv(lv) = 1.0_pr
                  case (4)
                    xv(lv) = 0.0_pr
                    yv(lv) = abs(v(q+1)/(v(q)-v(q+1)))
                end select
  
                if (lv==1 .or. xv(lv)/=xv(lv-1) .or. yv(lv)/=yv(lv-1)) then
                  lv = lv+1
                end if
              end if
            end do ! in q
            if (lu>2 .and. lv>2) then
              den = xv(2)*(yu(1)-yu(2))+xv(1)*(yu(2)-yu(1))+ &
                   (xu(1)-xu(2))*(yv(1)-yv(2))
              if (den==0.0_pr)  then
                write (*,*) i, j, k, xu(1), yu(1), xu(2), yu(2), &
                            xv(1), yv(1), xv(2), yv(2)
                print*, 'ZERO DENOM IN linelength'
                den = den + 0.0000001_pr
              end if
              xx = (xu(1)*(xv(2)*(yv(1)-yu(2))+xv(1)*(yu(2)-yv(2)))+ &
                    xu(2)*(xv(2)*(yu(1)-yv(1))+xv(1)*(yv(2)-yu(1))))/den
              yy = (xv(2)*(yu(1)-yu(2))*yv(1)+xu(1)*yu(2)*yv(1)- &
                    xv(1)*yu(1)*yv(2)-xu(1)*yu(2)*yv(2)+xv(1)*yu(2)*yv(2)+ &
                    xu(2)*yu(1)*(yv(2)-yv(1)))/den
  
              if (xx>=0.0_pr .and. xx<=1.0_pr .and. yy>=0.0_pr .and. yy<=1.0_pr) then
  ! found zero inside the square
                select case (n)
                  case (1)
                    x(m) = xx
                    y(m) = yy
                    z(m) = 0.0_pr
                  case (2)
                    x(m) = xx
                    y(m) = yy
                    z(m) = 1.0_pr
                  case (3)
                    x(m) = 1.0_pr
                    y(m) = yy
                    z(m) = xx
                  case (4)
                    x(m) = 0.0_pr
                    y(m) = yy
                    z(m) = xx
                  case (5)
                    x(m) = xx
                    y(m) = 1.0_pr
                    z(m) = yy
                  case (6)
                    x(m) = xx
                    y(m) = 0.0_pr
                    z(m) = yy
                end select
  
                !if (k==41 .and. i==41 .and. j==41) then
                !  print*, 'together=', n, m, den, x(m), y(m), z(m)
                !end if
                m=m+1
              end if
            end if
          end do nloop
          if (m>1) then      ! found zero in at least two sides
            a = 1.0_pr !scale
            if (x(1)==x(2) .and. (x(1)==0.0_pr .or. x(1)==1.0_pr) .or. &
                y(1)==y(2) .and. (y(1)==0.0_pr .or. y(1)==1.0_pr) .or. &
                y(1)==y(2) .and. (y(1)==0.0_pr .or. y(1)==1.0_pr)) then
              a = a*0.5_pr
              if ((x(1)==x(2) .and. (x(1)==0.0_pr .or. x(1)==1.0_pr) .and. &
                  (y(1)==y(2) .and. (y(1)==0.0_pr .or. y(1)==1.0_pr) .or. &
                   y(1)==y(2) .and. (y(1)==0.0_pr .or. y(1)==1.0_pr))) .or. &
                  (y(1)==y(2) .and. (y(1)==0.0_pr .or. y(1)==1.0_pr) .and. &
                  (z(1)==z(2) .and. (z(1)==0.0_pr .or. z(1)==1.0_pr)))) then
                a=a*0.5_pr
              end if
            end if
            linelength = linelength+a*sqrt((x(1)-x(2))**2+(y(1)-y(2))**2+ &
                                           (z(1)-z(2))**2)
          end if
          a = 1.0_pr
        end do iloop
      end do
    end do
    linelength = linelength*dx

    return
  end function linelength

! ***************************************************************************

  function imprint_vortex_line(in_var)
    ! Imprint a vortex line on the wave function.  This is used mainly for
    ! solving case 4 when in imaginary time, to constantly imprint the phase to
    ! construct a vortex line.
    use ic, only : vortex_line
    use parameters
    implicit none

    complex (pr), dimension(0:nx1,js:je,ks:ke) :: imprint_vortex_line
    complex (pr), dimension(0:nx1,js:je,ks:ke), intent(in) :: in_var

    imprint_vortex_line = abs(in_var) * vortex_line(vl1)

    return
  end function imprint_vortex_line

end module variables
