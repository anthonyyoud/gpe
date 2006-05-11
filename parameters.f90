module parameters
  ! Parameters to set
  implicit none
  save

  include 'mpif.h'

  integer,      parameter :: nyprocs      = 2
  integer,      parameter :: nzprocs      = 2
  integer,      parameter :: nx           = 64
  integer,      parameter :: ny           = 64
  integer,      parameter :: nz           = 128
  complex                 :: time_step    = (0.0,-0.0001)
  real,         parameter :: end_time     = 0.5
  real,         parameter :: xr           = 32.0
  real,         parameter :: yr           = 32.0
  real,         parameter :: zr           = 64.0
  ! bcs = 1 for periodic, 2 for reflective
  integer,      parameter :: bcs          = 2
  ! order = 2 for 2nd order derivatives, 4 for 4th order derivatives
  integer,      parameter :: order        = 4
  integer,      parameter :: save_rate    = 1
  real,         parameter :: save_rate2   = 1.0
  logical,      parameter :: save_contour = .true.
  logical,      parameter :: save_3d      = .true.
  logical,      parameter :: save_zeros   = .false.
  logical,      parameter :: restart      = .false.
  logical                 :: real_time    = .true.
  character(*), parameter :: scheme       = 'rk_adaptive'

  ! Parameters for adaptive time stepping
  real, parameter :: eps              = 1e-6
  real, parameter :: safety           = 0.9
  real, parameter :: dt_decrease      = -0.25
  real, parameter :: dt_increase      = -0.20
  real            :: errcon

  ! Vortex line parameters **************************************************
  !
  type :: line_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: amp         ! amplitude of a disturbance of the vortex line
    real :: ll          ! wavelength of the above disturbance
    real :: sgn         ! sign of the argument of the line
  end type line_param

  type (line_param), parameter :: vl1 = line_param( 0.0, 3.0, 0.1,33.0, 1.0)
  type (line_param), parameter :: vl2 = line_param(-3.0, 3.0,-0.0,33.0,-1.0)
  type (line_param), parameter :: vl3 = line_param( 0.0,-3.0,-0.1,33.0,-1.0)
  type (line_param), parameter :: vl4 = line_param(-3.0,-3.0,-0.0,33.0, 1.0)
  !  
  ! *************************************************************************

  ! Vortex ring parameters **************************************************
  !
  type :: ring_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: r0          ! radius
  end type ring_param

  type (ring_param), parameter :: vr1 = ring_param(10.0,17.0,10.0)
  !type (ring_param), parameter :: vr1 = ring_param(5.0,0.0,0.0)
  type (ring_param), parameter :: vr2 = ring_param(-10.0,0.0,7.0)
  !type (ring_param), parameter :: vr2 = ring_param(-5.0,0.0,0.0)
  !
  ! *************************************************************************

  ! Parameters that don't need changing
  integer, parameter :: nx1         = nx-1
  integer, parameter :: ny1         = ny-1
  integer, parameter :: nz1         = nz-1
  integer, parameter :: nprocs      = nyprocs*nzprocs
  integer            :: end_proc    = 0
  integer            :: myrank, myranky, myrankz
  integer            :: jsta, jend, jlen, yprev, ynext
  integer            :: ksta, kksta, kend, kkend, klen, zprev, znext
  integer            :: ierr
  
  integer,              dimension(MPI_STATUS_SIZE) :: istatus
  integer,              dimension(0:nzprocs-1)     :: jdisp, kklen
  integer,              dimension(0:nyprocs-1)     :: kdisp, jjlen
  complex, allocatable, dimension(:,:,:)           :: works1, works2, &
                                                      workr1, workr2
  character(7)       :: proc_dir = 'proc**/'
  logical            :: first_write = .true.
  complex            :: dt
  real               :: t           = 0
  real               :: im_t        = 0
  real,    parameter :: pi          = 3.14159265358979
  real,    parameter :: xl          = -xr
  real,    parameter :: yl          = -yr
  real,    parameter :: zl          = -zr
  real,    parameter :: dx          = (xr-xl)/nx
  real,    parameter :: dy          = (yr-yl)/ny
  real,    parameter :: dz          = (zr-zl)/nz
  real,    parameter :: dx2         = dx**2
  real,    parameter :: dy2         = dy**2
  real,    parameter :: dz2         = dz**2
  complex, parameter :: eye         = (0.0,1.0)
  logical            :: switched    = .false.
  integer            :: p
  
end module parameters
