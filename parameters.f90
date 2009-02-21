! $Id: parameters.f90,v 1.64 2009-02-21 14:10:19 youd Exp $
!----------------------------------------------------------------------------

module parameters
  ! Parameters to set
  implicit none
  save

  include 'mpif.h'

  logical, parameter :: pp_filtered_surface = .false.
  integer, parameter :: nlines              = 41
  integer, parameter :: nfilter             = 1
  real,    parameter :: fscale              = 1.0

  integer,      parameter :: nyprocs      = 4
  integer,      parameter :: nzprocs      = 4
  integer,      parameter :: nx           = 256
  integer,      parameter :: ny           = 256
  integer,      parameter :: nz           = 256
  real,         parameter :: tau          = 0.001
  real,         parameter :: end_time     = 500.0
  real,         parameter :: xr           = 128.0
  real,         parameter :: yr           = 128.0
  real,         parameter :: zr           = 128.0
  real,         parameter :: Urhs         = 0.0 !0.35
  real,         parameter :: diss_amp     = 0.0 !0.005
  real,         parameter :: scal         = 1.0 !0.64315009229562
  real,         parameter :: nv           = 0.5
  real,         parameter :: enerv        = 0.75
  ! see bottom of solve.f90 for possible values
  integer,      parameter :: eqn_to_solve = 1
  ! bcs = 1 for periodic, 2 for reflective
  integer,      parameter :: bcs          = 2
  ! order = 2 for 2nd order derivatives, 4 for 4th order derivatives
  integer,      parameter :: order        = 4
  integer,      parameter :: save_rate    = 50
  real,         parameter :: save_rate2   = 10.0
  real,         parameter :: save_rate3   = 10.0
  real,         parameter :: p_save       = 10.0
  logical,      parameter :: save_contour = .true.
  logical,      parameter :: save_3d      = .true.
  logical,      parameter :: save_filter  = .false.
  logical,      parameter :: save_average = .false.
  logical,      parameter :: save_spectrum= .false.
  logical,      parameter :: save_ll      = .true.
  logical,      parameter :: save_zeros   = .false.
  logical,      parameter :: restart      = .false.
  logical,      parameter :: saved_restart= .false.
  logical                 :: real_time    = .true.
  logical                 :: diagnostic   = .false.
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
    real         :: x0  ! x position
    real         :: y0  ! y position
    real         :: z0  ! z position
    real         :: amp ! amplitude of a disturbance of the vortex line
    real         :: ll  ! wavelength of the above disturbance
    real         :: sgn ! sign of the argument of the line
    character(1) :: dir ! direction in which the line should extend
  end type line_param

  !type (line_param), parameter :: &
  !  vl1 = line_param(   0.0,   0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl2 = line_param(   0.0,-128.0, 0.0, 0.0, 33.0,-1.0, 'z')
  !type (line_param), parameter :: &
  !  vl3 = line_param(   0.0, 128.0, 0.0, 0.0, 33.0,-1.0, 'z')
  !type (line_param), parameter :: &
  !  vl4 = line_param(-128.0,   0.0, 0.0, 0.0, 33.0,-1.0, 'z')
  !type (line_param), parameter :: &
  !  vl5 = line_param(-128.0,-128.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl6 = line_param(-128.0, 128.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl7 = line_param( 128.0,   0.0, 0.0, 0.0, 33.0,-1.0, 'z')
  !type (line_param), parameter :: &
  !  vl8 = line_param( 128.0,-128.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl9 = line_param( 128.0, 128.0, 0.0, 0.0, 33.0, 1.0, 'z')

  type (line_param), parameter :: &
    vl1 = line_param(-4.0, 0.0, 0.0, 0.1, 33.0, 1.0, 'y')
  type (line_param), parameter :: &
    vl2 = line_param( 4.0, 0.0, 0.0, 0.1, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl3 = line_param( 4.0, 4.0, 0.0, 0.1, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl4 = line_param( 4.0,-4.0, 0.0, 0.1, 33.0, 1.0, 'z')

  ! 2 bundles of 5 vortex lines
  !type (line_param), parameter :: &
  !  vl1 = line_param( -4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl2 = line_param(-12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl3 = line_param(-20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl4 = line_param(-12.0, 0.0, 8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl5 = line_param(-12.0, 0.0,-8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl6 = line_param(  4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl7 = line_param( 12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl8 = line_param( 20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl9 = line_param( 12.0, 8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl10 = line_param( 12.0,-8.0, 0.0, 0.0, 33.0, 1.0, 'z')

  ! 2 bundles of 7 vortex lines
  !type (line_param), parameter :: &
  !  vl1 = line_param( -4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl2 = line_param(-12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl3 = line_param(-20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl4 = line_param( -8.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl5 = line_param(-16.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl6 = line_param( -8.0, 0.0,-6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl7 = line_param(-16.0, 0.0,-6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl8 = line_param(  4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl9 = line_param( 12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl10 = line_param( 20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl11 = line_param(  8.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl12 = line_param( 16.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl13 = line_param(  8.0,-6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl14 = line_param( 16.0,-6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')

  ! 2 bundles of 9 vortex lines
  !type (line_param), parameter :: &
  !  vl1 = line_param( -4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl2 = line_param(-12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl3 = line_param(-20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl4 = line_param( -4.0, 0.0, 8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl5 = line_param(-12.0, 0.0, 8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl6 = line_param(-20.0, 0.0, 8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl7 = line_param( -4.0, 0.0,-8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl8 = line_param(-12.0, 0.0,-8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl9 = line_param(-20.0, 0.0,-8.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl10 = line_param(  4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl11 = line_param( 12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl12 = line_param( 20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl13 = line_param(  4.0, 8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl14 = line_param( 12.0, 8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl15 = line_param( 20.0, 8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl16 = line_param(  4.0,-8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl17 = line_param( 12.0,-8.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl18 = line_param( 20.0,-8.0, 0.0, 0.0, 33.0, 1.0, 'z')

  ! 2 bundles of 19 vortex lines
  !type (line_param), parameter :: &
  !  vl1 = line_param( -4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl2 = line_param(-12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl3 = line_param(-20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl4 = line_param(-28.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl5 = line_param(-36.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl6 = line_param( -8.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl7 = line_param(-16.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl8 = line_param(-24.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl9 = line_param(-32.0, 0.0, 6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl10 = line_param( -8.0, 0.0, -6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl11 = line_param(-16.0, 0.0, -6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl12 = line_param(-24.0, 0.0, -6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl13 = line_param(-32.0, 0.0, -6.92820323, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl14 = line_param(-12.0, 0.0, 13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl15 = line_param(-20.0, 0.0, 13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl16 = line_param(-28.0, 0.0, 13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl17 = line_param(-12.0, 0.0, -13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl18 = line_param(-20.0, 0.0, -13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl19 = line_param(-28.0, 0.0, -13.85640646, 0.0, 33.0, 1.0, 'y')
  !type (line_param), parameter :: &
  !  vl20 = line_param( 4.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl21 = line_param(12.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl22 = line_param(20.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl23 = line_param(28.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl24 = line_param(36.0, 0.0, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl25 = line_param( 8.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl26 = line_param(16.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl27 = line_param(24.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl28 = line_param(32.0, 6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl29 = line_param( 8.0, -6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl30 = line_param(16.0, -6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl31 = line_param(24.0, -6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl32 = line_param(32.0, -6.92820323, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl33 = line_param(12.0, 13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl34 = line_param(20.0, 13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl35 = line_param(28.0, 13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl36 = line_param(12.0, -13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl37 = line_param(20.0, -13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !type (line_param), parameter :: &
  !  vl38 = line_param(28.0, -13.85640646, 0.0, 0.0, 33.0, 1.0, 'z')
  !  
  ! *************************************************************************

  ! Vortex ring parameters **************************************************
  !
  type :: ring_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: z0          ! z position
    real :: r0          ! radius
    real :: dir         ! Propagation direction (+/-1)
  end type ring_param

  type (ring_param), parameter :: vr1 = ring_param(0.0, 0.0, 0.0, 16.0, -1.0)
  type (ring_param), parameter :: vr2 = ring_param(-60.0, 0.0, 0.0, 16.0, -1.0)
  type (ring_param), parameter :: vr3 = ring_param(60.0, 0.0, 0.0, 16.0, -1.0)
  type (ring_param), parameter :: vr4 = ring_param(0.0, -60.0, 0.0, 16.0, -1.0)
  type (ring_param), parameter :: vr5 = ring_param(0.0, 60.0, 0.0, 16.0, -1.0)
  type (ring_param), parameter :: vr6 = ring_param(0.0, 0.0, -60.0, 16.0, -1.0)
  type (ring_param), parameter :: vr7 = ring_param(0.0, 0.0, 60.0, 16.0, -1.0)
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
  real,    allocatable, dimension(:,:,:)           :: ave
  character(7)       :: proc_dir = 'proc**/'
  character(15)      :: end_state_file = 'end_state**.dat'
  character(24)      :: filt_end_state_file = 'end_state_filtered**.dat'
  logical            :: first_write = .true.
  complex            :: dt
  real               :: t           = 0.0
  real               :: im_t        = 0.0
  real               :: kc2         = 0.0
  real               :: comp_amp    = 0.0
  real, dimension(3) :: maxvar      = 0.0
  real, dimension(3) :: minvar      = 0.0
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
  integer            :: snapshots   = 1
  
end module parameters
