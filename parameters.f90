! $Id$
!----------------------------------------------------------------------------

module parameters
  ! Parameters to set
  implicit none
  save

  include 'mpif.h'

  ! Parameterised real kind
  !integer, parameter :: pr = selected_real_kind(6,37)
  integer, parameter :: pr = selected_real_kind(15,307)

  logical, parameter :: pp_filtered_surface = .false.
  integer, parameter :: nlines              = 3
  integer, parameter :: nfilter             = 1
  real (pr),    parameter :: fscale         = 1.0_pr

  integer,      parameter :: nyprocs       = 2
  integer,      parameter :: nzprocs       = 2
  integer,      parameter :: nx            = 120
  integer,      parameter :: ny            = 120
  integer,      parameter :: nz            = 120
  real (pr),    parameter :: tau           = 0.0001_pr
  real (pr),    parameter :: end_time      = 10.0_pr
  real (pr),    parameter :: xr            = 12.0_pr
  real (pr),    parameter :: yr            = 12.0_pr
  real (pr),    parameter :: zr            = 12.0_pr
  real (pr),    parameter :: Urhs          = 0.0_pr !0.35_pr
  real (pr),    parameter :: diss_amp      = 0.0_pr !0.005_pr
  real (pr),    parameter :: scal          = 1.0_pr !0.64315009229562_pr
  real (pr),    parameter :: nv            = 0.75_pr
  real (pr),    parameter :: enerv         = 0.75_pr
  real (pr),    parameter :: g             = 5.0_pr
  real (pr),    parameter :: mu            = 8.0_pr !83.28581966_pr
  real (pr),    parameter :: nn            = 1.0_pr !1e5_pr
  real (pr), dimension(3),  parameter :: omega = (/1.0_pr, 1.0_pr, 2.0_pr/)
  ! see bottom of solve.f90 for possible values
  integer,      parameter :: eqn_to_solve  = 4
  ! bcs = 1 for periodic, 2 for reflective
  integer,      parameter :: bcs           = 2
  ! order = 2 for 2nd order derivatives, 4 for 4th order derivatives
  integer,      parameter :: order         = 2
  integer,      parameter :: nbins         = 256
  integer,      parameter :: save_rate     = 10
  real (pr),    parameter :: save_rate2    = 0.1_pr
  real (pr),    parameter :: save_rate3    = 10.0_pr
  real (pr),    parameter :: p_save        = 10.0_pr
  logical,      parameter :: save_contour  = .false.
  logical,      parameter :: save_3d       = .true.
  logical,      parameter :: save_filter   = .false.
  logical,      parameter :: save_average  = .false.
  logical,      parameter :: save_spectrum = .false.
  logical,      parameter :: save_pdf      = .false.
  logical,      parameter :: save_vcf      = .false.
  logical,      parameter :: save_ll       = .false.
  logical,      parameter :: save_zeros    = .false.
  logical,      parameter :: restart       = .false.
  logical,      parameter :: saved_restart = .false.
  logical,      parameter :: renorm        = .false.
  logical,      parameter :: stop_imag     = .true.
  logical                 :: real_time     = .false.
  logical                 :: diagnostic    = .false.
  character(*), parameter :: scheme        = 'rk2'

  ! Parameters for adaptive time stepping
  real (pr), parameter :: eps              = 1e-6_pr
  real (pr), parameter :: safety           = 0.9_pr
  real (pr), parameter :: dt_decrease      = -0.25_pr
  real (pr), parameter :: dt_increase      = -0.20_pr
  real (pr)            :: errcon

  ! Vortex line parameters **************************************************
  !
  type :: line_param
    real (pr)    :: x0   ! x position
    real (pr)    :: y0   ! y position
    real (pr)    :: z0   ! z position
    real (pr)    :: amp1 ! amplitude of a disturbance of the vortex line (dir1)
    real (pr)    :: amp2 ! amplitude of a disturbance of the vortex line (dir2)
    real (pr)    :: ll   ! wavelength of the above disturbance
    real (pr)    :: sgn  ! sign of the argument of the line
    character(1) :: dir  ! direction in which the line should extend
    logical      :: imprint_phase ! imprint phase only (no vortex core)
  end type line_param

  type (line_param), parameter :: &
    vl1 = line_param(-1.5_pr, 0.0_pr, 0.0_pr, &
      0.0_pr, 0.0_pr, 0.0_pr, 1.0_pr, 'z', .true.)
  !  
  ! *************************************************************************

  ! Vortex ring parameters **************************************************
  !
  type :: ring_param
    real (pr) :: x0      ! x position
    real (pr) :: y0      ! y position
    real (pr) :: z0      ! z position
    real (pr) :: r0      ! radius
    real (pr) :: amp     ! Amplitude of a planar disturbance
    integer   :: mm      ! Wavenumber of planar disturbance
    real (pr) :: r1      ! Radius of helical disturbance
    integer   :: kk      ! Wavenumber of helical disturbance
    real (pr) :: dir     ! Propagation in x-direction (+/-1)
  end type ring_param

  type (ring_param), parameter :: &
    vr1 = ring_param(-48.0_pr, 0.0_pr, 0.0_pr, 25.0_pr, &
      0.0_pr, 5, 0.0_pr, 10, -1.0_pr)
  !
  ! *************************************************************************

  ! Parameters that don't need changing
  integer, parameter :: nx1 = nx-1
  integer, parameter :: ny1 = ny-1
  integer, parameter :: nz1 = nz-1
  integer :: nprocs
  integer :: end_proc = 0
  integer :: myrank, myranky, myrankz
  integer :: js, je, jlen, yprev, ynext
  integer :: ks, kks, ke, kke, klen, zprev, znext
  integer :: ierr
  
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  integer, dimension(0:nzprocs-1) :: jdisp, kklen
  integer, dimension(0:nyprocs-1) :: kdisp, jjlen
  complex (pr), allocatable, dimension(:,:,:) :: works1, works2, workr1, workr2
  real (pr), allocatable, dimension(:,:,:) :: ave
  character(7) :: proc_dir = 'proc**/'
  character(15) :: end_state_file = 'end_state**.dat'
  character(24) :: filt_end_state_file = 'end_state_filtered**.dat'
  complex (pr) :: dt
  real (pr) :: t = 0.0_pr
  real (pr) :: im_t = 0.0_pr
  real (pr) :: kc2 = 0.0_pr
  real (pr) :: comp_amp = 0.0_pr
  real (pr), dimension(4) :: maxvar = 0.0_pr
  real (pr), dimension(4) :: minvar = 0.0_pr
  real (pr), parameter :: pi = 3.1415926535897932384626433832795_pr
  real (pr), parameter :: xl = -xr
  real (pr), parameter :: yl = -yr
  real (pr), parameter :: zl = -zr
  real (pr), parameter :: dx = (xr-xl)/nx1
  real (pr), parameter :: dy = (yr-yl)/ny1
  real (pr), parameter :: dz = (zr-zl)/nz1
  real (pr), parameter :: dx2 = dx**2
  real (pr), parameter :: dy2 = dy**2
  real (pr), parameter :: dz2 = dz**2
  complex (pr), parameter :: eye = (0.0_pr,1.0_pr)
  integer :: p
  integer :: snapshots   = 1
  integer :: gpe_mpi_real
  integer :: gpe_mpi_2real
  integer :: gpe_mpi_complex
  integer :: gpe_mpi_2complex
  
end module parameters
