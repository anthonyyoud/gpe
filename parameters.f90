! $Id$
!----------------------------------------------------------------------------

module parameters
  ! Parameters amd global variables.  Nothing here needs to be changed.  Edit
  ! run.in and parameters.in instead.
  implicit none
  save

  include 'mpif.h'
  include 'parameters.in'

  logical :: pp_filtered_surface
  integer :: nlines
  integer :: nfilter
  real (pr) :: fscale

  real (pr) :: tau
  real (pr) :: end_time
  real (pr) :: xr
  real (pr) :: yr
  real (pr) :: zr
  real (pr) :: Urhs
  real (pr) :: diss_amp
  real (pr) :: scal
  real (pr) :: nv
  real (pr) :: enerv
  real (pr) :: g
  real (pr) :: mu
  real (pr) :: nn
  real (pr) :: omx, omy, omz
  integer :: eqn_to_solve
  integer :: bcs
  integer :: order
  integer :: nbins
  integer :: save_rate
  real (pr) :: save_rate2
  real (pr) :: save_rate3
  real (pr) :: p_save
  logical  :: save_contour
  logical  :: save_3d
  logical  :: save_filter
  logical  :: save_average
  logical  :: save_spectrum
  logical  :: save_pdf
  logical  :: save_vcf
  logical  :: save_ll
  logical  :: save_zeros
  logical  :: restart
  logical  :: saved_restart
  logical  :: renorm
  logical  :: imprint_vl
  logical  :: stop_imag
  logical :: real_time
  logical :: diagnostic
  character(5) :: scheme

  namelist /run_params/ tau, end_time, xr, yr, zr, scheme, eqn_to_solve, &
    bcs, order, restart, saved_restart, renorm, imprint_vl, stop_imag, &
    real_time

  namelist /eqn_params/ Urhs, diss_amp, scal, nv, enerv, g, mu, nn, &
    omx, omy, omz

  namelist /io_params/ save_rate, save_rate2, save_rate3, p_save, &
    save_contour, save_3d, save_filter, save_average, save_spectrum, &
    save_pdf, save_vcf, save_ll, save_zeros

  namelist /misc_params/ nbins, diagnostic, pp_filtered_surface, nlines, &
    nfilter, fscale

  ! Parameters for adaptive time stepping
  real (pr), parameter :: eps              = 1e-6_pr
  real (pr), parameter :: safety           = 0.9_pr
  real (pr), parameter :: dt_decrease      = -0.25_pr
  real (pr), parameter :: dt_increase      = -0.20_pr
  real (pr)            :: errcon

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
  real (pr) :: xl
  real (pr) :: yl
  real (pr) :: zl
  real (pr) :: dx
  real (pr) :: dy
  real (pr) :: dz
  real (pr) :: dx2
  real (pr) :: dy2
  real (pr) :: dz2
  complex (pr), parameter :: eye = (0.0_pr,1.0_pr)
  integer :: p
  integer :: snapshots   = 1
  integer :: gpe_mpi_real
  integer :: gpe_mpi_2real
  integer :: gpe_mpi_complex
  integer :: gpe_mpi_2complex

end module parameters
