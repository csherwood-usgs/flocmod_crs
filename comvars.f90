MODULE comvars
  IMPLICIT NONE
  SAVE
! l_testcase - if .TRUE. sets G(t) to values from lab experiment
  logical, parameter :: l_testcase = .TRUE.
!  integer, parameter :: rsh = SELECTED_REAL_KIND(p=6,r=37)
  integer, parameter :: rsh = SELECTED_REAL_KIND(p=13,r=200)

! not sure all of these can be set to one
  integer, parameter :: limin =1 
  integer, parameter :: liminp2 =1 
  integer, parameter :: limax =1 
  integer, parameter :: limaxm1 =1 
  integer, parameter :: ljmin =1
  integer, parameter :: ljmax =1
  integer, parameter :: ljminp2 =1
  integer, parameter :: ljmaxm1 =1
  integer, parameter :: kmax =1

  integer, parameter :: imin =1 
  integer, parameter :: imax =1 
  integer, parameter :: jmin =1
  integer, parameter :: jmax =1

  integer, parameter :: nv_mud = 15  !number of mud (floc) classes
  integer, parameter :: nv_tot = 15  !total number of particle classes
  integer, parameter :: imud1 = 1    !index to first mud (floc) class
  integer, parameter :: nvpc = 15    !total number of particle classes

! lchain - length of characters in filename
! rsh = size of REAL in kind() statement and _rsh for CPP
! rlg = not used
! pi = 3.14159...
! epsilon = truncation error used in absolute difference reference
! grav = 9.81 (m2/s) grav. accel 
  real(kind=rsh), parameter :: pi = 3.14159_rsh
  real(kind=rsh), parameter :: grav = 9.81_rsh
  real(kind=rsh), parameter :: rhoref = 1030._rsh
  real(kind=rsh), parameter :: epsilon = 1e-8
  real(kind=rsh) :: eps ! this is initialized as a tiny number

! Input:
! diam_sed - sediment diameters (m)
! ros - sediment density (kg/m3)
! f_ws - sediment settling velocity (m/s) (Stokes estimate)
  real(kind=rsh),dimension(nv_mud) :: diam_sed
  real(kind=rsh),dimension(nv_mud) :: ros
  real(kind=rsh),dimension(nv_mud) :: f_ws

  real(kind=rsh),dimension(nv_mud,kmax,imax,jmax) :: cv_wat ! input concentration array

  real(kind=rsh) :: dt
  real(kind=rsh) :: t

  LOGICAL,PUBLIC        :: l_ASH,l_ADS,l_COLLFRAG
  INTEGER,PUBLIC        :: f_ero_iv
  REAL(KIND=rsh),PUBLIC :: f_nf, f_dp0,f_alpha,f_beta,f_nb_frag,f_fter,f_dmax,f_ater, &
       f_ero_frac,f_ero_nbfrag,f_mneg_param,f_dt,f_collfragparam

  !ALLOCATE(f_diam(1:nv_mud))     ! floc diameter
  !ALLOCATE(f_vol(1:nv_mud))      ! floc volume
  !ALLOCATE(f_rho(1:nv_mud))      ! floc density
  !ALLOCATE(f_mass(0:nv_mud+1))     ! floc mass
  !ALLOCATE(f_cv(1:nv_mud))       ! extracted from cv_wat(:,k,j,j) mass concentration for every mud variables

  ! agregation kernels
  !ALLOCATE(f_coll_prob_sh(1:nv_mud,1:nv_mud)) !  shear agregation collision probability
  !ALLOCATE(f_coll_prob_ds(1:nv_mud,1:nv_mud)) ! differential settling collision probability
  !ALLOCATE(f_g1_sh(1:nv_mud,1:nv_mud,1:nv_mud)) ! shear agregation gain term
  !ALLOCATE(f_g1_ds(1:nv_mud,1:nv_mud,1:nv_mud)) ! differential settling agregation gain term
  !ALLOCATE(f_l1_sh(1:nv_mud,1:nv_mud)) ! shear agregation loss term
  !ALLOCATE(f_l1_ds(1:nv_mud,1:nv_mud)) ! differential settling agregation loss term  
  !ALLOCATE(f_g3(1:nv_mud,1:nv_mud)) ! fragmentation gain term     
  !ALLOCATE(f_l3(1:nv_mud)) ! fragmentation loss term
  !ALLOCATE(f_g4(1:nv_mud,1:nv_mud,1:nv_mud)) ! Collision fragmentation gain term  
  !ALLOCATE(f_l4(1:nv_mud,1:nv_mud)) ! Collision fragmentation loss term  

  REAL(KIND=rsh),DIMENSION(nv_mud),PUBLIC          :: f_diam,f_vol,f_rho,f_cv,f_l3
  REAL(KIND=rsh),DIMENSION(0:nv_mud+1),PUBLIC      :: f_mass
  REAL(KIND=rsh),DIMENSION(nv_mud,nv_mud),PUBLIC   :: f_coll_prob_sh, f_coll_prob_ds
  REAL(KIND=rsh),DIMENSION(nv_mud,nv_mud),PUBLIC   :: f_l1_sh, f_l1_ds, f_g3, f_l4
  REAL(KIND=rsh),DIMENSION(nv_mud,nv_mud,nv_mud),PUBLIC  :: f_g1_sh,f_g1_ds,f_g4
  REAL(KIND=rsh),DIMENSION(kmax,imax,jmax),PUBLIC :: f_d50,f_d90,f_d10,f_davg,f_gval,f_dtmin

END MODULE
