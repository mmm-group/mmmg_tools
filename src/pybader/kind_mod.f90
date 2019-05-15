MODULE kind_mod

!----------------------------------------------------------------------!
!                                                                      !
! kind_mod:  module for storing public parameters and types.           !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! contents:                                                            !
!                                                                      !
!   26 :: q1                                                           !
!   27 :: q2                                                           ! 
!   28 :: pi                                                           !
!   31 :: bader_obj                                                    !
!   45 :: charge_obj                                                   !
!   55 :: ions_obj                                                     !
!   63 :: options_obj                                                  !
!                                                                      !
!----------------------------------------------------------------------!

  IMPLICIT NONE

  PUBLIC

  ! Public parameters
  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)
  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)
  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2    

  ! Types
  TYPE bader_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: volpos_lat, volpos_car
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: volpos_dir
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: volchg, iondist, ionchg
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: minsurfdist, ionvol
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: volnum, known
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: path
    INTEGER,ALLOCATABLE,DIMENSION(:) :: nnion, vol
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: mask
    REAL(q2) :: stepsize, tol
    REAL(q2) :: vacchg, vacvol
    INTEGER nvols, pnum, bnum, pdim, bdim, refine_edge_itrs, mnum
  END TYPE

  TYPE :: charge_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
    REAL(q2),DIMENSION(3,3) :: lat2car, car2lat
    REAL(q2),DIMENSION(-1:1,-1:1,-1:1) :: lat_dist, lat_i_dist
    REAL(q2),DIMENSION(3) :: org_lat, org_car
    REAL(q2),DIMENSION(3) :: i_npts
    INTEGER,DIMENSION(3) :: npts
    INTEGER :: nrho
  END TYPE

  TYPE :: ions_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir,r_lat
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ion_chg
    REAL(q2),DIMENSION(3,3) :: lattice,dir2car,car2dir
    INTEGER :: nions
    REAL(q2) :: scalefactor
  END TYPE
  
  TYPE :: options_obj
    REAL(q2) :: badertol, stepsize, vacval
    INTEGER :: bader_opt, bader_offgrid = 0, bader_ongrid = 1
    INTEGER :: bader_neargrid = 2, bader_weight = 3
    INTEGER :: quit_opt, quit_max = 0, quit_known = 1
    INTEGER :: refine_edge_itrs
    LOGICAL :: bader_flag, voronoi, dipole, ldos_flag
    LOGICAL :: vac_flag, critpoints
  END TYPE options_obj

END MODULE kind_mod

