MODULE bader_wrap

!----------------------------------------------------------------------!
!                                                                      !
! bader_wrap:  python wrapper for streamlined version of the Bader     !
!              charge density analysis program.                        !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! source info:  http://theory.cm.utexas.edu/henkelman/code/bader/      !
!                                                                      !
! Based on algorithms described in the following publications:         !
!                                                                      !
!    A fast and robust algorithm for Bader decomposition of            !
!    charge density                                                    !
!    G. Henkelman, A. Arnaldsson, and H. Jonsson                       !
!    Comput. Mater. Sci. 36, 254-360 (2006).                           !
!                                                                      !
!    An improved grid-based algorithm for Bader charge allocation      !
!    E. Sanville, S. Kenny, R. Smith, and G. Henkelman                 !
!    J. Comput. Chem. 28, 899-908 (2007).                              !
!                                                                      !
!    A grid-based Bader analysis algorithm without lattice bias        !
!    W. Tang, E. Sanville, and G. Henkelman                            !
!    J. Phys.: Condens. Matter 21, 084204 (2009)                       !
!                                                                      !
!    Accurate and efficient algorithm for Bader charge integration     !
!    M. Yu and D. Trinkle                                              !
!    J. Chem. Phys. 134, 064111 (2011).                                !
!                                                                      !
!----------------------------------------------------------------------!

  USE kind_mod
  USE interface_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bader_run

  CONTAINS

  SUBROUTINE bader_run(nx,ny,nz,lat,ni,rho,rrho,ions,arg,callback,mask)

  !--------------------------------------------------------------------!
  ! bader_run:  python collector passing to wrapper of Bader function. !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! args:                                                              !
  !                                                                    !
  ! nx,ny,nz : grid dimensions.                                        !
  ! lattice  : 3x3 array of lattice vectors.                           !
  ! ni       : total number of ions.                                   !
  ! rho      : (nx,ny,nz) shaped array of charge density.              !
  ! rrho     : (nx,ny,nz) shaped array of reference charge density.    !
  ! ions     : (num_ion,3) array of site fractional coords.            !
  ! arg      : character string of command line flags.                 !
  ! callback : python function to capture output.                      !
  !                                                                    !
  !--------------------------------------------------------------------!

    INTEGER :: nx,ny,nz,ni
    REAL(q2) :: lat(3,3),rho(nx,ny,nz),rrho(nx,ny,nz),ions(ni,3)
    CHARACTER(LEN=128) :: arg

    EXTERNAL callback
    EXTERNAL mask

    INTENT(in) :: nx,ny,nz,lat,ni,rho,rrho,ions,arg

    CALL main_wrap(nx,ny,nz,lat,ni,rho,rrho,ions,arg,callback,mask)

  RETURN
  END SUBROUTINE bader_run

END MODULE bader_wrap
