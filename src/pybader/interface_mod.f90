MODULE interface_mod

!----------------------------------------------------------------------!
!                                                                      !
! interface_mod:  python wrapper for streamlined version of the Bader  !
!                 charge density analysis program.                     !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! contents:                                                            !
!                                                                      !
!   31 :: main_wrap                                                    !
!  112 :: output                                                       ! 
!  142 :: parse_args                                                   !
!  315 :: build_charge                                                 !
!                                                                      !
!----------------------------------------------------------------------!

  USE kind_mod
  USE matrix_mod
  USE charge_mod
  USE bader_mod
  USE weight_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: main_wrap

  CONTAINS

  SUBROUTINE main_wrap(nx,ny,nz,lat,ni,rho,rrho,sites,arg,callback,mask)

  !--------------------------------------------------------------------!
  ! main_wrap:  main Bader function.                                   !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! args:                                                              !
  !                                                                    !
  ! nx,ny,nz : grid dimensions.                                        !
  ! lattice  : 3x3 array of lattice vectors.                           !
  ! i_n      : total number of ions.                                   !
  ! rho      : (nx,ny,nz) shaped array of charge density.              !
  ! rrho     : (nx,ny,nz) shaped array of reference charge density.    !
  ! sites    : (num_ion,3) array of site fractional coords.            !
  ! arg      : character string of command line flags.                 !
  ! callback : python callback function to export output.              !
  ! mask     : python callback function to export masks.               !
  !                                                                    !
  ! returns:                                                           !
  !                                                                    !
  ! ba    : (nx,ny,nz) mask for each bader volume and ion volume.      !
  ! bc    : array of charge for each bader volume.                     !
  ! bi    : array of the corresponding ion for each Bader volume.      !
  ! ic    : array of the charge for each ion.                          !
  ! iv    : array of the volume for each ion.                          !
  ! id    : array of the distance of the maxima to ion site.           !
  ! imsd  : array of the distance of the maxima to the volume surface. !
  ! vc    : charge of the vacuum.                                      !
  ! vv    : volume of the vacuum.                                      !
  !                                                                    !
  !--------------------------------------------------------------------!

    TYPE(options_obj) :: opts
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chgval,chgref
    TYPE(bader_obj) :: bdr

    INTEGER :: nx,ny,nz,ni,npts(3)
    REAL(q2) :: lat(3,3),rho(nx,ny,nz),rrho(nx,ny,nz),sites(ni,3)
    CHARACTER(LEN=128) :: arg

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: bp
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: bc,bd,ic,iv,imsd
    INTEGER,ALLOCATABLE,DIMENSION(:) :: bi
    REAL(q2) :: vc,vv

    EXTERNAL callback
    EXTERNAL mask
    
    INTENT(in) :: nx,ny,nz,lat,ni,rho,sites,arg

    npts = (/nx,ny,nz/)

    CALL parse_args(arg,opts)
    CALL build_charge(lat,ni,npts,rho,sites,ions,chgval)
    chgref = chgval
    chgref%rho(:,:,:) = rrho(:,:,:)

    IF (opts%bader_flag) THEN
      IF (opts%bader_opt == opts%bader_weight) THEN
        CALL bader_weight_calc(bdr,ions,chgval,chgref,opts,mask)
      ELSE
        CALL bader_calc(bdr,ions,chgval,chgref,opts,mask)
      ENDIF
      CALL bader_mindist(bdr,ions,chgval)
    ENDIF

!-- Do I want to add these? callback is already pretty full------------!

!    IF (opts%critpoints) CALL critpoint_find(bdr,chgval,opts,ions)
!    IF (opts%dipole) CALL multipole_calc(bdr,ions,chgval,chgref,opts)
!    IF (opts%voronoi) CALL voronoi(vor,ions,chgval)

!----------------------------------------------------------------------!

    CALL output(bdr,ni,bc,bi,bd,bp,ic,iv,imsd,vc,vv)
    CALL callback(bdr%nvols,ni,bc,bi,bd,bp,ic,iv,imsd,vc,vv)

  RETURN
  END SUBROUTINE main_wrap

  SUBROUTINE output(bdr,ni,bc,bi,bd,bp,ic,iv,imsd,vc,vv)

  !--------------------------------------------------------------------!
  ! output:  prepare output for callback function as f2py can't handle !
  !          derived data types.                                       !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: bc, bd, ic, iv, imsd
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: bp
    INTEGER,ALLOCATABLE,DIMENSION(:) :: bi
    INTEGER :: ni
    REAL(q2) :: vc, vv

    ALLOCATE(bc(bdr%nvols),bi(bdr%nvols),bd(bdr%nvols),bp(bdr%nvols,3))
    ALLOCATE(ic(ni), iv(ni), imsd(ni))

    bc = bdr%volchg
    bi = bdr%nnion
    bd = bdr%iondist
    bp = bdr%volpos_dir
    ic = bdr%ionchg
    iv = bdr%ionvol
    imsd = bdr%minsurfdist
    vc = bdr%vacchg
    vv = bdr%vacvol

  RETURN
  END SUBROUTINE output

  SUBROUTINE parse_args(args,opts)

  !--------------------------------------------------------------------!
  ! parse_args:  parse passed flags and parameters and save them to    !
  !              opts(options_obj)                                     !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! args:                                                              !
  !                                                                    !
  ! arg      :  character string of passed command line flags.         !
  !                                                                    !
  !--------------------------------------------------------------------!

    TYPE(options_obj) :: opts
    CHARACTER(LEN=128) :: args
    INTEGER :: p,ip,it

    ! Default values
    opts%vac_flag = .FALSE.
    opts%vacval = 1E-3
    opts%bader_opt = opts%bader_neargrid
    opts%quit_opt = opts%quit_known
    opts%refine_edge_itrs = -1
    opts%bader_flag = .TRUE.
    opts%voronoi = .FALSE.
    opts%dipole = .FALSE.
    opts%ldos_flag = .FALSE.
    opts%badertol = 1E-3
    opts%stepsize = 0.0_q2
    opts%critpoints = .FALSE.
    opts%badermasks = .FALSE.

    p = 1
    DO 
      ip = index(args(p:),'-')
      IF (ip == 0) EXIT
      ip = p+ip-1
      ! Find critical points
      IF (args(ip:ip+2) == '-cp') THEN
        opts%critpoints = .TRUE.

      ! Vacuum options
      ELSEIF (args(ip:ip+3) == '-vac') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'AUTO') /=0 &
          & .OR. index(args(ip:it),'auto') /= 0) THEN
          opts%vac_flag = .TRUE.
        ELSEIF (index(args(ip:it), 'OFF') /=0 &
              & .OR. index(args(ip:it),'off') /= 0) THEN
          opts%vac_flag = .FALSE.
        ELSE
          READ(args(ip+4:it),*) opts%vacval
          opts%vac_flag = .TRUE.
        END IF

      ! Bader options
      ELSEIF (args(ip:ip+1) == '-b') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'OFFGRID') /=0 &
          & .OR. index(args(ip:it),'offgrid') /= 0) THEN
          opts%bader_opt = opts%bader_offgrid
        ELSEIF (index(args(ip:it), 'ONGRID') /=0 &
              & .OR. index(args(ip:it),'ongrid') /= 0) THEN
          opts%bader_opt = opts%bader_ongrid
          opts%refine_edge_itrs=0
        ELSEIF (index(args(ip:it), 'NEARGRID') /=0 &
              & .OR. index(args(ip:it),'neargrid') /= 0) THEN
          opts%bader_opt = opts%bader_neargrid
        ELSEIF (index(args(ip:it), 'WEIGHT') /=0 &
              & .OR. index(args(ip:it),'weight') /= 0) THEN
          opts%bader_opt = opts%bader_weight
        END IF

      ! Quit options
      ELSEIF (args(ip:ip+1) == '-m') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'MAX') /=0 &
          & .OR. index(args(ip:it),'max') /= 0) THEN
          opts%quit_opt = opts%quit_max
        ELSEIF (index(args(ip:it), 'KNOWN') /=0 &
              & .OR. index(args(ip:it),'known') /= 0) THEN
          opts%quit_opt = opts%quit_known
        END IF

      ! Calculation options
      ELSEIF (args(ip:ip+1) == '-c') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'ALL') /=0 &
          & .OR. index(args(ip:it),'all') /= 0) THEN
          opts%bader_flag = .TRUE.
          opts%voronoi = .TRUE.
          opts%dipole = .TRUE.
          opts%ldos_flag = .TRUE.
        ELSEIF (index(args(ip:it), 'BADER') /=0 &
              & .OR. index(args(ip:it),'bader') /= 0) THEN
          opts%bader_flag = .TRUE.
        ELSEIF (index(args(ip:it), 'VORONOI') /=0 &
              & .OR. index(args(ip:it),'voronoi') /= 0) THEN
          opts%voronoi = .TRUE.
        ELSEIF (index(args(ip:it), 'DIPOLE') /=0 &
              & .OR. index(args(ip:it),'dipole') /= 0) THEN
          opts%dipole = .TRUE.
        ELSEIF (index(args(ip:it), 'LDOS') /=0 &
              & .OR. index(args(ip:it),'ldos') /= 0) THEN
          opts%ldos_flag = .TRUE.
        END IF
        
      ELSEIF (args(ip:ip+1) == '-n') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'ALL') /=0 &
          & .OR. index(args(ip:it),'all') /= 0) THEN
          opts%bader_flag = .FALSE.
          opts%voronoi = .FALSE.
          opts%dipole = .FALSE.
          opts%ldos_flag = .FALSE.
        ELSEIF (index(args(ip:it), 'BADER') /=0 &
              & .OR. index(args(ip:it),'bader') /= 0) THEN
          opts%bader_flag = .FALSE.
        ELSEIF (index(args(ip:it), 'VORONOI') /=0 &
              & .OR. index(args(ip:it),'voronoi') /= 0) THEN
          opts%voronoi = .FALSE.
        ELSEIF (index(args(ip:it), 'DIPOLE') /=0 &
              & .OR. index(args(ip:it),'dipole') /= 0) THEN
          opts%dipole = .FALSE.
        ELSEIF (index(args(ip:it), 'LDOS') /=0 &
              & .OR. index(args(ip:it),'ldos') /= 0) THEN
          opts%ldos_flag = .FALSE.
      END IF

      ! Bader tolerance
      ELSEIF (args(ip:ip+1) == '-t') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        READ(args(ip+2:it),*) opts%badertol

      ! Bader masks
      ELSEIF (args(ip:ip+1) == '-p') THEN
        opts%badermasks = .TRUE.

      ! Refine edge iterations  -- change this to a flag once working
      ELSEIF (args(ip:ip+1) == '-r') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        IF (index(args(ip:it), 'AUTO') /=0 &
          & .OR. index(args(ip:it),'auto') /= 0) THEN
          opts%refine_edge_itrs=-1
        ELSE
          READ(args(ip+2:it-1),'(I16)') opts%refine_edge_itrs
          IF (ip+2 == it-1) THEN
            ip = it
            it = index(args(ip+1:),'-')+ip
            IF (it == ip) it = lnblnk(args) + 1
            READ(args(ip-1:it-1),'(I16)') opts%refine_edge_itrs
          END IF
        END IF

      ! Step size
      ELSEIF (args(ip:ip+1) == '-s') THEN
        it = index(args(ip+1:),'-')+ip
        IF (it == ip) it = lnblnk(args) + 1
        READ(args(ip+2:it),*) opts%stepsize
      END IF
      p = ip+1
    END DO
    
  RETURN
  END SUBROUTINE parse_args

  SUBROUTINE build_charge(lattice,numion,npts,rho,r_dir,ions,chg)

  !--------------------------------------------------------------------!
  ! build_charge:  create charge and ions objects.                     !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! args:                                                              !
  !                                                                    !
  ! lattice  :  character string of passed command line flags.         !
  ! numion   :  total number of ions.                                  !
  ! npts     :  1x3 array of the grid points.                          !
  ! rho      :  (nx,ny,nz) array of the charge density.                !
  ! r_dir    :  (numion,3) array of ion fractional coords.             !
  ! ions     :  ions_obj type.                                         !
  ! chg      :  charge_obj type.                                       !
  !                                                                    !
  !--------------------------------------------------------------------!

    REAL(q2),DIMENSION(3,3) :: lattice   
    INTEGER :: numion
    REAL(q2) :: r_dir(:,:)
    INTEGER,DIMENSION(3) :: npts
    REAL(q2) :: rho(:,:,:)

    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3) :: dlat,dcar
    INTEGER :: i,d1,d2,d3

    ions%scalefactor = 1._q2
    ions%lattice = lattice
    ions%nions = numion
    ions%dir2car = TRANSPOSE(ions%lattice)
    ions%car2dir = inverse(ions%dir2car)
    ALLOCATE(ions%r_dir(ions%nions,3))
    ALLOCATE(ions%r_car(ions%nions,3))

    DO i = 1, ions%nions
      ions%r_dir(i,:) = r_dir(i,:)
      ions%r_car(i,:) = MATMUL(ions%dir2car, ions%r_dir(i,:))
    END DO
 
    chg%npts(:) = npts(:)
    chg%nrho = PRODUCT(chg%npts(:))
    chg%i_npts = 1.0_q2/REAL(chg%npts,q2)
    ALLOCATE(chg%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    chg%rho(:,:,:) = rho(:,:,:)
    DO i=1,3
      chg%lat2car(:,i) = ions%dir2car(:,i)/REAL(chg%npts(i),q2)
    END DO
    chg%car2lat = inverse(chg%lat2car)

    ! origin of the lattice is at chg(1,1,1)
    chg%org_lat = (/1._q2,1._q2,1._q2/)
    chg%org_car = (/0._q2,0._q2,0._q2/)

    ! ion positions in grid points
    ALLOCATE(ions%r_lat(ions%nions,3))
    DO i = 1,ions%nions
      ions%r_lat(i,:) = MATMUL(chg%car2lat, ions%r_car(i,:))
      ions%r_lat(i,:) = ions%r_lat(i,:) + chg%org_lat
      CALL pbc_r_lat(ions%r_lat(i,:), chg%npts)
    END DO

    ! distance between neighboring points
    DO d1 = -1, 1
      dlat(1) = REAL(d1,q2)
      DO d2 = -1, 1
        dlat(2) = REAL(d2,q2)
        DO d3 = -1, 1
          dlat(3) = REAL(d3,q2)
          dcar = MATMUL(chg%lat2car, dlat)
          chg%lat_dist(d1,d2,d3) = SQRT(SUM(dcar*dcar))
          IF ((d1 == 0).AND.(d2 == 0).AND.(d3 == 0)) THEN
            chg%lat_i_dist(d1,d2,d3) = 0._q2
          ELSE
            chg%lat_i_dist(d1,d2,d3) = 1._q2/chg%lat_dist(d1,d2,d3)
          END IF
        END DO
      END DO
    END DO

  RETURN
  END SUBROUTINE build_charge

END MODULE interface_mod
