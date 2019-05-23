MODULE bader_mod

!----------------------------------------------------------------------!
!                                                                      !
! bader_mod:  module for analyzing the charge with the Bader atom in   !
!             molecules approach.                                      !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! contents:                                                            !
!                                                                      !
!    51 ::  bader_calc                                                 !
!   279 ::  max_offgrid                                                !
!   328 ::  step_offgrid                                               !
!   349 ::  max_ongrid                                                 !
!   380 ::  step_ongrid                                                !
!   418 ::  max_neargrid                                               !
!   451 ::  step_neargrid                                              !
!   503 ::  refine_edge                                                !
!   668 ::  assign_chg2atom                                            !
!   707 ::  bader_mindist                                              !
!   774 ::  cal_atomic_vol                                             !
!   812 ::  volnum_val                                                 !
!   845 ::  known_volnum                                               !
!   889 ::  assign_surrounding_pts                                     !
!   935 ::  known_volnum_ongrid                                        !
!   967 ::  reassign_volnum_ongrid                                     !
!   993 ::  is_vol_edge                                                !
!  1027 ::  is_atm_edge                                                !
!  1064 ::  is_vol_neighbor                                            !
!  1102 ::  reallocate_volpos                                          !
!  1125 ::  reallocate_path                                            !
!  1148 ::  reallocate_mask                                            !
!                                                                      !
!----------------------------------------------------------------------!

  USE kind_mod
  USE matrix_mod
  USE charge_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: bader_obj
  PUBLIC :: bader_calc, bader_mindist
  PUBLIC :: assign_chg2atom, cal_atomic_vol
  PUBLIC :: reallocate_volpos, reallocate_mask

  CONTAINS

  SUBROUTINE bader_calc(bdr,ions,chgval,chgtemp,opts,mask)

  !--------------------------------------------------------------------!
  ! bader_calc:  calculate the Bader volumes and integrate to give the !
  !              total charge in each volume.                          !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chgval
    TYPE(options_obj) :: opts

    INTEGER,DIMENSION(3) :: p, ptemp
    INTEGER :: n1, n2, n3, i, path_volnum, tenths_done
    INTEGER :: cr, count_max, t1, t2, c1, c2
    INTEGER :: ref_itrs = 1
    REAL(q2),DIMENSION(3) :: voxlen
    REAL(q2) :: vol
    TYPE(charge_obj) :: chgtemp
    EXTERNAL mask

    CALL SYSTEM_CLOCK(t1,cr,count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING BADER CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ! copy bader variable from opts
    bdr%tol = opts%badertol
    bdr%stepsize = opts%stepsize
    IF (opts%stepsize == 0) THEN  ! check for transpose error
      DO i=1,3
        voxlen(i) = SQRT(SUM(chgval%lat2car(:,i)*chgval%lat2car(:,i)))
      END DO
      bdr%stepsize = MINVAL(voxlen(:))
    END IF

    bdr%bdim = 64
    bdr%pdim = 64
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3)) ! will be expanded as needed
    ALLOCATE(bdr%path(bdr%pdim,3))
    ALLOCATE(bdr%volnum(chgval%npts(1),chgval%npts(2),chgval%npts(3)))
    ALLOCATE(bdr%known(chgval%npts(1),chgval%npts(2),chgval%npts(3)))
    bdr%volnum = 0
    bdr%known = 0
    bdr%bnum = 0
    bdr%nvols = 0  ! true number of Bader volumes

    ! find vacuum points, get the charge and volume
    bdr%vacchg = 0.0_q2
    bdr%vacvol = 0.0_q2

    vol = matrix_volume(ions%lattice)
    IF (opts%vac_flag) THEN
      DO n1=1,chgval%npts(1)
        DO n2=1,chgval%npts(2)
          DO n3=1,chgval%npts(3)
            IF (ABS(rho_val(chgtemp,n1,n2,n3)/vol) <= opts%vacval) THEN
               bdr%volnum(n1,n2,n3) = -1
               bdr%vacchg = bdr%vacchg + chgval%rho(n1,n2,n3)
               bdr%vacvol = bdr%vacvol + 1
            END IF
          END DO
        END DO
      END DO
    END IF
    bdr%vacchg = bdr%vacchg/REAL(chgval%nrho,q2)
    bdr%vacvol = bdr%vacvol*vol/chgval%nrho

    tenths_done = 0
    DO n1=1,chgval%npts(1)
      IF ((n1*10/chgval%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chgval%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2=1,chgval%npts(2)
        DO n3=1,chgval%npts(3)
          p = (/n1,n2,n3/)
          IF (bdr%volnum(p(1),p(2),p(3)) == 0) THEN
            IF (opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chgtemp,p)
            ELSEIF (opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chgtemp,p)
            ELSE 
              CALL max_neargrid(bdr,chgtemp,opts,p)
            END IF
            path_volnum = bdr%volnum(p(1),p(2),p(3))

            ! final point has not been assigned, assign new maximum
            IF (path_volnum == 0) THEN
              IF (bdr%bnum >= bdr%bdim) THEN
                CALL reallocate_volpos(bdr,bdr%bdim*2)
              END IF
              bdr%bnum = bdr%bnum+1
              path_volnum = bdr%bnum
              bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
            END IF

            ! assign all points along the trajectory
            DO i = 1,bdr%pnum
              ptemp = (/bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)/)
              IF(bdr%volnum(ptemp(1),ptemp(2),ptemp(3)) /= -1) THEN
                bdr%volnum(ptemp(1),ptemp(2),ptemp(3)) = path_volnum
              END IF
              IF (opts%bader_opt /= opts%bader_ongrid) THEN
                IF (bdr%known(ptemp(1),ptemp(2),ptemp(3)) /= 2) THEN
                  bdr%known(ptemp(1),ptemp(2),ptemp(3)) = 0
                END IF
              END IF
              IF (opts%quit_opt==opts%quit_known &
                & .AND. opts%bader_opt/=opts%bader_ongrid) THEN
                CALL assign_surrounding_pts(bdr,chgtemp,ptemp)
              END IF
            END DO

          END IF
        END DO
      END DO
    END DO
    WRITE(*,*) ''

    IF (opts%vac_flag) THEN
      DO n1 = 1, chgval%npts(1)
        DO n2 = 1, chgval%npts(2)
          DO n3 = 1, chgval%npts(3)
            IF(bdr%volnum(n1,n2,n3) == -1) THEN 
              bdr%volnum(n1,n2,n3) = bdr%bnum + 1
            END IF
          END DO
        END DO
      END DO
    END IF

    ! make a variable which can be changed to indicate convergence
    bdr%refine_edge_itrs = opts%refine_edge_itrs

    IF(opts%refine_edge_itrs < 0) THEN
      WRITE(*,'(/,2x,A)') 'REFINING AUTOMATICALLY'
      DO
        WRITE(*,'(2x,A,I2)') 'ITERATION:',ref_itrs
        CALL refine_edge(bdr,chgtemp,opts,ref_itrs)
        IF (bdr%refine_edge_itrs == 0) EXIT
        ref_itrs = ref_itrs + 1
      END DO
    ELSEIF (opts%refine_edge_itrs > 0) THEN
      WRITE(*,'(/,2x,A)') 'REFINING EDGE'
      DO i=1,opts%refine_edge_itrs
        WRITE(*,'(2x,A,I2)') 'ITERATION:',i
        CALL refine_edge(bdr,chgtemp,opts,i)
      END DO
    ENDIF

    ! total number of bader volumes is now known
    bdr%nvols = bdr%bnum
    CALL reallocate_volpos(bdr, bdr%nvols)
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3))
    ALLOCATE(bdr%volpos_car(bdr%nvols,3))
    DO i=1,bdr%nvols
      bdr%volpos_dir(i,:) = lat2dir(chgtemp, bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:) = lat2car(chgtemp, bdr%volpos_lat(i,:))
    END DO

    ! Sum up the charge included in each volume
    ALLOCATE(bdr%volchg(bdr%nvols))
    ALLOCATE(bdr%vol(bdr%nvols))
    p = chgval%npts
    bdr%volchg = 0._q2
    bdr%vol = 0
    DO n1 = 1, chgval%npts(1)
      DO n2 = 1, chgval%npts(2)
        DO n3 = 1, chgval%npts(3)
          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE
          bdr%volchg(bdr%volnum(n1,n2,n3)) = &
            bdr%volchg(bdr%volnum(n1,n2,n3)) + chgval%rho(n1,n2,n3)
          bdr%vol(bdr%volnum(n1,n2,n3)) = &
            bdr%vol(bdr%volnum(n1,n2,n3)) + 1 
        END DO
      END DO
    END DO
    bdr%volchg = bdr%volchg/REAL(chgval%nrho,q2)

    IF (opts%badermasks) THEN
      WRITE(*,'(/,2x,A)') 'CALCULATING BADER MASKS'
      WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
      WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
      tenths_done = 0
      DO i = 1, bdr%nvols
        IF ((i*10/bdr%nvols) > tenths_done) THEN
          tenths_done = (i*10/bdr%nvols)
          WRITE(*,'(A,$)') '**'
        END IF
        ALLOCATE(bdr%mask(2,bdr%vol(i)))
        bdr%mnum = bdr%vol(i)
        c1 = 0
        c2 = 0
        DO n1 = 1, chgval%npts(1)
          DO n2 = 1, chgval%npts(2)
            DO n3 = 1, chgval%npts(3)
              IF (bdr%volnum(n1,n2,n3) == i) THEN
                c2 = c2 + 1
                IF (bdr%mnum < c2) CALL reallocate_mask(bdr,2*bdr%mnum)
                bdr%mask(1,c2) = REAL(c1,q2)
                bdr%mask(2,c2) = 1._q2
              END IF
              c1 = c1 + 1
            END DO
          END DO
        END DO
        bdr%mnum = c2
        CALL reallocate_mask(bdr,c2)
        CALL mask(bdr%vol(i),bdr%mask)
        DEALLOCATE(bdr%mask)
      END DO
    END IF

    ALLOCATE(bdr%nnion(bdr%nvols))
    ALLOCATE(bdr%iondist(bdr%nvols))
    ALLOCATE(bdr%ionchg(ions%nions))
    CALL assign_chg2atom(bdr,ions)
    CALL cal_atomic_vol(bdr,ions,chgval)
    DEALLOCATE(bdr%path)

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F10.2,1A1)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),'s'

  RETURN
  END SUBROUTINE bader_calc

  SUBROUTINE max_offgrid(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! max_offgrid:  from the point p, do a maximization in the charge    !
  !               density.                                             !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3) :: pm
    REAL(q2) :: cp,cm 
    REAL(q2),DIMENSION(3) :: r

    bdr%pnum = 1
    bdr%path(bdr%pnum,:) = p
    IF (is_max(chg,p)) THEN
      write(*,*) '   max found (init)'
      RETURN
    END IF
    r = REAL(p,q2)

    DO
      cp = rho_val(chg,p(1),p(2),p(3))
      CALL step_offgrid(bdr,chg,r)
      pm = to_lat(chg,r)
      cm = rho_val(chg,pm(1),pm(2),pm(3))

      IF (cm < cp) EXIT
      p = pm

      ! if the point is new, add it to the path
      IF (.NOT.ALL(p(:) == bdr%path(bdr%pnum,:))) THEN
        IF (bdr%pnum >= bdr%pdim) THEN
          CALL reallocate_path(bdr,bdr%pdim*2)
        ENDIF
        bdr%pnum = bdr%pnum+1
        CALL pbc(p,chg%npts)
        bdr%path(bdr%pnum,:) = p
      END IF

      ! quit if this is a maximum or a known volume number
      IF (is_max(chg,p) .OR. known_volnum(bdr,chg,r) /= 0) EXIT

    END DO
    
  RETURN
  END SUBROUTINE max_offgrid

  SUBROUTINE step_offgrid(bdr,chg,r)

  !--------------------------------------------------------------------!
  ! step_offgrid:  step a distance of StepSize along rho_grad.         !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(INOUT) :: r

    REAL(q2),DIMENSION(3) :: grad,dr_car,dr_lat
    REAL(q2) :: rho

    grad = rho_grad(chg,r,rho)
    dr_car = grad*bdr%stepsize/SQRT(SUM(grad*grad))
    dr_lat = MATMUL(chg%car2lat,dr_car)
    r = r+dr_lat

  RETURN
  END SUBROUTINE step_offgrid

  SUBROUTINE max_ongrid(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! max_ongrid:  from the point p do a maximization on the charge      !
  !              density grid and assign the maximum found to the      !
  !              volnum array.                                         !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    bdr%pnum = 1
    bdr%path(bdr%pnum,:) = p
    DO
      CALL step_ongrid(chg,p)
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum = bdr%pnum+1
      CALL pbc(p,chg%npts)
      bdr%path(bdr%pnum,:) = p
      IF (bdr%volnum(p(1),p(2),p(3)) > 0) EXIT
    END DO

  RETURN
  END SUBROUTINE max_ongrid

  SUBROUTINE step_ongrid(chg,p)

  !--------------------------------------------------------------------!
  ! step_ongrid:  do a single iteration of a maximization on the       !
  !               charge density grid from the point (px,py,pz).       !
  !               Return a logical indicating if the current point is  !
  !               a charge density maximum.                            !
  !--------------------------------------------------------------------!

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    REAL(q2) :: rho_max,rho_tmp,rho_ctr
    INTEGER,DIMENSION(3) :: pt,pm
    INTEGER :: d1,d2,d3

    pm = p
    rho_ctr = rho_val(chg,p(1),p(2),p(3))
    rho_max = rho_ctr
    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p+(/d1,d2,d3/)
          rho_tmp = rho_val(chg,pt(1),pt(2),pt(3))
          rho_tmp = rho_ctr+(rho_tmp-rho_ctr)*chg%lat_i_dist(d1,d2,d3)
          IF (rho_tmp > rho_max) THEN
            rho_max = rho_tmp
            pm = pt
          END IF
        END DO
      END DO
    END DO
    CALL pbc(pm,chg%npts)
    p = pm

  RETURN
  END SUBROUTINE step_ongrid

  SUBROUTINE max_neargrid(bdr,chg,opts,p)

  !--------------------------------------------------------------------!
  ! max_neargrid:  from the point p do a maximization on the charge    !
  !                density grid and assign the maximum found to the    !
  !                volnum array.                                       !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

    INTEGER,DIMENSION(3),INTENT(INOUT) :: p

    bdr%pnum=1
    bdr%path(bdr%pnum,:)=p
    DO
      CALL step_neargrid(bdr,chg,p)
      ! if we didn't move, we're at a maximum
      IF (ALL(p == bdr%path(bdr%pnum,:))) EXIT
      ! otherwise, add point to path
      IF (bdr%pnum >= bdr%pdim) THEN
        CALL reallocate_path(bdr,bdr%pdim*2)
      ENDIF
      bdr%pnum = bdr%pnum+1
      bdr%path(bdr%pnum,:) = p
      IF (opts%quit_opt==opts%quit_known &
        & .AND. bdr%known(p(1),p(2),p(3))==2) EXIT
    END DO

  RETURN
  END SUBROUTINE max_neargrid

  SUBROUTINE step_neargrid(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! step_neargrid:  do a single iteration of a maximization on the     !
  !                 charge density grid from the point (px,py,pz).     !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3) :: pm
    REAL(q2),DIMENSION(3) :: gradrl, dr=(/0._q2,0._q2,0._q2/)
    REAL(q2) :: coeff
    SAVE dr

    IF (bdr%pnum == 1) THEN
      dr = (/0._q2,0._q2,0._q2/)
    END IF

    gradrl = rho_grad_dir(chg,p) 

    IF (MAXVAL(ABS(gradrl)) < 1E-30) THEN
      IF (is_max(chg,p)) THEN
        dr = (/0._q2,0._q2,0._q2/)
        RETURN
      ELSE
        pm = p
        CALL step_ongrid(chg,pm)
        dr = (/0._q2,0._q2,0._q2/)
      END IF 
    ELSE
      coeff = 1._q2/MAXVAL(ABS(gradrl))
      gradrl = coeff*gradrl
      pm = p + ANINT(gradrl)
      dr = dr + gradrl - ANINT(gradrl)
      pm = pm + ANINT(dr)
      dr = dr - ANINT(dr)
    END IF
    bdr%known(p(1),p(2),p(3)) = 1

    CALL pbc(pm,chg%npts)
    IF (bdr%known(pm(1),pm(2),pm(3)) == 1) THEN
      pm = p
      CALL step_ongrid(chg,pm)
      dr = (/0._q2,0._q2,0._q2/)
    END IF

    p = pm

  RETURN
  END SUBROUTINE step_neargrid

  SUBROUTINE refine_edge(bdr,chg,opts,ref_itrs)

  !--------------------------------------------------------------------!
  ! refine_edge:  refine the grid points on the edge of the Bader      !
  !               volumes.                                             !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts
    INTEGER :: ref_itrs

    INTEGER,DIMENSION(3) :: p,pt
    INTEGER :: n1,n2,n3,path_volnum,bvolnum,i
    INTEGER :: num_edge,num_reassign,num_check
    INTEGER :: d1,d2,d3,p1,p2,p3

     IF(opts%refine_edge_itrs/=-1 .OR. ref_itrs==1) THEN
       num_edge = 0
       DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            ! change for calculating the vacuum volume
            IF (bdr%volnum(n1,n2,n3) == bdr%bnum+1) CYCLE

            IF (is_vol_edge(bdr,chg,p) .AND. (.NOT.is_max(chg,p))) THEN
              num_edge = num_edge + 1
              bdr%volnum(p(1),p(2),p(3)) = -bdr%volnum(p(1),p(2),p(3))
              IF (opts%quit_opt == opts%quit_known) THEN
                bdr%known(p(1),p(2),p(3)) = 0
                CALL reassign_volnum_ongrid(bdr,chg,p)
              END IF 
            END IF
          END DO
        END DO
      END DO
      WRITE(*,'(2x,A,6x,1I8)') 'EDGE POINTS:',num_edge
    END IF

    IF(opts%refine_edge_itrs==-1 .AND. ref_itrs>1) THEN
      num_check=0
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            ! change for calculating the vacuum volume
            IF (bdr%volnum(n1,n2,n3)==bdr%bnum+1) CYCLE

            IF (bdr%volnum(n1,n2,n3) < 0 &
              & .AND. bdr%known(n1,n2,n3) /=-1) THEN
              DO d1 = -1,1
               DO d2 = -1,1
                DO d3 = -1,1
                  pt = p + (/d1,d2,d3/)
                  CALL pbc(pt,chg%npts)
                  p1 = pt(1)
                  p2 = pt(2)
                  p3 = pt(3)
                  ! change for calculating the vacuum volume
                  IF (bdr%volnum(p1,p2,p3) == bdr%bnum+1) CYCLE
                  IF (.NOT.is_max(chg,pt)) THEN 
                    IF (bdr%volnum(p1,p2,p3) > 0) THEN
                      bdr%volnum(p1,p2,p3) = -bdr%volnum(p1,p2,p3)
                      bdr%known(p1,p2,p3) = -1
                      num_check=num_check+1
                    ELSE IF (bdr%volnum(p1,p2,p3) < 0 &
                           & .AND. bdr%known(p1,p2,p3) == 0) THEN
                      bdr%known(p1,p2,p3) = -2
                      num_check = num_check + 1
                    END IF
                  END IF
                END DO
               END DO
              END DO
              num_check = num_check - 1
              IF (bdr%known(p1,p2,p3) /= -2) THEN
                bdr%volnum(n1,n2,n3) = ABS(bdr%volnum(n1,n2,n3))
              END IF
            END IF

          END DO
        END DO
      END DO
      WRITE(*,'(2x,A,3x,1I8)') 'CHECKED POINTS:', num_check

      ! make the surrounding points unkown
      DO n1 = 1,chg%npts(1)
        DO n2 = 1,chg%npts(2)
          DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            bvolnum = bdr%volnum(n1,n2,n3)

            IF (bvolnum < 0) THEN
              DO d1 = -1,1
               DO d2 = -1,1
                DO d3 = -1,1
                  pt = p + (/d1,d2,d3/)
                  CALL pbc(pt,chg%npts)
                  p1 = pt(1)
                  p2 = pt(2)
                  p3 = pt(3)
                  IF(bdr%known(p1,p2,p3) == 2) bdr%known(p1,p2,p3) = 0
                END DO
               END DO
              END DO
            END IF

          END DO
        END DO
      END DO

    END IF

    num_reassign = 0
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          p = (/n1,n2,n3/)
          bvolnum = bdr%volnum(n1,n2,n3)
          IF (bvolnum<0) THEN
            IF (opts%bader_opt == opts%bader_offgrid) THEN
              CALL max_offgrid(bdr,chg,p)
            ELSEIF (opts%bader_opt == opts%bader_ongrid) THEN
              CALL max_ongrid(bdr,chg,p)
            ELSE
              CALL max_neargrid(bdr,chg,opts,p)
            END IF
            path_volnum = bdr%volnum(p(1),p(2),p(3))
            IF (path_volnum<0 .OR. path_volnum>bdr%bnum) THEN
              WRITE(*,*) 'ERROR: new maxima in edge refinement'
            END IF
            bdr%volnum(n1,n2,n3) = path_volnum
            IF (ABS(bvolnum) /= path_volnum) THEN
              num_reassign = num_reassign + 1
              IF (opts%refine_edge_itrs == -1 &
                & .OR. opts%refine_edge_itrs==-3) THEN 
                bdr%volnum(n1,n2,n3) = -path_volnum
              END IF
            END IF
            DO i = 1,bdr%pnum
              pt = (/bdr%path(i,1),bdr%path(i,2),bdr%path(i,3)/)
              IF (bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
                bdr%known(pt(1),pt(2),pt(3)) = 0
              END IF
            END DO 
          END IF
        END DO
      END DO
    END DO
 
    WRITE(*,'(2x,A,1I8)') 'REASSIGNED POINTS:',num_reassign

    ! flag to indicate that we are done refining
    IF ((opts%refine_edge_itrs == -1 &
      & .OR. opts%refine_edge_itrs==-3) .AND. num_reassign==0) THEN
      bdr%refine_edge_itrs = 0
    END IF
    IF (opts%refine_edge_itrs==-2 .AND. num_reassign==0) THEN
      bdr%refine_edge_itrs = 0
    END IF

  RETURN
  END SUBROUTINE refine_edge

  SUBROUTINE assign_chg2atom(bdr,ions)

  !--------------------------------------------------------------------!
  ! assign_chg2atom:  assign an element of charge to a Bader atom.     !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions

    REAL(q2),DIMENSION(3) :: dv, v
    REAL(q2) :: dsq, dminsq
    INTEGER :: i, j, dindex

    bdr%ionchg = 0._q2

    DO i = 1,bdr%nvols
      dv = bdr%volpos_dir(i,:) - ions%r_dir(1,:)
      v = MATMUL(ions%dir2car, dv)
      CALL dpbc_car(ions, v)
      dminsq = DOT_PRODUCT(v, v)
      dindex = 1
      DO j = 2,ions%nions
        dv = bdr%volpos_dir(i,:) - ions%r_dir(j,:)
        v = MATMUL(ions%dir2car, dv)
        CALL dpbc_car(ions, v)
        dsq = DOT_PRODUCT(v, v)
        IF (dsq < dminsq) THEN
          dminsq = dsq
          dindex = j
        END IF
      END DO
      bdr%iondist(i) = SQRT(dminsq)
      bdr%nnion(i) = dindex
      bdr%ionchg(dindex) = bdr%ionchg(dindex) + bdr%volchg(i)
    END DO

  RETURN
  END SUBROUTINE assign_chg2atom

  SUBROUTINE bader_mindist(bdr,ions,chg)

  !--------------------------------------------------------------------!
  ! bader_mindist:  find the minimum distance from the surface of each !
  !                 volume to each ion.                                !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: v, dv_dir, dv_car
    INTEGER,DIMENSION(3) :: p
    REAL :: dist
    INTEGER :: i, atom, n1, n2, n3
    INTEGER :: cr, count_max, t1, t2, tenths_done
    
    CALL SYSTEM_CLOCK(t1, cr, count_max)

 
    WRITE(*,'(/,2x,A)') 'CALCULATING MINIMUM DISTANCES TO ATOMS'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ! store the minimum distance and the vector
    ALLOCATE(bdr%minsurfdist(ions%nions))
    bdr%minsurfdist = 0._q2

    tenths_done = 0
    DO n1 = 1,chg%npts(1)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          p = (/n1,n2,n3/)
          ! change for calculating the vacuum volume
          IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE

          ! if an edge cell is it the closest to the atom so far
          IF(is_atm_edge(bdr,chg,p,atom)) THEN 
            v = REAL((/n1,n2,n3/),q2)
            dv_dir = (v-chg%org_lat)/REAL(chg%npts,q2) 
            dv_dir = dv_dir - ions%r_dir(atom,:)
            dv_car = MATMUL(ions%dir2car, dv_dir)
            CALL dpbc_car(ions, dv_car)
            dist = DOT_PRODUCT(dv_car, dv_car)
            IF ((bdr%minsurfdist(atom) == 0.0_q2) .OR.  &
            &   (bdr%minsurfdist(atom) > dist)) THEN
              bdr%minsurfdist(atom) = dist
            END IF
          END IF
        END DO
      END DO
    END DO

    DO i = 1,ions%nions
      bdr%minsurfdist(i) = SQRT(bdr%minsurfdist(i))
    END DO

    CALL SYSTEM_CLOCK(t2,cr,count_max)
    WRITE(*,'(/,1A12,1F10.2,1A1)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2), 's'

  RETURN
  END SUBROUTINE bader_mindist

  SUBROUTINE cal_atomic_vol(bdr,ions,chg)

  !--------------------------------------------------------------------!
  ! cal_atomic_vol:  integrate the atomic volume for each atom.        !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2):: vol
    INTEGER :: n1, n2, n3, i, atom

    ALLOCATE(bdr%ionvol(ions%nions))
    bdr%ionvol = 0.0

    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
          ! change to calculate vacuum volume
          IF (bdr%volnum(n1,n2,n3) /= bdr%nvols + 1) THEN
            atom = bdr%nnion(bdr%volnum(n1,n2,n3))
            bdr%ionvol(atom) = bdr%ionvol(atom) + 1
          END IF
        END DO
      END DO
    END DO

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho

    DO i = 1,ions%nions
      bdr%ionvol(i) = bdr%ionvol(i)*vol
    END DO

  RETURN
  END SUBROUTINE cal_atomic_vol

  FUNCTION volnum_val(bdr,chg,p1,p2,p3)

  !--------------------------------------------------------------------!
  ! volnum_val:  return the density at the point (p1,p2,p3) taking     !
  !              into account the boundary conditions. This function   !
  !              is used to address points outside the charge density  !
  !              array without a bunch of if statements at the place   !
  !              the value is needed.                                  !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,INTENT(IN) :: p1, p2, p3
    INTEGER :: i, volnum_val
    INTEGER,DIMENSION(3) :: p

    p = (/p1,p2,p3/)
    DO i = 1,3
      DO
        IF (p(i) >= 1) EXIT
        p(i) = p(i) + chg%npts(i)
      END DO
      DO
        IF (p(i) <= chg%npts(i)) EXIT
        p(i) = p(i) - chg%npts(i)
      END DO
    END DO

    volnum_val = bdr%volnum(p(1),p(2),p(3))

  RETURN
  END FUNCTION volnum_val

  FUNCTION known_volnum(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! known_volnum:  return number of the associated bader volnum if all !
  !                surrounding !grid points are known to be associated !
  !                with the same bader volnum.                         !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    INTEGER :: known_volnum

    INTEGER :: volnum, d1, d2, d3, p1, p2, p3,pt(3)
    LOGICAL :: first_flag

    known_volnum = 0
    first_flag = .TRUE.

    p1 = FLOOR(p(1))
    p2 = FLOOR(p(2))
    p3 = FLOOR(p(3))
    pt = (/p1,p2,p3/)
    DO d1 = 0,1
      p1 = pt(1)+d1
      DO d2 = 0,1
        p2 = pt(2)+d2
        DO d3 = 0,1
          p3 = pt(3)+d3
          IF (first_flag) THEN
            volnum = volnum_val(bdr,chg,p1,p2,p3)
            IF (volnum <= 0) RETURN
            first_flag = .FALSE.
          ELSE
            IF (volnum /= volnum_val(bdr,chg,p1,p2,p3)) RETURN
          END IF
        END DO
      END DO
    END DO
    known_volnum = volnum

  RETURN
  END FUNCTION known_volnum

  SUBROUTINE assign_surrounding_pts(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! assign_surrounding_pts:  check the surrounding points of p to see  !
  !                          if their volnum is known.                 !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) :: pt

    pt = p + (/1,0,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN 
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/-1,0,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,1,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,-1,0/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,0,1/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF
    pt = p + (/0,0,-1/)
    CALL pbc(pt,chg%npts)
    IF(bdr%known(pt(1),pt(2),pt(3)) /= 2) THEN
      CALL known_volnum_ongrid(bdr,chg,pt)
    END IF

  RETURN
  END SUBROUTINE assign_surrounding_pts

  SUBROUTINE known_volnum_ongrid(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! known_volnum_ongrid:  return number of the associated bader volnum !
  !                       if nearest grid points are known to be       !
  !                       associated with the same bader volnum.       !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER :: volnum, p1, p2, p3

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)     

    volnum = volnum_val(bdr,chg,p1,p2,p3)
    IF(volnum <= 0) RETURN

    IF (volnum_val(bdr,chg,p1,p2,p3+1) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2,p3-1) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2+1,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1,p2-1,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1+1,p2,p3) /= volnum) RETURN
    IF (volnum_val(bdr,chg,p1-1,p2,p3) /= volnum) RETURN

    bdr%known(p1,p2,p3) = 2

  RETURN
  END SUBROUTINE known_volnum_ongrid

  SUBROUTINE reassign_volnum_ongrid(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! reassign_volnum_ongrid:  reassign the surrounding points of a edge !
  !                          point as unknown points.                  !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) :: pt
    INTEGER :: d1, d2, d3

    DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          bdr%known(pt(1),pt(2),pt(3)) = 0
        END DO
      END DO
    END DO

  RETURN
  END SUBROUTINE reassign_volnum_ongrid

  FUNCTION is_vol_edge(bdr,chg,p)

  !--------------------------------------------------------------------!
  ! is_vol_edge:  return .true. if the grid point is on the edge of a  !
  !               Bader volume.                                        !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_vol_edge

    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) ::pt
    INTEGER :: d1,d2,d3,volnum,volnbr

    volnum = bdr%volnum(p(1),p(2),p(3))
    is_vol_edge = .FALSE.
    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          volnbr = bdr%volnum(pt(1),pt(2),pt(3))
          IF (ABS(volnbr) /= ABS(volnum)) THEN
            is_vol_edge = .TRUE.
            EXIT neighborloop  
          END IF
        END DO
      END DO
    END DO neighborloop

  RETURN
  END FUNCTION is_vol_edge

  FUNCTION is_atm_edge(bdr,chg,p,atom)

  !--------------------------------------------------------------------!
  ! is_atm_edge:  return .true. if the grid point is on the edge of a  !
  !               Bader atom.                                          !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_atm_edge

    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) ::pt
    INTEGER :: d1, d2, d3, atmnbr
    INTEGER,INTENT(INOUT) ::atom 
    
    atom = bdr%nnion(bdr%volnum(p(1),p(2),p(3)))
    is_atm_edge = .FALSE.
    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          ! check if vacuum
          IF (bdr%volnum(pt(1),pt(2),pt(3)) == bdr%bnum + 1) CYCLE 
          atmnbr = bdr%nnion(bdr%volnum(pt(1),pt(2),pt(3)))
          IF (atmnbr /= atom) THEN
            is_atm_edge = .TRUE.
            EXIT neighborloop
          END IF
        END DO
      END DO
    END DO neighborloop
    
    RETURN
    END FUNCTION is_atm_edge

  FUNCTION is_vol_neighbor(bdr,chg,p,vol)

  !--------------------------------------------------------------------!
  ! is_vol_neighbor:  return .true. if the grid point neighboring the  !
  !                   bader volume.                                    !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    LOGICAL :: is_vol_neighbor

    INTEGER,DIMENSION(3),INTENT(IN) :: p
    INTEGER,DIMENSION(3) :: pt
    INTEGER :: d1, d2, d3, volneighbor
    INTEGER,INTENT(IN) :: vol

    is_vol_neighbor = .FALSE.

    ! only find neighbors, not points in the volume
    IF (bdr%volnum(pt(1),pt(2),pt(3)) == vol) RETURN

    neighborloop: DO d1 = -1,1
      DO d2 = -1,1
        DO d3 = -1,1
          pt = p + (/d1,d2,d3/)
          CALL pbc(pt,chg%npts)
          volneighbor = bdr%volnum(pt(1),pt(2),pt(3))
          IF (volneighbor == vol) THEN
            is_vol_neighbor = .TRUE.
            EXIT neighborloop
          END IF
        END DO
      END DO
    END DO neighborloop

    RETURN
    END FUNCTION is_vol_neighbor

  SUBROUTINE reallocate_volpos(bdr,newsize)

  !--------------------------------------------------------------------!
  ! reallocate_volpos:  extend volpos array.                           !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmpvolpos

    IF(newsize<bdr%bnum) write(*,*) 'Error: new volpos length too small'

    ALLOCATE(tmpvolpos(bdr%bdim,3))
    tmpvolpos = bdr%volpos_lat
    DEALLOCATE(bdr%volpos_lat)
    bdr%bdim = newsize
    ALLOCATE(bdr%volpos_lat(bdr%bdim,3))
    bdr%volpos_lat(1:bdr%bnum,:) = tmpvolpos(1:bdr%bnum,:)
    DEALLOCATE(tmpvolpos)

  END SUBROUTINE reallocate_volpos

  SUBROUTINE reallocate_path(bdr,newsize)

  !--------------------------------------------------------------------!
  ! reallocate_path:  extend path array.                               !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize

    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: tmppath

    IF(newsize < bdr%pnum) write(*,*) 'Error: new path length too small'

    ALLOCATE(tmppath(bdr%pdim,3))
    tmppath = bdr%path
    DEALLOCATE(bdr%path)
    bdr%pdim = newsize
    ALLOCATE(bdr%path(bdr%pdim,3))
    bdr%path(1:bdr%pnum,:) = tmppath(1:bdr%pnum,:)
    DEALLOCATE(tmppath)

  END SUBROUTINE reallocate_path

  SUBROUTINE reallocate_mask(bdr,newsize)

  !--------------------------------------------------------------------!
  ! reallocate_mask:  extend mask array.                               !
  !--------------------------------------------------------------------!

    TYPE(bader_obj) :: bdr
    INTEGER :: newsize
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: tmp
    
    IF(newsize < bdr%mnum) WRITE(*,*) 'Error: new mask length too small'

    ALLOCATE(tmp(2,bdr%mnum))
    tmp = bdr%mask
    DEALLOCATE(bdr%mask)
    ALLOCATE(bdr%mask(2,newsize))
    bdr%mask(:,1:bdr%mnum) = tmp(:,1:bdr%mnum)
    DEALLOCATE(tmp)
    bdr%mnum = newsize

  END SUBROUTINE

END MODULE bader_mod
