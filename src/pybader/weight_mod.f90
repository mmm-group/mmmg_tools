MODULE weight_mod

!----------------------------------------------------------------------!
!                                                                      !
! weight_mod:  module implementing the Yu and Trinkle weight method    !
!              [JCP 134, 064111 (2011)].                               !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! contents:                                                            !
!                                                                      !
!   43 :: bader_weight_calc                                            !
!  307 :: ws_voronoi                                                   ! 
!  427 :: incell                                                       !
!  496 :: sort_weight                                                  !
!  538 :: sort_vert                                                    !
!                                                                      !
!----------------------------------------------------------------------!

  USE kind_mod
  USE matrix_mod
  USE bader_mod
  USE charge_mod 
  IMPLICIT NONE
  PRIVATE

  TYPE weight_obj
    REAL(q2) :: rho
    INTEGER, DIMENSION(3) :: pos
  END TYPE

  TYPE rvert_obj
    REAL(q2), DIMENSION(3) :: r
    REAL(q2) :: phi
  END TYPE

  PUBLIC :: weight_obj
  PUBLIC :: bader_weight_calc
  PUBLIC :: sort_weight

  CONTAINS

  SUBROUTINE bader_weight_calc(bdr, ions, chgval, chgref, opts, mask)

  !--------------------------------------------------------------------!
  ! bader_weight_calc:  calculate the Bader volumes and integrate to   !
  !                     give the total charge in each volume.          !
  !--------------------------------------------------------------------!


    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chgval, chgref
    TYPE(options_obj) :: opts
    TYPE(weight_obj),ALLOCATABlE,DIMENSION(:) :: chgList
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: prob
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: alpha, t, w, bdrvol
    REAL(q2),DIMENSION(3,3) :: cell
    REAL(q2) :: vol, tsum, tw, maxProb
    INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: indList
    INTEGER,ALLOCATABLE,DIMENSION(:,:) :: vect, neigh
    INTEGER,ALLOCATABLE,DIMENSION(:) :: numbelow, basin, above
    INTEGER,DIMENSION(3) :: p, ngrid, pt, ptt
    INTEGER :: nPts, i, n, n1, n2, n3, numVect, nv, nVac, ion
    INTEGER :: t1, t2, cr, cm, nabove, na, nb, tbasin, m, bv, c1
    INTEGER :: tenths_done
    LOGICAL :: boundary

    EXTERNAL mask

    ngrid(:) = chgval%npts(:)

    bdr%nvols = 0
    bdr%bnum = 0

    DO i=1,3
      cell(i,:) = ions%lattice(i,:)/chgval%npts(i)
    END DO

    CALL ws_voronoi(cell, numVect, vect, alpha)
    CALL SYSTEM_CLOCK(t1, cr, cm)

    nPts = chgref%npts(1)*chgref%npts(2)*chgref%npts(3)
    ALLOCATE (numbelow(nPts))
    ALLOCATE (w(nPts))
    ALLOCATE (neigh(nPts, numVect))
    ALLOCATE (prob(nPts, numVect))
    ALLOCATE (basin(nPts))
    ALLOCATE (chgList(nPts))
    ALLOCATE (bdr%volnum(chgref%npts(1),chgref%npts(2),chgref%npts(3)))
    ALLOCATE (indList(chgref%npts(1), chgref%npts(2), chgref%npts(3)))
    bdr%bdim = 64 ! will expand as needed
    ALLOCATE (bdr%volpos_lat(bdr%bdim, 3))

    ! count vacuum points
    bdr%vacchg = 0.0_q2
    nVac = 0
    vol = matrix_volume(ions%lattice)
    IF (opts%vac_flag) THEN
    DO n1 = 1, chgref%npts(1)
      DO n2 = 1, chgref%npts(2)
        DO n3 = 1, chgref%npts(3)
          IF (rho_val(chgref,n1,n2,n3)/vol <= opts%vacval) THEN
             bdr%volnum(n1,n2,n3) = -1
             bdr%vacchg = bdr%vacchg + chgval%rho(n1,n2,n3)
             nVac = nVac + 1
          END IF
        END DO
      END DO
    END DO
    END IF
    bdr%vacchg = bdr%vacchg/REAL(chgval%nrho,q2)
    bdr%vacvol = nVac*vol/chgval%nrho

    ! assign points to chgList
    n = 1
    DO n1 = 1, chgref%npts(1)
      DO n2 = 1, chgref%npts(2)
        DO n3 = 1, chgref%npts(3)
          chgList(n)%rho = chgref%rho(n1,n2,n3)
          chgList(n)%pos = (/n1,n2,n3/)
          n = n + 1
        END DO
      END DO
    END DO

    WRITE(*,'(/,2x,A,$)') 'SORTING CHARGE VALUES ... '
    CALL sort_weight(nPts, chgList) ! max value first

    DO n = 1, nPts
      pt = chgList(n)%pos
      indList(pt(1),pt(2),pt(3)) = n
    END DO

    WRITE(*,'(A)') 'DONE'

    WRITE(*,'(2x,A,$)') 'CALCULATING FLUX ... '
    ALLOCATE (t(numVect))
    ALLOCATE (above(numVect))
    ! first loop, deal with all interior points
    DO n = 1, (nPts - nVac)
      basin(n) = 0
      numbelow(n) = 0
      nabove = 0
      tsum = 0
      DO nv = 1, numVect
        p = chgList(n)%pos + vect(nv,:)
        CALL pbc(p, chgref%npts)
        m = indList(p(1), p(2), p(3))
        IF (m < n) THEN ! point p has higher rho
          nabove = nabove + 1
          above(nabove) = m 
          t(nabove) = alpha(nv)*(chgList(m)%rho - chgList(n)%rho)
          tsum = tsum + t(nabove)
        END IF
      END DO
      IF (nabove == 0) THEN ! maxima
        pt = chgList(n)%pos
        bdr%bnum = bdr%bnum + 1
        bdr%nvols = bdr%nvols + 1
        basin(n) = bdr%nvols
        bdr%volnum(pt(1),pt(2),pt(3)) = bdr%nvols 
        IF (bdr%bnum >= bdr%bdim) THEN
          CALL reallocate_volpos(bdr, bdr%bdim*2)
        END IF
        bdr%volpos_lat(bdr%bnum,:) = REAL(p,q2)
        CYCLE
      END IF
      ! else, either an interior point, or a boundary point
      tbasin = basin(above(1))
      boundary = .FALSE.
      DO na = 1, nabove
        IF (basin(above(na))/=tbasin .OR. tbasin==0) boundary = .TRUE.
      END DO 
      IF (boundary) THEN ! boundary
        basin(n) = 0
        maxProb = 0
        DO na = 1, nabove
          pt = chgList(n)%pos
          m = above(na)
          ptt =chgList(m)%pos
          numbelow(m) = numbelow(m) + 1
          neigh(m, numbelow(m)) = n
          prob(m, numbelow(m)) = t(na) / tsum
          IF (prob(m,numbelow(m)) > maxProb) THEN
            maxProb = prob(m,numbelow(m))
            bdr%volnum(pt(1),pt(2),pt(3)) = &
              bdr%volnum(ptt(1),ptt(2),ptt(3))
          END IF
        END DO
      ELSE ! interior
        basin(n) = tbasin
        pt = chgList(n)%pos
        bdr%volnum(pt(1),pt(2),pt(3)) = tbasin
      END IF
    END DO
    DEALLOCATE(t)
    DEALLOCATE(above)

    ! restore chglist rho to values from chgval
    DO n = 1, (nPts - nVac)
      pt = chgList(n)%pos
      chgList(n)%rho = chgval%rho(pt(1),pt(2),pt(3))
    END DO
    WRITE(*,'(A)') 'DONE'

    WRITE(*,'(/,2x,A)') 'INTEGRATING CHARGES'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'
    tenths_done = 0

    ALLOCATE (bdr%volchg(bdr%nvols))
    ALLOCATE (bdrvol(bdr%nvols))
    ALLOCATE (bdr%vol(bdr%nvols))
    bdr%volchg = 0
    bdrvol = 0
    bdr%vol = 0

    DO bv = 1, bdr%nvols
      IF ((bv*10/bdr%nvols) > tenths_done) THEN
        tenths_done = (bv*10/bdr%nvols)
        WRITE(*,'(A,$)') '**'
      END IF
      DO n = 1, nPts
        IF (basin(n) == bv) THEN
          w(n) = 1
          bdr%vol(bv) = bdr%vol(bv) + 1
        ELSE
          w(n) = 0
        END IF
      END DO
      c1 = 0
      ALLOCATE(bdr%mask(2,bdr%vol(bv)))
      bdr%mask = 0._q2
      bdr%mnum = bdr%vol(bv)
      DO n = 1, nPts
        tw = w(n)
        IF (tw /= 0) THEN
          c1 = c1 + 1
          DO nb = 1, numbelow(n)
            w(neigh(n, nb)) = w(neigh(n, nb)) + prob(n, nb)*tw
          END DO
          p = chgList(n)%pos
          p = (/p(1)-1,p(2)-1,p(3)-1/)
          bdr%volchg(bv) = bdr%volchg(bv) + tw * chgList(n)%rho
          bdrvol(bv) = bdrvol(bv) + tw
          i=(((p(1)*chgval%npts(2))+p(2))*chgval%npts(3))+p(3)
          IF (bdr%mnum < c1) CALL reallocate_mask(bdr,2*bdr%mnum)
          bdr%mask(1,c1) = REAL(i,q2)
          bdr%mask(2,c1) = REAL(tw,q2)
        END IF
      END DO
      bdr%mnum = c1
      CALL reallocate_mask(bdr,c1)
      CALL mask(c1,bdr%mask)
      DEALLOCATE(bdr%mask)
    END DO
    bdr%volchg = bdr%volchg / REAL(chgval%nrho,q2)

    CALL SYSTEM_CLOCK(t2, cr, cm)
    WRITE(*,'(/,1A12,1F10.2,1A1)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),'s'

    ! reassign vacuum points from -1 to bdr%bnum + 1
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

    ! check that everything got assigned
    DO n = 1, nPts
      pt = chgList(n)%pos
      IF (bdr%volnum(pt(1),pt(2),pt(3)) == 0) THEN
        write(*,*) n,chgList(n)%pos(:),chgList(n)%rho
        STOP "some volume not assigned"
      END IF
    END DO

    ALLOCATE(bdr%nnion(bdr%nvols))
    ALLOCATE(bdr%iondist(bdr%nvols))
    ALLOCATE(bdr%ionchg(ions%nions))
    ALLOCATE(bdr%volpos_dir(bdr%nvols,3),bdr%volpos_car(bdr%nvols,3))

    DO i = 1, bdr%nvols
      bdr%volpos_dir(i,:) = lat2dir(chgref, bdr%volpos_lat(i,:))
      bdr%volpos_car(i,:) = lat2car(chgref, bdr%volpos_lat(i,:))
    END DO

    CALL assign_chg2atom(bdr, ions)

    ! Calculate atomic volumes from bader volumes

    ALLOCATE(bdr%ionvol(ions%nions))
    bdr%ionvol = 0.0

    DO bv = 1, bdr%nvols
      ion = bdr%nnion(bv)
      bdr%ionvol(ion) = bdr%ionvol(ion) + bdrvol(bv)
    END DO

    vol = matrix_volume(ions%lattice)
    vol = vol/chgref%nrho

    DO ion = 1, ions%nions
      bdr%ionvol(ion) = bdr%ionvol(ion)*vol
    END DO


    DEALLOCATE (numbelow, w, neigh, prob, bdrvol)
    DEALLOCATE (chgList, indList, basin)

  END SUBROUTINE bader_weight_calc

  SUBROUTINE ws_voronoi(cell, numVect, vect, alpha)

  !--------------------------------------------------------------------!
  ! source:  adapted from ws_voronoi.H                                 !
  ! author:  D. Trinkle                                                !
  ! date:    2010 December 27                                          !
  ! purpose: Determines the prefactors for computation of flux in      !
  !          Wigner-Seitz grid cells, based on the Voronoi             !
  !          decomposition.                                            !
  !--------------------------------------------------------------------!
  !                                                                    !
  ! Construct a list of the neighboring vectors that define the        !
  ! Wigner-Seitz cell and compute the "alpha" factors needed for flux; !
  ! you multiply the difference in densities by alpha, and use this to !
  ! compute the transition probabilities.                              !
  !                                                                    !
  !--------------------------------------------------------------------!

    REAL(q2), DIMENSION(3,3) :: cell
    INTEGER, INTENT(INOUT) :: numVect
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vect
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alpha

    TYPE(rvert_obj), ALLOCATABLE, DIMENSION(:) :: rVert
    REAL(q2), ALLOCATABLE, DIMENSION(:,:) :: R, Rtmp
    REAL(q2), ALLOCATABLE, DIMENSION(:) :: alphtmp, Rmag, Rmagtmp
    REAL(q2), DIMENSION(3,3) :: Rdot, Radj
    REAL(q2), DIMENSION(3) :: R2, Rx, Ry, nv
    REAL(q2) :: detR, tol, rdRn
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nVect, nVecttmp
    INTEGER :: nv1, nv2, nv3, nA, nB, n, nvi
    INTEGER :: nRange, numVert, maxVert, neigh, numNeigh, maxNeigh
    LOGICAL :: zeroarea

    tol = 1E-8
    nRange = 3
    maxNeigh = (2*nRange + 1)**3 - 1

    ALLOCATE (nVect(maxNeigh,3), nVecttmp(maxNeigh,3))
    ALLOCATE (R(maxNeigh,3), Rtmp(maxNeigh,3))
    ALLOCATE (Rmag(maxNeigh), Rmagtmp(maxNeigh))

    neigh = 0
    DO nv1 = -nRange, nRange
      DO nv2 = -nRange,nRange
        DO nv3 = -nRange,nRange
          nv = (/nv1,nv2,nv3/)
          IF (ALL(nv == 0)) CYCLE
          neigh = neigh + 1
          nVect(neigh,:) = nv
          R(neigh,:) = MATMUL(cell, nv)
          Rmag(neigh) = SUM(R(neigh,:)*R(neigh,:))*0.5_q2
        END DO
      END DO
    END DO

    ! find the number of neighboring vectors in the WS cell, numNeigh
    numNeigh = 0
    DO neigh = maxNeigh, 1, -1
      ! check to see if R/2 is within the WS cell
      IF( incell(R(neigh,:)*0.5_q2, maxNeigh, R, Rmag, tol) ) THEN
        numNeigh = numNeigh + 1
        Rtmp(numNeigh,:) = R(neigh,:)
        Rmagtmp(numNeigh) = Rmag(neigh)
        nVecttmp(numNeigh,:) = nVect(neigh,:)
      END IF
    END DO

    DEALLOCATE(R, Rmag, nVect)
    ALLOCATE(R(numNeigh,3), nVect(numNeigh,3), Rmag(numNeigh))
    DO neigh = 1, numNeigh
      R(NumNeigh+1 - neigh,:) = Rtmp(neigh,:)
      Rmag(NumNeigh+1 - neigh) = Rmagtmp(neigh)
      nVect(NumNeigh+1 - neigh,:) = nVecttmp(neigh,:)
    END DO
    DEALLOCATE(Rtmp, Rmagtmp, nVecttmp)

    ! next step is to find all of the vertex points
    maxVert = (numNeigh-2)*(numNeigh-4)
    ALLOCATE (rVert(maxVert))
    ALLOCATE (alphtmp(numNeigh))
    numVect = numNeigh
    DO neigh = 1, numNeigh
      numVert = 1
      Rdot(1,:) = R(neigh,:)
      R2(1) = Rmag(neigh)
      DO nA = 1, numNeigh
        Rdot(2,:) = R(nA,:)
        R2(2) = Rmag(nA)
        DO nB = nA + 1, numNeigh
          Rdot(3,:) = R(nB,:)
          R2(3) = Rmag(nB)
          detR = determinant(Rdot) 
          Radj = adjoint(Rdot)

          IF (ABS(detR) >= tol) THEN
            rVert(numVert)%r = MATMUL(Radj, R2)/detR
            ! check if this vertex is in the WS cell
            IF ( incell(rVert(numVert)%r, numNeigh, R, Rmag, tol) ) THEN
                ! inside the cell
              numVert = numVert + 1
            END IF
          END IF
        END DO
      END DO

      !check to make sure none of the vertices correspond to R/2:
      zeroarea = .False.
      DO n = 1, numVert
       zeroarea = (ABS(SUM(rVert(n)%r(:)**2) - 0.5*Rmag(neigh)) < tol)
      END DO
      IF (zeroarea .OR. numVert == 0) THEN
        alphtmp(neigh) = 0
        numVect = numVect - 1
        CYCLE
      END IF

      ! Now we have a list of all the vertices for the polygon
      ! defining the facet along the direction R[n].
      ! Last step is to sort the list in terms of a winding angle around
      ! R[n].  To do that, we define rx and ry which are perpendicular
      ! to R[n], normalized, and right-handed: ry = R x rx, so that
      ! rx x ry points along R[n].

      Rx = rvert(1)%r
      rdRn = SUM(rx(:)*R(neigh,:)) / SUM(R(neigh,:)**2)
      Rx(:) = Rx(:) - rdRn*R(neigh,:)
      rdRn = SQRT(SUM(Rx(:)**2))
      Rx = Rx/rdRn
      Ry = cross_product(R(neigh,:), Rx)
      rdRn = SQRT(SUM(Ry(:)**2))
      Ry = Ry/rdRn
      DO nvi = 1, numVert - 1
        rvert(nvi)%phi = ATAN2(SUM(rvert(nvi)%r(:)*Ry(:)), &
          SUM(rvert(nvi)%r(:)*Rx(:)))
      END DO
      CALL sort_vert(numVert - 1, rvert)
      alphtmp(neigh) = 0
      DO nvi = 1, numVert - 1
        alphtmp(neigh) = alphtmp(neigh) + &
                         triple_product(rvert(nvi)%r, &
                         rvert(MOD(nvi,(numVert-1))+1)%r, R(neigh,:))
      END DO
      alphtmp(neigh) = alphtmp(neigh)*0.25/Rmag(neigh)
      IF (ABS(alphtmp(neigh)) < tol) THEN
        alphtmp(neigh) = 0
        numVect = numVect - 1
      END IF
    END DO 

    ! assign the vertex array with the known number of verticies
    ALLOCATE (vect(numVect,3))
    ALLOCATE (alpha(numVect))

    nvi = 1
    DO n = 1, numNeigh
      IF (alphtmp(n) /= 0 ) THEN
        vect(nvi,:) = nVect(n,:)
        alpha(nvi) = alphtmp(n)
        nvi = nvi + 1
      END IF
    END DO 

  END SUBROUTINE ws_voronoi

  FUNCTION incell(r, numNeigh, rNeigh, rmag, tol)

  !--------------------------------------------------------------------!
  ! incell:  is a r inside the WS cell defined by the vectors R        !
  !          rmag = R.R/2.                                             !
  !--------------------------------------------------------------------!

    REAL(q2), INTENT(IN), DIMENSION(:) :: r, rmag
    REAL(q2), INTENT(IN), DIMENSION(:,:) :: rNeigh
    LOGICAL incell
    INTEGER n, numNeigh
    REAL(q2) tol

    incell = .TRUE.
    DO n = 1, numNeigh
      IF (DOT_PRODUCT(r(:), rNeigh(n,:)) > (rmag(n) + tol)) THEN
        incell = .FALSE.
        EXIT
      END IF
    END DO

    RETURN
  END FUNCTION

  SUBROUTINE sort_weight(array_size, weightList)

  !--------------------------------------------------------------------!
  ! sort_weight:  sort a charge list.                                  !
  !--------------------------------------------------------------------!

    INTEGER, INTENT(in) :: array_size
    TYPE(weight_obj),INTENT(inout),DIMENSION(array_size) :: weightList
    INCLUDE "qsort_inline.inc"
    
    CONTAINS
      SUBROUTINE init()
      END SUBROUTINE init
      LOGICAL &
  
      FUNCTION less_than(a,b)
        INTEGER, INTENT(IN) :: a,b
        IF ( weightList(a)%rho == weightList(b)%rho ) then
          less_than = a < b
        ELSE
          less_than = weightList(a)%rho > weightList(b)%rho ! max first
        END IF
      END FUNCTION less_than
  
      SUBROUTINE swap(a,b)
        INTEGER, INTENT(IN) :: a,b
        TYPE(weight_obj) :: hold
        hold = weightList(a)
        weightList(a) = weightList(b)
        weightList(b) = hold
      END SUBROUTINE swap
  
      ! circular shift-right by one:
      SUBROUTINE rshift(left,right)
        INTEGER, INTENT(in) :: left, right
        TYPE(weight_obj) :: hold
        hold = weightList(right)
        weightList(left+1:right) = weightList(left:right-1)
        weightList(left) = hold
      END SUBROUTINE rshift
  END SUBROUTINE sort_weight

  SUBROUTINE sort_vert(array_size, vertList)

  !--------------------------------------------------------------------!
  ! sort_vert:  sort a vertex list.                                    !
  !--------------------------------------------------------------------!

    INTEGER, INTENT(in) :: array_size
    TYPE(rvert_obj), INTENT(inout), DIMENSION(array_size) :: vertList
    INCLUDE "qsort_inline.inc"

    CONTAINS
      SUBROUTINE init()
      END SUBROUTINE init
      LOGICAL &

      FUNCTION less_than(a,b)
        INTEGER, INTENT(IN) :: a,b
        IF ( vertList(a)%phi == vertList(b)%phi ) THEN
          less_than = a < b
        ELSE
          less_than = vertList(a)%phi < vertList(b)%phi
        END IF
      END FUNCTION less_than

      SUBROUTINE swap(a,b)
        INTEGER, INTENT(IN) :: a,b
        TYPE(rvert_obj) :: hold
        hold = vertList(a)
        vertList(a) = vertList(b)
        vertList(b) = hold
      END SUBROUTINE swap

      ! circular shift-right by one:
      SUBROUTINE rshift(left,right)
        INTEGER, INTENT(IN) :: left, right
        TYPE(rvert_obj) :: hold
        hold = vertList(right)
        vertList(left+1:right) = vertList(left:right-1)
        vertList(left) = hold
      END SUBROUTINE rshift
  END SUBROUTINE sort_vert

END MODULE
