MODULE matrix_mod

!----------------------------------------------------------------------!
!                                                                      !
! matrix_mod:  module for performing matrix operatoins.                !
!                                                                      !
!----------------------------------------------------------------------!
!                                                                      !
! contents:                                                            !
!                                                                      !
!    28 :: inverse                                                     !
!    45 :: adjoint                                                     ! 
!    69 :: matrix_volume                                               !
!    87 :: triple_product                                              !
!   104 :: cross_product                                               !
!   120 :: determinant                                                 !
!                                                                      !
!----------------------------------------------------------------------!

  USE kind_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: inverse, adjoint, matrix_volume, triple_product
  PUBLIC :: cross_product, determinant
  CONTAINS

  FUNCTION inverse(A)

  !--------------------------------------------------------------------!
  ! inverse:  return the inverse of A(3,3).                            !
  !--------------------------------------------------------------------!

    REAL(q2), INTENT(IN), DIMENSION(3,3) :: A
    REAL(q2), DIMENSION(3,3) :: inverse
    REAL(q2) :: det

    det = determinant(A)
    IF (det == 0) STOP 'Divide by zero in matrix inverse'
    inverse = adjoint(A)/det

  RETURN
  END FUNCTION inverse

  FUNCTION adjoint(A)

  !--------------------------------------------------------------------!
  ! adjoint:  return the adjoint of A(3,3).                            !
  !--------------------------------------------------------------------!

    REAL(q2), INTENT(IN), DIMENSION(3,3) :: A
    REAL(q2), DIMENSION(3,3) :: adjoint

    adjoint(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
    adjoint(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    adjoint(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)

    adjoint(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    adjoint(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    adjoint(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)

    adjoint(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    adjoint(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    adjoint(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  RETURN
  END FUNCTION adjoint

  FUNCTION matrix_volume(h)

  !--------------------------------------------------------------------!
  ! matrix_volume:  function returning the triple product of the       !
  !                 lattice vectors.                                   !
  !--------------------------------------------------------------------!

    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: matrix_volume

    matrix_volume = h(1,1)*(h(2,2)*h(3,3) - h(2,3)*h(3,2))  &
    &             - h(1,2)*(h(2,1)*h(3,3) - h(3,1)*h(2,3))  &
    &             + h(1,3)*(h(2,1)*h(3,2) - h(3,1)*h(2,2))
    matrix_volume = ABS(matrix_volume)

  RETURN
  END FUNCTION matrix_volume
  
  FUNCTION triple_product(a,b,c)

  !--------------------------------------------------------------------!
  ! triple_product:  returns volume of parallelopiped defined by       !
  !                  vectors a, b & c.                                 !
  !--------------------------------------------------------------------!

    REAL(q2),INTENT(IN),DIMENSION(:) :: a,b,c
    REAL(q2) :: triple_product

    triple_product = c(1)*(a(2)*b(3) - a(3)*b(2))  & 
                   + c(2)*(a(3)*b(1) - a(1)*b(3))  &
                   + c(3)*(a(1)*b(2) - a(2)*b(1))

    RETURN
  END FUNCTION
 
  FUNCTION cross_product(A,B)

  !--------------------------------------------------------------------!
  ! cross_product:  calculate the cross product of two vectors.        !
  !--------------------------------------------------------------------!

    REAL(q2),INTENT(IN),DIMENSION(3) :: A,B
    REAL(q2), DIMENSION(3) :: cross_product

    cross_product(1) = A(2)*B(3) - A(3)*B(2)
    cross_product(2) = A(3)*B(1) - A(1)*B(3)
    cross_product(3) = A(1)*B(2) - A(2)*B(1)

    RETURN
  END FUNCTION

  FUNCTION determinant(A)

  !--------------------------------------------------------------------!
  ! determinant:  returns determinant of a 3x3 matrix.                 ! 
  !--------------------------------------------------------------------!

    REAL(q2),INTENT(IN),DIMENSION(3,3) :: A
    REAL(q2) :: determinant

    determinant = A(1,1)*A(2,2)*A(3,3) &
                - A(1,1)*A(2,3)*A(3,2) &
                - A(1,2)*A(2,1)*A(3,3) &
                + A(1,2)*A(2,3)*A(3,1) &
                + A(1,3)*A(2,1)*A(3,2) &
                - A(1,3)*A(2,2)*A(3,1)

    RETURN
    END FUNCTION

END MODULE matrix_mod

