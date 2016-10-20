!  procs.f03
!
!  Author: David N Alpert (david.alpert@gmail.com)
!  Advisor: Dr. Dong Qian (dong.qian@uc.edu)
!  University of Cincinnati, College of Engineering and Applied Science
!  School of Dynamic Systems, Mechanical Engineering Program
!  Computer Aided Engineering Research Laboratory
!
!  Record of Revisions
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   02/22/2010     D.N.A.     Original Code adapted from MATLAB and 1D FORTRAN code
!   11/04/2010     D.N.A.     Added explicit interfaces for Lapack subroutines
!   11/21/2011     D.N.A.     Add in CULA version of inverse and determinant for matrices of large order
!   02/22/2012     D.N.A.     Removal of cula_check_status subroutine which is now located in cula MODULE
!
!  All units should be converted to SI units.
!
!  The procs MODULE is to designed to contain general subroutines and functions, such as shape functions, Gauss quadrature, vector
!  norm, matrix trace, norm, condition number, determinant, and inverse.


MODULE PROCS

USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
USE interface_definitions ! Access explicit interface for LAPACK subroutines

IMPLICIT NONE

! Generic Interface for functions determinant, inverse, and norm

INTERFACE determinant
     MODULE PROCEDURE sdeterminant
     MODULE PROCEDURE ddeterminant
END INTERFACE determinant
INTERFACE inverse
     MODULE PROCEDURE sinverse
     MODULE PROCEDURE dinverse
END INTERFACE inverse
INTERFACE condition
     MODULE PROCEDURE scondition
     MODULE PROCEDURE dcondition
END INTERFACE
INTERFACE norm
     MODULE PROCEDURE snorm
     MODULE PROCEDURE dnorm
END INTERFACE norm
INTERFACE trace
     MODULE PROCEDURE strace
     MODULE PROCEDURE dtrace
END INTERFACE trace

CONTAINS

!---------------------------------------------------------------

FUNCTION strace(M)
!
! Calculate trace of square matrix M
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/28/2011     D.N.A.     Initial code
!

IMPLICIT NONE

REAL(kind=SGL) :: strace
REAL(kind=SGL), DIMENSION(:,:), INTENT(IN) :: M
INTEGER :: ii, N

N = SIZE(M,1)

IF( N /= SIZE(M,2) ) THEN
     WRITE(*,*) 'Error: Only a square matrix M may be used to calculate the trace(M).'
     STOP
END IF
strace = 0._SGL
IF (N == 0) THEN
     WRITE(*,*) 'Warning: Matrix M has size 0.  Function returning with TRACE = 0.'
ELSE
     DO ii = 1, N
          strace = strace + M(ii,ii)
     END DO
END IF
END FUNCTION strace

!---------------------------------------------------------------

FUNCTION dtrace(M)
!
! Calculate trace of square matrix M
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/28/2011     D.N.A.     Initial code
!

IMPLICIT NONE

REAL(kind=DBL) :: dtrace
REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) :: M
INTEGER :: ii, N

N = SIZE(M,1)

IF( N /= SIZE(M,2) ) THEN
     WRITE(*,*) 'Error: Only a square matrix M may be used to calculate the trace(M).'
     STOP
END IF
dtrace = 0._DBL
IF (N == 0) THEN
     WRITE(*,*) 'Warning: Matrix M has size 0.  Function returning with TRACE = 0.'
ELSE
     DO ii = 1, N
          dtrace = dtrace + M(ii,ii)
     END DO
END IF
END FUNCTION dtrace

!---------------------------------------------------------------

FUNCTION snorm(V, P)
!
! Calculate P-norm for real vector V. P defaults to 2.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/07/2011     D.N.A.     Initial code
!

IMPLICIT NONE

REAL(kind=SGL) :: snorm
INTEGER, INTENT(IN), OPTIONAL :: P
REAL(kind=SGL), DIMENSION(:), INTENT(IN) :: V
INTEGER :: ii, N, pp

IF(PRESENT(P)) THEN
     pp = P
ELSE
     pp = 2
END IF
N = SIZE(V,1)

IF (pp < 1) THEN
     WRITE(*,*) 'Error: Only a positive integer P may be used to find a P-norm of a vector.', pp, ' is an invalid input.'
     STOP
END IF
snorm = 0._SGL
IF (N == 0) THEN
     WRITE(*,*) 'Warning: Vector V has length 0.  Function returning with NORM = 0.'
ELSE IF (pp == 1) THEN
     snorm = SUM(ABS(V))
ELSE IF (pp == 2) THEN
     snorm = SQRT(DOT_PRODUCT(V,V))
ELSE
     DO ii = 1, N
          snorm = snorm + ABS(V(ii)) ** pp
     END DO
     snorm = snorm ** (1._SGL / REAL(pp, SGL))
END IF
END FUNCTION snorm

!---------------------------------------------------------------

FUNCTION dnorm(V, P)
!
! Calculate P-norm for real vector V. P defaults to 2.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/07/2011     D.N.A.     Initial code
!


IMPLICIT NONE

REAL(kind=DBL) :: dnorm
INTEGER, INTENT(IN), OPTIONAL :: P
REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: V
INTEGER :: ii, N, pp

IF(PRESENT(P)) THEN
     pp = P
ELSE
     pp = 2
END IF
N = SIZE(V,1)

IF (pp < 1) THEN
     WRITE(*,*) 'Error: Only a positive integer P may be used to find a P-norm of a vector.', pp, ' is an invalid input.'
     STOP
END IF
dnorm = 0._DBL
IF (N == 0) THEN
     WRITE(*,*) 'Warning: Vector V has length 0.  Function returning with NORM = 0.'
ELSE IF (pp == 1) THEN
     dnorm = SUM(ABS(V))
ELSE IF (pp == 2) THEN
     dnorm = SQRT(DOT_PRODUCT(V,V))
ELSE
     DO ii = 1, N
          dnorm = dnorm + ABS(V(ii)) ** pp
     END DO
     dnorm = dnorm ** (1._DBL / REAL(pp, DBL))
END IF

END FUNCTION dnorm

!---------------------------------------------------------------

FUNCTION sinverse( A, n )
!
! Computes inverse of a square matrix
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   02/23/2010     D.N.A.     Original Code
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! Order of matrix mat
REAL(kind=SGL), DIMENSION(n,n) :: sinverse ! Inverse of matrix mat
REAL(kind=SGL), DIMENSION(n,n), INTENT(IN) :: A ! Matrix to calculate inverse of
INTEGER :: errstat ! Error signal
INTEGER, DIMENSION(n) :: pivot ! Pivot vector
REAL(kind=SGL), DIMENSION(n*64) :: working ! Used for Lapack's GETRI subroutine
REAL(kind=SGL) :: det ! Determinant of matrix mat

det = sdeterminant(A)
IF(ABS(det) > 10._SGL * EPSILON(1._SGL)) THEN
     IF(n == 1) THEN
          sinverse(1,1) = A(1,1) ** (-1)
     ELSE IF(n == 2) THEN
          sinverse(1,1) =  A(2,2)
          sinverse(1,2) = -A(1,2)
          sinverse(2,1) = -A(2,1)
          sinverse(2,2) =  A(1,1)
          sinverse = sinverse / det
     ELSE IF(n == 3) THEN
          sinverse(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
          sinverse(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
          sinverse(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
          sinverse(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
          sinverse(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
          sinverse(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)
          sinverse(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
          sinverse(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
          sinverse(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)
          sinverse = sinverse / det
     ELSE
          sinverse = A
          CALL getrf( n, n, sinverse, n, pivot, errstat ) ! LU factor the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,A,I6)') 'There was an error calculating the LU factorization of the matrix.',&
               & ' The error code is ', errstat
          CALL getri( n, sinverse, n, pivot, working, n*64, errstat ) ! Calculate inverse of the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,A,I6)') 'There was an error calculating the inverse of the matrix. The error code is ',&
               & errstat
     END IF
ELSE IF (det == 0._SGL) THEN
     WRITE(*,*) 'Error: The matrix is singular, Its determinant is 0. The inverse cannot be calculated.'
     STOP
ELSE
     WRITE(*,*) 'Error: The matrix is close to singular.  Its determinant is ', det, '. A different approach should be sought.'
     STOP
END IF
END FUNCTION sinverse

!---------------------------------------------------------------

FUNCTION dinverse( A, n )
!
! Computes inverse of a square matrix
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   02/23/2010     D.N.A.     Original Code
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! Order of matrix mat
REAL(kind=DBL), DIMENSION(n,n) :: dinverse ! Inverse of matrix mat
REAL(kind=DBL), DIMENSION(n,n), INTENT(IN) :: A ! Matrix to calculate inverse of
INTEGER :: errstat ! Error signal
INTEGER, DIMENSION(n) :: pivot ! Pivot vector
REAL(kind=DBL), DIMENSION(n*64) :: working ! Used for Lapack's GETRI subroutine
REAL(kind=DBL) :: det ! Determinant of matrix mat

det = ddeterminant(A)
IF(ABS(det) > 10._DBL * EPSILON(1._DBL)) THEN
     IF(n == 1) THEN
          dinverse(1,1) = 1._SGL / A(1,1)
     ELSE IF(n == 2) THEN
          dinverse(1,1) =  A(2,2)
          dinverse(1,2) = -A(1,2)
          dinverse(2,1) = -A(2,1)
          dinverse(2,2) =  A(1,1)
          dinverse = dinverse / det
     ELSE IF(n == 3) THEN
          dinverse(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
          dinverse(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
          dinverse(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
          dinverse(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
          dinverse(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
          dinverse(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)
          dinverse(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
          dinverse(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
          dinverse(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)
          dinverse = dinverse / det
     ELSE
          dinverse = A
          CALL getrf( n, n, dinverse, n, pivot, errstat ) ! LU factor the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,A,I6)') 'There was an error calculating the LU factorization of the matrix.',&
               & ' The error code is ', errstat
          CALL getri( n, dinverse, n, pivot, working, n*64, errstat ) ! Calculate inverse of the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,A,I6)') 'There was an error calculating the inverse of the matrix. The error code is ',&
               & errstat
     END IF
ELSE IF (det == 0._DBL) THEN
     WRITE(*,*) 'Error: The matrix is singular, Its determinant is 0. The inverse cannot be calculated.'
     STOP
ELSE
     WRITE(*,*) 'Error: The matrix is close to singular.  Its determinant is ', det, '. A different approach should be sought.'
     STOP
END IF
END FUNCTION dinverse


!---------------------------------------------------------------

FUNCTION Scondition( A, n )
!
! Computes condition number of a square matrix using the equation: C(A) = norm(A*A^{-1})
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/18/2012     D.N.A.     Original Code
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! Order of matrix mat
REAL(kind=SGL) :: Scondition ! Condition number of matrix mat
REAL(kind=SGL), DIMENSION(n,n), INTENT(IN) :: A ! Matrix to calculate condition of
INTEGER :: ii
REAL(kind=SGL) :: normc
REAL(kind=SGL), DIMENSION(n, n) :: c

Scondition = 0._SGL
DO ii = 1, n
     Scondition = Scondition + DOT_PRODUCT(A(:,ii), A(:,ii))
END DO
Scondition = SQRT(Scondition)
c = Sinverse(A, n)
normc = 0._SGL
DO ii = 1, n
     normc = normc + DOT_PRODUCT(c(:,ii), c(:,ii))
END DO
Scondition = Scondition * SQRT(normc)
END FUNCTION Scondition

!---------------------------------------------------------------

FUNCTION Dcondition( A, n )
!
! Computes condition number of a square matrix using the equation: C(A) = norm(A*A^{-1})
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   01/18/2012     D.N.A.     Original Code
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! Order of matrix mat
REAL(kind=DBL) :: Dcondition ! Condition number of matrix mat
REAL(kind=DBL), DIMENSION(n,n), INTENT(IN) :: A ! Matrix to calculate condition of
INTEGER :: ii
REAL(kind=DBL) :: normc
REAL(kind=DBL), DIMENSION(n, n) :: c

Dcondition = 0._DBL
DO ii = 1, n
     Dcondition = Dcondition + DOT_PRODUCT(A(:,ii), A(:,ii))
END DO
Dcondition = SQRT(Dcondition)
c = Dinverse(A, n)
normc = 0._DBL
DO ii = 1, n
     normc = normc + DOT_PRODUCT(c(:,ii), c(:,ii))
END DO
Dcondition = Dcondition * SQRT(normc)
END FUNCTION Dcondition

!---------------------------------------------------------------


!---------------------------------------------------------------

FUNCTION sdeterminant(A)
!
! Computes determinant of a square matrix
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   02/23/2010     D.N.A.     Original Code
!   12/17/2010     D.N.A.     Removed n as input
!
IMPLICIT NONE

REAL(kind=SGL), DIMENSION(:,:), INTENT(IN) :: A ! Matrix to calculate determinant of
REAL(kind=SGL) :: sdeterminant ! Determinant of square matrix A
INTEGER :: errstat, ii, n ! Error signal, counter, Order of matrix mat
INTEGER, DIMENSION(:), ALLOCATABLE :: pivot ! Pivot vector
REAL(kind=SGL), DIMENSION(:,:), ALLOCATABLE :: mat
CHARACTER(len=120) :: errmesg ! Error message
 
n = SIZE(A,1) ! Order of matrix Af
IF( n == SIZE(A,2)) THEN ! Make sure matrix is square
     IF(n == 1) THEN
          sdeterminant = A(1,1)
     ELSE IF(n == 2) THEN
          sdeterminant = A(1,1) * A(2,2) - A(1,2) * A(2,1)
     ELSE IF(n == 3) THEN
          sdeterminant =   A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2)&
               &         - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)
     ELSE
          ALLOCATE(mat(n,n))
          !ALLOCATE(mat(n,n), pivot(n), STAT=errstat, ERRMSG = errmesg) ! Allocate memory
          !IF(errstat /= 0) WRITE(*,*) 'Error allocating Lapack inputs for GETRF in determinant. Code: ', errstat, &
          !     & ', Message: ', errmesg
          mat = A ! A cannot be changed because it's intent is IN
          CALL getrf( n, n, mat, n, pivot, errstat ) ! LU factor for the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,I6)') &
               &'There was an error calculating the LU factorization of the matrix.  The error code is ', errstat ! Error signal
          sdeterminant = 1.0_SGL ! Initialize determinant
          DO ii = 1, n ! Loop over matrix diagonal
               IF(pivot(ii) /= ii) THEN ! Conditional if permutation occured at diagonal
                    sdeterminant = -1._SGL * sdeterminant * mat(ii,ii) ! Update determinant
               ELSE ! Conditional if permuation did not occur at diagonal
                    sdeterminant = sdeterminant * mat(ii,ii) ! Update determinant
               END IF ! End conditional if permutation occured at diagonal
          END DO ! End loop over matrix diagonal
          DEALLOCATE(mat)
          !DEALLOCATE(mat, pivot, STAT=errstat, ERRMSG = errmesg) ! Deallocate memory
          !IF(errstat /= 0) WRITE(*,*) 'Error deallocating Lapack inputs for GETRF in determinant. Code: ', errstat, &
          !     & ', Message: ', errmesg
     END IF
ELSE
     WRITE(*,*) 'Error: Determinant function only accepts square matrices as input, your matrix is size ', n, ' x ', SIZE(A,2)
     STOP
END IF

END FUNCTION sdeterminant

!---------------------------------------------------------------

FUNCTION ddeterminant(A)
!
! Computes determinant of a square matrix
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   02/23/2010     D.N.A.     Original Code
!   12/17/2010     D.N.A.     Removed n as input
!
IMPLICIT NONE

REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) :: A ! Matrix to calculate determinant of
REAL(kind=DBL) :: ddeterminant ! Determinant of a square matrix A
INTEGER :: errstat, ii, n ! Error signal, counter, Order of matrix mat
INTEGER, DIMENSION(:), ALLOCATABLE :: pivot ! Pivot vector
REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: mat
CHARACTER(len=120) :: errmesg ! Error message

n = SIZE(A,1) ! Order of matrix A
IF( n == SIZE(A,2)) THEN ! Make sure matrix is square
     IF(n == 1) THEN
          ddeterminant = A(1,1)
     ELSE IF(n == 2) THEN
          ddeterminant = A(1,1) * A(2,2) - A(1,2) * A(2,1)
     ELSE IF(n == 3) THEN
          ddeterminant =   A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2)&
               &         - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) - A(1,3)*A(2,2)*A(3,1)
     ELSE
          ALLOCATE(mat(n,n))
          !ALLOCATE(mat(n,n), pivot(n), STAT=errstat, ERRMSG = errmesg) ! Allocate memory
          !IF(errstat /= 0) WRITE(*,*) 'Error allocating Lapack inputs for GETRF in determinant. Code: ', errstat, &
          !     & ', Message: ', errmesg
          mat = A ! A cannot be changed because it's intent is IN
          CALL getrf( n, n, mat, n, pivot, errstat ) ! LU factor for the Jacobian
          IF (errstat /=0 ) WRITE(*,'(1X,A,I6)') &
               &'There was an error calculating the LU factorization of the matrix.  The error code is ', errstat ! Error signal
          ddeterminant = 1.0_DBL ! Initialize determinant
          DO ii = 1, n ! Loop over matrix diagonal
               IF(pivot(ii) /= ii) THEN ! Conditional if permutation occured at diagonal
                    ddeterminant = -1._DBL * ddeterminant * mat(ii,ii) ! Update determinant
               ELSE ! Conditional if permuation did not occur at diagonal
                    ddeterminant = ddeterminant * mat(ii,ii) ! Update determinant
               END IF ! End conditional if permutation occured at diagonal
          END DO ! End loop over matrix diagonal
          DEALLOCATE(mat)
          !DEALLOCATE(mat, pivot, STAT=errstat, ERRMSG = errmesg) ! Deallocate memory
          !IF(errstat /= 0) WRITE(*,*) 'Error deallocating Lapack inputs for GETRF in determinant. Code: ', errstat, &
          !     & ', Message: ', errmesg
     END IF
ELSE
     WRITE(*,*) 'Error: Determinant function only accepts square matrices as input, your matrix is size ', n, ' x ', SIZE(A,2)
     STOP
END IF

END FUNCTION ddeterminant


!---------------------------------------------------------------

SUBROUTINE gaussquadrature( n, xi, wout )
!
! Calculate points xi and weights w for Gauss Legendre quadrature of order n. This subroutine is limited to orders of 1 through 20.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/02/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Changed from equations to values, Added KIND parameters
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n ! Order of quadrature
REAL(kind=REKIND), DIMENSION(n), INTENT(OUT) :: xi ! Gauss points
REAL(kind=REKIND), DIMENSION(n), OPTIONAL, INTENT(OUT) :: wout ! Gauss weights
REAL(kind=REKIND), DIMENSION(n) :: w
IF ( n == 1 ) THEN !Order of 1
     xi = [ 0._REKIND ] !Gauss point
     w  = [ 2._REKIND ] !Gauss weight
ELSE IF ( n == 2 ) THEN !Order of 2
     xi = [ -0.577350269189626_REKIND, 0.577350269189626_REKIND ] !Gauss points
     w  = [ 1._REKIND, 1._REKIND ] !Gauss weights
ELSE IF ( n == 3 ) THEN !Order of 3
     xi = [ -0.774596669241483_REKIND, 0._REKIND, 0.774596669241483_REKIND ] !Gauss points
     w  = [ 0.555555555555556_REKIND, 0.888888888888889_REKIND, 0.555555555555556_REKIND ] !Gauss weights
ELSE IF ( n == 4 ) THEN !Order of 4
     xi = [ - 0.861136311594053_REKIND, -0.339981043584856_REKIND, 0.339981043584856_REKIND,  0.861136311594053_REKIND ]
     w  = [ 0.347854845137454_REKIND, 0.652145154862546_REKIND, 0.652145154862546_REKIND, 0.347854845137454_REKIND ] !Gauss weights
ELSE IF ( n == 5 ) THEN !Order of 5
     xi = [ -0.906179845938664_REKIND, -0.538469310105683_REKIND, 0._REKIND, 0.538469310105683_REKIND, 0.906179845938664_REKIND ]
     w  = [ 0.236926885056189_REKIND, 0.478628670499366_REKIND, 0.568888888888889_REKIND, 0.478628670499366_REKIND, &
          & 0.236926885056189_REKIND ]
ELSE IF ( n == 6 ) THEN !Order of 6
     xi = [ -0.9324695142031521_REKIND, -0.6612093864662645_REKIND, -0.2386191860831969_REKIND, 0.2386191860831969_REKIND, &
          & 0.6612093864662645_REKIND, 0.9324695142031521_REKIND ]
     w  = [ 0.1713244923791704_REKIND, 0.3607615730481386_REKIND, 0.4679139345726910_REKIND, 0.4679139345726910_REKIND, &
          & 0.3607615730481386_REKIND, 0.1713244923791704_REKIND ]
ELSE IF ( n == 7 ) THEN !Order of 7
     xi = [ -0.9491079123427585_REKIND, -0.7415311855993945_REKIND, -0.4058451513773972_REKIND, 0._REKIND, &
          & 0.4058451513773972_REKIND, 0.7415311855993945_REKIND, 0.9491079123427585_REKIND ]
     w  = [ 0.1294849661688697_REKIND, 0.2797053914892766_REKIND, 0.3818300505051189_REKIND, 0.4179591836734694_REKIND, &
          & 0.3818300505051189_REKIND, 0.2797053914892766_REKIND, 0.1294849661688697_REKIND ]
ELSE IF ( n == 8 ) THEN !Order of 8
     xi = [ -0.9602898564975363_REKIND, -0.7966664774136267_REKIND, -0.5255324099163290_REKIND, -0.1834346424956498_REKIND, &
          & 0.1834346424956498_REKIND, 0.5255324099163290_REKIND, 0.7966664774136267_REKIND, 0.9602898564975363_REKIND ]
     w  = [ 0.1012285362903763_REKIND, 0.2223810344533745_REKIND, 0.3137066458778873_REKIND, 0.3626837833783620_REKIND, &
          & 0.3626837833783620_REKIND, 0.3137066458778873_REKIND, 0.2223810344533745_REKIND, 0.1012285362903763_REKIND ]
ELSE IF ( n == 9 ) THEN !Order of 9
     xi = [ -0.9681602395076261_REKIND, -0.8360311073266358_REKIND, -0.6133714327005904_REKIND, -0.3242534234038089_REKIND, &
          & 0._REKIND, 0.3242534234038089_REKIND, 0.6133714327005904_REKIND, 0.8360311073266358_REKIND, 0.9681602395076261_REKIND ]
     w  = [ 0.0812743883615744_REKIND, 0.1806481606948574_REKIND, 0.2606106964029354_REKIND, 0.3123470770400029_REKIND, &
          & 0.3302393550012598_REKIND, 0.3123470770400029_REKIND, 0.2606106964029354_REKIND, 0.1806481606948574_REKIND, &
          & 0.0812743883615744_REKIND ]
ELSE IF ( n == 10 ) THEN !Order of 10
     xi = [ -0.9739065285171717_REKIND, -0.8650633666889845_REKIND, -0.6794095682990244_REKIND, -0.4333953941292472_REKIND, &
          & -0.1488743389816312_REKIND, 0.1488743389816312_REKIND, 0.4333953941292472_REKIND, 0.6794095682990244_REKIND, &
          & 0.8650633666889845_REKIND, 0.9739065285171717_REKIND ]
     w  = [ 0.0666713443086881_REKIND, 0.1494513491505806_REKIND, 0.2190863625159820_REKIND, 0.2692667193099963_REKIND, &
          & 0.2955242247147529_REKIND, 0.2955242247147529_REKIND, 0.2692667193099963_REKIND, 0.2190863625159820_REKIND, &
          & 0.1494513491505806_REKIND, 0.0666713443086881_REKIND ]
ELSE IF ( n == 11 ) THEN !Order of 11
     xi = [ -0.9782286581460570_REKIND, -0.8870625997680953_REKIND, -0.7301520055740494_REKIND, -0.5190961292068118_REKIND, &
          & -0.2695431559523450_REKIND, 0._REKIND, 0.2695431559523450_REKIND, 0.5190961292068118_REKIND, 0.7301520055740494_REKIND,&
          & 0.8870625997680953_REKIND, 0.9782286581460570_REKIND ]
     w  = [ 0.0556685671161737_REKIND, 0.1255803694649046_REKIND, 0.1862902109277343_REKIND, 0.2331937645919905_REKIND, &
          & 0.2628045445102467_REKIND, 0.2729250867779006_REKIND, 0.2628045445102467_REKIND, 0.2331937645919905_REKIND, &
          & 0.1862902109277343_REKIND, 0.1255803694649046_REKIND, 0.0556685671161737_REKIND ]
ELSE IF ( n == 12 ) THEN !Order of 12
     xi = [ -0.9815606342467192_REKIND, -0.9041172563704749_REKIND, -0.7699026741943047_REKIND, -0.5873179542866175_REKIND, &
          & -0.3678314989981802_REKIND, -0.1252334085114689_REKIND, 0.1252334085114689_REKIND, 0.3678314989981802_REKIND, &
          & 0.5873179542866175_REKIND, 0.7699026741943047_REKIND, 0.9041172563704749_REKIND, 0.9815606342467192_REKIND ]
     w  = [ 0.0471753363865118_REKIND, 0.1069393259953184_REKIND, 0.1600783285433462_REKIND, 0.2031674267230659_REKIND, &
          & 0.2334925365383548_REKIND, 0.2491470458134028_REKIND, 0.2334925365383548_REKIND, 0.2031674267230659_REKIND, &
          & 0.1600783285433462_REKIND, 0.1069393259953184_REKIND, 0.0471753363865118_REKIND ]
ELSE IF ( n == 13 ) THEN !Order of 13
     xi = [ -0.9841830547185881_REKIND, -0.9175983992229779_REKIND, -0.8015780907333099_REKIND, -0.6423493394403402_REKIND, &
          & -0.4484927510364469_REKIND, -0.2304583159551348_REKIND, 0._REKIND, 0.2304583159551348_REKIND, &
          & 0.4484927510364469_REKIND, 0.6423493394403402_REKIND, 0.8015780907333099_REKIND, 0.9175983992229779_REKIND, &
          & 0.9841830547185881_REKIND ]
     w  = [ 0.0404840047653159_REKIND, 0.0921214998377285_REKIND, 0.1388735102197872_REKIND, 0.1781459807619457_REKIND, &
          & 0.2078160475368885_REKIND, 0.2262831802628972_REKIND, 0.2325515532308739_REKIND, 0.2262831802628972_REKIND, &
          & 0.2078160475368885_REKIND, 0.1781459807619457_REKIND, 0.1388735102197872_REKIND, 0.0921214998377285_REKIND, &
          & 0.0404840047653159_REKIND ]
ELSE IF ( n == 14 ) THEN !Order of 14
     xi = [ -0.9862838086968123_REKIND, -0.9284348836635735_REKIND, -0.8272013150697650_REKIND, -0.6872929048116855_REKIND, &
          & -0.5152486363581541_REKIND, -0.3191123689278897_REKIND, -0.1080549487073437_REKIND, 0.1080549487073437_REKIND, &
          & 0.3191123689278897_REKIND, 0.5152486363581541_REKIND, 0.6872929048116855_REKIND, 0.8272013150697650_REKIND, &
          & 0.9284348836635735_REKIND, 0.9862838086968123_REKIND ]
     w  = [ 0.0351194603317519_REKIND, 0.0801580871597602_REKIND, 0.1215185706879032_REKIND, 0.1572031671581935_REKIND, &
          & 0.1855383974779378_REKIND, 0.2051984637212956_REKIND, 0.2152638534631578_REKIND, 0.2152638534631578_REKIND, &
          & 0.2051984637212956_REKIND, 0.1855383974779378_REKIND, 0.1572031671581935_REKIND, 0.1215185706879032_REKIND, &
          & 0.0801580871597602_REKIND, 0.0351194603317519_REKIND ]
ELSE IF ( n == 15 ) THEN !Order of 15
     xi = [ -0.9879925180204854_REKIND, -0.9372733924007060_REKIND, -0.8482065834104272_REKIND, -0.7244177313601701_REKIND, &
          & -0.5709721726085388_REKIND, -0.3941513470775634_REKIND, -0.2011940939974345_REKIND, 0._REKIND, &
          & 0.2011940939974345_REKIND, 0.3941513470775634_REKIND, 0.5709721726085388_REKIND, 0.7244177313601701_REKIND, &
          & 0.8482065834104272_REKIND, 0.9372733924007060_REKIND, 0.9879925180204854_REKIND ]
     w  = [ 0.0307532419961173_REKIND, 0.0703660474881081_REKIND, 0.1071592204671719_REKIND, 0.1395706779261543_REKIND, &
          & 0.1662692058169939_REKIND, 0.1861610000155622_REKIND, 0.1984314853271116_REKIND, 0.2025782419255613_REKIND, &
          & 0.1984314853271116_REKIND, 0.1861610000155622_REKIND, 0.1662692058169939_REKIND, 0.1395706779261543_REKIND, &
          & 0.1071592204671719_REKIND, 0.0703660474881081_REKIND, 0.0307532419961173_REKIND ]
ELSE IF ( n == 16 ) THEN !Order of 16
     xi = [ -0.9894009349916499_REKIND, -0.9445750230732326_REKIND, -0.8656312023878318_REKIND, -0.7554044083550030_REKIND, &
          & -0.6178762444026438_REKIND, -0.4580167776572274_REKIND, -0.2816035507792589_REKIND, -0.0950125098376374_REKIND, &
          & 0.0950125098376374_REKIND, 0.2816035507792589_REKIND, 0.4580167776572274_REKIND, 0.6178762444026438_REKIND, &
          & 0.7554044083550030_REKIND, 0.8656312023878318_REKIND, 0.9445750230732326_REKIND, 0.9894009349916499_REKIND ]
     w  = [ 0.0271524594117541_REKIND, 0.0622535239386479_REKIND, 0.0951585116824928_REKIND, 0.1246289712555339_REKIND, &
          & 0.1495959888165767_REKIND, 0.1691565193950025_REKIND, 0.1826034150449236_REKIND, 0.1894506104550685_REKIND, &
          & 0.1894506104550685_REKIND, 0.1826034150449236_REKIND, 0.1691565193950025_REKIND, 0.1495959888165767_REKIND, &
          & 0.1246289712555339_REKIND, 0.0951585116824928_REKIND, 0.0622535239386479_REKIND, 0.0271524594117541_REKIND ]
ELSE IF ( n == 17 ) THEN !Order of 17
     xi = [ -0.9905754753144174_REKIND, -0.9506755217687678_REKIND, -0.8802391537269859_REKIND, -0.7815140038968014_REKIND, &
          & -0.6576711592166907_REKIND, -0.5126905370864769_REKIND, -0.3512317634538763_REKIND, -0.1784841814958479_REKIND, &
          & 0._REKIND, 0.1784841814958479_REKIND, 0.3512317634538763_REKIND, 0.5126905370864769_REKIND, 0.6576711592166907_REKIND, &
          & 0.7815140038968014_REKIND, 0.8802391537269859_REKIND, 0.9506755217687678_REKIND, 0.9905754753144174_REKIND ]
     w  = [ 0.0241483028685479_REKIND, 0.0554595293739872_REKIND, 0.0850361483171792_REKIND, 0.1118838471934040_REKIND, &
          & 0.1351363684685255_REKIND, 0.1540457610768103_REKIND, 0.1680041021564500_REKIND, 0.1765627053669926_REKIND, &
          & 0.1794464703562065_REKIND, 0.1765627053669926_REKIND, 0.1680041021564500_REKIND, 0.1540457610768103_REKIND, &
          & 0.1351363684685255_REKIND, 0.1118838471934040_REKIND, 0.0850361483171792_REKIND, 0.0554595293739872_REKIND, &
          & 0.0241483028685479_REKIND ]
ELSE IF ( n == 18 ) THEN !Order of 18
     xi = [ -0.9915651684209309_REKIND, -0.9558239495713977_REKIND, -0.8926024664975557_REKIND, -0.8037049589725231_REKIND, &
          & -0.6916870430603532_REKIND, -0.5597708310739475_REKIND, -0.4117511614628426_REKIND, -0.2518862256915055_REKIND, &
          & -0.0847750130417353_REKIND, 0.0847750130417353_REKIND, 0.2518862256915055_REKIND, 0.4117511614628426_REKIND, &
          & 0.5597708310739475_REKIND, 0.6916870430603532_REKIND, 0.8037049589725231_REKIND, 0.8926024664975557_REKIND, &
          & 0.9558239495713977_REKIND, 0.9915651684209309_REKIND ]
     w  = [ 0.0216160135264833_REKIND, 0.0497145488949698_REKIND, 0.0764257302548891_REKIND, 0.1009420441062872_REKIND, &
          & 0.1225552067114785_REKIND, 0.1406429146706507_REKIND, 0.1546846751262652_REKIND, 0.1642764837458327_REKIND, &
          & 0.1691423829631436_REKIND, 0.1691423829631436_REKIND, 0.1642764837458327_REKIND, 0.1546846751262652_REKIND, &
          & 0.1406429146706507_REKIND, 0.1225552067114785_REKIND, 0.1009420441062872_REKIND, 0.0764257302548891_REKIND, &
          & 0.0497145488949698_REKIND, 0.0216160135264833_REKIND ]
ELSE IF ( n == 19 ) THEN !Order of 19
     xi = [ -0.9924068438435844_REKIND, -0.9602081521348300_REKIND, -0.9031559036148179_REKIND, -0.8227146565371428_REKIND, &
          & -0.7209661773352294_REKIND, -0.6005453046616810_REKIND, -0.4645707413759609_REKIND, -0.3165640999636298_REKIND, &
          & -0.1603586456402254_REKIND, 0._REKIND, 0.1603586456402254_REKIND, 0.3165640999636298_REKIND, 0.4645707413759609_REKIND,&
          & 0.6005453046616810_REKIND, 0.7209661773352294_REKIND, 0.8227146565371428_REKIND, 0.9031559036148179_REKIND, &
          & 0.9602081521348300_REKIND, 0.9924068438435844_REKIND ]
     w  = [ 0.0194617882297265_REKIND, 0.0448142267656996_REKIND, 0.0690445427376412_REKIND, 0.0914900216224500_REKIND, &
          & 0.1115666455473340_REKIND, 0.1287539625393362_REKIND, 0.1426067021736066_REKIND, 0.1527660420658597_REKIND, &
          & 0.1589688433939543_REKIND, 0.1610544498487837_REKIND, 0.1589688433939543_REKIND, 0.1527660420658597_REKIND, &
          & 0.1426067021736066_REKIND, 0.1287539625393362_REKIND, 0.1115666455473340_REKIND, 0.0914900216224500_REKIND, &
          & 0.0690445427376412_REKIND, 0.0448142267656996_REKIND, 0.0194617882297265_REKIND ]
ELSE IF ( n == 20 ) THEN !Order of 20
     xi = [ -0.9931285991850949_REKIND, -0.9639719272779138_REKIND, -0.9639719272779138_REKIND, -0.9122344282513259_REKIND, &
          & -0.8391169718222188_REKIND, -0.7463319064601508_REKIND, -0.6360536807265150_REKIND, -0.5108670019508271_REKIND, &
          & -0.3737060887154195_REKIND, -0.2277858511416451_REKIND, -0.0765265211334973_REKIND, 0.0765265211334973_REKIND, &
          & 0.2277858511416451_REKIND, 0.3737060887154195_REKIND, 0.5108670019508271_REKIND, 0.6360536807265150_REKIND, &
          & 0.7463319064601508_REKIND, 0.8391169718222188_REKIND, 0.9122344282513259_REKIND, 0.9639719272779138_REKIND, &
          & 0.9639719272779138_REKIND, 0.9931285991850949_REKIND ]
     w  = [ 0.0176140071391521_REKIND, 0.0406014298003869_REKIND, 0.0626720483341091_REKIND, 0.0832767415767048_REKIND, &
          & 0.1019301198172404_REKIND, 0.1181945319615184_REKIND, 0.1316886384491766_REKIND, 0.1420961093183820_REKIND, &
          & 0.1491729864726037_REKIND, 0.1527533871307258_REKIND, 0.1527533871307258_REKIND, 0.1491729864726037_REKIND, &
          & 0.1420961093183820_REKIND, 0.1316886384491766_REKIND, 0.1181945319615184_REKIND, 0.1019301198172404_REKIND, &
          & 0.0832767415767048_REKIND, 0.0626720483341091_REKIND, 0.0406014298003869_REKIND, 0.0176140071391521_REKIND ]
ELSE !Order out of range
     WRITE (*,*) 'The Gauss quadrature order is limited to integers between of 1 and 20' !Error message
END IF !End conditional
wexist: IF (PRESENT(wout)) THEN !Conditional if Gauss weight output is present
     wout=w !Copy Gauss weight
END IF wexist !End conditional if w output is present

END SUBROUTINE gaussquadrature

!---------------------------------------------------------------

SUBROUTINE shape1x(ns, x, n)
!
! Calculate array containing the full Lagrangian shape function (shapefunc) at point x of order n.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/15/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Added KIND parameters
!   01/17/2011     D.N.A.     Added higher orders and fixed order 5 elements 2&3
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! Point to evaluate shape function at
REAL(kind=REKIND), DIMENSION(n), INTENT(OUT) :: ns ! Shape function array of order n at point x
INTEGER :: ii, jj
REAL(kind=REKIND), DIMENSION(n) :: iso

order: IF ( n == 2 ) THEN ! Order 2
     ns = [0.5_REKIND - 0.5_REKIND * x, 0.5_REKIND + 0.5_REKIND * x ] ! [ N1 N2 ]
ELSE IF ( n == 3 ) THEN order ! Order 3
     ns = [ 0.5_REKIND * (-x + x**2), 1._REKIND - x**2, 0.5_REKIND * (x + x**2) ] ! [N1 N2 N3]
ELSE IF ( n == 4 ) THEN order ! Order 4
     ns( 1 ) = 0.0625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) ! N1
     ns( 2 ) = 0.5625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) ! N2
     ns( 3 ) = 0.5625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) ! N3
     ns( 4 ) = 0.0625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) ! N4
ELSE IF ( n == 5 ) THEN order ! Order 5
     ns(1) = (x - x ** 2 - 4._REKIND * x ** 3 + 4._REKIND * x ** 4 ) / 6._REKIND ! N1
     ns(2) = 4._REKIND * ( -x + 2._REKIND * x ** 2 + x ** 3 - 2._REKIND * x ** 4 ) / 3._REKIND ! N2
     ns(3) = 1._REKIND - 5._REKIND * x ** 2 + 4._REKIND * x ** 4 ! N3
     ns(4) = 4._REKIND * ( x + 2._REKIND * x ** 2 - x ** 3 - 2._REKIND * x ** 4 ) / 3._REKIND ! N4
     ns(5) = (-x - x ** 2 + 4._REKIND * x ** 3 + 4._REKIND * x ** 4 ) / 6._REKIND ! N5
ELSE IF ( n > 5 ) THEN order ! Order > 5
     iso(1) = -1._REKIND - 2._REKIND / REAL(n-1, REKIND)
     DO ii = 2, n ! Loop over nodes
          iso(ii) = iso(ii-1) + 2._REKIND / REAL(n-1, REKIND) ! Isoparametric nodal coordinate, dimensionless
     END DO ! End loop over nodes
     ns = 1._REKIND ! Initialize N
     DO ii = 1, n ! Loop over nodes
          DO jj = 1, n ! Loop over nodes
               IF( jj /= ii) THEN ! Check if denominator is 0
                    ns(ii) = ns(ii) * (iso(jj) - x) / (iso(jj) - iso(ii)) ! N
               END IF ! End check if denominator is 0
          END DO ! End loop over nodes
     END DO ! End loop over nodes
END IF order ! End conditional
END SUBROUTINE shape1x

!---------------------------------------------------------------

SUBROUTINE shape1dx(dndx, x, n)
!
! Calculate array containing the first derivative of the full Lagrangian shape function (dndx) at x of order n
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/02/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Added KIND parameters
!   01/17/2011     D.N.A.     Added higher orders and fixed order 5 elements 2&3
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x ! Point to evaluate shape function derivative at
REAL(kind=REKIND), DIMENSION(n), INTENT(OUT) :: dndx ! Shape function derivative array of order n at point x
INTEGER :: ii, jj
REAL(kind=REKIND) :: eps
REAL(kind=REKIND), DIMENSION(n) :: iso, dndxm

IF ( n == 2 ) THEN ! Order 2
     dndx = [-0.5_REKIND, 0.5_REKIND ] ! [ dN1/dx dN2/x ]
ELSE IF ( n == 3 ) THEN ! Order 3
     dndx = [ -0.5_REKIND + x, -x - x, 0.5_REKIND + x ] ! [ dN1/dx dN2/dx dN3/dx ]
ELSE IF ( n == 4 ) THEN ! Order 4
     dndx( 1 ) = 0.0625_REKIND * (  1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) ! dN1/dx
     dndx( 2 ) = 0.5625_REKIND * ( -3._REKIND -  2._REKIND * x +  9._REKIND * x ** 2 ) ! dN2/dx
     dndx( 3 ) = 0.5625_REKIND * (  3._REKIND -  2._REKIND * x -  9._REKIND * x ** 2 ) ! dN3/dx
     dndx( 4 ) = 0.0625_REKIND * ( -1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) ! dN4/dx
ELSE IF ( n == 5 ) THEN  ! Order 5
     dndx(1) = (1._REKIND - 2._REKIND * x - 12._REKIND * x ** 2 + 16._REKIND * x ** 3 ) / 6._REKIND ! dN1/dx
     dndx(2) = 4._REKIND * ( -1._REKIND + 4._REKIND * x + 3._REKIND * x ** 2 - 8._REKIND * x ** 3 ) / 3._REKIND ! dN2/dx
     dndx(3) = -10._REKIND * x + 16._REKIND * x ** 3 ! dN3/dx
     dndx(4) = 4._REKIND * ( 1._REKIND + 4._REKIND * x - 3._REKIND * x ** 2 - 8._REKIND * x ** 3 ) / 3._REKIND ! dN4/dx
     dndx(5) = (-1._REKIND - 2._REKIND * x + 12._REKIND * x ** 2 + 16._REKIND * x ** 3 ) / 6._REKIND ! dN5/dx
ELSE IF( n > 5 ) THEN ! Order > 5
     eps = SQRT(EPSILON(1._REKIND)) ! Differential
     iso(1) = -1._REKIND - 2._REKIND / REAL(n-1, REKIND)
     DO ii = 2, n ! Loop over nodes
          iso(ii) = iso(ii-1) + 2._REKIND / REAL(n-1, REKIND) ! Isoparametric nodal coordinate, dimensionless
     END DO ! End loop over nodes
     dndxm = 1._REKIND ! Initialize dNm
     dndx  = 1._REKIND ! Initialize dN
     DO ii = 1, n ! Loop over nodes
          DO jj = 1, n ! Loop over nodes
               IF( jj /= ii) THEN ! Check if denominator is 0
                    dndxm(ii) = dndxm(ii) * (iso(jj) - (x-eps)) / (iso(jj) - iso(ii)) ! dNm
                    dndx(ii)  = dndx(ii)  * (iso(jj) - (x+eps)) / (iso(jj) - iso(ii)) ! dN
               END IF ! End check if denominator is 0
          END DO ! End loop over nodes
     END DO ! End loop over nodes
     dndx = (dndx - dndxm) / (eps + eps)
END IF ! End conditional
END SUBROUTINE shape1dx

!---------------------------------------------------------------

SUBROUTINE shape1ddx(ddndxdx, x, n)
!
! Calculates array containing the second derivative of the full Lagrangian shape function (ddndxdx) at x of order n
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/02/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Added KIND parameters
!   01/17/2011     D.N.A.     Added higher orders and fixed order 5 elements 2&3
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL(kind=REKIND), INTENT(IN) :: x
REAL(kind=REKIND), DIMENSION(n), INTENT(OUT) :: ddndxdx
REAL(kind=REKIND) :: eps
REAL(kind=REKIND), DIMENSION(n) :: ddndxdxm

IF ( n == 2 ) THEN
     ddndxdx = [ 0._REKIND, 0._REKIND ]
ELSE IF ( n == 3 ) THEN
     ddndxdx = [ 1._REKIND, -2._REKIND, 1._REKIND ]
ELSE IF ( n == 4 ) THEN
     ddndxdx( 1 ) =  1.125_REKIND -  3.375_REKIND * x
     ddndxdx( 2 ) = -1.125_REKIND + 10.125_REKIND * x
     ddndxdx( 3 ) = -1.125_REKIND - 10.125_REKIND * x
     ddndxdx( 4 ) =  1.125_REKIND +  3.375_REKIND * x
ELSE IF ( n == 5 ) THEN
     ddndxdx(1) = (-1._REKIND  - 12._REKIND * x + 24._REKIND * x ** 2 ) / 3._REKIND
     ddndxdx(2) = ( 16._REKIND + 24._REKIND * x - 96._REKIND * x ** 2 ) / 3._REKIND
     ddndxdx(3) = -10._REKIND + 48._REKIND * x ** 2
     ddndxdx(4) = ( 16._REKIND - 24._REKIND * x - 96._REKIND * x ** 2 ) / 3._REKIND
     ddndxdx(5) = (-1._REKIND  + 12._REKIND * x + 24._REKIND * x ** 2 ) / 3._REKIND
ELSE IF ( n > 5) THEN
     eps = 10._REKIND * SQRT(EPSILON(1._REKIND))
     CALL shape1dx(ddndxdxm, x-eps, n)
     CALL shape1dx(ddndxdx,  x+eps, n)
     ddndxdx = (ddndxdx - ddndxdxm) / (eps + eps)
END IF
END SUBROUTINE shape1ddx

!---------------------------------------------------------------

SUBROUTINE shape2x(nout, x, y, n)
!
! Calculates array containing the full Lagrangian shape function (nout) at coordinate (x,y) of order n. This subroutine is
! limited to orders of 2 - 4 (for Q4, Q9, and Q16 elements).
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/15/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Added KIND parameters
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! X-coordinate of point to evaluate shape function at
REAL(kind=REKIND), INTENT(IN) :: y  ! Y-coordinate of point to evaluate shape function at
REAL(kind=REKIND), DIMENSION(n ** 2) :: ns ! Shape function array of order n at point (x,y)
REAL(kind=REKIND), INTENT(OUT), DIMENSION(2,2*(n**2)) :: nout ! 2D Shape functions
INTEGER :: ii, jj ! Counter

order: IF ( n == 2 ) THEN ! Order 2
     ! Shape:
     ! 4 -- 3
     ! |    |
     ! 1 -- 2
     ns( 1 ) = 0.25_REKIND * ( 1._REKIND - x ) * ( 1._REKIND - y ) ! N1
     ns( 2 ) = 0.25_REKIND * ( x + 1._REKIND ) * ( 1._REKIND - y ) ! N2
     ns( 3 ) = 0.25_REKIND * ( x + 1._REKIND ) * ( y + 1._REKIND ) ! N3
     ns( 4 ) = 0.25_REKIND * ( 1._REKIND - x ) * ( y + 1._REKIND ) ! N4
ELSE IF ( n == 3 ) THEN order ! Order 3
     ! Shape:
     ! 4 - 7 - 3
     ! |   |   |
     ! 8 - 9 - 6
     ! |   |   |
     ! 1 - 5 - 2
     ns( 1 ) = 0.25_REKIND * ( -x + x ** 2 ) * ( -y + y ** 2 ) ! N1
     ns( 2 ) = 0.25_REKIND * (  x + x ** 2 ) * ( -y + y ** 2 ) !N2
     ns( 3 ) = 0.25_REKIND * (  x + x ** 2 ) * (  y + y ** 2 ) !N3
     ns( 4 ) = 0.25_REKIND * ( -x + x ** 2 ) * (  y + y ** 2 ) ! N4
     ns( 5 ) = 0.5_REKIND * (1._REKIND - x ** 2) * ( -y + y ** 2 ) ! N5
     ns( 6 ) = 0.5_REKIND * ( x + x ** 2 ) * (1._REKIND - y ** 2) !N6
     ns( 7 ) = 0.5_REKIND * (1._REKIND - x ** 2) * ( y + y ** 2 ) ! N7
     ns( 8 ) = 0.5_REKIND * ( -x + x ** 2 ) * (1._REKIND - y ** 2) ! N8
     ns( 9 ) = (1._REKIND - x ** 2) * (1._REKIND - y ** 2) ! N9
ELSE IF (n == 4 ) THEN order ! Order 4
     ! Shape:
     ! 4  - 10 - 9  - 3
     ! |    |    |    |
     ! 11 - 16 - 15 - 8
     ! |    |    |    |
     ! 12 - 13 - 14 - 7
     ! |    |    |    |f
     ! 1  - 5  - 6  - 2
     ns(  1 ) = 0.00390625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! N1
     ns(  2 ) = 0.00390625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! N2
     ns(  3 ) = 0.00390625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! N3
     ns(  4 ) = 0.00390625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! N4
     ns(  5 ) = 0.03515625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! N5
     ns(  6 ) = 0.03515625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! N6
     ns(  7 ) = 0.03515625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! N7
     ns(  8 ) = 0.03515625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! N8
     ns(  9 ) = 0.03515625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! N9
     ns( 10 ) = 0.03515625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! N10
     ns( 11 ) = 0.03515625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! N11
     ns( 12 ) = 0.03515625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! N12
     ns( 13 ) = 0.31640625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! N13
     ns( 14 ) = 0.31640625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! N14
     ns( 15 ) = 0.31640625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! N15
     ns( 16 ) = 0.31640625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! N16
END IF order ! End Conditional
nout = 0._REKIND
jj = 0
nodel: DO ii = 1, n ** 2 ! Loop over spatial nodes
     jj = jj + 2
     nout(1, jj-1) = ns(ii) ! x-direction shape functions
     nout(2, jj)   = ns(ii) ! y-direction shape functions
END DO nodel ! End loop over spatial nodes

END SUBROUTINE shape2x

!---------------------------------------------------------------

SUBROUTINE shape2dx(dns, dnxy, x, y, n)
!
! Calculates array containing the first derivatives of the full Lagrangian shape function (dns and dnxy) at coordinate (x,y) of
! order n. This subroutine is limited to orders of 2 - 4 (for Q4, Q9, and Q16 elements). The outputs are used to calculate the 2D
! Jacobian and the strain displacement matrix.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   05/02/2009     D.N.A.     Original Code
!   07/07/2009     D.N.A.     Added KIND parameters
!

IMPLICIT NONE

INTEGER :: ii, jj, kk ! Counter
INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! X-coordinate of point to evaluate shape function at
REAL(kind=REKIND), INTENT(IN) :: y  ! Y-coordinate of point to evaluate shape function at
REAL(kind=REKIND), INTENT(OUT), DIMENSION(2,n**2) :: dns ! 2D Shape functions
REAL(kind=REKIND), INTENT(OUT), DIMENSION(4,2*(n**2)) :: dnxy ! 2D Shape functions
order: IF ( n == 2 ) THEN ! Order 2
     ! Shape:
     ! 4 -- 3
     ! |    |
     ! 1 -- 2
     dns( 1, 1 ) = -0.25_REKIND + 0.25_REKIND * y ! dN1dx
     dns( 1, 2 ) =  0.25_REKIND - 0.25_REKIND * y ! dN2dx
     dns( 1, 3 ) =  0.25_REKIND + 0.25_REKIND * y ! dN3dx
     dns( 1, 4 ) = -0.25_REKIND - 0.25_REKIND * y ! dN4dx
     dns( 2, 1 ) = -0.25_REKIND + 0.25_REKIND * x ! dN1dy
     dns( 2, 2 ) = -0.25_REKIND - 0.25_REKIND * x ! dN2dy
     dns( 2, 3 ) =  0.25_REKIND + 0.25_REKIND * x ! dN3dy
     dns( 2, 4 ) =  0.25_REKIND - 0.25_REKIND * x ! dN4dy
ELSE IF ( n == 3 ) THEN order ! Order 3
     ! Shape:
     ! 4 - 7 - 3
     ! |   |   |
     ! 8 - 9 - 6
     ! |   |   |
     ! 1 - 5 - 2
     dns( 1, 1 ) = 0.25_REKIND * ( -1._REKIND + x + x ) * ( -y + y ** 2 ) ! N1
     dns( 1, 2 ) = 0.25_REKIND * ( 1._REKIND + x + x ) * ( -y + y ** 2 ) !N2
     dns( 1, 3 ) = 0.25_REKIND * ( 1._REKIND + x + x ) * (  y + y ** 2 ) !N3
     dns( 1, 4 ) = 0.25_REKIND * ( -1._REKIND + x + x ) * (  y + y ** 2 ) ! N4
     dns( 1, 5 ) = 0.5_REKIND * (-2._REKIND * x ) * ( -y + y ** 2 ) ! N5
     dns( 1, 6 ) = 0.5_REKIND * ( 1._REKIND + x + x ) * (1._REKIND - y ** 2) !N6
     dns( 1, 7 ) = 0.5_REKIND * (-2._REKIND * x ) * ( y + y ** 2 ) ! N7
     dns( 1, 8 ) = 0.5_REKIND * ( -1._REKIND + x + x ) * (1._REKIND - y ** 2) ! N8
     dns( 1, 9 ) = -x - x + x*y*y + x*y*y ! N9
     dns( 2, 1 ) = 0.25_REKIND * ( -x + x ** 2 ) * ( -1._REKIND + y + y ) ! N1
     dns( 2, 2 ) = 0.25_REKIND * (  x + x ** 2 ) * ( -1._REKIND + y + y ) ! N2
     dns( 2, 3 ) = 0.25_REKIND * (  x + x ** 2 ) * (  1._REKIND + y + y ) ! N3
     dns( 2, 4 ) = 0.25_REKIND * ( -x + x ** 2 ) * (  1._REKIND + y + y ) ! N4
     dns( 2, 5 ) = 0.5_REKIND * (1._REKIND - x ** 2) * ( -1._REKIND + y + y ) ! N5
     dns( 2, 6 ) = -x*y - y*x*x !N6
     dns( 2, 7 ) = 0.5_REKIND * (1._REKIND - x ** 2) * (  1._REKIND + y + y ) ! N7
     dns( 2, 8 ) = x*y - y*x*x ! N8
     dns( 2, 9 ) = -y - y + y*x*x + y*x*x ! N9
ELSE IF (n == 4 ) THEN order ! Order 4
     ! Shape:
     ! 4  - 10 - 9  - 3
     ! |    |    |    |
     ! 11 - 16 - 15 - 8
     ! |    |    |    |
     ! 12 - 13 - 14 - 7
     ! |    |    |    |
     ! 1  - 5  - 6  - 2
     dns( 1,  1 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! dN1dx
     dns( 1,  2 ) = 0.00390625_REKIND * ( -1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! dN2dx
     dns( 1,  3 ) = 0.00390625_REKIND * ( -1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! dN3dx
     dns( 1,  4 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! dN4dx
     dns( 1,  5 ) = 0.03515625_REKIND * ( -3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! dN5dx
     dns( 1,  6 ) = 0.03515625_REKIND * ( 3._REKIND  - 2._REKIND * x - 9._REKIND * x ** 2 ) &
               & * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 ) ! dN6dx
     dns( 1,  7 ) = 0.03515625_REKIND * ( -1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! dN7dx
     dns( 1,  8 ) = 0.03515625_REKIND * ( -1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! dN8dx
     dns( 1,  9 ) = 0.03515625_REKIND * ( 3._REKIND  - 2._REKIND * x - 9._REKIND * x ** 2 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! dN9dx
     dns( 1, 10 ) = 0.03515625_REKIND * ( -3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
               & * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 ) ! dN10dx
     dns( 1, 11 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! dN11dx
     dns( 1, 12 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! dN12dx
     dns( 1, 13 ) = 0.31640625_REKIND * ( -3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! dN13dx
     dns( 1, 14 ) = 0.31640625_REKIND * ( 3._REKIND  - 2._REKIND * x - 9._REKIND * x ** 2 ) &
               & * (  1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 ) ! dN14dx
     dns( 1, 15 ) = 0.31640625_REKIND * ( 3._REKIND  - 2._REKIND * x - 9._REKIND * x ** 2 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! dN15dx
     dns( 1, 16 ) = 0.31640625_REKIND * ( -3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
               & * (  1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 ) ! dN16dx
     dns( 2,  1 ) = 0.00390625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2 ) ! dN1dy
     dns( 2,  2 ) = 0.00390625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2 ) ! dN2dy
     dns( 2,  3 ) = 0.00390625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2 ) ! dN3dy
     dns( 2,  4 ) = 0.00390625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( -1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2 ) ! dN4dy
     dns( 2,  5 ) = 0.03515625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2 ) ! dN5dy
     dns( 2,  6 ) = 0.03515625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2 ) ! dN6dy
     dns( 2,  7 ) = 0.03515625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( -3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2 ) ! dN7dy
     dns( 2,  8 ) = 0.03515625_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) &
               & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2 ) ! dN8dy
     dns( 2,  9 ) = 0.03515625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2 ) ! dN9dy
     dns( 2, 10 ) = 0.03515625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( -1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2 ) ! dN10dy
     dns( 2, 11 ) = 0.03515625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2 ) ! dN11dy
     dns( 2, 12 ) = 0.03515625_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) &
               & * ( -3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2 ) ! dN12dy
     dns( 2, 13 ) = 0.31640625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( -3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2 ) ! dN13dy
     dns( 2, 14 ) = 0.31640625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( -3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2 ) ! dN14dy
     dns( 2, 15 ) = 0.31640625_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) &
               & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2 ) ! dN15dy
     dns( 2, 16 ) = 0.31640625_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) &
               & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2 ) ! dN16dy
END IF order ! End conditional

! Initialize x/y direction shape functions
dnxy = 0._REKIND
jj = 0
kk = -1
nodel: DO ii = 1, n ** 2 ! Loop over nodes
     jj = jj + 2
     kk = kk + 2
     dnxy(1, kk) = dns(1, ii)
     dnxy(2, kk) = dns(2, ii)
     dnxy(3, jj) = dns(1, ii)
     dnxy(4, jj) = dns(2, ii)
! modified by Rui at May 27 2015
! refer to 'Robert D, et al. Concepts and application of FEA: isoparametric elements'
END DO nodel ! End loop over nodes

END SUBROUTINE shape2dx

!---------------------------------------------------------------

SUBROUTINE shape2dxx(dnxx, dnyy, x, y, n)
!
! Calculates array containing the second derivatives of the full Lagrangian shape function (nout) at coordinate (x,y) of order n.
! This subroutine is limited to orders of 2 - 4 (for Q4, Q9, and Q16 elements). The outputs are used to calculate the 2D Jacobian
! and the strain displacement matrix.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   10/26/2011     D.N.A.     Original Code
!

IMPLICIT NONE

INTEGER :: ii, jj, kk ! Counter
INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! X-coordinate of point to evaluate shape function at
REAL(kind=REKIND), INTENT(IN) :: y  ! Y-coordinate of point to evaluate shape function at
REAL(kind=REKIND), DIMENSION(n ** 2) :: ddndxdx, ddndxdy, ddndydx, ddndydy ! Shape function array of order n at point (x,y)
REAL(kind=REKIND), INTENT(OUT), DIMENSION(3,2*(n**2)) :: dnxx, dnyy ! 2D Shape functions

order: IF ( n == 2 ) THEN ! Order 2
     ! Shape:
     ! 4 -- 3
     ! |    |
     ! 1 -- 2
     ddndxdx( 1 ) =  0._REKIND   ! N1,xx
     ddndxdx( 2 ) =  0._REKIND   ! N2,xx
     ddndxdx( 3 ) =  0._REKIND   ! N3,xx
     ddndxdx( 4 ) =  0._REKIND   ! N4,xx
     ddndxdy( 1 ) =  0.25_REKIND ! N1,xy
     ddndxdy( 2 ) = -0.25_REKIND ! N2,xy
     ddndxdy( 3 ) =  0.25_REKIND ! N3,xy
     ddndxdy( 4 ) = -0.25_REKIND ! N4,xy
     ddndydy( 1 ) =  0._REKIND   ! N1,yy
     ddndydy( 2 ) =  0._REKIND   ! N2,yy
     ddndydy( 3 ) =  0._REKIND   ! N3,yy
     ddndydy( 4 ) =  0._REKIND   ! N4,yy
     ddndydx( 1 ) =  0.25_REKIND ! N1,yx
     ddndydx( 2 ) = -0.25_REKIND ! N2,yx
     ddndydx( 3 ) =  0.25_REKIND ! N3,yx
     ddndydx( 4 ) = -0.25_REKIND ! N4,yx
ELSE IF ( n == 3 ) THEN order ! Order 3
     ! Shape:
     ! 4 - 7 - 3
     ! |   |   |
     ! 8 - 9 - 6
     ! |   |   |
     ! 1 - 5 - 2
     ddndxdx( 1 ) = 0.5_REKIND * ( y ** 2 - y ) ! N1,xx
     ddndxdx( 2 ) = 0.5_REKIND * ( y ** 2 - y ) ! N2,xx
     ddndxdx( 3 ) = 0.5_REKIND * ( y ** 2 + y ) ! N3,xx
     ddndxdx( 4 ) = 0.5_REKIND * ( y ** 2 + y ) ! N4,xx
     ddndxdx( 5 ) = y - y ** 2 ! N5,xx
     ddndxdx( 6 ) = 1._REKIND - y ** 2 ! N6,xx
     ddndxdx( 7 ) = -y - y ** 2 ! N7,xx
     ddndxdx( 8 ) = 1._REKIND - y ** 2 ! N8,xx
     ddndxdx( 9 ) = -2._REKIND + 2._REKIND * y ** 2 ! N9,xx
     ddndxdy( 1 ) =  0.25_REKIND - 0.5_REKIND * ( x + y ) + x * y ! N1,xy
     ddndxdy( 2 ) = -0.25_REKIND - 0.5_REKIND * ( x - y ) + x * y ! N2,xy
     ddndxdy( 3 ) =  0.25_REKIND + 0.5_REKIND * ( x + y ) + x * y ! N3,xy
     ddndxdy( 4 ) = -0.25_REKIND + 0.5_REKIND * ( x - y ) + x * y ! N4,xy
     ddndxdy( 5 ) =  x - 2._REKIND * x * y ! N5,xy
     ddndxdy( 6 ) = -y - 2._REKIND * x * y ! N6,xy
     ddndxdy( 7 ) = -x - 2._REKIND * x * y ! N7,xy
     ddndxdy( 8 ) =  y - 2._REKIND * x * y ! N8,xy
     ddndxdy( 9 ) = 4._REKIND * x * y ! N9,xy
     ddndydy( 1 ) = 0.5_REKIND * ( x ** 2 - x ) ! N1,yy
     ddndydy( 2 ) = 0.5_REKIND * ( x ** 2 + x ) ! N2,yy
     ddndydy( 3 ) = 0.5_REKIND * ( x ** 2 + x ) ! N3,yy
     ddndydy( 4 ) = 0.5_REKIND * ( x ** 2 - x ) ! N4,yy
     ddndydy( 5 ) = 1._REKIND - x ** 2 ! N5,yy
     ddndydy( 6 ) = -x - x ** 2 ! N6,yy
     ddndydy( 7 ) = 1._REKIND - x ** 2 ! N7,yy
     ddndydy( 8 ) = x - x ** 2  ! N8,yy
     ddndydy( 9 ) = -2._REKIND + 2._REKIND * x ** 2 ! N9,yy
     ddndydx( 1 ) =  0.25_REKIND - 0.5_REKIND * ( x + y ) + x * y ! N1,yx
     ddndydx( 2 ) = -0.25_REKIND - 0.5_REKIND * ( x - y ) + x * y ! N2,yx
     ddndydx( 3 ) =  0.25_REKIND + 0.5_REKIND * ( x + y ) + x * y ! N3,yx
     ddndydx( 4 ) = -0.25_REKIND + 0.5_REKIND * ( x - y ) + x * y ! N4,yx
     ddndydx( 5 ) =  x - 2._REKIND * x * y ! N5,yx
     ddndydx( 6 ) = -y - 2._REKIND * x * y ! N6,yx
     ddndydx( 7 ) = -x - 2._REKIND * x * y ! N7,yx
     ddndydx( 8 ) =  y - 2._REKIND * x * y ! N8,yx
     ddndydx( 9 ) = 4._REKIND * x * y ! N9,yx
ELSE IF (n == 4 ) THEN order ! Order 4
     ! Shape:
     ! 4  - 10 - 9  - 3
     ! |    |    |    |
     ! 11 - 16 - 15 - 8
     ! |    |    |    |
     ! 12 - 13 - 14 - 7
     ! |    |    |    |
     ! 1  - 5  - 6  - 2
     ! N_xx:
     ddndxdx(  1 ) = 0.0703125_REKIND * ( 1._REKIND - 3._REKIND * x ) * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 )
     ddndxdx(  2 ) = 0.0703125_REKIND * ( 1._REKIND + 3._REKIND * x ) * ( -1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 )
     ddndxdx(  3 ) = 0.0703125_REKIND * ( 1._REKIND + 3._REKIND * x ) * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 )
     ddndxdx(  4 ) = 0.0703125_REKIND * ( 1._REKIND - 3._REKIND * x ) * ( -1._REKIND - y + 9._REKIND * y ** 2 + 9._REKIND * y ** 3 )
     ddndxdx(  5 ) = 0.0703125_REKIND * ( 1._REKIND - 9._REKIND * x ) * ( 1._REKIND - y - 9._REKIND * y ** 2 + 9._REKIND * y ** 3 )
     ddndxdx(  6 ) = 0.0703125_REKIND * ( 1._REKIND + 9._REKIND * x ) * ( 1._REKIND - y - 9._REKIND * y ** 2 + 9._REKIND * y ** 3 )
     ddndxdx(  7 ) = 0.6328125_REKIND * ( 1._REKIND + 3._REKIND * x ) * ( 1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 )
     ddndxdx(  8 ) = 0.6328125_REKIND * ( 1._REKIND + 3._REKIND * x ) * ( 1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 )
     ddndxdx(  9 ) = 0.0703125_REKIND * ( 1._REKIND + 9._REKIND * x ) * ( 1._REKIND + y + 9._REKIND * y ** 2 - 9._REKIND * y ** 3 )
     ddndxdx( 10 ) = 0.0703125_REKIND * ( 1._REKIND - 9._REKIND * x ) * ( 1._REKIND + y - 9._REKIND * y ** 2 - 9._REKIND * y ** 3 )
     ddndxdx( 11 ) = 0.6328125_REKIND * ( 1._REKIND - 3._REKIND * x ) * ( 1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 )
     ddndxdx( 12 ) = 0.6328125_REKIND * ( 1._REKIND - 3._REKIND * x ) * ( 1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 )
     ddndxdx( 13 ) = 0.6328125_REKIND * ( -1._REKIND + 9._REKIND * x ) * ( 1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 )
     ddndxdx( 14 ) = 0.6328125_REKIND * ( -1._REKIND - 9._REKIND * x ) * ( 1._REKIND - 3._REKIND * y - y ** 2 + 3._REKIND * y ** 3 )
     ddndxdx( 15 ) = 0.6328125_REKIND * ( -1._REKIND - 9._REKIND * x ) * ( 1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 )
     ddndxdx( 16 ) = 0.6328125_REKIND * ( -1._REKIND + 9._REKIND * x ) * ( 1._REKIND + 3._REKIND * y - y ** 2 - 3._REKIND * y ** 3 )
     ddndxdy(  1 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N1,xy
     ddndxdy(  2 ) = 0.00390625_REKIND * (-1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N2,xy
     ddndxdy(  3 ) = 0.00390625_REKIND * ( 1._REKIND - 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND - 18._REKIND * y - 27._REKIND * y ** 2) ! N3,xy
     ddndxdy(  4 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * (-1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2) ! N4,xy
     ddndxdy(  5 ) = 0.03515625_REKIND * (-3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N5,xy
     ddndxdy(  6 ) = 0.03515625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N6,xy
     ddndxdy(  7 ) = 0.03515625_REKIND * ( 1._REKIND - 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND + 2._REKIND * y - 9._REKIND * y ** 2) ! N7,xy
     ddndxdy(  8 ) = 0.03515625_REKIND * (-1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N8,xy
     ddndxdy(  9 ) = 0.03515625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * (-1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2) ! N9,xy
     ddndxdy( 10 ) = 0.03515625_REKIND * ( 3._REKIND + 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND  - 18._REKIND * y - 27._REKIND * y ** 2) ! N10,xy
     ddndxdy( 11 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N11,xy
     ddndxdy( 12 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * (-3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2) ! N12,xy
     ddndxdy( 13 ) = 0.31640625_REKIND * ( 3._REKIND + 2._REKIND * x - 9._REKIND * x ** 2) &
          & * ( 3._REKIND + 2._REKIND * y - 9._REKIND * y ** 2) ! N13,xy
     ddndxdy( 14 ) = 0.31640625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2) &
          & * (-3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2) ! N14,xy
     ddndxdy( 15 ) = 0.31640625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N15,xy
     ddndxdy( 16 ) = 0.31640625_REKIND * (-3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N16,xy

     ! N_yy
     ddndydy(  1 ) = 0.0703125_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) * ( 1._REKIND - 3._REKIND * y )
     ddndydy(  2 ) = 0.0703125_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) * ( 1._REKIND - 3._REKIND * y )
     ddndydy(  3 ) = 0.0703125_REKIND * ( -1._REKIND - x + 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) * ( 1._REKIND + 3._REKIND * y )
     ddndydy(  4 ) = 0.0703125_REKIND * ( -1._REKIND + x + 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) * ( 1._REKIND + 3._REKIND * y )
     ddndydy(  5 ) = 0.6328125_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) * ( 1._REKIND - 3._REKIND * y )
     ddndydy(  6 ) = 0.6328125_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) * ( 1._REKIND - 3._REKIND * y )
     ddndydy(  7 ) = 0.0703125_REKIND * (  1._REKIND + x - 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) * ( 1._REKIND - 9._REKIND * y )
     ddndydy(  8 ) = 0.0703125_REKIND * (  1._REKIND + x - 9._REKIND * x ** 2 - 9._REKIND * x ** 3 ) * ( 1._REKIND + 9._REKIND * y )
     ddndydy(  9 ) = 0.6328125_REKIND * (  1._REKIND + 3._REKIND * x - x ** 2 - 3._REKIND * x ** 3 ) * ( 1._REKIND + 3._REKIND * y )
     ddndydy( 10 ) = 0.6328125_REKIND * (  1._REKIND - 3._REKIND * x - x ** 2 + 3._REKIND * x ** 3 ) * ( 1._REKIND + 3._REKIND * y )
     ddndydy( 11 ) = 0.0703125_REKIND * (  1._REKIND - x - 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) * ( 1._REKIND + 9._REKIND * y )
     ddndydy( 12 ) = 0.0703125_REKIND * (  1._REKIND - x - 9._REKIND * x ** 2 + 9._REKIND * x ** 3 ) * ( 1._REKIND - 9._REKIND * y )
     ddndydy( 13 ) = 0.6328125_REKIND * ( -1._REKIND + 3._REKIND * x + x ** 2 - 3._REKIND * x ** 3 ) * ( 1._REKIND - 9._REKIND * y )
     ddndydy( 14 ) = 0.6328125_REKIND * ( -1._REKIND - 3._REKIND * x + x ** 2 + 3._REKIND * x ** 3 ) * ( 1._REKIND - 9._REKIND * y )
     ddndydy( 15 ) = 0.6328125_REKIND * ( -1._REKIND - 3._REKIND * x + x ** 2 + 3._REKIND * x ** 3 ) * ( 1._REKIND + 9._REKIND * y )
     ddndydy( 16 ) = 0.6328125_REKIND * ( -1._REKIND + 3._REKIND * x + x ** 2 - 3._REKIND * x ** 3 ) * ( 1._REKIND + 9._REKIND * y )

     ddndydx(  1 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N1,yx
     ddndydx(  2 ) = 0.00390625_REKIND * (-1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N2,yx
     ddndydx(  3 ) = 0.00390625_REKIND * ( 1._REKIND - 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 1._REKIND - 18._REKIND * y - 27._REKIND * y ** 2) ! N3,yx
     ddndydx(  4 ) = 0.00390625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * (-1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2) ! N4,yx
     ddndydx(  5 ) = 0.03515625_REKIND * (-3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N5,yx
     ddndydx(  6 ) = 0.03515625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND + 18._REKIND * y - 27._REKIND * y ** 2) ! N6,yx
     ddndydx(  7 ) = 0.03515625_REKIND * ( 1._REKIND - 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND + 2._REKIND * y - 9._REKIND * y ** 2) ! N7,yx
     ddndydx(  8 ) = 0.03515625_REKIND * (-1._REKIND + 18._REKIND * x + 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N8,yx
     ddndydx(  9 ) = 0.03515625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * (-1._REKIND + 18._REKIND * y + 27._REKIND * y ** 2) ! N9,yx
     ddndydx( 10 ) = 0.03515625_REKIND * ( 3._REKIND + 2._REKIND * x - 9._REKIND * x ** 2 ) &
          & * ( 1._REKIND  - 18._REKIND * y - 27._REKIND * y ** 2) ! N10,yx
     ddndydx( 11 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N11,yx
     ddndydx( 12 ) = 0.03515625_REKIND * ( 1._REKIND + 18._REKIND * x - 27._REKIND * x ** 2 ) &
          & * (-3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2) ! N12,yx
     ddndydx( 13 ) = 0.31640625_REKIND * ( 3._REKIND + 2._REKIND * x - 9._REKIND * x ** 2) &
          & * ( 3._REKIND + 2._REKIND * y - 9._REKIND * y ** 2) ! N13,yx
     ddndydx( 14 ) = 0.31640625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2) &
          & * (-3._REKIND - 2._REKIND * y + 9._REKIND * y ** 2) ! N14,yx
     ddndydx( 15 ) = 0.31640625_REKIND * ( 3._REKIND - 2._REKIND * x - 9._REKIND * x ** 2) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N15,yx
     ddndydx( 16 ) = 0.31640625_REKIND * (-3._REKIND - 2._REKIND * x + 9._REKIND * x ** 2) &
          & * ( 3._REKIND - 2._REKIND * y - 9._REKIND * y ** 2) ! N16,yx
END IF order ! End conditional
dnxx = 0._REKIND
dnyy = 0._REKIND
jj = 0
kk = -1
nodel: DO ii = 1, n ** 2 ! Loop over spatial nodes
     jj = jj + 2
     kk = kk + 2
     dnxx(1,kk) = 0.5_REKIND  * ddndxdx(ii) ! x-direction shape functions
     dnxx(2,jj) = 0.5_REKIND  * ddndxdy(ii) ! y-direction shape functions
     dnxx(3,kk) = ddndxdx(ii) + ddndydy(ii) ! x-direction shape functions
     dnyy(1,jj) = 0.5_REKIND  * ddndydy(ii) ! y-direction shape functions
     dnyy(2,kk) = 0.5_REKIND  * ddndydx(ii) ! y-direction shape functions
     dnyy(3,jj) = ddndxdx(ii) + ddndydy(ii) ! y-direction shape functions
END DO nodel ! End loop over spatial nodes

END SUBROUTINE shape2dxx

!---------------------------------------------------------------

SUBROUTINE shape3x(nout, x, y, z, n)
!
! Calculates array containing the full Lagrangian shape function (nout) at 
! coordinate (x,y,z) of order n. This subroutine is limited to 8-node and 
! 20-node brick elements in 3D
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   09/18/2015   Rui Zhang     Original Code
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! X-coordinate
REAL(kind=REKIND), INTENT(IN) :: y  ! Y-coordinate
REAL(kind=REKIND), INTENT(IN) :: z  ! Z-coordinate
REAL(kind=REKIND), DIMENSION(n**3) :: ns ! Shape function array of order n at point (x,y,z)
REAL(kind=REKIND), INTENT(OUT), DIMENSION(3,3*(n**3)) :: nout ! 3D Shape functions
INTEGER :: ii, jj ! Counter

order: IF ( n == 2 ) THEN ! 8-node brick element
     ns(1) = 0.125_REKIND*(1._REKIND-x)*(1._REKIND-y)*(1._REKIND+z) ! N1
     ns(2) = 0.125_REKIND*(1._REKIND-x)*(1._REKIND-y)*(1._REKIND-z) ! N2
     ns(3) = 0.125_REKIND*(1._REKIND-x)*(1._REKIND+y)*(1._REKIND-z) ! N3
     ns(4) = 0.125_REKIND*(1._REKIND-x)*(1._REKIND+y)*(1._REKIND+z) ! N4
     ns(5) = 0.125_REKIND*(1._REKIND+x)*(1._REKIND-y)*(1._REKIND+z) ! N5
     ns(6) = 0.125_REKIND*(1._REKIND+x)*(1._REKIND-y)*(1._REKIND-z) ! N6
     ns(7) = 0.125_REKIND*(1._REKIND+x)*(1._REKIND+y)*(1._REKIND-z) ! N7
     ns(8) = 0.125_REKIND*(1._REKIND+x)*(1._REKIND+y)*(1._REKIND+z) ! N8
ELSE IF ( n == 3 ) THEN order ! Order 3

END IF order ! End Conditional
nout = 0._REKIND
jj = 0
nodel: DO ii = 1, n ** 3 ! Loop over spatial nodes
     jj = jj + 3
     nout(1, jj-2) = ns(ii) ! x-direction shape functions
     nout(2, jj-1) = ns(ii) ! y-direction shape functions
     nout(3, jj)   = ns(ii) ! z-direction shape functions
END DO nodel ! End loop over spatial nodes

END SUBROUTINE shape3x

!---------------------------------------------------------------

SUBROUTINE shape3dx(dns, dnxyz, x, y, z, n)
!
! Calculates array containing the first derivatives of the full Lagrangian 
! shape function (dns and dnxyz) at coordinate (x,y,z) of order n. This 
! subroutine is limited 8-node and 20-node brick elements. The outputs 
! are used to calculate the 3D Jacobian and the strain displacement matrix.
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   09/18/2015   Rui Zhang     Original Code
!

IMPLICIT NONE

INTEGER :: ii, jj, kk, ll ! Counter
INTEGER, INTENT(IN) :: n  ! Order shape function
REAL(kind=REKIND), INTENT(IN) :: x  ! X-coordinate
REAL(kind=REKIND), INTENT(IN) :: y  ! Y-coordinate
REAL(kind=REKIND), INTENT(IN) :: z  ! Z-coordinate
REAL(kind=REKIND), INTENT(OUT), DIMENSION(3,n**3) :: dns ! 3D Shape functions
REAL(kind=REKIND), INTENT(OUT), DIMENSION(9,3*(n**3)) :: dnxyz ! 3D Shape functions
order: IF ( n == 2 ) THEN ! 8-node brick element
     dns(1, 1) = -0.125_REKIND*(1._REKIND-y)*(1._REKIND+z) ! dN1dx
     dns(1, 2) = -0.125_REKIND*(1._REKIND-y)*(1._REKIND-z) ! dN2dx
     dns(1, 3) = -0.125_REKIND*(1._REKIND+y)*(1._REKIND-z) ! dN3dx
     dns(1, 4) = -0.125_REKIND*(1._REKIND+y)*(1._REKIND+z) ! dN4dx
     dns(1, 5) =  0.125_REKIND*(1._REKIND-y)*(1._REKIND+z) ! dN5dx
     dns(1, 6) =  0.125_REKIND*(1._REKIND-y)*(1._REKIND-z) ! dN6dx
     dns(1, 7) =  0.125_REKIND*(1._REKIND+y)*(1._REKIND-z) ! dN7dx
     dns(1, 8) =  0.125_REKIND*(1._REKIND+y)*(1._REKIND+z) ! dN8dx
     dns(2, 1) = -0.125_REKIND*(1._REKIND-x)*(1._REKIND+z) ! dN1dy
     dns(2, 2) = -0.125_REKIND*(1._REKIND-x)*(1._REKIND-z) ! dN2dy
     dns(2, 3) =  0.125_REKIND*(1._REKIND-x)*(1._REKIND-z) ! dN3dy
     dns(2, 4) =  0.125_REKIND*(1._REKIND-x)*(1._REKIND+z) ! dN4dy
     dns(2, 5) = -0.125_REKIND*(1._REKIND+x)*(1._REKIND+z) ! dN5dy
     dns(2, 6) = -0.125_REKIND*(1._REKIND+x)*(1._REKIND-z) ! dN6dy
     dns(2, 7) =  0.125_REKIND*(1._REKIND+x)*(1._REKIND-z) ! dN7dy
     dns(2, 8) =  0.125_REKIND*(1._REKIND+x)*(1._REKIND+z) ! dN8dy
     dns(3, 1) =  0.125_REKIND*(1._REKIND-x)*(1._REKIND-y) ! dN1dz
     dns(3, 2) = -0.125_REKIND*(1._REKIND-x)*(1._REKIND-y) ! dN2dz
     dns(3, 3) = -0.125_REKIND*(1._REKIND-x)*(1._REKIND+y) ! dN3dz
     dns(3, 4) =  0.125_REKIND*(1._REKIND-x)*(1._REKIND+y) ! dN4dz
     dns(3, 5) =  0.125_REKIND*(1._REKIND+x)*(1._REKIND-y) ! dN5dz
     dns(3, 6) = -0.125_REKIND*(1._REKIND+x)*(1._REKIND-y) ! dN6dz
     dns(3, 7) = -0.125_REKIND*(1._REKIND+x)*(1._REKIND+y) ! dN7dz
     dns(3, 8) =  0.125_REKIND*(1._REKIND+x)*(1._REKIND+y) ! dN8dz
ELSE IF ( n == 3 ) THEN order ! 20-node brick element
     
END IF order ! End conditional

! Initialize x/y direction shape functions
dnxyz = 0._REKIND
jj = 0
kk = -1
ll = -2
nodel: DO ii = 1, n ** 3 ! Loop over nodes
     jj = jj + 3
     kk = kk + 3
     ll = ll + 3
     dnxyz(1, ll) = dns(1, ii)
     dnxyz(2, ll) = dns(2, ii)
     dnxyz(3, ll) = dns(3, ii)
     dnxyz(4, kk) = dns(1, ii)
     dnxyz(5, kk) = dns(2, ii)
     dnxyz(6, kk) = dns(3, ii)
     dnxyz(7, jj) = dns(1, ii)
     dnxyz(8, jj) = dns(2, ii)
     dnxyz(9, jj) = dns(3, ii)
END DO nodel ! End loop over nodes

END SUBROUTINE shape3dx

!---------------------------------------------------------------

SUBROUTINE feassembly(globmat, elemat, nodes, n, ndof, ndofn)
!
! Assembles elemental spatial matrix into global spatial matrix
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   09/19/2015   Rui Zhang    Add 3D assembly
!   05/02/2009     D.N.A.     Original Code
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: n, ndofn ! Nodes per element, DOF per node
INTEGER, INTENT(IN) :: ndof ! total global DOF
INTEGER, DIMENSION(n), INTENT(IN) :: nodes ! Nodal connectivity of element
REAL(kind=REKIND), DIMENSION(n * ndofn, n * ndofn), INTENT(IN) :: elemat ! Elemental matrix
REAL(kind=REKIND), DIMENSION(ndof, ndof), INTENT(INOUT) :: globmat ! Global matrix

INTEGER :: ii, jj, kk, mm ! Elemental DOF, indices
dimcheck: IF(ndofn == 1) THEN ! 1D problem
     DO ii = 1, n ! Loop over elemental DOFs, local matrix column index
          kk = nodes(ii) ! Global matrix column index
          DO jj = 1, n ! Loop over elemental DOFs, local matrix row index
               mm = nodes(jj) ! Global matrix row index
               globmat(mm, kk) = globmat(mm, kk) + elemat(jj, ii) ! Add local matrix to global matrix
          END DO ! End loop over elemental DOFs, local matrix column index
     END DO ! End loop over elemental DOFs, local matrix row index
ELSE IF(ndofn == 2) THEN dimcheck ! 2D problem
     DO ii = 1, n ! Loop over elemental DOFs, local matrix column index
          kk = nodes(ii) + nodes(ii) ! Global matrix column index
          DO jj = 1, n ! Loop over elemental DOFs, local matrix row index
               mm = nodes(jj) + nodes(jj) ! Global matrix row index
               ! Add local matrix to global matrix
               globmat(mm-1:mm, kk-1:kk) = globmat(mm-1:mm, kk-1:kk) + elemat(2*jj-1:2*jj, 2*ii-1:2*ii)
          END DO ! End loop over elemental DOFs, local matrix column index
     END DO ! End loop over elemental DOFs, local matrix row index
ELSE IF(ndofn == 3) THEN dimcheck ! 3D problem
     DO ii = 1, n ! Loop over elemental DOFs, local matrix column index
          kk = nodes(ii)*3 ! Global matrix column index
          DO jj = 1, n ! Loop over elemental DOFs, local matrix row index
               mm = nodes(jj)*3! Global matrix row index
               ! Add local matrix to global matrix
               globmat(mm-2:mm, kk-2:kk) = globmat(mm-2:mm, kk-2:kk) + elemat(3*jj-2:3*jj, 3*ii-2:3*ii)
          END DO ! End loop over elemental DOFs, local matrix column index
     END DO ! End loop over elemental DOFs, local matrix row index
ELSE dimcheck    
END IF dimcheck ! conditional dimension check

END SUBROUTINE feassembly

!---------------------------------------------------------------

SUBROUTINE spatialfext(ndof, ndofn, neh, nel, nn, nnl, nnsd, nnsdmo, nquads, wqs, xiqs, phia, xcoord, f_ext)
!
! Calculates and assembles spatial component of the external force (surface traction)
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   03/06/2012     D.N.A.     Moved code from main.f03

IMPLICIT NONE

INTEGER, INTENT(IN) :: ndof, ndofn, neh, nel, nn, nnl, nnsd, nnsdmo, nquads
REAL(kind=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs, xiqs
REAL(kind=REKIND), DIMENSION(2,4), INTENT(IN) :: phia
REAL(kind=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: xcoord
REAL(kind=REKIND), DIMENSION(ndof), INTENT(OUT) :: f_ext
INTEGER :: ii, jj, kk, ll
REAL(kind=REKIND) :: f_ebottomx, f_ebottomy, f_eleftx, f_elefty, f_erightx, f_erighty, f_etopx, f_etopy, jacbottom, jacleft, &
     & jacright, jactop, temp1
REAL(kind=REKIND), DIMENSION(nnsd) :: phiabottomx, phiabottomy, phialeftx, phialefty, phiarightx, phiarighty, phiatopx, phiatopy, &
     & tempns1, xcoordebottom, xcoordeleft, xcoorderight, xcoordetop
REAL(kind=REKIND), DIMENSION(nnsd, nquads) :: dns, ns

DO ii = 1, nquads ! Loop over quadrature points
     CALL shape1x(tempns1, xiqs(ii), nnsd) ! Calculate 1D shape function at Gauss point
     ns(:,ii) = tempns1 ! 1D shape function at all Gauss points
     CALL shape1dx(tempns1, xiqs(ii), nnsd) ! Calculate 1D shape function derivative at Gauss point
     dns(:,ii) = tempns1 ! 1D shape function derivative at all shape function
END DO ! End loop over quadrature points
! Line element traction amplitudes, in N/m
phiabottomx = phia(1,1)
phiabottomy = phia(2,1)
phiatopx    = phia(1,3)
phiatopy    = phia(2,3)
phiarightx  = phia(1,2)
phiarighty  = phia(2,2)
phialeftx   = phia(1,4)
phialefty   = phia(2,4)
f_ext = 0._REKIND !Initialize force vector to 0 N
ll = nnl+nnl
temp1 = REAL(nnsdmo + nnsdmo, REKIND)
DO ii = nnsdmo, nel*nnsdmo, nnsdmo ! Loop over elements axially
     ! Elemental spatial force amplitude, in N
     f_ebottomx = 0._REKIND
     f_ebottomy = 0._REKIND
     f_etopx    = 0._REKIND
     f_etopy    = 0._REKIND
     ! Elemental nodal coordinate vectors, in m
     xcoordebottom = xcoord( ii - nnsd + 2 : ii + 1, 1 )
     xcoordetop    = xcoord( nn - nnl + ii - nnsd + 2 : nn-nnl + ii + 1, 1 )
     DO jj = 1, nquads ! Loop over quadrature points
          jacbottom = DOT_PRODUCT(xcoordebottom,dns(:,jj)) ! Jacobian for bottom line element
          jactop = DOT_PRODUCT(xcoordetop,dns(:,jj)) ! Jacobian for top line element
          ! Elemental external force amplitude vector, in N
          f_ebottomx = f_ebottomx + DOT_PRODUCT(ns(:,jj),phiabottomx) * jacbottom * wqs(jj) ! For bottom line element
          f_ebottomy = f_ebottomy + DOT_PRODUCT(ns(:,jj),phiabottomy) * jacbottom * wqs(jj) ! For bottom line element
          f_etopx = f_etopx + DOT_PRODUCT(ns(:,jj),phiatopx) * jactop * wqs(jj) ! For top line element
          f_etopy = f_etopy + DOT_PRODUCT(ns(:,jj),phiatopy) * jactop * wqs(jj) ! For top line element
     END DO ! End loop over quadrature points
     DO jj = 1 + ii + ii - nnsdmo - nnsdmo, ii + ii, 2 ! Loop over nodes per line element for force vector assembly
          kk = ndof - ll
          ! External force amplitude vector, in N
          f_ext(jj:jj+3) = f_ext(jj:jj+3) + [f_ebottomx, f_ebottomy, f_ebottomx, f_ebottomy] / temp1
          f_ext(kk+jj:kk+jj+3) = f_ext(kk+jj:kk+jj+3) + [f_etopx, f_etopy, f_etopx, f_etopy] / temp1
     END DO ! End loop over nodes per line element for force vector assembly

END DO ! End loop over elements axially
DO ii = 0, (neh-1)*nnsdmo*nnl, nnsdmo*nnl ! Loop over elements radially
     ! Elemental spatial force amplitude, in N
     f_erightx = 0._REKIND
     f_erighty = 0._REKIND
     f_eleftx  = 0._REKIND
     f_elefty  = 0._REKIND
     ! Elemental nodal coordinate vectors, in m
     xcoorderight = xcoord( ii+nnl : ii + nnsd*nnl : nnl, 2 )
     xcoordeleft  = xcoord( ii + 1 : ii + nnsdmo*nnl + 1 : nnl, 2 )
     DO jj = 1, nquads ! Loop over quadrature points
          jacright = DOT_PRODUCT(xcoorderight, dns(:,jj)) ! Jacobian for bottom line element
          jacleft  = DOT_PRODUCT(xcoordeleft,  dns(:,jj)) ! Jacobian for top line element
          ! Elemental external force amplitude vector updates, in N
          f_erightx = f_erightx + DOT_PRODUCT(ns(:,jj),phiarightx) * jacright * wqs(jj) ! For bottom line element
          f_erighty = f_erighty + DOT_PRODUCT(ns(:,jj),phiarighty) * jacright * wqs(jj) ! For bottom line element
          f_eleftx  = f_eleftx  + DOT_PRODUCT(ns(:,jj),phialeftx)  * jacleft  * wqs(jj) ! For top line element
          f_elefty  = f_elefty  + DOT_PRODUCT(ns(:,jj),phialefty)  * jacleft  * wqs(jj) ! For top line element
     END DO ! End loop over quadrature points
     DO jj = 2*ii + ll-1, 2*ii + nnsdmo*ll-1, ll ! Loop over nodes per line element for force vector assembly
          ! External force amplitude vector, in N
          f_ext(jj      : jj+3   ) = f_ext(jj   : jj+3   ) + [f_erightx, f_erighty, f_eleftx, f_elefty] / temp1
          f_ext(jj+ll   : jj+ll+1) = f_ext(jj+ll   : jj+ll+1) + [f_erightx, f_erighty] / temp1
          f_ext(jj-ll+2 : jj-ll+3) = f_ext(jj-ll+2 : jj-ll+3) + [f_eleftx,  f_elefty]  / temp1
     END DO ! End loop over nodes per line element for force vector assembly
END DO ! End loop over elements radially

END SUBROUTINE spatialfext

!---------------------------------------------------------------

SUBROUTINE sinIC1d(ampic, le, natfreq, nt, nx, tcoord, xcoord, fmt1)
!
! Calculates 1D analytical displacement history for a straight fixed fixed bar with no external force, no initial velocity, and
! initial displacement proportional to the first harmonic of the bar
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   03/06/2012     D.N.A.     Moved code from main.f03

IMPLICIT NONE
INTEGER, INTENT(IN) :: nt, nx
REAL(kind=REKIND), INTENT(IN) :: ampic, le, natfreq
REAL(kind=REKIND), DIMENSION(nt), INTENT(IN) :: tcoord
REAL(kind=REKIND), DIMENSION(nx), INTENT(IN) :: xcoord
CHARACTER(len=100), INTENT(IN) :: fmt1
INTEGER :: ii, errstat
REAL(kind=REKIND), DIMENSION(nt) :: omegat
REAL(kind=REKIND), DIMENSION(nx) :: dispa, xnorm
CHARACTER(len=100) :: errmesg

! Open file analyticaldisphist.dat to save analytical history to.
OPEN(UNIT=41, FILE='analyticaldisphist.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=errstat, IOMSG=errmesg)
IF (errstat /= 0) WRITE(*,*) 'Error opening errordata.dat. Code: ', errstat, '. Message: ', errmesg  ! Error message
xnorm = xcoord / le ! Normalized Lagrangian coordinates multiplied by maximum displacement, in m
DO ii = 1, nx ! Loop over spatial points
     xnorm(ii) = ampic * SIN(PI * xnorm(ii)) ! Calculate initial displacement, in m
END DO ! End loop over spatial points
omegat = natfreq * tcoord ! Phase angle of displacement wave at each time point, in radians
DO ii = 1, nt ! Loop over time points
     dispa = xnorm * COS(omegat(ii)) ! Displacement at current time point, in m
     WRITE(UNIT=41, FMT=TRIM(fmt1), IOSTAT=errstat, IOMSG=errmesg) dispa ! Save analytical displacement history to file
     IF(errstat/=0) WRITE(*,*) 'Error writing analytical disp history to file. Code: ', errstat, '. Message: ', errmesg 
END DO ! End loop over time points
ENDFILE(UNIT=41, IOSTAT=errstat, IOMSG=errmesg) ! Write end of file record
IF (errstat /= 0) WRITE(*,*) 'Error ending analyticaldisphist.dat. Code: ', errstat, '. Message: ', errmesg  ! Error message
REWIND(UNIT=41, IOSTAT=errstat, IOMSG=errmesg) ! Write end of file record
IF (errstat /= 0) WRITE(*,*) 'Error rewinding analyticaldisphist.dat. Code: ', errstat, '. Message: ', errmesg  ! Error message
CLOSE (UNIT=41, STATUS='KEEP', IOSTAT=errstat, IOMSG=errmesg) ! Close file
IF (errstat /= 0) WRITE(*,*) 'Error closing analyticaldisphist.dat. Code: ', errstat, '. Message: ', errmesg  ! Error message

END SUBROUTINE sinIC1d

!---------------------------------------------------------------

SUBROUTINE stressdispplot(dt, fmt1, ndof, ndofn, ne, neh, nel, nframes, nn, nnl, nnsd, nnsdmo, nnse, nnte, nstep, sd, stiff, &
     & xconn2, xcoord2, xmag, ymag, zmag)
!
! Create PLT files of the displacement and stress history for input into TECPLOT
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   04/06/2010     D.N.A.     Start work on integrating directly into STFEA framework
!   10/14/2010     D.N.A.     Add algorithm to produce algorithm that will take advantage of independent stress algorithm
!   01/28/2011     D.N.A.     Added Kirchhoff material model, changed neohookeanmodel module to nlconstitutivemodel
!   03/13/2012     D.N.A.     Adapted subroutine from archive to fit into current implementation. Removed nonlinear materials.

! NOTE: Once the Fortran 2008 subroutine 'EXECUTE_COMMAND_LINE' becomes common it can be used to ZIP and RM the PLT files after
!       they are created. As of March 13 2012, only the latest versions of GNU's Fortran compiler has implemented this subroutine.

INTEGER, INTENT(IN) :: ndof, ndofn, ne, neh, nel, nframes, nn, nnl, nnsd, nnsdmo, nnse, nnte, nstep
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn2
REAL(kind=REKIND), INTENT(IN) :: dt
REAL(kind=REKIND), INTENT(IN), OPTIONAL :: xmag, ymag, zmag
REAL(kind=REKIND), DIMENSION(3,3), INTENT(IN) :: stiff
REAL(kind=REKIND), DIMENSION(3,4), INTENT(IN) :: sd
REAL(kind=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: xcoord2
CHARACTER(len=100), INTENT(IN) :: fmt1
INTEGER :: errstat, framenum, ii, jj, kk, ll, mm, nnhs, nnls, pp, slabnum
INTEGER, DIMENSION(ne, 4) :: xconn
INTEGER, DIMENSION(nnse) :: nodes
REAL(kind=REKIND) :: isot, tframe, tframeincrement, tstep, xmag2, ymag2, zmag2
REAL(kind=REKIND), DIMENSION(3) :: sigmat
REAL(kind=REKIND), DIMENSION(4) :: sigmae
REAL(kind=REKIND), DIMENSION(nnte) :: ntframe
REAL(kind=REKIND), DIMENSION(nnse+nnse) :: dframe
REAL(kind=REKIND), DIMENSION(2, 2) :: gam, jac
REAL(kind=REKIND), DIMENSION(2, 4) :: dplot
REAL(kind=REKIND), DIMENSION(2, nnse) :: tempdns
REAL(kind=REKIND), DIMENSION(2, nnse+nnse) :: tempns2
REAL(kind=REKIND), DIMENSION(3, 4) :: stiffsd
REAL(kind=REKIND), DIMENSION(4, 4) :: gamexpand
REAL(kind=REKIND), DIMENSION(4, nnse+nnse) :: bmat, tempdnxy
REAL(kind=REKIND), DIMENSION(8, nnse) :: dns
REAL(kind=REKIND), DIMENSION(8, nnse+nnse) :: ns
REAL(kind=REKIND), DIMENSION(16, nnse+nnse) :: dnxy
REAL(kind=REKIND), DIMENSION(nnse, 2) :: ncoorde
REAL(kind=REKIND), DIMENSION(nnse+nnse, nnte) :: dslab
REAL(kind=REKIND), DIMENSION(ndof, nnte) :: disphistory
REAL(kind=REKIND), ALLOCATABLE, DIMENSION(:) :: estress, weight, xdisp, ydisp, zdisp
REAL(kind=REKIND), ALLOCATABLE, DIMENSION(:,:) :: displaced, xcoord
CHARACTER(len=100) :: decimalc, dtc, errmesg, fmt2, fmt3, iwidthc, nec, nnc, widthc

WRITE(widthc,*) PRECISION(1._REKIND) + 7
WRITE(iwidthc,*) RANGE(1)
WRITE(decimalc,*) PRECISION(1._REKIND)
iwidthc  = ADJUSTL(iwidthc)
widthc   = ADJUSTL(widthc)
decimalc = ADJUSTL(decimalc)
WRITE(fmt2,*) '(1X,4(1X,1ES'//TRIM(widthc)//'.'//TRIM(decimalc)//'))'
WRITE(fmt3,*) '(1X,4(1X,1I'//TRIM(iwidthc)//'))'
fmt2 = ADJUSTL(fmt2)
fmt3 = ADJUSTL(fmt3)

! Set defaults for magnitude factors to multply displacement by
IF(PRESENT(xmag)) THEN
     xmag2 = xmag
ELSE
     xmag2 = 1._REKIND
END IF
IF(PRESENT(ymag)) THEN
     ymag2 = ymag
ELSE
     ymag2 = 1._REKIND
END IF
IF(PRESENT(zmag)) THEN
     zmag2 = zmag
ELSE
     zmag2 = 1._REKIND
END IF

nnls = nel + 1 ! Number of plotting nodes lengthwise
nnhs = neh + 1 ! Number of plotting nodes heightwise
tstep = REAL(nstep + nstep,REKIND)/REAL(nframes-1,REKIND) ! Isoparametric coordinate step of time slabs between frames
tframeincrement = REAL(nstep, REKIND)/REAL(nframes-1,REKIND)*dt ! Increment of time between frames, in s
stiffsd = MATMUL(stiff, sd) ! Multiply stiffness matrix and isoparametric strain displacement relationship outside of loops

! Conversion of # plotting nodes and # elments to CHAR
WRITE(nnc, FMT=*) nnhs*nnls
WRITE(nec, FMT=*) ne
nnc = ADJUSTL(nnc)
nec = ADJUSTL(nec)

ALLOCATE(displaced(nnls*nnhs,4), estress(nnhs*nnls), weight(nnhs*nnls), xcoord(nnls*nnhs,3), xdisp(nnhs*nnls), ydisp(nnhs*nnls), &
     & zdisp(nnhs*nnls))!, STAT=errstat, ERRMSG=errmesg)
!IF(errstat/=0) WRITE(*,*) 'Error allocating variables Location A. Code: ', errstat, '. Message: ', errmesg

xconn = 0 ! Intialize xconn to 0
ll = 0 - nel ! Initialize index
DO ii = 0, neh*nnls - nnls, nnls ! Loop over elements heightwise
     ll = ll + nel ! Iterate index
     DO jj = 1, nel ! Loop over elements lengthwise
          kk = ll + jj ! Plotting element #
          xconn(kk,1) = ii + jj ! Lower left node # of plotting element
          xconn(kk,2) = ii + jj + 1 ! Lower right node # of plotting element
          xconn(kk,3) = ii + jj + 1 + nnls ! Upper right node # of plotting element
          xconn(kk,4) = ii + jj + nnls ! Upper left node # of plotting element
     END DO ! End loop over elements lengthwise
END DO ! End loop over elements heightwise

xcoord = 0._REKIND !Initialize x-coordinate matrix to 0. (sets z-axis to 0. for 2D problems)
jj = 0 ! Initialize index
kk = nnl - nnsdmo*nnl ! Initialize index
DO ii = 1, nnhs ! Loop over plotting nodes heightwise
     jj = jj + nnls ! Iterate index
     kk = kk + nnsdmo*nnl ! Iterate index
     xcoord( jj - nnls + 1 : jj, 1:2 ) = xcoord2( kk - nnl+1 : kk : nnsdmo, 1:2 ) ! Plotting nodal coordinates, in m
END DO ! End loop over plotting nodes heigthwise
! Spatial shape functions
kk = 0 ! Initialize counter
DO ii = -1, 1, 2 ! Loop over plotting nodes per element side
     DO jj = -1, 1, 2 ! Loop over plotting nodes per element side
          kk = kk + 2 ! Iterate counter
          CALL shape2dx(tempdns, tempdnxy, REAL(ii,REKIND), REAL(jj,REKIND), nnsd) ! Substitute node in dN/dx
          dns(kk-1 : kk, :) = tempdns ! Assign tempdns to ordered matrix of shape functions
          dnxy(kk+kk-3 : kk+kk, :) = tempdnxy ! Assign tempdnxy to ordered matrix of shape functions
          CALL shape2x(tempns2, REAL(ii,REKIND), REAL(jj,REKIND), nnsd) ! Substitute node in N
          ns(kk-1 : kk,:) = tempns2 ! Assign tempns2 to ordered matrix of shape functions
     END DO ! End loop over plotting nodes per element side
END DO ! End loop over plotting nodes per element side

OPEN(UNIT=41, FILE='disphist.dat', STATUS='OLD', ACTION='READ', POSITION='REWIND', IOSTAT=errstat, IOMSG=errmesg) ! Open file
IF (errstat /= 0) WRITE(*,*) 'Error opening disphist.dat. Code: ', errstat, '. Message: ', errmesg ! Error message

framenum = 0 ! Intialize counter for number of frames
slabnum = 1 ! Initialize counter for number of slabs read from disphist
WRITE(*,FMT='(1X,"Plot Loop ",I5,"/",I5)') framenum+1, nframes
DO ii = 1, nnte ! Loop over time nodes per slab
     READ(UNIT=41, FMT=TRIM(fmt1), IOSTAT=errstat, IOMSG=errmesg) disphistory(:, ii) ! Read file
     IF(errstat/=0) WRITE(*,*) 'Error reading disphist at frame ', framenum,' slab ', slabnum, ' node ', ii, '. Code: ', errstat, &
          & '. Message: ', errmesg ! Error message
END DO ! End loop over time nodes per slab

isot = -1._REKIND ! Isoparametric coordinate of initial time/frame
CALL shape1x(ntframe, isot, nnte) ! Calculate temporal shape function at initial time/frame 0.0
estress = 0.0_REKIND ! Initialize estress to 0.0
weight  = 0.0_REKIND ! Initialize weight to 0.0
xdisp   = 0.0_REKIND ! Initialize x displacement to 0.0
ydisp   = 0.0_REKIND ! Initialize x displacement to 0.0
zdisp   = 0.0_REKIND ! Initialize z displacement to 0.0
DO ii = 1, ne ! Loop over elements
     nodes = xconn2(ii, :) ! Nodes for current element in loop  
     ncoorde = xcoord2(nodes, :) !Coordinates for nodes in current element in loop
     kk = 0 ! Initialize index
     DO jj = 1, nnse ! Loop over nodes per spatial element
          kk = kk + 2 ! Iterate index
          ll = 2*nodes(jj) ! Initialize index
          dslab(kk - 1 : kk, :) = disphistory(ll - 1 : ll, :) ! Nodal displacement of slab, in m
     END DO ! End loop over nodes per spatial element
     dframe = MATMUL(dslab, ntframe) ! Displacement of current element in current frame
     ll = 0 ! Initialize index
     mm = 0 ! Initialize index
     DO jj = 1, 2 ! Loop over nodes of plotting elements lengthwise
          DO kk = 1, 2 ! Loop over nodes of plotting elements heightwise
               ll = ll + 2 ! Iterate index
               mm = mm + 1 ! Iterate index
               jac = MATMUL(dns(ll-1:ll, :), ncoorde) ! Deformation matrix
               gam = inverse(jac, 2) ! Gamma = inverse of Deformation matrix
               gamexpand = 0._REKIND ! Intialize gamexpand
               gamexpand(1:2, 1:2) = gam ! Assign gam to gamexpand
               gamexpand(3:4, 3:4) = gam ! Assign gam to gamexpand
               bmat = MATMUL(gamexpand, dnxy(ll + ll - 3: ll + ll, :)) ! strain displacement matrix
               sigmat = MATMUL(stiffsd, MATMUL(bmat, dframe)) ! Stress in Pa at time of frame and at node (ll/2)
               sigmae(mm) = SQRT(sigmat(1)**2 - sigmat(1)*sigmat(2) + sigmat(2)**2 + 3._REKIND*sigmat(3)**2)
               dplot(:,mm) = MATMUL(ns(ll-1:ll, :), dframe)
          END DO ! End loop over nodes of plotting elements heightwise
     END DO ! End loop over nodes of plotting elements lengthwise
     DO jj = 1, 4 ! Loop over plotting nodes in element
          ! Effective stress at plotting nodes (unweighted)
          estress(xconn(ii,jj)) = estress(xconn(ii,jj)) + sigmae(jj)
          ! Effective displacement at plotting nodes (unweighted)
          xdisp(xconn(ii,jj)) = xdisp(xconn(ii,jj)) + dplot(1, jj)
          ydisp(xconn(ii,jj)) = ydisp(xconn(ii,jj)) + dplot(2, jj)
          ! Increment weight for eact time node that appears in the connectivity
          weight(xconn(ii,jj)) = weight(xconn(ii,jj)) + 1._REKIND
     END DO ! End loop ove plotting nodes in element
END DO ! End loop over elements
DO ii = 1, nnhs*nnls ! Loop over plotting nodes
     xdisp(ii) = xdisp(ii) / weight(ii) ! X displacement, in m
     ydisp(ii) = ydisp(ii) / weight(ii) ! X displacement, in m
     estress(ii) = estress(ii) / weight(ii) ! Effective stress, in Pa
END DO ! End loop over plotting nodes

tframe = 0._REKIND ! Time of current frame (initial frame)
WRITE(dtc,*) tframe ! Write time of frame
dtc = ADJUSTL(dtc) ! Adjust dtc to be flush left

OPEN (UNIT=42, FILE='stressdisp.plt', STATUS='REPLACE', ACTION='WRITE', IOSTAT=errstat, IOMSG=errmesg)! Open file
IF (errstat /= 0) WRITE(*,*) 'Error opening stressdisp.plt. Code: ', errstat, '. Message: ', errmesg ! Error message
WRITE(UNIT=42,FMT='(1X,A,17(/,A))',IOSTAT=errstat, IOMSG=errmesg) 'TITLE="2D Fixed Fixed Beam"', 'VARIABLES = "X"', '"Y"', '"Z"', &
     & '"EffectiveStress"', 'TEXT', 'CS=FRAME', 'X=50.0,Y=75.0', 'C=BLACK', 'S=GLOBAL', 'HU=POINT', 'LS=1 AN=MIDCENTER', &
     & 'BXM=20 LT=0.1 BXO=BLACK BXF=WHITE', 'F=HELV-BOLD', 'H=14 A=0', 'MFC=""', 'CLIPPING=CLIPTOVIEWPORT', &
     & 'T="Time &(SOLUTIONTIME%.4f) s"'
IF (errstat /= 0) WRITE(*,*) 'Error writing initial preamble to output file. Code: ', errstat, '. Message: ', errmesg
WRITE(UNIT=42,FMT='(1X,A,3(/,A))',IOSTAT=errstat, IOMSG=errmesg) 'ZONE STRANDID=1, SOLUTIONTIME='//TRIM(dtc), &
     & 'Nodes='//TRIM(nnc)//', Elements='//TRIM(nec)//', ZONETYPE=FEQuadrilateral', 'DATAPACKING=POINT', &
     & 'DT=(SINGLE SINGLE SINGLE SINGLE )'
IF (errstat /= 0) WRITE(*,*) 'Error writing preamble to output file at ', 0, '. Code: ', errstat, '. Message: ', errmesg
 


! displaced is a matrix order (number of plotting nodes, 4).
! The 1st column is the x coordinate
! The 2nd column is the y coordinate
! The 3rd column is the z coordinate
! The 4th column is the effective stress
displaced(:,1)=xcoord(:,1) + xmag2 * xdisp
displaced(:,2)=xcoord(:,2) + ymag2 * ydisp
displaced(:,3)=xcoord(:,3) + zmag2 * zdisp
displaced(:,4)=estress

! Write out displaced to file
DO ii = 1, nnhs*nnls ! Loop over plotting nodes
     WRITE(UNIT=42,FMT=TRIM(fmt2),IOSTAT=errstat, IOMSG=errmesg) displaced(ii,:)
     IF (errstat /= 0) WRITE(*,*) 'Error writing displaced to output file. Code: ', errstat, '. Message: ', errmesg
END DO ! End loop over plotting nodes

! Write out nodal connectivity matrix to file
DO ii = 1, ne ! Loop over plotting elements
     WRITE(UNIT=42,FMT=TRIM(fmt3),IOSTAT=errstat, IOMSG=errmesg) xconn(ii,:)
     IF (errstat /= 0) WRITE(*,*) 'Error writing xconn to output file. Code: ', errstat, '. Message: ', errmesg ! Error message
END DO ! End loop over plotting elements

DO framenum = 1, nframes - 1 ! Loop over frames
     WRITE(*,FMT='(1X,"Plot Loop ",I5,"/",I5)') framenum+1, nframes
     isot = isot + tstep ! Iterate isoparametric time frame coordinate
     isotcheck: DO ! Check if isot is in current slab
          IF(isot <= 1._REKIND) EXIT isotcheck ! Current frame in current slab
          IF((framenum == nframes - 1).AND.(slabnum==nstep)) THEN ! Check if last frame
               isot = 1._REKIND ! Set isoparametric time coordinate to end of slab
               EXIT isotcheck ! Continue calculated stress
          END IF ! End check if last frame
          isot = isot - 2._REKIND ! Reset isot by one slab (isoparametric length of 2.0)
          slabnum = slabnum + 1 ! Iterate number of slabs read from disphist
          DO ii = 1, nnte ! Loop over time nodes per slab
               READ(UNIT=41, FMT=TRIM(fmt1), IOSTAT=errstat, IOMSG=errmesg) disphistory(:, ii) ! Read file
               IF(errstat/=0) WRITE(*,*) 'Error reading disphist at frame ', framenum,' slab ', slabnum, ' node ', ii, '. Code: ', &
                    & errstat, '. Message: ', errmesg ! Error message
          END DO ! End loop over time nodes per slab
     END DO isotcheck ! End check if isot is in current slab

     CALL shape1x(ntframe, isot, nnte) ! Calculate temporal shape function at initial time/frame 0.0
     estress = 0.0_REKIND ! Initialize estress to 0.0
     weight  = 0.0_REKIND ! Initialize weight to 0.0
     xdisp   = 0.0_REKIND ! Initialize x displacement to 0.0
     ydisp   = 0.0_REKIND ! Initialize x displacement to 0.0
     zdisp   = 0.0_REKIND ! Initialize z displacement to 0.0
     DO ii = 1, ne ! Loop over elements
          nodes = xconn2(ii, :) ! Nodes for current element in loop  
          ncoorde = xcoord2(nodes, :) !Coordinates for nodes in current element in loop
          kk = 0 ! Initialize index
          DO jj = 1, nnse ! Loop over nodes per spatial element
               kk = kk + 2 ! Iterate index
               ll = 2*nodes(jj) ! Initialize index
               dslab(kk - 1 : kk, :) = disphistory(ll - 1 : ll, :) ! Nodal displacement of slab, in m
          END DO ! End loop over nodes per spatial element
          dframe = MATMUL(dslab, ntframe) ! Displacement of current element in current frame
          ll = 0 ! Initialize index
          mm = 0 ! Initialize index
          DO jj = 1, 2 ! Loop over nodes of plotting elements lengthwise
               DO kk = 1, 2 ! Loop over nodes of plotting elements heightwise
                    ll = ll + 2 ! Iterate index
                    mm = mm + 1 ! Iterate index
                    jac = MATMUL(dns(ll-1:ll, :), ncoorde) ! Deformation matrix
                    gam = inverse(jac, 2) ! Gamma = inverse of Deformation matrix
                    gamexpand = 0._REKIND ! Intialize gamexpand
                    gamexpand(1:2, 1:2) = gam ! Assign gam to gamexpand
                    gamexpand(3:4, 3:4) = gam ! Assign gam to gamexpand
                    bmat = MATMUL(gamexpand, dnxy(ll + ll - 3: ll + ll, :)) ! strain displacement matrix
                    sigmat = MATMUL(stiffsd, MATMUL(bmat, dframe)) ! Stress in Pa at time of frame and at node (ll/2)
                    sigmae(mm) = SQRT(sigmat(1)**2 - sigmat(1)*sigmat(2) + sigmat(2)**2 + 3._REKIND*sigmat(3)**2)
                    dplot(:,mm) = MATMUL(ns(ll-1:ll, :), dframe)
               END DO ! End loop over nodes of plotting elements heightwise
          END DO ! End loop over nodes of plotting elements lengthwise
          DO jj = 1, 4 ! Loop over plotting nodes in element
               ! Effective stress at plotting nodes (unweighted)
               estress(xconn(ii,jj)) = estress(xconn(ii,jj)) + sigmae(jj)
               ! Effective displacement at plotting nodes (unweighted)
               xdisp(xconn(ii,jj)) = xdisp(xconn(ii,jj)) + dplot(1, jj)
               ydisp(xconn(ii,jj)) = ydisp(xconn(ii,jj)) + dplot(2, jj)
               ! Increment weight for eact time node that appears in the connectivity
               weight(xconn(ii,jj)) = weight(xconn(ii,jj)) + 1._REKIND
          END DO ! End loop ove plotting nodes in element
     END DO ! End loop over elements
     DO ii = 1, nnhs*nnls ! Loop over plotting nodes
          xdisp(ii) = xdisp(ii) / weight(ii) ! X displacement, in m
          ydisp(ii) = ydisp(ii) / weight(ii) ! X displacement, in m
          estress(ii) = estress(ii) / weight(ii) ! Effective stress, in Pa
     END DO ! End loop over plotting nodes

     ! displaced is a matrix order (number of plotting nodes, 4).
     ! The 1st column is the x coordinate at frame jj
     ! The 2nd column is the y coordinate at frame jj
     ! The 3rd column is the z coordinate at frame jj
     ! The 4th column is the effective stress at frame jj
     displaced(:,1)=xcoord(:,1) + xmag2 * xdisp
     displaced(:,2)=xcoord(:,2) + ymag2 * ydisp
     displaced(:,3)=xcoord(:,3) + zmag2 * zdisp
     displaced(:,4)=estress

     IF((isot == 1._REKIND).AND.(framenum /= nframes - 1)) THEN ! Condition if isoparametric coordinate at time jump
          isot = isot - 2._REKIND ! Reset isot by one slab (isoparametric length of 2.0)
          slabnum = slabnum + 1
          DO ii = 1, nnte ! Loop over time nodes per slab
               READ(UNIT=41, FMT=TRIM(fmt1), IOSTAT=errstat, IOMSG=errmesg) disphistory(:, ii) ! Read file
               IF(errstat/=0) WRITE(*,*) 'Error reading disphist at frame ', framenum,' after jump now in slab ', slabnum, &
                    & ' node ', ii, '. Code: ', errstat, '. Message: ', errmesg ! Error message
          END DO ! End loop over time nodes per slab

          CALL shape1x(ntframe, isot, nnte) ! Calculate temporal shape function at initial time/frame 0.0
          estress = 0.0_REKIND ! Initialize estress to 0.0
          weight  = 0.0_REKIND ! Initialize weight to 0.0
          xdisp   = 0.0_REKIND ! Initialize x displacement to 0.0
          xdisp   = 0.0_REKIND ! Initialize x displacement to 0.0
          zdisp   = 0.0_REKIND ! Initialize z displacement to 0.0
          DO ii = 1, ne ! Loop over elements
               nodes = xconn2(ii, :) ! Nodes for current element in loop  
               ncoorde = xcoord2(nodes, :) !Coordinates for nodes in current element in loop
               kk = 0 ! Initialize index
               DO jj = 1, nnse ! Loop over nodes per spatial element
                    kk = kk + 2 ! Iterate index
                    ll = 2*nodes(jj) ! Initialize index
                    dslab(kk - 1 : kk, :) = disphistory(ll - 1 : ll, :) ! Nodal displacement of slab, in m
               END DO ! End loop over nodes per spatial element
               dframe = MATMUL(dslab, ntframe) ! Displacement of current element in current frame
               ll = 0 ! Initialize index
               mm = 0 ! Initialize index
               DO jj = 1, 2 ! Loop over nodes of plotting elements lengthwise
                    DO kk = 1, 2 ! Loop over nodes of plotting elements heightwise
                         ll = ll + 2 ! Iterate index
                         mm = mm + 1 ! Iterate index
                         jac = MATMUL(dns(ll-1:ll, :), ncoorde) ! Deformation matrix
                         gam = inverse(jac, 2) ! Gamma = inverse of Deformation matrix
                         gamexpand = 0._REKIND ! Intialize gamexpand
                         gamexpand(1:2, 1:2) = gam ! Assign gam to gamexpand
                         gamexpand(3:4, 3:4) = gam ! Assign gam to gamexpand
                         bmat = MATMUL(gamexpand, dnxy(ll + ll - 3: ll + ll, :)) ! strain displacement matrix
                         sigmat = MATMUL(stiffsd, MATMUL(bmat, dframe)) ! Stress in Pa at time of frame and at node (ll/2)
                         sigmae(mm) = SQRT(sigmat(1)**2 - sigmat(1)*sigmat(2) + sigmat(2)**2 + 3._REKIND*sigmat(3)**2)
                         dplot(:,mm) = MATMUL(ns(ll-1:ll, :), dframe)
                    END DO ! End loop over nodes of plotting elements heightwise
               END DO ! End loop over nodes of plotting elements lengthwise
               DO jj = 1, 4 ! Loop over plotting nodes in element
                    ! Effective stress at plotting nodes (unweighted)
                    estress(xconn(ii,jj)) = estress(xconn(ii,jj)) + sigmae(jj)
                    ! Effective displacement at plotting nodes (unweighted)
                    xdisp(xconn(ii,jj)) = xdisp(xconn(ii,jj)) + dplot(1, jj)
                    ydisp(xconn(ii,jj)) = ydisp(xconn(ii,jj)) + dplot(2, jj)
                    ! Increment weight for eact time node that appears in the connectivity
                    weight(xconn(ii,jj)) = weight(xconn(ii,jj)) + 1._REKIND
               END DO ! End loop ove plotting nodes in element
          END DO ! End loop over elements
          DO ii = 1, nnhs*nnls ! Loop over plotting nodes
               xdisp(ii) = xdisp(ii) / weight(ii) ! X displacement, in m
               ydisp(ii) = ydisp(ii) / weight(ii) ! X displacement, in m
               estress(ii) = estress(ii) / weight(ii) ! Effective stress, in Pa
          END DO ! End loop over plotting nodes
          ! Add displaced after jump to displaced before jump:
          ! displaced is a matrix order (number of plotting nodes, 4).
          ! The 1st column is the x coordinate at frame jj
          ! The 2nd column is the y coordinate at frame jj
          ! The 3rd column is the z coordinate at frame jj
          ! The 4th column is the effective stress at frame jj
          displaced(:,1) = displaced(:,1) + xcoord(:,1) + xmag2 * xdisp
          displaced(:,2) = displaced(:,2) + xcoord(:,2) + ymag2 * ydisp
          displaced(:,3) = displaced(:,3) + xcoord(:,3) + zmag2 * zdisp
          displaced(:,4) = displaced(:,4) + estress
          displaced = 0.5_REKIND * displaced ! Average displaced from both sides of temporal jump
     END IF ! End check if frame at temporal jump
     
     ! Output the initial data for the PLT file:
     tframe = tframe + tframeincrement ! Time of current frame
     WRITE(dtc,*) tframe ! Write time of frame
     dtc = ADJUSTL(dtc) ! Adjust dtc to be flush left
     WRITE(UNIT=42,FMT='(1X,A,3(/,A))',IOSTAT=errstat, IOMSG=errmesg) 'ZONE STRANDID=1, SOLUTIONTIME='//TRIM(dtc), &
          & 'Nodes='//TRIM(nnc)//', Elements='//TRIM(nec)//', ZONETYPE=FEQuadrilateral', 'DATAPACKING=POINT', &
          & 'DT=(SINGLE SINGLE SINGLE SINGLE )'
     IF (errstat /= 0) WRITE(*,*) 'Error writing preamble to output file at ', framenum, '. Code: ', errstat, '. Message: ', errmesg

     ! Write out displaced to file
     DO ii = 1, nnhs*nnls ! Loop over plotting nodes
          WRITE(UNIT=42,FMT=TRIM(fmt2),IOSTAT=errstat, IOMSG=errmesg) displaced(ii,:)
          IF (errstat /= 0) WRITE(*,*) 'Error writing displaced to output file. Code: ', errstat, '. Message: ', errmesg
     END DO ! End loop over plotting nodes

     ! Write out nodal connectivity matrix to file
     DO ii = 1, ne ! Loop over plotting elements
          WRITE(UNIT=42,FMT=TRIM(fmt3),IOSTAT=errstat, IOMSG=errmesg) xconn(ii,:)
          IF (errstat /= 0) WRITE(*,*) 'Error writing xconn to output file. Code: ', errstat, '. Message: ', errmesg ! Error message
     END DO ! End loop over plotting elements

END DO ! End loop over frames

DEALLOCATE(displaced, estress, weight, xcoord, xdisp, ydisp, zdisp)!, STAT=errstat, ERRMSG=errmesg)
!IF(errstat/=0) WRITE(*,*) 'Error deallocating variables Location B. Code: ', errstat, '. Message: ', errmesg

ENDFILE(UNIT=42, IOSTAT=errstat, IOMSG=errmesg)
IF (errstat /= 0) WRITE(*,*) 'Error ending output file. Code: ', errstat, '. Message: ', errmesg ! Error message
REWIND(UNIT=42, IOSTAT=errstat, IOMSG=errmesg) ! Rewind file
IF (errstat /= 0) WRITE(*,*) 'Error rewinding output file. Code: ', errstat, '. Message: ', errmesg  ! Error message
CLOSE(UNIT=42,STATUS='KEEP',IOSTAT=errstat,IOMSG=errmesg) ! Close file
IF (errstat /= 0) WRITE(*,*) 'Error closing output file. Code: ', errstat, '. Message: ', errmesg ! Error message

REWIND(UNIT=41, IOSTAT=errstat, IOMSG=errmesg) ! Rewind file
IF (errstat /= 0) WRITE(*,*) 'Error rewinding diphist file. Code: ', errstat, '. Message: ', errmesg  ! Error message
CLOSE(UNIT=41,STATUS='KEEP',IOSTAT=errstat,IOMSG=errmesg) ! Close file
IF (errstat /= 0) WRITE(*,*) 'Error closing disphist file. Code: ', errstat, '. Message: ', errmesg ! Error message

END SUBROUTINE stressdispplot

SUBROUTINE timeelapsed(timea, timeb)
INTEGER, DIMENSION(8), INTENT(IN) :: timea, timeb
INTEGER, DIMENSION(8) :: timec
timec = timeb - timea
IF(timec(8) < 0) THEN ! Negative milliseconds
     timec(7) = timec(7) - 1
     timec(8) = timec(8) + 1000
END IF
IF(timec(7) < 0) THEN ! Negative seconds
     timec(6) = timec(6) - 1
     timec(7) = timec(7) + 60
END IF
IF(timec(6) < 0) THEN ! Negative minutes
     timec(5) = timec(5) - 1
     timec(6) = timec(6) + 60
END IF
IF(timec(5) < 0) THEN ! Negative hours
     timec(3) = timec(3) - 1
     timec(5) = timec(5) + 24
END IF
IF(timec(3) < 0) THEN
     timec(2) = timec(2) - 1
     timec(3) = timec(3) + 31
END IF
IF(timec(2) < 0) THEN
     WRITE(*,*) 'User Error: STOP working on New Years Eve, go celebrate!'
! Thank you David  -- Rui
END IF

!WRITE(*,FMT='(1X,"Time = ",1I2," days, ", 1I2, " hours, ", 1I2, " minutes, ", 1F6.3, " seconds")') &
!     & timec(3), timec(5), timec(6), REAL(timec(7)) + 0.001*REAL(timec(8))
WRITE(*,FMT='(1X,"Time = ", 1I2, " minutes, ", 1F6.3, " seconds")') &
     & timec(6), REAL(timec(7)) + 0.001*REAL(timec(8))

END SUBROUTINE timeelapsed

END MODULE PROCS

