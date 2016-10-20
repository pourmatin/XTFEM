!  interface_definitions.f03
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
!   11/04/2010     D.N.A.     Original Code to use old Lapack subroutines
!   10/21/2011     D.N.A.     Expanded interfaces to include SGL and DBL COMPLEX subroutines and CULA_GESV
!
!  All units should be converted to SI units.
!
!  Purpose: The interface_definitions MODULE allows for an explicit interface to the old implicit Lapack subroutines


MODULE INTERFACE_DEFINITIONS

USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
IMPLICIT NONE

!---------------------------------------------------------------

INTERFACE GESV

     SUBROUTINE SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     REAL(kind=SGL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if INFO = 0, the solution matrix X
     INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE SGESV

     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     REAL(kind=DBL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if INFO = 0, the solution matrix X
     INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE DGESV

     SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     COMPLEX(kind=SGL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if INFO = 0, the solution matrix X
     INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE CGESV

     SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     COMPLEX(kind=DBL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if INFO = 0, the solution matrix X
     INTEGER, DIMENSION(N), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE ZGESV

     MODULE PROCEDURE sgesv95
     MODULE PROCEDURE dgesv95
     MODULE PROCEDURE cgesv95
     MODULE PROCEDURE zgesv95

END INTERFACE GESV

!---------------------------------------------------------------

INTERFACE GETRF

     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER, DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE DGETRF
     
     SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER, DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE SGETRF

     SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER, DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE ZGETRF

     SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER, DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END SUBROUTINE CGETRF

     MODULE PROCEDURE dgetrf95
     MODULE PROCEDURE sgetrf95
     MODULE PROCEDURE zgetrf95
     MODULE PROCEDURE cgetrf95

END INTERFACE GETRF

!---------------------------------------------------------------

INTERFACE GETRI

     SUBROUTINE SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, LDA, LWORK
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize
     !         returned by ILAENV
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER, DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     REAL(kind=SGL), DIMENSION(MAX(1,LWORK)), INTENT(OUT) :: WORK
     ! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
     END SUBROUTINE SGETRI

     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
     USE kinds
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, LDA, LWORK
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize
     !         returned by ILAENV
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     REAL(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by DGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER, DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from DGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     REAL(kind=DBL), DIMENSION(MAX(1,LWORK)), INTENT(OUT) :: WORK
     ! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
     END SUBROUTINE DGETRI

     SUBROUTINE CGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
     USE kinds
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, LDA, LWORK
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize
     !         returned by ILAENV
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=SGL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by DGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER, DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from DGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     COMPLEX(kind=SGL), DIMENSION(MAX(1,LWORK)), INTENT(OUT) :: WORK
     ! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
     END SUBROUTINE CGETRI

     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
     USE kinds
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: N, LDA, LWORK
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize
     !         returned by ILAENV
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so
     !           the solution could not be computed
     COMPLEX(kind=DBL), DIMENSION(LDA, N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by DGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER, DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from DGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     COMPLEX(kind=DBL), DIMENSION(MAX(1,LWORK)), INTENT(OUT) :: WORK
     ! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
     END SUBROUTINE ZGETRI

     MODULE PROCEDURE sgetri95
     MODULE PROCEDURE dgetri95
     MODULE PROCEDURE cgetri95
     MODULE PROCEDURE zgetri95

END INTERFACE GETRI

!---------------------------------------------------------------

INTERFACE GEMM

     SUBROUTINE SGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     REAL(kind=SGL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     REAL(kind=SGL), DIMENSION(LDA, *), INTENT(IN) :: A
     REAL(kind=SGL), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     REAL(kind=SGL), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE SGEMM

     SUBROUTINE DGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     REAL(kind=DBL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     REAL(kind=DBL), DIMENSION(LDA, *), INTENT(IN) :: A
     REAL(kind=DBL), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     REAL(kind=DBL), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE DGEMM

     SUBROUTINE CGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     COMPLEX(kind=SGL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**H
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**H
     COMPLEX(kind=SGL), DIMENSION(LDA, *), INTENT(IN) :: A
     COMPLEX(kind=SGL), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     COMPLEX(kind=SGL), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE CGEMM

     SUBROUTINE ZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     COMPLEX(kind=DBL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**H
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**H
     COMPLEX(kind=DBL), DIMENSION(LDA, *), INTENT(IN) :: A
     COMPLEX(kind=DBL), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     COMPLEX(kind=DBL), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE ZGEMM

     MODULE PROCEDURE sgemm95
     MODULE PROCEDURE dgemm95
     MODULE PROCEDURE cgemm95
     MODULE PROCEDURE zgemm95

END INTERFACE GEMM

!---------------------------------------------------------------

INTERFACE GEMV

     SUBROUTINE SGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     REAL(kind=SGL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     REAL(kind=SGL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     REAL(kind=SGL), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     REAL(kind=SGL), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE SGEMV

     SUBROUTINE DGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     REAL(kind=DBL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     REAL(kind=DBL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     REAL(kind=DBL), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     REAL(kind=DBL), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE DGEMV

     SUBROUTINE CGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     COMPLEX(kind=SGL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**H * x + beta * y
     COMPLEX(kind=SGL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     COMPLEX(kind=SGL), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     COMPLEX(kind=SGL), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE CGEMV

     SUBROUTINE ZGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
     USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     COMPLEX(kind=DBL), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(len=1), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**H * x + beta * y
     COMPLEX(kind=DBL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     COMPLEX(kind=DBL), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     COMPLEX(kind=DBL), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE ZGEMV

     MODULE PROCEDURE sgemv95
     MODULE PROCEDURE dgemv95
     MODULE PROCEDURE cgemv95
     MODULE PROCEDURE zgemv95

END INTERFACE GEMV

!---------------------------------------------------------------

INTERFACE TRTRS

     SUBROUTINE STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
     USE kinds
     IMPLICIT NONE
     CHARACTER(len=1), INTENT(IN) :: UPLO, TRANS, DIAG
     ! UPLO = U: A is upper triangular
     ! UPLO = L: A is lower triangular
     ! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose)
     !                                                 T: A**T * X = B (Transpose)
     !                                                 C: A**H * X = B (Conjugate transpose = Transpose)
     ! DIAG = N: A is a non-unit triangular
     !      = U: A is unit triangular
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = Order of matrix A. N >= 0.
     ! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
     REAL(kind=SGL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! Triangular matrix A.
     ! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
     !              is not referenced.
     ! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
     !              is not referenced.
     ! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
     REAL(kind=SGL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, right hand side matrix.
     ! On exit, if INFO=0, solution matrix X
     END SUBROUTINE STRTRS

     SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
     USE kinds
     IMPLICIT NONE
     CHARACTER(len=1), INTENT(IN) :: UPLO, TRANS, DIAG
     ! UPLO = U: A is upper triangular
     ! UPLO = L: A is lower triangular
     ! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose)
     !                                                 T: A**T * X = B (Transpose)
     !                                                 C: A**H * X = B (Conjugate transpose = Transpose)
     ! DIAG = N: A is a non-unit triangular
     !      = U: A is unit triangular
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = Order of matrix A. N >= 0.
     ! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
     REAL(kind=DBL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! Triangular matrix A.
     ! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
     !              is not referenced.
     ! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
     !              is not referenced.
     ! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
     REAL(kind=DBL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, right hand side matrix.
     ! On exit, if INFO=0, solution matrix X
     END SUBROUTINE DTRTRS

     SUBROUTINE CTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
     USE kinds
     IMPLICIT NONE
     CHARACTER(len=1), INTENT(IN) :: UPLO, TRANS, DIAG
     ! UPLO = U: A is upper triangular
     ! UPLO = L: A is lower triangular
     ! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose)
     !                                                 T: A**T * X = B (Transpose)
     !                                                 C: A**H * X = B (Conjugate transpose)
     ! DIAG = N: A is a non-unit triangular
     !      = U: A is unit triangular
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = Order of matrix A. N >= 0.
     ! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
     COMPLEX(kind=SGL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! Triangular matrix A.
     ! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
     !              is not referenced.
     ! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
     !              is not referenced.
     ! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
     COMPLEX(kind=SGL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, right hand side matrix.
     ! On exit, if INFO=0, solution matrix X
     END SUBROUTINE CTRTRS

     SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
     USE kinds
     IMPLICIT NONE
     CHARACTER(len=1), INTENT(IN) :: UPLO, TRANS, DIAG
     ! UPLO = U: A is upper triangular
     ! UPLO = L: A is lower triangular
     ! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose)
     !                                                 T: A**T * X = B (Transpose)
     !                                                 C: A**H * X = B (Conjugate transpose)
     ! DIAG = N: A is a non-unit triangular
     !      = U: A is unit triangular
     INTEGER, INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = Order of matrix A. N >= 0.
     ! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     INTEGER, INTENT(OUT) :: INFO
     ! INFO = 0: successful exit
     !      < 0: if INFO = -i, i-th argument had an illegal value
     !      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
     COMPLEX(kind=DBL), DIMENSION(LDA, N), INTENT(IN) :: A
     ! Triangular matrix A.
     ! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
     !              is not referenced.
     ! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
     !              is not referenced.
     ! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
     COMPLEX(kind=DBL), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, right hand side matrix.
     ! On exit, if INFO=0, solution matrix X
     END SUBROUTINE ZTRTRS

     MODULE PROCEDURE strtrs95
     MODULE PROCEDURE dtrtrs95
     MODULE PROCEDURE ctrtrs95
     MODULE PROCEDURE ztrtrs95

END INTERFACE TRTRS

!---------------------------------------------------------------

INTERFACE cula_gesv

     FUNCTION cula_Sgesv(N, NRHS, A, LDA, IPIV, B, LDB) BIND(C, NAME='CULA_SGESV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Sgesv
     INTEGER(C_INT), INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     REAL(C_FLOAT), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     REAL(C_FLOAT), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if culaNoError is returned, the solution matrix X
     INTEGER(C_INT), DIMENSION(N), INTENT(OUT) :: IPIV
     END FUNCTION cula_Sgesv
     
     FUNCTION cula_Dgesv(N, NRHS, A, LDA, IPIV, B, LDB) BIND(C, NAME='CULA_DGESV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Dgesv
     INTEGER(C_INT), INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     REAL(C_DOUBLE), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     REAL(C_DOUBLE), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if culaNoError is returned, the solution matrix X
     INTEGER(C_INT), DIMENSION(N), INTENT(OUT) :: IPIV
     END FUNCTION cula_Dgesv

     FUNCTION cula_Cgesv(N, NRHS, A, LDA, IPIV, B, LDB) BIND(C, NAME='CULA_CGESV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Cgesv
     INTEGER(C_INT), INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if culaNoError is returned, the solution matrix X
     INTEGER(C_INT), DIMENSION(N), INTENT(OUT) :: IPIV
     END FUNCTION cula_Cgesv
     
     FUNCTION cula_Zgesv(N, NRHS, A, LDA, IPIV, B, LDB) BIND(C, NAME='CULA_ZGESV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Zgesv
     INTEGER(C_INT), INTENT(IN) :: N, NRHS, LDA, LDB
     ! N = # of linear equations (order of matrix A). N >= 0.
     ! NRHS = # of right hand sides (columns of B). NRHS >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     ! LDB = Leading dimension of array B. LDB >= max(1,N).
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDB, NRHS), INTENT(INOUT) :: B
     ! On entry, the right hand side matrix
     ! On exit, if culaNoError is returned, the solution matrix X
     INTEGER(C_INT), DIMENSION(N), INTENT(OUT) :: IPIV
     END FUNCTION cula_Zgesv

END INTERFACE cula_gesv

!---------------------------------------------------------------

INTERFACE cula_getrf

     FUNCTION cula_Sgetrf(M, N, A, LDA, IPIV) BIND(C, NAME='CULA_SGETRF')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Sgetrf
     INTEGER(C_INT), INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     REAL(C_FLOAT), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER(C_INT), DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Sgetrf
     
     FUNCTION cula_Dgetrf(M, N, A, LDA, IPIV) BIND(C, NAME='CULA_DGETRF')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Dgetrf
     INTEGER(C_INT), INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     REAL(C_DOUBLE), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER(C_INT), DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Dgetrf

     FUNCTION cula_Cgetrf(M, N, A, LDA, IPIV) BIND(C, NAME='CULA_CGETRF')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Cgetrf
     INTEGER(C_INT), INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER(C_INT), DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Cgetrf

     FUNCTION cula_Zgetrf(M, N, A, LDA, IPIV) BIND(C, NAME='CULA_ZGETRF')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Zgetrf
     INTEGER(C_INT), INTENT(IN) :: M, N, LDA
     ! M = # rows of matrix A. M >= 0.
     ! N = # columns of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, MxN coefficient matrix A
     ! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
     INTEGER(C_INT), DIMENSION(MIN(M,N)), INTENT(OUT) :: IPIV
     ! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Zgetrf

END INTERFACE cula_getrf

!---------------------------------------------------------------
INTERFACE cula_getri

     FUNCTION cula_Sgetri(N, A, LDA, IPIV) BIND(C, NAME='CULA_SGETRI')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Sgetri
     INTEGER(C_INT), INTENT(IN) :: N, LDA
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     REAL(C_FLOAT), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER(C_INT), DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Sgetri

     FUNCTION cula_Dgetri(N, A, LDA, IPIV) BIND(C, NAME='CULA_DGETRI')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Dgetri
     INTEGER(C_INT), INTENT(IN) :: N, LDA
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     REAL(C_DOUBLE), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER(C_INT), DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Dgetri

     FUNCTION cula_Cgetri(N, A, LDA, IPIV) BIND(C, NAME='CULA_CGETRI')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Cgetri
     INTEGER(C_INT), INTENT(IN) :: N, LDA
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER(C_INT), DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Cgetri

     FUNCTION cula_Zgetri(N, A, LDA, IPIV) BIND(C, NAME='CULA_ZGETRI')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_Zgetri
     INTEGER(C_INT), INTENT(IN) :: N, LDA
     ! N = Order of matrix A. N >= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,N).
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDA,N), INTENT(INOUT) :: A
     ! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
     ! On exit, if info = 0, the inverse of the original matrix A
     INTEGER(C_INT), DIMENSION(N), INTENT(IN) :: IPIV
     ! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
     END FUNCTION cula_Zgetri

END INTERFACE cula_getri

!---------------------------------------------------------------

INTERFACE cula_gemm

     SUBROUTINE cula_Sgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) BIND(C, NAME='CULA_SGEMM')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     REAL(C_FLOAT), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     REAL(C_FLOAT), DIMENSION(LDA, *), INTENT(IN) :: A
     REAL(C_FLOAT), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     REAL(C_FLOAT), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE cula_Sgemm

     SUBROUTINE cula_Dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) BIND(C, NAME='CULA_DGEMM')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     REAL(C_DOUBLE), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     REAL(C_DOUBLE), DIMENSION(LDA, *), INTENT(IN) :: A
     REAL(C_DOUBLE), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     REAL(C_DOUBLE), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE cula_Dgemm

     SUBROUTINE cula_Cgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) BIND(C, NAME='CULA_CGEMM')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDA, *), INTENT(IN) :: A
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE cula_Cgemm

     SUBROUTINE cula_Zgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) BIND(C, NAME='CULA_ZGEMM')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: K, LDA, LDB, LDC, M, N
     ! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
     ! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
     ! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
     ! LDC = Leading dimension of array C. LDC >= max(1,M).
     ! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
     ! N = # of columns of matrix op( B) and of matrix C. N >= 0.
     COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then C need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANSA, TRANSB
     ! TRANSA = Form of op( A ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( A ) = A
     !        = 'T' or 't', op( A ) = A**T
     !        = 'C' or 'c', op( A ) = A**T
     ! TRANSB = Form of op( B ) to be used in matrix multiplication:
     !        = 'N' or 'n', op( B ) = B
     !        = 'T' or 't', op( B ) = B**T
     !        = 'C' or 'c', op( B ) = B**T
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDA, *), INTENT(IN) :: A
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDB, *), INTENT(IN) :: B
     ! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
     ! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
     ! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
     ! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDC, N), INTENT(INOUT) :: C
     ! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
     ! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
     END SUBROUTINE cula_Zgemm

END INTERFACE cula_gemm

!---------------------------------------------------------------

INTERFACE cula_gemv

     SUBROUTINE cula_Sgemv( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) BIND(C, NAME='CULA_SGEMV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     REAL(C_FLOAT), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     REAL(C_FLOAT), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     REAL(C_FLOAT), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     REAL(C_FLOAT), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE cula_Sgemv

     SUBROUTINE cula_Dgemv( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) BIND(C, NAME='CULA_DGEMV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     REAL(C_DOUBLE), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     REAL(C_DOUBLE), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     REAL(C_DOUBLE), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE cula_Dgemv

     SUBROUTINE cula_Cgemv( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) BIND(C, NAME='CULA_CGEMV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     COMPLEX(C_FLOAT_COMPLEX), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     COMPLEX(C_FLOAT_COMPLEX), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE cula_Cgemv

     SUBROUTINE cula_Zgemv( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) BIND(C, NAME='CULA_ZGEMV')
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT), INTENT(IN) :: INCX, INCY, LDA, M, N
     ! INCX = Increment for elements of X. INCX /= 0.
     ! INCY = Increment for elements of Y. INCY /= 0.
     ! LDA = Leading dimension of array A. LDA >= max(1,M).
     ! M = # of rows of matrix A. M >= 0.
     ! N = # of columns of matrix A. N >= 0.
     COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: ALPHA, BETA
     ! ALPHA = scalar alpha.
     ! BETA = scalar beta. When BETA is 0 then Y need not be set on input.
     CHARACTER(kind=C_CHAR), INTENT(IN) :: TRANS
     ! TRANS = Operation to be performed:
     !        = 'N' or 'n', y = alpha * A * x + beta * y
     !        = 'T' or 't', y = alpha * A**T * x + beta * y
     !        = 'C' or 'c', y = alpha * A**T * x + beta * y
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(LDA, N), INTENT(IN) :: A
     ! The leading M by N part contains the matrix of coefficients A.
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*), INTENT(IN) :: X
     ! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
     ! DIMENSION is at least (1+(M-1)*abs(INCX)).
     COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(*), INTENT(INOUT) :: Y
     ! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
     ! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).
     END SUBROUTINE cula_Zgemv

END INTERFACE cula_gemv

!---------------------------------------------------------------

INTERFACE
     FUNCTION cula_initialize()
     USE ISO_C_BINDING
     IMPLICIT NONE
     INTEGER(C_INT) :: cula_initialize
     END FUNCTION cula_initialize
     
     SUBROUTINE cula_shutdown()
     USE ISO_C_BINDING
     IMPLICIT NONE
     END SUBROUTINE cula_shutdown
END INTERFACE

!---------------------------------------------------------------

CONTAINS

!---------------------------------------------------------------

SUBROUTINE dgesv95(A, B, ipiv, info)
IMPLICIT NONE
INTEGER, INTENT(OUT), OPTIONAL :: info
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: info2 ! Dummy variable for info
INTEGER :: n, nrhs, lda, ldb
! N = # of linear equations (order of matrix A). N >= 0.
! NRHS = # of right hand sides (columns of B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
REAL(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
REAL(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: B
! On entry, the right hand side matrix
! On exit, if INFO = 0, the solution matrix X
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: ipiv
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv2
CHARACTER(len=120) :: errmesg
lda  = SIZE(A,1)
n    = SIZE(A,2)
ldb  = SIZE(B,1)
nrhs = SIZE(B,2)
ALLOCATE(ipiv2(n))
!ALLOCATE(ipiv2(n), STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not ALLOCATE correctly in DGESV95. The error message is: ', errmesg
CALL dgesv( n, nrhs, A, lda, ipiv2, B, ldb, info2 )
IF(PRESENT(ipiv)) ipiv = ipiv2
IF(PRESENT(info)) info = info2
DEALLOCATE(ipiv2)
!DEALLOCATE(ipiv2, STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not DEALLOCATE correctly in DGESV95. The error message is: ', errmesg
END SUBROUTINE dgesv95

SUBROUTINE sgesv95(A, B, ipiv, info)
IMPLICIT NONE
INTEGER, INTENT(OUT), OPTIONAL :: info
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: info2 ! Dummy variable for info
INTEGER :: n, nrhs, lda, ldb
! N = # of linear equations (order of matrix A). N >= 0.
! NRHS = # of right hand sides (columns of B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
REAL(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
REAL(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: B
! On entry, the right hand side matrix
! On exit, if INFO = 0, the solution matrix X
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: ipiv
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv2
CHARACTER(len=120) :: errmesg
lda  = SIZE(A,1)
n    = SIZE(A,2)
ldb  = SIZE(B,1)
nrhs = SIZE(B,2)
ALLOCATE(ipiv2(n))
!ALLOCATE(ipiv2(n),STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not ALLOCATE correctly in SGESV95. The error message is: ', errmesg
CALL sgesv( n, nrhs, A, lda, ipiv2, B, ldb, info2 )
IF(PRESENT(ipiv)) ipiv = ipiv2
IF(PRESENT(info)) info = info2
DEALLOCATE(ipiv2)
!DEALLOCATE(ipiv2,  STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not DEALLOCATE correctly in SGESV95. The error message is: ', errmesg
END SUBROUTINE sgesv95

SUBROUTINE zgesv95(A, B, ipiv, info)
IMPLICIT NONE
INTEGER, INTENT(OUT), OPTIONAL :: info
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: info2 ! Dummy variable for info
INTEGER :: n, nrhs, lda, ldb
! N = # of linear equations (order of matrix A). N >= 0.
! NRHS = # of right hand sides (columns of B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
COMPLEX(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
COMPLEX(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: B
! On entry, the right hand side matrix
! On exit, if INFO = 0, the solution matrix X
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: ipiv
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv2
CHARACTER(len=120) :: errmesg
lda  = SIZE(A,1)
n    = SIZE(A,2)
ldb  = SIZE(B,1)
nrhs = SIZE(B,2)
ALLOCATE(ipiv2(n))
!ALLOCATE(ipiv2(n), STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not ALLOCATE correctly in DGESV95. The error message is: ', errmesg
CALL zgesv( n, nrhs, A, lda, ipiv2, B, ldb, info2 )
IF(PRESENT(ipiv)) ipiv = ipiv2
IF(PRESENT(info)) info = info2
DEALLOCATE(ipiv2)
!DEALLOCATE(ipiv2, STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not DEALLOCATE correctly in DGESV95. The error message is: ', errmesg
END SUBROUTINE zgesv95

SUBROUTINE cgesv95(A, B, ipiv, info)
IMPLICIT NONE
INTEGER, INTENT(OUT), OPTIONAL :: info
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: info2 ! Dummy variable for info
INTEGER :: n, nrhs, lda, ldb
! N = # of linear equations (order of matrix A). N >= 0.
! NRHS = # of right hand sides (columns of B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
COMPLEX(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
COMPLEX(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: B
! On entry, the right hand side matrix
! On exit, if INFO = 0, the solution matrix X
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: ipiv
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv2
CHARACTER(len=120) :: errmesg
lda  = SIZE(A,1)
n    = SIZE(A,2)
ldb  = SIZE(B,1)
nrhs = SIZE(B,2)
ALLOCATE(ipiv2(n))
!ALLOCATE(ipiv2(n), STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not ALLOCATE correctly in DGESV95. The error message is: ', errmesg
CALL cgesv( n, nrhs, A, lda, ipiv2, B, ldb, info2 )
IF(PRESENT(ipiv)) ipiv = ipiv2
IF(PRESENT(info)) info = info2
DEALLOCATE(ipiv2)
!DEALLOCATE(ipiv2, STAT = info2, ERRMSG = errmesg)
!IF(info2 /= 0) WRITE(*,*) 'IPIV and B did not DEALLOCATE correctly in DGESV95. The error message is: ', errmesg
END SUBROUTINE cgesv95

!---------------------------------------------------------------

SUBROUTINE dgetrf95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: M, N, LDA
! M = # rows of matrix A. M >= 0.
! N = # columns of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
REAL(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, MxN coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: IPIV
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV2
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
M   = LDA
ALLOCATE(IPIV2(MIN(M,N)))
!ALLOCATE(IPIV2(MIN(M,N)), STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not ALLOCATE correctly in DGETRF95. The error message is: ', errmesg
CALL dgetrf( M, N, A, LDA, IPIV2, INFO2 )
IF(PRESENT(IPIV)) IPIV = IPIV2
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(IPIV2)
!DEALLOCATE(IPIV2, STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not DEALLOCATE correctly in DGETRF95. The error message is: ', errmesg
END SUBROUTINE dgetrf95
     
SUBROUTINE sgetrf95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: M, N, LDA
! M = # rows of matrix A. M >= 0.
! N = # columns of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
REAL(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, MxN coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: IPIV
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV2
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
M   = LDA
ALLOCATE(IPIV2(MIN(M,N)))
!ALLOCATE(IPIV2(MIN(M,N)), STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not ALLOCATE correctly in SGETRF95. The error message is: ', errmesg
CALL sgetrf( M, N, A, LDA, IPIV2, INFO2 )
IF(PRESENT(IPIV)) IPIV = IPIV2
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(IPIV2)
!DEALLOCATE(IPIV2, STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not DEALLOCATE correctly in SGETRF95. The error message is: ', errmesg
END SUBROUTINE sgetrf95

SUBROUTINE zgetrf95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: M, N, LDA
! M = # rows of matrix A. M >= 0.
! N = # columns of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
COMPLEX(kind=DBL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, MxN coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: IPIV
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV2
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
M   = LDA
ALLOCATE(IPIV2(MIN(M,N)))
!ALLOCATE(IPIV2(MIN(M,N)), STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not ALLOCATE correctly in ZGETRF95. The error message is: ', errmesg
CALL zgetrf( M, N, A, LDA, IPIV2, INFO2 )
IF(PRESENT(IPIV)) IPIV = IPIV2
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(IPIV2)
!DEALLOCATE(IPIV2, STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not DEALLOCATE correctly in ZGETRF95. The error message is: ', errmesg
END SUBROUTINE zgetrf95

SUBROUTINE cgetrf95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: M, N, LDA
! M = # rows of matrix A. M >= 0.
! N = # columns of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
COMPLEX(kind=SGL), DIMENSION(:,:), INTENT(INOUT) :: A
! On entry, MxN coefficient matrix A
! On exit, the factors L and U from the factorization A = P*L*U; the unit diagonal elements of L are not stored
INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: IPIV
! Indices that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV2
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
M   = LDA
ALLOCATE(IPIV2(MIN(M,N)))
!ALLOCATE(IPIV2(MIN(M,N)), STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not ALLOCATE correctly in ZGETRF95. The error message is: ', errmesg
CALL cgetrf( M, N, A, LDA, IPIV2, INFO2 )
IF(PRESENT(IPIV)) IPIV = IPIV2
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(IPIV2)
!DEALLOCATE(IPIV2, STAT = INFO2, ERRMSG = errmesg)
!IF(INFO2 /= 0) WRITE(*,*) 'IPIV did not DEALLOCATE correctly in ZGETRF95. The error message is: ', errmesg
END SUBROUTINE cgetrf95

!---------------------------------------------------------------

SUBROUTINE sgemm95( A, B, C, TRANSA, TRANSB, ALPHA, BETA )
IMPLICIT NONE
INTEGER :: K, LDA, LDB, LDC, M, N
! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
! LDC = Leading dimension of array C. LDC >= max(1,M).
! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
! N = # of columns of matrix op( B) and of matrix C. N >= 0.
REAL(kind=SGL), INTENT(IN), OPTIONAL :: ALPHA, BETA
! ALPHA = scalar alpha. Default value = 1.0
! BETA = scalar beta. When BETA is 0 then C need not be set on input. Default value = 0.0
REAL(kind=SGL) :: ALPHA2, BETA2
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANSA, TRANSB
! TRANSA = Form of op( A ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( A ) = A
!        = 'T' or 't', op( A ) = A**T
!        = 'C' or 'c', op( A ) = A**T
! TRANSB = Form of op( B ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( B ) = B
!        = 'T' or 't', op( B ) = B**T
!        = 'C' or 'c', op( B ) = B**T
CHARACTER(len=1) :: TRANSA2, TRANSB2
REAL(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A, B
! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
REAL(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: C
! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
IF(PRESENT(TRANSA)) THEN
     TRANSA2 = TRANSA
ELSE
     TRANSA2 = 'N'
END IF
IF(PRESENT(TRANSB)) THEN
     TRANSB2 = TRANSB
ELSE
     TRANSB2 = 'N'
END IF
IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = 1.0_SGL
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = 0.0_SGL
END IF
LDA = SIZE(A,1)
LDB = SIZE(B,1)
LDC = SIZE(C,1)
N   = SIZE(C,2)
IF((TRANSA2 == 'N').OR.(TRANSA2 == 'n')) THEN
     K = SIZE(A,2)
     M = LDA
ELSE
     M = SIZE(A,2)
     K = LDA
END IF
CALL sgemm( TRANSA2, TRANSB2, M, N, K, ALPHA2, A, LDA, B, LDB, BETA2, C, LDC )
END SUBROUTINE sgemm95

SUBROUTINE dgemm95( A, B, C, TRANSA, TRANSB, ALPHA, BETA )
IMPLICIT NONE
INTEGER :: K, LDA, LDB, LDC, M, N
! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
! LDC = Leading dimension of array C. LDC >= max(1,M).
! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
! N = # of columns of matrix op( B) and of matrix C. N >= 0.
REAL(kind=DBL), INTENT(IN), OPTIONAL :: ALPHA, BETA
! ALPHA = scalar alpha. Default value = 1.0
! BETA = scalar beta. When BETA is 0 then C need not be set on input. Default value = 0.0
REAL(kind=DBL) :: ALPHA2, BETA2
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANSA, TRANSB
! TRANSA = Form of op( A ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( A ) = A
!        = 'T' or 't', op( A ) = A**T
!        = 'C' or 'c', op( A ) = A**T
! TRANSB = Form of op( B ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( B ) = B
!        = 'T' or 't', op( B ) = B**T
!        = 'C' or 'c', op( B ) = B**T
CHARACTER(len=1) :: TRANSA2, TRANSB2
REAL(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A, B
! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
REAL(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: C
! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
IF(PRESENT(TRANSA)) THEN
     TRANSA2 = TRANSA
ELSE
     TRANSA2 = 'N'
END IF
IF(PRESENT(TRANSB)) THEN
     TRANSB2 = TRANSB
ELSE
     TRANSB2 = 'N'
END IF
IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = 1.0_DBL
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = 0.0_DBL
END IF
LDA = SIZE(A,1)
LDB = SIZE(B,1)
LDC = SIZE(C,1)
N   = SIZE(C,2)
IF((TRANSA2 == 'N').OR.(TRANSA2 == 'n')) THEN
     K = SIZE(A,2)
     M = LDA
ELSE
     M = SIZE(A,2)
     K = LDA
END IF
CALL dgemm( TRANSA2, TRANSB2, M, N, K, ALPHA2, A, LDA, B, LDB, BETA2, C, LDC )
END SUBROUTINE dgemm95

SUBROUTINE cgemm95( A, B, C, TRANSA, TRANSB, ALPHA, BETA )
IMPLICIT NONE
INTEGER :: K, LDA, LDB, LDC, M, N
! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
! LDC = Leading dimension of array C. LDC >= max(1,M).
! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
! N = # of columns of matrix op( B) and of matrix C. N >= 0.
COMPLEX(kind=SGL), INTENT(IN), OPTIONAL :: ALPHA, BETA
! ALPHA = scalar alpha. Default value = 1.0
! BETA = scalar beta. When BETA is 0 then C need not be set on input. Default value = 0.0
COMPLEX(kind=SGL) :: ALPHA2, BETA2
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANSA, TRANSB
! TRANSA = Form of op( A ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( A ) = A
!        = 'T' or 't', op( A ) = A**T
!        = 'C' or 'c', op( A ) = A**T
! TRANSB = Form of op( B ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( B ) = B
!        = 'T' or 't', op( B ) = B**T
!        = 'C' or 'c', op( B ) = B**T
CHARACTER(len=1) :: TRANSA2, TRANSB2
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A, B
! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: C
! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
! On exit, C is overwritten by the M byd N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
IF(PRESENT(TRANSA)) THEN
     TRANSA2 = TRANSA
ELSE
     TRANSA2 = 'N'
END IF
IF(PRESENT(TRANSB)) THEN
     TRANSB2 = TRANSB
ELSE
     TRANSB2 = 'N'
END IF
IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = ( 1.0_SGL, 0.0_SGL )
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = ( 0.0_SGL, 0.0_SGL )
END IF
LDA = SIZE(A,1)
LDB = SIZE(B,1)
LDC = SIZE(C,1)
N   = SIZE(C,2)
IF((TRANSA2 == 'N').OR.(TRANSA2 == 'n')) THEN
     K = SIZE(A,2)
     M = LDA
ELSE
     M = SIZE(A,2)
     K = LDA
END IF
CALL cgemm( TRANSA2, TRANSB2, M, N, K, ALPHA2, A, LDA, B, LDB, BETA2, C, LDC )
END SUBROUTINE cgemm95

SUBROUTINE zgemm95( A, B, C, TRANSA, TRANSB, ALPHA, BETA )
IMPLICIT NONE
INTEGER :: K, LDA, LDB, LDC, M, N
! K = # of columns of matrix op( A ) and # of rows of matrix op( B ). K >= 0.
! LDA = Leading dimension of array A. When TRANSA == 'N' or 'n', LDA >= max(1,M), otherwise LDA >= max(1,K).
! LDB = Leading dimension of array B. When TRANSB == 'N' or 'n', LDB >= max(1,K), otherwise LDB >= max(1,N).
! LDC = Leading dimension of array C. LDC >= max(1,M).
! M = # of rows of matrix op( A ) and of matrix C. M >= 0.
! N = # of columns of matrix op( B) and of matrix C. N >= 0.
COMPLEX(kind=DBL), INTENT(IN), OPTIONAL :: ALPHA, BETA
! ALPHA = scalar alpha. Default value = 1.0
! BETA = scalar beta. When BETA is 0 then C need not be set on input. Default value = 0.0
COMPLEX(kind=DBL) :: ALPHA2, BETA2
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANSA, TRANSB
! TRANSA = Form of op( A ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( A ) = A
!        = 'T' or 't', op( A ) = A**T
!        = 'C' or 'c', op( A ) = A**T
! TRANSB = Form of op( B ) to be used in matrix multiplication. Default value = 'N'
!        = 'N' or 'n', op( B ) = B
!        = 'T' or 't', op( B ) = B**T
!        = 'C' or 'c', op( B ) = B**T
CHARACTER(len=1) :: TRANSA2, TRANSB2
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A, B
! With TRANSA == 'N' or 'n', A is of DIMENSION (LDA, K) and the leading M by K part contains the matrix A, otherwise A is of
! DIMENSION(LDA, M) and the leading K by M part contains the matrix A.
! With TRANSB == 'N' or 'n', B is of DIMENSION (LDB, N) and the leading K by N part contains the matrix B, otherwise B is of
! DIMENSION(LDB, K) and the leading N by K part contains the matrix B.
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: C
! On entry, the leading M by N part contains the matrix C, except when BETA is 0, in which case C need not be set.
! On exit, C is overwritten by the M by N matrix ( ALPHA * op( A ) * op( B ) + BETA * C).
IF(PRESENT(TRANSA)) THEN
     TRANSA2 = TRANSA
ELSE
     TRANSA2 = 'N'
END IF
IF(PRESENT(TRANSB)) THEN
     TRANSB2 = TRANSB
ELSE
     TRANSB2 = 'N'
END IF
IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = ( 1.0_DBL, 0.0_DBL )
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = ( 0.0_DBL, 0.0_DBL )
END IF
LDA = SIZE(A,1)
LDB = SIZE(B,1)
LDC = SIZE(C,1)
N   = SIZE(C,2)
IF((TRANSA2 == 'N').OR.(TRANSA2 == 'n')) THEN
     K = SIZE(A,2)
     M = LDA
ELSE
     M = SIZE(A,2)
     K = LDA
END IF
CALL zgemm( TRANSA2, TRANSB2, M, N, K, ALPHA2, A, LDA, B, LDB, BETA2, C, LDC )
END SUBROUTINE zgemm95

!---------------------------------------------------------------

SUBROUTINE sgemv95( A, X, Y, ALPHA, BETA, TRANS )
USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
IMPLICIT NONE
INTEGER :: M, N
! INCX = Increment for elements of X. INCX /= 0. Set to 1 for Fortran95 implementation.
! INCY = Increment for elements of Y. INCY /= 0. Set to 1 for Fortran95 implementation.
! LDA = Leading dimension of array A. LDA >= max(1,M). Replaced by just using M.
! M = # of rows of matrix A. M >= 0.
! N = # of columns of matrix A. N >= 0.
REAL(kind=SGL), INTENT(IN), OPTIONAL :: ALPHA, BETA
REAL(kind=SGL) :: ALPHA2, BETA2
! ALPHA = scalar alpha. The default value is 1.
! BETA = scalar beta. When BETA is 0 then Y need not be set on input. The default value is 0.
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANS
CHARACTER(len=1) :: TRANS2
! TRANS = Operation to be performed:
!        = 'N' or 'n', y = alpha * A * x + beta * y. This is the default value.
!        = 'T' or 't', y = alpha * A**T * x + beta * y
!        = 'C' or 'c', y = alpha * A**T * x + beta * y
REAL(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A
! The leading M by N part contains the matrix of coefficients A.
REAL(kind=SGL), DIMENSION(:), INTENT(IN) :: X
! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
! DIMENSION is at least (1+(M-1)*abs(INCX)).
REAL(kind=SGL), DIMENSION(:), INTENT(INOUT) :: Y
! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).

IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = 1._SGL
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = 0._SGL
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
M = SIZE(A,1)
N = SIZE(A,2)
CALL sgemv( TRANS2, M, N, ALPHA2, A, M, X, 1, BETA2, Y, 1 )
END SUBROUTINE sgemv95

SUBROUTINE dgemv95( A, X, Y, ALPHA, BETA, TRANS )
USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
IMPLICIT NONE
INTEGER :: M, N
! INCX = Increment for elements of X. INCX /= 0. Set to 1 for Fortran95 implementation.
! INCY = Increment for elements of Y. INCY /= 0. Set to 1 for Fortran95 implementation.
! LDA = Leading dimension of array A. LDA >= max(1,M). Replaced by just using M.
! M = # of rows of matrix A. M >= 0.
! N = # of columns of matrix A. N >= 0.
REAL(kind=DBL), INTENT(IN), OPTIONAL :: ALPHA, BETA
REAL(kind=DBL) :: ALPHA2, BETA2
! ALPHA = scalar alpha. The default value is 1.
! BETA = scalar beta. When BETA is 0 then Y need not be set on input. The default value is 0.
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANS
CHARACTER(len=1) :: TRANS2
! TRANS = Operation to be performed:
!        = 'N' or 'n', y = alpha * A * x + beta * y. This is the default value.
!        = 'T' or 't', y = alpha * A**T * x + beta * y
!        = 'C' or 'c', y = alpha * A**T * x + beta * y
REAL(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A
! The leading M by N part contains the matrix of coefficients A.
REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: X
! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
! DIMENSION is at least (1+(M-1)*abs(INCX)).
REAL(kind=DBL), DIMENSION(:), INTENT(INOUT) :: Y
! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).

IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = 1._DBL
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = 0._DBL
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
M = SIZE(A,1)
N = SIZE(A,2)
CALL dgemv( TRANS2, M, N, ALPHA2, A, M, X, 1, BETA2, Y, 1 )
END SUBROUTINE dgemv95

SUBROUTINE cgemv95( A, X, Y, ALPHA, BETA, TRANS )
USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
IMPLICIT NONE
INTEGER :: M, N
! INCX = Increment for elements of X. INCX /= 0. Set to 1 for Fortran95 implementation.
! INCY = Increment for elements of Y. INCY /= 0. Set to 1 for Fortran95 implementation.
! LDA = Leading dimension of array A. LDA >= max(1,M). Replaced by just using M.
! M = # of rows of matrix A. M >= 0.
! N = # of columns of matrix A. N >= 0.
COMPLEX(kind=SGL), INTENT(IN), OPTIONAL :: ALPHA, BETA
COMPLEX(kind=SGL) :: ALPHA2, BETA2
! ALPHA = scalar alpha. The default value is 1.
! BETA = scalar beta. When BETA is 0 then Y need not be set on input. The default value is 0.
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANS
CHARACTER(len=1) :: TRANS2
! TRANS = Operation to be performed:
!        = 'N' or 'n', y = alpha * A * x + beta * y. This is the default value.
!        = 'T' or 't', y = alpha * A**T * x + beta * y
!        = 'C' or 'c', y = alpha * A**H * x + beta * y
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A
! The leading M by N part contains the matrix of coefficients A.
COMPLEX(kind=SGL), DIMENSION(:), INTENT(IN) :: X
! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
! DIMENSION is at least (1+(M-1)*abs(INCX)).
COMPLEX(kind=SGL), DIMENSION(:), INTENT(INOUT) :: Y
! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).

IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = ( 1.0_SGL, 0.0_SGL )
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = ( 0.0_SGL, 0.0_SGL )
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
M = SIZE(A,1)
N = SIZE(A,2)
CALL cgemv( TRANS2, M, N, ALPHA2, A, M, X, 1, BETA2, Y, 1 )
END SUBROUTINE cgemv95

SUBROUTINE zgemv95( A, X, Y, ALPHA, BETA, TRANS )
USE kinds ! Access kinds module to specifiy kind of variables (single and double reals)
IMPLICIT NONE
INTEGER :: M, N
! INCX = Increment for elements of X. INCX /= 0. Set to 1 for Fortran95 implementation.
! INCY = Increment for elements of Y. INCY /= 0. Set to 1 for Fortran95 implementation.
! LDA = Leading dimension of array A. LDA >= max(1,M). Replaced by just using M.
! M = # of rows of matrix A. M >= 0.
! N = # of columns of matrix A. N >= 0.
COMPLEX(kind=DBL), INTENT(IN), OPTIONAL :: ALPHA, BETA
COMPLEX(kind=DBL) :: ALPHA2, BETA2
! ALPHA = scalar alpha. The default value is 1.
! BETA = scalar beta. When BETA is 0 then Y need not be set on input. The default value is 0.
CHARACTER(len=1), INTENT(IN), OPTIONAL :: TRANS
CHARACTER(len=1) :: TRANS2
! TRANS = Operation to be performed:
!        = 'N' or 'n', y = alpha * A * x + beta * y. This is the default value.
!        = 'T' or 't', y = alpha * A**T * x + beta * y
!        = 'C' or 'c', y = alpha * A**H * x + beta * y
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A
! The leading M by N part contains the matrix of coefficients A.
COMPLEX(kind=DBL), DIMENSION(:), INTENT(IN) :: X
! X = Incremented array X containing vector x. When TRANS == 'N' or 'n', DIMENSION is at least (1+(n-1)*abs(INCX)). Otherwise
! DIMENSION is at least (1+(M-1)*abs(INCX)).
COMPLEX(kind=DBL), DIMENSION(:), INTENT(INOUT) :: Y
! Before entry with with BETA /= 0 the incremented array Y contains the vector y. On exit, Y is overwritten by updated vector
! y. When TRANS = 'N' or 'n', DIMENSION is at least (1+(M-1)*abs(INCY)). Otherwise, DIMENSION is at least (1+(N-1)*abs(INCY)).

IF(PRESENT(ALPHA)) THEN
     ALPHA2 = ALPHA
ELSE
     ALPHA2 = ( 1.0_DBL, 0.0_DBL )
END IF
IF(PRESENT(BETA)) THEN
     BETA2 = BETA
ELSE
     BETA2 = ( 0.0_DBL, 0.0_DBL )
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
M = SIZE(A,1)
N = SIZE(A,2)
CALL zgemv( TRANS2, M, N, ALPHA2, A, M, X, 1, BETA2, Y, 1 )
END SUBROUTINE zgemv95

!---------------------------------------------------------------

SUBROUTINE sgetri95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: N, LDA, LWORK
! N = Order of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize returned
!         by ILAENV (64)
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
REAL(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: A
! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
! On exit, if info = 0, the inverse of the original matrix A
INTEGER, DIMENSION(:), INTENT(IN) :: IPIV
! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
REAL(kind=SGL), DIMENSION(:), ALLOCATABLE :: WORK
! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
ALLOCATE(WORK(MAX(1,N*64)))
!ALLOCATE(WORK(MAX(1,N*64)), STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not ALLOCATE correctly in SGETRI95. The error message is: ', errmesg
CALL sgetri( N, A, LDA, IPIV, WORK, LWORK, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(WORK)
!DEALLOCATE(WORK, STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not DEALLOCATE correctly in SGETRI95. The error message is: ', errmesg
END SUBROUTINE sgetri95

SUBROUTINE dgetri95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: N, LDA, LWORK
! N = Order of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize returned
!         by ILAENV (64)
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
REAL(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: A
! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
! On exit, if info = 0, the inverse of the original matrix A
INTEGER, DIMENSION(:), INTENT(IN) :: IPIV
! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: WORK
! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
ALLOCATE(WORK(MAX(1,N*64)))
!ALLOCATE(WORK(MAX(1,N*64)), STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not ALLOCATE correctly in SGETRI95. The error message is: ', errmesg
CALL dgetri( N, A, LDA, IPIV, WORK, LWORK, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(WORK)
!DEALLOCATE(WORK, STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not DEALLOCATE correctly in SGETRI95. The error message is: ', errmesg
END SUBROUTINE dgetri95

SUBROUTINE cgetri95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: N, LDA, LWORK
! N = Order of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize returned
!         by ILAENV (64)
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: A
! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
! On exit, if info = 0, the inverse of the original matrix A
INTEGER, DIMENSION(:), INTENT(IN) :: IPIV
! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
COMPLEX(kind=SGL), DIMENSION(:), ALLOCATABLE :: WORK
! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
ALLOCATE(WORK(MAX(1,N*64)))
!ALLOCATE(WORK(MAX(1,N*64)), STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not ALLOCATE correctly in SGETRI95. The error message is: ', errmesg
CALL cgetri( N, A, LDA, IPIV, WORK, LWORK, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(WORK)
!DEALLOCATE(WORK, STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not DEALLOCATE correctly in SGETRI95. The error message is: ', errmesg
END SUBROUTINE cgetri95

SUBROUTINE zgetri95( A, IPIV, INFO )
IMPLICIT NONE
INTEGER :: N, LDA, LWORK
! N = Order of matrix A. N >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,M).
! LWORK = Dimension of array WORK. LWORK >= max(1,n). For optimal performance LWORK >= N*NB where NB is optimal blocksize returned
!         by ILAENV (64)
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, U(i,i) is exactly zero. The factorization has been completed but the factor U is exactly singular so the
!           solution could not be computed
INTEGER :: INFO2
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: A
! On entry, the factors L and U from the factorization A = P*L*U as computed by SGETRF
! On exit, if info = 0, the inverse of the original matrix A
INTEGER, DIMENSION(:), INTENT(IN) :: IPIV
! Indices from SGETRF that define the permutation matrix P; row i of matrix was interchanged with row IPIV(i)
COMPLEX(kind=DBL), DIMENSION(:), ALLOCATABLE :: WORK
! Workspace. On exit, if INFO = 0, WORK(1) returns the optimal LWORK
CHARACTER(len=120) :: errmesg
LDA = SIZE(A,1)
N   = SIZE(A,2)
ALLOCATE(WORK(MAX(1,N*64)))
!ALLOCATE(WORK(MAX(1,N*64)), STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not ALLOCATE correctly in SGETRI95. The error message is: ', errmesg
CALL zgetri( N, A, LDA, IPIV, WORK, LWORK, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
DEALLOCATE(WORK)
!DEALLOCATE(WORK, STAT=info2, ERRMSG=errmesg)
!IF(info2 /= 0) WRITE(*,*) 'WORK did not DEALLOCATE correctly in SGETRI95. The error message is: ', errmesg
END SUBROUTINE zgetri95

!---------------------------------------------------------------

SUBROUTINE strtrs95( A, B, UPLO, TRANS, DIAG, INFO )
IMPLICIT NONE
CHARACTER(len=1), INTENT(IN), OPTIONAL :: UPLO, TRANS, DIAG
! UPLO = U: A is upper triangular. This is the default value.
! UPLO = L: A is lower triangular
! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose). This is the default value.
!                                                 T: A**T * X = B (Transpose)
!                                                 C: A**H * X = B (Conjugate transpose = Transpose)
! DIAG = N: A is a non-unit triangular. This is the default value.
!      = U: A is unit triangular
CHARACTER(len=1) :: UPLO2, TRANS2, DIAG2
INTEGER :: N, NRHS, LDA, LDB, INFO2
! N = Order of matrix A. N >= 0.
! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
REAL(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A
! Triangular matrix A.
! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
!              is not referenced.
! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
!              is not referenced.
! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
REAL(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: B
! On entry, right hand side matrix.
! On exit, if INFO=0, solution matrix X

IF(PRESENT(UPLO)) THEN
     UPLO2 = UPLO
ELSE
     UPLO2 = 'U'
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
IF(PRESENT(DIAG)) THEN
     DIAG2 = DIAG
ELSE
     DIAG2 = 'N'
END IF
LDA  = SIZE(A,1)
N    = SIZE(A,2)
LDB  = SIZE(B,1)
NRHS = SIZE(B,2)
CALL strtrs( UPLO2, TRANS2, DIAG2, N, NRHS, A, LDA, B, LDB, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
END SUBROUTINE strtrs95

SUBROUTINE dtrtrs95( A, B, UPLO, TRANS, DIAG, INFO )
IMPLICIT NONE
CHARACTER(len=1), INTENT(IN), OPTIONAL :: UPLO, TRANS, DIAG
! UPLO = U: A is upper triangular. This is the default value.
! UPLO = L: A is lower triangular
! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose). This is the default value.
!                                                 T: A**T * X = B (Transpose)
!                                                 C: A**H * X = B (Conjugate transpose = Transpose)
! DIAG = N: A is a non-unit triangular. This is the default value.
!      = U: A is unit triangular
CHARACTER(len=1) :: UPLO2, TRANS2, DIAG2
INTEGER :: N, NRHS, LDA, LDB, INFO2
! N = Order of matrix A. N >= 0.
! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
REAL(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A
! Triangular matrix A.
! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
!              is not referenced.
! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
!              is not referenced.
! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
REAL(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: B
! On entry, right hand side matrix.
! On exit, if INFO=0, solution matrix X

IF(PRESENT(UPLO)) THEN
     UPLO2 = UPLO
ELSE
     UPLO2 = 'U'
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
IF(PRESENT(DIAG)) THEN
     DIAG2 = DIAG
ELSE
     DIAG2 = 'N'
END IF
LDA  = SIZE(A,1)
N    = SIZE(A,2)
LDB  = SIZE(B,1)
NRHS = SIZE(B,2)
CALL dtrtrs( UPLO2, TRANS2, DIAG2, N, NRHS, A, LDA, B, LDB, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
END SUBROUTINE dtrtrs95

SUBROUTINE ctrtrs95( A, B, UPLO, TRANS, DIAG, INFO )
IMPLICIT NONE
CHARACTER(len=1), INTENT(IN), OPTIONAL :: UPLO, TRANS, DIAG
! UPLO = U: A is upper triangular. This is the default value.
! UPLO = L: A is lower triangular
! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose). This is the default value.
!                                                 T: A**T * X = B (Transpose)
!                                                 C: A**H * X = B (Conjugate transpose = Transpose)
! DIAG = N: A is a non-unit triangular. This is the default value.
!      = U: A is unit triangular
CHARACTER(len=1) :: UPLO2, TRANS2, DIAG2
INTEGER :: N, NRHS, LDA, LDB, INFO2
! N = Order of matrix A. N >= 0.
! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(IN) :: A
! Triangular matrix A.
! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
!              is not referenced.
! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
!              is not referenced.
! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
COMPLEX(kind=SGL), DIMENSION(:, :), INTENT(INOUT) :: B
! On entry, right hand side matrix.
! On exit, if INFO=0, solution matrix X

IF(PRESENT(UPLO)) THEN
     UPLO2 = UPLO
ELSE
     UPLO2 = 'U'
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
IF(PRESENT(DIAG)) THEN
     DIAG2 = DIAG
ELSE
     DIAG2 = 'N'
END IF
LDA  = SIZE(A,1)
N    = SIZE(A,2)
LDB  = SIZE(B,1)
NRHS = SIZE(B,2)
CALL ctrtrs( UPLO2, TRANS2, DIAG2, N, NRHS, A, LDA, B, LDB, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
END SUBROUTINE ctrtrs95

SUBROUTINE ztrtrs95( A, B, UPLO, TRANS, DIAG, INFO )
IMPLICIT NONE
CHARACTER(len=1), INTENT(IN), OPTIONAL :: UPLO, TRANS, DIAG
! UPLO = U: A is upper triangular. This is the default value.
! UPLO = L: A is lower triangular
! TRANS = Specifies form of system of equations = N: A    * X = B (No transpose). This is the default value.
!                                                 T: A**T * X = B (Transpose)
!                                                 C: A**H * X = B (Conjugate transpose = Transpose)
! DIAG = N: A is a non-unit triangular. This is the default value.
!      = U: A is unit triangular
CHARACTER(len=1) :: UPLO2, TRANS2, DIAG2
INTEGER :: N, NRHS, LDA, LDB, INFO2
! N = Order of matrix A. N >= 0.
! NRHS = # right hand sides (columns of matrix B). NRHS >= 0.
! LDA = Leading dimension of array A. LDA >= max(1,N).
! LDB = Leading dimension of array B. LDB >= max(1,N).
INTEGER, INTENT(OUT), OPTIONAL :: INFO
! INFO = 0: successful exit
!      < 0: if INFO = -i, i-th argument had an illegal value
!      > 0: if INFO = i, A(i,i) is zero. A is singular and solutions X have not been computed
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(IN) :: A
! Triangular matrix A.
! If UPLO = U, leading NxN upper triangular part of A contains upper triangular matrix and strictly lower triangular part of A
!              is not referenced.
! If UPLO = L, leading NxN lower triangular part of A contains lower triangular matrix and strictly upper triangular part of A
!              is not referenced.
! If DIAG = U, diagonal elements of A are also not referenced and are assumed to be 1.
COMPLEX(kind=DBL), DIMENSION(:, :), INTENT(INOUT) :: B
! On entry, right hand side matrix.
! On exit, if INFO=0, solution matrix X

IF(PRESENT(UPLO)) THEN
     UPLO2 = UPLO
ELSE
     UPLO2 = 'U'
END IF
IF(PRESENT(TRANS)) THEN
     TRANS2 = TRANS
ELSE
     TRANS2 = 'N'
END IF
IF(PRESENT(DIAG)) THEN
     DIAG2 = DIAG
ELSE
     DIAG2 = 'N'
END IF
LDA  = SIZE(A,1)
N    = SIZE(A,2)
LDB  = SIZE(B,1)
NRHS = SIZE(B,2)
CALL ztrtrs( UPLO2, TRANS2, DIAG2, N, NRHS, A, LDA, B, LDB, INFO2 )
IF(PRESENT(INFO)) INFO = INFO2
END SUBROUTINE ztrtrs95

!---------------------------------------------------------------

END MODULE INTERFACE_DEFINITIONS

