! Timestamp: Fri Dec 11 16:13:35 CST 2015 
!_______________________________________________________________________________
! 
! Subourines for Kronecker product in space-time finite element method
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
! ---Subroutines---
!    Name               Description
! ------------------  ------------------------------------------------------
!  kronspmv             Compute y := (A*K + B*M) x
!  kronspgemv           Compute y := alpha (A*K + B*M) x + beta y
!  kronlusol            Solve (A*LU) y = x
!  kronpgmres           Solve (A*K + B*M) x = y by PGMRES
!  inverse              Inverse a general square matrix by LAPACK
! 
!  Note that * is Kronecker product


SUBROUTINE kronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, x, y)
!-------------------------------------------------------------------------------
! Compute y = (A*K + B*M)x, note that * is Kronecker product 
!
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
!  ---Entries---
!    Name     Type      Size        Description
!  -------- --------- ----------- --------------------------------
!     a      REAL      (q,q)        A matrix
!     b      REAL      (q,q)        B matrix
!     q      INTEGER     1          Dimension of A and B
!     ki     INTEGER    p+1         Row pointer of K in CSR format
!     kj     INTEGER     *          Colume indices of K in CSR format
!     kv     REAL        *          Nonzero values of K in CSR format
!     mi     INTEGER    p+1         Row pointer of M in CSR format
!     mj     INTEGER     *          Colume indices of M in CSR format
!     mv     REAL        *          Nonzero values of M in CSR format
!     p      INTEGER     1          Dimension of K and M
!     x      REAL       p*q         Input vector
!     y      REAL       p*q         Output vector

!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: p, q
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*)
REAL(KIND=8), INTENT(IN) :: kv(*), mv(*)
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a, b
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
INTEGER :: i, j
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u
ALLOCATE(u(p*q))
y = 0.d0
! K part
u = 0.d0
! Perform spMV first
!DO i = 1, 6
!    DO j= 1, 6
!        PRINT *, 'a(',i, ',',j,')=', a(i,j), b(i,j)
!    ENDDO
!ENDDO
DO i = 1, q
    CALL amux(p, x(1+(i-1)*p:i*p), u(1+(i-1)*p:i*p), kv, kj, ki)
ENDDO
! Perform Vector addition
DO i = 1, q
    DO j = 1, q
        y(1+(i-1)*p:i*p) = y(1+(i-1)*p:i*p) + a(i,j)*u(1+(j-1)*p:j*p)
    ENDDO
ENDDO
u = 0.d0
! Perform spMV first
DO i = 1, q
    CALL amux(p, x(1+(i-1)*p:i*p), u(1+(i-1)*p:i*p), mv, mj, mi) 
ENDDO
! Perform Vector addition
DO i = 1, q
    DO j = 1, q
        y(1+(i-1)*p:i*p) = y(1+(i-1)*p:i*p) + b(i,j)*u(1+(j-1)*p:j*p)
    ENDDO
ENDDO
DEALLOCATE(u)
RETURN
END SUBROUTINE kronspmv


SUBROUTINE kronspgemv(a, b, alpha, beta, q, ki, kj, kv, mi, mj, mv, p, x, y)
!-------------------------------------------------------------------------------
! Compute y := alpha*(A*K + B*M)*x + beta*y, note that * is Kronecker product 
!
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
!  ---Entries---
!    Name     Type      Size        Description
!  -------- --------- ----------- --------------------------------
!     a      REAL      (q,q)        A matrix
!     b      REAL      (q,q)        B matrix
!    alpha   REAL        1        
!    beta    REAL        1
!     q      INTEGER     1          Dimension of A and B
!     ki     INTEGER    p+1         Row pointer of K in CSR format
!     kj     INTEGER     *          Colume indices of K in CSR format
!     kv     REAL        *          Nonzero values of K in CSR format
!     mi     INTEGER    p+1         Row pointer of M in CSR format
!     mj     INTEGER     *          Colume indices of M in CSR format
!     mv     REAL        *          Nonzero values of M in CSR format
!     p      INTEGER     1          Dimension of K and M
!     x      REAL       p*q         Input vector
!     y      REAL       p*q         Output vector

!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: p, q
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*)
REAL(KIND=8), INTENT(IN) :: kv(*), mv(*), alpha, beta
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a, b
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
REAL(KIND=8), DIMENSION(p*q) :: z
z = 0.d0
CALL kronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, x, z)
y = alpha*z + beta*y
RETURN
END SUBROUTINE kronspgemv


SUBROUTINE kronlusol(a, q, alu, jlu, ju, p, x, y)
! Solve (A*LU)y = x for y, note that * is Kronecker product 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
!  ---Entries---
!    Name     Type      Size        Description
!  -------- --------- ----------- --------------------------------
!     a      REAL      (q,q)        A matrix
!     q      INTEGER     1          Dimension of A
!     alu    REAL        *          MSR format
!     jlu    INTEGER     *                      of
!     ju     INTEGER     *                          L and U
!     p      INTEGER     1          Dimension of L and U
!     x      REAL       p*q         Input vector
!     y      REAL       p*q         Output vector
! 
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: p, q
REAL(KIND=8), INTENT(IN) :: alu(*)
INTEGER, INTENT(IN) :: jlu(*), ju(*)
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
INTEGER :: i, j
REAL(KIND=8), DIMENSION(q,q) :: inva
REAL(KIND=8), DIMENSION(p*q) :: w
! Compute inverse of A matrix
CALL inverse(q, a, inva)
! Perform LUSOL first
w = 0.d0
!!$omp parallel num_threads(6)
!!$omp do
DO i = 1, q
    CALL lusol(p, x(1+(i-1)*p:i*p), w(1+(i-1)*p:i*p), alu, jlu, ju)
ENDDO
!!$omp end do
!!$omp end parallel
! Perform Vector addition
y = 0.d0
DO i = 1, q
    DO j = 1, q
        y(1+(i-1)*p:i*p) = y(1+(i-1)*p:i*p) + inva(i,j)*w(1+(j-1)*p:j*p)
    ENDDO
ENDDO
RETURN
END SUBROUTINE kronlusol


SUBROUTINE kronpgmres(a, b, q, ki, kj, kv, mi, mj, mv, p, alu, jlu, ju, rhs, &
    &  sol, im, eps, maxits, iout, ierr)
!-------------------------------------------------------------------------------
! Solve (A*K + B*M) x = y, note that * is Kronecker product 
!
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
!  ---Entries---
!    Name     Type      Size        Description
!  -------- --------- ----------- --------------------------------
!     a      REAL      (q,q)        A matrix
!     b      REAL      (q,q)        B matrix
!     q      INTEGER     1          Dimension of A and B
!     ki     INTEGER    p+1         Row pointer of K in CSR format
!     kj     INTEGER     *          Colume indices of K in CSR format
!     kv     REAL        *          Nonzero values of K in CSR format
!     mi     INTEGER    p+1         Row pointer of M in CSR format
!     mj     INTEGER     *          Colume indices of M in CSR format
!     mv     REAL        *          Nonzero values of M in CSR format
!     p      INTEGER     1          Order of K and M
!     alu    REAL        *          MSR format
!     jlu    INTEGER     *                      of
!     ju     INTEGER     *                          L and U
!     rhs    REAL       p*q         Input vector
!     sol    REAL       p*q         Output vector
!     im     INTEGER     1          Size of Krylov subspace
!     eps    REAL        1          Tolerance for stopping criterion
!  maxits    INTEGER     1          Maximum number of iterations allowed 
!    iout    INTEGER     1          Output file unit number
!    ierr    INTEGER     1          Error message with the following meaning
!                           0 --> successful return
!                           1 --> convergence not achieved in itmax iterations
!                          -1 --> the initial guess seems to be the exact
!                                 solution (initial residual computed was zero)
!
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: p, q
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a, b
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*), jlu(*), ju(*)
REAL(KIND=8), INTENT(IN) :: kv(*), mv(*), alu(*)
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: rhs
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: sol
INTEGER, INTENT(IN) :: im, maxits, iout
INTEGER, INTENT(OUT) :: ierr
REAL(KIND=8), INTENT(IN) :: eps

! ---BLAS functions---
REAL(KIND=8) :: dnrm2, ddot
! ---Internal variables---
INTEGER :: i, i1, ii, j, jj, k, k1, n, n1, its
REAL(KIND=8) :: epsmac, ro0, ro, t, eps1, gam
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: c, s, rs
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hh, vv
!
n = p*q
ALLOCATE(hh(im+1,im), c(im), s(im), rs(im+1), vv(n,im+1))
vv = 0.d0
epsmac = 1.d-16
!-------------------------------------------------------------
! Arnoldi size should not exceed kmax=500 in this version.
! to reset modify paramter kmax accordingly.
!-------------------------------------------------------------
n1 = n + 1
its = 0
ro0 = dnrm2(n, rhs, 1)
!-------------------------------------------------------------
! outer loop starts here
!-------------- compute initial residual vector --------------
CALL kronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, sol, vv(1:n,1))
DO j = 1, n
    vv(j,1) = rhs(j) - vv(j,1)
ENDDO
!-------------------------------------------------------------
20  ro = dnrm2(n, vv(1:n,1), 1)
IF ((iout .GE. 0) .AND. (its .EQ. 0)) THEN
    WRITE(iout, 199) its, ro/ro0
ENDIF
IF (ro .EQ. 0.0d0) GOTO 999
t = 1.0d0/ro
DO j = 1, n
    vv(j,1) = vv(j,1)*t
ENDDO
IF (its .EQ. 0) eps1 = eps*ro
! ** initialize 1-st term  of rhs of hessenberg system..
rs(1) = ro
i = 0
4   i = i + 1
its = its + 1
i1 = i + 1
CALL kronlusol(a, q, alu, jlu, ju, p, vv(1:n,i), rhs)
CALL kronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, rhs, vv(1:n,i1))
!-----------------------------------------
!     modified Gram - Schmidt
!-----------------------------------------
DO j = 1, i
    t = ddot(n, vv(1,j), 1, vv(1,i1), 1)
    hh(j,i) = t
    CALL daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
ENDDO
t = dnrm2(n, vv(1:n,i1), 1)
hh(i1,i) = t
IF (t .EQ. 0.0d0) GOTO 58
t = 1.0d0/t
DO  k = 1, n
    vv(k,i1) = vv(k,i1)*t
ENDDO
! done with modified gram schimd and arnoldi step..
! now  update factorization of hh
58  IF (i .EQ. 1) GOTO 121
!--------perfrom previous transformations  on i-th column of h
DO k = 2, i
    k1 = k-1
    t = hh(k1,i)
    hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
    hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
ENDDO
121 gam = DSQRT(hh(i,i)**2 + hh(i1,i)**2)
! if gamma is zero then any small value will do...
! will affect only residual estimate
IF (gam .EQ. 0.0d0) gam = epsmac
! get  next plane rotation
 c(i) = hh(i,i)/gam
s(i) = hh(i1,i)/gam
rs(i1) = -s(i)*rs(i)
rs(i) =  c(i)*rs(i)
! detrermine residual norm and test for convergence-
hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
ro = DABS(rs(i1))
IF (iout .GT. 0) WRITE(iout, 199) its, ro/ro0
IF ((i .LT. im) .AND. (ro .GT. eps1))  GOTO 4
! now compute solution. first solve upper triangular system.
rs(i) = rs(i)/hh(i,i)
DO ii = 2, i
    k = i - ii + 1
    k1 = k + 1
    t = rs(k)
    DO j = k1, i
        t = t - hh(k,j)*rs(j)
    ENDDO
    rs(k) = t/hh(k,k)
ENDDO
! form linear combination of v(*,i)'s to get solution
t = rs(1)
DO k = 1, n
    rhs(k) = vv(k,1)*t
ENDDO
DO j = 2, i
    t = rs(j)
    DO k = 1, n
        rhs(k) = rhs(k) + t*vv(k,j)
    ENDDO
ENDDO
! call preconditioner.
CALL kronlusol(a, q, alu, jlu, ju, p, rhs, rhs)
DO k = 1, n
    sol(k) = sol(k) + rhs(k)
ENDDO
! restart outer loop  when necessary
IF (ro .LE. eps1) GOTO 990
IF (its .GE. maxits) GOTO 991
! else compute residual vector and continue..
DO j = 1, i
    jj = i1 - j + 1
    rs(jj-1) = -s(jj-1)*rs(jj)
    rs(jj) = c(jj-1)*rs(jj)
ENDDO
DO j = 1, i1
    t = rs(j)
    IF (j .EQ. 1)  t = t-1.0d0
    CALL daxpy (n, t, vv(1,j), 1,  vv, 1)
ENDDO
199 FORMAT('   its =', i4, ' rel. res. norm =', d20.6)
! restart outer loop.
GOTO 20
990 ierr = 0
GOTO 1000
991 ierr = 1
GOTO 1000
999 ierr = -1
1000 CONTINUE
DEALLOCATE(hh, c, s, rs)
RETURN
END SUBROUTINE kronpgmres


SUBROUTINE inverse(n, A, Ainv)
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
IMPLICIT NONE
REAL(KIND=8), DIMENSION(n,n), INTENT(IN) :: A
REAL(KIND=8), DIMENSION(n,n) :: Ainv
REAL(KIND=8), DIMENSION(n) :: work  ! work array for LAPACK
INTEGER, DIMENSION(n) :: ipiv   ! pivot indices
INTEGER :: n, info
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
CALL DGETRF(n, n, Ainv, n, ipiv, info)
IF (info /= 0) THEN
   STOP 'Matrix is numerically singular!'
ENDIF
! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
CALL DGETRI(n, Ainv, n, ipiv, work, n, info)
IF (info /= 0) THEN
   STOP 'Matrix inversion failed!'
ENDIF
RETURN
END SUBROUTINE inverse
