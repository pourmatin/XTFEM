! Timestamp: Thu Mar 10 14:20:29 CST 2016 
!_______________________________________________________________________________
! 
! Subourines for Kronecker product in space-time finite element method
!   extended with MUMPS (MPI-based direct solver)
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Mar., 2016
! 
! ---Subroutines---
!    Name               Description
! ------------------  ------------------------------------------------------
!  mpipgmres           Solve (A*K + B*M) x = y by PGMRES
!  mpilusol            Solve (A*LU) y = x
! 
!  Note that * is Kronecker product

SUBROUTINE mpikronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, x, y)
!-------------------------------------------------------------------------------
! Compute y = (A*K + B*M)x, note that * is Kronecker product 
!
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2016
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
!
!! Variable declaration
IMPLICIT NONE
INCLUDE 'mpif.h'
! ---External variables---
INTEGER, INTENT(IN) :: p, q
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*)
REAL(KIND=8), INTENT(IN) :: kv(*), mv(*)
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a, b
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
INTEGER :: i, j, ROOT, MPI_IERR
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: z
ALLOCATE(z(p*q))
z = 0.d0
ROOT = 0
! Perform kronspmv on all threads
CALL kronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, x, z)
! Collect results on host thread
DO i = 1, q
    CALL MPI_REDUCE(z((i-1)*p+1:i*p), y((i-1)*p+1:i*p), p, MPI_REAL8, &
        & MPI_SUM, ROOT, MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'i = ', i
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        STOP
    ENDIF
ENDDO
!CALL MPI_REDUCE(z, y, p*q, MPI_REAL8, MPI_SUM, ROOT, MPI_COMM_WORLD, MPI_IERR)
!IF (MPI_IERR .NE. 0) THEN
!    PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
!    STOP
!ENDIF
DEALLOCATE(z)
RETURN
END SUBROUTINE mpikronspmv


SUBROUTINE mpikronspgemv(a, b, alpha, beta, q, ki, kj, kv, mi, mj, mv, p, x, y)
!-------------------------------------------------------------------------------
! Compute y := alpha*(A*K + B*M)*x + beta*y, note that * is Kronecker product 
!
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Apr., 2016
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
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: p, q
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*)
REAL(KIND=REKIND), INTENT(IN) :: kv(*), mv(*), alpha, beta
REAL(KIND=REKIND), DIMENSION(q,q), INTENT(IN) :: a, b
REAL(KIND=REKIND), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=REKIND), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: z
ALLOCATE(z(p*q))
z = 0.d0
CALL mpikronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, x, z)
y = alpha*z + beta*y
DEALLOCATE(z)
RETURN
END SUBROUTINE mpikronspgemv

SUBROUTINE mpilusol(a, p, q, mumps_par, x, y)
! Solve (A*LU)y = x for y, note that * is Kronecker product 
! Developed by Rui Zhang (rui.zhang4@utdallas.edu), Dec., 2015
! 
!  ---Entries---
!    Name     Type      Size        Description
!  -------- --------- ----------- --------------------------------
!     a      REAL      (q,q)        A matrix
!     p      INTEGER     1          Dimension of K and M
!     q      INTEGER     1          Dimension of A
!  mumps_par STRUC       1          MUMPS
!     x      REAL       p*q         Input vector
!     y      REAL       p*q         Output vector
!
!! Variable declaration
IMPLICIT NONE
INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
! ---External variables--- 
INTEGER, INTENT(IN) :: p, q
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a
TYPE (DMUMPS_STRUC), INTENT(INOUT) :: mumps_par
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: x
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: y
! ---Internal variables---
INTEGER :: i, j, MPI_IERR, ROOT
REAL(KIND=8), DIMENSION(q,q) :: inva
ROOT = 0
IF (mumps_par%MYID .EQ. 0) THEN
    ! Compute inverse of A matrix
    CALL inverse(q, a, inva)
    ! Perform LUSOL first

    DO i = 1, p*q
        mumps_par%RHS(i) = x(i)
    ENDDO
ENDIF
mumps_par%JOB = 3
!CALL MPI_BARRIER(mumps_par%COMM, MPI_IERR)
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
    WRITE(*,*) 'ERROR: MUMPS, JOB = 3'
    WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
    STOP
ENDIF
IF (mumps_par%MYID .EQ. 0) THEN
    ! Perform Vector addition
    y = 0.d0
    DO i = 1, q
        DO j = 1, q
            y(1+(i-1)*p:i*p) = y(1+(i-1)*p:i*p) + & 
                & inva(i,j)*mumps_par%RHS(1+(j-1)*p:j*p)
        ENDDO
    ENDDO
ENDIF
CALL MPI_BCAST(y, p*q, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    STOP
ENDIF
RETURN
END SUBROUTINE mpilusol


SUBROUTINE mpipgmres(a, b, q, ki, kj, kv, mi, mj, mv, p, mumps_par, rhs, &
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
!  mumps_par STRUC       1          MUMPS
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
INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
! ---External variables---
INTEGER, INTENT(IN) :: p, q
REAL(KIND=8), DIMENSION(q,q), INTENT(IN) :: a, b
INTEGER, INTENT(IN) :: ki(*), mi(*), kj(*), mj(*)
REAL(KIND=8), INTENT(IN) :: kv(*), mv(*)
TYPE (DMUMPS_STRUC), INTENT(INOUT) :: mumps_par
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: rhs
REAL(KIND=8), DIMENSION(p*q), INTENT(INOUT) :: sol
INTEGER, INTENT(IN) :: im, maxits, iout
INTEGER, INTENT(OUT) :: ierr
REAL(KIND=8), INTENT(IN) :: eps
! ---BLAS functions---
REAL(KIND=8) :: dnrm2, ddot
! ---Internal variables---
INTEGER :: i, i1, ii, j, jj, k, k1, n, n1, its, MPI_IERR, stat
REAL(KIND=8) :: epsmac, ro0, ro, t, eps1, gam
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: c, s, rs, v0
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hh, vv
199 FORMAT('   its =', i4, ' rel. res. norm =', d20.6)
n = p*q
ALLOCATE(v0(n))
v0 = 0.d0
IF (mumps_par%MYID .EQ. 0) THEN 
    ALLOCATE(hh(im+1,im), c(im), s(im), rs(im+1), vv(n,im+1))
    vv = 0.d0
    epsmac = 1.d-16
    ro0 = dnrm2(n, rhs, 1)
ENDIF
!-------------- compute initial residual vector --------------
!DO i = 1, 6
!    DO j= 1, 6
!        PRINT *, 'a(',i, ',',j,')=', a(i,j), b(i,j)
!    ENDDO
!ENDDO
CALL mpikronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, sol, v0)
IF (mumps_par%MYID .EQ. 0) THEN
    vv(1:n,1) = v0
    DO j = 1, n
        vv(j,1) = rhs(j) - vv(j,1)
    ENDDO
ENDIF
! Outer loop
stat = 0
its = 0
DO WHILE (its .LT. maxits)
    IF (mumps_par%MYID .EQ. 0) THEN
        ro = dnrm2(n, vv(1:n,1), 1)
        IF ((iout .GE. 0) .AND. (its .EQ. 0)) THEN
            WRITE(iout, 199) its, ro/ro0
        ENDIF
        IF (ro .EQ. 0.0d0) THEN
            ierr = -1
            stat = 1
        ENDIF
    ENDIF
    CALL MPI_BCAST(stat, 1, MPI_INTEGER, 0, mumps_par%COMM, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_BCAST error = ',  MPI_IERR
        STOP
    ENDIF
    IF (stat .NE. 0) GOTO 1000
    IF (mumps_par%MYID .EQ. 0) THEN
        t = 1.0d0/ro
        DO j = 1, n
            vv(j,1) = vv(j,1)*t
        ENDDO
        IF (its .EQ. 0) eps1 = eps*ro
        ! ** initialize 1-st term  of rhs of hessenberg system
        rs(1) = ro
    ENDIF
    ! Inner loop
    DO i = 1, im
        its = its + 1
        i1 = i + 1
        IF (mumps_par%MYID .EQ. 0) THEN
            v0 = vv(1:n,i)
        ENDIF
        CALL mpilusol(a, p, q, mumps_par, v0, rhs)
        CALL mpikronspmv(a, b, q, ki, kj, kv, mi, mj, mv, p, rhs, v0)
        IF (mumps_par%MYID .EQ. 0) THEN
            vv(1:n,i1) = v0
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
            IF (t .NE. 0.0d0) THEN
                t = 1.0d0/t
                DO  k = 1, n
                    vv(k,i1) = vv(k,i1)*t
                ENDDO
            ENDIF
            ! done with modified gram schimd and arnoldi step
            ! now  update factorization of hh
            IF (i .NE. 1) THEN
                !--------perfrom previous transformations  on i-th column of h
                DO k = 2, i
                    k1 = k-1
                    t = hh(k1,i)
                    hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
                    hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
                ENDDO
            ENDIF
            gam = DSQRT(hh(i,i)**2 + hh(i1,i)**2)
            ! if gamma is zero then any small value will do
            ! will affect only residual estimate
            IF (gam .EQ. 0.0d0) gam = epsmac
            ! get next plane rotation
            c(i) = hh(i,i)/gam
            s(i) = hh(i1,i)/gam
            rs(i1) = -s(i)*rs(i)
            rs(i) =  c(i)*rs(i)
            ! detrermine residual norm and test for convergence-
            hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
            ro = DABS(rs(i1))
!PRINT *, 'itr:', its, 'ro= ', ro
            IF (iout .GT. 0) WRITE(iout, 199) its, ro/ro0
            IF (ro .LE. eps1)  THEN
                stat = 1
            ENDIF
        ENDIF
        CALL MPI_BCAST(stat, 1, MPI_INTEGER, 0, mumps_par%COMM, MPI_IERR)
        IF (MPI_IERR .NE. 0) THEN
            PRINT *, 'MPI_BCAST error = ',  MPI_IERR
            STOP
        ENDIF
        IF (stat .NE. 0)  THEN
            ! now compute solution. 
            IF (mumps_par%MYID .EQ. 0) THEN
                ! first solve upper triangular system.
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
            ENDIF
            ! call preconditioner.
            CALL mpilusol(a, p, q, mumps_par, rhs, rhs)
            IF (mumps_par%MYID .EQ. 0) THEN
                DO k = 1, n
                    sol(k) = sol(k) + rhs(k)
                ENDDO
                ierr = 0
            ENDIF
            ! exit
            GOTO 1000
        ENDIF
    ENDDO ! inner loop
    IF (mumps_par%MYID .EQ. 0) THEN
        DO j = 1, i-1
            jj = i1 - j + 1
            rs(jj-1) = -s(jj-1)*rs(jj)
            rs(jj) = c(jj-1)*rs(jj)
        ENDDO
        DO j = 1, i1-1
            t = rs(j)
            IF (j .EQ. 1)  t = t-1.0d0
            CALL daxpy (n, t, vv(1,j), 1,  vv, 1)
        ENDDO
    ENDIF
ENDDO ! outer loop
ierr = 1
1000 CONTINUE
CALL MPI_BCAST(ierr, 1, MPI_INTEGER, 0, mumps_par%COMM, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    STOP
ENDIF
IF (mumps_par%MYID .EQ. 0) DEALLOCATE(hh, c, s, rs, vv)
DEALLOCATE(v0)
RETURN
END SUBROUTINE mpipgmres
