SUBROUTINE mumps_solver(mumps_par, ROOT, MPI_ID, MPI_P, ndof, ndofst, nnte, ki, kj, kv, job)
USE kinds
USE interface_definitions
!$ use omp_lib
!-------------------------------------------------------------------------------
! Declare variables
!-------------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'

INTEGER, INTENT(IN) :: nnte, ndof, ndofst, job
INTEGER, INTENT(IN) :: MPI_P, MPI_ID, ROOT
INTEGER, INTENT(IN) :: ki(*), kj(*)
REAL(KIND=8), INTENT(IN) :: kv(*)
TYPE (DMUMPS_STRUC), INTENT(INOUT) :: mumps_par
! Counters
INTEGER :: ii, jj, kk, MPI_IERR
!-------------------------------------------------------------------------------
! MUMPS: initialize, JOB = -1
!-------------------------------------------------------------------------------
IF (job .EQ. 0) THEN
   mumps_par%COMM = MPI_COMM_WORLD
   mumps_par%JOB = -1
   mumps_par%SYM = 2 ! 0 unsymmetric, 1 SPD, 2 general symmetric
   mumps_par%PAR = 1 ! 0 host not work, 1 host work
   CALL DMUMPS(mumps_par)
   IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
           & "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
           & "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
           !GOTO 1500
   ENDIF
   IF (MPI_ID .EQ. ROOT) THEN
       OPEN(UNIT=120, FILE='MUMPS_errmsg.log', STATUS='unknown')
       OPEN(UNIT=125, FILE='MUMPS_digmsg.log', STATUS='unknown')
       OPEN(UNIT=130, FILE='MUMPS_glbmsg.log', STATUS='unknown')
   ENDIF
   ! output stream for error messages
   mumps_par%ICNTL(1) = 120 ! <=0 suppressed
   ! output stream for diagnostic printing, statistics, warning messages
   mumps_par%ICNTL(2) = 0 ! <=0 suppressed
   ! output stream for global information
   mumps_par%ICNTL(3) = 130 ! <=0 suppressed
   ! level of printing for error, warning and diagnostic messages
   mumps_par%ICNTL(4) = 2
   ! using METIS for ordering
   mumps_par%ICNTL(7) = 5
   ! percentage increase in the estimated working space
   mumps_par%ICNTL(14) = 60
   mumps_par%N = ndof
ENDIF
kk = 0
DO ii = 1, ndof
    DO jj = ki(ii), ki(ii+1)-1
        IF (ii .GE. kj(jj)) THEN
            kk = kk + 1
        ENDIF
    ENDDO
ENDDO
IF (MPI_P .EQ. 1) THEN
!-------------------------------------------------------------------------------
! MUMPS: define problem on host thread
!-------------------------------------------------------------------------------
     mumps_par%NZ = kk
     IF (job .EQ. 0) Then
         ALLOCATE(mumps_par%RHS(ndofst))
     ELSE
         DEALLOCATE(mumps_par%IRN, mumps_par%JCN, mumps_par%A)
     ENDIF
     ALLOCATE(mumps_par%IRN(mumps_par%NZ))
     ALLOCATE(mumps_par%JCN(mumps_par%NZ))
     ALLOCATE(mumps_par%A(mumps_par%NZ))
     kk = 1
     DO ii = 1, ndof
         DO jj = ki(ii), ki(ii+1)-1
             IF (ii .GE. kj(jj)) THEN
                 mumps_par%IRN(kk) = ii
                 mumps_par%JCN(kk) = kj(jj)
                 mumps_par%A(kk) = kv(jj)
                 kk = kk + 1
             ENDIF
         ENDDO
     ENDDO
     IF (job .EQ. 0) Then
         mumps_par%NRHS = nnte         ! Multiple RHS
         mumps_par%LRHS = mumps_par%N    ! >= N
     ENDIF
ELSE
!-------------------------------------------------------------------------------
! MUMPS: define problem on all threads
!-------------------------------------------------------------------------------
    ! distribution of input matrix
    IF (job .EQ. 0) Then
        mumps_par%ICNTL(18) = 3
    ELSE
        DEALLOCATE(mumps_par%IRN_loc, mumps_par%JCN_loc, mumps_par%A_loc)
    ENDIF
    mumps_par%NZ_loc = kk
    ALLOCATE(mumps_par%IRN_loc(mumps_par%NZ_loc))
    ALLOCATE(mumps_par%JCN_loc(mumps_par%NZ_loc))
    ALLOCATE(mumps_par%A_loc(mumps_par%NZ_loc))
    kk = 1
    DO ii = 1, ndof
        DO jj = ki(ii), ki(ii+1)-1
            IF (ii .GE. kj(jj)) THEN
                mumps_par%IRN_loc(kk) = ii
                mumps_par%JCN_loc(kk) = kj(jj)
                mumps_par%A_loc(kk) = kv(jj)
                kk = kk + 1
            ENDIF
        ENDDO
    ENDDO
    IF (MPI_ID .EQ. ROOT .AND. job .EQ. 0) THEN
        ALLOCATE(mumps_par%RHS(ndofst))
        mumps_par%NRHS = nnte         ! Multiple RHS
        mumps_par%LRHS = mumps_par%N    ! >= N
    ENDIF
ENDIF
!-------------------------------------------------------------------------------
! MUMPS: analysis, JOB = 1
!-------------------------------------------------------------------------------
mumps_par%JOB = 1
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
    WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
    !GOTO 1500
ENDIF
!-------------------------------------------------------------------------------
! MUMPS: factorize, JOB = 2
!-------------------------------------------------------------------------------
mumps_par%JOB = 2
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
    WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
    !GOTO 1500
ENDIF
RETURN
ENDSUBROUTINE mumps_solver
