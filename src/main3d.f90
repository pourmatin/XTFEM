! Timestamp: Wed Oct 14 10:00:56 CDT 2015 
!_______________________________________________________________________________
!
! Space-time finite element method for elastodynamics in three-dimensions (3D) 
! Developed by Rui Zhang @ UT Dallas, Richardson, TX, Sep., 2015
! E-mail: rui.zhang4@utdallas.edu (USA), ruizhang@mail.nwpu.edu.cn (China)
! Supervised by Prof. Dong Qian (UTD, USA) and Prof. Lihua Wen (NPU, China)
!_______________________________________________________________________________
PROGRAM main3d
USE kinds
USE interface_definitions
!$ use omp_lib
!-------------------------------------------------------------------------------
! Declare variables
!-------------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'mpif.h'
INCLUDE 'dmumps_struc.h'
! Parameters
REAL(KIND=REKIND) :: kappa, cp, rho, le, area, freq
INTEGER :: ne, nn, ndofn, nnsd, nnte, nquads, nquadt, ndof, ndofst, nnse, &
    & ngp, ngpe, ncmp, n1, n2
REAL(KIND=REKIND) :: dt, dt0, q0
INTEGER :: nstep, nip, nipc, nout
INTEGER :: errstat
CHARACTER(LEN=80) :: job, inp, dat, temp, errmesg
! Mesh
INTEGER, DIMENSION(:,:), ALLOCATABLE :: xconn
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: xcoord
! Shape functions
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: wqs
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: ns
REAL(KIND=REKIND), DIMENSION(:,:,:), ALLOCATABLE :: dns
! Spatial matrices
INTEGER :: knz, mnz
INTEGER, DIMENSION(:), POINTER :: ki, kj, mi, mj, kit, kjt, mit, mjt, ki0, kj0, mi0, mj0
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: mv, kv, mvt, kvt, mv0, kv0
! Space-time coefficients
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: tcoorde
REAL(KIND=REKIND), DIMENSION(3,3) :: ck, cm, cktn, cmtn, cktn1, cmtn1
! Boundary conditions
REAL(KIND=REKIND) :: ntele, ftotal, fnode, T0
INTEGER, DIMENSION(:), ALLOCATABLE :: fixedbcst, fixedbcst0 
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: rext, nodeshare
! Linear system solver
INTEGER :: im, maxits, iout, iout2, ierr
REAL(KIND=REKIND) :: eps
! Time-loop
REAL(KIND=REKIND) :: angvel
REAL(KIND=REKIND), DIMENSION(3) :: fextc
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: dispst, fextst, fBound, fInit
! Two-scale damage model
!REAL(KIND=REKIND) :: Dc
!REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: p, D, ws, del_po, sig_D!
!REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: strain_tgo, strain_pl, sig_eff, sig_back
! Status flag
!INTEGER :: nef, ndf, nefe, itfail, itfail0
!INTEGER, DIMENSION(:), ALLOCATABLE :: necnt, necnt0, gp_bc, gpsta, elesta, nefail, ndfail
! Counters
INTEGER :: i, j, k, l, m, ii, jj, kk, ll, mm, pp
! Performance statistics
REAL(KIND=REKIND) :: mem(10), tcpu(10, 3), tusr(10, 3)
! MPI
CHARACTER(LEN=80) :: mpi_job
INTEGER :: MPI_IERR, MPI_P, MPI_ID, ROOT
INTEGER :: mpi_ne, mpi_ngp, fid
INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_ele, mpi_gp
! MUMPS
TYPE (DMUMPS_STRUC) mumps_par
!-------------------------------------------------------------------------------
! Initialize MPI
!-------------------------------------------------------------------------------
CALL MPI_Init(MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_Init error = ',  MPI_IERR
    GOTO 1500
ENDIF
CALL MPI_Comm_size(MPI_COMM_WORLD, MPI_P, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_Comm_size error = ',  MPI_IERR
    GOTO 1500
ENDIF
CALL MPI_Comm_rank(MPI_COMM_WORLD, MPI_ID, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_Comm_rank error = ',  MPI_IERR
    GOTO 1500
ENDIF
ROOT = 0 ! root thread
IF (MPI_P .GT. 1) THEN
    IF (MPI_ID .EQ. ROOT) THEN
        PRINT *, ''
        WRITE(*,'(A,I4)'), ' *** MPI enabled, processes:', MPI_P
        PRINT *, ''
    ENDIF 
ENDIF
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(10,1))
    !$ tusr(10,1) = omp_get_wtime()
    PRINT *, '----------------------------------------------------------'
    PRINT *, '           Extended space-time FEM (XTFEM) 3D'
    PRINT *, '----------------------------------------------------------'
ENDIF

!-------------------------------------------------------------------------------
! Initialize KOKKOS
!-------------------------------------------------------------------------------
!CALL kokkos_init()
!-------------------------------------------------------------------------------
! Set parameters
!-------------------------------------------------------------------------------
! Material, isotropic elastic only
kappa   = .1_REKIND                 ! heat conductivity W/mm.K
cp      = 452._REKIND                ! specific heat capacity J/Kg.K
rho     = .00000786_REKIND           ! Mass density Kg/mm^3 
! Damage model
!Dc      = 0.3_REKIND    ! Damage criteria
! Input file
job     = 'sent'
inp     = TRIM(job)//'.inp'
dat     = TRIM(job)//'.dat'
WRITE(temp,'(i8)') MPI_ID
mpi_job = TRIM(job)//'_'//TRIM(ADJUSTL(temp))
! Geometry
le      = 130._REKIND    ! Height of the plate (loading side)
! Mesh
ndofn   = 1             ! Number of spatial degrees of freedom per node
nnsd    = 2              ! Order of spatial shape functions
nnte    = 3               ! Order of temporal shape functions
nquadt  = 2             ! Number of quadrature points temporally
nquads  = 2            ! Number of quadrature points spatially
nnse    = 8              ! Number of nodes per spatial element
ngpe    = 8              ! Number of Gauss points per spatial element
! Time step
dt0     = 0.05_REKIND   ! Step time
nstep   = 500            ! Number of steps
n1      = 10           ! Initial time step (cycle)
n2      = 1            ! Time step after crack initiation (cycle)
nipc    = 256           ! Number of temporal interpolation points per cycle
! Applid heat flux
q0      = 10._REKIND   ! uniformly distributed heat flux W/mm^2
freq    = .02_REKIND    ! Frequency, in Hz
angvel  = 2._REKIND*pi*freq
area    = 10._REKIND    ! Traction area mm^2
T0      = 10._REKIND    ! Boundary temprature
! Linear system solver
im      = 10000            ! GMRES: size of Krylov subspace 
eps     = 1.d-10         ! GMRES: tolerance for stopping criterion
maxits  = 50000            ! GMRES: maximum number of iterations allowed 
iout    = 40            ! GMRES: solution info output file number
! Post-processing
nout    = 1             ! Output frequency

IF (MPI_ID .EQ. ROOT) THEN
!-------------------------------------------------------------------------------
! Read input file
!-------------------------------------------------------------------------------
    PRINT *,''
    PRINT *, '>>> Reading input file'
    OPEN(UNIT=100, FILE=TRIM(ADJUSTL(inp)), STATUS='old', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        WRITE(*,*) 'Error file: ', TRIM(ADJUSTL(inp))
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, MPI_IERR)
    ENDIF
    READ(100, *) nn ! Read the number of nodes 
    ALLOCATE(xcoord(nn, 3))
    DO i = 1, nn
        READ(100, *) ii, xcoord(ii, :) ! Read node coordinates
    ENDDO
    READ(100, *) ne ! Read the number of elements 
    ALLOCATE(xconn(ne,nnse))
    DO i = 1, ne
        READ(100, *) ii, xconn(ii, 1:nnse) ! Read element connectivities
    ENDDO
    CLOSE(100)
ENDIF
! MPI: broadcast elements and nodes to all threads
CALL MPI_BCAST(ne, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, MPI_IERR) 
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    GOTO 1500
ENDIF
CALL MPI_BCAST(nn, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    GOTO 1500
ENDIF
ndof    = ndofn*nn      ! Number of spatial DOFs
ndofst  = nnte*ndof    ! Number of space-time DOFs
ngp     = ngpe*ne       ! Number of spatial Gauss points
IF (MPI_ID .NE. ROOT) THEN
    ALLOCATE(xcoord(nn, 3), xconn(ne,nnse))
ENDIF
CALL MPI_BCAST(xcoord, nn*3, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    GOTO 1500
ENDIF
CALL MPI_BCAST(xconn, ne*nnse, MPI_INTEGER, ROOT, MPI_COMM_WORLD, MPI_IERR)
IF (MPI_IERR .NE. 0) THEN
    PRINT *, 'MPI_BCAST error = ',  MPI_IERR
    GOTO 1500
ENDIF
! MPI: divide elements and distribute to all threads
CALL mpigroupsize(ne, MPI_P, MPI_ID, mpi_ne)
ALLOCATE(mpi_ele(mpi_ne))
CALL mpisetele(ne, MPI_P, MPI_ID, mpi_ne, mpi_ele)
mpi_ngp = ngpe*mpi_ne
ALLOCATE(mpi_gp(mpi_ngp))
CALL mpisetgp(ngpe, mpi_ne, mpi_ele, mpi_ngp, mpi_gp)
IF (MPI_ID .EQ. ROOT) THEN
    WRITE(*,100) '* Number of spatial elements', ne
    WRITE(*,100) '* Number of spatial nodes   ', nn
    WRITE(*,100) '* Number of spatial DOFs    ', ndof
    WRITE(*,100) '* Number of space-time DOFs ', ndofst
    ! Open data file
    OPEN(UNIT=110, FILE=TRIM(ADJUSTL(dat)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        WRITE(*,*) 'Error file: ', TRIM(ADJUSTL(dat))
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, MPI_IERR)
    ENDIF
    ! Number of MPI threads
    WRITE(110,'(i8)') MPI_P
    ! Dimension
    WRITE(110,'(a8)') '3D'
    ! Number of nodes
    WRITE(110,'(i8)') nn
    ! Nodal coordinates
    DO ii = 1, nn
        WRITE(110,'(3(E20.11E3,2x))') (xcoord(ii,jj),jj=1,3)
    ENDDO
    ! Element type
    WRITE(110,'(a8)') 'C3D8'
    ! Number of elements
    WRITE(110,'(i8)') ne
    ! Element connectivity
    DO ii = 1, ne
        WRITE(110,'(1x,8i10)') (xconn(ii,jj),jj=1,nnse)
    ENDDO
    CLOSE(110)
    CALL initfile(job, MPI_P)
ENDIF
temp = TRIM(mpi_job)//'.ele'
fid = 400 + MPI_ID
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), STATUS='unknown', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    WRITE(*,*) 'Error file: ', TRIM(ADJUSTL(temp))
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, MPI_IERR)
ENDIF
WRITE(fid,'(i8)') mpi_ne
DO ii = 1, mpi_ne
    WRITE(fid,'(i8)') mpi_ele(ii)
ENDDO
CLOSE(fid)
temp = TRIM(mpi_job)//'.gp'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), STATUS='unknown', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    WRITE(*,*) 'Error file: ', TRIM(ADJUSTL(temp))
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, MPI_IERR)
ENDIF
WRITE(fid,'(i8)') mpi_ngp
DO ii = 1, mpi_ngp
    WRITE(fid,'(i8)') mpi_gp(ii)
ENDDO
CLOSE(fid)

!-------------------------------------------------------------------------------
! Set Gauss Points and shape functions
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(1,1))
    !$ tusr(1,1) = omp_get_wtime()
    PRINT *,''
    PRINT *, '>>> Assembling spatial matrices'
ENDIF
ALLOCATE(dns(ngpe, nnse, 3), wqs(nquads), ns(ngpe, nnse))
CALL setshapefunc(nnsd, nnse, nquads, ngpe, wqs, ns, dns)
!-------------------------------------------------------------------------------
! Formulate spatial matrices (spatial integration)
!-------------------------------------------------------------------------------
knz = mpi_ne*(nnse*ndofn)**2
mnz = knz
ALLOCATE(kit(MAX(knz,ndof+1)), kjt(knz), kvt(knz), mit(MAX(mnz,ndof+1)), &
    & mjt(mnz), mvt(mnz))
 CALL mpiglbmtx(mpi_ne, mpi_ele, ne, nn, ndofn, nnse, nquads, ngpe, &
    & xconn, xcoord, kappa, cp, rho, wqs, ns, dns, kit(1:knz), kjt, kvt, knz, mit(1:mnz), mjt, mvt, mnz)
 CALL tricsr(kit, kjt, kvt, ndof, knz, MAX(knz,ndof+1))
 CALL tricsr(mit, mjt, mvt, ndof, mnz, MAX(mnz,ndof+1))
ALLOCATE(ki(ndof+1), kj(knz), kv(knz), mi(ndof+1), mj(mnz), mv(mnz))
ALLOCATE(ki0(ndof+1), kj0(knz), kv0(knz), mi0(ndof+1), mj0(mnz), mv0(mnz))
ki = kit(1:(ndof+1))
kj = kjt(1:knz)
kv = kvt(1:knz)
mi = mit(1:(ndof+1))
mj = mjt(1:mnz)
mv = mvt(1:mnz)
DEALLOCATE(ns, kit, kjt, kvt, mit, mjt, mvt)

IF (MPI_ID .EQ. ROOT) THEN
    !temp = 'stiff'
    !CALL outputcsrmat(ki, kj, kv, ndof, knz, temp)
    !temp = 'mass'
    !CALL outputcsrmat(mi, mj, mv, ndof, mnz, temp)
    CALL CPU_TIME(tcpu(1,2))
    !$ tusr(1,2) = omp_get_wtime()
    mem(1) = (knz*16+mnz*16)/32.0**4.0
    tcpu(1,3) = tcpu(1,2) - tcpu(1,1)
    tusr(1,3) = tusr(1,2) - tusr(1,1)
    WRITE(*,100) '* NNZ of spatial K matrix   ', knz
    WRITE(*,100) '* NNZ of spatial M matrix   ', mnz
    WRITE(*,105) '* Spatial assembly  (sec)   ', tcpu(1,3)
    WRITE(*,105) '* Spatial matrices   (MB)   ', mem(1)
ENDIF
!-------------------------------------------------------------------------------
! Set boundary conditions
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(2,1))
    !$ tusr(2,1) = omp_get_wtime()
    !PRINT *,''
    !PRINT *, '>>> Setting boundary conditions'
ENDIF
ALLOCATE(fixedbcst(ndofst), fixedbcst0(ndofst))
fixedbcst = 0
fixedbcst0 = 0
DO ii = 1, nn
    IF (xcoord(ii,3) .EQ. 0._REKIND .OR. xcoord(ii,3) .EQ. le*10) THEN ! y = 0
        DO jj = 1, nnte
            DO kk = 1, ndofn
                fixedbcst( (ii-1)*ndofn+(jj-1)*ndof+kk ) = 1
                fixedbcst0( (ii-1)*ndofn+(jj-1)*ndof+kk ) = T0
            ENDDO
        ENDDO
    ENDIF
ENDDO

!Save the original K and M before applying BC
ki0 = ki
kj0 = kj
kv0 = kv
mi0 = mi
mj0 = mj
mv0 = mv

! Apply BC to K and M matrices
CALL csrapplybc(ki, kj, kv, ndof, knz, fixedbcst(1:ndof))
CALL csrapplybc(mi, mj, mv, ndof, mnz, fixedbcst(1:ndof))
ALLOCATE(rext(ndof), nodeshare(nn))
nodeshare = 0._REKIND ! number of elements that a node shared by
ntele = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        nodeshare(xconn(ii,jj)) = nodeshare(xconn(ii,jj)) + 1._REKIND
    ENDDO
ENDDO
DO ii = 1, nn
    IF (xcoord(ii, 3) .EQ. le) THEN
        ntele = ntele + nodeshare(ii)
    ENDIF
ENDDO
ftotal = q0*area
fnode = ftotal/ntele
rext = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        kk = xconn(ii,jj)
        IF (xcoord(kk, 3) .EQ. le) THEN ! z = 130
            rext(kk*ndofn-0) = rext(kk*ndofn-0) + fnode
        ENDIF
    ENDDO
ENDDO
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(2,2))
    !$ tusr(2,2) = omp_get_wtime()
    tcpu(2,3) = tcpu(2,2) - tcpu(2,1)
    tusr(2,3) = tusr(2,2) - tusr(2,1)
    !WRITE(*,105) '* Setting BCs  (sec)        ', tcpu(2,3)
ENDIF
!-------------------------------------------------------------------------------
! MUMPS: initialize, JOB = -1
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(3,1))
    !$ tusr(3,1) = omp_get_wtime()
    PRINT *,''
    PRINT *, '>>> MUMPS analysis and factorize'
ENDIF
mumps_par%COMM = MPI_COMM_WORLD
mumps_par%JOB = -1
mumps_par%SYM = 2 ! 0 unsymmetric, 1 SPD, 2 general symmetric
mumps_par%PAR = 1 ! 0 host not work, 1 host work
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
    WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        & "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        & "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
    GOTO 1500
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
! Out-of-core mode
!mumps_par%ICNTL(22) = 1 ! Out of core
!mumps_par%OOC_TMPDIR = 'ooc'
!mumps_par%OOC_PREFIX = 'ooc'

IF (MPI_ID .EQ. ROOT) THEN
    !PRINT *, 'MUMPS initialized'
ENDIF
IF (MPI_P .EQ. 1) THEN
    IF (MPI_ID .EQ. ROOT) THEN
!-------------------------------------------------------------------------------
! MUMPS: define problem on host thread
!-------------------------------------------------------------------------------
        mumps_par%N = ndof
        kk = 0
        DO ii = 1, ndof
            DO jj = ki(ii), ki(ii+1)-1
                IF (ii .GE. kj(jj)) THEN
                    kk = kk + 1
                ENDIF
            ENDDO
        ENDDO
        mumps_par%NZ = kk
        ALLOCATE(mumps_par%IRN(mumps_par%NZ))
        ALLOCATE(mumps_par%JCN(mumps_par%NZ))
        ALLOCATE(mumps_par%A(mumps_par%NZ))
        ALLOCATE(mumps_par%RHS(ndofst))
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
        mumps_par%NRHS = nnte         ! Multiple RHS
        mumps_par%LRHS = mumps_par%N    ! >= N
    ENDIF
ELSE
!-------------------------------------------------------------------------------
! MUMPS: define problem on all threads
!-------------------------------------------------------------------------------
    ! distribution of input matrix
    mumps_par%ICNTL(18) = 3
    mumps_par%N = ndof
    kk = 0
    DO ii = 1, ndof
        DO jj = ki(ii), ki(ii+1)-1
            IF (ii .GE. kj(jj)) THEN
                kk = kk + 1
            ENDIF
        ENDDO
    ENDDO
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
    IF (MPI_ID .EQ. ROOT) THEN
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
    GOTO 1500
ENDIF
IF (MPI_ID .EQ. ROOT) THEN
    !PRINT *, 'MUMPS analysis ... done'
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
    GOTO 1500
ENDIF
IF (MPI_ID .EQ. ROOT) THEN
    !PRINT *, 'MUMPS factorize ... done'
    CALL CPU_TIME(tcpu(3,2))
    !$ tusr(3,2) = omp_get_wtime()
    tcpu(3,3) = tcpu(3,2) - tcpu(3,1)
    tusr(3,3) = tusr(3,2) - tusr(3,1)
    WRITE(*,105) '* Time elapsed (sec)        ', tcpu(3,3)
ENDIF

!-------------------------------------------------------------------------------
! Time-loop
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    PRINT *,''
    PRINT *,'>>> Starting time-loop'
    OPEN(UNIT=iout, FILE='gmres.log', STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    iout2 = 500
    OPEN(UNIT=iout2, FILE='running.log', STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        WRITE(*,*) 'Error file: gmres.log'
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, MPI_IERR)
    ENDIF
ENDIF

ALLOCATE(dispst(ndofst), fextst(ndofst), fInit(ndofst), fBound(ndofst), tcoorde(nnte))
dt          = dt0*n1        ! initial time step
nip         = nipc*n1       ! intiial nubmer of temporal interpolation points
dispst      = 0._REKIND     ! Solution vector
tcoorde     = 0._REKIND
fInit       = 0._REKIND
fBound      = 0._REKIND

! Assiagn initial boundary condition
DO ii = 1, nn
    DO jj = 1, nnte
        DO kk = 1, ndofn
            !dispst( (ii-1)*ndofn+(jj-1)*ndof+kk ) = dispst( (ii-1)*ndofn+(jj-1)*ndof+kk ) + xcoord(ii,3)*(10-xcoord(ii,3))
        ENDDO
    ENDDO
ENDDO


! Get S-T matrix coefficients
ck = 0._REKIND
cm = 0._REKIND
cktn = 0._REKIND
cmtn = 0._REKIND
cktn1 = 0._REKIND
cmtn1 = 0._REKIND
CALL stcoeff(1, dt, angvel, ck, cm, cktn, cmtn, cktn1, cmtn1)
CALL mpikronspgemv(ck+cktn, cm+cmtn, 1._REKIND, 1._REKIND, nnte, ki0, kj0, &
        & kv0, mi0, mj0, mv0, ndof, dispst, fInit)
dispst = fixedbcst0
CALL mpikronspgemv(-ck-cktn, -cm-cmtn, 1._REKIND, 1._REKIND, nnte, ki0, kj0, &
        & kv0, mi0, mj0, mv0, ndof, dispst, fBound)

DO ii = 1, nstep ! Loop over space time slabs   
    ! Temporal nodal coordinates of iith S-T slab
    DO jj = 1, nnte
        tcoorde(jj) = tcoorde(nnte) + (jj-1)*dt/(nnte-1)
    ENDDO
    ! Calculate temporal component of external force vectors
    IF (ii .EQ. 1) THEN
       fextst = fInit
    ELSE
       fextst = 0._REKIND
    ENDIF
    CALL getfextc(ii, angvel, dt, fextc)
    kk = 0
    ll = 1 - ndof
    DO jj = 1, nnte
        kk = kk + ndof
        ll = ll + ndof
        ! Space time external force analagous vector
        fextst(ll:kk) = fextst(ll:kk) + fextc(jj)*rext
    ENDDO 
    ! Use KRONSPGEMV perform spgemv
    CALL mpikronspgemv(cktn1, cmtn1, 1._REKIND, 1._REKIND, nnte, ki, kj, &
        & kv, mi, mj, mv, ndof, dispst, fextst)
    fextst = fextst + fBound
    DO jj = 1, ndofst
        IF (fixedbcst(jj) .NE. 0) fextst(jj) = 0._REKIND
        !PRINT *, 'ID = ', jj, 'fextst = ', fextst(jj)
    ENDDO
    !PRINT *, 'ID = ', MPI_ID, 'sum(fextst) = ', sum(dabs(fextst))
    dispst = 0._REKIND
    CALL mpipgmres(ck+cktn, cm+cmtn, nnte, ki, kj, kv, mi, mj, mv, ndof, &
        &  mumps_par, fextst, dispst, im, eps, maxits, iout, ierr)
    DO jj = 1, ndofst
        IF (fixedbcst(jj) .NE. 0) dispst(jj) = fixedbcst0(jj)
        !PRINT *, 'ID = ', jj, 'disp = ', dispst(jj)
    ENDDO    
    !PRINT *, 'ID = ', MPI_ID, 'sum(dispst) = ', sum(dabs(dispst))
    IF (ierr .NE. 0) THEN
        PRINT *, 'Error: ', ierr
        GOTO 1000
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_Barrier error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    
    CALL MPI_BCAST(dispst, ndofst, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_BCAST error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    IF ((MOD(ii, nout) .EQ. 0)) THEN
        CALL mpipostprocess(job, mpi_ne, nn, ndofn, nnse, nnte, ngpe, &
            & nip, nipc, xconn(mpi_ele,:), tcoorde, dispst, MPI_ID)
    ENDIF
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Output results ... done, step = ', ii
        IF ((MOD(ii, 50) .EQ. 1) .OR. (ii .EQ. 1))  THEN
            WRITE(*,*) ''
            WRITE(*,110) 'STEP', 'TIME'
        ENDIF
        !DO jj = 1, nn
        !    IF (xcoord(jj,1) .EQ. 0._REKIND .AND. xcoord(jj,2) .EQ. 0._REKIND) THEN
        !        PRINT *, 'z= ', xcoord(jj,3),'t= ', tcoorde(2), 'disp = ', dispst(jj+nn)
        !    ENDIF
        !ENDDO
        WRITE(*,115) ii, tcoorde(nnte)
    ENDIF

ENDDO
! Clear memory and exit
DEALLOCATE(xcoord, xconn, dispst, fextst, fBound, fInit, tcoorde, dns, wqs)
DEALLOCATE(fixedbcst, fixedbcst0, rext)
DEALLOCATE(ki, kj, kv, mi, mj, mv)
DEALLOCATE(ki0, kj0, kv0, mi0, mj0, mv0)
DEALLOCATE(mpi_ele, mpi_gp)
IF (MPI_P .EQ. 1) THEN
    DEALLOCATE(mumps_par%IRN, mumps_par%JCN, mumps_par%A)
ELSE
    DEALLOCATE(mumps_par%IRN_loc, mumps_par%JCN_loc, mumps_par%A_loc)
ENDIF
IF (MPI_ID .EQ. ROOT) THEN
    DEALLOCATE(mumps_par%RHS)
ENDIF 
1000 CONTINUE
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(10,2))
    !$ tusr(10,2) = omp_get_wtime()
    tcpu(10,3) = tcpu(10,2) - tcpu(10,1)
    tusr(10,3) = tusr(10,2) - tusr(10,1)
    PRINT *, ''
    print *, '>>> Solution is done'
    WRITE(*,105) '* Total CPU time    (sec)   ', tcpu(10,3)
    WRITE(*,105) '* Total USR time    (sec)   ', tusr(10,3)
    PRINT *, ''
    PRINT *, '----------------------------------------------------------'
    PRINT *, ''
ENDIF
!  Destroy the MUMPS instance
mumps_par%JOB = -2
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL DMUMPS(mumps_par)
IF (mumps_par%INFOG(1).LT.0) THEN
    WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
        &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
    GOTO 1500
END IF
IF (MPI_ID .EQ. ROOT) THEN
    CLOSE(iout)
    CLOSE(120)
    CLOSE(125)
    CLOSE(130)
    CLOSE(iout2)
ENDIF
1500 CONTINUE
IF (MPI_ID .EQ. ROOT) THEN
    PRINT *, '>>> Program exit'
    PRINT *, ''
ENDIF
!CALL kokkos_finish()
CALL MPI_Finalize (MPI_IERR)
! FORMATS
100   FORMAT(A32, 8X, I16)
105   FORMAT(A32, 8X, F16.2)
110   FORMAT(A8,2X,7(A8, 2X))
115   FORMAT(I8,2X,F8.1,2X,6(F8.2,2X))
END PROGRAM main3d
