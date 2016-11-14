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
REAL(KIND=REKIND) :: ym, nu, rho, le, area, freq
INTEGER :: ne, nn, ndofn, nnsd, nnte, nquads, nquadt, ndof, ndofst, nnse, &
    & ngp, ngpe, ncmp, n1, n2
REAL(KIND=REKIND) :: dt, dt0, p0
INTEGER :: nstep, nip, nipc, nout
INTEGER :: errstat
CHARACTER(LEN=80) :: job, inp, dat, temp, errmesg
! Mesh
INTEGER, DIMENSION(:,:), ALLOCATABLE :: xconn
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: xcoord
! Shape functions
REAL(KIND=REKIND), DIMENSION(6,6) :: stiff
REAL(KIND=REKIND), DIMENSION(6,9) :: sd
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: wqs
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: dns, ns, dnxyz
! Spatial matrices
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nx, Ny, Nz
INTEGER :: knz, mnz
INTEGER, DIMENSION(:), POINTER :: ki, kj, mi, mj, kit, kjt, mit, mjt
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: mv, kv, mvt, kvt
! Space-time coefficients
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: tcoorde
REAL(KIND=REKIND), DIMENSION(6,6) :: ck, cm, cktn, cmtn, cktn1, cmtn1
! Boundary conditions
REAL(KIND=REKIND) :: ntele, ftotal, fnode
INTEGER, DIMENSION(:), ALLOCATABLE :: fixedbcst, fixedbcst0 
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: rext, nodeshare
! Linear system solver
INTEGER :: im, maxits, iout, iout2, ierr
REAL(KIND=REKIND) :: eps
! Time-loop
REAL(KIND=REKIND) :: angvel, dmax, wsmax
REAL(KIND=REKIND), DIMENSION(6) :: fextc
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: dispst, fextst
! Two-scale damage model
REAL(KIND=REKIND) :: Dc
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: p, D, ws, del_po, sig_D
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: strain_tgo, strain_pl, &
    & sig_eff, sig_back
! Status flag
INTEGER :: nef, ndf, nefe, itfail, itfail0
INTEGER, DIMENSION(:), ALLOCATABLE :: necnt, necnt0, gp_bc, gpsta, elesta, &
    & nefail, ndfail
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
CALL kokkos_init()
!-------------------------------------------------------------------------------
! Set parameters
!-------------------------------------------------------------------------------
! Material, isotropic elastic only
ym      = 1.97e5        ! Young's modulus
nu      = 0.3_REKIND    ! Poisson's ratio
rho     = 7.860e-9      ! Mass density
! Damage model
Dc      = 0.3_REKIND    ! Damage criteria
! Input file
job     = 'sent'
inp     = TRIM(job)//'.inp'
dat     = TRIM(job)//'.dat'
WRITE(temp,'(i8)') MPI_ID
mpi_job = TRIM(job)//'_'//TRIM(ADJUSTL(temp))
! Geometry
le      = 130._REKIND    ! Height of the plate (loading side)
! Mesh
ndofn   = 3             ! Number of spatial degrees of freedom per node
nnsd    = 2             ! Order of spatial shape functions
nnte    = 3             ! Order of temporal shape functions
nquadt  = 2             ! Number of quadrature points temporally
nquads  = 2             ! Number of quadrature points spatially
ncmp    = 6             ! Number of components in stress/strain tensor
nnse    = nnsd**ndofn   ! Number of nodes per spatial element
ngpe    = nquads**ndofn ! Number of Gauss points per spatial element
! Time step
dt0     = 0.05_REKIND   ! Step time
nstep   = 30000            ! Number of steps
n1      = 100           ! Initial time step (cycle)
n2      = 5            ! Time step after crack initiation (cycle)
nipc    = 256           ! Number of temporal interpolation points per cycle
! Applid load
p0      = 70.0_REKIND   ! Pressure amplitude (Traction load)
freq    = 20._REKIND    ! Frequency, in Hz
angvel  = 2._REKIND*pi*freq
area    = 664._REKIND    ! Traction area
! Linear system solver
im      = 10            ! GMRES: size of Krylov subspace 
eps     = 1.d-8         ! GMRES: tolerance for stopping criterion
maxits  = 50            ! GMRES: maximum number of iterations allowed 
iout    = 40            ! GMRES: solution info output file number
! Post-processing
nout    = 100             ! Output frequency

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
    ALLOCATE(xcoord(nn, ndofn))
    DO i = 1, nn
        READ(100, *) ii, xcoord(ii, 1:ndofn) ! Read node coordinates
    ENDDO
    READ(100, *) ne ! Read the number of elements 
    ALLOCATE(xconn(ne,nnse))
    DO i = 1, ne
        READ(100, *) ii, xconn(ii, 1:nnse) ! Read element connectivities
    ENDDO
    CLOSE(100)
    !print *, 'Number of Groups = ', MPI_P
    !print *, 'Size of Group = ', ceiling(REAL(ne)/REAL(MPI_P))
    !j = ceiling(REAL(ne)/REAL(MPI_P))
    !DO i = 1, MPI_P
    !    k = i*j
    !    if (k.gt.ne) k = ne
    !    PRINT *, 'GROUP',i,'Elements from', (i-1)*j+1, 'to', k
    !ENDDO
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
ndofst  = 2*nnte*ndof   ! Number of space-time DOFs
ngp     = ngpe*ne       ! Number of spatial Gauss points
IF (MPI_ID .NE. ROOT) THEN
    ALLOCATE(xcoord(nn, ndofn), xconn(ne,nnse))
ENDIF
CALL MPI_BCAST(xcoord, nn*ndofn, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
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
! Calculate the number of elements to which that a node belongs, it will be
! used to determine when a node is disconnected with others
ALLOCATE(necnt(nn), necnt0(nn))
necnt = 0
DO i = 1, mpi_ne
    k = mpi_ele(i)
    DO j = 1, nnse
        necnt(xconn(k,j)) = necnt(xconn(k,j)) + 1
    ENDDO
ENDDO
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
        WRITE(110,'(3(E20.11E3,2x))') (xcoord(ii,jj),jj=1,ndofn)
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

!PRINT *, 'ID = ', MPI_ID, 'from ', mpi_ele(1), 'to', mpi_ele(mpi_ne)
!PRINT *, 'ID = ', MPI_ID, 'from ', mpi_gp(1), 'to', mpi_gp(mpi_ngp)
!CALL pause_exe()
!-------------------------------------------------------------------------------
! Set Gauss Points and shape functions
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(1,1))
    !$ tusr(1,1) = omp_get_wtime()
    PRINT *,''
    PRINT *, '>>> Assembling spatial matrices'
ENDIF
ALLOCATE(dns(ndofn*ngpe, nnse), dnxyz(ndofn**2*ngpe, ndofn*nnse), &
    & wqs(nquads), ns(ndofn*ngpe, ndofn*nnse))
CALL setshapefunc(ndofn, nnsd, nnse, nquads, ngpe, wqs, ns, dns, dnxyz)
!-------------------------------------------------------------------------------
! Formulate spatial matrices (spatial integration)
!-------------------------------------------------------------------------------
knz = mpi_ne*(nnse*ndofn)**2
mnz = knz
ALLOCATE(Nx(nnse, mpi_ngp), Ny(nnse, mpi_ngp), Nz(nnse, mpi_ngp), &
    & kit(MAX(knz,ndof+1)), kjt(knz), kvt(knz), mit(MAX(mnz,ndof+1)), &
    & mjt(mnz), mvt(mnz))
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
 CALL mpiglbmtx(mpi_ne, mpi_ele, ne, nn, ndofn, nnse, nquads, ngpe, ncmp, &
    & xconn, xcoord, ym, nu, rho, wqs, ns, dns, dnxyz, sd, stiff, Nx, Ny, &
    & Nz, kit(1:knz), kjt, kvt, knz, mit(1:mnz), mjt, mvt, mnz)
 CALL tricsr(kit, kjt, kvt, ndof, knz, MAX(knz,ndof+1))
 CALL tricsr(mit, mjt, mvt, ndof, mnz, MAX(mnz,ndof+1))
!PRINT *, 'ID = ', MPI_ID, 'knz = ', kit(ndof+1)
!GOTO 1500
ALLOCATE(ki(ndof+1), kj(knz), kv(knz), mi(ndof+1), mj(mnz), mv(mnz))
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

!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!PRINT *, 'ID = ', MPI_ID, 'sum(mv) = ', sum(mv)
!CALL pause_exe()
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
DO ii = 1, nn
    IF (xcoord(ii,3) .EQ. 0._REKIND) THEN ! y = 0
        DO jj = 1, 2*nnte
            fixedbcst((ii-1)*ndofn+(jj-1)*ndof+1 : ii*ndofn+(jj-1)*ndof) &
                & = [1, 1, 1]
        ENDDO
    ENDIF
ENDDO
fixedbcst0 = fixedbcst
! Apply BC to K and M matrices
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL csrapplybc(ki, kj, kv, ndof, knz, fixedbcst(1:ndof))
CALL csrapplybc(mi, mj, mv, ndof, mnz, fixedbcst(1:ndof))
!temp = 'stiff_bc'//trim(char(mpi_id+49))
!CALL outputcsrmat(ki, kj, kv, ndof, knz, temp)
!call pause_exe()
!temp = 'mass_bc_'//trim(char(mpi_id+49))
!CALL outputcsrmat(mi, mj, mv, ndof, mnz, temp)
! Traction load (uniformly distributed constant)
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
ftotal = p0*area
fnode = ftotal/ntele
rext = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        kk = xconn(ii,jj)
        IF (xcoord(kk, 3) .EQ. le) THEN ! y = 10
            rext(kk*3-0) = rext(kk*3-0) + fnode
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
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!PRINT *, 'ID = ', MPI_ID, 'sum(kv) = ', sum(kv(1:knz))
!CALL pause_exe()
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
        mumps_par%NRHS = 2*nnte         ! Multiple RHS
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
        mumps_par%NRHS = 2*nnte         ! Multiple RHS
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
ALLOCATE(dispst(ndofst), fextst(ndofst), tcoorde(nnte))
ALLOCATE(p(mpi_ngp), D(mpi_ngp), ws(mpi_ngp), del_po(mpi_ngp), &
    & sig_D(mpi_ngp), strain_tgo(6,mpi_ngp), strain_pl(6,mpi_ngp), &
    & sig_eff(6,mpi_ngp), sig_back(6,mpi_ngp), elesta(mpi_ne), &
    & gpsta(mpi_ngp), nefail(mpi_ne), ndfail(nn))
dt          = dt0*n1        ! initial time step
nip         = nipc*n1       ! intiial nubmer of temporal interpolation points
itfail      = 0             ! > 0 if element(s) failed in current step, else 0
gpsta       = 1             ! Gauss point status, 0 = failed, 1 = active
elesta      = 1             ! Element status, 0 = failed, 1 = active
dispst      = 0._REKIND     ! Solution vector
p           = 0._REKIND     ! P (??) of last step
D           = 0._REKIND     ! Damage of last step
ws          = 0._REKIND     ! Stored energy of last step
del_po      = 0._REKIND     ! Delta_P of last step
sig_D       = 0._REKIND     ! Sig_D (??) of last step
strain_tgo  = 0._REKIND     ! Total global strain of last step
strain_pl   = 0._REKIND     ! Total plastic strain of last step
sig_eff     = 0._REKIND     ! Effective stress of last step
sig_back    = 0._REKIND     ! Backstress of last step
tcoorde     = 0._REKIND
DO ii = 1, nstep ! Loop over space time slabs   
    IF (MPI_ID .EQ. ROOT) THEN
        CALL CPU_TIME(tcpu(6,1))
        !$ tusr(6,1) = omp_get_wtime()
        WRITE(iout2,*) 'Check damage ... start, step = ', ii
    ENDIF
    ! Check damage first
    itfail = 0
    itfail0 = 0
    ndf = 0
    nef = 0
    ndfail = 0
    nefail = 0
    necnt0 = 0
    IF (ANY(D(1:mpi_ngp) .GE. Dc)) THEN 
        itfail0 = 1
        ! Check element and Gauss point status
        CALL chkelesta(nn, mpi_ne, nnse, ngpe, Dc, D, elesta, gpsta, nef, &
            & nefail)
    ENDIF
    ! Broadcast itfail to all MPI threads by MPI_ALLREDUCE
    CALL MPI_ALLREDUCE(itfail0, itfail, 1, MPI_INTEGER, MPI_SUM, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_ALLREDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    !PRINT *, 'ID = ', MPI_ID, 'ndf = ', ndf, 'nef = ', nef
    IF (itfail .NE. 0) THEN     
        ! Delete failed elements
        ! 1st step, delete the failed elements
        CALL delelement(nn, mpi_ne, ndofn, nnse, ngpe, ncmp, nquads, &
            & xconn(mpi_ele,:), xcoord, dns, dnxyz, sd, stiff, wqs, nef, &
            & nefail, ki, kj, kv, ndof, knz)
        !temp = 'stiff_bc'//trim(char(mpi_id+49))
        !CALL outputcsrmat(ki, kj, kv, ndof, knz, temp)
        ! Check connectivity status
        CALL chkcnsta(nn, mpi_ne, nnse, xconn(mpi_ele,:), nef, nefail, necnt)
        ! Collect necnt all together  
        CALL MPI_ALLREDUCE(necnt, necnt0, nn, MPI_INTEGER, MPI_SUM, &
            & MPI_COMM_WORLD, MPI_IERR)
        IF (MPI_IERR .NE. 0) THEN
            PRINT *, 'MPI_ALLREDUCE error = ',  MPI_IERR
            GOTO 1500
        ENDIF
        ! Check node status
        CALL chkndsta(nn, mpi_ne, nnse, xconn(mpi_ele,:), nef, nefail, necnt0, &
            & ndf, ndfail)
        ! 2nd step, add diagonal entries to K for disconnected nodes
        CALL adddiagonal(nnse, ndofn, ndf, ndfail, ki, kj, kv, ndof, knz)
        !temp = 'stiff_bc'//trim(char(mpi_id+49))
        !CALL outputcsrmat(ki, kj, kv, ndof, knz, temp)
        !call pause_exe()
        ! 3rd step, constrain the disconnected nodes
        !PRINT *, 'ID = ', MPI_ID, 'nef = ', nef, 'ef = ', mpi_ele(nefail(1:nef))
        IF (ndf.GT.0) THEN
        !PRINT *, 'ID = ', MPI_ID, 'ndf = ', ndf, 'df = ', ndfail(1:ndf)
            DO i = 1, ndf
                j = ndfail(i)
                DO k = 1, 2*nnte
                    fixedbcst0((j-1)*ndofn+(k-1)*ndof+1 : j*ndofn+(k-1)*ndof) &
                        & = [1, 1, 1]
                ENDDO
            ENDDO
        ENDIF
        ! update bc
        CALL MPI_ALLREDUCE(fixedbcst0, fixedbcst, ndofst, MPI_INTEGER, &
            & MPI_SUM, MPI_COMM_WORLD, MPI_IERR)
        IF (MPI_IERR .NE. 0) THEN
            PRINT *, 'MPI_ALLREDUCE error = ',  MPI_IERR
            GOTO 1500
        ENDIF
        fixedbcst0 = fixedbcst
        CALL csrapplybc(mi, mj, mv, ndof, mnz, fixedbcst(1:ndof))
        CALL csrapplybc(ki, kj, kv, ndof, knz, fixedbcst(1:ndof))
        ! 4th step, generate new preconditioner by MUMPS
        ! 4.1 copy K matrix
        k = 0
        DO i = 1, ndof
            DO j = ki(i), ki(i+1)-1
                IF (i .GE. kj(j)) THEN
                    k = k + 1
                ENDIF
            ENDDO
        ENDDO
        IF (MPI_P .EQ. 1) THEN
            mumps_par%NZ = k
            !WRITE(*,*) 'NZ = ', K
            DEALLOCATE(mumps_par%IRN, mumps_par%JCN, mumps_par%A)
            ALLOCATE(mumps_par%IRN(mumps_par%NZ))
            ALLOCATE(mumps_par%JCN(mumps_par%NZ))
            ALLOCATE(mumps_par%A(mumps_par%NZ))
            k = 1
            DO i = 1, ndof
                DO j = ki(i), ki(i+1)-1
                    IF (i .GE. kj(j)) THEN
                        mumps_par%IRN(k) = i
                        mumps_par%JCN(k) = kj(j)
                        mumps_par%A(k) = kv(j)
                        k = k + 1
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            mumps_par%NZ_loc = k
            !WRITE(*,*) 'NZ = ', K
            DEALLOCATE(mumps_par%IRN_loc, mumps_par%JCN_loc, mumps_par%A_loc)
            ALLOCATE(mumps_par%IRN_loc(mumps_par%NZ_loc))
            ALLOCATE(mumps_par%JCN_loc(mumps_par%NZ_loc))
            ALLOCATE(mumps_par%A_loc(mumps_par%NZ_loc))
            k = 1
            DO i = 1, ndof
                DO j = ki(i), ki(i+1)-1
                    IF (i .GE. kj(j)) THEN
                        mumps_par%IRN_loc(k) = i
                        mumps_par%JCN_loc(k) = kj(j)
                        mumps_par%A_loc(k) = kv(j)
                        k = k + 1
                    ENDIF
                ENDDO
            ENDDO 
        ENDIF
	! 4.2 analysis
        mumps_par%JOB = 1
        !CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
        CALL DMUMPS(mumps_par)
        IF (mumps_par%INFOG(1).LT.0) THEN
            WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
                &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
            GOTO 1500
        ENDIF
        ! 4.3 factorize
        mumps_par%JOB = 2
        !CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
        CALL DMUMPS(mumps_par)
        IF (mumps_par%INFOG(1).LT.0) THEN
            WRITE(*,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
                &   "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
                &   "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
            GOTO 1500
        ENDIF
        ! Reduce time step
        dt = dt0*n2
        nip = nipc*n2
    ENDIF
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Check damage ... done, step = ', ii
        CALL CPU_TIME(tcpu(6,2))
        !$ tusr(6,2) = omp_get_wtime()
    ENDIF
    !CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    !IF (MPI_IERR .NE. 0) THEN
    !    PRINT *, 'MPI_Barrier error = ',  MPI_IERR
    !    GOTO 1500
    !ENDIF
    ! Temporal nodal coordinates of iith S-T slab
    DO jj = 1, nnte
        tcoorde(jj) = tcoorde(nnte) + (jj-1)*dt/(nnte-1)
    ENDDO
    ! Calculate temporal component of external force vectors
    fextc = 0._REKIND
    fextst = 0._REKIND
    CALL getfextc(ii, angvel, dt, fextc)
    kk = 0
    ll = 1 - ndof
    DO jj = 1, 2*nnte
        kk = kk + ndof
        ll = ll + ndof
        ! Space time external force analagous vector
        fextst(ll:kk) = fextc(jj)*rext
    ENDDO 
    ! Get S-T matrix coefficients
    ck = 0._REKIND
    cm = 0._REKIND
    cktn = 0._REKIND
    cmtn = 0._REKIND
    cktn1 = 0._REKIND
    cmtn1 = 0._REKIND     
    CALL stcoeff(ii, dt, angvel, ck, cm, cktn, cmtn, cktn1, cmtn1)
    ! Use KRONSPGEMV perform spgemv
    CALL mpikronspgemv(cktn1, cmtn1, 1._REKIND, 1._REKIND, 2*nnte, ki, kj, &
        & kv, mi, mj, mv, ndof, dispst, fextst)
    !PRINT *, 'ID = ', MPI_ID, 'sum(fextst) = ', sum(dabs(fextst))
    DO jj = 1, ndofst
        IF (fixedbcst(jj) .NE. 0) fextst(jj) = 0._REKIND 
    ENDDO
    dispst = 0._REKIND
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Solve linear system ... start, step = ', ii
        CALL CPU_TIME(tcpu(4,1))
        !$ tusr(4,1) = omp_get_wtime()
    ENDIF
    !CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    !IF (MPI_IERR .NE. 0) THEN
    !    PRINT *, 'MPI_Barrier error = ',  MPI_IERR
    !    GOTO 1500
    !ENDIF

    CALL mpipgmres(ck+cktn, cm+cmtn, 2*nnte, ki, kj, kv, mi, mj, mv, ndof, &
        &  mumps_par, fextst, dispst, im, eps, maxits, iout, ierr)
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Solve linear system ... done, step = ', ii
        CALL CPU_TIME(tcpu(4,2))
        !$ tusr(4,2) = omp_get_wtime()
        IF (ierr/=0) THEN
            WRITE(*,*) 'GMRES error: ', ierr
        ENDIF
    ENDIF
    IF (ierr .NE. 0) GOTO 1000
    IF (MPI_ID .EQ. ROOT) THEN
        DO jj = 1, ndofst
            IF (fixedbcst(jj) .NE. 0) dispst(jj) = 0._REKIND 
        ENDDO
        !temp = 'dispst'
        !CALL outputrealvec(dispst, ndofst, temp)
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
    IF (MPI_ID .EQ. ROOT) THEN
        ! post-process
        CALL CPU_TIME(tcpu(5,1))
        !$ tusr(5,1) = omp_get_wtime()
        WRITE(iout2,*) 'Solve damage ... start, step = ', ii
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    CALL kokkos_damage(MPI_ID, job, mpi_ne, nn, ndofn, nnse, nnte, ngpe, &
            & ncmp, nip, xconn(mpi_ele,:), Nx, Ny, Nz, ym, nu, angvel, &
            & tcoorde, dispst, gpsta, strain_tgo, strain_pl, sig_eff, &
            & sig_back, p, D, ws, del_po, sig_D)
    
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_Barrier error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    !PRINT *, 'ID = ', MPI_ID, 'max(D) = ', maxval(D)
    dmax = 0._REKIND
    wsmax = 0._REKIND
    CALL MPI_REDUCE(MAXVAL(D), dmax, 1, MPI_REAL8, MPI_MAX, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL MPI_REDUCE(MAXVAL(ws), wsmax, 1, MPI_REAL8, MPI_MAX, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Solve damage ... done, step = ', ii
        CALL CPU_TIME(tcpu(5,2))
        !$ tusr(5,2) = omp_get_wtime()
        CALL CPU_TIME(tcpu(7,1))
        !$ tusr(7,1) = omp_get_wtime()
        WRITE(iout2,*) 'Output results ... start, step = ', ii
    ENDIF
    IF ((MOD(ii, nout) .EQ. 0) .OR. (itfail .NE. 0)) THEN
        CALL mpipostprocess(job, mpi_ne, nn, ndofn, nnse, nnte, ngpe, ncmp, &
            & nip, nipc, xconn(mpi_ele,:), Nx, Ny, Nz, ym, nu, angvel, &
            & tcoorde, dispst, gpsta, elesta, D, ws, MPI_ID)
    ENDIF
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Output results ... done, step = ', ii
        CALL CPU_TIME(tcpu(7,2))
        !$ tusr(7,2) = omp_get_wtime()
        tcpu(4,3) = tcpu(4,2) - tcpu(4,1)
        tusr(4,3) = tusr(4,2) - tusr(4,1)
        tcpu(5,3) = tcpu(5,2) - tcpu(5,1)
        tusr(5,3) = tusr(5,2) - tusr(5,1)
        tcpu(6,3) = tcpu(6,2) - tcpu(6,1)
        tusr(6,3) = tusr(6,2) - tusr(6,1)
        tcpu(7,3) = tcpu(6,2) - tcpu(6,1)
        tusr(7,3) = tusr(7,2) - tusr(7,1) 
        IF ((MOD(ii, 50) .EQ. 1) .OR. (ii .EQ. 1))  THEN
            WRITE(*,*) ''
            WRITE(*,110) 'STEP', 'TIME', 'SOLVER', 'DAMAGE', 'RMELE', & 
                & 'POSTPR', 'max(Ws)', 'max(D)'
        ENDIF
        WRITE(*,115) ii, tcoorde(nnte), tusr(4:7,3), wsmax, dmax
    ENDIF

    !CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    !IF (MPI_IERR .NE. 0) THEN
    !    PRINT *, 'MPI_Barrier error = ',  MPI_IERR
    !    GOTO 1500
    !ENDIF
ENDDO
! Clear memory and exit
DEALLOCATE(xcoord, xconn, dispst, fextst, tcoorde, dns, wqs)
DEALLOCATE(fixedbcst, fixedbcst0, rext)
DEALLOCATE(ki, kj, kv, mi, mj, mv)
DEALLOCATE(mpi_ele, mpi_gp)
DEALLOCATE(p, D, ws, del_po, sig_D, strain_tgo, strain_pl, sig_eff, sig_back, &
    & elesta, gpsta, necnt, necnt0, nefail, ndfail)
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
CALL kokkos_finish()
CALL MPI_Finalize (MPI_IERR)
! FORMATS
100   FORMAT(A32, 8X, I16)
105   FORMAT(A32, 8X, F16.2)
110   FORMAT(A8,2X,7(A8, 2X))
115   FORMAT(I8,2X,F8.1,2X,6(F8.2,2X))
END PROGRAM main3d
