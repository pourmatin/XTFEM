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
REAL(KIND=REKIND) :: kappa, cp, rho, ym0, nu0, le, area, alpha
INTEGER :: ne, nn, nnsd, nnte, nquads, nquadt, nnse, ngp, ngpe, ncmp, n1, n2
REAL(KIND=REKIND) :: dt, dt0
INTEGER :: nstep, nip, nipc, nout
INTEGER :: errstat
CHARACTER(LEN=80) :: job, inp, dat, temp, errmesg
! Mesh
INTEGER, DIMENSION(:,:), ALLOCATABLE :: xconn
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: xcoord
! Shape functions
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: wqs
! Space-time coefficients
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: tcoorde
! Linear system solver
INTEGER :: im, maxits, iout, iout2, ierr
REAL(KIND=REKIND) :: eps
! Counters
INTEGER :: i, j, k, l, m, ii, jj, kk, ll, mm, pp
! Performance statistics
REAL(KIND=REKIND) :: mem(10), tcpu(10, 3), tusr(10, 3)
! MPI
CHARACTER(LEN=80) :: mpi_job
INTEGER :: MPI_IERR, MPI_P, MPI_ID, ROOT
INTEGER :: mpi_ne, mpi_ngp, fid
INTEGER, DIMENSION(:), ALLOCATABLE :: mpi_ele, mpi_gp

!Mechanical Parameters
INTEGER :: ndofn_m, ndof_m, ndofst_m
REAL(KIND=REKIND) :: p0, freq_m
! Shape functions
REAL(KIND=REKIND), DIMENSION(:,:,:), ALLOCATABLE :: stiff
REAL(KIND=REKIND), DIMENSION(6,9) :: sd
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: dns_m, ns_m, dnxyz
! Spatial matrices
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nx, Ny, Nz
! Spatial matrices
INTEGER :: knz_m, mnz_m
INTEGER, DIMENSION(:), POINTER :: ki_m, kj_m, mi_m, mj_m, kit_m, kjt_m, mit_m, mjt_m
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: mv_m, kv_m, mvt_m, kvt_m
! Space-time coefficients
REAL(KIND=REKIND), DIMENSION(6,6) :: ck_m, cm_m, cktn_m, cmtn_m, cktn1_m, cmtn1_m
! Boundary conditions
REAL(KIND=REKIND) :: ntele_m, ftotal_m, fnode_m
INTEGER, DIMENSION(:), ALLOCATABLE :: fixedbcst_m, fixedbcst0_m
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: rext_m, nodeshare_m
! Linear system solver
INTEGER :: im_m, maxits_m
REAL(KIND=REKIND) :: eps_m
! Time-loop
REAL(KIND=REKIND) :: angvel_m, dmax, wsmax, Umin, Umax, Tmin, Tmax
REAL(KIND=REKIND), DIMENSION(6) :: fextc_m, fint_m
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: dispst, fextst_m, ym, nu
! Two-scale damage model
REAL(KIND=REKIND) :: Dc
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: p, D, ws, del_po, sig_D
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: strain_tgo, strain_pl, &
    & sig_eff, sig_back
! Status flag
INTEGER :: nef, ndf, nefe, itfail, itfail0
INTEGER, DIMENSION(:), ALLOCATABLE :: necnt, necnt0, gp_bc, gpsta, elesta, &
    & nefail, ndfail
! MUMPS
TYPE (DMUMPS_STRUC) mumps_par_m

! Temperature Parameters
INTEGER :: ndofn_t, ndof_t, ndofst_t
REAL(KIND=REKIND) :: q0, freq_t
! Shape functions
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: ns_t
REAL(KIND=REKIND), DIMENSION(:,:,:), ALLOCATABLE :: dns_t
! Spatial matrices
INTEGER :: knz_t, mnz_t
INTEGER, DIMENSION(:), POINTER :: ki_t, kj_t, mi_t, mj_t, kit_t, kjt_t, mit_t, mjt_t, ki0, kj0, mi0, mj0
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: mv_t, kv_t, mvt_t, kvt_t, mv0, kv0
! Space-time coefficients
REAL(KIND=REKIND), DIMENSION(3,3) :: ck_t, cm_t, cktn_t, cmtn_t, cktn1_t, cmtn1_t
! Boundary conditions
REAL(KIND=REKIND) :: ntele_t, ftotal_t, fnode_t, T0, Tref
INTEGER, DIMENSION(:), ALLOCATABLE :: fixedbcst_t, fixedbcst0_t
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: rext_t, nodeshare_t
! Thermal strain
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE ::vv, ve
! Linear system solver
INTEGER :: im_t, maxits_t
REAL(KIND=REKIND) :: eps_t
! Time-loop
REAL(KIND=REKIND), DIMENSION(3) :: fextc_t
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: tempst, fextst_t, fBound, fInit
! MUMPS
TYPE (DMUMPS_STRUC) mumps_par_t
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
CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
!-------------------------------------------------------------------------------
! Initialize KOKKOS
!-------------------------------------------------------------------------------
CALL kokkos_init()
!-------------------------------------------------------------------------------
! Set parameters
!-------------------------------------------------------------------------------
! Material, isotropic elastic only
ym0     = 1.97e5                     ! Initial Young's modulus N/mm^2 (MPa)
nu0     = 0.3_REKIND                 ! Initial Poisson's ratio
rho     = 7.860e-6                   ! Mass density Kg/mm^3
kappa   = .1_REKIND                  ! heat conductivity W/mm.T
cp      = 452._REKIND                ! specific heat capacity J/Kg.T
alpha   = 1.2e-5                     ! thermal expansion coefficient 1/T
! Damage model
Dc      = 0.3_REKIND    ! Damage criteria
! Input file
job     = 'sent'
inp     = TRIM(job)//'.inp'
dat     = TRIM(job)//'.dat'
WRITE(temp,'(i8)') MPI_ID
mpi_job = TRIM(job)//'_'//TRIM(ADJUSTL(temp))
! Geometry
le      = 10._REKIND    ! Height of the plate (loading side)
! Mesh
ndofn_t = 1             ! Number of spatial degrees of freedom per node
ndofn_m = 3
nnsd    = 2              ! Order of spatial shape functions
nnte    = 3               ! Order of temporal shape functions
nquadt  = 2             ! Number of quadrature points temporally
nquads  = 2            ! Number of quadrature points spatially
ncmp    = 6             ! Number of components in stress/strain tensor
nnse    = 8              ! Number of nodes per spatial element
ngpe    = 8              ! Number of Gauss points per spatial element
! Time step
dt0     = 0.05_REKIND   ! Step time
nstep   = 300000            ! Number of steps
n1      = 10           ! Initial time step (cycle)
n2      = 1            ! Time step after crack initiation (cycle)
nipc    = 256           ! Number of temporal interpolation points per cycle
! Applid load
p0      = 70._REKIND       ! Pressure amplitude (Traction load)
freq_m    = 20._REKIND    ! Frequency, in Hz
angvel_m  = 2._REKIND*pi*freq_m
! Applid heat flux
q0      = -0.3_REKIND   ! uniformly distributed heat flux W/mm^2
freq_t    = 0.0_REKIND    ! Frequency, in Hz
area    = 50._REKIND    ! Traction area mm^2
T0      = 0._REKIND      ! Boundary temprature
Tref    = 0._REKIND      ! Reference tempereture
! Linear system solver
im_m      = 100                      ! GMRES: size of Krylov subspace
eps_m     = 1.d-10               ! GMRES: tolerance for stopping criterion
maxits_m  = 500                   ! GMRES: maximum number of iterations allowed
iout    = 40                          ! GMRES: solution info output file number
! Linear system solver
im_t      = 100            ! GMRES: size of Krylov subspace 
eps_t     = 1.d-10         ! GMRES: tolerance for stopping criterion
maxits_t  = 500            ! GMRES: maximum number of iterations allowed 
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
CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
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
ndof_m    = ndofn_m*nn          ! Number of spatial DOFs for the Mechanical part
ndofst_m  = 2*nnte*ndof_m     ! Number of space-time DOFs for the Mechanical part
ndof_t    = ndofn_t*nn              ! Number of spatial DOFs for the Temperature part
ndofst_t  = nnte*ndof_t            ! Number of space-time DOFs for the Temperature part
ngp     = ngpe*ne       ! Number of spatial Gauss points
IF (MPI_ID .NE. ROOT) THEN
    ALLOCATE(xcoord(nn, 3), xconn(ne,nnse))
ENDIF
CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
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
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
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
    WRITE(*,100) '* Number of spatial DOFs    ', ndof_m
    WRITE(*,100) '* Number of space-time DOFs ', ndofst_m
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
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
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
CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
ALLOCATE(dns_t(ngpe, nnse, 3), wqs(nquads), ns_t(ngpe, nnse))
CALL setshapefunc_t(nnsd, nnse, nquads, ngpe, wqs, ns_t, dns_t)
ALLOCATE(dns_m(ndofn_m*ngpe, nnse), dnxyz(ndofn_m**2*ngpe, ndofn_m*nnse), ns_m(ndofn_m*ngpe, ndofn_m*nnse))
CALL setshapefunc_m(ndofn_m, nnsd, nnse, nquads, ngpe, wqs, ns_m, dns_m, dnxyz)
!-------------------------------------------------------------------------------
! Formulate spatial matrices for Temperature part (spatial integration)
!-------------------------------------------------------------------------------
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
knz_t = mpi_ne*(nnse*ndofn_t)**2
mnz_t = knz_t
ALLOCATE(kit_t(MAX(knz_t,ndof_t+1)), kjt_t(knz_t), kvt_t(knz_t), mit_t(MAX(mnz_t,ndof_t+1)), &
    & mjt_t(mnz_t), mvt_t(mnz_t))
 CALL mpiglbmtx_t(mpi_ne, mpi_ele, ne, nn, ndofn_t, nnse, nquads, ngpe, &
    & xconn, xcoord, kappa, cp, rho, wqs, ns_t, dns_t, kit_t(1:knz_t), &
    & kjt_t, kvt_t, knz_t, mit_t(1:mnz_t), mjt_t, mvt_t, mnz_t)
 CALL tricsr(kit_t, kjt_t, kvt_t, ndof_t, knz_t, MAX(knz_t,ndof_t+1))
 CALL tricsr(mit_t, mjt_t, mvt_t, ndof_t, mnz_t, MAX(mnz_t,ndof_t+1))
ALLOCATE(ki_t(ndof_t+1), kj_t(knz_t), kv_t(knz_t), mi_t(ndof_t+1), mj_t(mnz_t), mv_t(mnz_t))
ALLOCATE(ki0(ndof_t+1), kj0(knz_t), kv0(knz_t), mi0(ndof_t+1), mj0(mnz_t), mv0(mnz_t))
ki_t = kit_t(1:(ndof_t+1))
kj_t = kjt_t(1:knz_t)
kv_t = kvt_t(1:knz_t)
mi_t = mit_t(1:(ndof_t+1))
mj_t = mjt_t(1:mnz_t)
mv_t = mvt_t(1:mnz_t)
!-------------------------------------------------------------------------------
! Formulate spatial matrices for Mechanical part (spatial integration)
!-------------------------------------------------------------------------------
knz_m = mpi_ne*(nnse*ndofn_m)**2
mnz_m = knz_m
ALLOCATE(ym(mpi_ne*ngpe), nu(mpi_ne*ngpe), stiff(ncmp,ncmp,mpi_ngp), tempst(ndofst_t)) 
ym     = ym0
nu     = nu0
tempst = Tref
ALLOCATE(Nx(nnse, mpi_ngp), Ny(nnse, mpi_ngp), Nz(nnse, mpi_ngp), &
    & kit_m(MAX(knz_m,ndof_m+1)), kjt_m(knz_m), kvt_m(knz_m), mit_m(MAX(mnz_m,ndof_m+1)), &
    & mjt_m(mnz_m), mvt_m(mnz_m), vv(nn*ndofn_m), ve(nn*ndofn_m))
CALL mpiglbmtx_m(mpi_ne, mpi_ele, ne, nn, ndofn_m, nnse, nquads, ngpe, ncmp, &
    & xconn, xcoord, ym, nu, ym0, nu0, rho, alpha, Tref, wqs, ns_t, ns_m, dns_t, dns_m, &
    & dnxyz, sd, tempst(2*nn+1:3*nn), stiff, Nx, Ny, Nz, &
    & kit_m(1:knz_m), kjt_m, kvt_m, knz_m, mit_m(1:mnz_m), mjt_m, mvt_m, mnz_m, ve)
CALL tricsr(kit_m, kjt_m, kvt_m, ndof_m, knz_m, MAX(knz_m,ndof_m+1))
CALL tricsr(mit_m, mjt_m, mvt_m, ndof_m, mnz_m, MAX(mnz_m,ndof_m+1))
ALLOCATE(ki_m(ndof_m+1), kj_m(knz_m), kv_m(knz_m), mi_m(ndof_m+1), mj_m(mnz_m), mv_m(mnz_m))
ki_m = kit_m(1:(ndof_m+1))
kj_m = kjt_m(1:knz_m)
kv_m = kvt_m(1:knz_m)
mi_m = mit_m(1:(ndof_m+1))
mj_m = mjt_m(1:mnz_m)
mv_m = mvt_m(1:mnz_m)

IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(1,2))
    !$ tusr(1,2) = omp_get_wtime()
    mem(1) = (knz_t*16+mnz_t*16)/32.0**4.0
    tcpu(1,3) = tcpu(1,2) - tcpu(1,1)
    tusr(1,3) = tusr(1,2) - tusr(1,1)
    WRITE(*,100) '* NNZ of spatial K matrix   ', knz_m
    WRITE(*,100) '* NNZ of spatial M matrix   ', mnz_m
    WRITE(*,105) '* Spatial assembly  (sec)   ', tcpu(1,3)
    WRITE(*,105) '* Spatial matrices   (MB)   ', mem(1)
ENDIF
!-------------------------------------------------------------------------------
! Set boundary conditions for the Temperature part
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(2,1))
    !$ tusr(2,1) = omp_get_wtime()
    !PRINT *,''
    !PRINT *, '>>> Setting boundary conditions'
ENDIF
CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
ALLOCATE(fixedbcst_t(ndofst_t), fixedbcst0_t(ndofst_t))
fixedbcst_t = 0
fixedbcst0_t = 0
DO ii = 1, nn
    IF (xcoord(ii,3) .EQ. 1._REKIND) THEN ! y = 0
        DO jj = 1, nnte
            DO kk = 1, ndofn_t
                fixedbcst_t( (ii-1)*ndofn_t+(jj-1)*ndof_t+kk ) = 1
                fixedbcst0_t( (ii-1)*ndofn_t+(jj-1)*ndof_t+kk ) = T0
            ENDDO
        ENDDO
    ENDIF
ENDDO

!Save the original K and M before applying BC
ki0 = ki_t
kj0 = kj_t
kv0 = kv_t
mi0 = mi_t
mj0 = mj_t
mv0 = mv_t

! Apply BC to K and M matrices
CALL csrapplybc(ki_t, kj_t, kv_t, ndof_t, knz_t, fixedbcst_t(1:ndof_t))
CALL csrapplybc(mi_t, mj_t, mv_t, ndof_t, mnz_t, fixedbcst_t(1:ndof_t))
ALLOCATE(rext_t(ndof_t), nodeshare_t(nn))
nodeshare_t = 0._REKIND ! number of elements that a node shared by
ntele_t = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        nodeshare_t(xconn(ii,jj)) = nodeshare_t(xconn(ii,jj)) + 1._REKIND
    ENDDO
ENDDO
DO ii = 1, nn
    IF (xcoord(ii, 3) .EQ. 2.) THEN
        ntele_t = ntele_t + nodeshare_t(ii)
    ENDIF
ENDDO
ftotal_t = q0*area
fnode_t = ftotal_t/ntele_t
rext_t = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        kk = xconn(ii,jj)
        IF (xcoord(kk, 3) .EQ. 2. .OR. xcoord(kk, 3) .EQ. 0.) THEN ! z = 130
            rext_t(kk*ndofn_t-0) = rext_t(kk*ndofn_t-0) + fnode_t
        ENDIF
    ENDDO
ENDDO
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!-------------------------------------------------------------------------------
! Set boundary conditions for the Mechanical part
!-------------------------------------------------------------------------------
ALLOCATE(fixedbcst_m(ndofst_m), fixedbcst0_m(ndofst_m))
fixedbcst_m = 0
DO ii = 1, nn
    IF (xcoord(ii,2) .EQ. 0._REKIND .OR. xcoord(ii,2) .EQ. le) THEN ! y = 0
        DO jj = 1, 2*nnte
            fixedbcst_m((ii-1)*ndofn_m+(jj-1)*ndof_m+1 : ii*ndofn_m+(jj-1)*ndof_m) &
                & = [0, 1, 0]
        ENDDO
    ENDIF

    IF (xcoord(ii,3) .EQ. 0._REKIND .AND. xcoord(ii,2) .EQ. 0._REKIND .AND. xcoord(ii,1) .EQ. 0._REKIND) THEN ! y = 0
        DO jj = 1, 2*nnte
            fixedbcst_m((ii-1)*ndofn_m+(jj-1)*ndof_m+1 : ii*ndofn_m+(jj-1)*ndof_m) &
                & = [1, 1, 1]
        ENDDO
    ELSEIF (xcoord(ii,3) .EQ. 0._REKIND .AND. xcoord(ii,2) .EQ. 0._REKIND .AND. xcoord(ii,1) .EQ. 5._REKIND) THEN ! y = 0
        DO jj = 1, 2*nnte
            fixedbcst_m((ii-1)*ndofn_m+(jj-1)*ndof_m+1 : ii*ndofn_m+(jj-1)*ndof_m) &
                & = [0, 1, 1]
        ENDDO
    ENDIF
ENDDO
fixedbcst0_m = fixedbcst_m
! Apply BC to K and M matrices
CALL csrapplybc(ki_m, kj_m, kv_m, ndof_m, knz_m, fixedbcst_m(1:ndof_m))
CALL csrapplybc(mi_m, mj_m, mv_m, ndof_m, mnz_m, fixedbcst_m(1:ndof_m))
ALLOCATE(rext_m(ndof_m), nodeshare_m(nn))
nodeshare_m = 0._REKIND ! number of elements that a node shared by
ntele_m = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        nodeshare_m(xconn(ii,jj)) = nodeshare_m(xconn(ii,jj)) + 1._REKIND
    ENDDO
ENDDO
DO ii = 1, nn
    IF (xcoord(ii, 3) .EQ. le) THEN
        ntele_m = ntele_m + nodeshare_m(ii)
    ENDIF
ENDDO
ftotal_m = p0*area
fnode_m = ftotal_m/ntele_m
rext_m = 0._REKIND
DO ii = 1, ne
    DO jj = 1, nnse
        kk = xconn(ii,jj)
        IF (xcoord(kk, 2) .EQ. le) THEN ! y = 10
            rext_m(kk*3-1) = rext_m(kk*3-1) + fnode_m
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
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!-------------------------------------------------------------------------------
! CALL MUMPS
!-------------------------------------------------------------------------------
IF (MPI_ID .EQ. ROOT) THEN
    CALL CPU_TIME(tcpu(3,1))
    !$ tusr(3,1) = omp_get_wtime()
    PRINT *,''
    PRINT *, '>>> MUMPS analysis and factorize'
ENDIF
CALL mumps_solver(mumps_par_t, ROOT, MPI_ID, MPI_P, ndof_t, ndofst_t, nnte, ki_t, kj_t, kv_t, 0)
CALL mumps_solver(mumps_par_m, ROOT, MPI_ID, MPI_P, ndof_m, ndofst_m, 2*nnte, ki_m, kj_m, kv_m, 0)
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
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
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
ALLOCATE(dispst(ndofst_m), fextst_t(ndofst_t), fextst_m(ndofst_m),  fInit(ndofst_t), fBound(ndofst_t), tcoorde(nnte))
ALLOCATE(p(mpi_ngp), D(mpi_ngp), ws(mpi_ngp), del_po(mpi_ngp), sig_D(mpi_ngp), strain_tgo(ncmp,mpi_ngp), ndfail(nn), &
        & strain_pl(ncmp,mpi_ngp), sig_eff(ncmp,mpi_ngp), sig_back(ncmp,mpi_ngp), elesta(mpi_ne), gpsta(mpi_ngp), nefail(mpi_ne))
dt          = dt0*n1        ! initial time step
nip         = nipc*n1       ! intiial nubmer of temporal interpolation points
itfail      = 0             ! > 0 if element(s) failed in current step, else 0
gpsta       = 1             ! Gauss point status, 0 = failed, 1 = active
elesta      = 1             ! Element status, 0 = failed, 1 = active
dispst      = 0._REKIND     ! Solution vector
tempst      = Tref          ! Solution vector
tcoorde     = 0._REKIND
fInit       = 0._REKIND
fBound      = 0._REKIND
p           = 0._REKIND     ! P (??) of last step
D           = 0._REKIND     ! Damage of last step
ws          = 0._REKIND     ! Stored energy of last step
del_po      = 0._REKIND     ! Delta_P of last step
sig_D       = 0._REKIND     ! Sig_D (??) of last step
strain_tgo  = 0._REKIND     ! Total global strain of last step
strain_pl   = 0._REKIND     ! Total plastic strain of last step
sig_eff     = 0._REKIND     ! Effective stress of last step
sig_back    = 0._REKIND     ! Backstress of last step

! Get S-T matrix coefficients
CALL stcoeff_t(1, dt, ck_t, cm_t, cktn_t, cmtn_t, cktn1_t, cmtn1_t)
CALL mpikronspgemv(ck_t+cktn_t, cm_t+cmtn_t, 1._REKIND, 1._REKIND, nnte, ki0, kj0, &
        & kv0, mi0, mj0, mv0, ndof_t, tempst, fInit)
tempst = fixedbcst0_t
CALL mpikronspgemv(-ck_t-cktn_t, -cm_t-cmtn_t, 1._REKIND, 1._REKIND, nnte, ki0, kj0, &
        & kv0, mi0, mj0, mv0, ndof_t, tempst, fBound)
DO ii = 1, nstep ! Loop over space time slabs 
CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    ! Temporal nodal coordinates of iith S-T slab
    DO jj = 1, nnte
        tcoorde(jj) = tcoorde(nnte) + (jj-1)*dt/(nnte-1)
    ENDDO
    !--------------------
    ! Check damage first
    !--------------------
    IF (MPI_ID .EQ. ROOT) THEN
        CALL CPU_TIME(tcpu(6,1))
        !$ tusr(6,1) = omp_get_wtime()
        WRITE(iout2,*) 'Check damage ... start, step = ', ii
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
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
        !EXIT
        ! Delete failed elements
        ! 1st step, delete the failed elements
        CALL delelement(nn, mpi_ne, ndofn_m, nnse, ngpe, ncmp, nquads, &
            & xconn(mpi_ele,:), xcoord, dns_m, dnxyz, sd, stiff, wqs, nef, &
            & nefail, ki_m, kj_m, kv_m, ndof_m, knz_m)
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
        CALL adddiagonal(nnse, ndofn_m, ndf, ndfail, ki_m, kj_m, kv_m, ndof_m, knz_m)
        ! 3rd step, constrain the disconnected nodes
        IF (ndf.GT.0) THEN
            DO i = 1, ndf
                j = ndfail(i)
                DO k = 1, 2*nnte
                    fixedbcst0_m((j-1)*ndofn_m+(k-1)*ndof_m+1 : j*ndofn_m+(k-1)*ndof_m) = [1, 1, 1]
                ENDDO
            ENDDO
        ENDIF
        ! update bc
        CALL MPI_ALLREDUCE(fixedbcst0_m, fixedbcst_m, ndofst_m, MPI_INTEGER, &
            & MPI_SUM, MPI_COMM_WORLD, MPI_IERR)
        IF (MPI_IERR .NE. 0) THEN
            PRINT *, 'MPI_ALLREDUCE error = ',  MPI_IERR
            GOTO 1500
        ENDIF
        fixedbcst0_m = fixedbcst_m
        CALL csrapplybc(mi_m, mj_m, mv_m, ndof_m, mnz_m, fixedbcst_m(1:ndof_m))
        CALL csrapplybc(ki_m, kj_m, kv_m, ndof_m, knz_m, fixedbcst_m(1:ndof_m))
        ! 4th step, generate new preconditioner by MUMPS
        CALL mumps_solver(mumps_par_m, ROOT, MPI_ID, MPI_P, ndof_m, ndofst_m, 2*nnte, ki_m, kj_m, kv_m, 1)
        ! Reduce time step
        dt = dt0*n2
        nip = nipc*n2
    ENDIF
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Check damage ... done, step = ', ii
        CALL CPU_TIME(tcpu(6,2))
        !$ tusr(6,2) = omp_get_wtime()
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    !--------------------------
    !   SOLVE TEMPERATURE PART
    !--------------------------
    ! Calculate temporal component of external force vectors
    IF (ii .EQ. 1) THEN
       fextst_t = fInit
    ELSE
       fextst_t = 0._REKIND
    ENDIF
    CALL getfextc_t(ii, dt, fextc_t)
    kk = 0
    ll = 1 - ndof_t
    DO jj = 1, nnte
        kk = kk + ndof_t
        ll = ll + ndof_t
        ! Space time external force analagous vector
        fextst_t(ll:kk) = fextst_t(ll:kk) + fextc_t(jj)*rext_t
    ENDDO 
    ! Use KRONSPGEMV perform spgemv
    CALL stcoeff_t(ii, dt, ck_t, cm_t, cktn_t, cmtn_t, cktn1_t, cmtn1_t)
    CALL mpikronspgemv(cktn1_t, cmtn1_t, 1._REKIND, 1._REKIND, nnte, ki_t, kj_t, &
        & kv_t, mi_t, mj_t, mv_t, ndof_t, tempst, fextst_t)
    fextst_t = fextst_t + fBound
    DO jj = 1, ndofst_t
        IF (fixedbcst_t(jj) .NE. 0) fextst_t(jj) = 0._REKIND
        !PRINT *, 'ID = ', jj, 'fextst = ', fextst(jj)
    ENDDO
    !PRINT *, 'ID = ', MPI_ID, 'sum(fextst) = ', sum(dabs(fextst_t))
    tempst = 0._REKIND
    CALL mpipgmres(ck_t+cktn_t, cm_t+cmtn_t, nnte, ki_t, kj_t, kv_t, mi_t, mj_t, mv_t, ndof_t, &
        &  mumps_par_t, fextst_t, tempst, im_t, eps_t, maxits_t, iout, ierr)
    DO jj = 1, ndofst_t
        IF (fixedbcst_t(jj) .NE. 0) tempst(jj) = fixedbcst0_t(jj)
        !PRINT *, 'ID = ', jj, 'disp = ', dispst(jj)
    ENDDO    
    IF (ierr .NE. 0) THEN
        PRINT *, 'Error: ', ierr
        GOTO 1000
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_Barrier error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    
    CALL MPI_BCAST(tempst, ndofst_t, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_BCAST error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    !-------------------------------------------------------------------
    ! Update the Material property
    !-------------------------------------------------------------------
    knz_m = mpi_ne*(nnse*ndofn_m)**2
    mnz_m = knz_m
    vv = 0._REKIND
    DEALLOCATE(ki_m, kj_m, kv_m, mi_m, mj_m, mv_m)
    DEALLOCATE(kit_m, kjt_m, kvt_m, mit_m, mjt_m, mvt_m)
    ALLOCATE(kit_m(MAX(knz_m,ndof_m+1)), kjt_m(knz_m), kvt_m(knz_m), mit_m(MAX(mnz_m,ndof_m+1)), &
            & mjt_m(mnz_m), mvt_m(mnz_m))
    ALLOCATE(ki_m(ndof_m+1), kj_m(knz_m), kv_m(knz_m), mi_m(ndof_m+1), mj_m(mnz_m), mv_m(mnz_m))
    CALL material_update(nnse, nquads, ndofst_t, nn, ne, ngpe, mpi_ne, mpi_ele, xconn, &
                & ns_t, tempst, Tref, ym0, nu0, ym, nu)
    CALL mpiglbmtx_m(mpi_ne, mpi_ele, ne, nn, ndofn_m, nnse, nquads, ngpe, ncmp, &
           & xconn, xcoord, ym, nu, ym0, nu0, rho, alpha, Tref, wqs, ns_t, ns_m, dns_t, dns_m, &
           & dnxyz, sd, tempst(2*nn+1:3*nn), stiff, Nx, Ny, Nz, &
           & kit_m(1:knz_m), kjt_m, kvt_m, knz_m, mit_m(1:mnz_m), mjt_m, mvt_m, mnz_m, ve)
    CALL MPI_ALLREDUCE(ve, vv, nn*ndofn_m, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_ALLREDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL tricsr(kit_m, kjt_m, kvt_m, ndof_m, knz_m, MAX(knz_m,ndof_m+1))
    CALL tricsr(mit_m, mjt_m, mvt_m, ndof_m, mnz_m, MAX(mnz_m,ndof_m+1))

    ki_m = kit_m(1:(ndof_m+1))
    kj_m = kjt_m(1:knz_m)
    kv_m = kvt_m(1:knz_m)
    mi_m = mit_m(1:(ndof_m+1))
    mj_m = mjt_m(1:mnz_m)
    mv_m = mvt_m(1:mnz_m)
    ! Apply BC to K and M matrices
    CALL csrapplybc(ki_m, kj_m, kv_m, ndof_m, knz_m, fixedbcst_m(1:ndof_m))
    CALL csrapplybc(mi_m, mj_m, mv_m, ndof_m, mnz_m, fixedbcst_m(1:ndof_m))
    CALL mumps_solver(mumps_par_m, ROOT, MPI_ID, MPI_P, ndof_m, ndofst_m, 2*nnte, ki_m, kj_m, kv_m, 1)
    !----------------------------
    !   SOLVE THE MECHANICAL PART
    !----------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    Tmin = 0._REKIND
    CALL MPI_REDUCE(MINVAL(tempst(2*nn+1:3*nn)), Tmin, 1, MPI_REAL8, MPI_MIN, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    Tmax = 0._REKIND
    CALL MPI_REDUCE(MAXVAL(tempst(2*nn+1:3*nn)), Tmax, 1, MPI_REAL8, MPI_MAX, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    ! Calculate temporal component of external force vectors
    fextc_m = 0._REKIND
    fextst_m = 0._REKIND
    CALL getfextc_m(ii, angvel_m, dt, fextc_m, fint_m)
    kk = 0
    ll = 1 - ndof_m
    DO jj = 1, 2*nnte
        kk = kk + ndof_m
        ll = ll + ndof_m
        ! Space time external force analagous vector
        fextst_m(ll:kk) = fextc_m(jj)*rext_m + fint_m(jj)*vv
    ENDDO 
    ! Get S-T matrix coefficients
    CALL stcoeff_m(ii, dt, angvel_m, ck_m, cm_m, cktn_m, cmtn_m, cktn1_m, cmtn1_m)
    ! Use KRONSPGEMV perform spgemv
    CALL mpikronspgemv(cktn1_m, cmtn1_m, 1._REKIND, 1._REKIND, 2*nnte, ki_m, kj_m, &
        & kv_m, mi_m, mj_m, mv_m, ndof_m, dispst, fextst_m)
    !PRINT *, 'ID = ', MPI_ID, 'sum(fextst) = ', sum(dabs(fextst))
    DO jj = 1, ndofst_m
        IF (fixedbcst_m(jj) .NE. 0) fextst_m(jj) = 0._REKIND 
    ENDDO
    dispst = 0._REKIND
    IF (MPI_ID .EQ. ROOT) THEN
        WRITE(iout2,*) 'Solve linear system ... start, step = ', ii
        CALL CPU_TIME(tcpu(4,1))
        !$ tusr(4,1) = omp_get_wtime()
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    CALL mpipgmres(ck_m+cktn_m, cm_m+cmtn_m, 2*nnte, ki_m, kj_m, kv_m, mi_m, mj_m, mv_m, ndof_m, &
        &  mumps_par_m, fextst_m, dispst, im_m, eps_m, maxits_m, iout, ierr)
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
        DO jj = 1, ndofst_m
            IF (fixedbcst_m(jj) .NE. 0) dispst(jj) = 0._REKIND 
        ENDDO
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_Barrier error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    
    CALL MPI_BCAST(dispst, ndofst_m, MPI_REAL8, ROOT, MPI_COMM_WORLD, MPI_IERR)
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
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    Umin = 0._REKIND
    CALL MPI_REDUCE(MINVAL(dispst(nn*2+1:nn*3)), Umin, 1, MPI_REAL8, MPI_MIN, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    Umax = 0._REKIND
    CALL MPI_REDUCE(MAXVAL(dispst(nn*2+1:nn*3)), Umax, 1, MPI_REAL8, MPI_MAX, ROOT, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD, MPI_IERR)
    CALL kokkos_damage(MPI_ID, job, mpi_ne, nn, ndofn_m, nnse, nnte, ngpe, &
            & ncmp, nip, xconn(mpi_ele,:), Nx, Ny, Nz, ym, nu, angvel_m, &
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
    CALL MPI_ALLREDUCE(MAXVAL(D), dmax, 1, MPI_REAL8, MPI_MAX, &
        & MPI_COMM_WORLD, MPI_IERR)
    IF (MPI_IERR .NE. 0) THEN
        PRINT *, 'MPI_REDUCE error = ',  MPI_IERR
        GOTO 1500
    ENDIF
    CALL MPI_ALLREDUCE(MAXVAL(ws), wsmax, 1, MPI_REAL8, MPI_MAX, &
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
        CALL mpipostprocess(job, mpi_ne, nn, ndofn_m, nnse, nnte, ngpe, ncmp, &
            & nip, nipc, xconn(mpi_ele,:), Nx, Ny, Nz, ym, nu, angvel_m, &
            & tcoorde, dispst, tempst, gpsta, elesta, D, ws, MPI_ID)
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
                & 'POSTPR', 'mas(Ws)', 'max(D)', 'min(T)', 'max(T)'
        ENDIF
        WRITE(*,115) ii, tcoorde(nnte), tusr(4:7,3), wsmax, dmax, Tmin, Tmax
    ENDIF
ENDDO
! Clear memory and exit
DEALLOCATE(xcoord, xconn, tempst, fextst_t, fextst_m, fBound, fInit, tcoorde, dns_t, dns_m, wqs, vv)
DEALLOCATE(fixedbcst_t, fixedbcst0_t, rext_t)
DEALLOCATE(ki_t, kj_t, kv_t, mi_t, mj_t, mv_t)
DEALLOCATE(ki_m, kj_m, kv_m, mi_m, mj_m, mv_m)
DEALLOCATE(ns_t, kit_t, kjt_t, kvt_t, mit_t, mjt_t, mvt_t)
DEALLOCATE(ns_m, kit_m, kjt_m, kvt_m, mit_m, mjt_m, mvt_m)
DEALLOCATE(ki0, kj0, kv0, mi0, mj0, mv0)
DEALLOCATE(mpi_ele, mpi_gp)
DEALLOCATE(p, D, ws, del_po, sig_D, strain_tgo, strain_pl, sig_eff, sig_back, &
    & elesta, gpsta, necnt, necnt0, nefail, ndfail)
IF (MPI_P .EQ. 1) THEN
    DEALLOCATE(mumps_par_t%IRN, mumps_par_t%JCN, mumps_par_t%A)
    DEALLOCATE(mumps_par_m%IRN, mumps_par_m%JCN, mumps_par_m%A)
ELSE
    DEALLOCATE(mumps_par_t%IRN_loc, mumps_par_t%JCN_loc, mumps_par_t%A_loc)
    DEALLOCATE(mumps_par_m%IRN_loc, mumps_par_m%JCN_loc, mumps_par_m%A_loc)
ENDIF
IF (MPI_ID .EQ. ROOT) THEN
    DEALLOCATE(mumps_par_t%RHS)
    DEALLOCATE(mumps_par_m%RHS)
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
mumps_par_t%JOB = -2
mumps_par_m%JOB = -2
!CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
CALL DMUMPS(mumps_par_t)
CALL DMUMPS(mumps_par_m)
IF (mumps_par_t%INFOG(1).LT.0) THEN
    WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
        &   "  mumps_par%INFOG(1)= ", mumps_par_t%INFOG(1), &
        &   "  mumps_par%INFOG(2)= ", mumps_par_t%INFOG(2) 
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
100   FORMAT(A35, 9X, I16)
105   FORMAT(A35, 9X, F16.2)
110   FORMAT(A8,2X,9(A8, 2X))
115   FORMAT(I8,2X,F8.1,2X,6(F8.2,2X),2(F10.2))
END PROGRAM main3d
