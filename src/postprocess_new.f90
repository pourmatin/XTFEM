! Timestamp: Thu Feb 18 20:03:05 CST 2016 
!_______________________________________________________________________________
!
! Post-processing codes for
!   'Extended space-time FEM for elastodynamics in three-dimensions'
!
! Developed by Rui Zhang @ UT Dallas, Richardson, TX, Oct, 2015
! E-mail: rui.zhang4@utdallas.edu (USA), ruizhang@mail.nwpu.edu.cn (China)
! Supervised by Prof. Dong Qian (UTD, USA) and Prof. Lihua Wen (NPU, China)
!_______________________________________________________________________________
SUBROUTINE postprocess(jobname, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, &
    & nipc, xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, gpsta, elsta, D, ws)
USE kinds
!! Variable declaration
IMPLICIT NONE
INCLUDE 'mpif.h'
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, nipc
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: ym, nu, w
REAL(KIND=REKIND), DIMENSION(nnte), INTENT(IN) :: tcoorde
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
INTEGER, DIMENSION(ngpe*ne), INTENT(IN) :: gpsta
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D, ws
! ---Internal variables---
CHARACTER(LEN=80) :: odb
INTEGER :: ngp, noel, npt, nout
REAL(KIND=REKIND) :: t, dt, coef
REAL(KIND=REKIND), DIMENSION(ncmp, ncmp) :: Dmat, Cmat
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nt, strain, stress, disp
REAL(KIND=REKIND), DIMENSION(nn, ndofn) :: dispnd
INTEGER :: i, j, ii, jj
!--- Set parameters & allocate memory ------------------------------------------
ngp = ngpe*ne ! number of total Gauss points
odb = TRIM(jobname)//'.dat'
ALLOCATE(Nt(2*nnte, nip), disp(nnse*ndofn, ngp), strain(ngp, ncmp), &
    & stress(ngp, ncmp))
!--- Calculate temporal shape function -----------------------------------------
dt = tcoorde(nnte) - tcoorde(1)
!$omp parallel num_threads(24)
!$omp do private(t)
DO i = 1, nip
    t = tcoorde(1) + (i-1)*dt/nip
    Nt(1,i) = 2._REKIND*(tcoorde(2)-t)*(tcoorde(3)-t)/dt**2
    Nt(2,i) = -4._REKIND*(tcoorde(1)-t)*(tcoorde(3)-t)/dt**2
    Nt(3,i) = 2._REKIND*(tcoorde(1)-t)*(tcoorde(2)-t)/dt**2
    Nt(4,i) = Nt(1,i)*(DSIN(w*t)-DSIN(w*tcoorde(1)))
    Nt(5,i) = Nt(2,i)*(DSIN(w*t)-DSIN(w*tcoorde(2)))
    Nt(6,i) = Nt(3,i)*(DSIN(w*t)-DSIN(w*tcoorde(3)))
ENDDO
!$omp end do
!$omp end parallel
!--- Calculate constitutive matrices (stiffness and compliance) ----------------
coef = ym/((1._REKIND+nu)*(1._REKIND-2._REKIND*nu))
Dmat = 0._REKIND
Cmat = 0._REKIND
DO i = 1, 3
    Dmat(i,i) = coef*(1._REKIND-nu)
    Dmat(i+3,i+3) = coef*(1._REKIND-2._REKIND*nu)/2._REKIND
    Cmat(i,i) = 1._REKIND/ym
    Cmat(i+3,i+3) = 2._REKIND*(1._REKIND+nu)/ym
ENDDO
DO i = 1, 3
    DO j = 1, 3
        IF (i.NE.j) THEN
            Dmat(i,j) = coef*nu
            Cmat(i,j) = -nu/ym
        ENDIF
    ENDDO
ENDDO
nout = nip - nipc/4*3 + 1
!--- Evaluate strain-stress ----------------------------------------------------
!$omp parallel num_threads(24)
!$omp do private(noel, npt)
DO ii = 1, ngp ! Loop over spatial integration points
    noel = (ii-1)/ngpe + 1   ! Element index
    npt = MOD(ii, ngpe)     ! Gauss point index
    IF (npt .EQ. 0) npt = ngpe
    ! Interpolate displacement
    CALL dispinterp(nnte, nnse, nn, ndofn, xconn(noel,:), sol, Nt(:,nout), &
        & disp(:, ii))
    ! Calculate strain at mesoscale
    CALL calstrain(nnse, ndofn, ncmp, disp(:,ii), Nx(:,ii), Ny(:,ii), &
            & Nz(:,ii), Dmat, strain(ii, :), stress(ii,:))
ENDDO
!$omp end do
!$omp end parallel
!--- Clean up ------------------------------------------------------------------
DO i = 1, ngp
    IF (gpsta(i).EQ.0) THEN ! Failed integration points
        strain(i,:) = 0._REKIND
        stress(i,:) = 0._REKIND
    ENDIF
ENDDO
!--- Write Odb -----------------------------------------------------------------
t = tcoorde(1) + (nout-1)*dt/nip
CALL dispinterpt(nn, ndofn, nnte, Nt(:, nout), dispnd, sol)
CALL writeodb(odb, ne, nn, ndofn, nnse, ngpe, ngp, ncmp, t, dispnd, strain, &
    & stress, D, Ws, elsta)
DEALLOCATE(Nt, disp, strain, stress)
RETURN
END SUBROUTINE postprocess


SUBROUTINE writeodb(filename, ne, nn, ndofn, nnse, ngpe, ngp, ncmp, time, &
    & disp, strain, stress, D, Ws, elsta)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, ngpe, ngp, ncmp
REAL(KIND=REKIND), INTENT(IN) :: time
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: disp
REAL(KIND=REKIND), DIMENSION(ne*ngpe, ncmp), INTENT(IN) :: strain, stress
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D, Ws
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
! ---Internal variables---
INTEGER :: fid, ii, jj
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg
fid = 10
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(filename)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
WRITE(fid, *) time ! time
DO ii = 1, nn ! displacement
    WRITE(fid, *) (disp(ii, jj), jj = 1, ndofn)
ENDDO
DO ii = 1, ne*ngpe ! strain
    WRITE(fid, *) (strain(ii, jj), jj = 1, ncmp)
ENDDO
DO ii = 1, ne*ngpe ! stress
    WRITE(fid, *) (stress(ii, jj), jj = 1, ncmp)
ENDDO
DO ii = 1, ne*ngpe ! damage
    WRITE(fid, *) D(ii)
ENDDO
DO ii = 1, ne*ngpe ! energy
    WRITE(fid, *) Ws(ii)
ENDDO
DO ii = 1, ne ! status
    WRITE(fid, *) elsta(ii)
ENDDO
CLOSE(fid)
RETURN
END SUBROUTINE writeodb


SUBROUTINE dispinterpt(nn, ndofn, nnte, Nt, disp, sol)
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ndofn, nnte
REAL(KIND=REKIND), DIMENSION(2*nnte), INTENT(IN) :: Nt
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(OUT) :: disp
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
! ---Internal variables---
INTEGER :: i, j, k, idx
disp = 0._REKIND
DO i = 1, 2*nnte
    DO j = 1, nn
        DO k = 1, ndofn
            idx = (i-1)*nn*ndofn + (j-1)*ndofn + k
            disp(j,k) = disp(j,k) + Nt(i)*sol(idx)
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE dispinterpt

SUBROUTINE dispinterp(nnte, nnse, nn, ndofn, nodes, sol, Nt, disp)
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nnte, nnse, nn, ndofn
INTEGER, DIMENSION(nnse), INTENT(IN) :: nodes
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
REAL(KIND=REKIND), DIMENSION(2*nnte), INTENT(IN) :: Nt
REAL(KIND=REKIND), DIMENSION(nnse*ndofn), INTENT(OUT) :: disp
! ---Internal variables---
INTEGER :: i, j, k, ndid, dofid
disp = 0._REKIND
DO i = 1, 2*nnte
    DO j = 1, nnse
        ndid = nodes(j)
        DO k = 1, ndofn
            dofid = (i-1)*nn*ndofn + (ndid-1)*ndofn + k
            disp((j-1)*ndofn+k) = disp((j-1)*ndofn+k) + Nt(i)*sol(dofid)
        ENDDO
    ENDDO
ENDDO
RETURN
END SUBROUTINE dispinterp


SUBROUTINE calstrain(nnse, ndofn, ncmp, disp, Nx, Ny, Nz, dmat, strain, stress)
! Calculate the strain at mesoscale for a given Gauss point
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nnse, ndofn, ncmp
REAL(KIND=REKIND), DIMENSION(nnse*ndofn), INTENT(IN) ::  disp
REAL(KIND=REKIND), DIMENSION(nnse), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), DIMENSION(ncmp,ncmp), INTENT(IN) :: dmat
REAL(KIND=REKIND), DIMENSION(ncmp), INTENT(OUT) :: strain, stress
! ---Internal variables---
INTEGER :: i, j
strain = 0._REKIND
stress = 0._REKIND
DO i = 1, nnse
    j = (i-1)*ndofn
    strain(1) = strain(1) + Nx(i)*disp(j+1)
    strain(2) = strain(2) + Ny(i)*disp(j+2)
    strain(3) = strain(3) + Nz(i)*disp(j+3)
    strain(4) = strain(4) + Ny(i)*disp(j+1) + Nx(i)*disp(j+2)
    strain(5) = strain(5) + Nz(i)*disp(j+2) + Ny(i)*disp(j+3)
    strain(6) = strain(6) + Nz(i)*disp(j+1) + Nx(i)*disp(j+3)
ENDDO
stress = MATMUL(dmat, strain)
RETURN
END SUBROUTINE calstrain

SUBROUTINE mpipostprocess(jobname, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, &
    & nip, nipc, xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, gpsta, elsta, &
    & D, ws, stress, strain, MPI_ID)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, nipc
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: ym, nu, w
REAL(KIND=REKIND), DIMENSION(nnte), INTENT(IN) :: tcoorde
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
INTEGER, DIMENSION(ngpe*ne), INTENT(IN) :: gpsta
INTEGER, DIMENSION(ncmp, ngpe*ne), INTENT(INOUT) :: stress, strain
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D, ws
INTEGER, INTENT(IN) :: MPI_ID
! ---Internal variables---
CHARACTER(LEN=80) :: odb, temp, mpi_odb
INTEGER :: ngp, noel, npt, nout, fid
REAL(KIND=REKIND) :: t, dt, coef
REAL(KIND=REKIND), DIMENSION(ncmp, ncmp) :: Dmat, Cmat
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nt, strainT, stressT, disp
REAL(KIND=REKIND), DIMENSION(nn, ndofn) :: dispnd
INTEGER :: i, j, ii, jj
!--- Set parameters & allocate memory ------------------------------------------
ngp = ngpe*ne ! number of total Gauss points
odb = TRIM(jobname)//'.dat'
WRITE(temp,'(i8)') MPI_ID
fid = 700 + MPI_ID
mpi_odb = TRIM(jobname)//'_'//TRIM(ADJUSTL(temp))
ALLOCATE(Nt(2*nnte, nip), disp(nnse*ndofn, ngp), strainT(ngp, ncmp), &
    & stressT(ngp, ncmp))
!--- Calculate temporal shape function -----------------------------------------
dt = tcoorde(nnte) - tcoorde(1)
!!$omp parallel num_threads(24)
!!$omp do private(t)
DO i = 1, nip
    t = tcoorde(1) + (i-1)*dt/nip
    Nt(1,i) = 2._REKIND*(tcoorde(2)-t)*(tcoorde(3)-t)/dt**2
    Nt(2,i) = -4._REKIND*(tcoorde(1)-t)*(tcoorde(3)-t)/dt**2
    Nt(3,i) = 2._REKIND*(tcoorde(1)-t)*(tcoorde(2)-t)/dt**2
    Nt(4,i) = Nt(1,i)*(DSIN(w*t)-DSIN(w*tcoorde(1)))
    Nt(5,i) = Nt(2,i)*(DSIN(w*t)-DSIN(w*tcoorde(2)))
    Nt(6,i) = Nt(3,i)*(DSIN(w*t)-DSIN(w*tcoorde(3)))
ENDDO
!!$omp end do
!!$omp end parallel
!--- Calculate constitutive matrices (stiffness and compliance) ----------------
coef = ym/((1._REKIND+nu)*(1._REKIND-2._REKIND*nu))
Dmat = 0._REKIND
Cmat = 0._REKIND
DO i = 1, 3
    Dmat(i,i) = coef*(1._REKIND-nu)
    Dmat(i+3,i+3) = coef*(1._REKIND-2._REKIND*nu)/2._REKIND
    Cmat(i,i) = 1._REKIND/ym
    Cmat(i+3,i+3) = 2._REKIND*(1._REKIND+nu)/ym
ENDDO
DO i = 1, 3
    DO j = 1, 3
        IF (i.NE.j) THEN
            Dmat(i,j) = coef*nu
            Cmat(i,j) = -nu/ym
        ENDIF
    ENDDO
ENDDO
nout = nip - nipc/4*3 + 1
!--- Evaluate strain-stress ----------------------------------------------------
!!$omp parallel num_threads(24)
!!$omp do private(noel, npt)
DO ii = 1, ngp ! Loop over spatial integration points
    noel = (ii-1)/ngpe + 1   ! Element index
    npt = MOD(ii, ngpe)     ! Gauss point index
    IF (npt .EQ. 0) npt = ngpe
    ! Interpolate displacement
    CALL dispinterp(nnte, nnse, nn, ndofn, xconn(noel,:), sol, Nt(:,nout), &
        & disp(:, ii))
    ! Calculate strain at mesoscale
    CALL calstrain(nnse, ndofn, ncmp, disp(:,ii), Nx(:,ii), Ny(:,ii), &
            & Nz(:,ii), Dmat, strainT(ii, :), stressT(ii,:))
ENDDO
!!$omp end do
!!$omp end parallel
!--- Clean up ------------------------------------------------------------------
DO i = 1, ngp
    IF (gpsta(i).EQ.0) THEN ! Failed integration points
        strainT(i,:) = 0._REKIND
        stressT(i,:) = 0._REKIND
    ENDIF
ENDDO
stress = TRANSPOSE(stressT)
strain = Transpose(strainT)
!--- Write Odb -----------------------------------------------------------------
t = tcoorde(1) + (nout-1)*dt/nip
IF (MPI_ID.EQ.0) THEN
    CALL dispinterpt(nn, ndofn, nnte, Nt(:, nout), dispnd, sol)
    CALL writetime(jobname, t)
    CALL writedisp(jobname, nn, ndofn, dispnd)
ENDIF
!CALL writeother(MPI_ID, fid, mpi_odb, ne, ngpe, ncmp, strain, stress, D, Ws, elsta)
DEALLOCATE(Nt, disp, stressT, strainT)
RETURN
END SUBROUTINE mpipostprocess


SUBROUTINE writetime(filename, time)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
REAL(KIND=REKIND), INTENT(IN) :: time
! ---Internal variables---
INTEGER :: fid
CHARACTER(LEN=80) :: temp
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg
temp = TRIM(filename)//'.time'
fid = 100
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
WRITE(fid, *) time ! time
CLOSE(fid)
RETURN
END SUBROUTINE writetime


SUBROUTINE writedisp(filename, nn, ndofn, disp)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: nn, ndofn
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: disp
! ---Internal variables---
INTEGER :: fid, ii, jj
CHARACTER(LEN=80) :: temp
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg
temp = TRIM(filename)//'.disp'
fid = 100
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, nn ! displacement
    WRITE(fid, *) (disp(ii, jj), jj = 1, ndofn)
ENDDO
CLOSE(fid)
RETURN
END SUBROUTINE writedisp


SUBROUTINE writeother(MPI_ID, fid, filename, ne, ngpe, ncmp, strain, stress, D, Ws, elsta)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: ne, ngpe, ncmp, fid, MPI_ID
REAL(KIND=REKIND), DIMENSION(ne*ngpe, ncmp), INTENT(IN) :: strain, stress
REAL(KIND=REKIND), DIMENSION(ne*ngpe), INTENT(IN) :: D, Ws
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
! ---Internal variables---
INTEGER :: ii, jj
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg
!integer(kind=MPI_OFFSET_KIND) disp
INTEGER :: real_size, thefile,disp, ierr
real_size = SIZEOF(D(1))

!strainT = TRANSPOSE(strain)
!stressT = TRANSPOSE(stress)

disp = real_size * ne * ngpe * ncmp
!temp = TRIM(filename)//'.strain'
!CALL MPI_File_open(MPI_COMM_WORLD, 'SENT.strain', MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
!                       MPI_INFO_NULL, thefile, ierr)
!CALL MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, & 
!                           MPI_INTEGER, 'native', & 
!                           MPI_INFO_NULL, ierr) 
!CALL MPI_FILE_WRITE(thefile, strainT, ne * ngpe * ncmp, MPI_INTEGER, & 
!                        MPI_STATUS_IGNORE, ierr) 
!CALL MPI_FILE_CLOSE(thefile, ierr)
!OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
!    & IOMSG=errmesg)
!IF (errstat /= 0) THEN
!    WRITE(*,*) 'Error code: ', errstat
!    WRITE(*,*) 'Error message: ', errmesg
!    STOP
!ENDIF
!DO ii = 1, ne*ngpe !strain
!    WRITE(fid, '(6(E16.8,2X))') (strain(ii, jj), jj = 1, ncmp)
!ENDDO
!CLOSE(fid)

!CALL MPI_File_open(MPI_COMM_WORLD, 'SENT.stress', MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!                       MPI_INFO_NULL, thefile, ierr)
!CALL MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, &
!                           MPI_INTEGER, 'native', &
!                           MPI_INFO_NULL, ierr)
!CALL MPI_FILE_WRITE(thefile, stressT, ne * ngpe * ncmp, MPI_INTEGER, &
!                        MPI_STATUS_IGNORE, ierr)
!CALL MPI_FILE_CLOSE(thefile, ierr)

!temp = TRIM(filename)//'.stress'
!OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
!    & IOMSG=errmesg)
!IF (errstat /= 0) THEN
!    WRITE(*,*) 'Error code: ', errstat
!    WRITE(*,*) 'Error message: ', errmesg
!    STOP
!ENDIF
!DO ii = 1, ne*ngpe ! stress
!    WRITE(fid, '(6(E16.8,2X))') (stress(ii, jj), jj = 1, ncmp)
!ENDDO
!CLOSE(fid)

disp = real_size * ne * ngpe
!CALL MPI_File_open(MPI_COMM_WORLD, 'SENT.damage', MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!                       MPI_INFO_NULL, thefile, ierr)
!CALL MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, &
!                           MPI_INTEGER, 'native', &
!                           MPI_INFO_NULL, ierr)
!CALL MPI_FILE_WRITE(thefile, D, ne * ngpe, MPI_INTEGER, &
!                        MPI_STATUS_IGNORE, ierr)
!CALL MPI_FILE_CLOSE(thefile, ierr)
!temp = TRIM(filename)//'.damage'
!OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
!    & IOMSG=errmesg)
!IF (errstat /= 0) THEN
!    WRITE(*,*) 'Error code: ', errstat
!    WRITE(*,*) 'Error message: ', errmesg
!    STOP
!ENDIF
!DO ii = 1, ne*ngpe ! damage
!    WRITE(fid, *) D(ii)
!ENDDO
!CLOSE(fid)

!CALL MPI_File_open(MPI_COMM_WORLD, 'SENT.energy', MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!                       MPI_INFO_NULL, thefile, ierr)
!CALL MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, &
!                           MPI_INTEGER, 'native', &
!                           MPI_INFO_NULL, ierr)
!CALL MPI_FILE_WRITE(thefile, Ws, ne * ngpe, MPI_INTEGER, &
!                        MPI_STATUS_IGNORE, ierr)
!CALL MPI_FILE_CLOSE(thefile, ierr)
!temp = TRIM(filename)//'.energy'
!OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
!    & IOMSG=errmesg)
!IF (errstat /= 0) THEN
!    WRITE(*,*) 'Error code: ', errstat
!    WRITE(*,*) 'Error message: ', errmesg
!    STOP
!ENDIF
!DO ii = 1, ne*ngpe ! energy
!    WRITE(fid, *) Ws(ii)
!ENDDO
!CLOSE(fid)

disp = real_size * ne
!CALL MPI_File_open(MPI_COMM_WORLD, 'SENT.status', MPI_MODE_WRONLY + MPI_MODE_CREATE, &
!                       MPI_INFO_NULL, thefile, ierr)
!CALL MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, &
!                           MPI_INTEGER, 'native', &
!                           MPI_INFO_NULL, ierr)
!CALL MPI_FILE_WRITE(thefile, elsta, ne, MPI_INTEGER, &
!                        MPI_STATUS_IGNORE, ierr)
!CALL MPI_FILE_CLOSE(thefile, ierr)

!temp = TRIM(filename)//'.status'
!OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
!    & IOMSG=errmesg)
!IF (errstat /= 0) THEN
!    WRITE(*,*) 'Error code: ', errstat
!    WRITE(*,*) 'Error message: ', errmesg
!    STOP
!ENDIF
!DO ii = 1, ne ! status
!    WRITE(fid, *) elsta(ii)
!ENDDO
!CLOSE(fid)
RETURN
END SUBROUTINE writeother


SUBROUTINE initfile(filename, MPI_P)
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: MPI_P
! ---Internal variables---
INTEGER :: fid, i, j
CHARACTER(LEN=80) :: temp, mpifn0, mpifn1
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg
fid = 100
temp = TRIM(filename)//'.time'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), STATUS='unknown', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
CLOSE(fid)
temp = TRIM(filename)//'.disp'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), STATUS='unknown', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
CLOSE(fid)
RETURN
END SUBROUTINE initfile
