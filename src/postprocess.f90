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
SUBROUTINE tempinterpt(nn, ndofn, nnte, Nt, temp, sol)
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ndofn, nnte
REAL(KIND=REKIND), DIMENSION(2*nnte), INTENT(IN) :: Nt
REAL(KIND=REKIND), DIMENSION(nn), INTENT(OUT) :: temp
REAL(KIND=REKIND), DIMENSION(nnte*nn*ndofn), INTENT(IN) :: sol
! ---Internal variables---
INTEGER :: i, j, k, idx
temp = 0._REKIND
DO i = 1, nnte
    DO j = 1, nn
       idx = (i-1)*nn + j
       temp(j) = temp(j) + Nt(i)*sol(idx)
    ENDDO
ENDDO
RETURN
END SUBROUTINE tempinterpt


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
    & nip, nipc, xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, tempst, gpsta, elsta, &
    & D, ws, MPI_ID)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, nipc
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: w
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: ym, nu
REAL(KIND=REKIND), DIMENSION(nnte), INTENT(IN) :: tcoorde
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
REAL(KIND=REKIND), DIMENSION(nnte*nn), INTENT(IN) :: tempst
INTEGER, DIMENSION(ngpe*ne), INTENT(IN) :: gpsta
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D, ws
INTEGER, INTENT(IN) :: MPI_ID
! ---Internal variables---
CHARACTER(LEN=80) :: odb, temp, mpi_odb
INTEGER :: ngp, noel, npt, nout, fid
REAL(KIND=REKIND) :: t, dt, coef
REAL(KIND=REKIND), DIMENSION(ncmp, ncmp) :: Dmat, Cmat
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nt, strain, stress, disp
REAL(KIND=REKIND), DIMENSION(nn, ndofn) :: dispnd
REAL(KIND=REKIND), DIMENSION(nn) :: tempnd
INTEGER :: i, j, ii, jj
!--- Set parameters & allocate memory ------------------------------------------
ngp = ngpe*ne ! number of total Gauss points
odb = TRIM(jobname)//'.dat'
WRITE(temp,'(i8)') MPI_ID
fid = 700 + MPI_ID
mpi_odb = TRIM(jobname)//'_'//TRIM(ADJUSTL(temp))
ALLOCATE(Nt(2*nnte, nip), disp(nnse*ndofn, ngp), strain(ngp, ncmp), stress(ngp, ncmp))
!--- Calculate temporal shape function -----------------------------------------
dt = tcoorde(nnte) - tcoorde(1)
nout = nip - nipc/4*3 + 1
DO i = 1, nip
    t = tcoorde(1) + (i-1)*dt/nip
    Nt(1,i) = 2._REKIND*(tcoorde(2)-t)*(tcoorde(3)-t)/dt**2
    Nt(2,i) = -4._REKIND*(tcoorde(1)-t)*(tcoorde(3)-t)/dt**2
    Nt(3,i) = 2._REKIND*(tcoorde(1)-t)*(tcoorde(2)-t)/dt**2
    Nt(4,i) = Nt(1,i)*(DSIN(w*t)-DSIN(w*tcoorde(1)))
    Nt(5,i) = Nt(2,i)*(DSIN(w*t)-DSIN(w*tcoorde(2)))
    Nt(6,i) = Nt(3,i)*(DSIN(w*t)-DSIN(w*tcoorde(3)))
ENDDO
!--- Evaluate strain-stress ----------------------------------------------------
DO ii = 1, ngp ! Loop over spatial integration points
    noel = (ii-1)/ngpe + 1   ! Element index
    npt = MOD(ii, ngpe)     ! Gauss point index
    !--- Calculate constitutive matrices (stiffness and compliance) ----------------
    coef = ym(ii)/((1._REKIND+nu(ii))*(1._REKIND-2._REKIND*nu(ii)))
    Dmat = 0._REKIND
    Cmat = 0._REKIND
    DO i = 1, 3
        Dmat(i,i) = coef*(1._REKIND-nu(ii))
        Dmat(i+3,i+3) = coef*(1._REKIND-2._REKIND*nu(ii))/2._REKIND
        Cmat(i,i) = 1._REKIND/ym(ii)
        Cmat(i+3,i+3) = 2._REKIND*(1._REKIND+nu(ii))/ym(ii)
    ENDDO
    DO i = 1, 3
        DO j = 1, 3
            IF (i.NE.j) THEN
                Dmat(i,j) = coef*nu(ii)
                Cmat(i,j) = -nu(ii)/ym(ii)
            ENDIF
        ENDDO
    ENDDO

    IF (npt .EQ. 0) npt = ngpe
    ! Interpolate displacement
    CALL dispinterp(nnte, nnse, nn, ndofn, xconn(noel,:), sol, Nt(:,nout), &
        & disp(:, ii))
    ! Calculate strain at mesoscale
    CALL calstrain(nnse, ndofn, ncmp, disp(:,ii), Nx(:,ii), Ny(:,ii), &
            & Nz(:,ii), Dmat, strain(ii, :), stress(ii,:))
    !If (xconn(noel,1) .EQ. 2513 .or. xconn(noel,2) .EQ. 2513 .or. &
    !   &xconn(noel,3) .EQ. 2513 .or. xconn(noel,4) .EQ. 2513 .or. &
    !   &xconn(noel,5) .EQ. 2513 .or. xconn(noel,6) .EQ. 2513 .or. &
    !   &xconn(noel,7) .EQ. 2513 .or. xconn(noel,8) .EQ. 2513 ) then
    !    PRINT *, ii, 'stressZZ =', stress(ii,3)
    !ENDIF
ENDDO
!--- Clean up ------------------------------------------------------------------
DO i = 1, ngp
    IF (gpsta(i).EQ.0) THEN ! Failed integration points
        strain(i,:) = 0._REKIND
        stress(i,:) = 0._REKIND
    ENDIF
ENDDO
!--- Write Odb -----------------------------------------------------------------
t = tcoorde(1) + (nout-1)*dt/nip
IF (MPI_ID.EQ.0) THEN
    CALL dispinterpt(nn, ndofn, nnte, Nt(:, nout), dispnd, sol)
    CALL tempinterpt(nn, 1, nnte, Nt(:, nout), tempnd, tempst)
    CALL writetime(jobname, t)
    CALL writedisp(jobname, nn, ndofn, dispnd)
    CALL writetemp(jobname, nn, 1, tempnd)
ENDIF
CALL writeother(fid, mpi_odb, ne, ngpe, ncmp, strain, stress, D, Ws, elsta)
DEALLOCATE(Nt, disp, strain, stress)
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


SUBROUTINE writetemp(filename, nn, ndofn, disp)
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
temp = TRIM(filename)//'.temp'
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
END SUBROUTINE writetemp


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


SUBROUTINE writeother(fid, filename, ne, ngpe, ncmp, strain, stress, D, Ws, elsta)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: filename
INTEGER, INTENT(IN) :: ne, ngpe, ncmp, fid
REAL(KIND=REKIND), DIMENSION(ne*ngpe, ncmp), INTENT(IN) :: strain, stress
REAL(KIND=REKIND), DIMENSION(ne*ngpe), INTENT(IN) :: D, Ws
INTEGER, DIMENSION(ne), INTENT(IN) :: elsta
! ---Internal variables---
INTEGER :: ii, jj
CHARACTER(LEN=80) :: temp
INTEGER :: errstat
CHARACTER(LEN=80) :: errmesg

temp = TRIM(filename)//'.strain'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, ne*ngpe ! strain
    WRITE(fid, '(6(E16.8,2X))') (strain(ii, jj), jj = 1, ncmp)
ENDDO
CLOSE(fid)

temp = TRIM(filename)//'.stress'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, ne*ngpe ! stress
    WRITE(fid, '(6(E16.8,2X))') (stress(ii, jj), jj = 1, ncmp)
ENDDO
CLOSE(fid)

temp = TRIM(filename)//'.damage'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, ne*ngpe ! damage
    WRITE(fid, *) D(ii)
ENDDO
CLOSE(fid)

temp = TRIM(filename)//'.energy'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, ne*ngpe ! energy
    WRITE(fid, *) Ws(ii)
ENDDO
CLOSE(fid)

temp = TRIM(filename)//'.status'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), POSITION='append', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
DO ii = 1, ne ! status
    WRITE(fid, *) elsta(ii)
ENDDO
CLOSE(fid)
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
temp = TRIM(filename)//'.temp'
OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(temp)), STATUS='unknown', IOSTAT=errstat, &
    & IOMSG=errmesg)
IF (errstat /= 0) THEN
    WRITE(*,*) 'Error code: ', errstat
    WRITE(*,*) 'Error message: ', errmesg
    STOP
ENDIF
CLOSE(fid)
DO i = 1, MPI_P
    j = i - 1
    WRITE(temp,'(i8)') j
    mpifn0 = TRIM(filename)//'_'//TRIM(ADJUSTL(temp))
    mpifn1 = TRIM(mpifn0)//'.strain'
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mpifn1)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        STOP
    ENDIF
    CLOSE(fid)
    mpifn1 = TRIM(mpifn0)//'.stress'
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mpifn1)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        STOP
    ENDIF
    CLOSE(fid)
    mpifn1 = TRIM(mpifn0)//'.damage'
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mpifn1)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        STOP
    ENDIF
    CLOSE(fid)
    mpifn1 = TRIM(mpifn0)//'.energy'
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mpifn1)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        STOP
    ENDIF
    CLOSE(fid)
    mpifn1 = TRIM(mpifn0)//'.status'
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mpifn1)), STATUS='unknown', IOSTAT=errstat, &
        & IOMSG=errmesg)
    IF (errstat /= 0) THEN
        WRITE(*,*) 'Error code: ', errstat
        WRITE(*,*) 'Error message: ', errmesg
        STOP
    ENDIF
    CLOSE(fid)
ENDDO
RETURN
END SUBROUTINE initfile
