SUBROUTINE cpu_damage(jobname, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, &
    & xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, gpsta, eps_old, eps_pl, &
    & sig_eff, sig_back, p, D, ws, del_po, sig_d)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: ym, nu, w
REAL(KIND=REKIND), DIMENSION(nnte), INTENT(IN) :: tcoorde
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
INTEGER, DIMENSION(ngpe*ne), INTENT(IN) :: gpsta
REAL(KIND=REKIND), DIMENSION(ncmp, ngpe*ne), INTENT(INOUT) :: eps_old, eps_pl, &
    & sig_eff, sig_back
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(INOUT) :: p, D, ws, del_po, sig_d
! ---Internal variables---
CHARACTER(LEN=80) :: odb
INTEGER :: ngp, noel
REAL(KIND=REKIND) :: t, dt, coef
REAL(KIND=REKIND), DIMENSION(ncmp, ncmp) :: Dmat, Cmat
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nt, strain, stress, disp
INTEGER :: i, j, k, ii, jj
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
!--- Evaluate two-scale damage model -------------------------------------------
!$omp parallel num_threads(24)
!$omp do private(noel)
DO ii = 1, ngp ! Loop over spatial integration points
    noel = (ii-1)/ngpe + 1   ! Element index
    DO jj = 1, nip ! Loop over temporal interpolation points
        ! Interpolate displacement
        CALL dispinterp(nnte, nnse, nn, ndofn, xconn(noel,:), sol, Nt(:,jj), &
            & disp(:, ii))            
        ! Calculate strain at mesoscale
        CALL calstrain(nnse, ndofn, ncmp, disp(:,ii), Nx(:,ii), Ny(:,ii), &
            & Nz(:,ii), Dmat, strain(ii, :), stress(ii,:))
        ! Evaluate damage
        CALL getdfinal(eps_old(:,ii), strain(ii, :), eps_pl(:,ii), p(ii), &
            & D(ii), ws(ii), sig_eff(:,ii), sig_back(:,ii), del_po(ii), &
            & sig_d(ii), Cmat, Dmat)
    ENDDO
ENDDO
!$omp end do
!$omp end parallel
!--- Clean up ------------------------------------------------------------------
DO i = 1, ngp
    IF (gpsta(i).EQ.0) THEN ! Failed integration points
        D(i) = 0._REKIND
        ws(i) = 0._REKIND
    ENDIF
ENDDO
DEALLOCATE(Nt, disp, strain, stress)
RETURN
END SUBROUTINE cpu_damage



SUBROUTINE cpu_damage_c(jobname, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip, &
    & xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, gpsta, eps_old, eps_pl, &
    & sig_eff, sig_back, p, D, ws, del_po, sig_d)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: ym, nu, w
REAL(KIND=REKIND), DIMENSION(nnte), INTENT(IN) :: tcoorde
REAL(KIND=REKIND), DIMENSION(2*nnte*nn*ndofn), INTENT(IN) :: sol
INTEGER, DIMENSION(ngpe*ne), INTENT(IN) :: gpsta
REAL(KIND=REKIND), DIMENSION(ncmp, ngpe*ne), INTENT(INOUT) :: eps_old, eps_pl, &
    & sig_eff, sig_back
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(INOUT) :: p, D, ws, del_po, sig_d
! ---Internal variables---
CHARACTER(LEN=80) :: odb
INTEGER :: ngp, noel, npt
REAL(KIND=REKIND) :: t, dt
REAL(KIND=REKIND), DIMENSION(:,:), ALLOCATABLE :: Nt
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: c_Nx, c_Ny, c_Nz, c_Nt, &
    & c_eps_old, c_eps_pl, c_sig_eff, c_sig_back
INTEGER, DIMENSION(:), ALLOCATABLE :: c_nodes
INTEGER :: i, j, k
!--- Set parameters & allocate memory ------------------------------------------
ngp = ngpe*ne ! number of total Gauss points
odb = TRIM(jobname)//'.dat'
ALLOCATE(Nt(2*nnte, nip), c_Nx(nnse*ngp), c_Ny(nnse*ngp), c_Nz(nnse*ngp), &
    & c_nodes(nnse*ngp), c_Nt(2*nnte*nip), c_eps_old(ncmp*ngp), &
    & c_eps_pl(ncmp*ngp), c_sig_eff(ncmp*ngp), c_sig_back(ncmp*ngp))
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
!--- Convert 2D or 3D array to 1D array ----------------------------------------
k = 1
DO i = 1, ngp
    noel = (i-1)/ngpe + 1
    DO j = 1, nnse
        c_nodes(k) = xconn(noel, j)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, ngp
    DO j = 1, nnse
        c_Nx(k) = Nx(j,i)
        c_Ny(k) = Ny(j,i)
        c_Nz(k) = Nz(j,i)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, nip
    DO j = 1, 2*nnte
        c_Nt(k) = Nt(j,i)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, ngp
    DO j = 1, ncmp
        c_eps_old(k) = eps_old(j,i)
        c_eps_pl(k) = eps_pl(j,i)
        c_sig_eff(k) = sig_eff(j,i)
        c_sig_back(k) = sig_back(j,i)
        k = k + 1
    ENDDO
ENDDO
!--- Call c host -------------------------------------------------------------
CALL cpu_host(%VAL(ne), %VAL(nn), %VAL(ndofn), %VAL(nnse), %VAL(nnte), &
    & % VAL(ngpe), %VAL(ncmp), %VAL(nip), c_nodes, c_Nx, c_Ny, c_NZ, c_Nt, &
    & sol, c_eps_old, c_eps_pl, c_sig_eff, c_sig_back, p, D, ws, del_po, sig_d)
!--- Convert 1D array to 2D array ----------------------------------------------
k = 1
DO i = 1, ngp
    DO j = 1, ncmp
        eps_old(j,i) = c_eps_old(k)
        eps_pl(j,i) = c_eps_pl(k)
        sig_eff(j,i) = c_sig_eff(k)
        sig_back(j,i) = c_sig_back(k)
        k = k + 1
    ENDDO
ENDDO
!--- Clean up ------------------------------------------------------------------
DO i = 1, ngp
    IF (gpsta(i).EQ.0) THEN ! Failed integration points
        D(i) = 0._REKIND
        ws(i) = 0._REKIND
    ENDIF
ENDDO
DEALLOCATE(Nt, c_nodes, c_Nx, c_Ny, c_Nz, c_Nt, c_eps_old, c_eps_pl, &
    & c_sig_eff, c_sig_back)
RETURN
END SUBROUTINE cpu_damage_c


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
