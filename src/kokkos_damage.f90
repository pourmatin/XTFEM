SUBROUTINE kokkos_damage(rank, jobname, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, &
    & nip, xconn, Nx, Ny, Nz, ym, nu, w, tcoorde, sol, gpsta, eps_old, &
    & eps_pl, sig_eff, sig_back, p, D, ws, del_po, sig_d)
USE kinds
IMPLICIT NONE
! ---External variables---
CHARACTER(LEN=80), INTENT(IN) :: jobname
INTEGER, INTENT(IN) :: rank, ne, nn, ndofn, nnse, nnte, ngpe, ncmp, nip
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=REKIND), INTENT(IN) :: w
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: ym, nu
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
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: gpu_Nx, gpu_Ny, gpu_Nz, &
    & gpu_Nt, gpu_eps_old, gpu_eps_pl, gpu_sig_eff, gpu_sig_back
INTEGER, DIMENSION(:), ALLOCATABLE :: gpu_nodes
INTEGER :: i, j, k
!--- Set parameters & allocate memory ------------------------------------------
ngp = ngpe*ne ! number of total Gauss points
odb = TRIM(jobname)//'.dat'
ALLOCATE(Nt(2*nnte, nip), gpu_Nx(nnse*ngp), gpu_Ny(nnse*ngp), &
    & gpu_Nz(nnse*ngp), gpu_nodes(nnse*ngp), gpu_Nt(2*nnte*nip), &
    & gpu_eps_old(ncmp*ngp), gpu_eps_pl(ncmp*ngp), gpu_sig_eff(ncmp*ngp), &
    & gpu_sig_back(ncmp*ngp))
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
        gpu_nodes(k) = xconn(noel, j)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, ngp
    DO j = 1, nnse
        gpu_Nx(k) = Nx(j,i)
        gpu_Ny(k) = Ny(j,i)
        gpu_Nz(k) = Nz(j,i)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, nip
    DO j = 1, 2*nnte
        gpu_Nt(k) = Nt(j,i)
        k = k + 1
    ENDDO
ENDDO
k = 1
DO i = 1, ngp
    DO j = 1, ncmp
        gpu_eps_old(k) = eps_old(j,i)
        gpu_eps_pl(k) = eps_pl(j,i)
        gpu_sig_eff(k) = sig_eff(j,i)
        gpu_sig_back(k) = sig_back(j,i)
        k = k + 1
    ENDDO
ENDDO
!--- Call KOKKOS host -------------------------------------------------------------
CALL kokkos_host(%VAL(rank), %VAL(ne), %VAL(nn), %VAL(ndofn), %VAL(nnse), &
    & %VAL(nnte), %VAL(ngpe), %VAL(ncmp), %VAL(nip), gpu_nodes, gpu_Nx, &
    & gpu_Ny, gpu_Nz, gpu_Nt, sol, ym, nu, gpu_eps_old, gpu_eps_pl, gpu_sig_eff, &
    & gpu_sig_back, p, D, ws, del_po, sig_d)
!--- Convert 1D array to 2D array ----------------------------------------------
k = 1
DO i = 1, ngp
    DO j = 1, ncmp
        eps_old(j,i) = gpu_eps_old(k)
        eps_pl(j,i) = gpu_eps_pl(k)
        sig_eff(j,i) = gpu_sig_eff(k)
        sig_back(j,i) = gpu_sig_back(k)
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
DEALLOCATE(Nt, gpu_nodes, gpu_Nx, gpu_Ny, gpu_Nz, gpu_Nt, gpu_eps_old, &
    & gpu_eps_pl, gpu_sig_eff, gpu_sig_back)

RETURN
END SUBROUTINE kokkos_damage
