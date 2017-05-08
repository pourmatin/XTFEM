! Fri Apr 01 19:48:42 CDT 2016
! Subroutines for MPI

SUBROUTINE mpigroupsize(ne, MPI_P, MPI_ID, mpi_ne)
! Divide elements into groups and distribute to all threads
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ne, MPI_P, MPI_ID
INTEGER, INTENT(OUT) :: mpi_ne
! ---Internal variables---
INTEGER :: n
n = CEILING(REAL(ne)/REAL(MPI_P))
IF (n*(MPI_P-1) .GT. ne) THEN
    n = FLOOR(REAL(ne)/REAL(MPI_P))
ENDIF
IF (MPI_ID .LT. (MPI_P-1)) THEN
    mpi_ne = n
ELSE
    mpi_ne = ne - (MPI_P-1)*n
ENDIF
RETURN
END SUBROUTINE mpigroupsize


SUBROUTINE mpisetele(ne, MPI_P, MPI_ID, mpi_ne, mpi_ele)
! Divide elements into groups and distribute to all threads
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ne, MPI_P, MPI_ID, mpi_ne
INTEGER, DIMENSION(mpi_ne), INTENT(OUT) :: mpi_ele
! ---Internal variables---
INTEGER :: i, n
n = CEILING(REAL(ne)/REAL(MPI_P))
IF (n*(MPI_P-1) .GT. ne) THEN
    n = FLOOR(REAL(ne)/REAL(MPI_P))
ENDIF
n = n*MPI_ID
DO i = 1, mpi_ne
    mpi_ele(i) = n+i
ENDDO
RETURN
END SUBROUTINE mpisetele


SUBROUTINE mpisetgp(ngpe, mpi_ne, mpi_ele, mpi_ngp, mpi_gp)
! Divide elements into groups and distribute to all threads
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: ngpe, mpi_ne, mpi_ngp
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: mpi_ele
INTEGER, DIMENSION(mpi_ngp), INTENT(OUT) :: mpi_gp
! ---Internal variables---
INTEGER :: i, j, k
DO i = 1, mpi_ne
    k = mpi_ele(i)
    DO j = 1, ngpe
        mpi_gp((i-1)*ngpe+j) = (k-1)*ngpe+j
    ENDDO
ENDDO
RETURN
END SUBROUTINE mpisetgp

SUBROUTINE mpiglbmtx_m(mpi_ne, ele, ne, nn, ndofn, nnse, nquads, ngpe, ncmp, &
    & xconn, xcoord, ym, nu, ym0, nu0, rho, alpha, Tref, wqs, ns_t, ns, dns_t, dns, dnxyz, &
    &sd, tempst, stiff, Nx, Ny, Nz, ki, kj, kv, knz, mi, mj, mv, mnz, vv)
! Assembly global matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: mpi_ne, ne, nn, ndofn, nnse, nquads, ngpe, ncmp
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: ele
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nn, 3), INTENT(IN) :: xcoord
REAL(KIND=REKIND), INTENT(IN) :: rho, alpha, Tref, ym0, nu0
REAL(KIND=REKIND), DIMENSION(ngpe*mpi_ne), INTENT(IN) :: ym, nu
REAL(KIND=REKIND), DIMENSION(nn), INTENT(IN) :: tempst
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns_t
REAL(KIND=REKIND), DIMENSION(ngpe, nnse, 3), INTENT(IN) :: dns_t
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(6,9), INTENT(OUT) :: sd
REAL(KIND=REKIND), DIMENSION(6,6,mpi_ne*ngpe), INTENT(OUT) :: stiff
REAL(KIND=REKIND), DIMENSION(nnse, ngpe*ne), INTENT(OUT) :: Nx, Ny, Nz
INTEGER, INTENT(IN) :: knz
REAL(KIND=REKIND), DIMENSION(knz), INTENT(OUT) :: kv
INTEGER, DIMENSION(knz), INTENT(OUT) :: ki, kj
INTEGER, INTENT(IN) :: mnz
INTEGER, DIMENSION(mnz), INTENT(OUT) :: mi, mj
REAL(KIND=REKIND), DIMENSION(mnz), INTENT(OUT) :: mv
REAL(KIND=REKIND), DIMENSION(nn*ndofn), INTENT(OUT) :: vv
! ---Internal variables---
INTEGER :: i, j, k, noel, ndofe, ii, jj
INTEGER, DIMENSION(nnse) :: nodes
INTEGER, DIMENSION(nnse*ndofn) :: dofs
REAL(KIND=REKIND), DIMENSION(ndofn,ngpe) :: stress_t
REAL(KIND=REKIND), DIMENSION(nnse, ndofn) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ke, me
REAL(KIND=REKIND), DIMENSION(nnse*ndofn) :: ve
ndofe = (nnse*ndofn)**2
vv = 0._REKIND
DO i = 1, mpi_ne
    noel = ele(i)
    nodes = xconn(noel,:)
    ncoorde = xcoord(nodes,:)
    DO j = 1, nnse
        DO k = 0, ndofn-1
            dofs(j*ndofn-k) = nodes(j)*ndofn - k
        ENDDO
    ENDDO
    DO j = 1, ngpe
        CALL setsdstiff(ym((i-1)*ngpe+j), nu((i-1)*ngpe+j), sd, &
             & stiff(:,:,((i-1)*ngpe+j)))
    ENDDO
    CALL elemtx_m(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, rho, sd, stiff(:,:,(i-1)*ngpe+1:i*ngpe), &
        & wqs, ns, dns, dnxyz, Nx(:, (i-1)*ngpe+1:i*ngpe), &
        & Ny(:, (i-1)*ngpe+1:i*ngpe), Nz(:, (i-1)*ngpe+1:i*ngpe), ke, me)
    CALL calc_thermo_stress(ngpe, nn, nnse, nquads, nodes, ns_t, tempst, &
        & Tref, alpha, stiff(:,:,(i-1)*ngpe+1:i*ngpe), stress_t)
    CALL elevec_m(ndofn, nnse, nquads, ngpe, ncoorde, wqs, stress_t, ns_t, dns_t, ve)
    DO j = 1, nnse*ndofn
        vv(dofs(j)) = vv(dofs(j)) + ve(j)
    ENDDO
    CALL assembly(ki((i-1)*ndofe+1:i*ndofe), kj((i-1)*ndofe+1:i*ndofe), & 
        & kv((i-1)*ndofe+1:i*ndofe), ke, dofs, nnse, ndofn)
    CALL assembly(mi((i-1)*ndofe+1:i*ndofe), mj((i-1)*ndofe+1:i*ndofe), &
        & mv((i-1)*ndofe+1:i*ndofe), me, dofs, nnse, ndofn)
ENDDO
RETURN
END SUBROUTINE mpiglbmtx_m

SUBROUTINE mpiglbmtx_t(mpi_ne, ele, ne, nn, ndofn, nnse, nquads, ngpe, &
    & xconn, xcoord, kappa, cp, rho, wqs, ns, dns, ki, kj, kv, knz, mi, mj, mv, mnz)
! Assembly global matrices
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: mpi_ne, ne, nn, ndofn, nnse, nquads, ngpe
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: ele
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nn, 3), INTENT(IN) :: xcoord
REAL(KIND=REKIND), INTENT(IN) :: kappa, cp, rho
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, ndofn*nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
INTEGER, INTENT(IN) :: knz
REAL(KIND=REKIND), DIMENSION(knz), INTENT(OUT) :: kv
INTEGER, DIMENSION(knz), INTENT(OUT) :: ki, kj
INTEGER, INTENT(IN) :: mnz
INTEGER, DIMENSION(mnz), INTENT(OUT) :: mi, mj
REAL(KIND=REKIND), DIMENSION(mnz), INTENT(OUT) :: mv
! ---Internal variables---
INTEGER :: i, j, k, noel, ndofe
INTEGER, DIMENSION(nnse) :: nodes
INTEGER, DIMENSION(nnse*ndofn) :: dofs
REAL(KIND=REKIND), DIMENSION(nnse, 3) :: ncoorde
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ke, me
ndofe = (nnse*ndofn)**2
!CALL setsdstiff(kappa, cp, sd, stiff)
DO i = 1, mpi_ne
    noel = ele(i)
    nodes = xconn(noel,:)
    ncoorde = xcoord(nodes,:)
    DO j = 1, nnse
        DO k = 0, ndofn-1
            dofs(j*ndofn-k) = nodes(j)*ndofn - k
        ENDDO
    ENDDO
    CALL elemtx_t(ndofn, nnse, nquads, ngpe, ncoorde, rho, kappa, cp, &
        & wqs, ns, dns, ke, me)
    CALL assembly(ki((i-1)*ndofe+1:i*ndofe), kj((i-1)*ndofe+1:i*ndofe), & 
        & kv((i-1)*ndofe+1:i*ndofe), ke, dofs, nnse, ndofn)
    CALL assembly(mi((i-1)*ndofe+1:i*ndofe), mj((i-1)*ndofe+1:i*ndofe), &
        & mv((i-1)*ndofe+1:i*ndofe), me, dofs, nnse, ndofn)
ENDDO
RETURN
END SUBROUTINE mpiglbmtx_t
