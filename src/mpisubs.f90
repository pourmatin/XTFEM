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
n = CEILING(REAL(ne)/REAL(MPI_P))*MPI_ID
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


SUBROUTINE mpiglbmtx(mpi_ne, ele, ne, nn, ndofn, nnse, nquads, ngpe, &
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
    CALL elemtx(ndofn, nnse, nquads, ngpe, ncoorde, rho, kappa, cp, &
        & wqs, ns, dns, ke, me)
!print *, 'ke:'
!    do, j=1,8
!        write(*,'(100g15.5)') ( ke(j,k), k=1,8 )
!    enddo
!print *, 'me:'
!    do, j=1,8
!        write(*,'(100g15.5)') ( me(j,k), k=1,8 )
!    enddo
!print *, "i= ",i
    CALL assembly(ki((i-1)*ndofe+1:i*ndofe), kj((i-1)*ndofe+1:i*ndofe), & 
        & kv((i-1)*ndofe+1:i*ndofe), ke, dofs, nnse, ndofn)
    CALL assembly(mi((i-1)*ndofe+1:i*ndofe), mj((i-1)*ndofe+1:i*ndofe), &
        & mv((i-1)*ndofe+1:i*ndofe), me, dofs, nnse, ndofn)
ENDDO
RETURN
END SUBROUTINE mpiglbmtx
