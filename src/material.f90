SUBROUTINE material_update(nnse, nquads, ndofnst, nn, ne, ngpe, mpi_ne, ele, xconn, &
                & ns, tempst, ym0, nu0, ym, nu)
USE kinds
!-------------------------------------------------------------------------------
! Declare variables
!-------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: nnse, mpi_ne, nquads, nn, ne, ngpe, ndofnst
REAL(KIND=8), INTENT(IN) :: ym0, nu0
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: ele
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=8), DIMENSION(ndofnst), INTENT(IN) :: tempst
REAL(KIND=8), DIMENSION(mpi_ne*ngpe), INTENT(OUT) :: ym, nu
! Counters
INTEGER :: i, j, k, l, ii, jj, idx, noel
REAL(KIND=8) :: temp
INTEGER, DIMENSION(nnse) :: nodes

DO ii = 1, mpi_ne
    l = 0
    noel = ele(ii)
    nodes = xconn(noel,:)
    DO i = 1, nquads
        DO j = 1, nquads
            DO k = 1, nquads
                l = l + 1
                temp = 0.
                DO jj = 1, nnse
                    idx = nn + nodes(jj)
                    temp = temp + ns(l, jj) * tempst(idx)
                ENDDO
                ym((ii-1)*ngpe+l) = ym0 - 1.6*temp-0.1*temp*temp
                nu((ii-1)*ngpe+l) = nu0 + 5.e-5*temp
            ENDDO
        ENDDO
    ENDDO
ENDDO
ENDSUBROUTINE material_update
