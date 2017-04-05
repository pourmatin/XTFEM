SUBROUTINE material_update(nnse, nquads, ne, ngpe, mpi_ne, ele, xconn, &
                & ns, tempst, ndofst, ym0, nu0, ym, nu)
!-------------------------------------------------------------------------------
! Declare variables
!-------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: nnse, mpi_ne, nquads, ne, ngpe
REAL(KIND=8), INTENT(IN) :: ym0, nu0
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: ele
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=8), INTENT(IN) :: tempst(*)
REAL(KIND=8), INTENT(OUT) :: ym(*), nu(*)
! Counters
INTEGER :: i, j, k, l, ii, jj, idx, noel
REAL(KIND=8) :: temp
INTEGER, DIMENSION(nnse) :: nodes

l = 0
DO ii = 1, mpi_ne
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
                ym(l) = ym0 - 1.6*temp-0.1*temp*temp
                nu(l) = nu0 + 5.e-5*temp
            ENDDO
        ENDDO
    ENDDO
ENDDO
ENDSUBROUTINE material_update