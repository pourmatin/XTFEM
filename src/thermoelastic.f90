SUBROUTINE material_update(nnse, nquads, ndofnst, nn, ne, ngpe, mpi_ne, ele, xconn, &
                & ns, tempst, Tref, ym0, nu0, ym, nu)
USE kinds
!-------------------------------------------------------------------------------
! Declare variables
!-------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: nnse, mpi_ne, nquads, nn, ne, ngpe, ndofnst
REAL(KIND=REKIND), INTENT(IN) :: ym0, nu0, Tref
REAL(KIND=REKIND), DIMENSION(ngpe, nnse), INTENT(IN) :: ns
INTEGER, DIMENSION(mpi_ne), INTENT(IN) :: ele
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(ndofnst), INTENT(IN) :: tempst
REAL(KIND=REKIND), DIMENSION(mpi_ne*ngpe), INTENT(OUT) :: ym, nu
! Counters
INTEGER :: i, j, k, l, ii, jj, idx, noel
REAL(KIND=REKIND) :: temp
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
                    idx = 2*nn + nodes(jj)
                    temp = temp + ns(l, jj) * tempst(idx)
                ENDDO
                ym((ii-1)*ngpe+l) = ym0 - 1.6*(temp-Tref)-0.1*(temp-Tref)*(temp-Tref)
                nu((ii-1)*ngpe+l) = nu0 + 5.e-5*(temp-Tref)
            ENDDO
        ENDDO
    ENDDO
ENDDO
RETURN
ENDSUBROUTINE material_update

SUBROUTINE calc_thermo_stress(ngpe, nn, nnse, nquads, nodes, &
           & ns, tempst, Tref, alpha, stiff, stress)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables--- 
INTEGER, INTENT(IN) :: nnse, nquads, nn, ngpe
REAL(KIND=REKIND), INTENT(IN) :: Tref, alpha
INTEGER, DIMENSION(nnse) :: nodes
REAL(KIND=REKIND), DIMENSION (ngpe, nnse), INTENT(IN) :: ns
REAL(KIND=REKIND), DIMENSION(nn), INTENT(IN) :: tempst
REAL(KIND=REKIND), DIMENSION(6,6,ngpe), INTENT(IN) :: stiff
REAL(KIND=REKIND), DIMENSION(3,ngpe), INTENT(OUT) :: stress
! Counters
INTEGER :: l, jj, idx
REAL(KIND=REKIND) :: strain, temp

stress = 0
DO l = 1, nquads**3
    temp = 0.
    DO jj = 1, nnse
        idx = nodes(jj)
        temp = temp + ns(l, jj) * tempst(idx)
    ENDDO
    strain = alpha * (temp - Tref)
    DO jj = 1, 3
        stress(jj,l) = (stiff(jj,1,l) + stiff(jj,2,l) + stiff(jj,3,l)) * strain
    ENDDO
ENDDO
RETURN
END SUBROUTINE calc_thermo_stress

