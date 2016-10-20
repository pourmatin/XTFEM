SUBROUTINE delelement(nn, ne, ndofn, nnse, ngpe, ncmp, nquads, xconn, xcoord, &
    & dns, dnxyz, sd, stiff, wqs, nef, nefail, ki, kj, kv, ndof, knz)
USE kinds
USE procs
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ne, ndofn, nnse, ngpe, ncmp, nquads, ndof
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), DIMENSION(nn, ndofn), INTENT(IN) :: xcoord
REAL(KIND=REKIND), DIMENSION(ndofn*ngpe, nnse), INTENT(IN) :: dns
REAL(KIND=REKIND), DIMENSION(ndofn**2*ngpe, ndofn*nnse), INTENT(IN) :: dnxyz
REAL(KIND=REKIND), DIMENSION(6, 9), INTENT(IN) :: sd
REAL(KIND=REKIND), DIMENSION(ncmp, ncmp), INTENT(IN) :: stiff
REAL(KIND=REKIND), DIMENSION(nquads), INTENT(IN) :: wqs
INTEGER, INTENT(IN) :: nef
INTEGER, DIMENSION(nef), INTENT(IN) :: nefail
INTEGER, DIMENSION(ndof+1), INTENT(INOUT) :: ki
INTEGER, DIMENSION(knz), INTENT(INOUT) :: kj
REAL(KIND=REKIND), DIMENSION(knz), INTENT(INOUT) :: kv
INTEGER, INTENT(INOUT) :: knz
! ---Internal variables---
INTEGER :: i, j, k, l, noel, ndofe
INTEGER, DIMENSION(nnse) :: nodes
REAL(KIND=REKIND), DIMENSION(nnse, ndofn) :: ncoorde
INTEGER, DIMENSION(nnse*ndofn) :: dofs
REAL(KIND=REKIND), DIMENSION(nnse*ndofn, nnse*ndofn) :: ke
INTEGER, DIMENSION(:), ALLOCATABLE :: ia, ja
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: a
INTEGER :: nnz
IF (nef .GT. 0) THEN
    nnz = nef*(nnse*ndofn)**2 
    ndofe = (nnse*ndofn)**2
    ! Allocate memory
    ALLOCATE(ia(MAX(ndof+1,nnz)), ja(nnz), a(nnz))
    DO k = 1, nef
        noel = nefail(k)
        ! Formulate elemental stiffness matrix
        nodes = xconn(noel,:)
        ncoorde = xcoord(nodes,:)
        DO i = 1, nnse
            DO j = 1, ndofn
                dofs((i-1)*ndofn+j) = (nodes(i)-1)*ndofn + j
            ENDDO
        ENDDO
        CALL elekmtx(ndofn, nnse, nquads, ngpe, ncmp, ncoorde, sd, stiff, wqs, &
            & dns, dnxyz, ke)
        ! Assemble elemental matrix to a temporary global matrix
        CALL assembly(ia((k-1)*ndofe+1:k*ndofe), ja((k-1)*ndofe+1:k*ndofe), &
            & a((k-1)*ndofe+1:k*ndofe), ke, dofs, nnse, ndofn)
    ENDDO
    CALL tricsr(ia, ja, a, ndof, nnz, MAX(nnz,ndof+1))
    ! Subtract the temporary matrix from the original global matrix
    CALL cs_add(ki, kj, kv, knz, ia, ja, a, nnz, ndof, 1.0D+00, -1.0D+00)
    ! Clean memory
    DEALLOCATE(ia, ja, a)
ENDIF
RETURN
END SUBROUTINE delelement


SUBROUTINE chkdmgsta(nn, ne, nnse, ngpe, xconn, Dc, D, elesta, gpsta, nef, &
    & nefail, ndf, ndfail, necnt)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ne, nnse, ngpe
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
REAL(KIND=REKIND), INTENT(IN) :: Dc
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D
INTEGER, DIMENSION(ne), INTENT(INOUT) :: elesta
INTEGER, DIMENSION(ngpe*ne), INTENT(INOUT) :: gpsta
INTEGER, INTENT(OUT) :: nef
INTEGER, DIMENSION(ne), INTENT(OUT) :: nefail
INTEGER, INTENT(OUT) :: ndf
INTEGER, DIMENSION(nn), INTENT(OUT) :: ndfail
INTEGER, DIMENSION(nn), INTENT(INOUT) :: necnt
! ---Internal variables---
INTEGER :: i, j, k
nef = 0
nefail = 0
ndf = 0
ndfail = 0
! Determine the status of elements and Gauss points
DO i = 1, ne
    IF (ANY(D(ngpe*(i-1)+1:ngpe*i) .GE. Dc)) then
        gpsta(ngpe*(i-1)+1:ngpe*i) = 0
        elesta(i) = 0
        nef = nef + 1
        nefail(nef) = i
    ENDIF
ENDDO
! Determine the status of nodes
DO i = 1, nef
    k = nefail(i)          
    DO j = 1, nnse
        necnt(xconn(k,j)) = necnt(xconn(k,j)) - 1
        IF (necnt(xconn(k,j)) .EQ. 0) THEN
            ndf = ndf + 1
            ndfail(ndf) = xconn(k,j)
        ENDIF
    ENDDO 
ENDDO
RETURN
END SUBROUTINE chkdmgsta


SUBROUTINE chkelesta(nn, ne, nnse, ngpe, Dc, D, elesta, gpsta, nef, &
    & nefail)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ne, nnse, ngpe
REAL(KIND=REKIND), INTENT(IN) :: Dc
REAL(KIND=REKIND), DIMENSION(ngpe*ne), INTENT(IN) :: D
INTEGER, DIMENSION(ne), INTENT(INOUT) :: elesta
INTEGER, DIMENSION(ngpe*ne), INTENT(INOUT) :: gpsta
INTEGER, INTENT(OUT) :: nef
INTEGER, DIMENSION(ne), INTENT(OUT) :: nefail
! ---Internal variables---
INTEGER :: i, j, k
nef = 0
nefail = 0
! Determine the status of elements and Gauss points
DO i = 1, ne
    IF (ANY(D(ngpe*(i-1)+1:ngpe*i) .GE. Dc)) then
        gpsta(ngpe*(i-1)+1:ngpe*i) = 0
        elesta(i) = 0
        nef = nef + 1
        nefail(nef) = i
    ENDIF
ENDDO
RETURN
END SUBROUTINE chkelesta


SUBROUTINE chkcnsta(nn, ne, nnse, xconn, nef, nefail, necnt)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ne, nnse
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
INTEGER, INTENT(IN) :: nef
INTEGER, DIMENSION(ne), INTENT(IN) :: nefail
INTEGER, DIMENSION(nn), INTENT(INOUT) :: necnt
! ---Internal variables---
INTEGER :: i, j, k
! Determine the status of nodes
IF (nef.GT.0) THEN
    DO i = 1, nef
        k = nefail(i)          
        DO j = 1, nnse
            necnt(xconn(k,j)) = necnt(xconn(k,j)) - 1
        ENDDO 
    ENDDO
ENDIF
RETURN
END SUBROUTINE chkcnsta


SUBROUTINE chkndsta(nn, ne, nnse, xconn, nef, nefail, necnt, ndf, ndfail)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nn, ne, nnse
INTEGER, DIMENSION(ne, nnse), INTENT(IN) :: xconn
INTEGER, INTENT(IN) :: nef
INTEGER, DIMENSION(ne), INTENT(IN) :: nefail
INTEGER, DIMENSION(nn), INTENT(IN) :: necnt
INTEGER, INTENT(OUT) :: ndf
INTEGER, DIMENSION(nn), INTENT(OUT) :: ndfail
! ---Internal variables---
INTEGER :: i, j, k
ndf = 0
ndfail = 0
IF (nef.GT.0) THEN
    ! Determine the status of nodes
    DO i = 1, nef
        k = nefail(i)
        DO j = 1, nnse
            IF (necnt(xconn(k,j)) .EQ. 0) THEN
                IF (ndf .LT. 1) THEN
                    ndf = ndf + 1
                    ndfail(ndf) = xconn(k,j)
                ELSE
                    IF (ANY(xconn(k,j).EQ.ndfail(1:ndf))) THEN
                    ELSE
                        ndf = ndf + 1
                        ndfail(ndf) = xconn(k,j)
                    ENDIF
                ENDIF
            ENDIF
        ENDDO 
    ENDDO
ENDIF
RETURN
END SUBROUTINE chkndsta


SUBROUTINE adddiagonal(nnse, ndofn, ndf, ndfail, ki, kj, kv, ndof, knz)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nnse, ndofn, ndf, ndof
INTEGER, DIMENSION(ndf), INTENT(IN) :: ndfail
INTEGER, DIMENSION(ndof+1), INTENT(INOUT) :: ki
INTEGER, DIMENSION(knz), INTENT(INOUT) :: kj
REAL(KIND=REKIND), DIMENSION(knz), INTENT(INOUT) :: kv
INTEGER, INTENT(INOUT) :: knz
! ---Internal variables---
INTEGER :: i, j, k, l
INTEGER, DIMENSION(:), ALLOCATABLE :: ia, ja
REAL(KIND=REKIND), DIMENSION(:), ALLOCATABLE :: aa
INTEGER :: nz
IF (ndf .GT. 0) THEN
    ! Allocate memory
    nz = ndof+1 ! this is not the actual nz, which is ndf*ndofn < ndof
    ALLOCATE(ia(nz), ja(nz), aa(nz))
    ia = 0
    ja = 0
    aa = 0._REKIND
    ! Formulate diagonal matrix
    l = 1
    DO i = 1, ndf
        j = ndfail(i)
        DO k = 1, ndofn
            ia(l) = (j-1)*ndofn + k
            ja(l) = (j-1)*ndofn + k
            aa(l) = 1._REKIND
            l = l + 1 
        ENDDO
    ENDDO
    nz = ndf*ndofn
    CALL tricsr(ia, ja, aa, ndof, nz)
    ! Add the diagonal matrix to stiffness matrix
    knz = knz + nz
    CALL cs_add(ki, kj, kv, knz, ia, ja, aa, nz, ndof, 1.0D+00, 1.0D+00)
    ! Clean memory
    DEALLOCATE(ia, ja, aa)
ENDIF
RETURN
END SUBROUTINE adddiagonal
