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
