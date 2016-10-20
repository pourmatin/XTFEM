! This file contains subroutines for debugging
! Developed by Rui @ UTD 2/22/2016


SUBROUTINE outputcsrmat(ia, ja, aa, n, nz, fn)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n, nz
INTEGER, DIMENSION(n+1), INTENT(IN) :: ia
INTEGER, DIMENSION(nz), INTENT(IN) :: ja
REAL(KIND=REKIND), DIMENSION(nz), INTENT(IN) :: aa
CHARACTER(LEN=80), INTENT(IN) :: fn
! ---Internal variables---
INTEGER :: i
OPEN(UNIT = 200, FILE = TRIM(fn)//'.rowp', STATUS = 'unknown')
OPEN(UNIT = 210, FILE = TRIM(fn)//'.cols', STATUS = 'unknown')
OPEN(UNIT = 220, FILE = TRIM(fn)//'.vals', STATUS = 'unknown')
DO i = 1, n+1
    WRITE(200,*) ia(i)
ENDDO
DO i = 1, nz
    WRITE(210,*) ja(i)
    WRITE(220,*) aa(i)
ENDDO
CLOSE(200)
CLOSE(210)
CLOSE(220)
RETURN
END SUBROUTINE outputcsrmat


SUBROUTINE outputrealvec(a, n, fn)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n
REAL(KIND=REKIND), DIMENSION(n), INTENT(IN) :: a
CHARACTER(LEN=80), INTENT(IN) :: fn
! ---Internal variables---
INTEGER :: i
OPEN(UNIT = 200, FILE = TRIM(fn)//'.realvec', STATUS = 'unknown')
DO i = 1, n
    WRITE(200,*) a(i)
ENDDO
CLOSE(200)
RETURN
END SUBROUTINE outputrealvec


SUBROUTINE outputintvec(a, n, fn)
USE kinds
!! Variable declaration
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: n
INTEGER, DIMENSION(n), INTENT(IN) :: a
CHARACTER(LEN=80), INTENT(IN) :: fn
! ---Internal variables---
INTEGER :: i
OPEN(UNIT = 200, FILE = TRIM(fn)//'.intvec', STATUS = 'unknown')
DO i = 1, n
    WRITE(200,*) a(i)
ENDDO
CLOSE(200)
RETURN
END SUBROUTINE outputintvec


SUBROUTINE pause_exe()
IMPLICIT NONE

print *, 'warning: program pause'
DO WHILE(.TRUE.)
ENDDO
CALL SLEEP(120)
print *, 'warning: time out, program exit'
STOP
END SUBROUTINE pause_exe
