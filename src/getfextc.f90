SUBROUTINE getfextc(nstep, w, dt, fextc)
USE kinds
IMPLICIT NONE
! ---External variables---
INTEGER, INTENT(IN) :: nstep
REAL(KIND=REKIND), INTENT(IN) :: w, dt
REAL(KIND=REKIND), DIMENSION(3), INTENT(OUT) :: fextc
! ---Internal variables---
REAL(KIND=REKIND) :: t1
t1=(nstep-1)*dt
fextc = 0._REKIND      
fextc(1) = -((3._REKIND*dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+4._REKIND*          &
    & SIN(t1*w)-4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
fextc(2) = (4._REKIND*(dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+2._REKIND*SIN(t1*w)  &
    & -2._REKIND*SIN((dt+t1)*w)))/(dt**2*w**2)
fextc(3) = -((dt*w*COS(t1*w)+3._REKIND*dt*w*COS((dt+t1)*w)+4*SIN(t1*w)-        &
    & 4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
RETURN
END SUBROUTINE getfextc
