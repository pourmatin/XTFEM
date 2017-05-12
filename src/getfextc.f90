SUBROUTINE getfextc_t(t1, dt, fextc)
USE kinds
IMPLICIT NONE
! ---External variables---
REAL(KIND=REKIND), INTENT(IN) :: dt, t1
REAL(KIND=REKIND), DIMENSION(3), INTENT(OUT) :: fextc
! ---Internal variables---
REAL(KIND=REKIND) :: w
w = 2*3.14159/24
fextc(1) = t1*dt/6
fextc(2) = dt*(2*t1+dt)/6
fextc(3) = dt*(t1+dt)/6
fextc(1) = -((3._REKIND*dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+4._REKIND*          &
    & SIN(t1*w)-4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
fextc(2) = (4._REKIND*(dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+2._REKIND*SIN(t1*w)  &
    & -2._REKIND*SIN((dt+t1)*w)))/(dt**2*w**2)
fextc(3) = -((dt*w*COS(t1*w)+3._REKIND*dt*w*COS((dt+t1)*w)+4*SIN(t1*w)-        &
    & 4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
RETURN
END SUBROUTINE getfextc_t

SUBROUTINE getfextc_m(t1, w, dt, fextc, fint)
USE kinds
IMPLICIT NONE
! ---External variables---
REAL(KIND=REKIND), INTENT(IN) :: w, dt, t1
REAL(KIND=REKIND), DIMENSION(6), INTENT(OUT) :: fextc, fint
! ---Internal variables---
fextc = 0._REKIND      
fextc(1) = -((3._REKIND*dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+4._REKIND*          &
    & SIN(t1*w)-4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
fextc(2) = (4._REKIND*(dt*w*COS(t1*w)+dt*w*COS((dt+t1)*w)+2._REKIND*SIN(t1*w)  &
    & -2._REKIND*SIN((dt+t1)*w)))/(dt**2*w**2)
fextc(3) = -((dt*w*COS(t1*w)+3._REKIND*dt*w*COS((dt+t1)*w)+4*SIN(t1*w)-        &
    & 4._REKIND*SIN((dt+t1)*w))/(dt**2*w**2))
fextc(4) = (16._REKIND-4._REKIND*dt**2*w**2+2._REKIND*(-7._REKIND+dt**2*w**2)* &
    & COS(2._REKIND*t1*w)-2._REKIND*COS(2._REKIND*(dt+t1)*w)+8._REKIND*dt*w*   &
    & COS((dt+t1)*w)*SIN(t1*w)+9._REKIND*dt*w*SIN(2._REKIND*t1*w)-32._REKIND*  &
    & SIN(t1*w)*SIN((dt+t1)*w)-dt*w*SIN(2._REKIND*(dt+t1)*w))/                 &
    & (8._REKIND*dt**2*w**2) 
fextc(5) = (-COS(2._REKIND*t1*w)+COS(2._REKIND*(dt+t1)*w)+8._REKIND*COS(((dt+  &
    & 4._REKIND*t1)*w)/2._REKIND)-8._REKIND*COS((3._REKIND*dt*w)/2._REKIND+    &
    & 2._REKIND*t1*w)+dt*w*(SIN(2._REKIND*t1*w)+SIN(2._REKIND*(dt+t1)*w)-      &
    & 4._REKIND*(SIN(((dt+4._REKIND*t1)*w)/2._REKIND)+SIN((3._REKIND*dt*w)/    &
    & 2._REKIND+2._REKIND*t1*w))))/(2._REKIND*dt**2*w**2)
fextc(6) = (-16._REKIND+4._REKIND*dt**2*w**2+2._REKIND*COS(2._REKIND*t1*w)-    &
    & 2._REKIND*(-7._REKIND+dt**2*w**2)*COS(2._REKIND*(dt+t1)*w)-dt*w*         &
    & SIN(2._REKIND*t1*w)+8._REKIND*(dt*w*COS(t1*w)+4._REKIND*SIN(t1*w))*      &
    & SIN((dt+t1)*w)+9._REKIND*dt*w*SIN(2._REKIND*(dt+t1)*w))/(8._REKIND*      &
    & dt**2*w**2)
fint(1) = -1._REKIND
fint(2) = 0._REKIND
fint(3) = 1._REKIND
fint(4) = 0._REKIND
fint(5) = 0._REKIND
fint(6) = -SIN(t1*w) + SIN(w*(t1 + dt))
RETURN
END SUBROUTINE getfextc_m
