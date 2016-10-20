SUBROUTINE stcoeff(t1, dt, w, k, m, ktn, mtn, ktn1, mtn1)
! Exact temporal shape function integration by Mathematica
!   k = stiffness matrix of the complete slab
!   m = mass matrix of the complete slab 
!   ktn = stiffness matrix at tn_1 (+)
!   mtn = mass matric at tn_1 (+)
!   ktn1 = stiffness matrix at tn_1 (-)
!   mtne1 = mass matrix at tn_1 (-)
!   w = frequency of the forcing function
!   dt = timestep
USE kinds
IMPLICIT NONE
! ---External variables---
REAL(KIND=REKIND), INTENT(IN) :: t1, dt, w
REAL(KIND=REKIND), DIMENSION(6, 6), INTENT(OUT) :: k, m, ktn, mtn, ktn1, mtn1
! ---Internal variables---
INTEGER :: i, j

! Coefficients for Stiffness matrix of the complete slab
k(1,1) = -0.5_REKIND
k(1,2) = -2._REKIND/3._REKIND
k(1,3) = 1._REKIND/6._REKIND
k(1,4) = (-6._REKIND*dt*w*(-12._REKIND+dt**2*w**2)*COS(t1*w)+24._REKIND*dt*w*  &
    & COS((dt+t1)*w)+(96._REKIND-26._REKIND*dt**2*w**2+dt**4*w**4)*SIN(t1*w)+  &
    & 2._REKIND*(-48._REKIND+dt**2*w**2)*SIN((dt+t1)*w))/(2._REKIND*dt**4*w**4)  
k(1,5) = (2._REKIND*(dt*w*(-84._REKIND*COS(t1*w)-60._REKIND*COS((dt+t1)*w)+    &
    & dt*w*(18._REKIND*SIN(t1*w)+dt**2*w**2*SIN((dt/2._REKIND+t1)*w)-          &
    & 6._REKIND*SIN((dt+t1)*w)))-144._REKIND*(SIN(t1*w)-SIN((dt+t1)*w))))/     &
    & (3._REKIND*dt**4*w**4)
k(1,6) = (120._REKIND*dt*w*COS(t1*w)-6._REKIND*dt*w*(-28._REKIND+dt**2*w**2)*  &
    & COS((dt+t1)*w)-18._REKIND*(-16._REKIND+dt**2*w**2)*SIN(t1*w)-            &
    & (288._REKIND-42._REKIND*dt**2*w**2+dt**4*w**4)*SIN((dt+t1)*w))/          &
    & (6._REKIND*dt**4*w**4)
	
k(2,1) = 2._REKIND/3._REKIND
k(2,2) = 0.0_REKIND
k(2,3) = -2._REKIND/3._REKIND
k(2,4) = (-2._REKIND*(-6._REKIND*dt*w*(-16._REKIND+dt**2*w**2)*COS(t1*w)+      &
    & 48._REKIND*dt*w*COS((dt+t1)*w)+(-24._REKIND+dt**2*w**2)*((-6._REKIND+    &
    & dt**2*w**2)*SIN(t1*w)+6._REKIND*SIN((dt+t1)*w))))/(3._REKIND*dt**4*w**4)
k(2,5) = (16._REKIND*(6._REKIND*dt*w*COS(t1*w)+6._REKIND*dt*w*COS((dt+t1)*w)-  &
    & (-12._REKIND+dt**2*w**2)*(SIN(t1*w)-SIN((dt+t1)*w))))/(dt**4*w**4)
k(2,6) = (2._REKIND*(-48._REKIND*dt*w*COS(t1*w)+6._REKIND*dt*w*(-16._REKIND    &
    & +dt**2*w**2) *COS((dt+t1)*w)+(-24._REKIND+dt**2*w**2)*(6._REKIND*        &
    & SIN(t1*w)+(-6._REKIND+dt**2*w**2)*SIN((dt+t1)*w))))/(3._REKIND*dt**4*    &
    & w**4)
     
k(3,1) = -1._REKIND/6._REKIND
k(3,2) = 2._REKIND/3._REKIND
k(3,3) = 1._REKIND/2._REKIND
k(3,4) = (-6._REKIND*dt*w*(-28._REKIND+dt**2*w**2)*COS(t1*w)+120._REKIND*dt*w* &
    & COS((dt+t1)*w)+(288._REKIND-42._REKIND*dt**2*w**2+dt**4*w**4)*SIN(t1*w)+ &
    & 18._REKIND*(-16._REKIND+dt**2*w**2)*SIN((dt+t1)*w))/(6._REKIND*dt**4     &
    & *w**4)
k(3,5) = (-2._REKIND*(144._REKIND*SIN(t1*w)+dt*w*(60._REKIND*COS(t1*w)+        &
    & 84._REKIND*COS((dt+t1)*w)-6._REKIND*dt*w*SIN(t1*w)+dt**3*w**3*SIN((dt    &
    & /2._REKIND+t1)*w))+18._REKIND*(-8._REKIND+dt**2*w**2)*SIN((dt+t1)*w)))/  &
    & (3._REKIND*dt**4*w**4) 
k(3,6) = (96._REKIND*(SIN(t1*w)-SIN((dt+t1)*w))+dt*w*(24._REKIND*COS(t1*w)+    &
    & (72._REKIND-6._REKIND*dt**2*w**2)*COS((dt+t1)*w)+dt*w*(-2._REKIND*       &
    & SIN(t1*w)+(26._REKIND-dt**2*w**2)*SIN((dt+t1)*w))))/(2._REKIND*dt**4*    &
    & w**4)

k(4,1) = -k(1,4)
k(4,2) = -k(2,4)
k(4,3) = -k(3,4)
k(4,4) = 0._REKIND
k(4,5) = (576._REKIND-72._REKIND*dt**2*w**2-4._REKIND*dt**4*w**4+6._REKIND*    &
    & (-96._REKIND+13._REKIND*dt**2*w**2)*COS(t1*w)**2+576._REKIND*            &
    & SIN(t1*w)**2-78._REKIND*dt**2*w**2*SIN(t1*w)**2+339._REKIND*dt*w*        &
    & SIN(2._REKIND*t1*w)-6._REKIND*dt**3*w**3*SIN(2._REKIND*t1*w)+ 96._REKIND &
    & *dt*w*COS((dt+t1)*w)*(5._REKIND*SIN(t1*w)-4._REKIND*                     &
    & SIN((dt/2._REKIND+t1)*w))+48*dt*w*(-16 + dt**2*w**2)*COS(t1*w)*          &
    & SIN((dt/2._REKIND+t1)*w)-1152._REKIND*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w) &
    & + 240._REKIND*dt**2*w**2*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w)-8._REKIND    &
    & *dt**4*w**4*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w)-1152._REKIND*SIN(t1*w)*   &
    & SIN((dt+t1)*w)+48._REKIND*dt**2*w**2*SIN(t1*w)*SIN((dt+t1)*w)+           &
    & 1152._REKIND*SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w)-48._REKIND*         &
    & dt**2*w**2*SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w)-3._REKIND*dt*w*       &
    & SIN(2._REKIND*(dt+t1)*w))/(12._REKIND*dt**4*w**4)    
k(4,6) = (-1152._REKIND+72._REKIND*dt**2*w**2+2._REKIND*dt**4*w**4+2._REKIND   &
    & *(576._REKIND-84._REKIND*dt**2*w**2+dt**4*w**4)*COS(dt*w)-6._REKIND*     &
    & (-96._REKIND+7._REKIND*dt**2*w**2)*COS(2._REKIND*t1*w)+576._REKIND*      &
    & COS(2._REKIND*(dt+t1)*w)-42._REKIND*dt**2*w**2*COS(2._REKIND*(dt+t1)*w)  &
    & -1152._REKIND*COS((dt+2._REKIND*t1)*w)+168._REKIND*dt**2*w**2*           &
    & COS((dt+2._REKIND*t1)*w)-2._REKIND*dt**4*w**4*COS((dt+2._REKIND*t1)*w)+  &
    & 672._REKIND*dt*w*SIN(dt*w)-24._REKIND*dt**3*w**3*SIN(dt*w)-246._REKIND   &
    & *dt*w*SIN(2._REKIND*t1*w)+3._REKIND*dt**3*w**3*SIN(2._REKIND*t1*w)+      &
    & 246._REKIND*dt*w*SIN(2._REKIND*(dt+t1)*w)-3._REKIND*dt**3*w**3*          &
    & SIN(2._REKIND*(dt+t1)*w))/(24._REKIND*dt**4*w**4)

k(5,1) = -k(1,5)
k(5,2) = -k(2,5)
k(5,3) = -k(3,5)
k(5,4) = -k(4,5)
k(5,5) = 0._REKIND
k(5,6) = -(-576._REKIND+72._REKIND*dt**2*w**2+4._REKIND*dt**4*w**4+            &
    & (576._REKIND-78._REKIND*dt**2*w**2)*COS((dt+t1)*w)**2-3._REKIND*dt*w*    &
    & SIN(2._REKIND*t1*w)-384._REKIND*dt*w*COS(t1*w)*SIN((dt/2._REKIND+t1)*w)+ &
    & 48._REKIND*dt*w*(-16._REKIND+dt**2*w**2)*COS((dt+t1)*w)*SIN((dt/         &
    & 2._REKIND+t1)*w)-1152._REKIND*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w)+        &
    & 48._REKIND*dt**2*w**2*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w)+480._REKIND     &
    & *dt*w*COS(t1*w)*SIN((dt+t1)*w)+1152._REKIND*SIN(t1*w)*SIN((dt+t1)*w)     &
    & -48._REKIND*dt**2*w**2*SIN(t1*w)*SIN((dt+t1)*w)+1152._REKIND*            &
    & SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w)-240._REKIND*dt**2*w**2*          &
    & SIN((dt/2._REKIND+ t1)*w)*SIN((dt+t1)*w)+8._REKIND*dt**4*w**4*           &
    & SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w)-576._REKIND*SIN((dt+t1)*w)**2+   &
    & 78._REKIND*dt**2*w**2*SIN((dt+t1)*w)**2+339._REKIND*dt*w*SIN(2._REKIND   &
    & *(dt+t1)*w)-6._REKIND*dt**3*w**3*SIN(2._REKIND*(dt+t1)*w))/(12._REKIND   &
    & *dt**4*w**4)

k(6,1) = -k(1,6)
k(6,2) = -k(2,6)
k(6,3) = -k(3,6)
k(6,4) = -k(4,6)
k(6,5) = -k(5,6)
k(6,6) = 0._REKIND

! Coefficients for Mass matrix of the complete slab
m(1,1) = -4._REKIND/dt**2
m(1,2) = 8._REKIND/dt**2
m(1,3) = -4._REKIND/dt**2
m(1,4) = (3._REKIND*dt*w*COS(t1*w)-SIN(t1*w)+SIN((dt+t1)*w))/dt**2
m(1,5) = (4._REKIND*(3._REKIND*SIN(t1*w)-2._REKIND*SIN((dt/2._REKIND+t1)*w)-   &
    & SIN((dt+t1)*w)))/dt**2
m(1,6) = (dt*w*COS((dt+t1)*w)-3._REKIND*SIN(t1*w)+3._REKIND*SIN((dt+t1)*w))    &
    & /dt**2
  
m(2,1) = 0._REKIND
m(2,2) = 0._REKIND
m(2,3) = 0._REKIND
m(2,4) = (-4._REKIND*(dt*w*COS(t1*w)-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
m(2,5) = (16._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
m(2,6) = (-4._REKIND*(dt*w*COS((dt+t1)*w)-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2

m(3,1) = 4._REKIND/dt**2
m(3,2) = -8._REKIND/dt**2
m(3,3) = 4._REKIND/dt**2
m(3,4) = (dt*w*COS(t1*w)-3._REKIND*SIN(t1*w)+3._REKIND*SIN((dt+t1)*w))/dt**2
m(3,5) = (4._REKIND*(SIN(t1*w)+2._REKIND*SIN((dt/2._REKIND+t1)*w)-3._REKIND*   &
    & SIN((dt+t1)*w)))/dt**2
m(3,6) = (3._REKIND*dt*w*COS((dt+t1)*w)-SIN(t1*w)+SIN((dt+t1)*w))/dt**2

m(4,1) = 0._REKIND
m(4,2) = 0._REKIND
m(4,3) = 0._REKIND
m(4,4) = -(-2._REKIND+dt**2*w**2+2._REKIND*COS(dt*w)+(1._REKIND+dt**2*w**2)*   &
    & COS(2._REKIND*t1*w)+COS(2._REKIND*(dt+t1)*w)-2._REKIND*COS((dt+2._REKIND &
    & *t1)*w))/(4._REKIND*dt**2)
m(4,5) = (-8._REKIND*dt*w+4._REKIND*dt**3*w**3+3._REKIND*SIN(2._REKIND*t1*w)+  &
    & 2._REKIND*dt*w*(4._REKIND*COS(dt*w)+5._REKIND*COS(2._REKIND*t1*w)+       &
    & 2._REKIND*COS(2._REKIND*(dt+t1)*w)-4._REKIND*COS((dt+2._REKIND*t1)*w)-   &
    & 3._REKIND*dt*w*SIN(2._REKIND*t1*w))-3._REKIND*SIN(2._REKIND*(dt+t1)*w))/ &
    & (4._REKIND*dt**3*w)
m(4,6) = (-16._REKIND*dt*w-2._REKIND*dt**3*w**3+8._REKIND*dt*w*COS(t1*w)**2-   &
    & 4._REKIND*dt*w*COS((dt+t1)*w)**2-8._REKIND*dt**2*w**2*COS((dt+t1)*w)*    &
    & SIN(t1*w)+16._REKIND*dt*w*SIN(t1*w)**2-6._REKIND*SIN(2._REKIND*t1*w)+    &
    & 3._REKIND*dt**2*w**2*SIN(2._REKIND*t1*w)+8._REKIND*dt*w*SIN(t1*w)*       &
    & SIN((dt+t1)*w)+4._REKIND*dt*w*SIN((dt+t1)*w)**2+6._REKIND*SIN(2._REKIND* &
    & (dt+t1)*w)+dt**2*w**2*SIN(2._REKIND*(dt+t1)*w))/(8._REKIND*dt**3*w)

m(5,1) = 0._REKIND
m(5,2) = 0._REKIND
m(5,3) = 0._REKIND	
m(5,4) = -(4._REKIND*dt**3*w**3+10._REKIND*dt*w*COS(t1*w)**2-4._REKIND*dt*w*   & 
    & COS((dt+t1)*w)**2-10._REKIND*dt*w*SIN(t1*w)**2+3._REKIND*SIN(2._REKIND   &
    & *t1*w)+2._REKIND*dt**2*w**2*SIN(2._REKIND*t1*w)-16._REKIND*dt**2*w**2    &
    & *COS(t1*w)*SIN((dt/2._REKIND+t1)*w)+16._REKIND*dt*w*SIN(t1*w)*SIN((dt/   &
    & 2._REKIND+t1)*w)-16._REKIND*dt*w*SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w) &
    & +4._REKIND*dt*w*SIN((dt+t1)*w)**2-3._REKIND*SIN(2._REKIND*(dt+t1)*w))/   &
    & (4._REKIND*dt**3*w)
m(5,5) = (4._REKIND*(COS(2._REKIND*t1*w)-COS(2._REKIND*(dt+t1)*w)+2._REKIND*   &
    & COS(((3*dt)/2._REKIND+2._REKIND*t1)*w)-2._REKIND*COS(((dt+4._REKIND*t1)* &
    & w)/2._REKIND)))/dt**2
m(5,6)=(4._REKIND*dt**3*w**3-4._REKIND*dt*w*COS(t1*w)**2+10._REKIND*dt*w*      &
    & COS((dt+t1)*w)**2+4._REKIND*dt*w*SIN(t1*w)**2+3._REKIND*SIN(2._REKIND*   &
    & t1*w)+16._REKIND*dt**2*w**2*COS((dt+t1)*w)*SIN((dt/2._REKIND+t1)*w)-     &
    & 16._REKIND*dt*w*SIN(t1*w)*SIN((dt/2._REKIND+t1)*w)+16._REKIND*dt*w*      &
    & SIN((dt/2._REKIND+t1)*w)*SIN((dt+t1)*w)-10._REKIND*dt*w*SIN((dt+t1)*w)   &
    & **2-3._REKIND*SIN(2._REKIND*(dt+t1)*w)-2._REKIND*dt**2*w**2*             &
    & SIN(2._REKIND*(dt+t1)*w))/(4._REKIND*dt**3*w)

m(6,1) = 0._REKIND
m(6,2) = 0._REKIND
m(6,3) = 0._REKIND
m(6,4) = (2._REKIND*dt*w+2._REKIND*dt**3*w**3-4._REKIND*dt*w*COS(dt*w)+        &
    & 4._REKIND*dt*w*COS(t1*w)**2+2._REKIND*dt*w*COS(2._REKIND*t1*w)+          &
    & 4._REKIND*dt*w*COS(2._REKIND*(dt+t1)*w)+4._REKIND*dt*w*                  &
    & COS((dt+2._REKIND*t1)*w)+6._REKIND*SIN(2._REKIND*t1*w)+dt**2*w**2*       &
    & SIN(2._REKIND*t1*w)-8._REKIND*dt**2*w**2*COS(t1*w)*SIN((dt+t1)*w)-       &
    & 6._REKIND*SIN(2._REKIND*(dt+t1)*w)+3._REKIND*dt**2*w**2*                 &
    & SIN(2._REKIND*(dt+t1)*w))/(8._REKIND*dt**3*w)
m(6,5) = (8._REKIND*dt*w-4._REKIND*dt**3*w**3-3._REKIND*SIN(2._REKIND*t1*w)+   &
    & 3._REKIND*SIN(2._REKIND*(dt+t1)*w)-2._REKIND*dt*w*(4._REKIND*COS(dt*w)+  &
    & 2._REKIND*COS(2*t1*w)+5._REKIND*COS(2._REKIND*(dt+t1)*w)-4._REKIND*      &
    & COS((dt+2._REKIND*t1)*w)+3._REKIND*dt*w*SIN(2._REKIND*(dt+t1)*w)))/      &
    & (4._REKIND*dt**3*w)
m(6,6) = (-2._REKIND+dt**2*w**2+2._REKIND*COS(dt*w)+COS(2._REKIND*t1*w)+       &
    & (1._REKIND+dt**2*w**2)*COS(2._REKIND*(dt+t1)*w)-2._REKIND*COS((dt+       &
    & 2._REKIND*t1)*w))/(4._REKIND*dt**2)

! Coefficients for the Stiffness matrix at tn_1 (+)
DO i = 1, 6
    DO j = 1, 6
        ktn(i,j) = 0._REKIND
    ENDDO
ENDDO
ktn(1,1) = 1._REKIND

! Coefficients for the Mass matrix at tn_1 (+) 
mtn(1,1) = 9._REKIND/dt**2
mtn(1,2) = -12._REKIND/dt**2
mtn(1,3) = 3._REKIND/dt**2
mtn(1,4) = (-3._REKIND*w*COS(t1*w))/dt
mtn(1,5) = (12._REKIND*(-SIN(t1*w)+SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn(1,6) = (3._REKIND*(SIN(t1*w)-SIN((dt+t1)*w)))/dt**2

mtn(2,1) = mtn(1,2)
mtn(2,2) = 16._REKIND/dt**2
mtn(2,3) = -4._REKIND/dt**2
mtn(2,4) = (4._REKIND*w*COS(t1*w))/dt
mtn(2,5) = (16._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn(2,6) = (4._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2

mtn(3,1) = mtn(1,3)
mtn(3,2) = mtn(2,3)
mtn(3,3) = 1._REKIND/dt**2
mtn(3,4) = -((w*COS(t1*w))/dt)
mtn(3,5) = (4._REKIND*(-SIN(t1*w)+SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn(3,6) = (SIN(t1*w)-SIN((dt+t1)*w))/dt**2

mtn(4,1) = mtn(1,4)
mtn(4,2) = mtn(2,4)
mtn(4,3) = mtn(3,4)
mtn(4,4) = w**2*COS(t1*w)**2
mtn(4,5) = (4._REKIND*w*COS(t1*w)*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))/dt
mtn(4,6) = (w*COS(t1*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt

mtn(5,1) = mtn(1,5)
mtn(5,2) = mtn(2,5)
mtn(5,3) = mtn(3,5)
mtn(5,4) = mtn(4,5)
mtn(5,5) = (16._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w))**2)/dt**2
mtn(5,6) = (4._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w))*(-SIN(t1*w)+        &
    & SIN((dt+t1)*w)))/dt**2

mtn(6,1) = mtn(1,6)
mtn(6,2) = mtn(2,6)
mtn(6,3) = mtn(3,6)
mtn(6,4) = mtn(4,6)
mtn(6,5) = mtn(5,6)
mtn(6,6) = (SIN(t1*w)-SIN((dt+t1)*w))**2/dt**2

! Coefficients for the Stiffness matrix at tn_1 (-) 
DO i=1,6
    DO j=1,6
        ktn1(i,j)=0._REKIND
    ENDDO
ENDDO
ktn1(1,3) = 1._REKIND

! Coefficients for the Mass matrix at tn_1 (-)
mtn1(1,1) = -3._REKIND/dt**2
mtn1(1,2) = 12._REKIND/dt**2
mtn1(1,3) = -9/dt**2
mtn1(1,4) = (-3._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(1,5) = (-12._REKIND*(SIN((dt/2._REKIND+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(1,6) = (-3._REKIND*w*COS((dt+t1)*w))/dt

mtn1(2,1) = 4._REKIND/dt**2
mtn1(2,2) = -16._REKIND/dt**2
mtn1(2,3) = 12._REKIND/dt**2
mtn1(2,4) = (4._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(2,5) = (16._REKIND*(SIN((dt/2._REKIND+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(2,6) = (4._REKIND*w*COS((dt+t1)*w))/dt

mtn1(3,1) = -1._REKIND/dt**2
mtn1(3,2) = 4._REKIND/dt**2
mtn1(3,3) = -3._REKIND/dt**2
mtn1(3,4) = -((-SIN(t1*w)+SIN((dt+t1)*w))/dt**2)
mtn1(3,5) = (-4._REKIND*(SIN((dt/2._REKIND+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(3,6) = -((w*COS((dt+t1)*w))/dt)

mtn1(4,1) = (w*COS(t1*w))/dt
mtn1(4,2) = (-4._REKIND*w*COS(t1*w))/dt
mtn1(4,3) = (3._REKIND*w*COS(t1*w))/dt
mtn1(4,4) = (w*COS(t1*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt
mtn1(4,5) = (4._REKIND*w*COS(t1*w)*(SIN((dt/2._REKIND+t1)*w)-                  &
    & SIN((dt+t1)*w)))/dt
mtn1(4,6) = w**2*COS(t1*w)*COS((dt+t1)*w)

mtn1(5,1) = (4._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn1(5,2) = (-16._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn1(5,3) = (12._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))/dt**2
mtn1(5,4) = (4._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w))*(-SIN(t1*w)+       &
    & SIN((dt+t1)*w)))/dt**2
mtn1(5,5) = (16._REKIND*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w))*(SIN((dt/         &
    & 2._REKIND+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(5,6) = (4._REKIND*w*COS((dt+t1)*w)*(SIN(t1*w)-SIN((dt/2._REKIND+t1)*w)))  &
    & /dt

mtn1(6,1) = (-SIN(t1*w)+SIN((dt+t1)*w))/dt**2
mtn1(6,2) = (-4._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(6,3) = (3._REKIND*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(6,4) = (-SIN(t1*w)+SIN((dt+t1)*w))**2/dt**2
mtn1(6,5) = (4._REKIND*(SIN((dt/2._REKIND+t1)*w)-SIN((dt+t1)*w))*(-SIN(t1*w)+  &
    & SIN((dt+t1)*w)))/dt**2
mtn1(6,6) = (w*COS((dt+t1)*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt
	
RETURN
END SUBROUTINE stcoeff
