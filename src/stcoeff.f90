SUBROUTINE stcoeff_t(t1, dt, w, k, m, ktn, mtn, ktn1, mtn1)
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
REAL(KIND=REKIND), DIMENSION(3, 3), INTENT(OUT) :: k, m, ktn, mtn, ktn1, mtn1
! ---Internal variables---
INTEGER :: i, j

! Coefficients for Stiffness matrix of the complete slab
k(1,1) = 2./15.*dt
k(1,2) = 1./15.*dt
k(1,3) = -1./30.*dt
	
k(2,1) = 1./15.*dt
k(2,2) = 8./15.*dt
k(2,3) = 1./15.*dt
     
k(3,1) = -1./30.*dt
k(3,2) = 1./15.*dt
k(3,3) = 2./15.*dt

! Coefficients for Mass matrix of the complete slab
m(1,1) = -.5
m(1,2) = 2./3.
m(1,3) = -1./6.
  
m(2,1) = -2./3.
m(2,2) = 0.
m(2,3) = 2./3.

m(3,1) = 1./6.
m(3,2) = -2./3.
m(3,3) = .5

! Coefficients for the Stiffness matrix at tn_1 (+)
ktn = 0.

! Coefficients for the Mass matrix at tn_1 (+)
mtn = 0.
mtn(1,1) = 1.

! Coefficients for the Stiffness matrix at tn_1 (-) 
ktn1 = 0.

! Coefficients for the Mass matrix at tn_1 (-)
mtn1 = 0.
mtn1(1,3) = 1.
	
RETURN
END SUBROUTINE stcoeff_t

SUBROUTINE stcoeff_m(t1, dt, w, k, m, ktn, mtn, ktn1, mtn1)
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
k(1,1) = -0.5
k(1,2) = -2./3.
k(1,3) = 1./6.
k(1,4) = (-6.*dt*w*(-12.+dt**2*w**2)*COS(t1*w)+24.*dt*w*  &
    & COS((dt+t1)*w)+(96.-26.*dt**2*w**2+dt**4*w**4)*SIN(t1*w)+  &
    & 2.*(-48.+dt**2*w**2)*SIN((dt+t1)*w))/(2.*dt**4*w**4)  
k(1,5) = (2.*(dt*w*(-84.*COS(t1*w)-60.*COS((dt+t1)*w)+    &
    & dt*w*(18.*SIN(t1*w)+dt**2*w**2*SIN((dt/2.+t1)*w)-          &
    & 6.*SIN((dt+t1)*w)))-144.*(SIN(t1*w)-SIN((dt+t1)*w))))/     &
    & (3.*dt**4*w**4)
k(1,6) = (120.*dt*w*COS(t1*w)-6.*dt*w*(-28.+dt**2*w**2)*  &
    & COS((dt+t1)*w)-18.*(-16.+dt**2*w**2)*SIN(t1*w)-            &
    & (288.-42.*dt**2*w**2+dt**4*w**4)*SIN((dt+t1)*w))/          &
    & (6.*dt**4*w**4)
	
k(2,1) = 2./3.
k(2,2) = 0.0
k(2,3) = -2./3.
k(2,4) = (-2.*(-6.*dt*w*(-16.+dt**2*w**2)*COS(t1*w)+      &
    & 48.*dt*w*COS((dt+t1)*w)+(-24.+dt**2*w**2)*((-6.+    &
    & dt**2*w**2)*SIN(t1*w)+6.*SIN((dt+t1)*w))))/(3.*dt**4*w**4)
k(2,5) = (16.*(6.*dt*w*COS(t1*w)+6.*dt*w*COS((dt+t1)*w)-  &
    & (-12.+dt**2*w**2)*(SIN(t1*w)-SIN((dt+t1)*w))))/(dt**4*w**4)
k(2,6) = (2.*(-48.*dt*w*COS(t1*w)+6.*dt*w*(-16.    &
    & +dt**2*w**2) *COS((dt+t1)*w)+(-24.+dt**2*w**2)*(6.*        &
    & SIN(t1*w)+(-6.+dt**2*w**2)*SIN((dt+t1)*w))))/(3.*dt**4*    &
    & w**4)
     
k(3,1) = -1./6.
k(3,2) = 2./3.
k(3,3) = 1./2.
k(3,4) = (-6.*dt*w*(-28.+dt**2*w**2)*COS(t1*w)+120.*dt*w* &
    & COS((dt+t1)*w)+(288.-42.*dt**2*w**2+dt**4*w**4)*SIN(t1*w)+ &
    & 18.*(-16.+dt**2*w**2)*SIN((dt+t1)*w))/(6.*dt**4     &
    & *w**4)
k(3,5) = (-2.*(144.*SIN(t1*w)+dt*w*(60.*COS(t1*w)+        &
    & 84.*COS((dt+t1)*w)-6.*dt*w*SIN(t1*w)+dt**3*w**3*SIN((dt    &
    & /2.+t1)*w))+18.*(-8.+dt**2*w**2)*SIN((dt+t1)*w)))/  &
    & (3.*dt**4*w**4) 
k(3,6) = (96.*(SIN(t1*w)-SIN((dt+t1)*w))+dt*w*(24.*COS(t1*w)+    &
    & (72.-6.*dt**2*w**2)*COS((dt+t1)*w)+dt*w*(-2.*       &
    & SIN(t1*w)+(26.-dt**2*w**2)*SIN((dt+t1)*w))))/(2.*dt**4*    &
    & w**4)
k(4,1) = -k(1,4)
k(4,2) = -k(2,4)
k(4,3) = -k(3,4)
k(4,4) = 0.
k(4,5) = (576.-72.*dt**2*w**2-4.*dt**4*w**4+6.*    &
    & (-96.+13.*dt**2*w**2)*COS(t1*w)**2+576.*            &
    & SIN(t1*w)**2-78.*dt**2*w**2*SIN(t1*w)**2+339.*dt*w*        &
    & SIN(2.*t1*w)-6.*dt**3*w**3*SIN(2.*t1*w)+ 96. &
    & *dt*w*COS((dt+t1)*w)*(5.*SIN(t1*w)-4.*                     &
    & SIN((dt/2.+t1)*w))+48*dt*w*(-16 + dt**2*w**2)*COS(t1*w)*          &
    & SIN((dt/2.+t1)*w)-1152.*SIN(t1*w)*SIN((dt/2.+t1)*w) &
    & + 240.*dt**2*w**2*SIN(t1*w)*SIN((dt/2.+t1)*w)-8.    &
    & *dt**4*w**4*SIN(t1*w)*SIN((dt/2.+t1)*w)-1152.*SIN(t1*w)*   &
    & SIN((dt+t1)*w)+48.*dt**2*w**2*SIN(t1*w)*SIN((dt+t1)*w)+           &
    & 1152.*SIN((dt/2.+t1)*w)*SIN((dt+t1)*w)-48.*         &
    & dt**2*w**2*SIN((dt/2.+t1)*w)*SIN((dt+t1)*w)-3.*dt*w*       &
    & SIN(2.*(dt+t1)*w))/(12.*dt**4*w**4)    
k(4,6) = (-1152.+72.*dt**2*w**2+2.*dt**4*w**4+2.   &
    & *(576.-84.*dt**2*w**2+dt**4*w**4)*COS(dt*w)-6.*     &
    & (-96.+7.*dt**2*w**2)*COS(2.*t1*w)+576.*      &
    & COS(2.*(dt+t1)*w)-42.*dt**2*w**2*COS(2.*(dt+t1)*w)  &
    & -1152.*COS((dt+2.*t1)*w)+168.*dt**2*w**2*           &
    & COS((dt+2.*t1)*w)-2.*dt**4*w**4*COS((dt+2.*t1)*w)+  &
    & 672.*dt*w*SIN(dt*w)-24.*dt**3*w**3*SIN(dt*w)-246.   &
    & *dt*w*SIN(2.*t1*w)+3.*dt**3*w**3*SIN(2.*t1*w)+      &
    & 246.*dt*w*SIN(2.*(dt+t1)*w)-3.*dt**3*w**3*          &
    & SIN(2.*(dt+t1)*w))/(24.*dt**4*w**4)

k(5,1) = -k(1,5)
k(5,2) = -k(2,5)
k(5,3) = -k(3,5)
k(5,4) = -k(4,5)
k(5,5) = 0.
k(5,6) = -(-576.+72.*dt**2*w**2+4.*dt**4*w**4+            &
    & (576.-78.*dt**2*w**2)*COS((dt+t1)*w)**2-3.*dt*w*    &
    & SIN(2.*t1*w)-384.*dt*w*COS(t1*w)*SIN((dt/2.+t1)*w)+ &
    & 48.*dt*w*(-16.+dt**2*w**2)*COS((dt+t1)*w)*SIN((dt/         &
    & 2.+t1)*w)-1152.*SIN(t1*w)*SIN((dt/2.+t1)*w)+        &
    & 48.*dt**2*w**2*SIN(t1*w)*SIN((dt/2.+t1)*w)+480.     &
    & *dt*w*COS(t1*w)*SIN((dt+t1)*w)+1152.*SIN(t1*w)*SIN((dt+t1)*w)     &
    & -48.*dt**2*w**2*SIN(t1*w)*SIN((dt+t1)*w)+1152.*            &
    & SIN((dt/2.+t1)*w)*SIN((dt+t1)*w)-240.*dt**2*w**2*          &
    & SIN((dt/2.+ t1)*w)*SIN((dt+t1)*w)+8.*dt**4*w**4*           &
    & SIN((dt/2.+t1)*w)*SIN((dt+t1)*w)-576.*SIN((dt+t1)*w)**2+   &
    & 78.*dt**2*w**2*SIN((dt+t1)*w)**2+339.*dt*w*SIN(2.   &
    & *(dt+t1)*w)-6.*dt**3*w**3*SIN(2.*(dt+t1)*w))/(12.   &
    & *dt**4*w**4)

k(6,1) = -k(1,6)
k(6,2) = -k(2,6)
k(6,3) = -k(3,6)
k(6,4) = -k(4,6)
k(6,5) = -k(5,6)
k(6,6) = 0.
! Coefficients for Mass matrix of the complete slab
m(1,1) = -4./dt**2
m(1,2) = 8./dt**2
m(1,3) = -4./dt**2
m(1,4) = (3.*dt*w*COS(t1*w)-SIN(t1*w)+SIN((dt+t1)*w))/dt**2
m(1,5) = (4.*(3.*SIN(t1*w)-2.*SIN((dt/2.+t1)*w)-   &
    & SIN((dt+t1)*w)))/dt**2
m(1,6) = (dt*w*COS((dt+t1)*w)-3.*SIN(t1*w)+3.*SIN((dt+t1)*w))    &
    & /dt**2
  
m(2,1) = 0.
m(2,2) = 0.
m(2,3) = 0.
m(2,4) = (-4.*(dt*w*COS(t1*w)-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
m(2,5) = (16.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
m(2,6) = (-4.*(dt*w*COS((dt+t1)*w)-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2

m(3,1) = 4./dt**2
m(3,2) = -8./dt**2
m(3,3) = 4./dt**2
m(3,4) = (dt*w*COS(t1*w)-3.*SIN(t1*w)+3.*SIN((dt+t1)*w))/dt**2
m(3,5) = (4.*(SIN(t1*w)+2.*SIN((dt/2.+t1)*w)-3.*   &
    & SIN((dt+t1)*w)))/dt**2
m(3,6) = (3.*dt*w*COS((dt+t1)*w)-SIN(t1*w)+SIN((dt+t1)*w))/dt**2
m(4,1) = 0.
m(4,2) = 0.
m(4,3) = 0.
m(4,4) = -(-2.+dt**2*w**2+2.*COS(dt*w)+(1.+dt**2*w**2)*   &
    & COS(2.*t1*w)+COS(2.*(dt+t1)*w)-2.*COS((dt+2. &
    & *t1)*w))/(4.*dt**2)
m(4,5) = (-8.*dt*w+4.*dt**3*w**3+3.*SIN(2.*t1*w)+  &
    & 2.*dt*w*(4.*COS(dt*w)+5.*COS(2.*t1*w)+       &
    & 2.*COS(2.*(dt+t1)*w)-4.*COS((dt+2.*t1)*w)-   &
    & 3.*dt*w*SIN(2.*t1*w))-3.*SIN(2.*(dt+t1)*w))/ &
    & (4.*dt**3*w)
m(4,6) = (-16.*dt*w-2.*dt**3*w**3+8.*dt*w*COS(t1*w)**2-   &
    & 4.*dt*w*COS((dt+t1)*w)**2-8.*dt**2*w**2*COS((dt+t1)*w)*    &
    & SIN(t1*w)+16.*dt*w*SIN(t1*w)**2-6.*SIN(2.*t1*w)+    &
    & 3.*dt**2*w**2*SIN(2.*t1*w)+8.*dt*w*SIN(t1*w)*       &
    & SIN((dt+t1)*w)+4.*dt*w*SIN((dt+t1)*w)**2+6.*SIN(2.* &
    & (dt+t1)*w)+dt**2*w**2*SIN(2.*(dt+t1)*w))/(8.*dt**3*w)

m(5,1) = 0.
m(5,2) = 0.
m(5,3) = 0.	
m(5,4) = -(4.*dt**3*w**3+10.*dt*w*COS(t1*w)**2-4.*dt*w*   & 
    & COS((dt+t1)*w)**2-10.*dt*w*SIN(t1*w)**2+3.*SIN(2.   &
    & *t1*w)+2.*dt**2*w**2*SIN(2.*t1*w)-16.*dt**2*w**2    &
    & *COS(t1*w)*SIN((dt/2.+t1)*w)+16.*dt*w*SIN(t1*w)*SIN((dt/   &
    & 2.+t1)*w)-16.*dt*w*SIN((dt/2.+t1)*w)*SIN((dt+t1)*w) &
    & +4.*dt*w*SIN((dt+t1)*w)**2-3.*SIN(2.*(dt+t1)*w))/   &
    & (4.*dt**3*w)
m(5,5) = (4.*(COS(2.*t1*w)-COS(2.*(dt+t1)*w)+2.*   &
    & COS(((3*dt)/2.+2.*t1)*w)-2.*COS(((dt+4.*t1)* &
    & w)/2.)))/dt**2
m(5,6)=(4.*dt**3*w**3-4.*dt*w*COS(t1*w)**2+10.*dt*w*      &
    & COS((dt+t1)*w)**2+4.*dt*w*SIN(t1*w)**2+3.*SIN(2.*   &
    & t1*w)+16.*dt**2*w**2*COS((dt+t1)*w)*SIN((dt/2.+t1)*w)-     &
    & 16.*dt*w*SIN(t1*w)*SIN((dt/2.+t1)*w)+16.*dt*w*      &
    & SIN((dt/2.+t1)*w)*SIN((dt+t1)*w)-10.*dt*w*SIN((dt+t1)*w)   &
    & **2-3.*SIN(2.*(dt+t1)*w)-2.*dt**2*w**2*             &
    & SIN(2.*(dt+t1)*w))/(4.*dt**3*w)

m(6,1) = 0.
m(6,2) = 0.
m(6,3) = 0.
m(6,4) = (2.*dt*w+2.*dt**3*w**3-4.*dt*w*COS(dt*w)+        &
    & 4.*dt*w*COS(t1*w)**2+2.*dt*w*COS(2.*t1*w)+          &
    & 4.*dt*w*COS(2.*(dt+t1)*w)+4.*dt*w*                  &
    & COS((dt+2.*t1)*w)+6.*SIN(2.*t1*w)+dt**2*w**2*       &
    & SIN(2.*t1*w)-8.*dt**2*w**2*COS(t1*w)*SIN((dt+t1)*w)-       &
    & 6.*SIN(2.*(dt+t1)*w)+3.*dt**2*w**2*                 &
    & SIN(2.*(dt+t1)*w))/(8.*dt**3*w)
m(6,5) = (8.*dt*w-4.*dt**3*w**3-3.*SIN(2.*t1*w)+   &
    & 3.*SIN(2.*(dt+t1)*w)-2.*dt*w*(4.*COS(dt*w)+  &
    & 2.*COS(2*t1*w)+5.*COS(2.*(dt+t1)*w)-4.*      &
    & COS((dt+2.*t1)*w)+3.*dt*w*SIN(2.*(dt+t1)*w)))/      &
    & (4.*dt**3*w)
m(6,6) = (-2.+dt**2*w**2+2.*COS(dt*w)+COS(2.*t1*w)+       &
    & (1.+dt**2*w**2)*COS(2.*(dt+t1)*w)-2.*COS((dt+       &
    & 2.*t1)*w))/(4.*dt**2)
! Coefficients for the Stiffness matrix at tn_1 (+)
ktn = 0.
ktn(1,1) = 1.

! Coefficients for the Mass matrix at tn_1 (+) 
mtn(1,1) = 9./dt**2
mtn(1,2) = -12./dt**2
mtn(1,3) = 3./dt**2
mtn(1,4) = (-3.*w*COS(t1*w))/dt
mtn(1,5) = (12.*(-SIN(t1*w)+SIN((dt/2.+t1)*w)))/dt**2
mtn(1,6) = (3.*(SIN(t1*w)-SIN((dt+t1)*w)))/dt**2

mtn(2,1) = mtn(1,2)
mtn(2,2) = 16./dt**2
mtn(2,3) = -4./dt**2
mtn(2,4) = (4.*w*COS(t1*w))/dt
mtn(2,5) = (16.*(SIN(t1*w)-SIN((dt/2.+t1)*w)))/dt**2
mtn(2,6) = (4.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2

mtn(3,1) = mtn(1,3)
mtn(3,2) = mtn(2,3)
mtn(3,3) = 1./dt**2
mtn(3,4) = -((w*COS(t1*w))/dt)
mtn(3,5) = (4.*(-SIN(t1*w)+SIN((dt/2.+t1)*w)))/dt**2
mtn(3,6) = (SIN(t1*w)-SIN((dt+t1)*w))/dt**2
mtn(4,1) = mtn(1,4)
mtn(4,2) = mtn(2,4)
mtn(4,3) = mtn(3,4)
mtn(4,4) = w**2*COS(t1*w)**2
mtn(4,5) = (4.*w*COS(t1*w)*(SIN(t1*w)-SIN((dt/2.+t1)*w)))/dt
mtn(4,6) = (w*COS(t1*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt

mtn(5,1) = mtn(1,5)
mtn(5,2) = mtn(2,5)
mtn(5,3) = mtn(3,5)
mtn(5,4) = mtn(4,5)
mtn(5,5) = (16.*(SIN(t1*w)-SIN((dt/2.+t1)*w))**2)/dt**2
mtn(5,6) = (4.*(SIN(t1*w)-SIN((dt/2.+t1)*w))*(-SIN(t1*w)+        &
    & SIN((dt+t1)*w)))/dt**2

mtn(6,1) = mtn(1,6)
mtn(6,2) = mtn(2,6)
mtn(6,3) = mtn(3,6)
mtn(6,4) = mtn(4,6)
mtn(6,5) = mtn(5,6)
mtn(6,6) = (SIN(t1*w)-SIN((dt+t1)*w))**2/dt**2
! Coefficients for the Stiffness matrix at tn_1 (-) 
ktn1=0.
ktn1(1,3) = 1.

! Coefficients for the Mass matrix at tn_1 (-)
mtn1(1,1) = -3./dt**2
mtn1(1,2) = 12./dt**2
mtn1(1,3) = -9/dt**2
mtn1(1,4) = (-3.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(1,5) = (-12.*(SIN((dt/2.+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(1,6) = (-3.*w*COS((dt+t1)*w))/dt

mtn1(2,1) = 4./dt**2
mtn1(2,2) = -16./dt**2
mtn1(2,3) = 12./dt**2
mtn1(2,4) = (4.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(2,5) = (16.*(SIN((dt/2.+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(2,6) = (4.*w*COS((dt+t1)*w))/dt

mtn1(3,1) = -1./dt**2
mtn1(3,2) = 4./dt**2
mtn1(3,3) = -3./dt**2
mtn1(3,4) = -((-SIN(t1*w)+SIN((dt+t1)*w))/dt**2)
mtn1(3,5) = (-4.*(SIN((dt/2.+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(3,6) = -((w*COS((dt+t1)*w))/dt)
mtn1(4,1) = (w*COS(t1*w))/dt
mtn1(4,2) = (-4.*w*COS(t1*w))/dt
mtn1(4,3) = (3.*w*COS(t1*w))/dt
mtn1(4,4) = (w*COS(t1*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt
mtn1(4,5) = (4.*w*COS(t1*w)*(SIN((dt/2.+t1)*w)-                  &
    & SIN((dt+t1)*w)))/dt
mtn1(4,6) = w**2*COS(t1*w)*COS((dt+t1)*w)

mtn1(5,1) = (4.*(SIN(t1*w)-SIN((dt/2.+t1)*w)))/dt**2
mtn1(5,2) = (-16.*(SIN(t1*w)-SIN((dt/2.+t1)*w)))/dt**2
mtn1(5,3) = (12.*(SIN(t1*w)-SIN((dt/2.+t1)*w)))/dt**2
mtn1(5,4) = (4.*(SIN(t1*w)-SIN((dt/2.+t1)*w))*(-SIN(t1*w)+       &
    & SIN((dt+t1)*w)))/dt**2
mtn1(5,5) = (16.*(SIN(t1*w)-SIN((dt/2.+t1)*w))*(SIN((dt/         &
    & 2.+t1)*w)-SIN((dt+t1)*w)))/dt**2
mtn1(5,6) = (4.*w*COS((dt+t1)*w)*(SIN(t1*w)-SIN((dt/2.+t1)*w)))  &
    & /dt

mtn1(6,1) = (-SIN(t1*w)+SIN((dt+t1)*w))/dt**2
mtn1(6,2) = (-4.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(6,3) = (3.*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt**2
mtn1(6,4) = (-SIN(t1*w)+SIN((dt+t1)*w))**2/dt**2
mtn1(6,5) = (4.*(SIN((dt/2.+t1)*w)-SIN((dt+t1)*w))*(-SIN(t1*w)+  &
    & SIN((dt+t1)*w)))/dt**2
mtn1(6,6) = (w*COS((dt+t1)*w)*(-SIN(t1*w)+SIN((dt+t1)*w)))/dt
RETURN
END SUBROUTINE stcoeff_m
