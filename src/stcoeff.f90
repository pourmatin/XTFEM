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
REAL(KIND=REKIND), DIMENSION(3, 3), INTENT(OUT) :: k, m, ktn, mtn, ktn1, mtn1
! ---Internal variables---
INTEGER :: i, j

! Coefficients for Stiffness matrix of the complete slab
k(1,1) = 2._REKIND/15._REKIND*dt
k(1,2) = 1._REKIND/15._REKIND*dt
k(1,3) = -1._REKIND/30._REKIND*dt
	
k(2,1) = 1._REKIND/15._REKIND*dt
k(2,2) = 8._REKIND/15._REKIND*dt
k(2,3) = 1._REKIND/15._REKIND*dt
     
k(3,1) = -1._REKIND/30._REKIND*dt
k(3,2) = 1._REKIND/15._REKIND*dt
k(3,3) = 2._REKIND/15._REKIND*dt

! Coefficients for Mass matrix of the complete slab
m(1,1) = -.5_REKIND
m(1,2) = 2._REKIND/3._REKIND
m(1,3) = -1._REKIND/6._REKIND
  
m(2,1) = -2._REKIND/3._REKIND
m(2,2) = 0._REKIND
m(2,3) = 2._REKIND/3._REKIND

m(3,1) = 1._REKIND/6._REKIND
m(3,2) = -2._REKIND/3._REKIND
m(3,3) = .5_REKIND

! Coefficients for the Stiffness matrix at tn_1 (+)
ktn = 0._REKIND

! Coefficients for the Mass matrix at tn_1 (+)
mtn = 0._REKIND
mtn(1,1) = 1._REKIND

! Coefficients for the Stiffness matrix at tn_1 (-) 
ktn1 = 0._REKIND

! Coefficients for the Mass matrix at tn_1 (-)
mtn1 = 0._REKIND
mtn1(1,3) = 1._REKIND
	
RETURN
END SUBROUTINE stcoeff
