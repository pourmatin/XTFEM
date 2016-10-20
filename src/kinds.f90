!  kinds.f03
!
!  Author: David N Alpert (david.alpert@gmail.com)
!  Advisor: Dr. Dong Qian (dong.qian@uc.edu)
!  University of Cincinnati, College of Engineering and Applied Science
!  School of Dynamic Systems, Mechanical Engineering Program
!  Computer Aided Engineering Research Laboratory
!
!  Record of Revisions
!
!      Date      Programmer   Description
!  ------------ ------------ -----------------------------------
!   07/06/2009     D.N.A.     Original program
!   11/04/2010     D.N.A.     Removed the use of LONG integers
!
!  Specifications for KIND parameters and PI parameter (to the 20th decimal place)

MODULE KINDS

IMPLICIT NONE

!---------------------------------------------------------------

!  Data dictionary

!---------------------------------------------------------------

!  REAL and COMPLEX KINDs
!
INTEGER, PARAMETER :: SGL = SELECTED_REAL_KIND(p=6, r=37)
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13, r=200)
INTEGER, PARAMETER :: REKIND = SELECTED_REAL_KIND(p=13, r=200)
!  INT KINDs
!
INTEGER, PARAMETER :: SHORT = SELECTED_INT_KIND(6)
INTEGER, PARAMETER :: LONG = SELECTED_INT_KIND(10)
INTEGER, PARAMETER :: INKIND = SELECTED_INT_KIND(6)

!  GEOMETRIC CONSTANTS
! PI to the 70th decimal place
REAL(kind=REKIND), PARAMETER :: PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164_REKIND

END MODULE KINDS
