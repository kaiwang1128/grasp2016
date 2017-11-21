      SUBROUTINE IINIT(N, A, X, INCX) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:40   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: A 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(OUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: XADDR, I 
!-----------------------------------------------
!
!     ==================================================================
!
!     PURPOSE ... INITIALIZES INTEGER VECTOR TO
!                 A CONSTANT VALUE 'A'
!
!     CREATED ... AUG. 06, 1992
!
!     ==================================================================
!
!     -------------
!     ... ARGUMENTS
!     -------------
!
!
!     -------------------
!     ... LOCAL VARIABLES
!     -------------------
!
!
!     ==================================================================
!
      IF (INCX == 1) THEN 
!
!         ----------------------------------
!         ... UNIT INCREMENT (STANDARD CASE)
!         ----------------------------------
!
         X(:N) = A 
!
      ELSE 
!
!         ----------------------
!         ... NON-UNIT INCREMENT
!         ----------------------
!
         XADDR = 1 
         IF (INCX < 0) XADDR = ((-N) + 1)*INCX + 1 
!
         X(XADDR:(N-1)*INCX+XADDR:INCX) = A 
!
      ENDIF 
!
      RETURN  
!
      END SUBROUTINE IINIT 
