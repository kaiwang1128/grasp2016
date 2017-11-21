      SUBROUTINE SCATTER(N, A, INDEX, B) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: INDEX(1) 
      REAL , INTENT(OUT) :: A(1) 
      REAL , INTENT(IN) :: B(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     This subroutine disperses array elements accessed via an
!         integer pointer from contiguous storage to the appropriate
!         location.
!
 
      DO I = 1, N 
         A(INDEX(I)) = B(I) 
      END DO 
 
      RETURN  
      END SUBROUTINE SCATTER 
