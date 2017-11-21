      SUBROUTINE GATHER(N, A, B, INDEX) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:40   2/12/04  
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
!     This subroutine collects array elements accessed via an
!         integer pointer to contiguous storage.
!
 
      A(:N) = B(INDEX(:N)) 
 
      RETURN  
      END SUBROUTINE GATHER 
