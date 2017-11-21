!
!     ---------------------------------------------------------------
!      I T R E X G 
!     ---------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      * 
!     Vanderbilt University,  Nashville             October 1996   * 
!
      INTEGER FUNCTION ITREXG(I1,I2,I3,I4,K)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: I1, I2, I3, I4
      INTEGER, INTENT(OUT) :: K
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J
!-----------------------------------------------
      J=MAX0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG=J
      RETURN
      END FUNCTION ITREXG
