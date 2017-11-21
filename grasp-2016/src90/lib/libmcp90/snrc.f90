!***********************************************************************
!                                                                      *
      SUBROUTINE SNRC (IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
!                                                                      *
!***********************************************************************
!
      IMPLICIT NONE
!      DIMENSION IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(OUT) :: ND1
      INTEGER, INTENT(OUT) :: ND2
      INTEGER, INTENT(OUT) :: NE1
      INTEGER, INTENT(OUT) :: NE2
      INTEGER, INTENT(OUT) :: IBRD
      INTEGER, INTENT(OUT) :: IBRE
      INTEGER, INTENT(IN) :: IS(4)
      INTEGER, INTENT(IN) :: KAPS(4)
      INTEGER, INTENT(IN) :: KS(4)
!-----------------------------------------------
      STOP 'SNRC: Error '
      END
