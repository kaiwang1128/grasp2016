!*******************************************************************
!                                                                  *
      SUBROUTINE EILE(JA,JB,JC,JAA,JBB,JCC)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 02  ------------   *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: JA, JB, JC 
      INTEGER, INTENT(OUT) :: JAA, JBB, JCC
!-----------------------------------------------
      JAA=JA
      JCC=JA
      IF(JAA.GT.JB)JAA=JB
      IF(JCC.LT.JB)JCC=JB
      IF(JAA.GT.JC)JAA=JC
      IF(JCC.LT.JC)JCC=JC
      IF((JA.GT.JAA).AND.(JA.LT.JCC))JBB=JA
      IF((JB.GT.JAA).AND.(JB.LT.JCC))JBB=JB
      IF((JC.GT.JAA).AND.(JC.LT.JCC))JBB=JC
      RETURN
      END SUBROUTINE EILE
