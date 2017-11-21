!*******************************************************************
!                                                                  *
      SUBROUTINE RECOP1(NS,JA1,KA,IRE,IAT,RECC)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE m_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE diaga1_I
      USE diaga3_I
      USE diaga5_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: NS, JA1, KA, IRE
      INTEGER, INTENT(OUT)      :: IAT
      REAL(DOUBLE), INTENT(OUT) :: RECC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: IJ1, ISKR, NPEELGG
      REAL(DOUBLE) :: S, RE
!-----------------------------------------------
      IAT=1
      IJ1=JLIST(JA1)
      S=DBLE(JJQ1(3,IJ1))
      RECC=ONE/DSQRT(S)
      IF(NPEEL.EQ.1 .AND. NS.EQ.-1)RETURN
      IF(NS .EQ. -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IAT=0
      IF(IRE.NE.0) THEN
        IF(KA.EQ.0) THEN
          IAT=1
          RETURN
        END IF
      END IF
      IAT=1
      IF(NPEELGG.EQ.1) RETURN
      IAT=0
      IF(NPEELGG.NE.2) THEN
        CALL DIAGA5(NPEELGG,JA1,2*KA,IRE,IAT,RE)
        RECC=RE*RECC
        IF(IAT.EQ.0) RETURN
        IF(JA1.EQ.NPEELGG) RETURN
        IAT=0
      END IF
      CALL DIAGA1(JA1,2*KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      RECC=RE*RECC
      IF(NPEELGG.EQ.2) RETURN
      ISKR=NPEELGG-JA1
      IF(JA1.EQ.1)ISKR=NPEELGG-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DIAGA3(JA1,NPEELGG,2*KA,IRE,IAT,RE)
      RECC=RE*RECC
      RETURN
      END SUBROUTINE RECOP1
