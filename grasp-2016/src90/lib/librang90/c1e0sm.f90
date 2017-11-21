!
!     -------------------------------------------------------------
!      C 1 E 0 S M
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!
      SUBROUTINE C1E0SM(Q,QM,C,CM,A)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), INTENT(IN)  :: Q, QM, C, CM
      REAL(DOUBLE), INTENT(OUT) :: A
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IS, IG
!-----------------------------------------------
      A=ZERO
      IIQ=TWO*Q+TENTH
      IIC=TWO*C+TENTH
      IF(ITTK(IIQ,IIC,2).EQ.0)RETURN
      IF(DABS(QM-CM).GT.EPS) RETURN 
      IF((Q+TENTH).LT.DABS(QM)) RETURN
      IF((C+TENTH).LT.DABS(CM)) RETURN
      IF(DABS(QM).LE.EPS) THEN
       IS=Q+C+ONE+TENTH
       IF((IS/2)*2.NE.IS) RETURN
      END IF
      IG=Q-C+TWO+TENTH
      IF(IG.LE.0) RETURN
      IF(IG.GT.3) RETURN
      IF (IG .EQ. 1) THEN
        A=DSQRT(((C+CM)*(C-CM))/((TWO*C-ONE)*C))
      ELSE IF (IG .EQ. 2) THEN
        A=CM/DSQRT(C*(C+ONE))
      ELSE IF (IG .EQ. 3) THEN
        A=-DSQRT(((C+CM+ONE)*(C-CM+ONE))/((C+ONE)*(TWO*C+THREE)))
      END IF
      RETURN
      END SUBROUTINE C1E0SM
