!
!     -------------------------------------------------------------
!      I Z A S 1
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                             March 1995   *
!
      INTEGER FUNCTION IZAS1(IB,QB,IK,QK)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN) :: IB, IK
      REAL(DOUBLE), INTENT(IN) :: QB, QK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IQB, IQK
!-----------------------------------------------
      IZAS1=0
      IQB=TWO*DABS(QB)+TENTH
      IF(IQB.GT.IB)RETURN
      IF(MOD(IB+IQB,2).NE.0)RETURN
      IQK=TWO*DABS(QK)+TENTH
      IF(IQK.GT.IK)RETURN
      IF(MOD(IK+IQK,2).NE.0)RETURN
      IZAS1=1
      RETURN
      END  FUNCTION IZAS1
