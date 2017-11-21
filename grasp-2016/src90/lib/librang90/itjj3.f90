!******************************************************************
!                                                                  *
      INTEGER FUNCTION ITJJ3(IK,ID,KG1,BK,BD,IBT,BT,ITP,ITG,IQ)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 08  --------------   *
!                                                                  *
!     FUNCTION CALLED: ITTK                                        *
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
      USE ribojj_C
      USE ribojj9_C
      USE ribojj11_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE mes_I
      USE ittk_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)                :: KG1, IQ
      INTEGER,      INTENT(OUT)               :: ITP, ITG
      INTEGER,      INTENT(IN),  DIMENSION(7) :: IK, ID
      INTEGER,      INTENT(OUT), DIMENSION(7) :: IBT
      REAL(DOUBLE), INTENT(IN),  DIMENSION(3) :: BK, BD
      REAL(DOUBLE), INTENT(OUT), DIMENSION(3) :: BT
!      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITK, ITD, ITP1, ITG1
!-----------------------------------------------
      ITJJ3=0
      IF(ID(3).GT.37) RETURN
      IF(ITTK(ID(6),IK(6),KG1).EQ.0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3).LT.9) THEN
        ITP1=IMPTJJ(ITK)
        ITP=IMPNJJ(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGTJJ(ITK)
        ITG=IMGNJJ(ITD)
      ELSEIF(ID(3).EQ.9) THEN
        IF(ITK.GT.300) THEN
          IF(ITD.LT.300) CALL MES(53)
          IF(ID(4).GT.2) CALL MES(13)
          IF(IK(4).GT.2) CALL MES(13)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPTJJ9(ITK)
          ITP=IMPNJJ9(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGTJJ9(ITK)
          ITG=IMGNJJ9(ITD)
        ELSE
          PRINT*, "ERROR in ITJJ3"
          STOP
        ENDIF
      ELSE
        IF(ID(4).GT.2) CALL MES(13)
        IF(IK(4).GT.2) CALL MES(13)
        ITP1=IMPTJJ11(ITK)
        ITP=IMPNJJ11(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGTJJ11(ITK)
        ITG=IMGNJJ11(ITD)
      ENDIF
      IF(ITG1.NE.ITG)RETURN
      ITJJ3=1
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      RETURN
      END FUNCTION ITJJ3
