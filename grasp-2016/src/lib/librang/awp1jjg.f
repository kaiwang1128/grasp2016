      SUBROUTINE AWP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,AW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      AW=ZERO
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK1=K1*2
      KK2=BK2+BK2+TENTH*BK2
      IQ2=QM2*TWO+QM2*TENTH
      IQ3=QM3*TWO+QM3*TENTH
      IQ=IQ2+IQ2
      IF(ID(3).GT.37) RETURN
      IF(ITJJ2(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ).EQ.0)RETURN
      IQM=TWO*DABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM.LE.IBT(7)) THEN
          IF(IXJTIK(IK(3),KK1,KK2,ID(6),IK(6),IBT(6)).NE.0) THEN
            IBT(1)=IT
            BT(2)=DBLE(IBT(6))/TWO
            BT(1)=DBLE(IBT(7))/TWO
            CALL A1JJ(IK,BK,IBT,BT,QM1,D1)
            IF(DABS(D1).GT.EPS) THEN
              CALL W1JJG(K1,QM2,QM3,IBT,BT,ID,BD,W)
              IF(DABS(W).GT.EPS) THEN
                D1=D1*W
                CALL SIXJ(IK(3),KK1,KK2,ID(6),IK(6),IBT(6),0,SI1)
                ENQP=ENQP+D1*SI1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
    1 CONTINUE
      AW=ENQP*DSQRT(DBLE((KK2+1)))
      IE=KK2+IK(6)+ID(6)
      IF(((IE/4)*4).NE.IE)AW=-AW
      RETURN
      END
