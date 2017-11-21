*
*     ------------------------------------------------------------------
*	C O P V E C
*     ------------------------------------------------------------------
*
      SUBROUTINE COPVEC(FROM,TO,NDIM)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION FROM(1),TO(1)
C
      DO 100 I=1,NDIM
       TO(I)=FROM(I)
  100 CONTINUE
C
      RETURN
      END
