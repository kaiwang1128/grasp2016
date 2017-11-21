      SUBROUTINE DSPMV(UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:35   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE lsame_I 
      !USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: INCY 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: AP(*) 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY 
      REAL(DOUBLE) :: TEMP1, TEMP2 
      logical::lsame
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSPMV  performs the matrix-vector operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix, supplied in packed form.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
!     .. Local Scalars ..
!     .. External Functions ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0 
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = 1 
      ELSE IF (N < 0) THEN 
         INFO = 2 
      ELSE IF (INCX == 0) THEN 
         INFO = 6 
      ELSE IF (INCY == 0) THEN 
         INFO = 9 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSPMV ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N==0 .OR. ALPHA==ZERO .AND. BETA==ONE) RETURN  
!
!     Set up the start points in  X  and  Y.
!
      IF (INCX > 0) THEN 
         KX = 1 
      ELSE 
         KX = 1 - (N - 1)*INCX 
      ENDIF 
      IF (INCY > 0) THEN 
         KY = 1 
      ELSE 
         KY = 1 - (N - 1)*INCY 
      ENDIF 
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
!     First form  y := beta*y.
!
      IF (BETA /= ONE) THEN 
         IF (INCY == 1) THEN 
            IF (BETA == ZERO) THEN 
               Y(:N) = ZERO 
            ELSE 
               Y(:N) = BETA*Y(:N) 
            ENDIF 
         ELSE 
            IY = KY 
            IF (BETA == ZERO) THEN 
               Y(IY:(N-1)*INCY+IY:INCY) = ZERO 
            ELSE 
               Y(IY:(N-1)*INCY+IY:INCY) = BETA*Y(IY:(N-1)*INCY+IY:INCY) 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ALPHA == ZERO) RETURN  
      KK = 1 
      IF (LSAME(UPLO,'U')) THEN 
!
!        Form  y  when AP contains the upper triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               K = KK 
               Y(:J-1) = Y(:J-1) + TEMP1*AP(K:J-2+K) 
               TEMP2 = TEMP2 + SUM(AP(K:J-2+K)*X(:J-1)) 
               Y(J) = Y(J) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2 
               KK = KK + J 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               IX = KX 
               IY = KY 
               Y(IY:(J-2)*INCY+IY:INCY) = Y(IY:(J-2)*INCY+IY:INCY) + TEMP1*AP(&
                  KK:J-2+KK) 
               TEMP2 = TEMP2 + SUM(AP(KK:J-2+KK)*X(IX:(J-2)*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + TEMP1*AP(KK+J-1) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
               KK = KK + J 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  y  when AP contains the lower triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               Y(J) = Y(J) + TEMP1*AP(KK) 
               K = KK + 1 
               Y(J+1:N) = Y(J+1:N) + TEMP1*AP(K:N-J-1+K) 
               TEMP2 = TEMP2 + SUM(AP(K:N-J-1+K)*X(J+1:N)) 
               Y(J) = Y(J) + ALPHA*TEMP2 
               KK = KK + (N - J + 1) 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               Y(JY) = Y(JY) + TEMP1*AP(KK) 
               IX = JX 
               IY = JY 
               Y(IY+INCY:(N-J)*INCY+IY:INCY) = Y(IY+INCY:(N-J)*INCY+IY:INCY) + &
                  TEMP1*AP(KK+1:N-J+KK) 
               TEMP2 = TEMP2 + SUM(AP(KK+1:N-J+KK)*X(IX+INCX:(N-J)*INCX+IX:INCX&
                  )) 
               Y(JY) = Y(JY) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
               KK = KK + (N - J + 1) 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSPMV .
!
      END SUBROUTINE DSPMV 
