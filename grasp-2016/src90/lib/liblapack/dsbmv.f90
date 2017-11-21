      SUBROUTINE DSBMV(UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:35   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: INCY 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: A(LDA,*) 
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
      INTEGER :: I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L 
      REAL(DOUBLE) :: TEMP1, TEMP2 
      logical::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX, MIN 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSBMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric band matrix, with k super-diagonals.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the band matrix A is being supplied as
!           follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  being supplied.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  being supplied.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry, K specifies the number of super-diagonals of the
!           matrix A. K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the symmetric matrix, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer the upper
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the symmetric matrix, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer the lower
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
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
!     .. Intrinsic Functions ..
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
      ELSE IF (K < 0) THEN 
         INFO = 3 
      ELSE IF (LDA < K + 1) THEN 
         INFO = 6 
      ELSE IF (INCX == 0) THEN 
         INFO = 8 
      ELSE IF (INCY == 0) THEN 
         INFO = 11 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSBMV ', INFO) 
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
!     Start the operations. In this version the elements of the array A
!     are accessed sequentially with one pass through A.
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
      IF (LSAME(UPLO,'U')) THEN 
!
!        Form  y  when upper triangle of A is stored.
!
         KPLUS1 = K + 1 
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               L = KPLUS1 - J 
               Y(MAX(1,J-K):J-1) = Y(MAX(1,J-K):J-1) + TEMP1*A(L+MAX(1,J-K):J-1&
                  +L,J) 
               TEMP2 = TEMP2 + SUM(A(L+MAX(1,J-K):J-1+L,J)*X(MAX(1,J-K):J-1)) 
               Y(J) = Y(J) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               IX = KX 
               IY = KY 
               L = KPLUS1 - J 
               Y(IY:(J-MAX(1,J-K)-1)*INCY+IY:INCY) = Y(IY:(J-MAX(1,J-K)-1)*INCY&
                  +IY:INCY) + TEMP1*A(L+MAX(1,J-K):J-1+L,J) 
               TEMP2 = TEMP2 + SUM(A(L+MAX(1,J-K):J-1+L,J)*X(IX:(J-MAX(1,J-K)-1&
                  )*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
               IF (J <= K) CYCLE  
               KX = KX + INCX 
               KY = KY + INCY 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  y  when lower triangle of A is stored.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               Y(J) = Y(J) + TEMP1*A(1,J) 
               L = 1 - J 
               Y(J+1:MIN(N,J+K)) = Y(J+1:MIN(N,J+K)) + TEMP1*A(L+J+1:MIN(N,J+K)&
                  +L,J) 
               TEMP2 = TEMP2 + SUM(A(L+J+1:MIN(N,J+K)+L,J)*X(J+1:MIN(N,J+K))) 
               Y(J) = Y(J) + ALPHA*TEMP2 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               Y(JY) = Y(JY) + TEMP1*A(1,J) 
               L = 1 - J 
               IX = JX 
               IY = JY 
               Y(IY+INCY:(MIN(N,J+K)-J)*INCY+IY:INCY) = Y(IY+INCY:(MIN(N,J+K)-J&
                  )*INCY+IY:INCY) + TEMP1*A(L+J+1:MIN(N,J+K)+L,J) 
               TEMP2 = TEMP2 + SUM(A(L+J+1:MIN(N,J+K)+L,J)*X(IX+INCX:(MIN(N,J+K&
                  )-J)*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSBMV .
!
      END SUBROUTINE DSBMV 
