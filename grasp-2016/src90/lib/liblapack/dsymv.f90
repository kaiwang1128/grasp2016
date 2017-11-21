      SUBROUTINE DSYMV(UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:37   2/12/04  
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
      INTEGER :: I, INFO, IX, IY, J, JX, JY, KX, KY 
      REAL(DOUBLE) :: TEMP1, TEMP2 
      logical::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSYMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
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
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
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
      ELSE IF (LDA < MAX(1,N)) THEN 
         INFO = 5 
      ELSE IF (INCX == 0) THEN 
         INFO = 7 
      ELSE IF (INCY == 0) THEN 
         INFO = 10 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSYMV ', INFO) 
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
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
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
!        Form  y  when A is stored in upper triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               Y(:J-1) = Y(:J-1) + TEMP1*A(:J-1,J) 
               TEMP2 = TEMP2 + SUM(A(:J-1,J)*X(:J-1)) 
               Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               IX = KX 
               IY = KY 
               Y(IY:(J-2)*INCY+IY:INCY) = Y(IY:(J-2)*INCY+IY:INCY) + TEMP1*A(:J&
                  -1,J) 
               TEMP2 = TEMP2 + SUM(A(:J-1,J)*X(IX:(J-2)*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + TEMP1*A(J,J) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  y  when A is stored in lower triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               TEMP1 = ALPHA*X(J) 
               TEMP2 = ZERO 
               Y(J) = Y(J) + TEMP1*A(J,J) 
               Y(J+1:N) = Y(J+1:N) + TEMP1*A(J+1:N,J) 
               TEMP2 = TEMP2 + SUM(A(J+1:N,J)*X(J+1:N)) 
               Y(J) = Y(J) + ALPHA*TEMP2 
            END DO 
         ELSE 
            JX = KX 
            JY = KY 
            DO J = 1, N 
               TEMP1 = ALPHA*X(JX) 
               TEMP2 = ZERO 
               Y(JY) = Y(JY) + TEMP1*A(J,J) 
               IX = JX 
               IY = JY 
               Y(IY+INCY:(N-J)*INCY+IY:INCY) = Y(IY+INCY:(N-J)*INCY+IY:INCY) + &
                  TEMP1*A(J+1:N,J) 
               TEMP2 = TEMP2 + SUM(A(J+1:N,J)*X(IX+INCX:(N-J)*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + ALPHA*TEMP2 
               JX = JX + INCX 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSYMV .
!
      END SUBROUTINE DSYMV 
