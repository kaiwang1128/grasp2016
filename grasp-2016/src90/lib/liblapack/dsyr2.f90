      SUBROUTINE DSYR2(UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:38   2/12/04  
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
      INTEGER , INTENT(IN) :: LDA 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(IN) :: Y(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: A(LDA,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
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
!  DSYR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.
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
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
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
      ELSE IF (INCX == 0) THEN 
         INFO = 5 
      ELSE IF (INCY == 0) THEN 
         INFO = 7 
      ELSE IF (LDA < MAX(1,N)) THEN 
         INFO = 9 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSYR2 ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N==0 .OR. ALPHA==ZERO) RETURN  
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
      IF (INCX/=1 .OR. INCY/=1) THEN 
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
         JX = KX 
         JY = KY 
      ENDIF 
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF (LSAME(UPLO,'U')) THEN 
!
!        Form  A  when A is stored in the upper triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               IF (X(J)==ZERO .AND. Y(J)==ZERO) CYCLE  
               TEMP1 = ALPHA*Y(J) 
               TEMP2 = ALPHA*X(J) 
               A(:J,J) = A(:J,J) + X(:J)*TEMP1 + Y(:J)*TEMP2 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (X(JX)/=ZERO .OR. Y(JY)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(JY) 
                  TEMP2 = ALPHA*X(JX) 
                  IX = KX 
                  IY = KY 
                  A(:J,J) = A(:J,J) + X(IX:(J-1)*INCX+IX:INCX)*TEMP1 + Y(IY:(J-&
                     1)*INCY+IY:INCY)*TEMP2 
               ENDIF 
               JX = JX + INCX 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  A  when A is stored in the lower triangle.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               IF (X(J)==ZERO .AND. Y(J)==ZERO) CYCLE  
               TEMP1 = ALPHA*Y(J) 
               TEMP2 = ALPHA*X(J) 
               A(J:N,J) = A(J:N,J) + X(J:N)*TEMP1 + Y(J:N)*TEMP2 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (X(JX)/=ZERO .OR. Y(JY)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(JY) 
                  TEMP2 = ALPHA*X(JX) 
                  IX = JX 
                  IY = JY 
                  A(J:N,J) = A(J:N,J) + X(IX:(N-J)*INCX+IX:INCX)*TEMP1 + Y(IY:(&
                     N-J)*INCY+IY:INCY)*TEMP2 
               ENDIF 
               JX = JX + INCX 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSYR2 .
!
      END SUBROUTINE DSYR2 
