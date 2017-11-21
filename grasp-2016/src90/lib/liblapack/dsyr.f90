      SUBROUTINE DSYR(UPLO, N, ALPHA, X, INCX, A, LDA) 
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
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: LDA 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: A(LDA,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, J, JX, KX 
      REAL(DOUBLE) :: TEMP 
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
!  DSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
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
      ELSE IF (LDA < MAX(1,N)) THEN 
         INFO = 7 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSYR  ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N==0 .OR. ALPHA==ZERO) RETURN  
!
!     Set the start point in X if the increment is not unity.
!
      IF (INCX <= 0) THEN 
         KX = 1 - (N - 1)*INCX 
      ELSE IF (INCX /= 1) THEN 
         KX = 1 
      ENDIF 
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF (LSAME(UPLO,'U')) THEN 
!
!        Form  A  when A is stored in upper triangle.
!
         IF (INCX == 1) THEN 
            DO J = 1, N 
               IF (X(J) == ZERO) CYCLE  
               TEMP = ALPHA*X(J) 
               A(:J,J) = A(:J,J) + X(:J)*TEMP 
            END DO 
         ELSE 
            JX = KX 
            DO J = 1, N 
               IF (X(JX) /= ZERO) THEN 
                  TEMP = ALPHA*X(JX) 
                  IX = KX 
                  A(:J,J) = A(:J,J) + X(IX:(J-1)*INCX+IX:INCX)*TEMP 
               ENDIF 
               JX = JX + INCX 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  A  when A is stored in lower triangle.
!
         IF (INCX == 1) THEN 
            DO J = 1, N 
               IF (X(J) == ZERO) CYCLE  
               TEMP = ALPHA*X(J) 
               A(J:N,J) = A(J:N,J) + X(J:N)*TEMP 
            END DO 
         ELSE 
            JX = KX 
            DO J = 1, N 
               IF (X(JX) /= ZERO) THEN 
                  TEMP = ALPHA*X(JX) 
                  IX = JX 
                  A(J:N,J) = A(J:N,J) + X(IX:(N-J)*INCX+IX:INCX)*TEMP 
               ENDIF 
               JX = JX + INCX 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSYR  .
!
      END SUBROUTINE DSYR 
