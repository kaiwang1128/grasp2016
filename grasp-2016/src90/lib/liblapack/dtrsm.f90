      SUBROUTINE DTRSM(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:39   2/12/04  
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
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: LDB 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      CHARACTER  :: SIDE 
      CHARACTER  :: UPLO 
      CHARACTER  :: TRANSA 
      CHARACTER  :: DIAG 
      REAL(DOUBLE) , INTENT(IN) :: A(LDA,*) 
      REAL(DOUBLE) , INTENT(INOUT) :: B(LDB,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, J, K, NROWA 
      REAL(DOUBLE) :: TEMP 
      LOGICAL :: LSIDE, NOUNIT, UPPER 
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
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!     .. External Subroutines ..
!     .. Intrinsic Functions ..
!     .. Local Scalars ..
!     .. Parameters ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE = LSAME(SIDE,'L') 
      IF (LSIDE) THEN 
         NROWA = M 
      ELSE 
         NROWA = N 
      ENDIF 
      NOUNIT = LSAME(DIAG,'N') 
      UPPER = LSAME(UPLO,'U') 
!
      INFO = 0 
      IF (.NOT.LSIDE .AND. .NOT.LSAME(SIDE,'R')) THEN 
         INFO = 1 
      ELSE IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = 2 
      ELSE IF (.NOT.LSAME(TRANSA,'N') .AND. .NOT.LSAME(TRANSA,'T') .AND. .NOT.&
            LSAME(TRANSA,'C')) THEN 
         INFO = 3 
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN 
         INFO = 4 
      ELSE IF (M < 0) THEN 
         INFO = 5 
      ELSE IF (N < 0) THEN 
         INFO = 6 
      ELSE IF (LDA < MAX(1,NROWA)) THEN 
         INFO = 9 
      ELSE IF (LDB < MAX(1,M)) THEN 
         INFO = 11 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DTRSM ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N == 0) RETURN  
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA == ZERO) THEN 
         B(:M,:N) = ZERO 
         RETURN  
      ENDIF 
!
!     Start the operations.
!
      IF (LSIDE) THEN 
         IF (LSAME(TRANSA,'N')) THEN 
!
!           Form  B := alpha*inv( A )*B.
!
            IF (UPPER) THEN 
               DO J = 1, N 
                  IF (ALPHA /= ONE) THEN 
                     B(:M,J) = ALPHA*B(:M,J) 
                  ENDIF 
                  DO K = M, 1, -1 
                     IF (B(K,J) == ZERO) CYCLE  
                     IF (NOUNIT) B(K,J) = B(K,J)/A(K,K) 
                     B(:K-1,J) = B(:K-1,J) - B(K,J)*A(:K-1,K) 
                  END DO 
               END DO 
            ELSE 
               DO J = 1, N 
                  IF (ALPHA /= ONE) THEN 
                     B(:M,J) = ALPHA*B(:M,J) 
                  ENDIF 
                  DO K = 1, M 
                     IF (B(K,J) == ZERO) CYCLE  
                     IF (NOUNIT) B(K,J) = B(K,J)/A(K,K) 
                     B(K+1:M,J) = B(K+1:M,J) - B(K,J)*A(K+1:M,K) 
                  END DO 
               END DO 
            ENDIF 
         ELSE 
!
!           Form  B := alpha*inv( A' )*B.
!
            IF (UPPER) THEN 
               DO J = 1, N 
                  DO I = 1, M 
                     TEMP = ALPHA*B(I,J) 
                     TEMP = TEMP - SUM(A(:I-1,I)*B(:I-1,J)) 
                     IF (NOUNIT) TEMP = TEMP/A(I,I) 
                     B(I,J) = TEMP 
                  END DO 
               END DO 
            ELSE 
               DO J = 1, N 
                  DO I = M, 1, -1 
                     TEMP = ALPHA*B(I,J) 
                     TEMP = TEMP - SUM(A(I+1:M,I)*B(I+1:M,J)) 
                     IF (NOUNIT) TEMP = TEMP/A(I,I) 
                     B(I,J) = TEMP 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      ELSE 
         IF (LSAME(TRANSA,'N')) THEN 
!
!           Form  B := alpha*B*inv( A ).
!
            IF (UPPER) THEN 
               DO J = 1, N 
                  IF (ALPHA /= ONE) THEN 
                     B(:M,J) = ALPHA*B(:M,J) 
                  ENDIF 
                  DO K = 1, J - 1 
                     IF (A(K,J) == ZERO) CYCLE  
                     B(:M,J) = B(:M,J) - A(K,J)*B(:M,K) 
                  END DO 
                  IF (.NOT.NOUNIT) CYCLE  
                  TEMP = ONE/A(J,J) 
                  B(:M,J) = TEMP*B(:M,J) 
               END DO 
            ELSE 
               DO J = N, 1, -1 
                  IF (ALPHA /= ONE) THEN 
                     B(:M,J) = ALPHA*B(:M,J) 
                  ENDIF 
                  DO K = J + 1, N 
                     IF (A(K,J) == ZERO) CYCLE  
                     B(:M,J) = B(:M,J) - A(K,J)*B(:M,K) 
                  END DO 
                  IF (.NOT.NOUNIT) CYCLE  
                  TEMP = ONE/A(J,J) 
                  B(:M,J) = TEMP*B(:M,J) 
               END DO 
            ENDIF 
         ELSE 
!
!           Form  B := alpha*B*inv( A' ).
!
            IF (UPPER) THEN 
               DO K = N, 1, -1 
                  IF (NOUNIT) THEN 
                     TEMP = ONE/A(K,K) 
                     B(:M,K) = TEMP*B(:M,K) 
                  ENDIF 
                  DO J = 1, K - 1 
                     IF (A(J,K) == ZERO) CYCLE  
                     TEMP = A(J,K) 
                     B(:M,J) = B(:M,J) - TEMP*B(:M,K) 
                  END DO 
                  IF (ALPHA == ONE) CYCLE  
                  B(:M,K) = ALPHA*B(:M,K) 
               END DO 
            ELSE 
               DO K = 1, N 
                  IF (NOUNIT) THEN 
                     TEMP = ONE/A(K,K) 
                     B(:M,K) = TEMP*B(:M,K) 
                  ENDIF 
                  DO J = K + 1, N 
                     IF (A(J,K) == ZERO) CYCLE  
                     TEMP = A(J,K) 
                     B(:M,J) = B(:M,J) - TEMP*B(:M,K) 
                  END DO 
                  IF (ALPHA == ONE) CYCLE  
                  B(:M,K) = ALPHA*B(:M,K) 
               END DO 
            ENDIF 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DTRSM .
!
      END SUBROUTINE DTRSM 
