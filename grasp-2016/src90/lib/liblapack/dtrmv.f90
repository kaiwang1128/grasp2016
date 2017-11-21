      SUBROUTINE DTRMV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX) 
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
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: INCX 
      CHARACTER  :: UPLO 
      CHARACTER  :: TRANS 
      CHARACTER  :: DIAG 
      REAL(DOUBLE) , INTENT(IN) :: A(LDA,*) 
      REAL(DOUBLE) , INTENT(INOUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, J, JX, KX 
      REAL(DOUBLE) :: TEMP 
      LOGICAL :: NOUNIT 
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
!  DTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
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
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
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
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.&
            LSAME(TRANS,'C')) THEN 
         INFO = 2 
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN 
         INFO = 3 
      ELSE IF (N < 0) THEN 
         INFO = 4 
      ELSE IF (LDA < MAX(1,N)) THEN 
         INFO = 6 
      ELSE IF (INCX == 0) THEN 
         INFO = 8 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DTRMV ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N == 0) RETURN  
!
      NOUNIT = LSAME(DIAG,'N') 
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF (INCX <= 0) THEN 
         KX = 1 - (N - 1)*INCX 
      ELSE IF (INCX /= 1) THEN 
         KX = 1 
      ENDIF 
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (LSAME(TRANS,'N')) THEN 
!
!        Form  x := A*x.
!
         IF (LSAME(UPLO,'U')) THEN 
            IF (INCX == 1) THEN 
               DO J = 1, N 
                  IF (X(J) == ZERO) CYCLE  
                  TEMP = X(J) 
                  X(:J-1) = X(:J-1) + TEMP*A(:J-1,J) 
                  IF (.NOT.NOUNIT) CYCLE  
                  X(J) = X(J)*A(J,J) 
               END DO 
            ELSE 
               JX = KX 
               DO J = 1, N 
                  IF (X(JX) /= ZERO) THEN 
                     TEMP = X(JX) 
                     IX = KX 
                     X(IX:(J-2)*INCX+IX:INCX) = X(IX:(J-2)*INCX+IX:INCX) + TEMP&
                        *A(:J-1,J) 
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J) 
                  ENDIF 
                  JX = JX + INCX 
               END DO 
            ENDIF 
         ELSE 
            IF (INCX == 1) THEN 
               DO J = N, 1, -1 
                  IF (X(J) == ZERO) CYCLE  
                  TEMP = X(J) 
                  X(N:1+J:(-1)) = X(N:1+J:(-1)) + TEMP*A(N:1+J:(-1),J) 
                  IF (.NOT.NOUNIT) CYCLE  
                  X(J) = X(J)*A(J,J) 
               END DO 
            ELSE 
               KX = KX + (N - 1)*INCX 
               JX = KX 
               DO J = N, 1, -1 
                  IF (X(JX) /= ZERO) THEN 
                     TEMP = X(JX) 
                     IX = KX 
                     X(IX:INCX*(J-N+1)+IX:(-INCX)) = X(IX:INCX*(J-N+1)+IX:(-&
                        INCX)) + TEMP*A(N:1+J:(-1),J) 
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J) 
                  ENDIF 
                  JX = JX - INCX 
               END DO 
            ENDIF 
         ENDIF 
      ELSE 
!
!        Form  x := A'*x.
!
         IF (LSAME(UPLO,'U')) THEN 
            IF (INCX == 1) THEN 
               DO J = N, 1, -1 
                  TEMP = X(J) 
                  IF (NOUNIT) TEMP = TEMP*A(J,J) 
                  TEMP = TEMP + SUM(A(J-1:1:(-1),J)*X(J-1:1:(-1))) 
                  X(J) = TEMP 
               END DO 
            ELSE 
               JX = KX + (N - 1)*INCX 
               DO J = N, 1, -1 
                  TEMP = X(JX) 
                  IX = JX 
                  IF (NOUNIT) TEMP = TEMP*A(J,J) 
                  TEMP = TEMP + SUM(A(J-1:1:(-1),J)*X(IX-INCX:INCX*(1-J)+IX:(-&
                     INCX))) 
                  X(JX) = TEMP 
                  JX = JX - INCX 
               END DO 
            ENDIF 
         ELSE 
            IF (INCX == 1) THEN 
               DO J = 1, N 
                  TEMP = X(J) 
                  IF (NOUNIT) TEMP = TEMP*A(J,J) 
                  TEMP = TEMP + SUM(A(J+1:N,J)*X(J+1:N)) 
                  X(J) = TEMP 
               END DO 
            ELSE 
               JX = KX 
               DO J = 1, N 
                  TEMP = X(JX) 
                  IX = JX 
                  IF (NOUNIT) TEMP = TEMP*A(J,J) 
                  TEMP = TEMP + SUM(A(J+1:N,J)*X(IX+INCX:(N-J)*INCX+IX:INCX)) 
                  X(JX) = TEMP 
                  JX = JX + INCX 
               END DO 
            ENDIF 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DTRMV .
!
      END SUBROUTINE DTRMV 
