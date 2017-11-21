      SUBROUTINE DTPMV(UPLO, TRANS, DIAG, N, AP, X, INCX) 
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
      INTEGER , INTENT(IN) :: INCX 
      CHARACTER  :: UPLO 
      CHARACTER  :: TRANS 
      CHARACTER  :: DIAG 
      REAL(DOUBLE) , INTENT(IN) :: AP(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, J, JX, K, KK, KX 
      REAL(DOUBLE) :: TEMP 
      LOGICAL :: NOUNIT 
      logical::lsame
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DTPMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix, supplied in packed form.
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
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
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
      ELSE IF (INCX == 0) THEN 
         INFO = 7 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DTPMV ', INFO) 
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
!     Start the operations. In this version the elements of AP are
!     accessed sequentially with one pass through AP.
!
      IF (LSAME(TRANS,'N')) THEN 
!
!        Form  x:= A*x.
!
         IF (LSAME(UPLO,'U')) THEN 
            KK = 1 
            IF (INCX == 1) THEN 
               DO J = 1, N 
                  IF (X(J) /= ZERO) THEN 
                     TEMP = X(J) 
                     K = KK 
                     X(:J-1) = X(:J-1) + TEMP*AP(K:J-2+K) 
                     IF (NOUNIT) X(J) = X(J)*AP(KK+J-1) 
                  ENDIF 
                  KK = KK + J 
               END DO 
            ELSE 
               JX = KX 
               DO J = 1, N 
                  IF (X(JX) /= ZERO) THEN 
                     TEMP = X(JX) 
                     IX = KX 
                     X(IX:(J-2)*INCX+IX:INCX) = X(IX:(J-2)*INCX+IX:INCX) + TEMP&
                        *AP(KK:J-2+KK) 
                     IF (NOUNIT) X(JX) = X(JX)*AP(KK+J-1) 
                  ENDIF 
                  JX = JX + INCX 
                  KK = KK + J 
               END DO 
            ENDIF 
         ELSE 
            KK = (N*(N + 1))/2 
            IF (INCX == 1) THEN 
               DO J = N, 1, -1 
                  IF (X(J) /= ZERO) THEN 
                     TEMP = X(J) 
                     K = KK 
                     X(N:1+J:(-1)) = X(N:1+J:(-1)) + TEMP*AP(K:K+1+J-N:(-1)) 
                     IF (NOUNIT) X(J) = X(J)*AP(KK-N+J) 
                  ENDIF 
                  KK = KK - (N - J + 1) 
               END DO 
            ELSE 
               KX = KX + (N - 1)*INCX 
               JX = KX 
               DO J = N, 1, -1 
                  IF (X(JX) /= ZERO) THEN 
                     TEMP = X(JX) 
                     IX = KX 
                     X(IX:INCX*(J-N+1)+IX:(-INCX)) = X(IX:INCX*(J-N+1)+IX:(-&
                        INCX)) + TEMP*AP(KK:KK+1+J-N:(-1)) 
                     IF (NOUNIT) X(JX) = X(JX)*AP(KK-N+J) 
                  ENDIF 
                  JX = JX - INCX 
                  KK = KK - (N - J + 1) 
               END DO 
            ENDIF 
         ENDIF 
      ELSE 
!
!        Form  x := A'*x.
!
         IF (LSAME(UPLO,'U')) THEN 
            KK = (N*(N + 1))/2 
            IF (INCX == 1) THEN 
               DO J = N, 1, -1 
                  TEMP = X(J) 
                  IF (NOUNIT) TEMP = TEMP*AP(KK) 
                  K = KK - 1 
                  TEMP = TEMP + SUM(AP(K:K+2-J:(-1))*X(J-1:1:(-1))) 
                  X(J) = TEMP 
                  KK = KK - J 
               END DO 
            ELSE 
               JX = KX + (N - 1)*INCX 
               DO J = N, 1, -1 
                  TEMP = X(JX) 
                  IX = JX 
                  IF (NOUNIT) TEMP = TEMP*AP(KK) 
                  TEMP = TEMP + SUM(AP(KK-1:KK+1-J:(-1))*X(IX-INCX:INCX*(1-J)+&
                     IX:(-INCX))) 
                  X(JX) = TEMP 
                  JX = JX - INCX 
                  KK = KK - J 
               END DO 
            ENDIF 
         ELSE 
            KK = 1 
            IF (INCX == 1) THEN 
               DO J = 1, N 
                  TEMP = X(J) 
                  IF (NOUNIT) TEMP = TEMP*AP(KK) 
                  K = KK + 1 
                  TEMP = TEMP + SUM(AP(K:N-J-1+K)*X(J+1:N)) 
                  X(J) = TEMP 
                  KK = KK + (N - J + 1) 
               END DO 
            ELSE 
               JX = KX 
               DO J = 1, N 
                  TEMP = X(JX) 
                  IX = JX 
                  IF (NOUNIT) TEMP = TEMP*AP(KK) 
                  TEMP = TEMP + SUM(AP(KK+1:N-J+KK)*X(IX+INCX:(N-J)*INCX+IX:&
                     INCX)) 
                  X(JX) = TEMP 
                  JX = JX + INCX 
                  KK = KK + (N - J + 1) 
               END DO 
            ENDIF 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DTPMV .
!
      END SUBROUTINE DTPMV 
