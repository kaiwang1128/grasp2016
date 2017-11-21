      SUBROUTINE DLASR(SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:33   2/12/04  
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
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: LDA 
      CHARACTER  :: SIDE 
      CHARACTER  :: PIVOT 
      CHARACTER  :: DIRECT 
      REAL(DOUBLE) , INTENT(IN) :: C(*) 
      REAL(DOUBLE) , INTENT(IN) :: S(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: A(LDA,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, J 
      REAL(DOUBLE) :: CTEMP, STEMP, TEMP 
      LOGICAL::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0 
      IF (.NOT.(LSAME(SIDE,'L') .OR. LSAME(SIDE,'R'))) THEN 
         INFO = 1 
      ELSE IF (.NOT.(LSAME(PIVOT,'V') .OR. LSAME(PIVOT,'T') .OR. LSAME(PIVOT,&
            'B'))) THEN 
         INFO = 2 
      ELSE IF (.NOT.(LSAME(DIRECT,'F') .OR. LSAME(DIRECT,'B'))) THEN 
         INFO = 3 
      ELSE IF (M < 0) THEN 
         INFO = 4 
      ELSE IF (N < 0) THEN 
         INFO = 5 
      ELSE IF (LDA < MAX(1,M)) THEN 
         INFO = 9 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DLASR ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (M==0 .OR. N==0) RETURN  
      IF (LSAME(SIDE,'L')) THEN 
!
!        Form  P * A
!
         IF (LSAME(PIVOT,'V')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 1, M - 1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J+1,I) 
                     A(J+1,I) = CTEMP*TEMP - STEMP*A(J,I) 
                     A(J,I) = STEMP*TEMP + CTEMP*A(J,I) 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = M - 1, 1, -1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J+1,I) 
                     A(J+1,I) = CTEMP*TEMP - STEMP*A(J,I) 
                     A(J,I) = STEMP*TEMP + CTEMP*A(J,I) 
                  END DO 
               END DO 
            ENDIF 
         ELSE IF (LSAME(PIVOT,'T')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 2, M 
                  CTEMP = C(J-1) 
                  STEMP = S(J-1) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J,I) 
                     A(J,I) = CTEMP*TEMP - STEMP*A(1,I) 
                     A(1,I) = STEMP*TEMP + CTEMP*A(1,I) 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = M, 2, -1 
                  CTEMP = C(J-1) 
                  STEMP = S(J-1) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J,I) 
                     A(J,I) = CTEMP*TEMP - STEMP*A(1,I) 
                     A(1,I) = STEMP*TEMP + CTEMP*A(1,I) 
                  END DO 
               END DO 
            ENDIF 
         ELSE IF (LSAME(PIVOT,'B')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 1, M - 1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J,I) 
                     A(J,I) = STEMP*A(M,I) + CTEMP*TEMP 
                     A(M,I) = CTEMP*A(M,I) - STEMP*TEMP 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = M - 1, 1, -1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, N 
                     TEMP = A(J,I) 
                     A(J,I) = STEMP*A(M,I) + CTEMP*TEMP 
                     A(M,I) = CTEMP*A(M,I) - STEMP*TEMP 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      ELSE IF (LSAME(SIDE,'R')) THEN 
!
!        Form A * P'
!
         IF (LSAME(PIVOT,'V')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 1, N - 1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J+1) 
                     A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J) 
                     A(I,J) = STEMP*TEMP + CTEMP*A(I,J) 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = N - 1, 1, -1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J+1) 
                     A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J) 
                     A(I,J) = STEMP*TEMP + CTEMP*A(I,J) 
                  END DO 
               END DO 
            ENDIF 
         ELSE IF (LSAME(PIVOT,'T')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 2, N 
                  CTEMP = C(J-1) 
                  STEMP = S(J-1) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J) 
                     A(I,J) = CTEMP*TEMP - STEMP*A(I,1) 
                     A(I,1) = STEMP*TEMP + CTEMP*A(I,1) 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = N, 2, -1 
                  CTEMP = C(J-1) 
                  STEMP = S(J-1) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J) 
                     A(I,J) = CTEMP*TEMP - STEMP*A(I,1) 
                     A(I,1) = STEMP*TEMP + CTEMP*A(I,1) 
                  END DO 
               END DO 
            ENDIF 
         ELSE IF (LSAME(PIVOT,'B')) THEN 
            IF (LSAME(DIRECT,'F')) THEN 
               DO J = 1, N - 1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J) 
                     A(I,J) = STEMP*A(I,N) + CTEMP*TEMP 
                     A(I,N) = CTEMP*A(I,N) - STEMP*TEMP 
                  END DO 
               END DO 
            ELSE IF (LSAME(DIRECT,'B')) THEN 
               DO J = N - 1, 1, -1 
                  CTEMP = C(J) 
                  STEMP = S(J) 
                  IF (CTEMP==ONE .AND. STEMP==ZERO) CYCLE  
                  DO I = 1, M 
                     TEMP = A(I,J) 
                     A(I,J) = STEMP*A(I,N) + CTEMP*TEMP 
                     A(I,N) = CTEMP*A(I,N) - STEMP*TEMP 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DLASR
!
      END SUBROUTINE DLASR 
