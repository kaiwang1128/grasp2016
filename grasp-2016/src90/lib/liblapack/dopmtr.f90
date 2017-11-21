      SUBROUTINE DOPMTR(SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:34   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USE dlarf_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER  :: LDC 
      INTEGER , INTENT(OUT) :: INFO 
      CHARACTER  :: SIDE 
      CHARACTER  :: UPLO 
      CHARACTER  :: TRANS 
      REAL(DOUBLE)  :: AP(*) 
      REAL(DOUBLE)  :: TAU(*) 
      REAL(DOUBLE)  :: C(LDC,*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, I1, I2, I3, IC, II, JC, MI, NI, NQ 
      REAL(DOUBLE) :: AII 
      LOGICAL :: FORWRD, LEFT, NOTRAN, UPPER 
      logical::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DOPMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSPTRD using packed
!  storage:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension
!                               (M*(M+1)/2) if SIDE = 'L'
!                               (N*(N+1)/2) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.  AP is modified by the routine but
!          restored on exit.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L'
!                                     or (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L'
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!     Test the input arguments
!
      INFO = 0 
      LEFT = LSAME(SIDE,'L') 
      NOTRAN = LSAME(TRANS,'N') 
      UPPER = LSAME(UPLO,'U') 
!
!     NQ is the order of Q
!
      IF (LEFT) THEN 
         NQ = M 
      ELSE 
         NQ = N 
      ENDIF 
      IF (.NOT.LEFT .AND. .NOT.LSAME(SIDE,'R')) THEN 
         INFO = -1 
      ELSE IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = -2 
      ELSE IF (.NOT.NOTRAN .AND. .NOT.LSAME(TRANS,'T')) THEN 
         INFO = -3 
      ELSE IF (M < 0) THEN 
         INFO = -4 
      ELSE IF (N < 0) THEN 
         INFO = -5 
      ELSE IF (LDC < MAX(1,M)) THEN 
         INFO = -9 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DOPMTR', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (M==0 .OR. N==0) RETURN  
!
      IF (UPPER) THEN 
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
         FORWRD = LEFT .AND. NOTRAN .OR. .NOT.LEFT .AND. .NOT.NOTRAN 
!
         IF (FORWRD) THEN 
            I1 = 1 
            I2 = NQ - 1 
            I3 = 1 
            II = 2 
         ELSE 
            I1 = NQ - 1 
            I2 = 1 
            I3 = -1 
            II = NQ*(NQ + 1)/2 - 1 
         ENDIF 
!
         IF (LEFT) THEN 
            NI = N 
         ELSE 
            MI = M 
         ENDIF 
!
         DO I = I1, I2, I3 
            IF (LEFT) THEN 
!
!              H(i) is applied to C(1:i,1:n)
!
               MI = I 
            ELSE 
!
!              H(i) is applied to C(1:m,1:i)
!
               NI = I 
            ENDIF 
!
!           Apply H(i)
!
            AII = AP(II) 
            AP(II) = ONE 
            CALL DLARF (SIDE, MI, NI, AP(II-I+1), 1, TAU(I), C, LDC, WORK) 
            AP(II) = AII 
!
            IF (FORWRD) THEN 
               II = II + I + 2 
            ELSE 
               II = II - I - 1 
            ENDIF 
         END DO 
      ELSE 
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
         FORWRD = LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN 
!
         IF (FORWRD) THEN 
            I1 = 1 
            I2 = NQ - 1 
            I3 = 1 
            II = 2 
         ELSE 
            I1 = NQ - 1 
            I2 = 1 
            I3 = -1 
            II = NQ*(NQ + 1)/2 - 1 
         ENDIF 
!
         IF (LEFT) THEN 
            NI = N 
            JC = 1 
         ELSE 
            MI = M 
            IC = 1 
         ENDIF 
!
         DO I = I1, I2, I3 
            AII = AP(II) 
            AP(II) = ONE 
            IF (LEFT) THEN 
!
!              H(i) is applied to C(i+1:m,1:n)
!
               MI = M - I 
               IC = I + 1 
            ELSE 
!
!              H(i) is applied to C(1:m,i+1:n)
!
               NI = N - I 
               JC = I + 1 
            ENDIF 
!
!           Apply H(i)
!
            CALL DLARF (SIDE, MI, NI, AP(II), 1, TAU(I), C(IC,JC), LDC, WORK) 
            AP(II) = AII 
!
            IF (FORWRD) THEN 
               II = II + NQ - I + 1 
            ELSE 
               II = II - NQ + I - 2 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
!
!     End of DOPMTR
!
      END SUBROUTINE DOPMTR 
