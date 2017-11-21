      SUBROUTINE DOPGTR(UPLO, N, AP, TAU, Q, LDQ, WORK, INFO) 
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
      !USE dorg2l_I 
      !USE dorg2r_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER  :: LDQ 
      INTEGER , INTENT(OUT) :: INFO 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: AP(*) 
      REAL(DOUBLE)  :: TAU(*) 
      REAL(DOUBLE)  :: Q(LDQ,*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IINFO, IJ, J 
      LOGICAL :: UPPER 
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
!  DOPGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  DSPTRD using packed storage:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
!          The N-by-N orthogonal matrix Q.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N-1)
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
      UPPER = LSAME(UPLO,'U') 
      IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = -1 
      ELSE IF (N < 0) THEN 
         INFO = -2 
      ELSE IF (LDQ < MAX(1,N)) THEN 
         INFO = -6 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DOPGTR', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N == 0) RETURN  
!
      IF (UPPER) THEN 
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
         IJ = 2 
         DO J = 1, N - 1 
            IF (J - 1 > 0) THEN 
               Q(:J-1,J) = AP(IJ:J-2+IJ) 
               IJ = J - 1 + IJ 
            ENDIF 
            IJ = IJ + 2 
            Q(N,J) = ZERO 
         END DO 
         Q(:N-1,N) = ZERO 
         Q(N,N) = ONE 
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL DORG2L (N - 1, N - 1, N - 1, Q, LDQ, TAU, WORK, IINFO) 
!
      ELSE 
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
         Q(1,1) = ONE 
         Q(2:N,1) = ZERO 
         IJ = 3 
         DO J = 2, N 
            Q(1,J) = ZERO 
            IF (N - J > 0) THEN 
               Q(J+1:N,J) = AP(IJ:N-J-1+IJ) 
               IJ = N - J + IJ 
            ENDIF 
            IJ = IJ + 2 
         END DO 
!
!           Generate Q(2:n,2:n)
!
         IF (N > 1) CALL DORG2R (N - 1, N - 1, N - 1, Q(2,2), LDQ, TAU, WORK, &
            IINFO) 
      ENDIF 
      RETURN  
!
!     End of DOPGTR
!
      END SUBROUTINE DOPGTR 
