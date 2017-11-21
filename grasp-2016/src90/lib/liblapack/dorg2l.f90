      SUBROUTINE DORG2L(M, N, K, A, LDA, TAU, WORK, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:34   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlarf_I 
      !USE dscal_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: K 
      INTEGER  :: LDA 
      INTEGER , INTENT(OUT) :: INFO 
      REAL(DOUBLE)  :: A(LDA,*) 
      REAL(DOUBLE)  :: TAU(*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, II, J, L 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DORG2L generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
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
      IF (M < 0) THEN 
         INFO = -1 
      ELSE IF (N<0 .OR. N>M) THEN 
         INFO = -2 
      ELSE IF (K<0 .OR. K>N) THEN 
         INFO = -3 
      ELSE IF (LDA < MAX(1,M)) THEN 
         INFO = -5 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DORG2L', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N <= 0) RETURN  
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      DO J = 1, N - K 
         A(:M,J) = ZERO 
         A(M-N+J,J) = ONE 
      END DO 
!
      DO I = 1, K 
         II = N - K + I 
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         A(M-N+II,II) = ONE 
         CALL DLARF ('Left', M - N + II, II - 1, A(1,II), 1, TAU(I), A, LDA, &
            WORK) 
         CALL DSCAL (M - N + II - 1, (-TAU(I)),A(1,II), 1) 
         A(M-N+II,II) = ONE - TAU(I) 
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         A(M-N+II+1:M,II) = ZERO 
      END DO 
      RETURN  
!
!     End of DORG2L
!
      END SUBROUTINE DORG2L 
