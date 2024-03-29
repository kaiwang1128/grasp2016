      SUBROUTINE DSPTRD(UPLO, N, AP, D, E, TAU, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:36   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE daxpy_I 
      !USE dlarfg_I 
      !USE dspmv_I 
      !USE dspr2_I 
      !!USE xerbla_I 
      !!USE lsame_I 
      !USE ddot_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: INFO 
      CHARACTER  :: UPLO 
      REAL(DOUBLE)  :: AP(*) 
      REAL(DOUBLE) , INTENT(OUT) :: D(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: E(*) 
      REAL(DOUBLE)  :: TAU(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: HALF = 1.0D0/2.0D0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, I1, I1I1, II 
      REAL(DOUBLE) :: ALPHA, TAUI 
      LOGICAL :: UPPER 
      logical::lsame
      real(kind(0.0d0))::ddot
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSPTRD reduces a real symmetric matrix A stored in packed form to
!  symmetric tridiagonal form T by an orthogonal similarity
!  transformation: Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0 
      UPPER = LSAME(UPLO,'U') 
      IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = -1 
      ELSE IF (N < 0) THEN 
         INFO = -2 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSPTRD', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N <= 0) RETURN  
!
      IF (UPPER) THEN 
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
         I1 = N*(N - 1)/2 + 1 
         DO I = N - 1, 1, -1 
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            CALL DLARFG (I, AP(I1+I-1), AP(I1), 1, TAUI) 
            E(I) = AP(I1+I-1) 
!
            IF (TAUI /= ZERO) THEN 
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               AP(I1+I-1) = ONE 
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
               CALL DSPMV (UPLO, I, TAUI, AP, AP(I1), 1, ZERO, TAU, 1) 
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT(I,TAU,1,AP(I1),1) 
               CALL DAXPY (I, ALPHA, AP(I1), 1, TAU, 1) 
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSPR2 (UPLO, I, (-ONE), AP(I1), 1, TAU, 1, AP) 
!
               AP(I1+I-1) = E(I) 
            ENDIF 
            D(I+1) = AP(I1+I) 
            TAU(I) = TAUI 
            I1 = I1 - I 
         END DO 
         D(1) = AP(1) 
      ELSE 
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
         II = 1 
         DO I = 1, N - 1 
            I1I1 = II + N - I + 1 
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            CALL DLARFG (N - I, AP(II+1), AP(II+2), 1, TAUI) 
            E(I) = AP(II+1) 
!
            IF (TAUI /= ZERO) THEN 
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               AP(II+1) = ONE 
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
               CALL DSPMV (UPLO, N - I, TAUI, AP(I1I1), AP(II+1), 1, ZERO, TAU(&
                  I), 1) 
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT(N - I,TAU(I),1,AP(II+1),1) 
               CALL DAXPY (N - I, ALPHA, AP(II+1), 1, TAU(I), 1) 
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSPR2 (UPLO, N - I, (-ONE), AP(II+1), 1, TAU(I), 1, AP(I1I1&
                  )) 
!
               AP(II+1) = E(I) 
            ENDIF 
            D(I) = AP(II) 
            TAU(I) = TAUI 
            II = I1I1 
         END DO 
         D(N) = AP(II) 
      ENDIF 
!
      RETURN  
!
!     End of DSPTRD
!
      END SUBROUTINE DSPTRD 
