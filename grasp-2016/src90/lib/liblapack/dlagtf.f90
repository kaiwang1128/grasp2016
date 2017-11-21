      SUBROUTINE DLAGTF(N, A, LAMBDA, B, C, TOL, D, IN, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:30   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlamch_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: INFO 
      REAL(DOUBLE) , INTENT(IN) :: LAMBDA 
      REAL(DOUBLE) , INTENT(IN) :: TOL 
      INTEGER , INTENT(INOUT) :: IN(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: A(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: B(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: C(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: D(*) 
      REAL(KIND(0.0D0)) :: dlamch
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K 
      REAL(DOUBLE) :: EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX 
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
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
!  DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n
!  tridiagonal matrix and lambda is a scalar, as
!
!     T - lambda*I = PLU,
!
!  where P is a permutation matrix, L is a unit lower tridiagonal matrix
!  with at most one non-zero sub-diagonal elements per column and U is
!  an upper triangular matrix with at most two non-zero super-diagonal
!  elements per column.
!
!  The factorization is obtained by Gaussian elimination with partial
!  pivoting and implicit row scaling.
!
!  The parameter LAMBDA is included in the routine so that DLAGTF may
!  be used, in conjunction with DLAGTS, to obtain eigenvectors of T by
!  inverse iteration.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix T.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, A must contain the diagonal elements of T.
!
!          On exit, A is overwritten by the n diagonal elements of the
!          upper triangular matrix U of the factorization of T.
!
!  LAMBDA  (input) DOUBLE PRECISION
!          On entry, the scalar lambda.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, B must contain the (n-1) super-diagonal elements of
!          T.
!
!          On exit, B is overwritten by the (n-1) super-diagonal
!          elements of the matrix U of the factorization of T.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, C must contain the (n-1) sub-diagonal elements of
!          T.
!
!          On exit, C is overwritten by the (n-1) sub-diagonal elements
!          of the matrix L of the factorization of T.
!
!  TOL     (input) DOUBLE PRECISION
!          On entry, a relative tolerance used to indicate whether or
!          not the matrix (T - lambda*I) is nearly singular. TOL should
!          normally be chose as approximately the largest relative error
!          in the elements of T. For example, if the elements of T are
!          correct to about 4 significant figures, then TOL should be
!          set to about 5*10**(-4). If TOL is supplied as less than eps,
!          where eps is the relative machine precision, then the value
!          eps is used in place of TOL.
!
!  D       (output) DOUBLE PRECISION array, dimension (N-2)
!          On exit, D is overwritten by the (n-2) second super-diagonal
!          elements of the matrix U of the factorization of T.
!
!  IN      (output) INTEGER array, dimension (N)
!          On exit, IN contains details of the permutation matrix P. If
!          an interchange occurred at the kth step of the elimination,
!          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)
!          returns the smallest positive integer j such that
!
!             abs( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,
!
!          where norm( A(j) ) denotes the sum of the absolute values of
!          the jth row of the matrix A. If no such j exists then IN(n)
!          is returned as zero. If IN(n) is returned as positive, then a
!          diagonal element of U is small, indicating that
!          (T - lambda*I) is singular or nearly singular,
!
!  INFO    (output)
!          = 0   : successful exit
!          .lt. 0: if INFO = -k, the kth argument had an illegal value
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      INFO = 0 
      IF (N < 0) THEN 
         INFO = -1 
         CALL XERBLA ('DLAGTF', (-INFO)) 
         RETURN  
      ENDIF 
!
      IF (N == 0) RETURN  
!
      A(1) = A(1) - LAMBDA 
      IN(N) = 0 
      IF (N == 1) THEN 
         IF (A(1) == ZERO) IN(1) = 1 
         RETURN  
      ENDIF 
!
      EPS = DLAMCH('Epsilon') 
!
      TL = MAX(TOL,EPS) 
      SCALE1 = ABS(A(1)) + ABS(B(1)) 
      DO K = 1, N - 1 
         A(K+1) = A(K+1) - LAMBDA 
         SCALE2 = ABS(C(K)) + ABS(A(K+1)) 
         IF (K < N - 1) SCALE2 = SCALE2 + ABS(B(K+1)) 
         IF (A(K) == ZERO) THEN 
            PIV1 = ZERO 
         ELSE 
            PIV1 = ABS(A(K))/SCALE1 
         ENDIF 
         IF (C(K) == ZERO) THEN 
            IN(K) = 0 
            PIV2 = ZERO 
            SCALE1 = SCALE2 
            IF (K < N - 1) D(K) = ZERO 
         ELSE 
            PIV2 = ABS(C(K))/SCALE2 
            IF (PIV2 <= PIV1) THEN 
               IN(K) = 0 
               SCALE1 = SCALE2 
               C(K) = C(K)/A(K) 
               A(K+1) = A(K+1) - C(K)*B(K) 
               IF (K < N - 1) D(K) = ZERO 
            ELSE 
               IN(K) = 1 
               MULT = A(K)/C(K) 
               A(K) = C(K) 
               TEMP = A(K+1) 
               A(K+1) = B(K) - MULT*TEMP 
               IF (K < N - 1) THEN 
                  D(K) = B(K+1) 
                  B(K+1) = -MULT*D(K) 
               ENDIF 
               B(K) = TEMP 
               C(K) = MULT 
            ENDIF 
         ENDIF 
         IF (MAX(PIV1,PIV2)>TL .OR. IN(N)/=0) CYCLE  
         IN(N) = K 
      END DO 
      IF (ABS(A(N))<=SCALE1*TL .AND. IN(N)==0) IN(N) = N 
!
      RETURN  
!
!     End of DLAGTF
!
      END SUBROUTINE DLAGTF 
