      SUBROUTINE DLAGTS(JOB, N, A, B, C, D, IN, Y, TOL, INFO) 
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
      INTEGER , INTENT(IN) :: JOB 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: INFO 
      REAL(DOUBLE) , INTENT(INOUT) :: TOL 
      INTEGER , INTENT(IN) :: IN(*) 
      REAL(DOUBLE) , INTENT(IN) :: A(*) 
      REAL(DOUBLE) , INTENT(IN) :: B(*) 
      REAL(DOUBLE) , INTENT(IN) :: C(*) 
      REAL(DOUBLE) , INTENT(IN) :: D(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: Y(*) 
      REAL(KIND(0.0D0)) :: dlamch
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K 
      REAL(DOUBLE) :: ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, SIGN 
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
!  DLAGTS may be used to solve one of the systems of equations
!
!     (T - lambda*I)*x = y   or   (T - lambda*I)'*x = y,
!
!  where T is an n by n tridiagonal matrix, for x, following the
!  factorization of (T - lambda*I) as
!
!     (T - lambda*I) = P*L*U ,
!
!  by routine DLAGTF. The choice of equation to be solved is
!  controlled by the argument JOB, and in each case there is an option
!  to perturb zero or very small diagonal elements of U, this option
!  being intended for use in applications such as inverse iteration.
!
!  Arguments
!  =========
!
!  JOB     (input) INTEGER
!          Specifies the job to be performed by DLAGTS as follows:
!          =  1: The equations  (T - lambda*I)x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -1: The equations  (T - lambda*I)x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!          =  2: The equations  (T - lambda*I)'x = y  are to be solved,
!                but diagonal elements of U are not to be perturbed.
!          = -2: The equations  (T - lambda*I)'x = y  are to be solved
!                and, if overflow would otherwise occur, the diagonal
!                elements of U are to be perturbed. See argument TOL
!                below.
!
!  N       (input) INTEGER
!          The order of the matrix T.
!
!  A       (input) DOUBLE PRECISION array, dimension (N)
!          On entry, A must contain the diagonal elements of U as
!          returned from DLAGTF.
!
!  B       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, B must contain the first super-diagonal elements of
!          U as returned from DLAGTF.
!
!  C       (input) DOUBLE PRECISION array, dimension (N-1)
!          On entry, C must contain the sub-diagonal elements of L as
!          returned from DLAGTF.
!
!  D       (input) DOUBLE PRECISION array, dimension (N-2)
!          On entry, D must contain the second super-diagonal elements
!          of U as returned from DLAGTF.
!
!  IN      (input) INTEGER array, dimension (N)
!          On entry, IN must contain details of the matrix P as returned
!          from DLAGTF.
!
!  Y       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the right hand side vector y.
!          On exit, Y is overwritten by the solution vector x.
!
!  TOL     (input/output) DOUBLE PRECISION
!          On entry, with  JOB .lt. 0, TOL should be the minimum
!          perturbation to be made to very small diagonal elements of U.
!          TOL should normally be chosen as about eps*norm(U), where eps
!          is the relative machine precision, but if TOL is supplied as
!          non-positive, then it is reset to eps*max( abs( u(i,j) ) ).
!          If  JOB .gt. 0  then TOL is not referenced.
!
!          On exit, TOL is changed as described above, only if TOL is
!          non-positive on entry. Otherwise TOL is unchanged.
!
!  INFO    (output) INTEGER
!          = 0   : successful exit
!          .lt. 0: if INFO = -i, the i-th argument had an illegal value
!          .gt. 0: overflow would occur when computing the INFO(th)
!                  element of the solution vector x. This can only occur
!                  when JOB is supplied as positive and either means
!                  that a diagonal element of U is very small, or that
!                  the elements of the right-hand side vector y are very
!                  large.
!
!  =====================================================================
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
      IF (ABS(JOB)>2 .OR. JOB==0) THEN 
         INFO = -1 
      ELSE IF (N < 0) THEN 
         INFO = -2 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DLAGTS', (-INFO)) 
         RETURN  
      ENDIF 
!
      IF (N == 0) RETURN  
!
      EPS = DLAMCH('Epsilon') 
      SFMIN = DLAMCH('Safe minimum') 
      BIGNUM = ONE/SFMIN 
!
      IF (JOB < 0) THEN 
         IF (TOL <= ZERO) THEN 
            TOL = ABS(A(1)) 
            IF (N > 1) TOL = MAX(MAX(TOL,ABS(A(2))),ABS(B(1))) 
            DO K = 3, N 
               TOL = MAX(MAX(MAX(TOL,ABS(A(K))),ABS(B(K-1))),ABS(D(K-2))) 
            END DO 
            TOL = TOL*EPS 
            IF (TOL == ZERO) TOL = EPS 
         ENDIF 
      ENDIF 
!
      IF (ABS(JOB) == 1) THEN 
         DO K = 2, N 
            IF (IN(K-1) == 0) THEN 
               Y(K) = Y(K) - C(K-1)*Y(K-1) 
            ELSE 
               TEMP = Y(K-1) 
               Y(K-1) = Y(K) 
               Y(K) = TEMP - C(K-1)*Y(K) 
            ENDIF 
         END DO 
         IF (JOB == 1) THEN 
            DO K = N, 1, -1 
               IF (K <= N - 2) THEN 
                  TEMP = Y(K) - B(K)*Y(K+1) - D(K)*Y(K+2) 
               ELSE IF (K == N - 1) THEN 
                  TEMP = Y(K) - B(K)*Y(K+1) 
               ELSE 
                  TEMP = Y(K) 
               ENDIF 
               AK = A(K) 
               ABSAK = ABS(AK) 
               IF (ABSAK < ONE) THEN 
                  IF (ABSAK < SFMIN) THEN 
                     IF (ABSAK==ZERO .OR. ABS(TEMP)*SFMIN>ABSAK) THEN 
                        INFO = K 
                        RETURN  
                     ELSE 
                        TEMP = TEMP*BIGNUM 
                        AK = AK*BIGNUM 
                     ENDIF 
                  ELSE IF (ABS(TEMP) > ABSAK*BIGNUM) THEN 
                     INFO = K 
                     RETURN  
                  ENDIF 
               ENDIF 
               Y(K) = TEMP/AK 
            END DO 
         ELSE 
            DO K = N, 1, -1 
               IF (K <= N - 2) THEN 
                  TEMP = Y(K) - B(K)*Y(K+1) - D(K)*Y(K+2) 
               ELSE IF (K == N - 1) THEN 
                  TEMP = Y(K) - B(K)*Y(K+1) 
               ELSE 
                  TEMP = Y(K) 
               ENDIF 
               AK = A(K) 
               PERT = SIGN(TOL,AK) 
   40          CONTINUE 
               ABSAK = ABS(AK) 
               IF (ABSAK < ONE) THEN 
                  IF (ABSAK < SFMIN) THEN 
                     IF (ABSAK==ZERO .OR. ABS(TEMP)*SFMIN>ABSAK) THEN 
                        AK = AK + PERT 
                        PERT = 2*PERT 
                        GO TO 40 
                     ELSE 
                        TEMP = TEMP*BIGNUM 
                        AK = AK*BIGNUM 
                     ENDIF 
                  ELSE IF (ABS(TEMP) > ABSAK*BIGNUM) THEN 
                     AK = AK + PERT 
                     PERT = 2*PERT 
                     GO TO 40 
                  ENDIF 
               ENDIF 
               Y(K) = TEMP/AK 
            END DO 
         ENDIF 
      ELSE 
!
!        Come to here if  JOB = 2 or -2
!
         IF (JOB == 2) THEN 
            DO K = 1, N 
               IF (K >= 3) THEN 
                  TEMP = Y(K) - B(K-1)*Y(K-1) - D(K-2)*Y(K-2) 
               ELSE IF (K == 2) THEN 
                  TEMP = Y(K) - B(K-1)*Y(K-1) 
               ELSE 
                  TEMP = Y(K) 
               ENDIF 
               AK = A(K) 
               ABSAK = ABS(AK) 
               IF (ABSAK < ONE) THEN 
                  IF (ABSAK < SFMIN) THEN 
                     IF (ABSAK==ZERO .OR. ABS(TEMP)*SFMIN>ABSAK) THEN 
                        INFO = K 
                        RETURN  
                     ELSE 
                        TEMP = TEMP*BIGNUM 
                        AK = AK*BIGNUM 
                     ENDIF 
                  ELSE IF (ABS(TEMP) > ABSAK*BIGNUM) THEN 
                     INFO = K 
                     RETURN  
                  ENDIF 
               ENDIF 
               Y(K) = TEMP/AK 
            END DO 
         ELSE 
            DO K = 1, N 
               IF (K >= 3) THEN 
                  TEMP = Y(K) - B(K-1)*Y(K-1) - D(K-2)*Y(K-2) 
               ELSE IF (K == 2) THEN 
                  TEMP = Y(K) - B(K-1)*Y(K-1) 
               ELSE 
                  TEMP = Y(K) 
               ENDIF 
               AK = A(K) 
               PERT = SIGN(TOL,AK) 
   70          CONTINUE 
               ABSAK = ABS(AK) 
               IF (ABSAK < ONE) THEN 
                  IF (ABSAK < SFMIN) THEN 
                     IF (ABSAK==ZERO .OR. ABS(TEMP)*SFMIN>ABSAK) THEN 
                        AK = AK + PERT 
                        PERT = 2*PERT 
                        GO TO 70 
                     ELSE 
                        TEMP = TEMP*BIGNUM 
                        AK = AK*BIGNUM 
                     ENDIF 
                  ELSE IF (ABS(TEMP) > ABSAK*BIGNUM) THEN 
                     AK = AK + PERT 
                     PERT = 2*PERT 
                     GO TO 70 
                  ENDIF 
               ENDIF 
               Y(K) = TEMP/AK 
            END DO 
         ENDIF 
!
         DO K = N, 2, -1 
            IF (IN(K-1) == 0) THEN 
               Y(K-1) = Y(K-1) - C(K-1)*Y(K) 
            ELSE 
               TEMP = Y(K-1) 
               Y(K-1) = Y(K) 
               Y(K) = TEMP - C(K-1)*Y(K) 
            ENDIF 
         END DO 
      ENDIF 
      RETURN  
!
!     End of DLAGTS
!
      END SUBROUTINE DLAGTS 
