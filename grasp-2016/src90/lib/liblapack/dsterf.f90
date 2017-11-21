      SUBROUTINE DSTERF(N, D, E, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:37   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlamch_I 
      !USE dlanst_I 
      !USE dlapy2_I 
      !USE dlae2_I 
      !USE dlascl_I 
      !USE dlasrt_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER  :: INFO 
      REAL(DOUBLE)  :: D(*) 
      REAL(DOUBLE)  :: E(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D0 
      REAL(DOUBLE), PARAMETER :: THREE = 3.0D0 
      INTEGER, PARAMETER :: MAXIT = 30 
      REAL(KIND(0.0D0)) :: dlapy2
      REAL(KIND(0.0D0)) :: dlamch
      REAL(KIND(0.0D0)) :: DLANST
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ISCALE, JTOT, L, L1, LEND, LENDM1, LENDP1, LENDSV, LM1, LSV&
         , M, MM1, NM1, NMAXIT 
      REAL(DOUBLE) :: ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, OLDGAM, P, R&
         , RT1, RT2, RTE, S, SAFMAX, SAFMIN, SIGMA, SSFMAX, SSFMIN, TST 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, SIGN, SQRT 
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
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
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
!     Test the input parameters.
!
      INFO = 0 
!
!     Quick return if possible
!
      IF (N < 0) THEN 
         INFO = -1 
         CALL XERBLA ('DSTERF', (-INFO)) 
         RETURN  
      ENDIF 
      IF (N <= 1) RETURN  
!
!     Determine the unit roundoff for this environment.
!
      EPS = DLAMCH('E') 
      EPS2 = EPS**2 
      SAFMIN = DLAMCH('S') 
      SAFMAX = ONE/SAFMIN 
      SSFMAX = SQRT(SAFMAX)/THREE 
      SSFMIN = SQRT(SAFMIN)/EPS2 
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      NMAXIT = N*MAXIT 
      SIGMA = ZERO 
      JTOT = 0 
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1 
      NM1 = N - 1 
!
   10 CONTINUE 
      IF (L1 > N) GO TO 170 
      IF (L1 > 1) E(L1-1) = ZERO 
      IF (L1 <= NM1) THEN 
         DO M = L1, NM1 
            TST = ABS(E(M)) 
            IF (TST == ZERO) GO TO 30 
            IF (TST > (SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))*EPS) CYCLE  
            E(M) = ZERO 
            GO TO 30 
         END DO 
      ENDIF 
      M = N 
!
   30 CONTINUE 
      L = L1 
      LSV = L 
      LEND = M 
      LENDSV = LEND 
      L1 = M + 1 
      IF (LEND == L) GO TO 10 
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST('I',LEND - L + 1,D(L),E(L)) 
      ISCALE = 0 
      IF (ANORM > SSFMAX) THEN 
         ISCALE = 1 
         CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND - L + 1, 1, D(L), N, INFO) 
         CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND - L, 1, E(L), N, INFO) 
      ELSE IF (ANORM < SSFMIN) THEN 
         ISCALE = 2 
         CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND - L + 1, 1, D(L), N, INFO) 
         CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND - L, 1, E(L), N, INFO) 
      ENDIF 
!
      E(L:LEND-1) = E(L:LEND-1)**2 
!
!     Choose between QL and QR iteration
!
      IF (ABS(D(LEND)) < ABS(D(L))) THEN 
         LEND = LSV 
         L = LENDSV 
      ENDIF 
!
      IF (LEND >= L) THEN 
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   50    CONTINUE 
         IF (L /= LEND) THEN 
            LENDM1 = LEND - 1 
            DO M = L, LENDM1 
               TST = ABS(E(M)) 
               IF (TST <= EPS2*ABS(D(M)*D(M+1))) GO TO 70 
            END DO 
         ENDIF 
!
         M = LEND 
!
   70    CONTINUE 
         IF (M < LEND) E(M) = ZERO 
         P = D(L) 
         IF (M == L) GO TO 90 
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF (M == L + 1) THEN 
            RTE = SQRT(E(L)) 
            CALL DLAE2 (D(L), RTE, D(L+1), RT1, RT2) 
            D(L) = RT1 
            D(L+1) = RT2 
            E(L) = ZERO 
            L = L + 2 
            IF (L <= LEND) GO TO 50 
            GO TO 150 
         ENDIF 
!
         IF (JTOT == NMAXIT) GO TO 150 
         JTOT = JTOT + 1 
!
!        Form shift.
!
         RTE = SQRT(E(L)) 
         SIGMA = (D(L+1)-P)/(TWO*RTE) 
         R = DLAPY2(SIGMA,ONE) 
         SIGMA = P - RTE/(SIGMA + SIGN(R,SIGMA)) 
!
         C = ONE 
         S = ZERO 
         GAMMA = D(M) - SIGMA 
         P = GAMMA*GAMMA 
!
!        Inner loop
!
         MM1 = M - 1 
         DO I = MM1, L, -1 
            BB = E(I) 
            R = P + BB 
            IF (I /= M - 1) E(I+1) = S*R 
            OLDC = C 
            C = P/R 
            S = BB/R 
            OLDGAM = GAMMA 
            ALPHA = D(I) 
            GAMMA = C*(ALPHA - SIGMA) - S*OLDGAM 
            D(I+1) = OLDGAM + (ALPHA - GAMMA) 
            IF (C /= ZERO) THEN 
               P = (GAMMA*GAMMA)/C 
            ELSE 
               P = OLDC*BB 
            ENDIF 
         END DO 
!
         E(L) = S*P 
         D(L) = SIGMA + GAMMA 
         GO TO 50 
!
!        Eigenvalue found.
!
   90    CONTINUE 
         D(L) = P 
!
         L = L + 1 
         IF (L <= LEND) GO TO 50 
         GO TO 150 
!
      ELSE 
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
  100    CONTINUE 
         IF (L /= LEND) THEN 
            LENDP1 = LEND + 1 
            DO M = L, LENDP1, -1 
               TST = ABS(E(M-1)) 
               IF (TST <= EPS2*ABS(D(M)*D(M-1))) GO TO 120 
            END DO 
         ENDIF 
!
         M = LEND 
!
  120    CONTINUE 
         IF (M > LEND) E(M-1) = ZERO 
         P = D(L) 
         IF (M == L) GO TO 140 
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF (M == L - 1) THEN 
            RTE = SQRT(E(L-1)) 
            CALL DLAE2 (D(L), RTE, D(L-1), RT1, RT2) 
            D(L) = RT1 
            D(L-1) = RT2 
            E(L-1) = ZERO 
            L = L - 2 
            IF (L >= LEND) GO TO 100 
            GO TO 150 
         ENDIF 
!
         IF (JTOT == NMAXIT) GO TO 150 
         JTOT = JTOT + 1 
!
!        Form shift.
!
         RTE = SQRT(E(L-1)) 
         SIGMA = (D(L-1)-P)/(TWO*RTE) 
         R = DLAPY2(SIGMA,ONE) 
         SIGMA = P - RTE/(SIGMA + SIGN(R,SIGMA)) 
!
         C = ONE 
         S = ZERO 
         GAMMA = D(M) - SIGMA 
         P = GAMMA*GAMMA 
!
!        Inner loop
!
         LM1 = L - 1 
         DO I = M, LM1 
            BB = E(I) 
            R = P + BB 
            IF (I /= M) E(I-1) = S*R 
            OLDC = C 
            C = P/R 
            S = BB/R 
            OLDGAM = GAMMA 
            ALPHA = D(I+1) 
            GAMMA = C*(ALPHA - SIGMA) - S*OLDGAM 
            D(I) = OLDGAM + (ALPHA - GAMMA) 
            IF (C /= ZERO) THEN 
               P = (GAMMA*GAMMA)/C 
            ELSE 
               P = OLDC*BB 
            ENDIF 
         END DO 
!
         E(LM1) = S*P 
         D(L) = SIGMA + GAMMA 
         GO TO 100 
!
!        Eigenvalue found.
!
  140    CONTINUE 
         D(L) = P 
!
         L = L - 1 
         IF (L >= LEND) GO TO 100 
         GO TO 150 
!
      ENDIF 
!
!     Undo scaling if necessary
!
  150 CONTINUE 
      IF (ISCALE == 1) CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV + 1&
         , 1, D(LSV), N, INFO) 
      IF (ISCALE == 2) CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV + 1&
         , 1, D(LSV), N, INFO) 
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF (JTOT == NMAXIT) THEN 
         INFO = INFO + COUNT(E(:N-1)/=ZERO) 
         RETURN  
      ENDIF 
      GO TO 10 
!
!     Sort eigenvalues in increasing order.
!
  170 CONTINUE 
      CALL DLASRT ('I', N, D, INFO) 
!
      RETURN  
!
!     End of DSTERF
!
      END SUBROUTINE DSTERF 
