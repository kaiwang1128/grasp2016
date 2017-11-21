      SUBROUTINE DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:36   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USE dlamch_I 
      !USE dlanst_I 
      !USE dlapy2_I 
      !USE dlae2_I 
      !USE dlaev2_I 
      !USE dlartg_I 
      !USE dlascl_I 
      !USE dlaset_I 
      !USE dlasr_I 
      !USE dlasrt_I 
      !USE dswap_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER  :: LDZ 
      INTEGER  :: INFO 
      CHARACTER  :: COMPZ 
      REAL(DOUBLE)  :: D(*) 
      REAL(DOUBLE)  :: E(*) 
      REAL(DOUBLE)  :: Z(LDZ,*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D0 
      REAL(DOUBLE), PARAMETER :: THREE = 3.0D0 
      INTEGER, PARAMETER :: MAXIT = 30 
      REAL(KIND(0.0D0))::dlapy2
      REAL(KIND(0.0D0)) :: dlamch
      logical::LSAME
      REAL(KIND(0.0D0)) :: DLANST
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, LENDM1, LENDP1&
         , LENDSV, LM1, LSV, M, MM, MM1, NM1, NMAXIT 
      REAL(DOUBLE) :: ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, S, SAFMAX, &
         SAFMIN, SSFMAX, SSFMIN, TST 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, SIGN, SQRT 
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
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
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
      IF (LSAME(COMPZ,'N')) THEN 
         ICOMPZ = 0 
      ELSE IF (LSAME(COMPZ,'V')) THEN 
         ICOMPZ = 1 
      ELSE IF (LSAME(COMPZ,'I')) THEN 
         ICOMPZ = 2 
      ELSE 
         ICOMPZ = -1 
      ENDIF 
      IF (ICOMPZ < 0) THEN 
         INFO = -1 
      ELSE IF (N < 0) THEN 
         INFO = -2 
      ELSE IF (LDZ<1 .OR. ICOMPZ>0 .AND. LDZ<MAX(1,N)) THEN 
         INFO = -6 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSTEQR', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N == 0) RETURN  
!
      IF (N == 1) THEN 
         IF (ICOMPZ == 2) Z(1,1) = ONE 
         RETURN  
      ENDIF 
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH('E') 
      EPS2 = EPS**2 
      SAFMIN = DLAMCH('S') 
      SAFMAX = ONE/SAFMIN 
      SSFMAX = SQRT(SAFMAX)/THREE 
      SSFMIN = SQRT(SAFMIN)/EPS2 
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF (ICOMPZ == 2) CALL DLASET ('Full', N, N, ZERO, ONE, Z, LDZ) 
!
      NMAXIT = N*MAXIT 
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
      IF (L1 > N) GO TO 160 
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
      IF (ANORM == ZERO) GO TO 10 
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
!     Choose between QL and QR iteration
!
      IF (ABS(D(LEND)) < ABS(D(L))) THEN 
         LEND = LSV 
         L = LENDSV 
      ENDIF 
!
      IF (LEND > L) THEN 
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE 
         IF (L /= LEND) THEN 
            LENDM1 = LEND - 1 
            DO M = L, LENDM1 
               TST = ABS(E(M))**2 
               IF (TST <= (EPS2*ABS(D(M)))*ABS(D(M+1)) + SAFMIN) GO TO 60 
            END DO 
         ENDIF 
!
         M = LEND 
!
   60    CONTINUE 
         IF (M < LEND) E(M) = ZERO 
         P = D(L) 
         IF (M == L) GO TO 80 
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF (M == L + 1) THEN 
            IF (ICOMPZ > 0) THEN 
               CALL DLAEV2 (D(L), E(L), D(L+1), RT1, RT2, C, S) 
               WORK(L) = C 
               WORK(N-1+L) = S 
               CALL DLASR ('R', 'V', 'B', N, 2, WORK(L), WORK(N-1+L), Z(1,L), &
                  LDZ) 
            ELSE 
               CALL DLAE2 (D(L), E(L), D(L+1), RT1, RT2) 
            ENDIF 
            D(L) = RT1 
            D(L+1) = RT2 
            E(L) = ZERO 
            L = L + 2 
            IF (L <= LEND) GO TO 40 
            GO TO 140 
         ENDIF 
!
         IF (JTOT == NMAXIT) GO TO 140 
         JTOT = JTOT + 1 
!
!        Form shift.
!
         G = (D(L+1)-P)/(TWO*E(L)) 
         R = DLAPY2(G,ONE) 
         G = D(M) - P + E(L)/(G + SIGN(R,G)) 
!
         S = ONE 
         C = ONE 
         P = ZERO 
!
!        Inner loop
!
         MM1 = M - 1 
         DO I = MM1, L, -1 
            F = S*E(I) 
            B = C*E(I) 
            CALL DLARTG (G, F, C, S, R) 
            IF (I /= M - 1) E(I+1) = R 
            G = D(I+1) - P 
            R = (D(I)-G)*S + TWO*C*B 
            P = S*R 
            D(I+1) = G + P 
            G = C*R - B 
!
!           If eigenvectors are desired, then save rotations.
!
            IF (ICOMPZ <= 0) CYCLE  
            WORK(I) = C 
            WORK(N-1+I) = -S 
!
         END DO 
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF (ICOMPZ > 0) THEN 
            MM = M - L + 1 
            CALL DLASR ('R', 'V', 'B', N, MM, WORK(L), WORK(N-1+L), Z(1,L), LDZ&
               ) 
         ENDIF 
!
         D(L) = D(L) - P 
         E(L) = G 
         GO TO 40 
!
!        Eigenvalue found.
!
   80    CONTINUE 
         D(L) = P 
!
         L = L + 1 
         IF (L <= LEND) GO TO 40 
         GO TO 140 
!
      ELSE 
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE 
         IF (L /= LEND) THEN 
            LENDP1 = LEND + 1 
            DO M = L, LENDP1, -1 
               TST = ABS(E(M-1))**2 
               IF (TST <= (EPS2*ABS(D(M)))*ABS(D(M-1)) + SAFMIN) GO TO 110 
            END DO 
         ENDIF 
!
         M = LEND 
!
  110    CONTINUE 
         IF (M > LEND) E(M-1) = ZERO 
         P = D(L) 
         IF (M == L) GO TO 130 
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF (M == L - 1) THEN 
            IF (ICOMPZ > 0) THEN 
               CALL DLAEV2 (D(L-1), E(L-1), D(L), RT1, RT2, C, S) 
               WORK(M) = C 
               WORK(N-1+M) = S 
               CALL DLASR ('R', 'V', 'F', N, 2, WORK(M), WORK(N-1+M), Z(1,L-1)&
                  , LDZ) 
            ELSE 
               CALL DLAE2 (D(L-1), E(L-1), D(L), RT1, RT2) 
            ENDIF 
            D(L-1) = RT1 
            D(L) = RT2 
            E(L-1) = ZERO 
            L = L - 2 
            IF (L >= LEND) GO TO 90 
            GO TO 140 
         ENDIF 
!
         IF (JTOT == NMAXIT) GO TO 140 
         JTOT = JTOT + 1 
!
!        Form shift.
!
         G = (D(L-1)-P)/(TWO*E(L-1)) 
         R = DLAPY2(G,ONE) 
         G = D(M) - P + E(L-1)/(G + SIGN(R,G)) 
!
         S = ONE 
         C = ONE 
         P = ZERO 
!
!        Inner loop
!
         LM1 = L - 1 
         DO I = M, LM1 
            F = S*E(I) 
            B = C*E(I) 
            CALL DLARTG (G, F, C, S, R) 
            IF (I /= M) E(I-1) = R 
            G = D(I) - P 
            R = (D(I+1)-G)*S + TWO*C*B 
            P = S*R 
            D(I) = G + P 
            G = C*R - B 
!
!           If eigenvectors are desired, then save rotations.
!
            IF (ICOMPZ <= 0) CYCLE  
            WORK(I) = C 
            WORK(N-1+I) = S 
!
         END DO 
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF (ICOMPZ > 0) THEN 
            MM = L - M + 1 
            CALL DLASR ('R', 'V', 'F', N, MM, WORK(M), WORK(N-1+M), Z(1,M), LDZ&
               ) 
         ENDIF 
!
         D(L) = D(L) - P 
         E(LM1) = G 
         GO TO 90 
!
!        Eigenvalue found.
!
  130    CONTINUE 
         D(L) = P 
!
         L = L - 1 
         IF (L >= LEND) GO TO 90 
         GO TO 140 
!
      ENDIF 
!
!     Undo scaling if necessary
!
  140 CONTINUE 
      IF (ISCALE == 1) THEN 
         CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV + 1, 1, D(LSV), N&
            , INFO) 
         CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV - LSV, 1, E(LSV), N, &
            INFO) 
      ELSE IF (ISCALE == 2) THEN 
         CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV + 1, 1, D(LSV), N&
            , INFO) 
         CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV - LSV, 1, E(LSV), N, &
            INFO) 
      ENDIF 
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF (JTOT < NMAXIT) GO TO 10 
      INFO = INFO + COUNT(E(:N-1)/=ZERO) 
      GO TO 190 
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE 
      IF (ICOMPZ == 0) THEN 
!
!        Use Quick Sort
!
         CALL DLASRT ('I', N, D, INFO) 
!
      ELSE 
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO II = 2, N 
            I = II - 1 
            K = I 
            P = D(I) 
            DO J = II, N 
               IF (D(J) >= P) CYCLE  
               K = J 
               P = D(J) 
            END DO 
            IF (K == I) CYCLE  
            D(K) = D(I) 
            D(I) = P 
            CALL DSWAP (N, Z(1,I), 1, Z(1,K), 1) 
         END DO 
      ENDIF 
!
  190 CONTINUE 
      RETURN  
!
!     End of DSTEQR
!
      END SUBROUTINE DSTEQR 
