      SUBROUTINE DSTEBZ(RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, M, &
         NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, INFO) 
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
      !USE ilaenv_I 
      !USE dlamch_I 
      !USE dlaebz_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER , INTENT(IN) :: IL 
      INTEGER , INTENT(IN) :: IU 
      INTEGER , INTENT(OUT) :: M 
      INTEGER , INTENT(OUT) :: NSPLIT 
      INTEGER , INTENT(OUT) :: INFO 
      REAL(DOUBLE) , INTENT(IN) :: VL 
      REAL(DOUBLE) , INTENT(IN) :: VU 
      REAL(DOUBLE) , INTENT(IN) :: ABSTOL 
      CHARACTER  :: RANGE 
      CHARACTER  :: ORDER 
      INTEGER  :: IBLOCK(*) 
      INTEGER , INTENT(INOUT) :: ISPLIT(*) 
      INTEGER  :: IWORK(*) 
      REAL(DOUBLE)  :: D(*) 
      REAL(DOUBLE)  :: E(*) 
      REAL(DOUBLE)  :: W(*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D0 
      REAL(DOUBLE), PARAMETER :: HALF = 1.0D0/TWO 
      REAL(DOUBLE), PARAMETER :: FUDGE = 2.0D0 
      REAL(DOUBLE), PARAMETER :: RELFAC = 2.0D0 
      REAL(KIND(0.0D0)) :: dlamch
      INTEGER :: ILAENV
      logical::lsame
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM, IN, IOFF, &
         IORDER, IOUT, IRANGE, ITMAX, ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, &
         NWL, NWU 
      INTEGER, DIMENSION(1) :: IDUMMA 
      REAL(DOUBLE) :: ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, TMP1, TMP2, &
         TNORM, ULP, WKILL, WL, WLU, WU, WUL 
      LOGICAL :: NCNVRG, TOOFEW 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, INT, LOG, MAX, MIN, SQRT 
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
!  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
!  matrix T.  The user may ask for all eigenvalues, all eigenvalues
!  in the half-open interval (VL, VU], or the IL-th through IU-th
!  eigenvalues.
!
!  To avoid overflow, the matrix must be scaled so that its
!  largest element is no greater than overflow**(1/2) *
!  underflow**(1/4) in absolute value, and for greatest
!  accuracy, it should not be much smaller than that.
!
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!  Matrix", Report CS41, Computer Science Dept., Stanford
!  University, July 21, 1966.
!
!  Arguments
!  =========
!
!  RANGE   (input) CHARACTER
!          = 'A': ("All")   all eigenvalues will be found.
!          = 'V': ("Value") all eigenvalues in the half-open interval
!                           (VL, VU] will be found.
!          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
!                           entire matrix) will be found.
!
!  ORDER   (input) CHARACTER
!          = 'B': ("By Block") the eigenvalues will be grouped by
!                              split-off block (see IBLOCK, ISPLIT) and
!                              ordered from smallest to largest within
!                              the block.
!          = 'E': ("Entire matrix")
!                              the eigenvalues for the entire matrix
!                              will be ordered from smallest to
!                              largest.
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix T.  N >= 0.
!
!  VL      (input) DOUBLE PRECISION
!  VU      (input) DOUBLE PRECISION
!          If RANGE='V', the lower and upper bounds of the interval to
!          be searched for eigenvalues.  Eigenvalues less than or equal
!          to VL, or greater than VU, will not be returned.  VL < VU.
!          Not referenced if RANGE = 'A' or 'I'.
!
!  IL      (input) INTEGER
!  IU      (input) INTEGER
!          If RANGE='I', the indices (in ascending order) of the
!          smallest and largest eigenvalues to be returned.
!          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!          Not referenced if RANGE = 'A' or 'V'.
!
!  ABSTOL  (input) DOUBLE PRECISION
!          The absolute tolerance for the eigenvalues.  An eigenvalue
!          (or cluster) is considered to be located if it has been
!          determined to lie in an interval whose width is ABSTOL or
!          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
!          will be used, where |T| means the 1-norm of T.
!
!          Eigenvalues will be computed most accurately when ABSTOL is
!          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) off-diagonal elements of the tridiagonal matrix T.
!
!  M       (output) INTEGER
!          The actual number of eigenvalues found. 0 <= M <= N.
!          (See also the description of INFO=2,3.)
!
!  NSPLIT  (output) INTEGER
!          The number of diagonal blocks in the matrix T.
!          1 <= NSPLIT <= N.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On exit, the first M elements of W will contain the
!          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  IBLOCK  (output) INTEGER array, dimension (N)
!          At each row/column j where E(j) is zero or small, the
!          matrix T is considered to split into a block diagonal
!          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
!          block (from 1 to the number of blocks) the eigenvalue W(i)
!          belongs.  (DSTEBZ may use the remaining N-M elements as
!          workspace.)
!
!  ISPLIT  (output) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to ISPLIT(1),
!          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!          etc., and the NSPLIT-th consists of rows/columns
!          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!          (Only the first NSPLIT elements will actually be used, but
!          since the user cannot know a priori what value NSPLIT will
!          have, N words must be reserved for ISPLIT.)
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
!
!  IWORK   (workspace) INTEGER array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  some or all of the eigenvalues failed to converge or
!                were not computed:
!                =1 or 3: Bisection failed to converge for some
!                        eigenvalues; these eigenvalues are flagged by a
!                        negative block number.  The effect is that the
!                        eigenvalues may not be as accurate as the
!                        absolute and relative tolerances.  This is
!                        generally caused by unexpectedly inaccurate
!                        arithmetic.
!                =2 or 3: RANGE='I' only: Not all of the eigenvalues
!                        IL:IU were found.
!                        Effect: M < IU+1-IL
!                        Cause:  non-monotonic arithmetic, causing the
!                                Sturm sequence to be non-monotonic.
!                        Cure:   recalculate, using RANGE='A', and pick
!                                out eigenvalues IL:IU.  In some cases,
!                                increasing the PARAMETER "FUDGE" may
!                                make things work.
!                = 4:    RANGE='I', and the Gershgorin interval
!                        initially used was too small.  No eigenvalues
!                        were computed.
!                        Probable cause: your machine has sloppy
!                                        floating-point arithmetic.
!                        Cure: Increase the PARAMETER "FUDGE",
!                              recompile, and try again.
!
!  Internal Parameters
!  ===================
!
!  RELFAC  DOUBLE PRECISION, default = 2.0e0
!          The relative tolerance.  An interval (a,b] lies within
!          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
!          where "ulp" is the machine precision (distance from 1 to
!          the next larger floating point number.)
!
!  FUDGE   DOUBLE PRECISION, default = 2
!          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
!          a value of 1 should work, but on machines with sloppy
!          arithmetic, this needs to be larger.  The default for
!          publicly released versions should be large enough to handle
!          the worst machine around.  Note that this has no effect
!          on accuracy of the solution.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      INFO = 0 
!
!     Decode RANGE
!
      IF (LSAME(RANGE,'A')) THEN 
         IRANGE = 1 
      ELSE IF (LSAME(RANGE,'V')) THEN 
         IRANGE = 2 
      ELSE IF (LSAME(RANGE,'I')) THEN 
         IRANGE = 3 
      ELSE 
         IRANGE = 0 
      ENDIF 
!
!     Decode ORDER
!
      IF (LSAME(ORDER,'B')) THEN 
         IORDER = 2 
      ELSE IF (LSAME(ORDER,'E')) THEN 
         IORDER = 1 
      ELSE 
         IORDER = 0 
      ENDIF 
!
!     Check for Errors
!
      IF (IRANGE <= 0) THEN 
         INFO = -1 
      ELSE IF (IORDER <= 0) THEN 
         INFO = -2 
      ELSE IF (N < 0) THEN 
         INFO = -3 
      ELSE IF (IRANGE==2 .AND. VL>=VU) THEN 
         INFO = -5 
      ELSE IF (IRANGE==3 .AND. (IL<1 .OR. IL>MAX(1,N))) THEN 
         INFO = -6 
      ELSE IF (IRANGE==3 .AND. (IU<MIN(N,IL) .OR. IU>N)) THEN 
         INFO = -7 
      ENDIF 
!
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSTEBZ', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Initialize error flags
!
      INFO = 0 
      NCNVRG = .FALSE. 
      TOOFEW = .FALSE. 
!
!     Quick return if possible
!
      M = 0 
      IF (N == 0) RETURN  
!
!     Simplifications:
!
      IF (IRANGE==3 .AND. IL==1 .AND. IU==N) IRANGE = 1 
!
!     Get machine constants
!     NB is the minimum vector length for vector bisection, or 0
!     if only scalar is to be done.
!
      SAFEMN = DLAMCH('S') 
      ULP = DLAMCH('P') 
      RTOLI = ULP*RELFAC 
      NB = ILAENV(1,'DSTEBZ',' ',N,-1,-1,-1) 
      IF (NB <= 1) NB = 0 
!
!     Special Case when N=1
!
      IF (N == 1) THEN 
         NSPLIT = 1 
         ISPLIT(1) = 1 
         IF (IRANGE==2 .AND. (VL>=D(1) .OR. VU<D(1))) THEN 
            M = 0 
         ELSE 
            W(1) = D(1) 
            IBLOCK(1) = 1 
            M = 1 
         ENDIF 
         RETURN  
      ENDIF 
!
!     Compute Splitting Points
!
      NSPLIT = 1 
      WORK(N) = ZERO 
      PIVMIN = ONE 
!
      DO J = 2, N 
         TMP1 = E(J-1)**2 
         IF (ABS(D(J)*D(J-1))*ULP**2 + SAFEMN > TMP1) THEN 
            ISPLIT(NSPLIT) = J - 1 
            NSPLIT = NSPLIT + 1 
            WORK(J-1) = ZERO 
         ELSE 
            WORK(J-1) = TMP1 
            PIVMIN = MAX(PIVMIN,TMP1) 
         ENDIF 
      END DO 
      ISPLIT(NSPLIT) = N 
      PIVMIN = PIVMIN*SAFEMN 
!
!     Compute Interval and ATOLI
!
      IF (IRANGE == 3) THEN 
!
!        RANGE='I': Compute the interval containing eigenvalues
!                   IL through IU.
!
!        Compute Gershgorin interval for entire (split) matrix
!        and use it as the initial interval
!
         GU = D(1) 
         GL = D(1) 
         TMP1 = ZERO 
!
         DO J = 1, N - 1 
            TMP2 = SQRT(WORK(J)) 
            GU = MAX(GU,D(J)+TMP1+TMP2) 
            GL = MIN(GL,D(J)-TMP1-TMP2) 
            TMP1 = TMP2 
         END DO 
!
         GU = MAX(GU,D(N)+TMP1) 
         GL = MIN(GL,D(N)-TMP1) 
         TNORM = MAX(ABS(GL),ABS(GU)) 
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN 
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN 
!
!        Compute Iteration parameters
!
         ITMAX = INT((LOG(TNORM + PIVMIN) - LOG(PIVMIN))/LOG(TWO)) + 2 
         IF (ABSTOL <= ZERO) THEN 
            ATOLI = ULP*TNORM 
         ELSE 
            ATOLI = ABSTOL 
         ENDIF 
!
         WORK(N+1) = GL 
         WORK(N+2) = GL 
         WORK(N+3) = GU 
         WORK(N+4) = GU 
         WORK(N+5) = GL 
         WORK(N+6) = GU 
         IWORK(1) = -1 
         IWORK(2) = -1 
         IWORK(3) = N + 1 
         IWORK(4) = N + 1 
         IWORK(5) = IL - 1 
         IWORK(6) = IU 
!
         CALL DLAEBZ (3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, WORK, &
            IWORK(5), WORK(N+1), WORK(N+5), IOUT, IWORK, W, IBLOCK, IINFO) 
!
         IF (IWORK(6) == IU) THEN 
            WL = WORK(N+1) 
            WLU = WORK(N+3) 
            NWL = IWORK(1) 
            WU = WORK(N+4) 
            WUL = WORK(N+2) 
            NWU = IWORK(4) 
         ELSE 
            WL = WORK(N+2) 
            WLU = WORK(N+4) 
            NWL = IWORK(2) 
            WU = WORK(N+3) 
            WUL = WORK(N+1) 
            NWU = IWORK(3) 
         ENDIF 
!
         IF (NWL<0 .OR. NWL>=N .OR. NWU<1 .OR. NWU>N) THEN 
            INFO = 4 
            RETURN  
         ENDIF 
      ELSE 
!
!        RANGE='A' or 'V' -- Set ATOLI
!
         TNORM = MAX(ABS(D(1))+ABS(E(1)),ABS(D(N))+ABS(E(N-1))) 
!
         DO J = 2, N - 1 
            TNORM = MAX(TNORM,ABS(D(J))+ABS(E(J-1))+ABS(E(J))) 
         END DO 
!
         IF (ABSTOL <= ZERO) THEN 
            ATOLI = ULP*TNORM 
         ELSE 
            ATOLI = ABSTOL 
         ENDIF 
!
         IF (IRANGE == 2) THEN 
            WL = VL 
            WU = VU 
         ENDIF 
      ENDIF 
!
!     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
!     NWL accumulates the number of eigenvalues .le. WL,
!     NWU accumulates the number of eigenvalues .le. WU
!
      M = 0 
      IEND = 0 
      INFO = 0 
      NWL = 0 
      NWU = 0 
!
      DO JB = 1, NSPLIT 
         IOFF = IEND 
         IBEGIN = IOFF + 1 
         IEND = ISPLIT(JB) 
         IN = IEND - IOFF 
!
         IF (IN == 1) THEN 
!
!           Special Case -- IN=1
!
            IF (IRANGE==1 .OR. WL>=D(IBEGIN)-PIVMIN) NWL = NWL + 1 
            IF (IRANGE==1 .OR. WU>=D(IBEGIN)-PIVMIN) NWU = NWU + 1 
            IF (IRANGE==1 .OR. WL<D(IBEGIN)-PIVMIN .AND. WU>=D(IBEGIN)-PIVMIN) &
               THEN 
               M = M + 1 
               W(M) = D(IBEGIN) 
               IBLOCK(M) = JB 
            ENDIF 
         ELSE 
!
!           General Case -- IN > 1
!
!           Compute Gershgorin Interval
!           and use it as the initial interval
!
            GU = D(IBEGIN) 
            GL = D(IBEGIN) 
            TMP1 = ZERO 
!
            DO J = IBEGIN, IEND - 1 
               TMP2 = ABS(E(J)) 
               GU = MAX(GU,D(J)+TMP1+TMP2) 
               GL = MIN(GL,D(J)-TMP1-TMP2) 
               TMP1 = TMP2 
            END DO 
!
            GU = MAX(GU,D(IEND)+TMP1) 
            GL = MIN(GL,D(IEND)-TMP1) 
            BNORM = MAX(ABS(GL),ABS(GU)) 
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN 
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN 
!
!           Compute ATOLI for the current submatrix
!
            IF (ABSTOL <= ZERO) THEN 
               ATOLI = ULP*MAX(ABS(GL),ABS(GU)) 
            ELSE 
               ATOLI = ABSTOL 
            ENDIF 
!
            IF (IRANGE > 1) THEN 
               IF (GU < WL) THEN 
                  NWL = NWL + IN 
                  NWU = NWU + IN 
                  CYCLE  
               ENDIF 
               GL = MAX(GL,WL) 
               GU = MIN(GU,WU) 
               IF (GL >= GU) CYCLE  
            ENDIF 
!
!           Set Up Initial Interval
!
            WORK(N+1) = GL 
            WORK(N+IN+1) = GU 
            CALL DLAEBZ (1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D(IBEGIN), &
               E(IBEGIN), WORK(IBEGIN), IDUMMA, WORK(N+1), WORK(N+2*IN+1), IM, &
               IWORK, W(M+1), IBLOCK(M+1), IINFO) 
!
            NWL = NWL + IWORK(1) 
            NWU = NWU + IWORK(IN+1) 
            IWOFF = M - IWORK(1) 
!
!           Compute Eigenvalues
!
            ITMAX = INT((LOG(GU - GL + PIVMIN) - LOG(PIVMIN))/LOG(TWO)) + 2 
            CALL DLAEBZ (2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D(&
               IBEGIN), E(IBEGIN), WORK(IBEGIN), IDUMMA, WORK(N+1), WORK(N+2*IN&
               +1), IOUT, IWORK, W(M+1), IBLOCK(M+1), IINFO) 
!
!           Copy Eigenvalues Into W and IBLOCK
!           Use -JB for block number for unconverged eigenvalues.
!
            DO J = 1, IOUT 
               TMP1 = HALF*(WORK(J+N)+WORK(J+IN+N)) 
!
!              Flag non-convergence.
!
               IF (J > IOUT - IINFO) THEN 
                  NCNVRG = .TRUE. 
                  IB = -JB 
               ELSE 
                  IB = JB 
               ENDIF 
               W(IWORK(J)+1+IWOFF:IWORK(J+IN)+IWOFF) = TMP1 
               IBLOCK(IWORK(J)+1+IWOFF:IWORK(J+IN)+IWOFF) = IB 
            END DO 
!
            M = M + IM 
         ENDIF 
      END DO 
!
!     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      IF (IRANGE == 3) THEN 
         IM = 0 
         IDISCL = IL - 1 - NWL 
         IDISCU = NWU - IU 
!
         IF (IDISCL>0 .OR. IDISCU>0) THEN 
            DO JE = 1, M 
               IF (W(JE)<=WLU .AND. IDISCL>0) THEN 
                  IDISCL = IDISCL - 1 
               ELSE IF (W(JE)>=WUL .AND. IDISCU>0) THEN 
                  IDISCU = IDISCU - 1 
               ELSE 
                  IM = IM + 1 
                  W(IM) = W(JE) 
                  IBLOCK(IM) = IBLOCK(JE) 
               ENDIF 
            END DO 
            M = IM 
         ENDIF 
         IF (IDISCL>0 .OR. IDISCU>0) THEN 
!
!           Code to deal with effects of bad arithmetic:
!           Some low eigenvalues to be discarded are not in (WL,WLU],
!           or high eigenvalues to be discarded are not in (WUL,WU]
!           so just kill off the smallest IDISCL/largest IDISCU
!           eigenvalues, by simply finding the smallest/largest
!           eigenvalue(s).
!
!           (If N(w) is monotone non-decreasing, this should never
!               happen.)
!
            IF (IDISCL > 0) THEN 
               WKILL = WU 
               DO JDISC = 1, IDISCL 
                  IW = 0 
                  DO JE = 1, M 
                     IF (.NOT.(IBLOCK(JE)/=0 .AND. (W(JE)<WKILL .OR. IW==0))) &
                        CYCLE  
                     IW = JE 
                     WKILL = W(JE) 
                  END DO 
                  IBLOCK(IW) = 0 
               END DO 
            ENDIF 
            IF (IDISCU > 0) THEN 
!
               WKILL = WL 
               DO JDISC = 1, IDISCU 
                  IW = 0 
                  DO JE = 1, M 
                     IF (.NOT.(IBLOCK(JE)/=0 .AND. (W(JE)>WKILL .OR. IW==0))) &
                        CYCLE  
                     IW = JE 
                     WKILL = W(JE) 
                  END DO 
                  IBLOCK(IW) = 0 
               END DO 
            ENDIF 
            IM = 0 
            DO JE = 1, M 
               IF (IBLOCK(JE) == 0) CYCLE  
               IM = IM + 1 
               W(IM) = W(JE) 
               IBLOCK(IM) = IBLOCK(JE) 
            END DO 
            M = IM 
         ENDIF 
         IF (IDISCL<0 .OR. IDISCU<0) TOOFEW = .TRUE. 
      ENDIF 
!
!     If ORDER='B', do nothing -- the eigenvalues are already sorted
!        by block.
!     If ORDER='E', sort the eigenvalues from smallest to largest
!
      IF (IORDER==1 .AND. NSPLIT>1) THEN 
         DO JE = 1, M - 1 
            IE = 0 
            TMP1 = W(JE) 
            DO J = JE + 1, M 
               IF (W(J) >= TMP1) CYCLE  
               IE = J 
               TMP1 = W(J) 
            END DO 
!
            IF (IE == 0) CYCLE  
            ITMP1 = IBLOCK(IE) 
            W(IE) = W(JE) 
            IBLOCK(IE) = IBLOCK(JE) 
            W(JE) = TMP1 
            IBLOCK(JE) = ITMP1 
         END DO 
      ENDIF 
!
      INFO = 0 
      IF (NCNVRG) INFO = INFO + 1 
      IF (TOOFEW) INFO = INFO + 2 
      RETURN  
!
!     End of DSTEBZ
!
      END SUBROUTINE DSTEBZ 
