      SUBROUTINE DSTEIN(N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK, IWORK, &
         IFAIL, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:36   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE idamax_I 
      !USE dasum_I 
      !USE ddot_I 
      !!USE dlamch_I 
      !USE dnrm2_I 
      !USE daxpy_I 
      !USE dcopy_I 
      !USE dlagtf_I 
      !USE dlagts_I 
      !USE dlarnv_I 
      !USE dscal_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: LDZ 
      INTEGER , INTENT(OUT) :: INFO 
      INTEGER , INTENT(IN) :: IBLOCK(*) 
      INTEGER , INTENT(IN) :: ISPLIT(*) 
      INTEGER  :: IWORK(*) 
      INTEGER , INTENT(OUT) :: IFAIL(*) 
      REAL(DOUBLE)  :: D(*) 
      REAL(DOUBLE)  :: E(*) 
      REAL(DOUBLE) , INTENT(IN) :: W(*) 
      REAL(DOUBLE)  :: Z(LDZ,*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: TEN = 1.0D+1 
      REAL(DOUBLE), PARAMETER :: ODM3 = 1.0D-3 
      REAL(DOUBLE), PARAMETER :: ODM1 = 1.0D-1 
      INTEGER, PARAMETER :: MAXITS = 5 
      INTEGER, PARAMETER :: EXTRA = 2 
      real(kind(0.0d0))::dnrm2
      REAL(KIND(0.0D0)) :: dlamch
      real(kind(0.0d0))::ddot
       real(kind(0.0d0))::dasum
      integer::idamax
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1, INDRV2, INDRV3, &
         INDRV4, INDRV5, ITS, J, J1, JBLK, JMAX, NBLK, NRMCHK 
      INTEGER, DIMENSION(4) :: ISEED 
      REAL(DOUBLE) :: DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL, SCL, SEP, &
         TOL, XJ, XJM, ZTR 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, SQRT 
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
!  DSTEIN computes the eigenvectors of a real symmetric tridiagonal
!  matrix T corresponding to specified eigenvalues, using inverse
!  iteration.
!
!  The maximum number of iterations allowed for each eigenvector is
!  specified by an internal parameter MAXITS (currently set to 5).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix T.
!
!  E       (input) DOUBLE PRECISION array, dimension (N)
!          The (n-1) subdiagonal elements of the tridiagonal matrix
!          T, in elements 1 to N-1.  E(N) need not be set.
!
!  M       (input) INTEGER
!          The number of eigenvectors to be found.  0 <= M <= N.
!
!  W       (input) DOUBLE PRECISION array, dimension (N)
!          The first M elements of W contain the eigenvalues for
!          which eigenvectors are to be computed.  The eigenvalues
!          should be grouped by split-off block and ordered from
!          smallest to largest within the block.  ( The output array
!          W from DSTEBZ with ORDER = 'B' is expected here. )
!
!  IBLOCK  (input) INTEGER array, dimension (N)
!          The submatrix indices associated with the corresponding
!          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
!          the first submatrix from the top, =2 if W(i) belongs to
!          the second submatrix, etc.  ( The output array IBLOCK
!          from DSTEBZ is expected here. )
!
!  ISPLIT  (input) INTEGER array, dimension (N)
!          The splitting points, at which T breaks up into submatrices.
!          The first submatrix consists of rows/columns 1 to
!          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
!          through ISPLIT( 2 ), etc.
!          ( The output array ISPLIT from DSTEBZ is expected here. )
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)
!          The computed eigenvectors.  The eigenvector associated
!          with the eigenvalue W(i) is stored in the i-th column of
!          Z.  Any vector which fails to converge is set to its current
!          iterate after MAXITS iterations.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  IFAIL   (output) INTEGER array, dimension (M)
!          On normal exit, all elements of IFAIL are zero.
!          If one or more eigenvectors fail to converge after
!          MAXITS iterations, then their indices are stored in
!          array IFAIL.
!
!  INFO    (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, then i eigenvectors failed to converge
!               in MAXITS iterations.  Their indices are stored in
!               array IFAIL.
!
!  Internal Parameters
!  ===================
!
!  MAXITS  INTEGER, default = 5
!          The maximum number of iterations performed.
!
!  EXTRA   INTEGER, default = 2
!          The number of iterations performed after norm growth
!          criterion is satisfied, should be at least 1.
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
!     Test the input parameters.
!
      INFO = 0 
      IFAIL(:M) = 0 
!
      IF (N < 0) THEN 
         INFO = -1 
      ELSE IF (M<0 .OR. M>N) THEN 
         INFO = -4 
      ELSE IF (LDZ < MAX(1,N)) THEN 
         INFO = -9 
      ELSE 
         DO J = 2, M 
            IF (IBLOCK(J) < IBLOCK(J-1)) THEN 
               INFO = -6 
               EXIT  
            ENDIF 
            IF (IBLOCK(J)/=IBLOCK(J-1) .OR. W(J)>=W(J-1)) CYCLE  
            INFO = -5 
            EXIT  
         END DO 
      ENDIF 
!
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSTEIN', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N==0 .OR. M==0) THEN 
         RETURN  
      ELSE IF (N == 1) THEN 
         Z(1,1) = ONE 
         RETURN  
      ENDIF 
!
!     Get machine constants.
!
      EPS = DLAMCH('Precision') 
!
!     Initialize seed for random number generator DLARNV.
!
      ISEED = 1 
!
!     Initialize pointers.
!
      INDRV1 = 0 
      INDRV2 = INDRV1 + N 
      INDRV3 = INDRV2 + N 
      INDRV4 = INDRV3 + N 
      INDRV5 = INDRV4 + N 
!
!     Compute eigenvectors of matrix blocks.
!
      J1 = 1 
      L160: DO NBLK = 1, IBLOCK(M) 
!
!        Find starting and ending indices of block nblk.
!
         IF (NBLK == 1) THEN 
            B1 = 1 
         ELSE 
            B1 = ISPLIT(NBLK-1) + 1 
         ENDIF 
         BN = ISPLIT(NBLK) 
         BLKSIZ = BN - B1 + 1 
         IF (BLKSIZ /= 1) THEN 
            GPIND = B1 
!
!        Compute reorthogonalization criterion and stopping criterion.
!
            ONENRM = ABS(D(B1)) + ABS(E(B1)) 
            ONENRM = MAX(ONENRM,ABS(D(BN))+ABS(E(BN-1))) 
            DO I = B1 + 1, BN - 1 
               ONENRM = MAX(ONENRM,ABS(D(I))+ABS(E(I-1))+ABS(E(I))) 
            END DO 
            ORTOL = ODM3*ONENRM 
!
            DTPCRT = SQRT(ODM1/BLKSIZ) 
         ENDIF 
!
!        Loop through eigenvalues of block nblk.
!
         JBLK = 0 
         DO J = J1, M 
            IF (IBLOCK(J) /= NBLK) THEN 
               J1 = J 
               CYCLE  L160 
            ENDIF 
            JBLK = JBLK + 1 
            XJ = W(J) 
!
!           Skip all the work if the block size is one.
!
            IF (BLKSIZ == 1) THEN 
               WORK(INDRV1+1) = ONE 
               GO TO 120 
            ENDIF 
!
!           If eigenvalues j and j-1 are too close, add a relatively
!           small perturbation.
!
            IF (JBLK > 1) THEN 
               EPS1 = ABS(EPS*XJ) 
               PERTOL = TEN*EPS1 
               SEP = XJ - XJM 
               IF (SEP < PERTOL) XJ = XJM + PERTOL 
            ENDIF 
!
            ITS = 0 
            NRMCHK = 0 
!
!           Get random starting vector.
!
            CALL DLARNV (2, ISEED, BLKSIZ, WORK(INDRV1+1)) 
!
!           Copy the matrix T so it won't be destroyed in factorization.
!
!jbdbx           print *, ' BLKSIZ = ', BLKSIZ
            CALL DCOPY (BLKSIZ, D(B1), 1, WORK(INDRV4+1), 1) 
            CALL DCOPY (BLKSIZ - 1, E(B1), 1, WORK(INDRV2+2), 1) 
            CALL DCOPY (BLKSIZ - 1, E(B1), 1, WORK(INDRV3+1), 1) 
!
!           Compute LU factors with partial pivoting  ( PT = LU )
!
            TOL = ZERO 
            CALL DLAGTF (BLKSIZ, WORK(INDRV4+1), XJ, WORK(INDRV2+2), WORK(&
               INDRV3+1), TOL, WORK(INDRV5+1), IWORK, IINFO) 
!
!           Update iteration count.
!
   70       CONTINUE 
            ITS = ITS + 1 
            IF (ITS <= MAXITS) THEN 
!
!           Normalize and scale the righthand side vector Pb.
!
               SCL = BLKSIZ*ONENRM*MAX(EPS,ABS(WORK(INDRV4+BLKSIZ)))/DASUM(&
                  BLKSIZ,WORK(INDRV1+1),1) 
               CALL DSCAL (BLKSIZ, SCL, WORK(INDRV1+1), 1) 
!
!           Solve the system LU = Pb.
!
               CALL DLAGTS ((-1), BLKSIZ, WORK(INDRV4+1), WORK(INDRV2+2), WORK(&
                  INDRV3+1), WORK(INDRV5+1), IWORK, WORK(INDRV1+1), TOL, IINFO) 
!
!           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!           close enough.
!
               IF (JBLK /= 1) THEN 
                  IF (ABS(XJ - XJM) > ORTOL) GPIND = J 
                  IF (GPIND /= J) THEN 
                     DO I = GPIND, J - 1 
                        ZTR = -DDOT(BLKSIZ,WORK(INDRV1+1),1,Z(B1,I),1) 
                        CALL DAXPY (BLKSIZ, ZTR, Z(B1,I), 1, WORK(INDRV1+1), 1) 
                     END DO 
                  ENDIF 
               ENDIF 
!
!           Check the infinity norm of the iterate.
!
               JMAX = IDAMAX(BLKSIZ,WORK(INDRV1+1),1) 
               NRM = ABS(WORK(INDRV1+JMAX)) 
!
!           Continue for additional iterations after norm reaches
!           stopping criterion.
!
               IF (NRM < DTPCRT) GO TO 70 
               NRMCHK = NRMCHK + 1 
               IF (NRMCHK < EXTRA + 1) GO TO 70 
!
            ELSE 
!
!           If stopping criterion was not satisfied, update info and
!           store eigenvector number in array ifail.
!
               INFO = INFO + 1 
               IFAIL(INFO) = J 
            ENDIF 
!
!           Accept iterate as jth eigenvector.
!
            SCL = ONE/DNRM2(BLKSIZ,WORK(INDRV1+1),1) 
            JMAX = IDAMAX(BLKSIZ,WORK(INDRV1+1),1) 
            IF (WORK(INDRV1+JMAX) < ZERO) SCL = -SCL 
            CALL DSCAL (BLKSIZ, SCL, WORK(INDRV1+1), 1) 
  120       CONTINUE 
            Z(:N,J) = ZERO 
            Z(B1:BLKSIZ-1+B1,J) = WORK(INDRV1+1:BLKSIZ+INDRV1) 
!
!           Save the shift to check eigenvalue spacing at next
!           iteration.
!
            XJM = XJ 
!
         END DO 
      END DO L160 
!
      RETURN  
!
!     End of DSTEIN
!
      END SUBROUTINE DSTEIN 
