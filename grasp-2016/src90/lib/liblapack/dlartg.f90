      SUBROUTINE DLARTG (F, G, CS, SN, R) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:32   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlamch_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: F 
      REAL(DOUBLE) , INTENT(IN) :: G 
      REAL(DOUBLE) , INTENT(OUT) :: CS 
      REAL(DOUBLE) , INTENT(OUT) :: SN 
      REAL(DOUBLE) , INTENT(OUT) :: R 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: COUNTerotation, I 
      REAL(DOUBLE) :: EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE 
      LOGICAL :: FIRST 
      REAL(KIND(0.0D0)) :: dlamch
      SAVE FIRST, SAFMIN, SAFMN2, SAFMX2 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, INT, LOG, MAX, SQRT 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.
!
!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.
!
!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Save statement ..
!     ..
!     .. Data statements ..
      DATA FIRST/ .TRUE./  
!     ..
!     .. Executable Statements ..
!
      IF (FIRST) THEN 
         FIRST = .FALSE. 
         SAFMIN = DLAMCH('S') 
         EPS = DLAMCH('E') 
         SAFMN2 = DLAMCH('B')**INT(LOG(SAFMIN/EPS)/LOG(DLAMCH('B'))/TWO) 
         SAFMX2 = ONE/SAFMN2 
      ENDIF 
      IF (G == ZERO) THEN 
         CS = ONE 
         SN = ZERO 
         R = F 
      ELSE IF (F == ZERO) THEN 
         CS = ZERO 
         SN = ONE 
         R = G 
      ELSE 
         F1 = F 
         G1 = G 
         SCALE = MAX(ABS(F1),ABS(G1)) 
         IF (SCALE >= SAFMX2) THEN 
            COUNTerotation = 0 
            COUNTerotation = COUNTerotation + 1 
            F1 = F1*SAFMN2 
            G1 = G1*SAFMN2 
            SCALE = MAX(ABS(F1),ABS(G1)) 
            DO WHILE(SCALE >= SAFMX2) 
               COUNTerotation = COUNTerotation + 1 
               F1 = F1*SAFMN2 
               G1 = G1*SAFMN2 
               SCALE = MAX(ABS(F1),ABS(G1)) 
            END DO 
            R = SQRT(F1**2 + G1**2) 
            CS = F1/R 
            SN = G1/R 
            DO I = 1, COUNTerotation 
               R = R*SAFMX2 
            END DO 
         ELSE IF (SCALE <= SAFMN2) THEN 
            COUNTerotation = 0 
            COUNTerotation = COUNTerotation + 1 
            F1 = F1*SAFMX2 
            G1 = G1*SAFMX2 
            SCALE = MAX(ABS(F1),ABS(G1)) 
            DO WHILE(SCALE <= SAFMN2) 
               COUNTerotation = COUNTerotation + 1 
               F1 = F1*SAFMX2 
               G1 = G1*SAFMX2 
               SCALE = MAX(ABS(F1),ABS(G1)) 
            END DO 
            R = SQRT(F1**2 + G1**2) 
            CS = F1/R 
            SN = G1/R 
            DO I = 1, COUNTerotation 
               R = R*SAFMN2 
            END DO 
         ELSE 
            R = SQRT(F1**2 + G1**2) 
            CS = F1/R 
            SN = G1/R 
         ENDIF 
         IF (ABS(F)>ABS(G) .AND. CS<ZERO) THEN 
            CS = -CS 
            SN = -SN 
            R = -R 
         ENDIF 
      ENDIF 
      RETURN  
!
!     End of DLARTG
!
      END SUBROUTINE DLARTG 
