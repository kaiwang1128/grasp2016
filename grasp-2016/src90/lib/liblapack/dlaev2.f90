      SUBROUTINE DLAEV2(A, B, C, RT1, RT2, CS1, SN1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:29   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: A 
      REAL(DOUBLE) , INTENT(IN) :: B 
      REAL(DOUBLE) , INTENT(IN) :: C 
      REAL(DOUBLE) , INTENT(OUT) :: RT1 
      REAL(DOUBLE) , INTENT(OUT) :: RT2 
      REAL(DOUBLE) , INTENT(OUT) :: CS1 
      REAL(DOUBLE) , INTENT(OUT) :: SN1 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: HALF = 0.5D0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: SGN1, SGN2 
      REAL(DOUBLE) :: AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, SQRT 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) DOUBLE PRECISION
!  SN1     (output) DOUBLE PRECISION
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C 
      DF = A - C 
      ADF = ABS(DF) 
      TB = B + B 
      AB = ABS(TB) 
      IF (ABS(A) > ABS(C)) THEN 
         ACMX = A 
         ACMN = C 
      ELSE 
         ACMX = C 
         ACMN = A 
      ENDIF 
      IF (ADF > AB) THEN 
         RT = ADF*SQRT(ONE + (AB/ADF)**2) 
      ELSE IF (ADF < AB) THEN 
         RT = AB*SQRT(ONE + (ADF/AB)**2) 
      ELSE 
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT(TWO) 
      ENDIF 
      IF (SM < ZERO) THEN 
         RT1 = HALF*(SM - RT) 
         SGN1 = -1 
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B 
      ELSE IF (SM > ZERO) THEN 
         RT1 = HALF*(SM + RT) 
         SGN1 = 1 
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B 
      ELSE 
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT 
         RT2 = -HALF*RT 
         SGN1 = 1 
      ENDIF 
!
!     Compute the eigenvector
!
      IF (DF >= ZERO) THEN 
         CS = DF + RT 
         SGN2 = 1 
      ELSE 
         CS = DF - RT 
         SGN2 = -1 
      ENDIF 
      ACS = ABS(CS) 
      IF (ACS > AB) THEN 
         CT = -TB/CS 
         SN1 = ONE/SQRT(ONE + CT*CT) 
         CS1 = CT*SN1 
      ELSE 
         IF (AB == ZERO) THEN 
            CS1 = ONE 
            SN1 = ZERO 
         ELSE 
            TN = -CS/TB 
            CS1 = ONE/SQRT(ONE + TN*TN) 
            SN1 = TN*CS1 
         ENDIF 
      ENDIF 
      IF (SGN1 == SGN2) THEN 
         TN = CS1 
         CS1 = -SN1 
         SN1 = TN 
      ENDIF 
      RETURN  
!
!     End of DLAEV2
!
      END SUBROUTINE DLAEV2 
