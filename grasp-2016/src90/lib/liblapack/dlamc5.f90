!
!***********************************************************************
!
      SUBROUTINE DLAMC5(BETA, P, EMIN, IEEE, EMAX, RMAX) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:31   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlamc3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: BETA 
      INTEGER , INTENT(IN) :: P 
      INTEGER , INTENT(IN) :: EMIN 
      INTEGER , INTENT(OUT) :: EMAX 
      REAL(DOUBLE) , INTENT(OUT) :: RMAX 
      LOGICAL , INTENT(IN) :: IEEE 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D0 
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP 
      REAL(DOUBLE) :: OLDY, RECBAS, Y, Z 
      real(double) :: dlamc3
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MOD 
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
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1 
      EXBITS = 1 
   10 CONTINUE 
      TRY = LEXP*2 
      IF (TRY <= (-EMIN)) THEN 
         LEXP = TRY 
         EXBITS = EXBITS + 1 
         GO TO 10 
      ENDIF 
      IF (LEXP == (-EMIN)) THEN 
         UEXP = LEXP 
      ELSE 
         UEXP = TRY 
         EXBITS = EXBITS + 1 
      ENDIF 
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF (UEXP + EMIN > (-LEXP) - EMIN) THEN 
         EXPSUM = 2*LEXP 
      ELSE 
         EXPSUM = 2*UEXP 
      ENDIF 
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1 
      NBITS = 1 + EXBITS + P 
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF (MOD(NBITS,2)==1 .AND. BETA==2) THEN 
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1 
      ENDIF 
!
      IF (IEEE) THEN 
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1 
      ENDIF 
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE/BETA 
      Z = BETA - ONE 
      Y = ZERO 
      DO I = 1, P 
         Z = Z*RECBAS 
         IF (Y < ONE) OLDY = Y 
         Y = DLAMC3(Y,Z) 
      END DO 
      IF (Y >= ONE) Y = OLDY 
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO I = 1, EMAX 
         Y = DLAMC3(Y*BETA,ZERO) 
      END DO 
!
      RMAX = Y 
      RETURN  
!
!     End of DLAMC5
!
      END SUBROUTINE DLAMC5 
