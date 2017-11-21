      SUBROUTINE DLARNV(IDIST, ISEED, N, X) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:32   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlaruv_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDIST 
      INTEGER , INTENT(IN) :: N 
      INTEGER  :: ISEED(4) 
      REAL(DOUBLE) , INTENT(OUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: TWO = 2.0D+0 
      INTEGER, PARAMETER :: LV = 128 
      REAL(DOUBLE), PARAMETER :: TWOPI = 6.2831853071795864769252867663D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IL, IL2, IV 
      REAL(DOUBLE), DIMENSION(LV) :: U 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC COS, LOG, MIN, SQRT 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
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
!  DLARNV returns a vector of n random real numbers from a uniform or
!  normal distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) DOUBLE PRECISION array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine DLARUV to generate random
!  real numbers from a uniform (0,1) distribution, in batches of up to
!  128 using vectorisable code. The Box-Muller method is used to
!  transform numbers from a uniform to a normal distribution.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      DO IV = 1, N, LV/2 
         IL = MIN(LV/2,N - IV + 1) 
         IF (IDIST == 3) THEN 
            IL2 = 2*IL 
         ELSE 
            IL2 = IL 
         ENDIF 
!
!        Call DLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)
!
         CALL DLARUV (ISEED, IL2, U) 
!
         SELECT CASE (IDIST)  
         CASE (1)  
!
!           Copy generated numbers
!
            X(IV:IL-1+IV) = U(:IL) 
         CASE (2)  
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            X(IV:IL-1+IV) = TWO*U(:IL) - ONE 
         CASE (3)  
!
!           Convert generated numbers to normal (0,1) distribution
!
            DO I = 1, IL 
               X(IV+I-1) = SQRT((-TWO*LOG(U(2*I-1))))*COS(TWOPI*U(2*I)) 
            END DO 
         END SELECT 
      END DO 
      RETURN  
!
!     End of DLARNV
!
      END SUBROUTINE DLARNV 
