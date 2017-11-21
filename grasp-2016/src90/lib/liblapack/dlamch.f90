      REAL(KIND(0.0D0)) FUNCTION DLAMCH (CMACH) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:31   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USE dlamc2_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: CMACH 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: BETA, IMAX, IMIN, IT 
      REAL(DOUBLE) :: BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, RND, &
         SFMIN, SMALL, T 
      LOGICAL :: FIRST, LRND 
      logical::lsame
      SAVE FIRST, BASE, EMAX, EMIN, EPS, PREC, RMAX, RMIN, RND, SFMIN, T 
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
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
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
         CALL DLAMC2 (BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX) 
         BASE = BETA 
         T = IT 
         IF (LRND) THEN 
            RND = ONE 
            EPS = BASE**(1 - IT)/2 
         ELSE 
            RND = ZERO 
            EPS = BASE**(1 - IT) 
         ENDIF 
         PREC = EPS*BASE 
         EMIN = IMIN 
         EMAX = IMAX 
         SFMIN = RMIN 
         SMALL = ONE/RMAX 
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
         IF (SMALL >= SFMIN) SFMIN = SMALL*(ONE + EPS) 
      ENDIF 
!
      IF (LSAME(CMACH,'E')) THEN 
         RMACH = EPS 
      ELSE IF (LSAME(CMACH,'S')) THEN 
         RMACH = SFMIN 
      ELSE IF (LSAME(CMACH,'B')) THEN 
         RMACH = BASE 
      ELSE IF (LSAME(CMACH,'P')) THEN 
         RMACH = PREC 
      ELSE IF (LSAME(CMACH,'N')) THEN 
         RMACH = T 
      ELSE IF (LSAME(CMACH,'R')) THEN 
         RMACH = RND 
      ELSE IF (LSAME(CMACH,'M')) THEN 
         RMACH = EMIN 
      ELSE IF (LSAME(CMACH,'U')) THEN 
         RMACH = RMIN 
      ELSE IF (LSAME(CMACH,'L')) THEN 
         RMACH = EMAX 
      ELSE IF (LSAME(CMACH,'O')) THEN 
         RMACH = RMAX 
      ENDIF 
!
      DLAMCH = RMACH 
      RETURN  
!
!     End of DLAMCH
!
      END FUNCTION DLAMCH 
