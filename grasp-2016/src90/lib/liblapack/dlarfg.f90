      SUBROUTINE DLARFG(N, ALPHA, X, INCX, TAU) 
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
      !USE dlapy2_I 
      !USE dnrm2_I 
      !USE dscal_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER  :: INCX 
      REAL(DOUBLE)  :: ALPHA 
      REAL(DOUBLE) , INTENT(OUT) :: TAU 
      REAL(DOUBLE)  :: X(*) 
      REAL(KIND(0.0D0))::dlamch, dlapy2
      real(kind(0.0d0))::dnrm2
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, KNT 
      REAL(DOUBLE) :: BETA, RSAFMN, SAFMIN, XNORM 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, SIGN 
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
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
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
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      IF (N <= 1) THEN 
         TAU = ZERO 
         RETURN  
      ENDIF 
!
      XNORM = DNRM2(N - 1,X,INCX) 
!
      IF (XNORM == ZERO) THEN 
!
!        H  =  I
!
         TAU = ZERO 
      ELSE 
!
!        general case
!
         BETA = -SIGN(DLAPY2(ALPHA,XNORM),ALPHA) 
         SAFMIN = DLAMCH('S')/DLAMCH('E') 
         IF (ABS(BETA) < SAFMIN) THEN 
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE/SAFMIN 
            KNT = 0 
            KNT = KNT + 1 
            CALL DSCAL (N - 1, RSAFMN, X, INCX) 
            BETA = BETA*RSAFMN 
            ALPHA = ALPHA*RSAFMN 
            DO WHILE(ABS(BETA) < SAFMIN) 
               KNT = KNT + 1 
               CALL DSCAL (N - 1, RSAFMN, X, INCX) 
               BETA = BETA*RSAFMN 
               ALPHA = ALPHA*RSAFMN 
            END DO 
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DNRM2(N - 1,X,INCX) 
            BETA = -SIGN(DLAPY2(ALPHA,XNORM),ALPHA) 
            TAU = (BETA - ALPHA)/BETA 
            CALL DSCAL (N - 1, ONE/(ALPHA - BETA), X, INCX) 
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA 
            DO J = 1, KNT 
               ALPHA = ALPHA*SAFMIN 
            END DO 
         ELSE 
            TAU = (BETA - ALPHA)/BETA 
            CALL DSCAL (N - 1, ONE/(ALPHA - BETA), X, INCX) 
            ALPHA = BETA 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DLARFG
!
      END SUBROUTINE DLARFG 
