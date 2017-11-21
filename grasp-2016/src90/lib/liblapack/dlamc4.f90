!
!***********************************************************************
!
      SUBROUTINE DLAMC4(EMIN, START, BASE) 
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
      INTEGER , INTENT(OUT) :: EMIN 
      INTEGER , INTENT(IN) :: BASE 
      REAL(DOUBLE) , INTENT(IN) :: START 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO 
      real(double)::dlamc3
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
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      A = START 
      ONE = 1 
      RBASE = ONE/BASE 
      ZERO = 0 
      EMIN = 1 
      B1 = DLAMC3(A*RBASE,ZERO) 
      C1 = A 
      C2 = A 
      D1 = A 
      D2 = A 
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE 
      IF (C1==A .AND. C2==A .AND. D1==A .AND. D2==A) THEN 
         EMIN = EMIN - 1 
         A = B1 
         B1 = DLAMC3(A/BASE,ZERO) 
         C1 = DLAMC3(B1*BASE,ZERO) 
         D1 = ZERO 
         IF (BASE > 0) THEN 
            D1 = (BASE - 1)*B1 + D1 + B1 
         ENDIF 
         B2 = DLAMC3(A*RBASE,ZERO) 
         C2 = DLAMC3(B2/RBASE,ZERO) 
         D2 = ZERO 
         IF (BASE > 0) THEN 
            D2 = (BASE - 1)*B2 + D2 + B2 
         ENDIF 
         GO TO 10 
      ENDIF 
!+    END WHILE
!
      RETURN  
!
!     End of DLAMC4
!
      END SUBROUTINE DLAMC4 
