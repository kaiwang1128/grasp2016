!
!***********************************************************************
!
      SUBROUTINE DLAMC1(BETA, T, RND, IEEE1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:30   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlamc3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: BETA 
      INTEGER , INTENT(OUT) :: T 
      LOGICAL , INTENT(OUT) :: RND 
      LOGICAL , INTENT(OUT) :: IEEE1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LBETA, LT 
      REAL(DOUBLE) :: A, B, C, F, ONE, QTR, SAVEC, T1, T2 
      LOGICAL :: FIRST, LIEEE1, LRND 
      real(double) :: dlamc3
      SAVE FIRST, LIEEE1, LRND, LBETA, LT 
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
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
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
         ONE = 1 
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1 
         C = 1 
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE 
         IF (C == ONE) THEN 
            A = 2*A 
            C = DLAMC3(A,ONE) 
            C = DLAMC3(C,(-A)) 
            GO TO 10 
         ENDIF 
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1 
         C = DLAMC3(A,B) 
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE 
         IF (C == A) THEN 
            B = 2*B 
            C = DLAMC3(A,B) 
            GO TO 20 
         ENDIF 
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE/4 
         SAVEC = C 
         C = DLAMC3(C,(-A)) 
         LBETA = C + QTR 
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA 
         F = DLAMC3(B/2,(-B/100)) 
         C = DLAMC3(F,A) 
         IF (C == A) THEN 
            LRND = .TRUE. 
         ELSE 
            LRND = .FALSE. 
         ENDIF 
         F = DLAMC3(B/2,B/100) 
         C = DLAMC3(F,A) 
         IF (LRND .AND. C==A) LRND = .FALSE. 
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = DLAMC3(B/2,A) 
         T2 = DLAMC3(B/2,SAVEC) 
         LIEEE1 = T1==A .AND. T2>SAVEC .AND. LRND 
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0 
         A = 1 
         C = 1 
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE 
         IF (C == ONE) THEN 
            LT = LT + 1 
            A = A*LBETA 
            C = DLAMC3(A,ONE) 
            C = DLAMC3(C,(-A)) 
            GO TO 30 
         ENDIF 
!+       END WHILE
!
      ENDIF 
!
      BETA = LBETA 
      T = LT 
      RND = LRND 
      IEEE1 = LIEEE1 
      RETURN  
!
!     End of DLAMC1
!
      END SUBROUTINE DLAMC1 
