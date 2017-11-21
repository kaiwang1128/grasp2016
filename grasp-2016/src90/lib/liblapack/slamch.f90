      REAL FUNCTION SLAMCH (CMACH) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USE slamc2_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: CMACH 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL, PARAMETER :: ONE = 1.0E+0 
      REAL, PARAMETER :: ZERO = 0.0E+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: BETA, IMAX, IMIN, IT 
      REAL :: BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, RND, SFMIN, SMALL&
         , T 
      LOGICAL :: FIRST, LRND 

      LOGICAL :: LSAME
      real :: slamc3

      SAVE FIRST, BASE, EMAX, EMIN, EPS, PREC, RMAX, RMIN, RND, SFMIN, T 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
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
         CALL SLAMC2 (BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX) 
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
      SLAMCH = RMACH 
      RETURN  
!
!     End of SLAMCH
!
      END FUNCTION SLAMCH 


!
!***********************************************************************
!
      SUBROUTINE SLAMC1(BETA, T, RND, IEEE1) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE slamc3_I 
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
      REAL :: A, B, C, F, ONE, QTR, SAVEC, T1, T2 
      LOGICAL :: FIRST, LIEEE1, LRND 
      real :: slamc3

      SAVE FIRST, LIEEE1, LRND, LBETA, LT 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
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
!        Throughout this routine  we use the function  SLAMC3  to ensure
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
            C = SLAMC3(A,ONE) 
            C = SLAMC3(C,(-A)) 
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
         C = SLAMC3(A,B) 
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE 
         IF (C == A) THEN 
            B = 2*B 
            C = SLAMC3(A,B) 
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
         C = SLAMC3(C,(-A)) 
         LBETA = C + QTR 
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA 
         F = SLAMC3(B/2,(-B/100)) 
         C = SLAMC3(F,A) 
         IF (C == A) THEN 
            LRND = .TRUE. 
         ELSE 
            LRND = .FALSE. 
         ENDIF 
         F = SLAMC3(B/2,B/100) 
         C = SLAMC3(F,A) 
         IF (LRND .AND. C==A) LRND = .FALSE. 
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3(B/2,A) 
         T2 = SLAMC3(B/2,SAVEC) 
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
            C = SLAMC3(A,ONE) 
            C = SLAMC3(C,(-A)) 
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
!     End of SLAMC1
!
      END SUBROUTINE SLAMC1 


!
!***********************************************************************
!
      SUBROUTINE SLAMC2(BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE slamc3_I 
      !USE slamc1_I 
      !USE slamc4_I 
      !USE slamc5_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: BETA 
      INTEGER , INTENT(OUT) :: T 
      INTEGER , INTENT(OUT) :: EMIN 
      INTEGER , INTENT(OUT) :: EMAX 
      REAL , INTENT(OUT) :: EPS 
      REAL , INTENT(OUT) :: RMIN 
      REAL , INTENT(OUT) :: RMAX 
      LOGICAL , INTENT(OUT) :: RND 
      real :: slamc3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, NGNMIN, NGPMIN 
      REAL :: A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, SIXTH, SMALL, &
         THIRD, TWO, ZERO 
      LOGICAL :: FIRST, IEEE, IWARN, LIEEE1, LRND 

      SAVE FIRST, IWARN, LBETA, LEMAX, LEMIN, LT, LEPS, LRMAX, LRMIN 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, MIN 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
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
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Save statement ..
!     ..
!     .. Data statements ..
      DATA FIRST/ .TRUE./  
      DATA IWARN/ .FALSE./  
!     ..
!     .. Executable Statements ..
!
      IF (FIRST) THEN 
         FIRST = .FALSE. 
         ZERO = 0 
         ONE = 1 
         TWO = 2 
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1 (LBETA, LT, LRND, LIEEE1) 
!
!        Start to find EPS.
!
         B = LBETA 
         A = B**(-LT) 
         LEPS = A 
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO/3 
         HALF = ONE/2 
         SIXTH = SLAMC3(B,(-HALF)) 
         THIRD = SLAMC3(SIXTH,SIXTH) 
         B = SLAMC3(THIRD,(-HALF)) 
         B = SLAMC3(B,SIXTH) 
         B = ABS(B) 
         B = AMAX1(LEPS,B) 
!
         LEPS = 1 
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE 
         IF (LEPS>B .AND. B>ZERO) THEN 
            LEPS = B 
            C = SLAMC3(HALF*LEPS,TWO**5*LEPS**2) 
            C = SLAMC3(HALF,(-C)) 
            B = SLAMC3(HALF,C) 
            C = SLAMC3(HALF,(-B)) 
            B = SLAMC3(HALF,C) 
            GO TO 10 
         ENDIF 
!+       END WHILE
!
         LEPS = AMIN1(A,LEPS) 
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE/LBETA 
         SMALL = ONE 
         DO I = 1, 3 
            SMALL = SLAMC3(SMALL*RBASE,ZERO) 
         END DO 
         A = SLAMC3(ONE,SMALL) 
         CALL SLAMC4 (NGPMIN, ONE, LBETA) 
         CALL SLAMC4 (NGNMIN, (-ONE), LBETA) 
         CALL SLAMC4 (GPMIN, A, LBETA) 
         CALL SLAMC4 (GNMIN, (-A), LBETA) 
         IEEE = .FALSE. 
!
         IF (NGPMIN==NGNMIN .AND. GPMIN==GNMIN) THEN 
            IF (NGPMIN == GPMIN) THEN 
               LEMIN = NGPMIN 
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF (GPMIN - NGPMIN == 3) THEN 
               LEMIN = NGPMIN - 1 + LT 
               IEEE = .TRUE. 
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE 
               LEMIN = MIN(NGPMIN,GPMIN) 
!            ( A guess; no known machine )
               IWARN = .TRUE. 
            ENDIF 
!
         ELSE IF (NGPMIN==GPMIN .AND. NGNMIN==GNMIN) THEN 
            IF (ABS(NGPMIN - NGNMIN) == 1) THEN 
               LEMIN = MAX(NGPMIN,NGNMIN) 
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE 
               LEMIN = MIN(NGPMIN,NGNMIN) 
!            ( A guess; no known machine )
               IWARN = .TRUE. 
            ENDIF 
!
         ELSE IF (ABS(NGPMIN - NGNMIN)==1 .AND. GPMIN==GNMIN) THEN 
            IF (GPMIN - MIN(NGPMIN,NGNMIN) == 3) THEN 
               LEMIN = MAX(NGPMIN,NGNMIN) - 1 + LT 
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE 
               LEMIN = MIN(NGPMIN,NGNMIN) 
!            ( A guess; no known machine )
               IWARN = .TRUE. 
            ENDIF 
!
         ELSE 
            LEMIN = MIN(MIN(MIN(NGPMIN,NGNMIN),GPMIN),GNMIN) 
!         ( A guess; no known machine )
            IWARN = .TRUE. 
         ENDIF 
!**
! Comment out this if block if EMIN is ok
         IF (IWARN) THEN 
            FIRST = .TRUE. 
            WRITE (6, FMT=9999) LEMIN 
         ENDIF 
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1 
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1 
         DO I = 1, 1 - LEMIN 
            LRMIN = SLAMC3(LRMIN*RBASE,ZERO) 
         END DO 
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5 (LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX) 
      ENDIF 
!
      BETA = LBETA 
      T = LT 
      RND = LRND 
      EPS = LEPS 
      EMIN = LEMIN 
      RMIN = LRMIN 
      EMAX = LEMAX 
      RMAX = LRMAX 
!
      RETURN  
!
 9999 FORMAT(/,/,' WARNING. The value EMIN may be incorrect:-','  EMIN = ',I8,/&
         ,' If, after inspection, the value EMIN looks',&
         ' acceptable please comment out ',/,&
         ' the IF block as marked within the code of routine',' SLAMC2,',/,&
         ' otherwise supply EMIN explicitly.',/) 
      RETURN  
!
!     End of SLAMC2
!
      END SUBROUTINE SLAMC2 


!
!***********************************************************************
!
      REAL FUNCTION SLAMC3 (A, B) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL , INTENT(IN) :: A 
      REAL , INTENT(IN) :: B 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B 
!
      RETURN  
!
!     End of SLAMC3
!
      END FUNCTION SLAMC3 


!
!***********************************************************************
!
      SUBROUTINE SLAMC4(EMIN, START, BASE) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE slamc3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: EMIN 
      INTEGER , INTENT(IN) :: BASE 
      REAL , INTENT(IN) :: START 
      real :: slamc3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL :: A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
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
      B1 = SLAMC3(A*RBASE,ZERO) 
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
         B1 = SLAMC3(A/BASE,ZERO) 
         C1 = SLAMC3(B1*BASE,ZERO) 
         D1 = ZERO 
         IF (BASE > 0) THEN 
            D1 = (BASE - 1)*B1 + D1 + B1 
         ENDIF 
         B2 = SLAMC3(A*RBASE,ZERO) 
         C2 = SLAMC3(B2/RBASE,ZERO) 
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
!     End of SLAMC4
!
      END SUBROUTINE SLAMC4 


!
!***********************************************************************
!
      SUBROUTINE SLAMC5(BETA, P, EMIN, IEEE, EMAX, RMAX) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:41   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE slamc3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: BETA 
      INTEGER , INTENT(IN) :: P 
      INTEGER , INTENT(IN) :: EMIN 
      INTEGER , INTENT(OUT) :: EMAX 
      REAL , INTENT(OUT) :: RMAX 
      LOGICAL , INTENT(IN) :: IEEE 
      real :: slamc3
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL, PARAMETER :: ZERO = 0.0E0 
      REAL, PARAMETER :: ONE = 1.0E0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP 
      REAL :: OLDY, RECBAS, Y, Z 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MOD 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
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
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
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
!  RMAX    (output) REAL
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
         Y = SLAMC3(Y,Z) 
      END DO 
      IF (Y >= ONE) Y = OLDY 
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO I = 1, EMAX 
         Y = SLAMC3(Y*BETA,ZERO) 
      END DO 
!
      RMAX = Y 
      RETURN  
!
!     End of SLAMC5
!
      END SUBROUTINE SLAMC5 
