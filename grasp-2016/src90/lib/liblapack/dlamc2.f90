!
!***********************************************************************
!
      SUBROUTINE DLAMC2(BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX) 
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
      !USE dlamc1_I 
      !USE dlamc4_I 
      !USE dlamc5_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: BETA 
      INTEGER , INTENT(OUT) :: T 
      INTEGER , INTENT(OUT) :: EMIN 
      INTEGER , INTENT(OUT) :: EMAX 
      REAL(DOUBLE) , INTENT(OUT) :: EPS 
      REAL(DOUBLE) , INTENT(OUT) :: RMIN 
      REAL(DOUBLE) , INTENT(OUT) :: RMAX 
      real(double) :: dlamc3
      LOGICAL , INTENT(OUT) :: RND 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, NGNMIN, NGPMIN 
      REAL(DOUBLE) :: A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, SIXTH, &
         SMALL, THIRD, TWO, ZERO 
      LOGICAL :: FIRST, IEEE, IWARN, LIEEE1, LRND 

      SAVE FIRST, IWARN, LBETA, LEMAX, LEMIN, LT, LEPS, LRMAX, LRMIN 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, MIN 
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
!  DLAMC2 determines the machine parameters specified in its argument
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
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
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
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1 (LBETA, LT, LRND, LIEEE1) 
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
         SIXTH = DLAMC3(B,(-HALF)) 
         THIRD = DLAMC3(SIXTH,SIXTH) 
         B = DLAMC3(THIRD,(-HALF)) 
         B = DLAMC3(B,SIXTH) 
         B = ABS(B) 
         B = DMAX1(LEPS,B) 
!
         LEPS = 1 
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE 
         IF (LEPS>B .AND. B>ZERO) THEN 
            LEPS = B 
            C = DLAMC3(HALF*LEPS,TWO**5*LEPS**2) 
            C = DLAMC3(HALF,(-C)) 
            B = DLAMC3(HALF,C) 
            C = DLAMC3(HALF,(-B)) 
            B = DLAMC3(HALF,C) 
            GO TO 10 
         ENDIF 
!+       END WHILE
!
         LEPS = MIN(A,LEPS) 
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
            SMALL = DLAMC3(SMALL*RBASE,ZERO) 
         END DO 
         A = DLAMC3(ONE,SMALL) 
         CALL DLAMC4 (NGPMIN, ONE, LBETA) 
         CALL DLAMC4 (NGNMIN, (-ONE), LBETA) 
         CALL DLAMC4 (GPMIN, A, LBETA) 
         CALL DLAMC4 (GNMIN, (-A), LBETA) 
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
!        in routine DLAMC1. A true IEEE machine should have both  things
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
            LRMIN = DLAMC3(LRMIN*RBASE,ZERO) 
         END DO 
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL DLAMC5 (LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX) 
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
         ' the IF block as marked within the code of routine',' DLAMC2,',/,&
         ' otherwise supply EMIN explicitly.',/) 
      RETURN  
!
!     End of DLAMC2
!
      END SUBROUTINE DLAMC2 
