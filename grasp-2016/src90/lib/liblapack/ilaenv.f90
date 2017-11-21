 
      INTEGER FUNCTION ILAENV (ISPEC, NAME, OPTS, N1, N2, N3, N4) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:40   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISPEC 
      INTEGER , INTENT(IN) :: N1 
      INTEGER , INTENT(IN) :: N2 
      INTEGER  :: N3 
      INTEGER , INTENT(IN) :: N4 
      CHARACTER , INTENT(IN) :: NAME*(*) 
      CHARACTER  :: OPTS*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IC, IZ, NB, NBMIN, NX 
      LOGICAL :: CNAME, SNAME 
      CHARACTER :: C1, C2*2, C4*2, C3*3, SUBNAM*6 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC CHAR, ICHAR, INT, MIN, REAL 
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
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      SELECT CASE (ISPEC)  
!
!     Invalid value for ISPEC
!
      CASE DEFAULT 
         ILAENV = -1 
         RETURN  
!
      CASE (1:3)  
         ILAENV = 1 
         SUBNAM = NAME 
         IC = ICHAR(SUBNAM(1:1)) 
         IZ = ICHAR('Z') 
         IF (IZ==90 .OR. IZ==122) THEN 
!
!        ASCII character set
!
            IF (IC>=97 .AND. IC<=122) THEN 
               SUBNAM(1:1) = CHAR(IC - 32) 
               DO I = 2, 6 
                  IC = ICHAR(SUBNAM(I:I)) 
                  IF (IC<97 .OR. IC>122) CYCLE  
                  SUBNAM(I:I) = CHAR(IC - 32) 
               END DO 
            ENDIF 
!
         ELSE IF (IZ==233 .OR. IZ==169) THEN 
!
!        EBCDIC character set
!
            IF (IC>=129 .AND. IC<=137 .OR. IC>=145 .AND. IC<=153 .OR. IC>=162&
                .AND. IC<=169) THEN 
               SUBNAM(1:1) = CHAR(IC + 64) 
               DO I = 2, 6 
                  IC = ICHAR(SUBNAM(I:I)) 
                  IF (.NOT.(IC>=129 .AND. IC<=137 .OR. IC>=145 .AND. IC<=153&
                      .OR. IC>=162 .AND. IC<=169)) CYCLE  
                  SUBNAM(I:I) = CHAR(IC + 64) 
               END DO 
            ENDIF 
!
         ELSE IF (IZ==218 .OR. IZ==250) THEN 
!
!        Prime machines:  ASCII+128
!
            IF (IC>=225 .AND. IC<=250) THEN 
               SUBNAM(1:1) = CHAR(IC - 32) 
               DO I = 2, 6 
                  IC = ICHAR(SUBNAM(I:I)) 
                  IF (IC<225 .OR. IC>250) CYCLE  
                  SUBNAM(I:I) = CHAR(IC - 32) 
               END DO 
            ENDIF 
         ENDIF 
!
         C1 = SUBNAM(1:1) 
         SNAME = C1=='S' .OR. C1=='D' 
         CNAME = C1=='C' .OR. C1=='Z' 
         IF (.NOT.(CNAME .OR. SNAME)) RETURN  
         C2 = SUBNAM(2:3) 
         C3 = SUBNAM(4:6) 
         C4 = C3(2:3) 
!
         SELECT CASE (ISPEC)  
!
         CASE DEFAULT 
            NB = 1 
!
            IF (C2 == 'GE') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ELSE IF (C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. C3=='QLF'&
                     ) THEN 
                  IF (SNAME) THEN 
                     NB = 32 
                  ELSE 
                     NB = 32 
                  ENDIF 
               ELSE IF (C3 == 'HRD') THEN 
                  IF (SNAME) THEN 
                     NB = 32 
                  ELSE 
                     NB = 32 
                  ENDIF 
               ELSE IF (C3 == 'BRD') THEN 
                  IF (SNAME) THEN 
                     NB = 32 
                  ELSE 
                     NB = 32 
                  ENDIF 
               ELSE IF (C3 == 'TRI') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'PO') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'SY') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ELSE IF (SNAME .AND. C3=='TRD') THEN 
                  NB = 1 
               ELSE IF (SNAME .AND. C3=='GST') THEN 
                  NB = 64 
               ENDIF 
            ELSE IF (CNAME .AND. C2=='HE') THEN 
               SELECT CASE (C3)  
               CASE ('TRF')  
                  NB = 64 
               CASE ('TRD')  
                  NB = 1 
               CASE ('GST')  
                  NB = 64 
               END SELECT 
            ELSE IF (SNAME .AND. C2=='OR') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NB = 32 
               ELSE IF (C3(1:1) == 'M') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NB = 32 
               ENDIF 
            ELSE IF (CNAME .AND. C2=='UN') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NB = 32 
               ELSE IF (C3(1:1) == 'M') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NB = 32 
               ENDIF 
            ELSE IF (C2 == 'GB') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     IF (N4 <= 64) THEN 
                        NB = 1 
                     ELSE 
                        NB = 32 
                     ENDIF 
                  ELSE 
                     IF (N4 <= 64) THEN 
                        NB = 1 
                     ELSE 
                        NB = 32 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'PB') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     IF (N2 <= 64) THEN 
                        NB = 1 
                     ELSE 
                        NB = 32 
                     ENDIF 
                  ELSE 
                     IF (N2 <= 64) THEN 
                        NB = 1 
                     ELSE 
                        NB = 32 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'TR') THEN 
               IF (C3 == 'TRI') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'LA') THEN 
               IF (C3 == 'UUM') THEN 
                  IF (SNAME) THEN 
                     NB = 64 
                  ELSE 
                     NB = 64 
                  ENDIF 
               ENDIF 
            ELSE IF (SNAME .AND. C2=='ST') THEN 
               IF (C3 == 'EBZ') NB = 1 
            ENDIF 
            ILAENV = NB 
            RETURN  
!
         CASE (2)  
            NBMIN = 2 
            IF (C2 == 'GE') THEN 
               IF (C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. C3=='QLF') THEN 
                  IF (SNAME) THEN 
                     NBMIN = 2 
                  ELSE 
                     NBMIN = 2 
                  ENDIF 
               ELSE IF (C3 == 'HRD') THEN 
                  IF (SNAME) THEN 
                     NBMIN = 2 
                  ELSE 
                     NBMIN = 2 
                  ENDIF 
               ELSE IF (C3 == 'BRD') THEN 
                  IF (SNAME) THEN 
                     NBMIN = 2 
                  ELSE 
                     NBMIN = 2 
                  ENDIF 
               ELSE IF (C3 == 'TRI') THEN 
                  IF (SNAME) THEN 
                     NBMIN = 2 
                  ELSE 
                     NBMIN = 2 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'SY') THEN 
               IF (C3 == 'TRF') THEN 
                  IF (SNAME) THEN 
                     NBMIN = 8 
                  ELSE 
                     NBMIN = 8 
                  ENDIF 
               ELSE IF (SNAME .AND. C3=='TRD') THEN 
                  NBMIN = 2 
               ENDIF 
            ELSE IF (CNAME .AND. C2=='HE') THEN 
               IF (C3 == 'TRD') NBMIN = 2 
            ELSE IF (SNAME .AND. C2=='OR') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NBMIN = 2 
               ELSE IF (C3(1:1) == 'M') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NBMIN = 2 
               ENDIF 
            ELSE IF (CNAME .AND. C2=='UN') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NBMIN = 2 
               ELSE IF (C3(1:1) == 'M') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NBMIN = 2 
               ENDIF 
            ENDIF 
            ILAENV = NBMIN 
            RETURN  
!
         CASE (3)  
            NX = 0 
            IF (C2 == 'GE') THEN 
               IF (C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. C3=='QLF') THEN 
                  IF (SNAME) THEN 
                     NX = 128 
                  ELSE 
                     NX = 128 
                  ENDIF 
               ELSE IF (C3 == 'HRD') THEN 
                  IF (SNAME) THEN 
                     NX = 128 
                  ELSE 
                     NX = 128 
                  ENDIF 
               ELSE IF (C3 == 'BRD') THEN 
                  IF (SNAME) THEN 
                     NX = 128 
                  ELSE 
                     NX = 128 
                  ENDIF 
               ENDIF 
            ELSE IF (C2 == 'SY') THEN 
               IF (SNAME .AND. C3=='TRD') NX = 1 
            ELSE IF (CNAME .AND. C2=='HE') THEN 
               IF (C3 == 'TRD') NX = 1 
            ELSE IF (SNAME .AND. C2=='OR') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NX = 128 
               ENDIF 
            ELSE IF (CNAME .AND. C2=='UN') THEN 
               IF (C3(1:1) == 'G') THEN 
                  IF (C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. C4=='QL' .OR. &
                     C4=='HR' .OR. C4=='TR' .OR. C4=='BR') NX = 128 
               ENDIF 
            ENDIF 
            ILAENV = NX 
            RETURN  
!
         CASE (4)  
            ILAENV = 6 
            RETURN  
!
         CASE (5)  
            ILAENV = 2 
            RETURN  
!
         CASE (6)  
            ILAENV = INT(REAL(MIN(N1,N2))*1.6E0) 
            RETURN  
!
         CASE (7)  
            ILAENV = 1 
            RETURN  
!
         CASE (8)  
            ILAENV = 50 
            RETURN  
         END SELECT 
      END SELECT 
!
!     End of ILAENV
!
      END FUNCTION ILAENV 
