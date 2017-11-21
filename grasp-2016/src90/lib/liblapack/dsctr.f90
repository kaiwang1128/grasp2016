      SUBROUTINE DSCTR(NZ, X, INDX, Y) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:35   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NZ 
      INTEGER , INTENT(IN) :: INDX(*) 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(OUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ==================================================================
!     ==================================================================
!     ====  DSCTR -- DOUBLE PRECISION SCATTER                       ====
!     ==================================================================
!     ==================================================================
!
!     PURPOSE
!     -------
!
!         DSCTR SCATTERS THE COMPONENTS OF
!             A SPARSE VECTOR  X  STORED IN COMPRESSED FORM  (X,INDX)
!         INTO
!             SPECIFIED COMPONENTS OF A DOUBLE PRECISION VECTOR  Y
!             IN FULL STORAGE FORM.
!
!         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX
!         ARE MODIFIED.  THE VALUES IN  INDX  MUST BE DISTINCT TO
!         ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
!
!         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
!         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
!         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME
!         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
!
!     ARGUMENTS
!     ---------
!
!     INPUT ...
!
!         NZ      INTEGER     NUMBER OF ELEMENTS TO BE SCATTERED FROM
!                             COMPRESSED FORM.
!         X       DOUBLE      ARRAY CONTAINING THE VALUES TO BE
!                             SCATTERED FROM COMPRESSED FORM INTO FULL
!                             STORAGE FORM.
!         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE VALUES
!                             TO BE SCATTERED FROM COMPRESSED FORM.
!                             IT IS ASSUMED THAT THE ELEMENTS IN  INDX
!                             ARE DISTINCT.
!
!     OUTPUT ...
!
!         Y       DOUBLE      ARRAY WHOSE ELEMENTS SPECIFIED BY  INDX
!                             HAVE BEEN SET TO THE CORRESPONDING
!                             ENTRIES OF  X.  ONLY THE ELEMENTS
!                             CORRESPONDING TO THE INDICES IN  INDX
!                             HAVE BEEN MODIFIED.
!
!     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
!
!     FORTRAN VERSION WRITTEN OCTOBER 1984
!     ROGER G GRIMES, BOEING COMPUTER SERVICES
!
!     ==================================================================
!
!     -------------
!     ... ARGUMENTS
!     -------------
!
!
!
!     -------------------
!     ... LOCAL VARIABLES
!     -------------------
!
!
!     ==================================================================
!
      IF (NZ <= 0) RETURN  
!
      DO I = 1, NZ 
         Y(INDX(I)) = X(I) 
      END DO 
!
      RETURN  
      END SUBROUTINE DSCTR 
