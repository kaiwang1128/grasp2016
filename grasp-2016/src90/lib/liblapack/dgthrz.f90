      SUBROUTINE DGTHRZ(NZ, Y, X, INDX) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:29   2/12/04  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NZ 
      INTEGER , INTENT(IN) :: INDX(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: Y(*) 
      REAL(DOUBLE) , INTENT(OUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ==================================================================
!     ==================================================================
!     ====  DGTHRZ -- DOUBLE PRECISION GATHER AND ZERO              ====
!     ==================================================================
!     ==================================================================
!
!     PURPOSE
!     -------
!
!         DGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM
!             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM
!         INTO
!             A DOUBLE PRECISION VECTOR  X  IN COMPRESSED FORM  (X,INDX).
!         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
!
!         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX
!         ARE REFERENCED OR MODIFIED.
!
!     ARGUMENTS
!     ---------
!
!     INPUT ...
!
!         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO
!                             COMPRESSED FORM.
!         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
!                             VALUES TO BE GATHERED INTO COMPRESSED FORM.
!
!     UPDATED ...
!
!         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE
!                             VECTOR  Y  IN FULL STORAGE FORM.  THE
!                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.
!                             ONLY THE ELEMENTS CORRESPONDING TO THE
!                             INDICES IN  INDX  HAVE BEEN ACCESSED.
!
!     OUTPUT ...
!
!         X       DOUBLE      ARRAY CONTAINING THE VALUES GATHERED INTO
!                             THE COMPRESSED FORM.
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
         X(I) = Y(INDX(I)) 
         Y(INDX(I)) = 0.0D0 
      END DO 
!
      RETURN  
      END SUBROUTINE DGTHRZ 
