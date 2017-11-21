      SUBROUTINE DGTHR(NZ, Y, X, INDX) 
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
      REAL(DOUBLE) , INTENT(IN) :: Y(*) 
      REAL(DOUBLE) , INTENT(OUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
!
!     ==================================================================
!     ==================================================================
!     ====  DGTHR -- DOUBLE PRECISION GATHER                        ====
!     ==================================================================
!     ==================================================================
!
!     PURPOSE
!     -------
!
!         DGTHR GATHERS THE SPECIFIED ELEMENTS FROM
!             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM
!         INTO
!             A DOUBLE PRECISION VECTOR  X  IN COMPRESSED FORM (X,INDX).
!
!         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX
!         ARE REFERENCED.
!
!     ARGUMENTS
!     ---------
!
!     INPUT ...
!
!         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO
!                             COMPRESSED FORM.
!         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE
!                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
!                             THE ELEMENTS CORRESPONDING TO THE INDICES
!                             IN  INDX  WILL BE ACCESSED.
!         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
!                             VALUES TO BE GATHERED INTO COMPRESSED FORM.
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
      X(:NZ) = Y(INDX(:NZ)) 
!
      RETURN  
      END SUBROUTINE DGTHR 
