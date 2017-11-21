      SUBROUTINE DROTI(NZ, X, INDX, Y, C, S) 
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
      REAL(DOUBLE) , INTENT(IN) :: C 
      REAL(DOUBLE) , INTENT(IN) :: S 
      INTEGER , INTENT(IN) :: INDX(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: X(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: TEMP 
!-----------------------------------------------
!
!     ==================================================================
!     ==================================================================
!     ====  DROTI  --  APPLY INDEXED DOUBLE PRECISION GIVENS        ====
!     ====             ROTATION                                     ====
!     ==================================================================
!     ==================================================================
!
!     -------------
!     ... ARGUMENTS
!     -------------
!
!     PURPOSE
!     -------
!
!         DROTI APPLIES A GIVENS ROTATION TO
!             A SPARSE VECTOR   X  STORED IN COMPRESSED FORM  (X,INDX)
!         AND
!             ANOTHER VECTOR  Y  IN FULL STORAGE FORM.
!
!         DROTI DOES NOT HANDLE FILL-IN IN  X  AND THEREFORE, IT IS
!         ASSUMED THAT ALL NONZERO COMPONENTS OF  Y  ARE LISTED IN
!         INDX.  ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN
!         INDX  ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST
!         BE DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
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
!         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
!         INDX    INTEGER     ARRAY CONTAINING THE INDICIES OF THE
!                             COMPRESSED FORM.  IT IS ASSUMED THAT
!                             THE ELEMENTS IN  INDX  ARE DISTINCT.
!         C,S     DOUBLE      THE TWO SCALARS DEFINING THE GIVENS
!                             ROTATION.
!
!     UPDATED ...
!
!         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE
!                             SPARSE VECTOR IN COMPRESSED FORM.
!         Y       DOUBLE      ARRAY WHICH CONTAINS THE VECTOR  Y
!                             IN FULL STORAGE FORM.  ONLY THE
!                             ELEMENTS WHOSE INDICIES ARE LISTED IN
!                             INDX  HAVE BEEN REFERENCED OR MODIFIED.
!
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
!
!     ==================================================================
!
      IF (NZ <= 0) RETURN  
!
      IF (C==1.0D0 .AND. S==0.0D0) RETURN  
!
      DO I = 1, NZ 
         TEMP = (-S*X(I)) + C*Y(INDX(I)) 
         X(I) = C*X(I) + S*Y(INDX(I)) 
         Y(INDX(I)) = TEMP 
      END DO 
!
      RETURN  
      END SUBROUTINE DROTI 
