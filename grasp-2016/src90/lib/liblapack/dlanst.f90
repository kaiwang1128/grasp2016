      REAL(KIND(0.0D0)) FUNCTION DLANST (NORM, N, D, E) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:31   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USE dlassq_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      CHARACTER  :: NORM 
      REAL(DOUBLE)  :: D(*) 
      REAL(DOUBLE)  :: E(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: ANORM, SCALE, SUM 
      logical::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ABS, MAX, SQRT 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  DLANST returns the value
!
!     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!          set to zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF (N <= 0) THEN 
         ANORM = ZERO 
      ELSE IF (LSAME(NORM,'M')) THEN 
!
!        Find max(abs(A(i,j))).
!
         ANORM = ABS(D(N)) 
         DO I = 1, N - 1 
            ANORM = MAX(ANORM,ABS(D(I))) 
            ANORM = MAX(ANORM,ABS(E(I))) 
         END DO 
      ELSE IF (LSAME(NORM,'O') .OR. NORM=='1' .OR. LSAME(NORM,'I')) THEN 
!
!        Find norm1(A).
!
         IF (N == 1) THEN 
            ANORM = ABS(D(1)) 
         ELSE 
            ANORM = MAX(ABS(D(1))+ABS(E(1)),ABS(E(N-1))+ABS(D(N))) 
            DO I = 2, N - 1 
               ANORM = MAX(ANORM,ABS(D(I))+ABS(E(I))+ABS(E(I-1))) 
            END DO 
         ENDIF 
      ELSE IF (LSAME(NORM,'F') .OR. LSAME(NORM,'E')) THEN 
!
!        Find normF(A).
!
         SCALE = ZERO 
         SUM = ONE 
         IF (N > 1) THEN 
            CALL DLASSQ (N - 1, E, 1, SCALE, SUM) 
            SUM = 2*SUM 
         ENDIF 
         CALL DLASSQ (N, D, 1, SCALE, SUM) 
         ANORM = SCALE*SQRT(SUM) 
      ENDIF 
!
      DLANST = ANORM 
      RETURN  
!
!     End of DLANST
!
      END FUNCTION DLANST 
