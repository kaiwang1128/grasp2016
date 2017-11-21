      REAL(KIND(0.0D0)) FUNCTION DLANSP (NORM, UPLO, N, AP, WORK) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:31   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dlassq_I 
      !!USE lsame_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      CHARACTER  :: NORM 
      CHARACTER  :: UPLO 
      REAL(DOUBLE)  :: AP(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K 
      REAL(DOUBLE) :: ABSA, SCALE, SUM, VALUE 
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
!     October 31, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLANSP  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A,  supplied in packed form.
!
!  Description
!  ===========
!
!  DLANSP returns the value
!
!     DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
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
!          Specifies the value to be returned in DLANSP as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is supplied.
!          = 'U':  Upper triangular part of A is supplied
!          = 'L':  Lower triangular part of A is supplied
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is
!          set to zero.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The upper or lower triangle of the symmetric matrix A, packed
!          columnwise in a linear array.  The j-th column of A is stored
!          in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF (N == 0) THEN 
         VALUE = ZERO 
      ELSE IF (LSAME(NORM,'M')) THEN 
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO 
         IF (LSAME(UPLO,'U')) THEN 
            K = 1 
            DO J = 1, N 
               DO I = K, K + J - 1 
                  VALUE = MAX(VALUE,ABS(AP(I))) 
               END DO 
               K = K + J 
            END DO 
         ELSE 
            K = 1 
            DO J = 1, N 
               DO I = K, K + N - J 
                  VALUE = MAX(VALUE,ABS(AP(I))) 
               END DO 
               K = K + N - J + 1 
            END DO 
         ENDIF 
      ELSE IF (LSAME(NORM,'I') .OR. LSAME(NORM,'O') .OR. NORM=='1') THEN 
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO 
         K = 1 
         IF (LSAME(UPLO,'U')) THEN 
            DO J = 1, N 
               SUM = ZERO 
               DO I = 1, J - 1 
                  ABSA = ABS(AP(K)) 
                  SUM = SUM + ABSA 
                  WORK(I) = WORK(I) + ABSA 
                  K = K + 1 
               END DO 
               WORK(J) = SUM + ABS(AP(K)) 
               K = K + 1 
            END DO 
            DO I = 1, N 
               VALUE = MAX(VALUE,WORK(I)) 
            END DO 
         ELSE 
            WORK(:N) = ZERO 
            DO J = 1, N 
               SUM = WORK(J) + ABS(AP(K)) 
               K = K + 1 
               DO I = J + 1, N 
                  ABSA = ABS(AP(K)) 
                  SUM = SUM + ABSA 
                  WORK(I) = WORK(I) + ABSA 
                  K = K + 1 
               END DO 
               VALUE = MAX(VALUE,SUM) 
            END DO 
         ENDIF 
      ELSE IF (LSAME(NORM,'F') .OR. LSAME(NORM,'E')) THEN 
!
!        Find normF(A).
!
         SCALE = ZERO 
         SUM = ONE 
         K = 2 
         IF (LSAME(UPLO,'U')) THEN 
            DO J = 2, N 
               CALL DLASSQ (J - 1, AP(K), 1, SCALE, SUM) 
               K = K + J 
            END DO 
         ELSE 
            DO J = 1, N - 1 
               CALL DLASSQ (N - J, AP(K), 1, SCALE, SUM) 
               K = K + N - J + 1 
            END DO 
         ENDIF 
         SUM = 2*SUM 
         K = 1 
         DO I = 1, N 
            IF (AP(K) /= ZERO) THEN 
               ABSA = ABS(AP(K)) 
               IF (SCALE < ABSA) THEN 
                  SUM = ONE + SUM*(SCALE/ABSA)**2 
                  SCALE = ABSA 
               ELSE 
                  SUM = SUM + (ABSA/SCALE)**2 
               ENDIF 
            ENDIF 
            IF (LSAME(UPLO,'U')) THEN 
               K = K + I + 1 
            ELSE 
               K = K + N - I + 1 
            ENDIF 
         END DO 
         VALUE = SCALE*SQRT(SUM) 
      ENDIF 
!
      DLANSP = VALUE 
      RETURN  
!
!     End of DLANSP
!
      END FUNCTION DLANSP 
