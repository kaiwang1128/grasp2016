      SUBROUTINE DSYR2K(UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:38   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE lsame_I 
      !USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: LDB 
      INTEGER , INTENT(IN) :: LDC 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: UPLO 
      CHARACTER  :: TRANS 
      REAL(DOUBLE) , INTENT(IN) :: A(LDA,*) 
      REAL(DOUBLE) , INTENT(IN) :: B(LDB,*) 
      REAL(DOUBLE) , INTENT(INOUT) :: C(LDC,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, J, L, NROWA 
      REAL(DOUBLE) :: TEMP1, TEMP2 
      LOGICAL :: UPPER 
      logical::lsame
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSYR2K  performs one of the symmetric rank 2k operations
!
!     C := alpha*A*B' + alpha*B*A' + beta*C,
!
!  or
!
!     C := alpha*A'*B + alpha*B'*A + beta*C,
!
!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
!  matrices in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
!                                        beta*C.
!
!              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns  of the  matrices  A and B,  and on  entry  with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrices  A and B.  K must be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  k by n  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDB must be at least  max( 1, n ), otherwise  LDB must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!     .. External Subroutines ..
!     .. Intrinsic Functions ..
!     .. Local Scalars ..
!     .. Parameters ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      IF (LSAME(TRANS,'N')) THEN 
         NROWA = N 
      ELSE 
         NROWA = K 
      ENDIF 
      UPPER = LSAME(UPLO,'U') 
!
      INFO = 0 
      IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = 1 
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.&
            LSAME(TRANS,'C')) THEN 
         INFO = 2 
      ELSE IF (N < 0) THEN 
         INFO = 3 
      ELSE IF (K < 0) THEN 
         INFO = 4 
      ELSE IF (LDA < MAX(1,NROWA)) THEN 
         INFO = 7 
      ELSE IF (LDB < MAX(1,NROWA)) THEN 
         INFO = 9 
      ELSE IF (LDC < MAX(1,N)) THEN 
         INFO = 12 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSYR2K', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N==0 .OR. (ALPHA==ZERO .OR. K==0) .AND. BETA==ONE) RETURN  
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA == ZERO) THEN 
         IF (UPPER) THEN 
            IF (BETA == ZERO) THEN 
               DO J = 1, N 
                  C(:J,J) = ZERO 
               END DO 
            ELSE 
               DO J = 1, N 
                  C(:J,J) = BETA*C(:J,J) 
               END DO 
            ENDIF 
         ELSE 
            IF (BETA == ZERO) THEN 
               DO J = 1, N 
                  C(J:N,J) = ZERO 
               END DO 
            ELSE 
               DO J = 1, N 
                  C(J:N,J) = BETA*C(J:N,J) 
               END DO 
            ENDIF 
         ENDIF 
         RETURN  
      ENDIF 
!
!     Start the operations.
!
      IF (LSAME(TRANS,'N')) THEN 
!
!        Form  C := alpha*A*B' + alpha*B*A' + C.
!
         IF (UPPER) THEN 
            DO J = 1, N 
               IF (BETA == ZERO) THEN 
                  C(:J,J) = ZERO 
               ELSE IF (BETA /= ONE) THEN 
                  C(:J,J) = BETA*C(:J,J) 
               ENDIF 
               DO L = 1, K 
                  IF (A(J,L)==ZERO .AND. B(J,L)==ZERO) CYCLE  
                  TEMP1 = ALPHA*B(J,L) 
                  TEMP2 = ALPHA*A(J,L) 
                  C(:J,J) = C(:J,J) + A(:J,L)*TEMP1 + B(:J,L)*TEMP2 
               END DO 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (BETA == ZERO) THEN 
                  C(J:N,J) = ZERO 
               ELSE IF (BETA /= ONE) THEN 
                  C(J:N,J) = BETA*C(J:N,J) 
               ENDIF 
               DO L = 1, K 
                  IF (A(J,L)==ZERO .AND. B(J,L)==ZERO) CYCLE  
                  TEMP1 = ALPHA*B(J,L) 
                  TEMP2 = ALPHA*A(J,L) 
                  C(J:N,J) = C(J:N,J) + A(J:N,L)*TEMP1 + B(J:N,L)*TEMP2 
               END DO 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  C := alpha*A'*B + alpha*B'*A + C.
!
         IF (UPPER) THEN 
            DO J = 1, N 
               DO I = 1, J 
                  TEMP1 = ZERO 
                  TEMP2 = ZERO 
                  TEMP1 = TEMP1 + SUM(A(:K,I)*B(:K,J)) 
                  TEMP2 = TEMP2 + SUM(B(:K,I)*A(:K,J)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2 
                  ELSE 
                     C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + ALPHA*TEMP2 
                  ENDIF 
               END DO 
            END DO 
         ELSE 
            DO J = 1, N 
               DO I = J, N 
                  TEMP1 = ZERO 
                  TEMP2 = ZERO 
                  TEMP1 = TEMP1 + SUM(A(:K,I)*B(:K,J)) 
                  TEMP2 = TEMP2 + SUM(B(:K,I)*A(:K,J)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = ALPHA*TEMP1 + ALPHA*TEMP2 
                  ELSE 
                     C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 + ALPHA*TEMP2 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSYR2K.
!
      END SUBROUTINE DSYR2K 
