      SUBROUTINE DSYMM(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:37   2/12/04  
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
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: LDB 
      INTEGER , INTENT(IN) :: LDC 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: SIDE 
      CHARACTER  :: UPLO 
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
      INTEGER :: I, INFO, J, K, NROWA 
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
!  DSYMM  performs one of the matrix-matrix operations
!
!     C := alpha*A*B + beta*C,
!
!  or
!
!     C := alpha*B*A + beta*C,
!
!  where alpha and beta are scalars,  A is a symmetric matrix and  B and
!  C are  m by n matrices.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE  specifies whether  the  symmetric matrix  A
!           appears on the  left or right  in the  operation as follows:
!
!              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
!
!              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of  the  symmetric  matrix   A  is  to  be
!           referenced as follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of the
!                                  symmetric matrix is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of the
!                                  symmetric matrix is to be referenced.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies the number of rows of the matrix  C.
!           M  must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix C.
!           N  must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
!           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading m by m upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  m by m  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced.
!           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
!           the array  A  must contain the  symmetric matrix,  such that
!           when  UPLO = 'U' or 'u', the leading n by n upper triangular
!           part of the array  A  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!           the leading  n by n  lower triangular part  of the  array  A
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  A  is not
!           referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry, the leading  m by n part of the array  B  must
!           contain the matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n updated
!           matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
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
!     Set NROWA as the number of rows of A.
!
      IF (LSAME(SIDE,'L')) THEN 
         NROWA = M 
      ELSE 
         NROWA = N 
      ENDIF 
      UPPER = LSAME(UPLO,'U') 
!
!     Test the input parameters.
!
      INFO = 0 
      IF (.NOT.LSAME(SIDE,'L') .AND. .NOT.LSAME(SIDE,'R')) THEN 
         INFO = 1 
      ELSE IF (.NOT.UPPER .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = 2 
      ELSE IF (M < 0) THEN 
         INFO = 3 
      ELSE IF (N < 0) THEN 
         INFO = 4 
      ELSE IF (LDA < MAX(1,NROWA)) THEN 
         INFO = 7 
      ELSE IF (LDB < MAX(1,M)) THEN 
         INFO = 9 
      ELSE IF (LDC < MAX(1,M)) THEN 
         INFO = 12 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSYMM ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (M==0 .OR. N==0 .OR. ALPHA==ZERO .AND. BETA==ONE) RETURN  
!
!     And when  alpha.eq.zero.
!
      IF (ALPHA == ZERO) THEN 
         IF (BETA == ZERO) THEN 
            C(:M,:N) = ZERO 
         ELSE 
            C(:M,:N) = BETA*C(:M,:N) 
         ENDIF 
         RETURN  
      ENDIF 
!
!     Start the operations.
!
      IF (LSAME(SIDE,'L')) THEN 
!
!        Form  C := alpha*A*B + beta*C.
!
         IF (UPPER) THEN 
            DO J = 1, N 
               DO I = 1, M 
                  TEMP1 = ALPHA*B(I,J) 
                  TEMP2 = ZERO 
                  C(:I-1,J) = C(:I-1,J) + TEMP1*A(:I-1,I) 
                  TEMP2 = TEMP2 + SUM(B(:I-1,J)*A(:I-1,I)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2 
                  ELSE 
                     C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2 
                  ENDIF 
               END DO 
            END DO 
         ELSE 
            DO J = 1, N 
               DO I = M, 1, -1 
                  TEMP1 = ALPHA*B(I,J) 
                  TEMP2 = ZERO 
                  C(I+1:M,J) = C(I+1:M,J) + TEMP1*A(I+1:M,I) 
                  TEMP2 = TEMP2 + SUM(B(I+1:M,J)*A(I+1:M,I)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = TEMP1*A(I,I) + ALPHA*TEMP2 
                  ELSE 
                     C(I,J) = BETA*C(I,J) + TEMP1*A(I,I) + ALPHA*TEMP2 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  C := alpha*B*A + beta*C.
!
         DO J = 1, N 
            TEMP1 = ALPHA*A(J,J) 
            IF (BETA == ZERO) THEN 
               C(:M,J) = TEMP1*B(:M,J) 
            ELSE 
               C(:M,J) = BETA*C(:M,J) + TEMP1*B(:M,J) 
            ENDIF 
            DO K = 1, J - 1 
               IF (UPPER) THEN 
                  TEMP1 = ALPHA*A(K,J) 
               ELSE 
                  TEMP1 = ALPHA*A(J,K) 
               ENDIF 
               C(:M,J) = C(:M,J) + TEMP1*B(:M,K) 
            END DO 
            DO K = J + 1, N 
               IF (UPPER) THEN 
                  TEMP1 = ALPHA*A(J,K) 
               ELSE 
                  TEMP1 = ALPHA*A(K,J) 
               ENDIF 
               C(:M,J) = C(:M,J) + TEMP1*B(:M,K) 
            END DO 
         END DO 
      ENDIF 
!
      RETURN  
!
!     End of DSYMM .
!
      END SUBROUTINE DSYMM 
