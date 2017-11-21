      SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C&
         , LDC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:28   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !USExerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: LDB 
      INTEGER , INTENT(IN) :: LDC 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: TRANSA 
      CHARACTER  :: TRANSB 
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
      INTEGER :: I, INFO, J, L, NCOLA, NROWA, NROWB 
      REAL(DOUBLE) :: TEMP 
      LOGICAL :: NOTA, NOTB 
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
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
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
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
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
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME(TRANSA,'N') 
      NOTB = LSAME(TRANSB,'N') 
      IF (NOTA) THEN 
         NROWA = M 
         NCOLA = K 
      ELSE 
         NROWA = K 
         NCOLA = M 
      ENDIF 
      IF (NOTB) THEN 
         NROWB = K 
      ELSE 
         NROWB = N 
      ENDIF 
!
!     Test the input parameters.
!
      INFO = 0 
      IF (.NOT.NOTA .AND. .NOT.LSAME(TRANSA,'C') .AND. .NOT.LSAME(TRANSA,'T')) &
         THEN 
         INFO = 1 
      ELSE IF (.NOT.NOTB .AND. .NOT.LSAME(TRANSB,'C') .AND. .NOT.LSAME(TRANSB,&
            'T')) THEN 
         INFO = 2 
      ELSE IF (M < 0) THEN 
         INFO = 3 
      ELSE IF (N < 0) THEN 
         INFO = 4 
      ELSE IF (K < 0) THEN 
         INFO = 5 
      ELSE IF (LDA < MAX(1,NROWA)) THEN 
         INFO = 8 
      ELSE IF (LDB < MAX(1,NROWB)) THEN 
         INFO = 10 
      ELSE IF (LDC < MAX(1,M)) THEN 
         INFO = 13 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DGEMM ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (M==0 .OR. N==0 .OR. (ALPHA==ZERO .OR. K==0) .AND. BETA==ONE) RETURN  
!
!     And if  alpha.eq.zero.
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
      IF (NOTB) THEN 
         IF (NOTA) THEN 
!
!           Form  C := alpha*A*B + beta*C.
!
            DO J = 1, N 
               IF (BETA == ZERO) THEN 
                  C(:M,J) = ZERO 
               ELSE IF (BETA /= ONE) THEN 
                  C(:M,J) = BETA*C(:M,J) 
               ENDIF 
               DO L = 1, K 
                  IF (B(L,J) == ZERO) CYCLE  
                  TEMP = ALPHA*B(L,J) 
                  C(:M,J) = C(:M,J) + TEMP*A(:M,L) 
               END DO 
            END DO 
         ELSE 
!
!           Form  C := alpha*A'*B + beta*C
!
            DO J = 1, N 
               DO I = 1, M 
                  TEMP = ZERO 
                  TEMP = TEMP + SUM(A(:K,I)*B(:K,J)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = ALPHA*TEMP 
                  ELSE 
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J) 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      ELSE 
         IF (NOTA) THEN 
!
!           Form  C := alpha*A*B' + beta*C
!
            DO J = 1, N 
               IF (BETA == ZERO) THEN 
                  C(:M,J) = ZERO 
               ELSE IF (BETA /= ONE) THEN 
                  C(:M,J) = BETA*C(:M,J) 
               ENDIF 
               DO L = 1, K 
                  IF (B(J,L) == ZERO) CYCLE  
                  TEMP = ALPHA*B(J,L) 
                  C(:M,J) = C(:M,J) + TEMP*A(:M,L) 
               END DO 
            END DO 
         ELSE 
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO J = 1, N 
               DO I = 1, M 
                  TEMP = ZERO 
                  TEMP = TEMP + SUM(A(:K,I)*B(J,:K)) 
                  IF (BETA == ZERO) THEN 
                     C(I,J) = ALPHA*TEMP 
                  ELSE 
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J) 
                  ENDIF 
               END DO 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DGEMM .
!
      END SUBROUTINE DGEMM 
