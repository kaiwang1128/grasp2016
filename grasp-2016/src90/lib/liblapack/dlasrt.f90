      SUBROUTINE DLASRT(ID, N, D, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:33   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !!USE lsame_I 
      !!USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: INFO 
      CHARACTER  :: ID 
      REAL(DOUBLE) , INTENT(INOUT) :: D(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: SELECT = 20 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: DIR, ENDD, I, J, START, STKPNT 
      INTEGER, DIMENSION(2,32) :: STACK 
      REAL(DOUBLE) :: D1, D2, D3, DMNMX, TMP 
      logical::lsame
!-----------------------------------------------
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
      INFO = 0 
      DIR = -1 
      IF (LSAME(ID,'D')) THEN 
         DIR = 0 
      ELSE IF (LSAME(ID,'I')) THEN 
         DIR = 1 
      ENDIF 
      IF (DIR == (-1)) THEN 
         INFO = -1 
      ELSE IF (N < 0) THEN 
         INFO = -2 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DLASRT', (-INFO)) 
         RETURN  
      ENDIF 
!
!     Quick return if possible
!
      IF (N <= 1) RETURN  
!
      STKPNT = 1 
      STACK(1,1) = 1 
      STACK(2,1) = N 
   10 CONTINUE 
      START = STACK(1,STKPNT) 
      ENDD = STACK(2,STKPNT) 
      STKPNT = STKPNT - 1 
      IF (ENDD - START<=SELECT .AND. ENDD-START>0) THEN 
!
!        Do Insertion sort on D( START:ENDD )
!
         IF (DIR == 0) THEN 
!
!           Sort into decreasing order
!
            L30: DO I = START + 1, ENDD 
               DO J = I, START + 1, -1 
                  IF (D(J) > D(J-1)) THEN 
                     DMNMX = D(J) 
                     D(J) = D(J-1) 
                     D(J-1) = DMNMX 
                  ELSE 
                     CYCLE  L30 
                  ENDIF 
               END DO 
            END DO L30 
!
         ELSE 
!
!           Sort into increasing order
!
            L50: DO I = START + 1, ENDD 
               DO J = I, START + 1, -1 
                  IF (D(J) < D(J-1)) THEN 
                     DMNMX = D(J) 
                     D(J) = D(J-1) 
                     D(J-1) = DMNMX 
                  ELSE 
                     CYCLE  L50 
                  ENDIF 
               END DO 
            END DO L50 
!
         ENDIF 
!
      ELSE IF (ENDD - START > SELECT) THEN 
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D(START) 
         D2 = D(ENDD) 
         I = (START + ENDD)/2 
         D3 = D(I) 
         IF (D1 < D2) THEN 
            IF (D3 < D1) THEN 
               DMNMX = D1 
            ELSE IF (D3 < D2) THEN 
               DMNMX = D3 
            ELSE 
               DMNMX = D2 
            ENDIF 
         ELSE 
            IF (D3 < D2) THEN 
               DMNMX = D2 
            ELSE IF (D3 < D1) THEN 
               DMNMX = D3 
            ELSE 
               DMNMX = D1 
            ENDIF 
         ENDIF 
!
         IF (DIR == 0) THEN 
!
!           Sort into decreasing order
!
            I = START - 1 
            J = ENDD + 1 
   60       CONTINUE 
            J = J - 1 
            DO WHILE(D(J) < DMNMX) 
               J = J - 1 
            END DO 
            I = I + 1 
            DO WHILE(D(I) > DMNMX) 
               I = I + 1 
            END DO 
            IF (I < J) THEN 
               TMP = D(I) 
               D(I) = D(J) 
               D(J) = TMP 
               GO TO 60 
            ENDIF 
            IF (J - START > ENDD - J - 1) THEN 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = START 
               STACK(2,STKPNT) = J 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = J + 1 
               STACK(2,STKPNT) = ENDD 
            ELSE 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = J + 1 
               STACK(2,STKPNT) = ENDD 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = START 
               STACK(2,STKPNT) = J 
            ENDIF 
         ELSE 
!
!           Sort into increasing order
!
            I = START - 1 
            J = ENDD + 1 
   90       CONTINUE 
            J = J - 1 
            DO WHILE(D(J) > DMNMX) 
               J = J - 1 
            END DO 
            I = I + 1 
            DO WHILE(D(I) < DMNMX) 
               I = I + 1 
            END DO 
            IF (I < J) THEN 
               TMP = D(I) 
               D(I) = D(J) 
               D(J) = TMP 
               GO TO 90 
            ENDIF 
            IF (J - START > ENDD - J - 1) THEN 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = START 
               STACK(2,STKPNT) = J 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = J + 1 
               STACK(2,STKPNT) = ENDD 
            ELSE 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = J + 1 
               STACK(2,STKPNT) = ENDD 
               STKPNT = STKPNT + 1 
               STACK(1,STKPNT) = START 
               STACK(2,STKPNT) = J 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (STKPNT > 0) GO TO 10 
      RETURN  
!
!     End of DLASRT
!
      END SUBROUTINE DLASRT 
