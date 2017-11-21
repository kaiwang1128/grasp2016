      SUBROUTINE DLARF(SIDE, M, N, V, INCV, TAU, C, LDC, WORK) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:32   2/12/04  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE dgemv_I 
      !USE dger_I 
      !!USE lsame_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: M 
      INTEGER  :: N 
      INTEGER  :: INCV 
      INTEGER  :: LDC 
      REAL(DOUBLE) , INTENT(IN) :: TAU 
      CHARACTER  :: SIDE 
      REAL(DOUBLE)  :: V(*) 
      REAL(DOUBLE)  :: C(LDC,*) 
      REAL(DOUBLE)  :: WORK(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
      logical::lsame
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      IF (LSAME(SIDE,'L')) THEN 
!
!        Form  H * C
!
         IF (TAU /= ZERO) THEN 
!
!           w := C' * v
!
            CALL DGEMV ('Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1) 
!
!           C := C - v * w'
!
            CALL DGER (M, N, (-TAU), V, INCV, WORK, 1, C, LDC) 
         ENDIF 
      ELSE 
!
!        Form  C * H
!
         IF (TAU /= ZERO) THEN 
!
!           w := C * v
!
            CALL DGEMV ('No transpose', M, N, ONE, C, LDC, V, INCV, ZERO, WORK&
               , 1) 
!
!           C := C - w * v'
!
            CALL DGER (M, N, (-TAU), WORK, 1, V, INCV, C, LDC) 
         ENDIF 
      ENDIF 
      RETURN  
!
!     End of DLARF
!
      END SUBROUTINE DLARF 
