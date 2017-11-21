!
!***********************************************************************
!                                                                      *
      MODULE orb_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE 
!...Created by Pacific-Sierra Research 77to90  4.3E  08:57:22  12/25/06  

!     Formerly orb10_C and orb10r_C
      INTEGER, PARAMETER :: NNNW = 214 
      CHARACTER(LEN=2), DIMENSION(NNNW) :: NH 
      CHARACTER(LEN=2), DIMENSION(NNNW) :: NHR 
!     Formerly orb1_C
      REAL(DOUBLE), DIMENSION(NNNW) :: E, GAMA 
!     Formerly orb2_C and orb2r_C
      INTEGER :: NCF, NW, NCFR, NWR 
      INTEGER(BYTE), DIMENSION(:,:), pointer :: IQA
      INTEGER, DIMENSION(:,:), pointer :: IQAR
!     Formerly orb4_C and orb4r_C
      INTEGER, DIMENSION(NNNW) :: NP, NAK, NPR, NAKR 
!     Formerly orb5_c
      INTEGER, DIMENSION(NNNW) :: NKL, NKJ 
      END MODULE orb_C 
