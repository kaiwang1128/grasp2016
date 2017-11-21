!
!***********************************************************************
!                                                                      *
      MODULE prnt_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
      INTEGER :: NVEC, NVECMX 
      INTEGER :: NVECFF, NVECMXFF 
      INTEGER :: NVECII, NVECMXII 
!     REAL(DOUBLE) :: PNIVEC 
!     REAL(DOUBLE) :: PNIVECFF 
      REAL(DOUBLE) :: PNIVECII 
      INTEGER, DIMENSION(:), pointer :: ivec, ivecff, ivecii
      END MODULE prnt_C 
