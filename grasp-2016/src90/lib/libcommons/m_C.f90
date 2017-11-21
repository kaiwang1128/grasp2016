!
!***********************************************************************
!                                                                      *
      MODULE m_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06  
      INTEGER, PARAMETER :: NNNW = 214 
      INTEGER, DIMENSION(NNNW) :: NQ1, NQ2 
      INTEGER, DIMENSION(NNNW) :: JJC1, JJC2
      INTEGER, DIMENSION(3,NNNW) :: JJQ1, JJQ2 
      INTEGER, DIMENSION(NNNW) :: JLIST, KLIST 
      INTEGER :: NPEEL, NCORE 
      END MODULE m_C 
