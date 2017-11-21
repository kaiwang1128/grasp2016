 
!***********************************************************************
!                                                                      *
      SUBROUTINE CHKPLT 
!***********************************************************************
!                                                                      *
!   This code checks for  consistent substitutions of plants between   *
!   GENMCP and the LIB92 subprograms.
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:32:22   1/ 3/07  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE LIB92P_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lodplt_I 
      USE convrt_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEYORB = 121 
      INTEGER, PARAMETER :: NNNW = 214
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENTH, NWP 
      CHARACTER :: RECORD*20 
!
!
!
!   Load COMMON/LIB92P/
!
      CALL LODPLT 
!
!   Consistently DOUBLEPRECISION or REAL?
!
      IF (LPLANT .NEQV. .TRUE.) THEN 
         IF (LPLANT) THEN 
            WRITE (6, *) 'Plant DB was set to .TRUE. in LIB92,' 
            WRITE (6, *) ' but to .FALSE. in GENMCP.' 
         ELSE 
            WRITE (6, *) 'Plant DB was set to .FALSE. in LIB92,' 
            WRITE (6, *) ' but to .TRUE. in GENMCP.' 
         ENDIF 
         STOP  
      ENDIF 
!
!   Consistent numerical plants? Check only those relevant
!   to GENMCP
!
      IF (NPLANT(1) /= KEYORB) THEN 
         CALL CONVRT (NPLANT(1), RECORD, LENTH) 
         WRITE (6, *) 'CHKPLT: Plant KEYORB has been set to ' 
         WRITE (6, *) ' '//RECORD(1:LENTH)//' in LIB92, but to KEYORB' 
         WRITE (6, *) ' in GENMCP.' 
         STOP  
      ENDIF 
!
      IF (NPLANT(3) /= NNNW) THEN 
         CALL CONVRT (NPLANT(3), RECORD, LENTH) 
         WRITE (6, *) 'CHKPLT: Plant NW has been set to ' 
         WRITE (6, *) ' '//RECORD(1:LENTH)//' in LIB92, but to NNNW' 
         WRITE (6, *) ' in GENMCP.' 
         STOP  
      ENDIF 
!
      IF (NPLANT(4) /= NNNW) THEN 
         CALL CONVRT (NPLANT(4), RECORD, LENTH) 
         WRITE (6, *) 'CHKPLT: Plant NWP has been set to ' 
         WRITE (6, *) ' '//RECORD(1:LENTH)//' in LIB92, but to NNNW' 
         WRITE (6, *) ' in GENMCP.' 
         STOP  
      ENDIF 
!
!      IF (MOD(NNNW,4) == 0) THEN 
!         NWP = NNNW/4 
!      ELSE 
!         NWP = NNNW/4 + 1 
!      ENDIF 
!
!      IF (NNNW /= NWP) THEN 
!         CALL CONVRT (NWP, RECORD, LENTH) 
!         WRITE (6, *) 'CHKPLT: Plant NWP should be set' 
!         WRITE (6, *) ' to '//RECORD(1:LENTH)//'.' 
!         STOP  
!      ENDIF 

      RETURN  
      END SUBROUTINE CHKPLT 
