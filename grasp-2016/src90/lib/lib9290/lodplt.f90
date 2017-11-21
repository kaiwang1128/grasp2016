!***********************************************************************
!                                                                      *
      SUBROUTINE LODPLT 
!                                                                      *
!   This subroutine loads   COMMON/LIB92P/    with the values of the   *
!   plants substituted in the LIB92 subprograms.                       *
!                                                                      *
!   Written by Farid A Parpia               Last update: 09 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:30:54   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE LIB92P_C 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NNNP = 590 
      INTEGER, PARAMETER :: NNN1 = NNNP + 10 
      INTEGER, PARAMETER :: NNNW = 214 
      INTEGER, PARAMETER :: KEYORB = 121 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP
!-----------------------------------------------
!
!
      LPLANT = .TRUE. 
!
      NPLANT(1) = KEYORB 
      NPLANT(2) = NNNP 
!
!   Consistency check
!
      IF (NNN1 /= NNNP + 10) THEN 
         WRITE (6, *) 'LODPLT: Plant N1 should be equal to ', NP + 10, &
            ' in the LIB92 subprograms.' 
         STOP  
      ENDIF 
!
      NPLANT(3) = NNNW 
      NPLANT(4) = NNNW 
!
      RETURN  
      END SUBROUTINE LODPLT 
