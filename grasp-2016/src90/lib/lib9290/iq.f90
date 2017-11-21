!***********************************************************************
!                                                                      *
      INTEGER FUNCTION IQ (ISUBSH, ICSF) 
!                                                                      *
!   IQ is the occupation of subshell ISUBSH in CSF  ICSF.              *
!                                                                      *
!   Call(s) to: [LIB92]: IUNPCK.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 30 Oct 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:38   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE ORB_C,       ONLY: IQA, NCF
      USE IOUNIT_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISUBSH 
      INTEGER              :: ICSF 
      INTEGER, PARAMETER :: NNNW = 214
!-----------------------------------------------
!
! IQA declared I*4 and COMMON/iounit/ added
! XHH 1997.02.12
!
! XHH added
!
!
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NNNW)) THEN
        IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCF)) THEN
!cff        IQ = IUNPCK (IQA(1,ICSF),ISUBSH)
!GG      IQ = IBITS(IQA((ISUBSH - 1)/4 + 1,ICSF),8*MOD(ISUBSH - 1,4),8) 
      iq = IQA(isubsh,icsf)      
        ELSE
           WRITE(istde,*) 'IQ: Argument ICSF is out of range.'
           STOP
        ENDIF
     ELSE
        WRITE(istde,*) 'IQ: Argument ISUBSH is out of range.'
        WRITE(istde,*) 'ISUBSH=',ISUBSH, '  NNNW=',NNNW
        STOP
     ENDIF
!
      RETURN  
      END FUNCTION IQ 
