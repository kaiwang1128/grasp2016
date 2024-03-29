************************************************************************
*                                                                      *
      FUNCTION JQSR (IWHICH,ISUBSH,ICSF)
*                                                                      *
*   JQSR is a subshell quantum number for subshell ISUBSH in configu-  *
*   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
*   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
! Packed array declared I*4, common/iounit/ added
! XHH 1997.02.12
Cww      INTEGER PJCUPR,PNTIQR
      POINTER (PJCUPR,JCUPRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      INTEGER*4 JQSAR
      POINTER (PNJQSR,JQSAR(NNNWP,3,1))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
      COMMON/iounit/istdi,istdo,istde
*
      IF ((IWHICH .GE. 1) .AND. (IWHICH .LE. 3)) THEN
         IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
            IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
               JQSR = IUNPCK (JQSAR(1,IWHICH,ICSF),ISUBSH)
            ELSE
               WRITE(istde,*) 'JQSR: Argument ICSF is out of range.'
               STOP
            ENDIF
         ELSE
            WRITE(istde,*) 'JQSR: Argument ISUBSH is out of range.'
            STOP
         ENDIF
      ELSE
         WRITE(istde,*) 'JQSR: Argument IWHICH is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
