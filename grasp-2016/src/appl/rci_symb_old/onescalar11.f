************************************************************************
*                                                                      *
      SUBROUTINE ONESCALAR11(IC,IR,NCOEC,ELEMNT)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 1                       *         
*                                                                      *
*   Written by Per J                                      April 2017   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INCLUDE 'parameters.def'
      
      POINTER (PNTRIQ,RIQDUMMY)

      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION TSHELL(NNNW)

      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ      
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-20)

      ATWINV = 1.D0/EMN

      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)

*   Accumulate the contribution from the one-body operators:
*   kinetic energy, electron-nucleus interaction; update the
*   angular integral counter
*
      IF (IA .NE. 0) THEN
        IF (IA .EQ. IB) THEN
          DO IA = 1,NW
            TCOEFF = DBLE(TSHELL(IA))
            IF (DABS (TCOEFF) .GT. CUTOFF) THEN
              NCOEC = NCOEC + 1
              CALL IABINT (IA, IA, TEGRAL)
              ELEMNT = ELEMNT + TEGRAL*TCOEFF
              IF (LNMS) THEN
                CALL KEINT (IA,IA,TEGRAL)
                ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
              ENDIF
              IF (LVP) THEN
                CALL VPINT (IA, IA, TEGRAL)
                ELEMNT = ELEMNT + TEGRAL*TCOEFF
              ENDIF
            ENDIF
          ENDDO
        ELSE
          TCOEFF = DBLE(TSHELL(1))
          IF (DABS (TCOEFF) .GT. CUTOFF) THEN
            NCOEC = NCOEC + 1
            CALL IABINT (IA, IB, TEGRAL)
            ELEMNT = ELEMNT + TEGRAL*TCOEFF
            IF (LNMS) THEN
              CALL KEINT (IA, IB, TEGRAL)
              ELEMNT = ELEMNT + TEGRAL*ATWINV*TCOEFF
            ENDIF
            IF (LVP) THEN
              CALL VPINT (IA, IB, TEGRAL)
              ELEMNT = ELEMNT + TEGRAL*TCOEFF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*
      RETURN
      END
