************************************************************************
*                                                                      *
      SUBROUTINE RKCO_GG11(IC,IR,INCOR,NCTEC,INC2,NMCBP,NCORE,
     :           ELSTO,ELEMNT)
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8 (A-H, O-Z)

      INCLUDE 'parameters.def'

      EXTERNAL BREID,CORD

      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))

      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS      

      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6      
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z 
*
      PARAMETER (CUTOFF = 1.0D-20)

      ATWINV = 1.D0/EMN

      IBUG1 = 0
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*
      NVCOEF = 0
*
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*
      DO 7 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + 1
          IF (LSMS) THEN
            IF (LABEL(5,I) .EQ. 1) THEN
              CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
              CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
              ELEMNT = ELEMNT - TGRL1*TGRL2*ATWINV*VCOEFF
            ENDIF
          ENDIF
          CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                 LABEL(3,I), LABEL(4,I),
     :                 LABEL(5,I), TEGRAL)
          ELEMNT = ELEMNT + TEGRAL*VCOEFF
        ENDIF
    7 CONTINUE
*
      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
        NVCOEF = 0
*
        CALL RKCO_GG (IC, IR, BREID, 1, 2)

        DO 8 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            IF (ITYPE .EQ. 1) THEN
              CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 2) THEN
              CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 3) THEN
              CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 4) THEN
              CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 5) THEN
              CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 6) THEN
              CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ENDIF 
            CONTR = COEFF(I)*TEGRAL
            IF (LABEL(6,I) .GT. 0) THEN
              ELEMNT = ELEMNT + CONTR
            ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
              NCORE = NCORE + 1
              ELSTO = ELSTO + CONTR
            ENDIF
          ENDIF
    8   CONTINUE
*
        IBUG1 = 0
* 
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
      ENDIF
!
      RETURN
      END
