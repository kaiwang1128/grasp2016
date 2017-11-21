************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK12(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN)
!      SUBROUTINE ONESCALAR11(IC,IR,NCOEC,ELEMNT)
!      SUBROUTINE RKCO_GG11(IC,IR,INCOR,NCTEC,INC2,NMCBP,NCORE,
!     :           ELSTO,ELEMNT)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 1 and 2                 *         
*                                                                      *
*   Written by Per Jönsson                                April 2017   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INCLUDE 'parameters.def'

      EXTERNAL BREID,CORD

      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))

      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      CHARACTER*2 NH
*
      INTEGER NTYPE(6,NCF)
      DIMENSION TSHELL(NNNW),EMTBLOCK(MAXSPAN,MAXSPAN)

      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6      
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ      
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)      
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-20)

      WRITE(*,*) 'In matrixblock12 IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0

*   Call onescalar        

      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
*   
*   Accumulate the contribution from the one-body operators:
*   kinetic energy, electron-nucleus interaction; update the
*   angular integral counter
*
*   For for type 1 and 2 there are only two possibilities.
*   IA = 0, the CSFs differ by two or more orbitals and no interaction
*   IA ~= 0 and IA ~= IB, the CSFs differ by one orbital       
*
      IF (IA .NE. 0) THEN

        NORB = NTYPE(4,IC) - NTYPE(3,IC) + 1 ! Number of symbolic orbs

*   Ensure that the indices are in `canonical' order

        IF (IA .GT. IB) THEN
          ISWAP = IB
          IB = IA
          IA = ISWAP
        ENDIF

        TCOEFF = DBLE(TSHELL(1))
        IF (DABS(TCOEFF) .GT. CUTOFF) THEN
          NCOEC = NCOEC + NTYPE(2,IC)
          DO I = 1,NORB
            IB = I + NTYPE(3,IC) - 1
            CALL IABINT(IA,IB,TEGRAL)
            EMTBLOCK(1,I) = EMTBLOCK(1,I) + TEGRAL*TCOEFF
            WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                   NP(IB),NH(IB),')'            
            IF (LNMS) THEN
              CALL KEINT(IA,IB,TEGRAL,TEGRAL)
              EMTBLOCK(1,I) = EMTBLOCK(1,I) + TEGRAL*ATWINV*TCOEFF
            ENDIF
            IF (LVP) THEN
              CALL VPINT(IA,IB,TEGRAL,TEGRAL)
              EMTBLOCK(1,I) = EMTBLOCK(1,I) + TEGRAL*TCOEFF
            ENDIF
          END DO
        ENDIF
      ENDIF

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
* Kai: I modify the code from this line.**********************************    
      DO 8 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC)
          DO 9 J = NORB, 1, -1
            IF(J .NE. NORB) THEN
              LABEL(4,I) = J + NTYPE(3,IC) - 1 
            END IF
            IF (LSMS) THEN
              IF (LABEL(5,I) .EQ. 1) THEN
                CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                EMTBLOCK(1,J) = EMTBLOCK(1,J) - TGRL1*TGRL2*ATWINV*VCOEFF
              ENDIF
            ENDIF
            CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
    
            WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                    NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                    NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                    NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                    NP(LABEL(4,I)),NH(LABEL(4,I)),')'
     
            EMTBLOCK(1,J) = EMTBLOCK(1,J) + TEGRAL*VCOEFF
    9     CONTINUE
        ENDIF
    8 CONTINUE
    
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

        DO 10 I = 1, NVCOEF
          CONTR = 0 ! Kai: For different symbolic orbital, 
                    !should I initialize the CONTR value?
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            NMCBP = NMCBP + NTYPE(2,IC)
            DO 11 J = NORB, 1, -1
              ITYPE = ABS (LABEL(6,I))
              IF(J .NE. NORB) THEN
                LABEL(4,I) = J + NTYPE(3,IC) - 1 
              END IF
              IF (ITYPE .EQ. 1) THEN
                CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 2) THEN
                CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 3) THEN
                CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 4) THEN
                CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 5) THEN
                CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 6) THEN
                CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ENDIF 
              CONTR = COEFF(I)*TEGRAL
              IF (LABEL(6,I) .GT. 0) THEN
                EMTBLOCK(1,J) = EMTBLOCK(1,J) + CONTR
              ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                NCORE = NCORE + NTYPE(2,IC)     
                ELSTO = ELSTO + CONTR  ! Kai: Should ELSTO be modified?
              ENDIF
   11     CONTINUE
        ENDIF
   10 CONTINUE
*
        IBUG1 = 0
* 
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
      ENDIF
!
  592 FORMAT (2I5,F15.9,1X,1A2,1I1,1A2,1A1,I1,1A2,1A1)
  593 FORMAT (2I5,F15.9,1X,1A,1I1,1A1,1I1,1A2,1I1,1A2,
     :        1A1,1I1,1A2,1I1,1A2,1A1) 

      RETURN
      END
