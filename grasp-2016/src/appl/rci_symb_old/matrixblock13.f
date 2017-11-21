************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK13(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
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

      WRITE(*,*) 'In matrixblock13, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0

*   Call onescalar to       


*   No contribution from the one electron operator

*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*
    
*   Here comes the off-diagonal part within the diagonal symbolic block.
      NVCOEF = 0
*
      NORBU = NTYPE(4,IC) - NTYPE(3,IC) + 1 ! Number of symbolic orbs
      NORBL = NTYPE(6,IC) - NTYPE(5,IC) + 1 ! Number of symbolic orbs

      WRITE(*,*)"Here comes the off-diagonal part within the diagonal 
     : symbolic block in matrixblock13"
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      MRANK = 0 
      DO 11 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          MRANK = MRANK + 1
          NCTEC = NCTEC + NTYPE(2,IC)
          M = NORBU * NORBL  
          DO 12 J = NORBL, 1, -1
            DO 13 K = NORBU, 1, -1      
              IF (J .LE. K) THEN
                VCOEFF = COEFF(I)
              ELSE
                VCOEFF = -COEFF(I)
              END IF
              M = M - 1 
              IF (J .EQ. NORBL .AND. K .EQ. NORBU) THEN
                 
              ELSE
                IF (MRANK .EQ. 1) THEN
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                  LABEL(4,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 2) THEN
                  LABEL(3,I) = K + NTYPE(5,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                ELSE
                  STOP
                END IF                  
              ENDIF
              IF (LSMS) THEN
                IF (LABEL(5,I) .EQ. 1) THEN
                  CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                  CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                  EMTBLOCK(M,1) = EMTBLOCK(M,1) - TGRL1*TGRL2*ATWINV*VCOEFF
                                               ! negative contributions
                ENDIF
              ENDIF
              CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
              EMTBLOCK(M,1) = EMTBLOCK(M,1) + TEGRAL*VCOEFF
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
   13       CONTINUE          
   12     CONTINUE
        ENDIF
   11 CONTINUE
      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN !Kai:What do these arugments mean?
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
        NVCOEF = 0
*       

      CALL RKCO_GG (IC, IR, BREID, 1, 2)
      MRANK = 0 

      DO 14 I = 1, NVCOEF
        IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
          NMCBP = NMCBP + NTYPE(2,IC)
          ITYPE = ABS (LABEL(6,I))
          M = NORBU * NORBL  
          DO 15 J = NORBL, 1, -1
            DO 16 K = NORBU, 1, -1      
              IF (J .LE. K) THEN
                VCOEFF = COEFF(I)
              ELSE
                VCOEFF = - COEFF(I)
              END IF
              M = M - 1 
              IF (J .EQ. NORBL .AND. K .EQ. NORBU) THEN
                 
              ELSE
                IF (MRANK .EQ. 1) THEN
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                  LABEL(4,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 2) THEN
                  LABEL(3,I) = K + NTYPE(5,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                ELSE
                  STOP
                END IF                  
              ENDIF
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
              CONTR = VCOEFF * TEGRAL
              IF (LABEL(6,I) .GT. 0) THEN
                EMTBLOCK(M,1) = EMTBLOCK(M,1) + CONTR
              ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
              NCORE = NCORE + NTYPE(2,IC)
              ELSTO = ELSTO + CONTR
            ENDIF
   16       CONTINUE          
   15     CONTINUE

          ENDIF
   14   CONTINUE
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
