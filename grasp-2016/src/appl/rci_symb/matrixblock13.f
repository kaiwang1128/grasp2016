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
*   Written by Per Jönsson & Kai Wang                         May 2017 *
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
      WRITE(*,*)"Here comes the off-diagonal part within the diagonal 
     : symbolic block in matrixblock13"

!         NVCOEF = 0
!         ELEMNT = 0
!         CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*
!         DO 7 I = 1, NVCOEF
!            VCOEFF = COEFF(I)
!            IF (DABS (VCOEFF) .GT. CUTOFF) THEN
!               NCTEC = NCTEC + 1
!               IF (LSMS) THEN
!                  IF (LABEL(5,I) .EQ. 1) THEN
!                     CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
!                     CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
!                     ELEMNT = ELEMNT - TGRL1*TGRL2*ATWINV*VCOEFF
!                  ENDIF
!               ENDIF
!               CALL RKINTC (LABEL(1,I), LABEL(2,I),
!     :                      LABEL(3,I), LABEL(4,I),
!     :                      LABEL(5,I), TEGRAL)

!               WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
!     :                    NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                    NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                    NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                    NP(LABEL(4,I)),NH(LABEL(4,I)),')'

!               ELEMNT = ELEMNT + TEGRAL*VCOEFF
!            ENDIF
!    7    CONTINUE

!           write(*,*)IC,IR,ELEMNT

*

      NVCOEF = 0
*
      NORBU = NTYPE(4,IC) - NTYPE(3,IC) + 1 ! Number of symbolic orbs
      NORBL = NTYPE(6,IC) - NTYPE(5,IC) + 1 ! Number of symbolic orbs
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)

      DO 11 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC)
          M = NORBU * NORBL + 1 
          MRANK = 0 
          DO 12 J = NORBL, 1, -1   
            DO 13 K = NORBU, 1, -1     
              M = M - 1 
              IF (J .EQ. NORBL .AND. K .EQ. NORBU) THEN
                              
              ELSE
                IF (MRANK .EQ. 1) THEN
                  LABEL(3,I) = K + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 2) THEN
                  LABEL(4,I) = K + NTYPE(3,IC) - 1
                  LABEL(3,I) = J + NTYPE(5,IC) - 1
                ELSE
                  WRITE(*,*)"STOP 1 IN MATRIXBLOCK13"
                  STOP
                END IF                  
              END IF
              IF (LSMS) THEN
                IF (LABEL(5,I) .EQ. 1) THEN
                  CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                  CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                  EMTBLOCK(1,M) = EMTBLOCK(1,M) - TGRL1*TGRL2*ATWINV*VCOEFF
                END IF
              END IF
              CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
              EMTBLOCK(1,M) = EMTBLOCK(1,M) + TEGRAL*VCOEFF
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
              IF(M .EQ. NORBU * NORBL)then
                IF(NP(LABEL(3,I)) .EQ. NP(NTYPE(4,IC)))THEN
                  MRANK = 1
                ELSE IF(NP(LABEL(4,I)) .EQ. NP(NTYPE(4,IC)))THEN
                  MRANK = 2
                ELSE
                  WRITE(*,*)'STOP IN 2 MATRIXBLOCK13.f'
                  STOP
                END IF
              END IF
   13       CONTINUE          
   12     CONTINUE
        END IF
   11 CONTINUE

        DO J = NORBU * NORBL, 1, -1 
           write(*,*)IC,IR,1,J,EMTBLOCK(1,J)
        END DO

      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN !Kai:What do these arugments mean?
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
        NVCOEF = 0
*       

      CALL RKCO_GG (IC, IR, BREID, 1, 2)

      DO 14 I = 1, NVCOEF
        IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
          MRANK = 0
          NMCBP = NMCBP + NTYPE(2,IC)
          ITYPE = ABS (LABEL(6,I))
          M = NORBU * NORBL + 1 
          DO 15 J = NORBL, 1, -1         !7d    NTYPE(5,IC)
            DO 16 K = NORBU, 1, -1       !9p    NTYPE(3,IC)
              M = M - 1 
              IF (J .EQ. NORBL .AND. K .EQ. NORBU) THEN
                 
              ELSE
                IF (MRANK .EQ. 1) THEN
                  LABEL(1,I) = J + NTYPE(5,IC) - 1
                  LABEL(2,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 2) THEN
                  LABEL(1,I) = J + NTYPE(5,IC) - 1
                  LABEL(3,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 3) THEN
                  LABEL(1,I) = J + NTYPE(5,IC) - 1
                  LABEL(4,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 4) THEN
                  LABEL(2,I) = J + NTYPE(5,IC) - 1
                  LABEL(1,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 5) THEN
                  LABEL(2,I) = J + NTYPE(5,IC) - 1
                  LABEL(3,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 6) THEN
                  LABEL(2,I) = J + NTYPE(5,IC) - 1
                  LABEL(4,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 7) THEN
                  LABEL(3,I) = J + NTYPE(5,IC) - 1
                  LABEL(1,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 8) THEN
                  LABEL(3,I) = J + NTYPE(5,IC) - 1
                  LABEL(2,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 9) THEN
                  LABEL(3,I) = J + NTYPE(5,IC) - 1
                  LABEL(4,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 10) THEN
                  LABEL(4,I) = J + NTYPE(5,IC) - 1
                  LABEL(1,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 11) THEN
                  LABEL(4,I) = J + NTYPE(5,IC) - 1
                  LABEL(2,I) = K + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 12) THEN
                  LABEL(4,I) = J + NTYPE(5,IC) - 1
                  LABEL(3,I) = K + NTYPE(3,IC) - 1
                ELSE
                   WRITE(*,*)'STOP IN 3 MATRIXBLOCK13.f'
                  STOP
                END IF                  
              ENDIF
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'2222'

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
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
              CONTR = COEFF(I) * TEGRAL
              IF (LABEL(6,I) .GT. 0) THEN
                EMTBLOCK(1,M) = EMTBLOCK(1,M) + CONTR
              ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                NCORE = NCORE + NTYPE(2,IC)
                ELSTO = ELSTO + CONTR
              END IF
              IF(M .EQ. NORBU * NORBL)then
                IF(LABEL(1,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 1
                ELSE IF(LABEL(1,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 2
                ELSE IF(LABEL(1,I) .EQ. NTYPE(6,IC).AND. 
     :             LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 3
                ELSE IF(LABEL(2,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 4
                ELSE IF(LABEL(2,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 5
                ELSE IF(LABEL(2,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 6
                ELSE IF(LABEL(3,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 7
                ELSE IF(LABEL(3,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 8
                ELSE IF(LABEL(3,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 9
                ELSE IF(LABEL(4,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 10
                ELSE IF(LABEL(4,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 11
                ELSE IF(LABEL(4,I) .EQ. NTYPE(6,IC) .AND. 
     :             LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 12
                ELSE
                  WRITE(*,*)'STOP IN 2 MATRIXBLOCK13.f'
                  STOP
                END IF
              END IF
   16       CONTINUE          
   15     CONTINUE
        END IF
   14 CONTINUE
   
        DO J = NORBU * NORBL, 1, -1 
           write(*,*)IC,IR,1,J,EMTBLOCK(1,J)
        END DO
   
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
  594 FORMAT (5I5,F15.9,1X,1A,1I1,1A1,1I2,1A2,1I2,1A2,
     :        1A1,1I2,1A2,1I2,1A2,1A1) 

      RETURN
      END
