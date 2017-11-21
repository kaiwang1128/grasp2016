************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK15(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 1 and 4                 *         
*                                                                      *
*   Written by Per Jönsson & Kai Wang                   September 2017 *
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

      WRITE(*,*) 'In matrixblock15, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0

*   Call onescalar to       
*   There are at least two different electrons 
*   No contribution from the one electron operator

*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*
      NVCOEF = 0

*
*   The line (NORBL = NTYPE(6,IC) - NTYPE(5,IC)) should be replaced by 
*   the line (NORBL = NTYPE(6,IC) - NTYPE(5,IC) + 1), when the value NTYPE(5,IC)
*   is corrected by Per. 
*   The lines (                  LABEL(3,I) = J + NTYPE(3,IC) - 1
!                  LABEL(4,I) = K + NTYPE(5,IC) - 1
*   ) and the lines (!                  LABEL(3,I) = K + NTYPE(5,IC) - 1
!                  LABEL(3,I) = K + NTYPE(5,IC)
*   ) should also be modified.
*
      NORBU = NTYPE(4,IC) - NTYPE(3,IC) + 1 ! Number of symbolic orbs
!      NORBL = NTYPE(6,IC) - NTYPE(5,IC) ! Number of symbolic orbs
!      NORBL = NTYPE(6,IC) - NTYPE(5,IC) + 1 ! Number of symbolic orbs

!      write(*,*)NORBU, NTYPE(4,IC),NTYPE(3,IC),NTYPE(4,IC) - NTYPE(3,IC) + 1 
!      write(*,*)NORBL, NTYPE(6,IC),NTYPE(5,IC),NTYPE(6,IC) - NTYPE(5,IC) + 1 
!      write(*,*)(NTYPE(I,IC),I=1,6) 
      
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      
*     There are two cases when IC and IR are of type 1 and 4.
*     The two symbolic orbits of IC of type4 are also at the places of LABEL(3,I) and LABEL(4,I).
*     Just chage the values of LABEL(3,I) and LABEL(4,I).
*     Case I:MRANK .EQ. 1, LABEL(3,I) = J + NTYPE(3,IC) - 1 and LABEL(4,I) = K + NTYPE(5,IC) - 1
*     Case 2:MRANK .EQ. 2, LABEL(4,I) = J + NTYPE(3,IC) - 1 and LABEL(3,I) = K + NTYPE(5,IC) - 1
 
!      DO  I = 1, NVCOEF
!        VCOEFF = COEFF(I)
!            CALL RKINTC (LABEL(1,I), LABEL(2,I),
!     :                   LABEL(3,I), LABEL(4,I),
!     :                   LABEL(5,I), TEGRAL)
!
!           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
!     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!      END DO
!      write(*,*)

      DO 11 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC)
          M = NTYPE(2,IC) + 1 
          MRANK = 0 
          DO 12 J = NORBU, 1, -1   
              M = M - 1 
              IF (J .EQ. NORBU .AND. K .EQ. NORBL) THEN
                              
              ELSE
                IF (MRANK .EQ. 1) THEN
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                ELSE
                  WRITE(*,*)"STOP 1 IN MATRIXBLOCK15"
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
              IF(M .EQ. NTYPE(2,IC))then
                IF(NP(LABEL(3,I)) .EQ. NP(NTYPE(4,IC)) .AND.
     :             NP(LABEL(4,I)) .EQ. NP(NTYPE(4,IC)))THEN
                  MRANK = 1
                ELSE
                  WRITE(*,*)'STOP IN 2 MATRIXBLOCK15.f'
                  STOP
                END IF
              END IF
   12     CONTINUE
        END IF
   11 CONTINUE

        DO J = NTYPE(2,IC), 1, -1 
           write(*,*)IC,IR,1,J,EMTBLOCK(1,J)
        END DO
!      goto 600
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
          M = NTYPE(2,IC) + 1  
          DO 15 J = NORBU, 1, -1         
              M = M - 1 
              IF (J .EQ. NORBU .AND. K .EQ. NORBL) THEN
                 
              ELSE
               SELECT CASE (MRANK)
                CASE (1)
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                CASE (2)
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                CASE (3)
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                CASE (4)
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                CASE (5)
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                CASE (6)
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                CASE (7)
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                CASE (8)
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                CASE DEFAULT
                  WRITE(*,*)'STOP 3 IN MATRIXBLOCK15'                  
                  STOP                                        
               END SELECT 
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
              IF(M .EQ. NTYPE(2,IC))then
                IF(LABEL(1,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 1
                ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC).AND. 
     :             LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 2
                ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 3
                ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 4
                ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 5
                ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 6
                ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 7
                ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC) .AND. 
     :             LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                  MRANK = 8
                ELSE
                  WRITE(*,*)'STOP IN 4 MATRIXBLOCK15'
                  STOP
                END IF
              END IF
   15     CONTINUE
        END IF
   14 CONTINUE
   
        DO J = NTYPE(2,IC), 1, -1 
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

  600 continue
      RETURN
      END
