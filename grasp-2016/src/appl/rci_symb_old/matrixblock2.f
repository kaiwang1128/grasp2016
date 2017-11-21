************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK2(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
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

      WRITE(*,*) 'In matrixblock2, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0

*   Call onescalar to       

      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)

*   Diagonal case IC = IR. Start with computing the contribution
*   ENONSYM that does not include the symbolic orbital. For
*   NTYPE = 2 there is only one symbolic orbital 
*
      ENONSYM = 0.D0

      DO IA = 1,NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
          IF (IA.GT.NORBGEN) THEN ! IA > NORBGEB --> IA symbolic orb
            TCOEFFSAVE = TCOEFF
            IASAVE = IA
          ELSE
            CALL IABINT(IA,IA,TEGRAL)
            ENONSYM = ENONSYM + TEGRAL*TCOEFF
            IF (LNMS) THEN
              CALL KEINT(IA,IA,TEGRAL)
              ENONSYM = ENONSYM + TEGRAL*ATWINV*TCOEFF
            ENDIF
            IF (LVP) THEN
              CALL VPINT(IA,IA,TEGRAL)
              ENONSYM = ENONSYM + TEGRAL*TCOEFF
            END IF
          END IF
        END IF
      ENDDO
* Kai: I modify the code from this line.**********************************    
*  To get the diagonal matrix elements, add contributions from each of
*  the symbolic orbitals to ENONSYM 
      NORB = NTYPE(4,IC) - NTYPE(3,IC) + 1 ! Number of symbolic orbs
      NCOEC = NCOEC + NORB
      DO I = 1,NORB
        IA = I + NTYPE(3,IC) - 1
        CALL IABINT(IA,IA,TEGRAL)
        WRITE(*,592) IC,IC,TCOEFFSAVE,'I(',NP(IA),NH(IA),',',
     :                   NP(IA),NH(IA),')'  
        EMTBLOCK(I,I) = ENONSYM + TCOEFFSAVE*TEGRAL
        IF (LNMS) THEN
          CALL KEINT(IA,IA,TEGRAL)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + TEGRAL*ATWINV*TCOEFF
        ENDIF                   
        IF (LVP) THEN
          CALL VPINT(IA,IA,TEGRAL)
          EMTBLOCK(I,I) = EMTBLOCK(I,I) + TEGRAL*TCOEFF
        END IF
        WRITE(*,*) 'IC,I',IC,I,EMTBLOCK(I,I)
      END DO

*   Here comes the off-diagonal part within the diagonal symbolic block.
      CALL ONESCALAR(IC,IC-1,IA,IB,TSHELL) 
 
      WRITE(*,592) IC,IC-1,DBLE(TSHELL(1)),'I(',NP(IA),NH(IA),',',
     :                   NP(IB),NH(IB),')'        

*   Ensure that the indices are in `canonical' order

      IF (IA .GT. IB) THEN
        ISWAP = IB
        IB = IA
        IA = ISWAP
      ENDIF

      TCOEFF = DBLE(TSHELL(1))
      IF (DABS(TCOEFF) .GT. CUTOFF) THEN
        NCOEC = NCOEC + NTYPE(2,IC)
        DO I = 1,NORB - 1
           DO J = I + 1,NORB 
             IA = I + NTYPE(3,IC) - 1
             IB = J + NTYPE(3,IC) - 1
             CALL IABINT(IA,IB,TEGRAL)
             EMTBLOCK(I,J) = EMTBLOCK(I,J) + TEGRAL*TCOEFF
             WRITE(*,592) IC,IC-1,TCOEFF,'I(',NP(IA),NH(IA),',',
     :          NP(IB),NH(IB),')'
             IF (LNMS) THEN
               CALL KEINT(IA,IB,TEGRAL,TEGRAL)
               EMTBLOCK(I,J) = EMTBLOCK(I,J) + TEGRAL*ATWINV*TCOEFF
             ENDIF
             IF (LVP) THEN
               CALL VPINT(IA,IB,TEGRAL,TEGRAL)
               EMTBLOCK(I,J) = EMTBLOCK(I,J) + TEGRAL*TCOEFF
             ENDIF
          END DO   
        END DO
      ENDIF

      IBUG1 = 0
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*

*   Diagonal case IC = IR. Start with computing the contribution
*   ENONSYM that does not include the symbolic orbital. 

      NVCOEF = 0
      ENONSYM = 0.D0
*
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*
      DO 7 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          IF (LABEL(1,I) .LT. NTYPE(3,IC) .AND. 
     :        LABEL(2,I) .LT. NTYPE(3,IC) .AND.
     :        LABEL(3,I) .LT. NTYPE(3,IC) .AND. 
     :	      LABEL(4,I) .LT. NTYPE(3,IC)) THEN 
            IF (LSMS) THEN
              IF (LABEL(5,I) .EQ. 1) THEN
                CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                ENONSYM = ENONSYM - TGRL1*TGRL2*ATWINV*VCOEFF
                                  ! negative contributions?  
              ENDIF
            ENDIF
            CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
            ENONSYM = ENONSYM + TEGRAL*VCOEFF
            WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!           write(*,*)ENONSYM
          
          END IF
        ENDIF
    7 CONTINUE

      DO 8 J = NORB, 1, -1      
        EMTBLOCK(J,J) = EMTBLOCK(J,J) + ENONSYM
    8 CONTINUE

*
*   Diagonal case IC = IR. compute the contribution
*   of the symbolic orbital. 


      DO 9 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          NCTEC = NCTEC + NTYPE(2,IC)  
          IF (LABEL(1,I) .GE. NTYPE(3,IC) .OR. 
     :        LABEL(2,I) .GE. NTYPE(3,IC) .OR.
     :        LABEL(3,I) .GE. NTYPE(3,IC) .OR. 
     :	      LABEL(4,I) .GE. NTYPE(3,IC)) THEN 
            DO 10 J = NORB, 1, -1      
              IF (J .NE. NORB) THEN
!                WRITE(*,*)"J=",J
                IF (LABEL(1,I) .GE. NTYPE(3,IC)) THEN
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                END IF
                IF (LABEL(2,I) .GE. NTYPE(3,IC)) THEN
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                END IF
                IF (LABEL(3,I) .GE. NTYPE(3,IC)) THEN
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                END IF
                IF (LABEL(4,I) .GE. NTYPE(3,IC)) THEN
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                END IF
              ENDIF
              IF (LSMS) THEN
                IF (LABEL(5,I) .EQ. 1) THEN
                  CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                  CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                  EMTBLOCK(J,J) = EMTBLOCK(J,J) - TGRL1*TGRL2*ATWINV*VCOEFF
                                               ! negative contributions
                ENDIF
              ENDIF
              CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
              EMTBLOCK(J,J) = EMTBLOCK(J,J) + TEGRAL*VCOEFF
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!             write(*,*)ENONSYM
          
   10       CONTINUE
          END IF
        ENDIF
    9 CONTINUE
    
*   Here comes the off-diagonal part within the diagonal symbolic block.
      NVCOEF = 0
      ENONSYM = 0.D0
*
      WRITE(*,*)"Here comes the off-diagonal part within the diagonal 
     : symbolic block in matrixblock22"
      CALL RKCO_GG (IC, IC-1, CORD, INCOR, 1)

      DO 11 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          DO 12 J = NORB, 2, -1
            DO 13 K = J - 1, 1, -1      
              IF (K .NE. (NORB - 1)) THEN
                IF (LABEL(5,I) .EQ. 0) THEN
                  LABEL(2,I) = K + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
!                END IF
!                IF (LABEL(5,I) .EQ. 1) THEN
                ELSE
                  LABEL(3,I) = K + NTYPE(3,IC) - 1
                  LABEL(4,I) = J + NTYPE(3,IC) - 1
                END IF
              ENDIF
!              write(*,*)"J=",J,"K=",K,"NORB=",NORB,"NORB-1=",NORB-1
!              write(*,*)(LABEL(M,I),M=1,6)
              
              IF (LSMS) THEN
                IF (LABEL(5,I) .EQ. 1) THEN
                  CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                  CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                  EMTBLOCK(J,K) = EMTBLOCK(J,K) - TGRL1*TGRL2*ATWINV*VCOEFF
                                               ! negative contributions
                ENDIF
              ENDIF
              CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)
              EMTBLOCK(J,K) = EMTBLOCK(J,K) + TEGRAL*VCOEFF
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!             write(*,*)ENONSYM
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

        DO 14 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
              NMCBP = NMCBP + NTYPE(2,IC)
              ITYPE = ABS (LABEL(6,I))
            DO 15 J = NORB, 2, -1
              DO 16 K = J - 1, 1, -1      
                IF (K .NE. (NORB - 1)) THEN
                  IF (LABEL(5,I) .EQ. 0) THEN
                    LABEL(2,I) = K + NTYPE(3,IC) - 1
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  ELSE
                    LABEL(3,I) = K + NTYPE(3,IC) - 1
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  END IF
                ENDIF
                IF (ITYPE .EQ. 1) THEN
                  CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ELSEIF (ITYPE .EQ. 2) THEN
                  CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ELSEIF (ITYPE .EQ. 3) THEN
                  CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ELSEIF (ITYPE .EQ. 4) THEN
                  CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ELSEIF (ITYPE .EQ. 5) THEN
                  CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ELSEIF (ITYPE .EQ. 6) THEN
                  CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                         LABEL(3,I), LABEL(4,I),
     :                         LABEL(5,I), TEGRAL)
                ENDIF 
                CONTR = COEFF(I)*TEGRAL
                IF (LABEL(6,I) .GT. 0) THEN
                  EMTBLOCK(J,K) = EMTBLOCK(J,K) + CONTR
                ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                  NCORE = NCORE + 1
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
