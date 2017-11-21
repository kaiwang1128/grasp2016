************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK23(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 2 and 3                 *         
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

      WRITE(*,*) 'In matrixblock23, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0

*   Call onescalar to       

*   Here comes the off-diagonal part within the diagonal symbolic block.
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL) 
 
      WRITE(*,592) IC,IR,DBLE(TSHELL(1)),'I(',NP(IA),NH(IA),',',
     :                   NP(IB),NH(IB),')'        

*   Ensure that the indices are in `canonical' order
      NORBUU = 0
      NORBUL = 0
      NORBL = 0
      IF (IA .GT. IB) THEN
        ISWAP = IB
        IB = IA
        IA = ISWAP
      ENDIF
*    In matrixblock 23, for symbolic CSF of Type 2, the total number of
*    symbolic CSF (and symbolic orbitals) is defined by NORBL. For symbolic 
*    CSF of Type 3, there are two kinds of symbolic orbitals. the total  
*    number of the first kind of symbolic orbitali defined by NORBUL. The  
*    total number of the sencond kind of symbolic orbital is defined by 
*    NORBUU. The total number of symbolic CSF of Type 3 is defined by 
*    N = NORBUL * NORBUU

      NORBUU = NTYPE(6,IC) - NTYPE(5,IC) + 1
      NORBUL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      NORBL = NTYPE(4,IR) - NTYPE(3,IR) + 1

      TCOEFF = DBLE(TSHELL(1))

*     2s [7-3]p (type 2) and [3-7]s [3-7]p
*     NORBL =5, the total number of symbolic CSF is 5.
*     NORBUL = 5, NORBL =5. The total number of symbolic CSF of Type 3 
*     is defined by N = NORBUL * NORBUU =125.
*     Among the 125*5 cases, one electron operator for some specific cases 
*     exists, such as, for the case---2s 7p and 3s 7p, but does not exist 
*     for cases, such as, for the case---2s 7p and 3s 6p. Therefore, when 
*     the second orbital ([7-3]p) of symbolic CSF of Type 2 equals to the
*     the first orbital ([3-7]s) or the second 
*     orbital ([7-3]p) of symbolic CSF of Type 3, one electron operator 
*     exists. This condition is defined by 
*     IF ((K+ NTYPE(5,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN or 
*     IF ((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN.
*     Among the 125*5 cases, one electron operator for 5*5 cases. 
*     When the second orbital of symbolic CSF of Type 2 equals to the
*     the first orbital ([3-7]s) or the second 
*     orbital ([7-3]p) of symbolic CSF of Type 3, the left orbital of 
*     of symbolic CSF of Type 3 is defined by IB = K + NTYPE(5,IC) - 1
*     or IB = J + NTYPE(3,IC) - 1. IB should be changed for different 
*     cases. IA is always the first orbital of symbolic CSF of Type 2.

      IF (DABS(TCOEFF) .GT. CUTOFF) THEN
        N = NORBUL * NORBUU 
        DO J = NORBUL, 1, -1
          DO K = NORBUU, 1, -1
            N = N -1
            DO M = NORBL, 1, -1              
              IF ((K+ NTYPE(5,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN
                NCOEC = NCOEC + 1
!              IA = I + NTYPE(3,IC) - 1
                IB = J + NTYPE(3,IC) - 1
                CALL IABINT(IA,IB,TEGRAL)
                EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*TCOEFF
                WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                       NP(IB),NH(IB),')'
                IF (LNMS) THEN
                  CALL KEINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*ATWINV*TCOEFF
                ENDIF
                IF (LVP) THEN
                  CALL VPINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*TCOEFF
                ENDIF
              END IF
              IF ((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN
                NCOEC = NCOEC + 1
                IB = K + NTYPE(5,IC) - 1
                CALL IABINT(IA,IB,TEGRAL)
                EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*TCOEFF
                WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                       NP(IB),NH(IB),')'
                IF (LNMS) THEN
                  CALL KEINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*ATWINV*TCOEFF
                ENDIF
                IF (LVP) THEN
                  CALL VPINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(N,M) = EMTBLOCK(N,M) + TEGRAL*TCOEFF
                ENDIF
              END IF
            END DO   
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

*   Here comes the off-diagonal part within the diagonal symbolic block.

*   The Rk has always two or four components
*   For 4s 4p and 2s 4p, The Rk has four components
*             1.000000000 R0(2s 4p ,4s 4p )
*            -0.333333333 R1(2s 4s ,4p 4p )
*             2.000000000 R0(1s 2s ,1s 4s )
*            -1.000000000 R0(1s 1s ,4s 2s )
*   For 4s 4p and 2s 3p, The Rk has two components
*             1.000000000 R0(2s 3p ,4s 4p )
*            -0.333333333 R1(2s 3p ,4p 4s )
*   If one orbital of CSF in type 3 is the same with one orbital of CSF 
*   in type 2, The Rk has four components. 
*   Otherwise, The Rk has two components.

*   For 4s 3p and 2s 3p
*             0.333333333 R1(2s 3p ,3p 4s )
*            -1.000000000 R0(2s 3p ,4s 3p )
*            -2.000000000 R0(1s 2s ,1s 4s )
*             1.000000000 R0(1s 1s ,2s 4s )
*     The absolute values of four components, with opposite in sign, 
*     are the same with the absolute values of four components 
*     of 4s 4p and 2s 4p.



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
     : symbolic block in matrixblock23"
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
