************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK22(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
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
     : ,TMPBLOCK(MAXSPAN,MAXSPAN)

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

      WRITE(*,*) 'In matrixblock22, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0
*   
*   Accumulate the contribution from the one-body operators:
*   kinetic energy, electron-nucleus interaction; update the
*   angular integral counter

*   Call onescalar to       
      TSHELL = 0.D0
      IF (NTYPE(3,IC).EQ.NTYPE(3,IR))THEN
        CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
        IF(DBLE(TSHELL(1)).NE.0.D0)THEN
          WRITE(*,*) 'ONESCALAR In matrixblock22 IC,IR',IC,IR

          WRITE(*,592) IC,IR,DBLE(TSHELL(1))
     :               ,'I(',NP(IA),NH(IA),',',NP(IB),NH(IB),')'
          DO J = NTYPE(2,IC), NTYPE(2,IC), -1
            DO K = NTYPE(2,IR), NTYPE(2,IR), -1
               TCOEFF = DBLE(TSHELL(1))
               IF (DABS (TCOEFF) .GT. CUTOFF) THEN
                  NCOEC = NCOEC + 1
                  CALL IABINT (IA, IB, TEGRAL)

                  WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                 NP(IB),NH(IB),')'

                        !------------------------
                  EMTBLOCK(K,J) = EMTBLOCK(K,J) + TEGRAL*TCOEFF
                  IF (LNMS) THEN
                     CALL KEINT (IA, IB, TEGRAL)
                        !------------------------
                     EMTBLOCK(K,J) = EMTBLOCK(K,J) + 
     :                TEGRAL*ATWINV*TCOEFF
                  ENDIF
                  IF (LVP) THEN
                     CALL VPINT (IA, IB, TEGRAL)
                        !------------------------
                     EMTBLOCK(K,J) = EMTBLOCK(K,J) + TEGRAL*TCOEFF
                  ENDIF
               END IF
            END DO
          END DO
          DO J = NTYPE(2,IC), 1, -1
            DO K = NTYPE(2,IR), 1, -1
              IF(J.EQ.K)THEN
              EMTBLOCK(K,J) = EMTBLOCK(NTYPE(2,IR),NTYPE(2,IC))
              END IF
            END DO
          END DO
        END IF
      END IF
       DO 52 K = NTYPE(2,IR), 1, -1
         DO 51 J = NTYPE(2,IC), 1, -1

           write(*,*)IC,IR,K,J,EMTBLOCK(K,J)
   51    continue
   52  continue    
      IBUG1 = 0
    
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*
      WRITE(*,*) 'RKCO_GG CORD In matrixblock22 IC,IR',IC,IR

      NVCOEF = 0
      MMIN=MIN(NTYPE(3,IC),NTYPE(3,IR))

      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      DO  I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
        
            CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)

           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
        END IF
      END DO

      DO  I = 1, NVCOEF
        MRANK = 0
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          DO J = NTYPE(2,IC), 1, -1
            DO K = NTYPE(2,IR), 1, -1
              IF(J .EQ. NTYPE(2,IC) .AND. K .EQ. NTYPE(2,IR))then
                NUMB = 0
                DO, MTYPE=1,4
                  IF(LABEL(MTYPE,I).LT.MMIN)THEN
                    NUMB = NUMB + 1                    
                  END IF
                END DO
                IF(NUMB.NE.2.AND.NUMB.NE.4)THEN
                  WRITE(*,*)NUMB
                  WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 1 IN MATRIXBLOCK22'
                  STOP
                END IF
              END IF
              IF(NUMB.EQ.4.AND.J.EQ.NTYPE(2,IC)
     :         .AND.K.EQ.NTYPE(2,IR))THEN
                MRANK = 1               
              ELSE IF(NUMB.EQ.2.AND.J.EQ.NTYPE(2,IC)
     :               .AND.K.EQ.NTYPE(2,IR))THEN
                IF(LABEL(1,I).EQ.NTYPE(4,IR))THEN
                  IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 2              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 3              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 4              
                  ELSE 
                    WRITE(*,*)'STOP 2 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 5              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 6              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 7              
                  ELSE 
                    WRITE(*,*)'STOP 3 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 8              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 9              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 10              
                  ELSE 
                    WRITE(*,*)'STOP 4 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 11              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 12              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 13              
                  ELSE 
                    WRITE(*,*)'STOP 5 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE 
                  WRITE(*,*)LABEL(1,I),LABEL(2,I),LABEL(3,I)
     :            ,LABEL(4,I),MRANK             
                  WRITE(*,*)'STOP 6 IN MATRIXBLOCK22'                  
                  STOP
                END IF
              END IF
              IF(MRANK.EQ.2)THEN
                LABEL(1,I)= K + NTYPE(3,IR) -1
                LABEL(2,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.3)THEN
                LABEL(1,I)= K + NTYPE(3,IR) -1
                LABEL(3,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.4)THEN
                LABEL(1,I)= K + NTYPE(3,IR) -1
                LABEL(4,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.5)THEN
                LABEL(2,I)= K + NTYPE(3,IR) -1
                LABEL(1,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.6)THEN
                LABEL(2,I)= K + NTYPE(3,IR) -1
                LABEL(3,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.7)THEN
                LABEL(2,I)= K + NTYPE(3,IR) -1
                LABEL(4,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.8)THEN
                LABEL(3,I)= K + NTYPE(3,IR) -1
                LABEL(1,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.9)THEN
                LABEL(3,I)= K + NTYPE(3,IR) -1
                LABEL(2,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.10)THEN
                LABEL(3,I)= K + NTYPE(3,IR) -1
                LABEL(4,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.11)THEN
                LABEL(4,I)= K + NTYPE(3,IR) -1
                LABEL(1,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.12)THEN
                LABEL(4,I)= K + NTYPE(3,IR) -1
                LABEL(2,I)= J + NTYPE(3,IC) -1
              ELSE IF(MRANK.EQ.13)THEN
                LABEL(4,I)= K + NTYPE(3,IR) -1
                LABEL(3,I)= J + NTYPE(3,IC) -1
              END IF              
              IF((MRANK.EQ.1.AND.K.EQ.J).OR.MRANK.GE.2)THEN
                IF (LSMS) THEN
                  IF (LABEL(5,I) .EQ. 1) THEN
                    CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                    CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                      EMTBLOCK(K,J) = EMTBLOCK(K,J) - 
     :                                TGRL1*TGRL2*ATWINV*VCOEFF
                     WRITE(*,*)IC,IR,MRANK,K,J,EMTBLOCK(K,J)
                  END IF
                END IF
                CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
                EMTBLOCK(K,J) = EMTBLOCK(K,J) + TEGRAL*VCOEFF
                WRITE(*,594)IC,IR,int(MRANK),K,J,VCOEFF,'R',
     :                   LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'

              END IF
              
            END DO
          END DO
        END IF
      END DO
      
      IF(NTYPE(3,IC).EQ.NTYPE(3,IR).AND.
     :   NTYPE(2,IC).EQ.NTYPE(2,IR))THEN
        TMPBLOCK=0.D0
        DO, J = NTYPE(2,IC), 1, -1
          DO, K = NTYPE(2,IR), 1, -1
            TMPBLOCK(J,K)=EMTBLOCK(K,J)
          END DO
        END DO
        EMTBLOCK=0.D0
        DO, J = NTYPE(2,IC), 1, -1
          DO, K = NTYPE(2,IR), 1, -1
            EMTBLOCK(K,J)=TMPBLOCK(K,J)
          END DO
        END DO

      END IF

       DO 53 J = NTYPE(2,IC), 1, -1
         DO 54 K = NTYPE(2,IR), 1, -1

           write(*,*)IC,IR,J,K,EMTBLOCK(K,J)
   54    continue
   53  continue    



    
      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
      
      WRITE(*,*) 'RKCO_GG BREID In matrixblock22, IC,IR',IC,IR
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
*     
      IF(NTYPE(4,IC).eq.NTYPE(4,IR).and.
     :   NTYPE(2,IC).eq.NTYPE(2,IR))THEN
       NVCOEF = 0      
       CALL RKCO_GG (IC, IR, BREID, 1, 2)
       DO  I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))        
                  IF (ITYPE .EQ. 1) THEN

                     CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 2) THEN
                     CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 3) THEN
                     CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 4) THEN
                     CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 5) THEN
                     CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 6) THEN
                     CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ENDIF 

           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
        END IF
      END DO

      DO  I = 1, NVCOEF
        MRANK = 0
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))
          DO J = NTYPE(2,IC), 1, -1
            DO K = NTYPE(2,IR), 1, -1
              IF(J .EQ. NTYPE(2,IC) .AND. K .EQ. NTYPE(2,IR))then
                NUMB = 0
                DO, MTYPE=1,4
                  IF(LABEL(MTYPE,I).LT.MMIN)THEN
                    NUMB = NUMB + 1                    
                  END IF
                END DO
                IF(NUMB.NE.2.AND.NUMB.NE.4)THEN
                  WRITE(*,*)NUMB
                  WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 11 IN MATRIXBLOCK22'
                  STOP
                END IF
              END IF
              IF(NUMB.EQ.4.AND.J.EQ.NTYPE(2,IC)
     :         .AND.K.EQ.NTYPE(2,IR))THEN
                MRANK = 1               
              ELSE IF(NUMB.EQ.2.AND.J.EQ.NTYPE(2,IC)
     :               .AND.K.EQ.NTYPE(2,IR))THEN
                IF(LABEL(1,I).EQ.NTYPE(4,IR))THEN
                  IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 2              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 3              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 4              
                  ELSE 
                    WRITE(*,*)'STOP 12 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 5              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 6              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 7              
                  ELSE 
                    WRITE(*,*)'STOP 13 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 8              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 9              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 10              
                  ELSE 
                    WRITE(*,*)'STOP 14 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 11              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 12              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 13              
                  ELSE 
                    WRITE(*,*)'STOP 15 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE 
                  WRITE(*,*)LABEL(1,I),LABEL(2,I),LABEL(3,I)
     :            ,LABEL(4,I),MRANK             
                  WRITE(*,*)'STOP 16 IN MATRIXBLOCK22'                  
                  STOP
                END IF
              END IF
              IF(J.EQ.K)THEN  
                IF(MRANK.EQ.2)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.3)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.4)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.5)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.6)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.7)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.8)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.9)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.10)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.11)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.12)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.13)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                END IF              
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
                    EMTBLOCK(K,J) = EMTBLOCK(K,J) + CONTR
                      WRITE(*,594)IC,IR,int(MRANK),K,J,VCOEFF,'R',
     :                     LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  ELSE
!                          ...It comes here only when ic=ir=1
!                             clue: rkco<-breid<-talk<-label(6,i)
                    NCORE = NCORE + 1
                    ELSTO = ELSTO + CONTR
                  END IF                
              END IF              
            END DO
          END DO
         END IF
       END DO
       
       NVCOEF = 0      
       CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
       DO  I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))        
                  IF (ITYPE .EQ. 1) THEN

                     CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 2) THEN
                     CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 3) THEN
                     CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 4) THEN
                     CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 5) THEN
                     CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 6) THEN
                     CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ENDIF 

           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
        END IF
      END DO

      DO  I = 1, NVCOEF
        MRANK = 0
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))
          DO J = NTYPE(2,IC), 1, -1
            DO K = NTYPE(2,IR), 1, -1
              IF(J .EQ. NTYPE(2,IC) .AND. K .EQ. NTYPE(2,IR))then
                NUMB = 0
                DO, MTYPE=1,4
                  IF(LABEL(MTYPE,I).LT.MMIN)THEN
                    NUMB = NUMB + 1                    
                  END IF
                END DO
                IF(NUMB.NE.2.AND.NUMB.NE.4)THEN
                  WRITE(*,*)NUMB
                  WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 11 IN MATRIXBLOCK22'
                  STOP
                END IF
              END IF
              IF(NUMB.EQ.4.AND.J.EQ.NTYPE(2,IC)
     :         .AND.K.EQ.NTYPE(2,IR))THEN
                MRANK = 1               
              ELSE IF(NUMB.EQ.2.AND.J.EQ.NTYPE(2,IC)
     :               .AND.K.EQ.NTYPE(2,IR))THEN
                IF(LABEL(1,I).EQ.NTYPE(4,IR)-1)THEN
                  LABEL(1,I)= NTYPE(4,IR)                  
                  IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 2              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 3              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 4              
                  ELSE 
                    WRITE(*,*)'STOP 12 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR)-1)THEN 
                  LABEL(2,I)= NTYPE(4,IR)                  
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 5              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 6              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 7              
                  ELSE 
                    WRITE(*,*)'STOP 13 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR)-1)THEN 
                  LABEL(3,I)= NTYPE(4,IR)                  
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 8              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 9              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 10              
                  ELSE 
                    WRITE(*,*)'STOP 14 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR)-1)THEN 
                  LABEL(4,I)= NTYPE(4,IR)                  
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 11              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 12              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 13              
                  ELSE 
                    WRITE(*,*)'STOP 15 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE 
                  WRITE(*,*)LABEL(1,I),LABEL(2,I),LABEL(3,I)
     :            ,LABEL(4,I),MRANK             
                  WRITE(*,*)'STOP 16 IN MATRIXBLOCK22'                  
                  STOP
                END IF
              END IF
              IF(J.NE.K)THEN  
                IF(MRANK.EQ.2)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.3)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.4)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.5)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.6)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.7)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.8)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.9)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.10)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.11)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.12)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.13)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                END IF              
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
                    EMTBLOCK(K,J) = EMTBLOCK(K,J) + CONTR
                      WRITE(*,594)IC,IR,int(MRANK),K,J,VCOEFF,'R',
     :                     LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  ELSE
!                          ...It comes here only when ic=ir=1
!                             clue: rkco<-breid<-talk<-label(6,i)
                    NCORE = NCORE + 1
                    ELSTO = ELSTO + CONTR
                  END IF                
              END IF              
            END DO
          END DO
         END IF
       END DO
      ELSE IF(NTYPE(4,IC).NE.NTYPE(4,IR).OR.
     :   NTYPE(2,IC).NE.NTYPE(2,IR))THEN
       NVCOEF = 0      
       CALL RKCO_GG (IC, IR, BREID, 1, 2)
       DO  I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))        
                  IF (ITYPE .EQ. 1) THEN

                     CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 2) THEN
                     CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 3) THEN
                     CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 4) THEN
                     CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 5) THEN
                     CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ELSEIF (ITYPE .EQ. 6) THEN
                     CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                            LABEL(3,I), LABEL(4,I),
     :                            LABEL(5,I), TEGRAL)
                  ENDIF 

           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
        END IF
      END DO

      DO  I = 1, NVCOEF
        MRANK = 0
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))
          DO J = NTYPE(2,IC), 1, -1
            DO K = NTYPE(2,IR), 1, -1
              IF(J .EQ. NTYPE(2,IC) .AND. K .EQ. NTYPE(2,IR))then
                NUMB = 0
                DO, MTYPE=1,4
                  IF(LABEL(MTYPE,I).LT.MMIN)THEN
                    NUMB = NUMB + 1                    
                  END IF
                END DO
                IF(NUMB.NE.2.AND.NUMB.NE.4)THEN
                  WRITE(*,*)NUMB
                  WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 11 IN MATRIXBLOCK22'
                  STOP
                END IF
              END IF
              IF(NUMB.EQ.4.AND.J.EQ.NTYPE(2,IC)
     :         .AND.K.EQ.NTYPE(2,IR))THEN
                MRANK = 1               
              ELSE IF(NUMB.EQ.2.AND.J.EQ.NTYPE(2,IC)
     :               .AND.K.EQ.NTYPE(2,IR))THEN
                IF(LABEL(1,I).EQ.NTYPE(4,IR))THEN
                  IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 2              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 3              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 4              
                  ELSE 
                    WRITE(*,*)'STOP 32 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 5              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 6              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 7              
                  ELSE 
                    WRITE(*,*)'STOP 33 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 8              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 9              
                  ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 10              
                  ELSE 
                    WRITE(*,*)'STOP 34 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN 
                  IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 11              
                  ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 12              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                    MRANK = 13              
                  ELSE 
                    WRITE(*,*)'STOP 35 IN MATRIXBLOCK22'                  
                    STOP
                  END IF
                ELSE 
                  WRITE(*,*)LABEL(1,I),LABEL(2,I),LABEL(3,I)
     :            ,LABEL(4,I),MRANK             
                  WRITE(*,*)'STOP 36 IN MATRIXBLOCK22'                  
                  STOP
                END IF
              END IF
!              IF(J.EQ.K)THEN  
                IF(MRANK.EQ.2)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.3)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.4)THEN
                  LABEL(1,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.5)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.6)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.7)THEN
                  LABEL(2,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.8)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.9)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.10)THEN
                  LABEL(3,I)= K + NTYPE(3,IR) -1
                  LABEL(4,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.11)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(1,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.12)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(2,I)= J + NTYPE(3,IC) -1
                ELSE IF(MRANK.EQ.13)THEN
                  LABEL(4,I)= K + NTYPE(3,IR) -1
                  LABEL(3,I)= J + NTYPE(3,IC) -1
                END IF              
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
                    EMTBLOCK(K,J) = EMTBLOCK(K,J) + CONTR
                      WRITE(*,594)IC,IR,int(MRANK),K,J,VCOEFF,'R',
     :                     LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  ELSE
!                          ...It comes here only when ic=ir=1
!                             clue: rkco<-breid<-talk<-label(6,i)
                    NCORE = NCORE + 1
                    ELSTO = ELSTO + CONTR
                  END IF                
!              END IF              
            END DO
          END DO
         END IF
       END DO
      ELSE
        WRITE(*,*)'STOP 36 IN MATRIXBLOC23'
        STOP        
      END IF 
       DO 55 J = NTYPE(2,IC), 1, -1
         DO 56 K = NTYPE(2,IR), 1, -1

           write(*,*)IC,IR,J,K,EMTBLOCK(K,J)
   56    continue
   55  continue    
            
      END IF

*   Diagonal case IC = IR. Start with computing the contribution
*   ENONSYM that does not include the symbolic orbital. For
*   NTYPE = 2 there is only one symbolic orbital 
*

  592 FORMAT (2I5,F15.9,1X,1A2,1I1,1A2,1A1,I1,1A2,1A1)
  593 FORMAT (2I5,F15.9,1X,1A,1I1,1A1,1I1,1A2,1I1,1A2,
     :        1A1,1I1,1A2,1I1,1A2,1A1) 
  594 FORMAT (5I5,F15.9,1X,1A,1I1,1A1,1I2,1A2,1I2,1A2,
     :        1A1,1I2,1A2,1I2,1A2,1A1) 

      RETURN
      END
