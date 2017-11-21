************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK3(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!      SUBROUTINE ONESCALAR11(IC,IR,NCOEC,ELEMNT)
!      SUBROUTINE RKCO_GG11(IC,IR,INCOR,NCTEC,INC2,NMCBP,NCORE,
!     :           ELSTO,ELEMNT)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 1 and 2                 *         
*                                                                      *
*   Written by Per Jönsson & Kai Wang                         June 2017*
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
      INTEGER I,J,K,L,M,N,O,NN,ITMP1,ITMP2
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
      WRITE(*,*) 'ONESCALAR In matrixblock3, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0
      TSHELL = 0.D0

*   Call onescalar to       
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
*
       DO IA =1,NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
          CALL IABINT(IA,IA,TEGRAL)
          WRITE(*,592) IC,IC,TCOEFF,'I(',
     :                 NP(IA),NH(IA),',',NP(IA),NH(IA),')'  

        END IF       
       END DO
       WRITE(*,*)
       NORBUU = NTYPE(6,IC) - NTYPE(5,IC) + 1
       NORBUL = NTYPE(4,IC) - NTYPE(3,IC) + 1
*      Computing the contribution from the one body operator,
*      when 1) two orbitals of one state of the matrix elements
*      are equal to two orbitals of the other state, respectively; (I.EQ.K.AND.J.EQ.M)
*      2) one orbital of one state of the matrix elements
*      is equal to one orbital of the other state, and IA is equal to NTYPE(4,IC);
*      3) one orbital of one state of the matrix elements
*      is equal to one orbital of the other state, and IA is equal to NTYPE(6,IC).
       DO IA = 1,NW
        TCOEFF = DBLE(TSHELL(IA))
        IF (DABS (TCOEFF) .GT. CUTOFF) THEN
         N=NTYPE(2,IC)+1          
         DO, I=NORBUL,1,-1
          DO, J=NORBUU,1,-1
           N=N-1
           L=N+1
           DO, K=I,1,-1 
            DO, M=J,1,-1
             L=L-1  
             IF((I.EQ.K.AND.J.EQ.M).OR.                   
     :        (I.NE.K.AND.J.EQ.M.AND.IA.EQ.NTYPE(4,IC)).OR.
     :        (I.EQ.K.AND.J.NE.M.AND.IA.EQ.NTYPE(6,IC)))THEN
              IF(IA.EQ.NTYPE(6,IC).AND.I.EQ.K)THEN
               ITMP1=J+NTYPE(5,IC)-1
               ITMP2=M+NTYPE(5,IC)-1  
              ELSE IF(IA.EQ.NTYPE(4,IC).AND.J.EQ.M)THEN
               ITMP1=I+NTYPE(3,IC)-1
               ITMP2=K+NTYPE(3,IC)-1
              ELSE
               ITMP1=IA
               ITMP2=IA               
              END IF
              WRITE(*,*) IC,IC,NP(I+NTYPE(3,IC)-1),NH(I+NTYPE(3,IC)-1),
     :                   ',',NP(J+NTYPE(5,IC)-1),NH(J+NTYPE(5,IC)-1),
     :                   ',',NP(K+NTYPE(3,IC)-1),NH(K+NTYPE(3,IC)-1),
     :                   ',',NP(M+NTYPE(5,IC)-1),NH(M+NTYPE(5,IC)-1)  
              WRITE(*,592) IC,IC,TCOEFF,'I(',NP(ITMP2),NH(ITMP2),',',
     :                    NP(ITMP1),NH(ITMP1),')'  
              CALL IABINT(ITMP2,ITMP1,TEGRAL)
              EMTBLOCK(L,N) = EMTBLOCK(L,N) + TEGRAL*TCOEFF
              IF (LNMS) THEN
               CALL KEINT(ITMP2,ITMP1,TEGRAL)
               EMTBLOCK(L,N) = EMTBLOCK(L,N) + TEGRAL*ATWINV*TCOEFF
              END IF
              IF (LVP) THEN
               CALL VPINT(ITMP2,ITMP1,TEGRAL)
               EMTBLOCK(L,N) = EMTBLOCK(L,N)+ TEGRAL*TCOEFF
              END IF
              write(*,*)IC,IR,N,L,EMTBLOCK(L,N) 
              WRITE(*,*)

             END IF
            END DO
           END DO
          END DO
         END DO
        END IF
       END DO
         WRITE(*,*)
         N=NTYPE(2,IC)+1          
         DO, I=NORBUL,1,-1
          DO, J=NORBUU,1,-1
           N=N-1
           L=N+1
           DO, K=I,1,-1 
            DO, M=J,1,-1
             L=L-1  
              write(*,*)IC,IR,N,L,EMTBLOCK(L,N) 
            END DO
           END DO
          END DO
         END DO
       
       
      WRITE(*,*) 'RKCO_GG CORD In matrixblock3, IC,IR',IC,IR
       
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*

*   Diagonal case IC = IR. Start with computing the contribution
*   ENONSYM that does not include the symbolic orbital. 

      NVCOEF = 0
      
*
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
*     
    
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
         N=NTYPE(2,IC)+1          
         DO, O=NORBUL,1,-1
          DO, J=NORBUU,1,-1
           N=N-1
           L=N+1
           DO, K=O,1,-1 
            DO, M=J,1,-1
              L=L-1
              IF(O .EQ. NORBUL .AND. J .EQ. NORBUU .AND.
     :           O.EQ.K.AND.J.EQ.M)then
                NUMB = 0
                DO, MTYPE=1,4
                  IF(LABEL(MTYPE,I).LT.NTYPE(3,IC))THEN
                    NUMB = NUMB + 1                    
                  END IF
                END DO
                IF(NUMB.NE.2.AND.NUMB.NE.4.AND.NUMB.NE.0)THEN
                  WRITE(*,*)NUMB
                  WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 1 IN MATRIXBLOCK33'
                  STOP
                END IF
                IF(NUMB.EQ.4)THEN
                  MRANK = 1               
                ELSE IF(NUMB.EQ.2)THEN
                  IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(4,IC))THEN
                      MRANK = 2              
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(4,IC))THEN
                      MRANK = 3      
                  ELSE IF(LABEL(2,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 4      
                  ELSE IF(LABEL(3,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 5      
                  ELSE 
                    WRITE(*,*)'STOP 2 IN MATRIXBLOCK33'                  
                    STOP
                  END IF
                ELSE IF(NUMB.EQ.0)THEN
                  IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(1,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(3,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 6      
                  ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(1,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(2,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 7      
                  ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(2,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 8      
                  ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :               LABEL(3,I).EQ.NTYPE(6,IC).AND.
     :               LABEL(4,I).EQ.NTYPE(6,IC))THEN
                      MRANK = 9      
                  ELSE 
                    WRITE(*,594)IC,IR,int(MRANK),K,J,
     :                   VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                    WRITE(*,*)'STOP 3 IN MATRIXBLOCK33'                  
                    STOP                      
                  END IF
                ELSE  
                  WRITE(*,594)IC,IR,int(MRANK),K,J,
     :                   VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  WRITE(*,*)'STOP 4 IN MATRIXBLOCK33'                  
                  STOP                                        
                END IF
              END IF
              SELECT CASE (MRANK)
                CASE (1)
                  GOTO 50
                CASE (2)
                  LABEL(2,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = O + NTYPE(3,IC) - 1                                  
                CASE (3)
                  LABEL(3,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = O + NTYPE(3,IC) - 1                                  
                CASE (4)
                  LABEL(2,I) = J + NTYPE(5,IC) - 1                                  
                  LABEL(4,I) = M + NTYPE(5,IC) - 1                                  
                CASE (5)
                  LABEL(3,I) = M + NTYPE(5,IC) - 1                                  
                  LABEL(4,I) = J + NTYPE(5,IC) - 1                                  
                CASE (6)
                  LABEL(2,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = O + NTYPE(3,IC) - 1                                  
                  LABEL(1,I) = M + NTYPE(5,IC) - 1                                  
                  LABEL(3,I) = J + NTYPE(5,IC) - 1                                  
                CASE (7)
                  LABEL(3,I) = O + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(1,I) = M + NTYPE(5,IC) - 1                                  
                  LABEL(2,I) = J + NTYPE(5,IC) - 1                                  
                CASE (8)
                  LABEL(1,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(3,I) = O + NTYPE(3,IC) - 1                                  
                  LABEL(2,I) = M + NTYPE(5,IC) - 1                                  
                  LABEL(4,I) = J + NTYPE(5,IC) - 1                                  
                CASE (9)
                  LABEL(1,I) = K + NTYPE(3,IC) - 1                                  
                  LABEL(2,I) = O + NTYPE(3,IC) - 1                                  
                  LABEL(3,I) = M + NTYPE(5,IC) - 1                                  
                  LABEL(4,I) = J + NTYPE(5,IC) - 1                                  
                CASE DEFAULT
                  WRITE(*,*)'STOP 5 IN MATRIXBLOCK33'                  
                  STOP                                        
              END SELECT              
   50        CONTINUE
             IF((O.EQ.K.AND.J.EQ.M).OR.                   
     :        (J.EQ.M.AND.(MRANK.EQ.2.OR.MRANK.EQ.3)).OR.
     :        (O.EQ.K.AND.(MRANK.EQ.4.OR.MRANK.EQ.5)).OR.
     :        (MRANK.GE.6.AND.MRANK.LE.9))THEN

                WRITE(*,594)IC,IR,int(MRANK),L,N,VCOEFF,'R',
     :                   LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                   WRITE(*,*)IC,IR,MRANK,L,N,EMTBLOCK(L,N),'111'                 
               IF (LSMS) THEN
                  IF (LABEL(5,I) .EQ. 1) THEN
                     CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                     write(*,*)LABEL(2,I), LABEL(4,I), TGRL2
                     CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                     write(*,*)LABEL(2,I), LABEL(4,I), TGRL2
                     EMTBLOCK(L,N) = EMTBLOCK(L,N) 
     :                              - TGRL1*TGRL2*ATWINV*VCOEFF
                     WRITE(*,*)IC,IR,MRANK,L,N,EMTBLOCK(L,N),'222'
                     WRITE(*,*)TGRL1,TGRL2,ATWINV,VCOEFF
                 END IF
               END IF
               CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
                WRITE(*,594)IC,IR,int(MRANK),L,N,VCOEFF,'R',
     :                   LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
               WRITE(*,*)IC,IR,MRANK,L,N,EMTBLOCK(L,N),'333'
               write(*,*)TEGRAL,VCOEFF
               EMTBLOCK(L,N) = EMTBLOCK(L,N) + TEGRAL*VCOEFF
               WRITE(*,*)IC,IR,MRANK,L,N,EMTBLOCK(L,N),'444'
             END IF 
                
            END DO
           END DO
          END DO
         END DO
        END IF
      END DO
      

      IBUG1 = 0
         WRITE(*,*)
         N=NTYPE(2,IC)+1          
         DO, I=NORBUL,1,-1
          DO, J=NORBUU,1,-1
           N=N-1
           L=N+1
           DO, K=I,1,-1 
            DO, M=J,1,-1
             L=L-1  
              write(*,*)IC,IR,N,L,EMTBLOCK(L,N) 
            END DO
           END DO
          END DO
         END DO
      goto 90

      IF (LTRANS .AND. (INC2.EQ.1)) THEN !Kai:What do these arugments mean?
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
*       
      WRITE(*,*)"Here comes the contribution from the two-electron 
     : transverse interaction operator in matrixblock2"

*   Diagonal case IC = IR. Start with computing the contribution
*   ENONSYM that does not include the symbolic orbital. For
*   NTYPE = 2 there is only one symbolic orbital 
*
      NVCOEF = 0
      COEFF = 0
      ENONSYM = 0.D0
*     
      CALL RKCO_GG (IC, IR, BREID, 1, 2)
      GOTO 600
*     enonsym is not used right now. To further rise up the speed of
*     caculations, it should be used.    
      DO 27 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))
          IF (LABEL(1,I) .LT. NTYPE(3,IC) .AND. 
     :        LABEL(2,I) .LT. NTYPE(3,IC) .AND.
     :        LABEL(3,I) .LT. NTYPE(3,IC) .AND. 
     :	      LABEL(4,I) .LT. NTYPE(3,IC)) THEN 
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
     :                    LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 6) THEN
              CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                     LABEL(3,I), LABEL(4,I),
     :                     LABEL(5,I), TEGRAL)
            ENDIF 
            WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
            CONTR = COEFF(I)*TEGRAL
            IF (LABEL(6,I) .GT. 0) THEN
              ENONSYM = ENONSYM + CONTR
            ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
!              NCORE = NCORE + 1
              ELSTO = ELSTO + CONTR
            ENDIF
          END IF
        ENDIF
   27 CONTINUE

      DO 28 J = NORB, 1, -1      
        EMTBLOCK(J,J) = EMTBLOCK(J,J) + ENONSYM
   28 CONTINUE

*
*   Diagonal case IC = IR. compute the contribution
*   of the symbolic orbital. 

  600 CONTINUE
      DO 29 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          ITYPE = ABS (LABEL(6,I))
!          IF (LABEL(1,I) .GE. NTYPE(3,IC) .OR. 
!     :        LABEL(2,I) .GE. NTYPE(3,IC) .OR.
!     :        LABEL(3,I) .GE. NTYPE(3,IC) .OR. 
!     :	      LABEL(4,I) .GE. NTYPE(3,IC)) THEN 
            DO 30 J = NORB, 1, -1      
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
     :                      LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 6) THEN
                CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ENDIF 
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
              CONTR = COEFF(I)*TEGRAL
              IF (LABEL(6,I) .GT. 0) THEN
                EMTBLOCK(J,J) = EMTBLOCK(J,J) + CONTR
              ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
!                NCORE = NCORE + 1
                ELSTO = ELSTO + CONTR
              ENDIF          
   30       CONTINUE
!          END IF
        END IF
   29 CONTINUE
       DO 62 J = NORB, 1, -1
         write(*,*)IC,IR,J,EMTBLOCK(J,J)
   62  continue     
*   Here comes the off-diagonal part within the diagonal symbolic block.
      NVCOEF = 0
      CALL RKCO_GG (IC, IC-1, BREID, 1, 2)
!      J = 6 
!      K = 7
!      write(*,*)'kkkkk',J,K
      DO 31 I = 1, NVCOEF
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          VCOEFF = COEFF(I)
          ITYPE = ABS (LABEL(6,I))
          DO 32 J = NORB-1, 1, -1
            DO 33 K = NORB, J+1, -1      
              IF (J .NE. (NORB - 1)) THEN
                DO, M = 1, 3
                  DO, N = M + 1, 4
                    IF (LABEL(M,I) .GE. NTYPE(3,IC) .AND. 
     :                  LABEL(N,I) .GE. NTYPE(3,IC)) THEN
                      IF (LABEL(M,I) .GT. LABEL(N,I)) THEN
                        LABEL(M,I) = K + NTYPE(3,IC) - 1
                        LABEL(N,I) = J + NTYPE(3,IC) - 1
                      ELSE IF (LABEL(M,I) .LT. LABEL(N,I)) THEN
                        LABEL(M,I) = J + NTYPE(3,IC) - 1
                        LABEL(N,I) = K + NTYPE(3,IC) - 1
                      ELSE
                        WRITE(*,*)"STOP 1 IN MATRIXBLOCK2"
                        STOP
                      END IF
                      GOTO 34
                    END IF
                  END DO
                END DO    
              END IF
   34         CONTINUE
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
     :                      LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ELSEIF (ITYPE .EQ. 6) THEN
                CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                       LABEL(3,I), LABEL(4,I),
     :                       LABEL(5,I), TEGRAL)
              ENDIF 
              WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                     NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                     NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                     NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                     NP(LABEL(4,I)),NH(LABEL(4,I)),')'
              CONTR = COEFF(I)*TEGRAL
              IF (LABEL(6,I) .GT. 0) THEN
!                write(*,*)IC,IR,J,K,EMTBLOCK(J,K),CONTR,'111111'
                EMTBLOCK(J,K) = EMTBLOCK(J,K) + CONTR
!                write(*,*)IC,IR,J,K,EMTBLOCK(J,K),CONTR,'222222'
              ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
!                NCORE = NCORE + 1
                ELSTO = ELSTO + CONTR
              ENDIF          

   33       CONTINUE          
   32     CONTINUE
        END IF
   31 CONTINUE   
   
        DO J = 1, NORB-1, 1
           DO K = J+1, NORB, 1 
           write(*,*)IC,IR,J,K,EMTBLOCK(J,K)
          END DO   
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
   90 CONTINUE
      RETURN
      END
