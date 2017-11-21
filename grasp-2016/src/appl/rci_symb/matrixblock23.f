************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK23(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :           NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
*                                                                      *
*   This subroutine calls onescalar and computes one electron          *
*   matrix elements when IC and IR are of type 2 and 3                 *         
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
      INTEGER NTYPE(6,NCF),NUMB
      DIMENSION TSHELL(NNNW),EMTBLOCK(MAXSPAN,MAXSPAN)
      DIMENSION BREITBLOCK(MAXSPAN,MAXSPAN),TMPBLOCK(MAXSPAN,MAXSPAN)

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
      BREITBLOCK = 0.D0
      TMPBLOCK = 0.D0
      TSHELL = 0.D0
*   Call onescalar to       
      WRITE(*,*) 'ONESCALAR In matrixblock23, IC,IR',IC,IR

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
        N = NORBUL * NORBUU + 1
        DO J = NORBUL, 1, -1
          DO K = NORBUU, 1, -1
            N = N -1
            DO M = NORBL, 1, -1              
              IF ((K+ NTYPE(5,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN
                NCOEC = NCOEC + 1
!              IA = I + NTYPE(3,IC) - 1
                IB = J + NTYPE(3,IC) - 1
                CALL IABINT(IA,IB,TEGRAL)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*TCOEFF
                WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                       NP(IB),NH(IB),')'
                IF (LNMS) THEN
                  CALL KEINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*ATWINV*TCOEFF
                ENDIF
                IF (LVP) THEN
                  CALL VPINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*TCOEFF
                ENDIF
              END IF
              IF ((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN
                NCOEC = NCOEC + 1
                IB = K + NTYPE(5,IC) - 1
                CALL IABINT(IA,IB,TEGRAL)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*TCOEFF
                WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                       NP(IB),NH(IB),')'
                IF (LNMS) THEN
                  CALL KEINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*ATWINV*TCOEFF
                ENDIF
                IF (LVP) THEN
                  CALL VPINT(IA,IB,TEGRAL,TEGRAL)
                  EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*TCOEFF
                ENDIF
              END IF
            END DO   
          END DO   
        END DO
      END IF

      DO M = NORBL, 1, -1
        DO N=NORBUL * NORBUU, 1, -1
           WRITE(*,*)M,N,EMTBLOCK(M,N)
        END DO
      END DO


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
*   Otherwise, The Rk has two components. MRANK=1-9

*   For 4s 3p and 2s 3p
*             0.333333333 R1(2s 3p ,3p 4s )
*            -1.000000000 R0(2s 3p ,4s 3p )
*            -2.000000000 R0(1s 2s ,1s 4s )
*             1.000000000 R0(1s 1s ,2s 4s )
*     The absolute values of four components, with opposite in sign, 
*     are the same with the absolute values of four components 
*     of 4s 4p and 2s 4p.
*     This comment is wrong. 
*     For MRANK=1-9, none orbital of CSF in type 3 is   
*     the same with one orbital of CSF in type 2.
*            -2.000000000 R0(1s 2s ,1s 4s )
*             1.000000000 R0(1s 1s ,2s 4s )
*     The values (-2.000000000) of above components should not be 
*     included in the calculations.
*     For MRANK=1-14, the orbitals (2s 3p ,3p 4s )
*     should be changed with different 
*     symbolic orbitals.
*      

      WRITE(*,*) 'RKCO_GG CORD In matrixblock23, IC,IR',IC,IR

        N = NORBUL * NORBUU + 1
        M = NORBL              
      NVCOEF = 0
      CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
      DO  I = 1, NVCOEF
        VCOEFF = COEFF(I)
            CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)

           WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
      END DO

*
      write(*,*)'111111'
      
        MMIN=MIN(NTYPE(3,IC),NTYPE(3,IR))
        write(*,*)MMIN
      write(*,*)'22222'

      DO 7 I = 1, NVCOEF
        VCOEFF = COEFF(I)
        MRANK = 0
        IF (DABS (VCOEFF) .GT. CUTOFF) THEN
          N = NORBUL * NORBUU + 1
          DO J = NORBUL, 1, -1
            DO K = NORBUU, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
!                WRITE(*,*)'MRANK',MRANK,'N',N,'M',M
                IF(MRANK .GE. 9 .OR.
     :           (N .EQ. NORBUL * NORBUU .AND. M .EQ. NORBL)
     :           .OR.((M + NTYPE(3,IR)).EQ.(J + NTYPE(3,IC)))
     :           .OR.((M + NTYPE(3,IR)).EQ.(K + NTYPE(5,IC))))THEN
                  VCOEFF = COEFF(I)
                ELSE
                  VCOEFF = 0.D0               
                END IF                  ! negative contributions?  

                IF (MRANK .EQ. 1) THEN
                  LABEL(1,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 2) THEN
                  LABEL(2,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 3) THEN
                  LABEL(3,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 4) THEN
                  LABEL(4,I) = K + NTYPE(5,IC) - 1
                ELSE IF (MRANK .EQ. 5) THEN
                  LABEL(1,I) = J + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 6) THEN
                  LABEL(2,I) = J + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 7) THEN
                  LABEL(3,I) = J + NTYPE(3,IC) - 1
                ELSE IF (MRANK .EQ. 8) THEN
                  LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                ELSE IF (MRANK .EQ. 9) THEN
                  LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
!                  write(*,*)LABEL(2,I),LABEL(3,I),LABEL(4,I)
!                  write(*,*)M,J,K
!                  write(*,*)NTYPE(3,IR),NTYPE(3,IC),NTYPE(5,IC)
!                  write(*,*)NP(LABEL(2,I)),NP(LABEL(3,I)),NP(LABEL(4,I))
!                  write(*,*)NH(LABEL(2,I)),NH(LABEL(3,I)),NH(LABEL(4,I))
                ELSE IF (MRANK .EQ. 10) THEN
                  LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                ELSE IF (MRANK .EQ. 11) THEN
                  LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                ELSE IF (MRANK .EQ. 12) THEN
                  LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                ELSE IF (MRANK .EQ. 13) THEN
                  LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                ELSE IF (MRANK .EQ. 14) THEN
                  LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                  LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                  LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
               END IF
            IF (LSMS) THEN
              IF (LABEL(5,I) .EQ. 1) THEN
                CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                EMTBLOCK(M,N) = EMTBLOCK(M,N) - TGRL1*TGRL2*ATWINV*VCOEFF
!                WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N),'11'
              ENDIF
            ENDIF

            CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                   LABEL(3,I), LABEL(4,I),
     :                   LABEL(5,I), TEGRAL)

            EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*VCOEFF
!            WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N),'22'

            WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),VCOEFF,'R',
     :                   LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                IF(N .EQ. NORBUL * NORBUU .AND. M .EQ. NORBL)then 
                  NUMB = 0
                  DO, MTYPE=1,4
                    IF(LABEL(MTYPE,I).GE.MMIN)THEN
                      NUMB = NUMB + 1
                    END IF
                  END DO
                  IF(NUMB .EQ. 1)THEN
                    IF( LABEL(1,I) .EQ. NTYPE(6,IC)) THEN
                      MRANK = 1
                    ELSE IF(LABEL(2,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 2
                    ELSE IF(LABEL(3,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 3
                    ELSE IF(LABEL(4,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 4
                    ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                      MRANK = 5
                    ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                      MRANK = 6
                    ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                      MRANK = 7
                    ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                      MRANK = 8
                    END IF
                  ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IR))THEN
                    IF(LABEL(3,I) .EQ. NTYPE(4,IC).AND. 
     :                 LABEL(4,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 9
                    ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC).AND. 
     :                     LABEL(3,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 10
                    END IF  
                  ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IR))THEN
                    IF(LABEL(2,I) .EQ. NTYPE(4,IC).AND. 
     :                 LABEL(4,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 11
                    ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC).AND. 
     :                      LABEL(2,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 12
                    END IF  
                  ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IR))THEN
                    IF(LABEL(2,I) .EQ. NTYPE(4,IC).AND. 
     :                 LABEL(3,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 13
                    ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC).AND. 
     :                      LABEL(2,I) .EQ. NTYPE(6,IC))THEN
                      MRANK = 14
                    END IF  
                  ELSE 
                      WRITE(*,*)'STOP 1 IN MATRIXBOLOCK23'
                      STOP  
                  END IF                                    
                END IF
              END DO   
            END DO   
          END DO
        END IF
    7 CONTINUE



         IF(NTYPE(4,IC).eq.NTYPE(4,IR).or.
     :    NTYPE(6,IC).eq.NTYPE(4,IR))THEN
            WRITE(*,*)NTYPE(4,IR),NTYPE(4,IC),NTYPE(6,IC)
            DO M = 1,NORBL,  1
             WRITE(*,*)(EMTBLOCK(M,N),N=1,NORBUL*NORBUU,1)
             write(*,*)
            END DO
            DO, i=0,(NORBUL * NORBUU)/NORBL-1
              TMPBLOCK(1:NORBL,1+i*NORBL:(i+1)*NORBL)=
     :        TRANSPOSE(EMTBLOCK(1:NORBL,1+i*NORBL:(i+1)*NORBL))
            END DO
            EMTBLOCK=0.D0
            EMTBLOCK(1:NORBL,1:NORBUL * NORBUU)
     :       =TMPBLOCK(1:NORBL,1:NORBUL * NORBUU)

         END IF
            DO M = 1, NORBL, 1
              WRITE(*,*)(EMTBLOCK(M,N),N=1,NORBUL*NORBUU,1)
              write(*,*)
            END DO


      WRITE(*,*) 'RKCO_GG BREID In matrixblock23, IC,IR',IC,IR




    
      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
*     
       IF(NTYPE(4,IC).EQ.NTYPE(4,IR).OR.
     :    NTYPE(6,IC).EQ.NTYPE(4,IR))THEN

         NVCOEF = 0
       
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
*
        DO 8 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            N = NORBUL * NORBUU + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
                 IF((J + NTYPE(3,IC)).EQ.(M + NTYPE(3,IR))
     :            .OR.(K + NTYPE(5,IC)).EQ.(M + NTYPE(3,IR)))THEN
                  IF (MRANK .EQ. 1) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 2) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 3) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 4) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 5) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 6) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 7) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 8) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 9) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 10) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 11) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 12) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 13) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 14) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 15) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 16) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 17) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 18) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 19) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 20) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 21) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 22) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 23) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 24) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 25) THEN
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 26) THEN
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 27) THEN
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 28) THEN
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 29) THEN
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 30) THEN
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 31) THEN
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 32) THEN
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE
                    IF(MRANK .NE. 0)THEN
                      WRITE(*,*)'STOP 10 IN MATRXIBLOCK23'
                      STOP  
                    END IF                                      
                  END IF
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
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'

                  CONTR = COEFF(I)*TEGRAL
                  IF (LABEL(6,I) .GT. 0) THEN
                     EMTBLOCK(M,N) = EMTBLOCK(M,N) + CONTR
                  ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     NCORE = NCORE + 1
                     ELSTO = ELSTO + CONTR
                  END IF
                  IF(N .EQ. NORBUL * NORBUU .AND. M .EQ. NORBL)then
                    NUMB = 0
                    DO, MTYPE=1,4
                      IF(LABEL(MTYPE,I).GE.MMIN)THEN
                        NUMB = NUMB + 1
                      END IF
                    END DO
                    IF(NUMB .EQ. 3)THEN
                      IF(LABEL(1,I).EQ.NTYPE(4,IR))THEN
                        IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 1
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 2
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 3
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 4
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 5
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 6
                        ELSE
                          WRITE(*,*)'STOP 3 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN                              
                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 7
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 8
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 9
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 10
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 11
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 12
                        ELSE
                          WRITE(*,*)'STOP 4 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN                              
                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 13
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 14
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 15
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 16
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 17
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 18
                        ELSE
                          WRITE(*,*)'STOP 4 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN                              
                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 19
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 20
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 21
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 22
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 23
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 24
                        ELSE
                          WRITE(*,*)'STOP 5 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE
                        WRITE(*,*)'STOP 6 IN MATRXIBLOCK23'
                        STOP                                                     
                      END IF
                    ELSE IF(NUMB .EQ. 1)THEN
                      IF(NTYPE(4,IR).EQ.NTYPE(4,IC))THEN
                        IF(LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 25
                        ELSE IF(LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 26
                        ELSE IF(LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 27
                        ELSE IF(LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 28
                        ELSE
                          WRITE(*,*)'STOP 7 IN MATRXIBLOCK23'
                          STOP
                        END IF                      
                      ELSE IF(NTYPE(4,IR).EQ.NTYPE(6,IC))THEN
                        IF(LABEL(1,I).EQ.NTYPE(4,IC))THEN
                          MRANK = 29
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC))THEN
                          MRANK = 30
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC))THEN
                          MRANK = 31
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC))THEN
                          MRANK = 32
                        ELSE
                          WRITE(*,*)'STOP 8 IN MATRXIBLOCK23'
                          STOP
                        END IF                      
                      ELSE
                        WRITE(*,*)'STOP 9 IN MATRXIBLOCK23'
                        STOP                      
                      END IF
                    ELSE
                      WRITE(*,*)'STOP 2 IN MATRXIBLOCK23'
                      STOP
                    END IF
                  END IF
                 END IF
             END DO   
            END DO   
           END DO

          END IF
    8   CONTINUE
        DO N=NORBUL * NORBUU, 1, -1
          DO M = NORBL, 1, -1
           WRITE(*,*)M,N,EMTBLOCK(M,N)
        END DO
      END DO
    
        NVCOEF = 0
        CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
*

        DO 9 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            N = NORBUL * NORBUU 
            M = NORBL
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
!                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
!     :                   ,'R',LABEL(5,I),'(',
!     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  IF(N .EQ. NORBUL * NORBUU .AND.M.EQ.NORBL)then
                    NUMB = 0
                    DO, MTYPE=1,4
                      IF(LABEL(MTYPE,I).GE.MMIN)THEN
                        NUMB = NUMB + 1
                      END IF
                    END DO
                    IF(NUMB .EQ. 3)THEN
                      IF(LABEL(1,I).EQ.NTYPE(4,IR)-1)THEN
                        LABEL(1,I)= NTYPE(4,IR)
                        IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 1
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 2
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 3
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 4
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 5
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 6
                        ELSE
                          WRITE(*,*)'STOP 13 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR)-1)THEN                              
                        LABEL(2,I)= NTYPE(4,IR)

                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 7
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 8
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 9
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 10
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 11
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 12
                        ELSE
                          WRITE(*,*)'STOP 14 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR)-1)THEN                              
                        LABEL(3,I)= NTYPE(4,IR)

                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 13
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 14
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 15
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 16
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 17
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 18
                        ELSE
                          WRITE(*,*)'STOP 20 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR)-1)THEN                              
                        LABEL(4,I)= NTYPE(4,IR)
                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 19
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 20
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 21
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 22
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 23
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 24
                        ELSE
                          WRITE(*,*)'STOP 15 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE
                        WRITE(*,*)'STOP 16 IN MATRXIBLOCK23'
                        STOP                                                     
                      END IF
                    ELSE
                      WRITE(*,*)'STOP 12 IN MATRXIBLOCK23'
                      STOP
                    END IF
                  END IF
                              
            N = NORBUL * NORBUU + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
                 IF((J + NTYPE(3,IC)).NE.(M + NTYPE(3,IR))
     :            .AND.(K + NTYPE(5,IC)).NE.(M + NTYPE(3,IR)))THEN
                  IF (MRANK .EQ. 1) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 2) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 3) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 4) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 5) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 6) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 7) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 8) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 9) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 10) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 11) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 12) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 13) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 14) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 15) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 16) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 17) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 18) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 19) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 20) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 21) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 22) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 23) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 24) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 25) THEN
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 26) THEN
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 27) THEN
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 28) THEN
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 29) THEN
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 30) THEN
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 31) THEN
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 32) THEN
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE
                    IF(MRANK .NE. 0)THEN
                      WRITE(*,*)'STOP 110 IN MATRXIBLOCK23'
                      STOP  
                    END IF                                      
                  END IF
!                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
!     :                   ,'R',LABEL(5,I),'(',
!     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!                  WRITE(*,*)'2222'
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
!                  WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N)
!                  write(*,*)'111'
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  CONTR = COEFF(I)*TEGRAL
                  IF (LABEL(6,I) .GT. 0) THEN
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + CONTR
!                    BREITBLOCK(M,N) = BREITBLOCK(M,N) + CONTR
!               WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),BREITBLOCK(M,N)
                  ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     NCORE = NCORE + 1
                     ELSTO = ELSTO + CONTR
                  END IF

                 END IF

             END DO   
            END DO   
           END DO

          END IF
    9   CONTINUE    
       ELSE IF(NTYPE(4,IC).NE.NTYPE(4,IR).AND.
     :    NTYPE(6,IC).NE.NTYPE(4,IR))THEN

         NVCOEF = 0
       
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
        DO 10 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            N = NORBUL * NORBUU 
            M = NORBL
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
!                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
!     :                   ,'R',LABEL(5,I),'(',
!     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  IF(N .EQ. NORBUL * NORBUU .AND.M.EQ.NORBL)then
                    NUMB = 0
                    DO, MTYPE=1,4
                      IF(LABEL(MTYPE,I).GE.MMIN)THEN
                        NUMB = NUMB + 1
                      END IF
                    END DO
                    IF(NUMB .EQ. 3)THEN
                      IF(LABEL(1,I).EQ.NTYPE(4,IR))THEN
                        LABEL(1,I)= NTYPE(4,IR)
                        IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 1
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 2
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 3
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 4
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 5
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 6
                        ELSE
                          WRITE(*,*)'STOP 23 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN                              
                        LABEL(2,I)= NTYPE(4,IR)

                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 7
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 8
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 9
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 10
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 11
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 12
                        ELSE
                          WRITE(*,*)'STOP 24 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN                              
                        LABEL(3,I)= NTYPE(4,IR)

                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 13
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 14
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 15
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(4,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 16
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 17
                        ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 18
                        ELSE
                          WRITE(*,*)'STOP 30 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN                              
                        LABEL(4,I)= NTYPE(4,IR)
                        IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                     LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 19
                        ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 20
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 21
                        ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(3,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 22
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(1,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 23
                        ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                          LABEL(2,I).EQ.NTYPE(6,IC))THEN
                          MRANK = 24
                        ELSE
                          WRITE(*,*)'STOP 25 IN MATRXIBLOCK23'
                          STOP                               
                        END IF
                      ELSE
                        WRITE(*,*)'STOP 26 IN MATRXIBLOCK23'
                        STOP                                                     
                      END IF
                    ELSE
                      WRITE(*,*)'STOP 22 IN MATRXIBLOCK23'
                      STOP
                    END IF
                  END IF
                              
            N = NORBUL * NORBUU + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
!                 IF((J + NTYPE(3,IC)).NE.(M + NTYPE(3,IR))
!     :            .AND.(K + NTYPE(5,IC)).NE.(M + NTYPE(3,IR)))THEN
                  IF (MRANK .EQ. 1) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 2) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 3) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 4) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 5) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 6) THEN
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 7) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 8) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 9) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 10) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 11) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 12) THEN
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 13) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 14) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 15) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 16) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 17) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 18) THEN
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 19) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 20) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 21) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 22) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 23) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 24) THEN
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 25) THEN
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 26) THEN
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 27) THEN
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 28) THEN
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 29) THEN
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 30) THEN
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 31) THEN
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE IF (MRANK .EQ. 32) THEN
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                  ELSE
                    IF(MRANK .NE. 0)THEN
                      WRITE(*,*)'STOP 210 IN MATRXIBLOCK23'
                      STOP  
                    END IF                                      
                  END IF
!                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
!     :                   ,'R',LABEL(5,I),'(',
!     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
!                  WRITE(*,*)'2222'
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
!                  WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N)
!                  write(*,*)'111'
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                  CONTR = COEFF(I)*TEGRAL
                  IF (LABEL(6,I) .GT. 0) THEN
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) + CONTR
!                    BREITBLOCK(M,N) = BREITBLOCK(M,N) + CONTR
!               WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),BREITBLOCK(M,N)
                  ELSE
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     NCORE = NCORE + 1
                     ELSTO = ELSTO + CONTR
                  END IF

!                 END IF

             END DO   
            END DO   
           END DO

          END IF
   10   CONTINUE    

       ELSE
         WRITE(*,*)'STOP 55 IN MATRIXBLOCK23'
         STOP
       END IF

       
        DO N=NORBUL * NORBUU, 1, -1
          DO M = NORBL, 1, -1
           WRITE(*,*)M,N,BREITBLOCK(M,N),
     :      EMTBLOCK(M,N),BREITBLOCK(M,N)+EMTBLOCK(M,N)
        END DO
      END DO

*
        IBUG1 = 0
* 
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
      ENDIF
!
  592 FORMAT (2I5,F15.9,1X,1A2,1I2,1A2,1A1,I2,1A2,1A1)
  593 FORMAT (2I5,F15.9,1X,1A,1I1,1A1,1I2,1A2,1I2,1A2,
     :        1A1,1I2,1A2,1I2,1A2,1A1) 
  594 FORMAT (5I5,F15.9,1X,1A,1I1,1A1,1I2,1A2,1I2,1A2,
     :        1A1,1I2,1A2,1I2,1A2,1A1) 
     
      RETURN
      END
