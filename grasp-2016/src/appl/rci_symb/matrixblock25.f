************************************************************************
*                                                                      *
      SUBROUTINE MATRIXBLOCK25(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
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
      INTEGER MRANKS
      WRITE(*,*) 'In matrixblock25, IC,IR',IC,IR

      ATWINV = 1.D0/EMN

      IBUG1 = 0

*   Set all matrix elements to zero

      EMTBLOCK = 0.D0
      BREITBLOCK = 0.D0
      TMPBLOCK = 0.D0
      TSHELL = 0.D0
      NORBUU = 0
      NORBUL = 0
      NORBL = 0
*   Call onescalar to       

      WRITE(*,*) 'ONESCALAR In matrixblock25, IC,IR',IC,IR
      CALL ONESCALAR(IC,IR,IA,IB,TSHELL) 
      WRITE(*,592) IC,IR,DBLE(TSHELL(1)),'I(',NP(IA),NH(IA),',',
     :                   NP(IB),NH(IB),')'        

*   Ensure that the indices are in `canonical' order
      IF (IA .GT. IB) THEN
        ISWAP = IB
        IB = IA
        IA = ISWAP
      ENDIF
      WRITE(*,592) IC,IR,DBLE(TSHELL(1)),'I(',NP(IA),NH(IA),',',
     :                   NP(IB),NH(IB),')'        
*    In matrixblock 25, for symbolic CSF of Type 2, the total number of
*    symbolic CSF (and symbolic orbitals) is defined by NORBL. For symbolic 
*    CSF of Type 5, there are two symbolic orbitals (the same kind). the total  
*    number of the first （second）kind of symbolic orbitali defined by NORBUL. 
*    The total number of symbolic CSF of Type 5 is defined by 
*    N = NTYPE(2,IC) = NORBUL

!      NORBUU = NTYPE(6,IC) - NTYPE(5,IC) + 1
      NORBUL = NTYPE(4,IC) - NTYPE(3,IC) + 1
      NORBL = NTYPE(4,IR) - NTYPE(3,IR) + 1
!
!      write(*,*)NTYPE(2,IC),NTYPE(2,IR)
!      write(*,*)NORBUL,NTYPE(4,IC),NTYPE(3,IC),NTYPE(4,IC)-NTYPE(3,IC)+1
!      write(*,*)NORBUU,NTYPE(6,IC),NTYPE(5,IC),NTYPE(6,IC)-NTYPE(5,IC)+1 
!      write(*,*)NORBL,NTYPE(4,IR),NTYPE(3,IR),NTYPE(4,IR)-NTYPE(3,IR)+1

      TCOEFF = DBLE(TSHELL(1))

*      1s ( 1)  2p ( 2)  9p ( 1)  
*          1/2        2      3/2  
*                      3/2      1-
*    
*      1s ( 1)  2p ( 1)  9p ( 2)    
*            1/2      3/2        0  
*                          1      1-
*       -1.000000000 I( 2p , 9p )
*    
*
*     NORBL =7, the total number of symbolic CSF is 7 [3-9]p.
*     NORBUL = NTYPE(2,IC) 7, The total number of symbolic CSF of Type 5 
*     Among the 7*7 cases, one electron operator for some specific cases 
*     exists, such as, for the case---1s ( 1)  2p ( 2)  9p ( 1)  and 
*     1s ( 1)  2p ( 1)  9p ( 2), but does not exist for cases, such as,
*     for the case---1s ( 1)  2p ( 2)  8p ( 1) and 
*     1s ( 1)  2p ( 1)  9p ( 2). Therefore, when 
*     the second orbital ([3-9]p) of symbolic CSF of Type 2 equals to the
*     the first (second) orbital ([3-9]p) of symbolic CSF of Type 5,  
*     one electron operator exists. This condition is defined by 
*     IF ((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN.
*     Among the 7*7 cases, one electron operator for 7 cases. 
*     When the second orbital of symbolic CSF of Type 2 equals to the
*     the first (senond) orbital ([3-9]p) or the second 
*     orbital ([3-11]s) of symbolic CSF of Type 4, IB should be changed 
*     by IB = M + NTYPE(3,IR) - 1 
*     cases. IA is always the first orbital of symbolic CSF of Type 2.

      IF (DABS(TCOEFF) .GT. CUTOFF) THEN
        N = NTYPE(2,IC) + 1
        DO J = NORBUL, 1, -1
          N = N -1
          DO M = NORBL, 1, -1              
            IF ((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1)) THEN
              NCOEC = NCOEC + 1
              IB = M + NTYPE(3,IR) - 1
              CALL IABINT(IA,IB,TEGRAL)
              EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*TCOEFF
              WRITE(*,592) IC,IR,TCOEFF,'I(',NP(IA),NH(IA),',',
     :                     NP(IB),NH(IB),')'
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
        N = NTYPE(2,IC) + 1
        DO J = NORBUL, 1, -1
          N = N -1
          DO M = NORBL, 1, -1              
            write(*,*)J,K,N,M,EMTBLOCK(M,N)
          END DO
        END DO
      END IF


      IBUG1 = 0
*
*   Accumulate the contributions from the two-electron
*   Coulomb operator and the mass polarisation; the latter
*   is computed first because the orbital indices may be
*   permuted by RKINTC
*

*   Here comes the off-diagonal part within the diagonal symbolic block.

*   The Rk has always two or four components
*   For 1s ( 1)  2p ( 1) 10s ( 1) 11s ( 1) and 1s ( 2)  2p ( 1) 11s ( 1)
*   , The Rk has four components
*             -1.118033989 R0( 1s  1s , 1s 10s )
*             -1.118033989 R0( 1s  2p ,10s  2p )
*             -1.118033989 R0( 1s 11s ,10s 11s )
*              1.118033989 R0( 1s 10s ,11s 11s )
*   For 1s ( 1)  2p ( 1) 10s ( 1) 11s ( 1) and 1s ( 2)  2p ( 1) 10s ( 1)
*   , The Rk has four components
*              1.118033989 R0( 1s  1s , 1s 11s )
*              1.118033989 R0( 1s  2p ,11s  2p )
*             -1.118033989 R0( 1s 11s ,10s 11s )
*              1.118033989 R0( 1s 10s ,11s 11s )
*     The absolute values of two components between 
*     1s ( 1)  2p ( 1) 10s ( 1) 11s ( 1) and 1s ( 2)  2p ( 1) 10s ( 1), 
*     with opposite in sign, 
*     are the same with the absolute values of two components between 
*     1s ( 1)  2p ( 1) 10s ( 1) 11s ( 1) and 1s ( 2)  2p ( 1) 11s ( 1) 

*   For 1s ( 1)  2p ( 1) 10s ( 1) 11s ( 1) and 1s ( 2)  2p ( 1) 9s ( 1)
*   , The Rk has two components
*             -1.118033989 R0( 1s 11s ,10s 11s )
*              1.118033989 R0( 1s 10s ,11s 11s )
*   If one orbital of CSF in type 4 is the same with one orbital of CSF 
*   in type 2, The Rk has four components. 
*   Otherwise, The Rk has two components. MRANK=1-9

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

      WRITE(*,*) 'RKCO_GG CORD In matrixblock25, IC,IR',IC,IR

      N = NTYPE(2,IC) + 1
      M = NORBL              
      IF(NTYPE(3,IC).EQ.NTYPE(3,IR))THEN
       NVCOEF = 0
       COEFF = 0       
       CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
       DO  I = 1, NVCOEF
         VCOEFF = COEFF(I)
         CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :          LABEL(3,I), LABEL(4,I),
     :          LABEL(5,I), TEGRAL)
       
         WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                NP(LABEL(4,I)),NH(LABEL(4,I)),')'
       END DO

       MMIN = MIN(NTYPE(3,IC),NTYPE(3,IR))
       DO 7 I = 1, NVCOEF
         VCOEFF = COEFF(I)
         MRANK = 0
         IF (DABS (VCOEFF) .GT. CUTOFF) THEN
           N = NTYPE(2,IC) + 1
           DO J = NORBUL, 1, -1
             N = N -1
             DO M = NORBL, 1, -1              
               IF(N .EQ. NTYPE(2,IC) .AND. M .EQ. NORBL)then 
                 NUMB = 0
                 DO, MTYPE=1,4
                   IF(LABEL(MTYPE,I).GE.MMIN)THEN
                     NUMB = NUMB + 1
                   END IF
                 END DO
                 IF(NUMB .EQ. 1)THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 1
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 2
                   ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 3
                   ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 4
                   END IF
                 ELSE IF(LABEL(1,I) .LT. MMIN)THEN
                     MRANK = 5
                 ELSE IF(LABEL(2,I) .LT. MMIN)THEN
                     MRANK = 6
                 ELSE IF(LABEL(3,I) .LT. MMIN)THEN
                     MRANK = 7
                 ELSE IF(LABEL(4,I) .LT. MMIN)THEN
                     MRANK = 8
                 ELSE 
                     WRITE(*,*)'STOP 2 IN MATRIXBOLOCK25'
                     STOP  
                 END IF                                    
               END IF
               IF((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1))THEN
                 SELECT CASE (MRANK)
                 CASE (0)
                   GOTO 50
                 CASE (1)
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                 CASE (2)
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                CASE (3)
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (4)
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (5)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (6)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (7)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (8)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 END SELECT              
   50            CONTINUE
                 
                 IF (LSMS) THEN
                   IF (LABEL(5,I) .EQ. 1) THEN
                    CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                    CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) - 
     :               TGRL1*TGRL2*ATWINV*VCOEFF
                   ENDIF
                 ENDIF
                 WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),VCOEFF,'R',
     :                        LABEL(5,I),'(',
     :                        NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                        NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                        NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                        NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                 
                 CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                        LABEL(3,I), LABEL(4,I),
     :                        LABEL(5,I), TEGRAL)
                 
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*VCOEFF
!                WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N),'22'
               END IF                    
             END DO   
           END DO
         END IF
    7  CONTINUE      
       DO M = 1, NORBL, 1
         WRITE(*,*)(EMTBLOCK(M,N),N=1,NTYPE(2,IC),1)
         write(*,*)
       END DO
       NVCOEF = 0
       COEFF = 0       
       CALL RKCO_GG (IC, IR-1, CORD, INCOR, 1)
       DO  I = 1, NVCOEF
         VCOEFF = COEFF(I)
         CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :          LABEL(3,I), LABEL(4,I),
     :          LABEL(5,I), TEGRAL)
       
         WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                NP(LABEL(4,I)),NH(LABEL(4,I)),')'
       END DO

       MMIN = MIN(NTYPE(3,IC),NTYPE(3,IR))
       DO 71 I = 1, NVCOEF
         VCOEFF = COEFF(I)
         MRANK = 0
         IF (DABS (VCOEFF) .GT. CUTOFF) THEN
           N = NTYPE(2,IC) + 1
           DO J = NORBUL, 1, -1
             N = N -1
             DO M = NORBL, 1, -1              
               IF(N .EQ. NTYPE(2,IC) .AND. M .EQ. NORBL)then 
                 NUMB = 0
                 DO, MTYPE=1,4
                   IF(LABEL(MTYPE,I).GE.MMIN)THEN
                     NUMB = NUMB + 1
                   END IF
                 END DO
                 IF(NUMB .EQ. 1)THEN
                   WRITE(*,*)'STOP 3 IN MATRIXBOLOCK25'
                   STOP  
                 ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IR)-1)THEN
                   IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 1
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 2
                   ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 3
                   END IF
                 ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IR)-1)THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 4
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 5
                   ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 6
                   END IF
                 ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IR)-1)THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 7
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 8
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 9
                   END IF
                 ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IR)-1)THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 10
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 11
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 12
                   END IF
                 ELSE 
                     WRITE(*,*)'STOP 4 IN MATRIXBOLOCK25'
                     STOP  
                 END IF                                    
               END IF
               IF((J+ NTYPE(3,IC) - 1) .NE. (M + NTYPE(3,IR) - 1))THEN
                 SELECT CASE (MRANK)
                 CASE (0)
                   GOTO 60
                 CASE (1)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (2)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                CASE (3)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (4)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (5)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (6)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (7)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                 CASE (8)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (9)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (10)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                 CASE (11)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (12)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 END SELECT              
   60            CONTINUE
                 
                 IF (LSMS) THEN
                   IF (LABEL(5,I) .EQ. 1) THEN
                    CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                    CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) - 
     :               TGRL1*TGRL2*ATWINV*VCOEFF
                   ENDIF
                 ENDIF
                 WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),VCOEFF,'R',
     :                        LABEL(5,I),'(',
     :                        NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                        NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                        NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                        NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                 
                 CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                        LABEL(3,I), LABEL(4,I),
     :                        LABEL(5,I), TEGRAL)
                 
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*VCOEFF
!                WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N),'22'
               END IF                    
             END DO   
           END DO
         END IF
   71  CONTINUE   
       DO M = 1, NORBL, 1
         WRITE(*,*)(EMTBLOCK(M,N),N=1,NTYPE(2,IC),1)
         write(*,*)
       END DO
       
      ELSE IF(NTYPE(3,IC).NE.NTYPE(3,IR))THEN    
       NVCOEF = 0
       COEFF = 0       
       CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
       DO  I = 1, NVCOEF
         VCOEFF = COEFF(I)
         CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :          LABEL(3,I), LABEL(4,I),
     :          LABEL(5,I), TEGRAL)
       
         WRITE(*,593) IC,IR,VCOEFF,'R',LABEL(5,I),'(',
     :                NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                NP(LABEL(4,I)),NH(LABEL(4,I)),')'
       END DO

       MMIN = MIN(NTYPE(3,IC),NTYPE(3,IR))
       DO 72 I = 1, NVCOEF
         VCOEFF = COEFF(I)
         MRANK = 0
         IF (DABS (VCOEFF) .GT. CUTOFF) THEN
           N = NTYPE(2,IC) + 1
           DO J = NORBUL, 1, -1
             N = N -1
             DO M = NORBL, 1, -1              
               IF(N .EQ. NTYPE(2,IC) .AND. M .EQ. NORBL)then 
                 NUMB = 0
                 DO, MTYPE=1,4
                   IF(LABEL(MTYPE,I).GE.MMIN)THEN
                     NUMB = NUMB + 1
                   END IF
                 END DO
                 IF(NUMB .EQ. 1)THEN
                   WRITE(*,*)'STOP 3 IN MATRIXBOLOCK25'
                   STOP  
                 ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IR))THEN
                   IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 1
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 2
                   ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 3
                   END IF
                 ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IR))THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 4
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 5
                   ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 6
                   END IF
                 ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IR))THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 7
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 8
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 9
                   END IF
                 ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IR))THEN
                   IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 10
                   ELSE IF(LABEL(1,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 11
                   ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC)
     :             .AND.LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                     MRANK = 12
                   END IF
                 ELSE 
                     WRITE(*,*)'STOP 4 IN MATRIXBOLOCK25'
                     STOP  
                 END IF                                    
               END IF
               IF((J+ NTYPE(3,IC) - 1) .NE. (M + NTYPE(3,IR) - 1))THEN
                 SELECT CASE (MRANK)
                 CASE (0)
                   GOTO 70
                 CASE (1)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (2)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                CASE (3)
                   LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (4)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (5)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (6)
                   LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (7)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                 CASE (8)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (9)
                   LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(4,I) = J + NTYPE(3,IC) - 1
                 CASE (10)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                 CASE (11)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(1,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 CASE (12)
                   LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                   LABEL(2,I) = J + NTYPE(3,IC) - 1
                   LABEL(3,I) = J + NTYPE(3,IC) - 1
                 END SELECT              
   70            CONTINUE
                 
                 IF (LSMS) THEN
                   IF (LABEL(5,I) .EQ. 1) THEN
                    CALL VINT (LABEL(1,I), LABEL(3,I), TGRL1)
                    CALL VINT (LABEL(2,I), LABEL(4,I), TGRL2)
                    EMTBLOCK(M,N) = EMTBLOCK(M,N) - 
     :               TGRL1*TGRL2*ATWINV*VCOEFF
                   ENDIF
                 ENDIF
                 WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),VCOEFF,'R',
     :                        LABEL(5,I),'(',
     :                        NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                        NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                        NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                        NP(LABEL(4,I)),NH(LABEL(4,I)),')'
                 
                 CALL RKINTC (LABEL(1,I), LABEL(2,I),
     :                        LABEL(3,I), LABEL(4,I),
     :                        LABEL(5,I), TEGRAL)
                 
                 EMTBLOCK(M,N) = EMTBLOCK(M,N) + TEGRAL*VCOEFF
!                WRITE(*,*)IC,IR,int(MRANK),int(N),int(M),EMTBLOCK(M,N),'22'
               END IF                    
             END DO   
           END DO
         END IF
   72  CONTINUE   
       DO M = 1, NORBL, 1
         WRITE(*,*)(EMTBLOCK(M,N),N=1,NTYPE(2,IC),1)
         write(*,*)
       END DO       
      END IF 
           





      WRITE(*,*) 'RKCO_GG BREID In matrixblock25, IC,IR',IC,IR
    
      IBUG1 = 0

      IF (LTRANS .AND. (INC2.EQ.1)) THEN 
*
*   Accumulate the contribution from the two-electron
*   transverse interaction operator
*
*     
       IF(NTYPE(4,IC).EQ.NTYPE(4,IR))THEN

        NVCOEF = 0
        COEFF =0 
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
*
        DO 8 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            N = NTYPE(2,IC) + 1
            IF (ITYPE .EQ. 1) THEN
               CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 2) THEN
               CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 3) THEN
               CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 4) THEN
               CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 5) THEN
               CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 6) THEN
               CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ENDIF
            NUMB = 0
            DO, MTYPE=1,4
              IF(LABEL(MTYPE,I).GE.MMIN)THEN
                NUMB = NUMB + 1
              END IF
            END DO
            IF(NUMB .EQ. 1)THEN
              IF(LABEL(1,I) .EQ. NTYPE(4,IC))THEN
                MRANK = 1
              ELSE IF(LABEL(2,I) .EQ. NTYPE(4,IC))THEN
                MRANK = 2
              ELSE IF(LABEL(3,I) .EQ. NTYPE(4,IC))THEN
                MRANK = 3
              ELSE IF(LABEL(4,I) .EQ. NTYPE(4,IC))THEN
                MRANK = 4
              END IF
            ELSE IF(LABEL(1,I) .LT. MMIN)THEN
                MRANK = 5
            ELSE IF(LABEL(2,I) .LT. MMIN)THEN
                MRANK = 6
            ELSE IF(LABEL(3,I) .LT. MMIN)THEN
                MRANK = 7
            ELSE IF(LABEL(4,I) .LT. MMIN)THEN
                MRANK = 8
            ELSE 
                WRITE(*,*)'STOP 5 IN MATRIXBOLOCK25'
                STOP  
            END IF                                    
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
            DO J = NORBUL, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
               IF((J+ NTYPE(3,IC) - 1) .EQ. (M + NTYPE(3,IR) - 1))THEN
                  SELECT CASE (MRANK)
                  CASE (:0)
                     WRITE(*,*)'STOP 6 IN MATRXIBLOCK25'
                     STOP
                  CASE (1)
                    LABEL(1,I) = J + NTYPE(3,IC) - 1
                  CASE (2)
                    LABEL(2,I) = J + NTYPE(3,IC) - 1
                  CASE (3)
                    LABEL(3,I) = J + NTYPE(3,IC) - 1
                  CASE (4)
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  CASE (5)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  CASE (6)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  CASE (7)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1
                    LABEL(4,I) = J + NTYPE(3,IC) - 1
                  CASE (8)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1
                    LABEL(3,I) = J + NTYPE(3,IC) - 1
                  END SELECT              
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
                END IF
              END DO   
            END DO
          END IF
    8   CONTINUE
            N = NTYPE(2,IC) + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, J, -1
              N = N -1
              DO M = NORBL, 1, -1              
                WRITE(*,*)M,N,J,K,M,EMTBLOCK(M,N)
             END DO   
            END DO   
           END DO

        NVCOEF = 0
        CALL RKCO_GG (IC, IR-1, BREID, 1, 2)
*

        DO 10 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            IF (ITYPE .EQ. 1) THEN

               CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 2) THEN
               CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 3) THEN
               CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 4) THEN
               CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 5) THEN
               CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 6) THEN
               CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ENDIF 
                  WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
     :                   ,'R',LABEL(5,I),'(',
     :                   NP(LABEL(1,I)),NH(LABEL(1,I)),
     :                   NP(LABEL(2,I)),NH(LABEL(2,I)),',',
     :                   NP(LABEL(3,I)),NH(LABEL(3,I)),
     :                   NP(LABEL(4,I)),NH(LABEL(4,I)),')'
            write(*,*)'111111'
            NUMB = 0
            DO, MTYPE=1,4
              IF(LABEL(MTYPE,I).GE.MMIN)THEN
                NUMB = NUMB + 1
              END IF
            END DO
            IF(NUMB .EQ. 1)THEN
              WRITE(*,*)'STOP 7 IN MATRXIBLOCK25'
              STOP
            ELSE IF(NUMB .EQ. 3)THEN
              IF(LABEL(1,I).EQ.NTYPE(4,IR)-1)THEN
                IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 1
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 2
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 3
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 4
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 5
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 6
                ELSE
                  WRITE(*,*)'STOP 8 IN MATRXIBLOCK25'
                  STOP                               
                END IF
              ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR)-1)THEN                              
                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 7
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 8
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 9
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 10
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 11
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 12
                ELSE
                  WRITE(*,*)'STOP 9 IN MATRXIBLOCK29'
                  STOP                               
                END IF
              ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR)-1)THEN                              
                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 13
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 14
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 15
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 16
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 17
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 18
                ELSE
                  WRITE(*,*)'STOP 10 IN MATRXIBLOCK25'
                  STOP                               
                END IF
              ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR)-1)THEN                              
                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 19
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 20
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 21
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 22
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 23
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 24
                ELSE
                  WRITE(*,*)'STOP 11 IN MATRXIBLOCK25'
                  STOP                               
                END IF
              END IF
            ELSE
              WRITE(*,*)'STOP 12 IN MATRXIBLOCK25'
              STOP
            END IF
                              
            N = NTYPE(2,IC) + 1
            DO J = NORBUL, 1, -1
              N = N -1
              DO M = NORBL, 1, -1              
                 IF((J + NTYPE(3,IC)).NE.(M + NTYPE(3,IR)))THEN
                  SELECT CASE (MRANK)
                  CASE (:0)
                    WRITE(*,*)'STOP 13 IN MATRXIBLOCK25'
                    STOP
                  CASE (1)                 
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (2)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (3)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (4)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (5)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (6)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (7)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (8)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (9)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (10)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (11)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (12) 
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (13)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (14)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (15)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (16)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (17)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (18)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (19)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (20)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (21)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (22)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (23)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (24)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (25:)
                    WRITE(*,*)'STOP 14 IN MATRXIBLOCK25'
                    STOP  
                  END SELECT              
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

          END IF
   10   CONTINUE  
            N = NTYPE(2,IC) + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, J, -1
              N = N -1
              DO M = NORBL, 1, -1              
                WRITE(*,*)M,N,J,K,M,EMTBLOCK(M,N)
             END DO   
            END DO   
           END DO
           
                 
       ELSE IF(NTYPE(4,IC).NE.NTYPE(4,IR).AND.
     :    NTYPE(6,IC).NE.NTYPE(4,IR))THEN
       goto 600

         NVCOEF = 0
       
        CALL RKCO_GG (IC, IR, BREID, 1, 2)
        DO 11 I = 1, NVCOEF
          IF (DABS (COEFF(I)) .GT. CUTOFF) THEN
            MRANK = 0
            NMCBP = NMCBP + 1
            ITYPE = ABS (LABEL(6,I))
            IF (ITYPE .EQ. 1) THEN

               CALL BRINT1 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 2) THEN
               CALL BRINT2 (LABEL(1,I), LABEL(2,I), 
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 3) THEN
               CALL BRINT3 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 4) THEN
               CALL BRINT4 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 5) THEN
               CALL BRINT5 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ELSEIF (ITYPE .EQ. 6) THEN
               CALL BRINT6 (LABEL(1,I), LABEL(2,I),
     :                      LABEL(3,I), LABEL(4,I),
     :                      LABEL(5,I), TEGRAL)
            ENDIF 
!            WRITE(*,594)IC,IR,int(MRANK),int(N),int(M),COEFF(I)
!     :             ,'R',LABEL(5,I),'(',
!     :             NP(LABEL(1,I)),NH(LABEL(1,I)),
!     :             NP(LABEL(2,I)),NH(LABEL(2,I)),',',
!     :             NP(LABEL(3,I)),NH(LABEL(3,I)),
!     :             NP(LABEL(4,I)),NH(LABEL(4,I)),')'

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
     :             LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 1
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 2
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 3
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 4
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 5
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 6
                ELSE
                  WRITE(*,*)'STOP 33 IN MATRXIBLOCK24'
                  STOP                               
                END IF
              ELSE IF(LABEL(2,I).EQ.NTYPE(4,IR))THEN                              
                LABEL(2,I)= NTYPE(4,IR)

                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 7
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 8
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 9
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 10
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 11
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 12
                ELSE
                  WRITE(*,*)'STOP 34 IN MATRXIBLOCK24'
                  STOP                               
                END IF
              ELSE IF(LABEL(3,I).EQ.NTYPE(4,IR))THEN                              
                LABEL(3,I)= NTYPE(4,IR)

                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 13
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 14
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 15
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(4,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 16
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 17
                ELSE IF(LABEL(4,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 18
                ELSE
                  WRITE(*,*)'STOP 35 IN MATRXIBLOCK24'
                  STOP                               
                END IF
              ELSE IF(LABEL(4,I).EQ.NTYPE(4,IR))THEN                              
                LABEL(4,I)= NTYPE(4,IR)
                IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :             LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 19
                ELSE IF(LABEL(1,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 20
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 21
                ELSE IF(LABEL(2,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(3,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 22
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(1,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 23
                ELSE IF(LABEL(3,I).EQ.NTYPE(4,IC).AND.
     :                  LABEL(2,I).EQ.NTYPE(6,IC))THEN
                  MRANK = 24
                ELSE
                  WRITE(*,*)'STOP 36 IN MATRXIBLOCK24'
                  STOP                               
                END IF
              ELSE
                WRITE(*,*)'STOP 37 IN MATRXIBLOCK24'
                STOP                                                     
              END IF
            ELSE
              WRITE(*,*)'STOP 38 IN MATRXIBLOCK24'
              STOP
            END IF
            N = NTYPE(2,IC) + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, J, -1
              N = N -1
              DO M = NORBL, 1, -1              
                  SELECT CASE (MRANK)
                  CASE (:0)
                    WRITE(*,*)'STOP 39 IN MATRXIBLOCK24'
                    STOP
                  CASE (1)                 
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (2)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (3)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (4)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (5)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (6)
                    LABEL(1,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (7)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (8)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (9)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (10)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (11)
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (12) 
                    LABEL(2,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (13)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (14)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (15)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (16)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(4,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (17)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (18)
                    LABEL(3,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(4,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (19)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (20)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(1,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (21)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (22)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(2,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(3,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (23)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(1,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (24)
                    LABEL(4,I) = M + NTYPE(3,IR) - 1                                  
                    LABEL(3,I) = J + NTYPE(3,IC) - 1                                  
                    LABEL(2,I) = K + NTYPE(5,IC) - 1                                  
                  CASE (25:)
                    WRITE(*,*)'STOP 40 IN MATRXIBLOCK24'
                    STOP  
                  END SELECT              
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
   11   CONTINUE    
            N = NTYPE(2,IC) + 1
            DO J = NORBUL, 1, -1
             DO K = NORBUU, J, -1
              N = N -1
              DO M = NORBL, 1, -1              
                WRITE(*,*)M,N,J,K,M,EMTBLOCK(M,N)
             END DO   
            END DO   
           END DO

       ELSE
         WRITE(*,*)'STOP 41 IN MATRIXBLOCK24'
         STOP
       END IF

       

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
  600 continue
     
      RETURN
      END
