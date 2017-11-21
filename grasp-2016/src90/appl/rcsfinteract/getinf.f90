!***********************************************************************
!                                                                      *
      SUBROUTINE GETINF 
!                                                                      *
!   Interactively determines data governing the generation of MCP co-  *
!   efficients.                                                        *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:32:22   1/ 3/07  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!GG      USE DEF1_C 
      USE DEF_C 
      USE DEFAULT_C 
      USE FOPARM_C 
!GG      USE MCPA_C 
      USE MCP_C 
!GG      USE MCPB_C 
!GG      USE ORB2_C 
      USE ORB_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE getyn_I 
      USE convrt_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LENTH, K 
      LOGICAL :: YES 
      CHARACTER :: CNUM*20 
!-----------------------------------------------
!
!
!   Determine the physical effects specifications
!
      IF (NDEF /= 0) THEN 
         WRITE (6, *) 'Generate MCP coefficients only for' 
         WRITE (6, *) ' diagonal matrix elements? (This is' 
         WRITE (6, *) ' appropriate to an (E)AL calculation):' 
         YES = GETYN() 
      ELSE 
         YES = .FALSE. 
      ENDIF 
      IF (YES) THEN 
         DIAG = .TRUE. 
         LFORDR = .FALSE. 
         ICCUT = 0 
      ELSE 
         DIAG = .FALSE. 
         IF (NDEF /= 0) THEN 
            WRITE (6, *) 'Treat contributions of some CSFs' 
            WRITE (6, *) ' as first-order perturbations?' 
            YES = GETYN() 
         ELSE 
            YES = .FALSE. 
         ENDIF 
         IF (YES) THEN 
            LFORDR = .TRUE. 
            WRITE (6, *) 'The contribution of CSFs' 
            WRITE (6, *) ' 1 -- ICCUT will be treated' 
            WRITE (6, *) ' variationally; the remainder' 
            WRITE (6, *) ' perturbatively; enter ICCUT:' 
    1       CONTINUE 
            READ (5, *) ICCUT 
            IF (ICCUT<=1 .OR. ICCUT>=NCF) THEN 
               CALL CONVRT (NCF, CNUM, LENTH) 
               WRITE (6, *) 'GETINF: ICCUT must be greater than 1' 
               WRITE (6, *) ' and less than '//CNUM(1:LENTH)//';' 
               WRITE (6, *) ' please reenter ICCUT:' 
               GO TO 1 
            ENDIF 
         ELSE 
            LFORDR = .FALSE. 
            ICCUT = 0 
         ENDIF 
      ENDIF 
!
!   Add the second and third records to the file headers; the
!   first is written in SETRES
!
      DO K = 30, 32 + KMAX 
         WRITE (K) NELEC, NCF, NW 
         WRITE (K) DIAG, ICCUT, LFORDR 
      END DO 
!
      RETURN  
      END SUBROUTINE GETINF 
