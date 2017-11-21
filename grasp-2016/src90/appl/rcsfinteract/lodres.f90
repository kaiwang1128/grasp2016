!***********************************************************************
!                                                                      *
      SUBROUTINE LODRES(IERR) 
!                                                                      *
!   Loads and checks certain data from the  .res  files.               *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETINF, LENGTH, LODRES, OPENFL.       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 07 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:32:22   1/ 3/07  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!GG      USE DEF1_C 
      USE DEF_C 
!CGG      USE MCPA_C 
      USE MCP_C 
      USE FOPARM_C 
!GG      USE MCPB_C 
!GG      USE ORB2_C 
      USE ORB_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: IERR 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NELECT, NCFT, NWT, K, ICCUTT 
      LOGICAL :: DIAGT, LFRDRT 
!-----------------------------------------------
!
!
      IERR = 0 
!
!   The second and third record in each file should be identical
!
      READ (30) NELECT, NCFT, NWT 
      IF (NELECT /= NELEC) THEN 
         WRITE (6, *) 'LODRES: The number of electrons does not' 
         WRITE (6, *) ' match that from the Configuration Symmetry' 
         WRITE (6, *) ' List File;' 
         IERR = IERR + 1 
      ENDIF 
      IF (NCFT /= NCF) THEN 
         WRITE (6, *) 'LODRES: The number of CSFs does not' 
         WRITE (6, *) ' match that from the Configuration' 
         WRITE (6, *) ' Symmetry List File;' 
         IERR = IERR + 1 
      ENDIF 
      IF (NWT /= NW) THEN 
         WRITE (6, *) 'LODRES: The number of subshells does not' 
         WRITE (6, *) ' match that from the Configuration' 
         WRITE (6, *) ' Symmetry List File;' 
         IERR = IERR + 1 
      ENDIF 
!
      READ (30) DIAG, ICCUT, LFORDR 
      IF (LFORDR) THEN 
         IF (ICCUT<1 .OR. ICCUT>NCF) THEN 
            WRITE (6, *) 'LODRES: The first-order calculation' 
            WRITE (6, *) ' determined from the GENMCP REStart' 
            WRITE (6, *) ' files is not appropriate to the' 
            WRITE (6, *) ' Configuration Symmetry List File;' 
            IERR = IERR + 1 
         ENDIF 
      ENDIF 
!
      DO K = 31, 32 + KMAX 
         READ (K) NELECT, NCFT, NWT 
         IF (NELEC/=NELECT .OR. NCF/=NCFT .OR. NWT/=NW) THEN 
            WRITE (6, *) 'LODRES: Inconsistent electron cloud data' 
            WRITE (6, *) ' detected on comparing GENMCP REStart' 
            WRITE (6, *) ' Files;' 
            IERR = IERR + 1 
         ENDIF 
         READ (K) DIAGT, ICCUTT, LFRDRT 
         IF (.NOT.((DIAG .NEQV. DIAGT) .OR. ICCUT/=ICCUTT .OR. (LFORDR .NEQV. &
            LFRDRT))) CYCLE  
         WRITE (6, *) 'LODRES: Inconsistent calculation type' 
         WRITE (6, *) ' detected on comparing GENMCP REStart' 
         WRITE (6, *) ' Files;' 
         IERR = IERR + 1 
      END DO 
!
      RETURN  
      END SUBROUTINE LODRES 
