************************************************************************
*                                                                      *
      SUBROUTINE MATELT (I1,K,I2,APART,GJPART,DGJPART)
*                                                                      *
*   This  routine computes  the angular part  of the  reduced matrix   *
*   elements of the magnetic and electric multipole operators as       *
*   well as for th operators associated with the g_j factor            *
*                                                                      *
*                                                                      *
*   Calls to: [LIB92]: CLRX                                            *
*                                                                      *
*   Written by Per Jonsson                Last revision: july 2015     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
c
cbieron include 'parameters.def'
c
      include 'parameters.def'
c
c      PARAMETER (NNNW = 120)
c
*
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
*
*   Set KAP1 and KAP2
*
      KAP1 = NAK(I1)
*
      IF (MOD (K,2) .EQ. 1) THEN
         KAP2 = -NAK(I2)
      ELSE
         KAP2 =  NAK(I2)
      ENDIF
*
*   Determine the l quantum numbers
*
      IF (KAP1 .GT. 0) THEN
         L1 =  KAP1
      ELSE
         L1 = -KAP1-1
      ENDIF
*
      IF (KAP2 .GT. 0) THEN
         L2 =  KAP2
      ELSE
         L2 = -KAP2-1
      ENDIF
*
      IF (MOD (L1+K+L2,2) .EQ. 0) THEN
*
*   Parity selection rule satisfied
*
*   Determine the phase factor
*
         IF (MOD (KAP1+1,2) .EQ. 0) THEN
            FASE =  1.0D 00
         ELSE
            FASE = -1.0D 00
         ENDIF
*
*   The other factor is \sqrt (2 j_2 + 1); since j = | \kappa | - 1/2,
*   we have 2 j_2 + 1 = 2 | \kappa |; the factor \sqrt (2 j_2 + 1)
*   has been accounted for in MCT
*
         OVLFAC = FASE * SQRT (DBLE (2*ABS (KAP2)))
*
         IF (MOD (K,2) .EQ. 1) THEN
*
*   These are for the magnetic multipole moments  and the two operators
*   of the g_j factor
*
            APART =   (KAP1+NAK(I2))
     :              * CLRX (KAP1,K,KAP2)
     :              * OVLFAC
*
            GJPART = APART
            DGJPART = -(KAP1+NAK(I2)-1.0D 00)
     :              * CLRX (KAP1,K,KAP2)
     :              * OVLFAC
         ELSE
*
*   These are for the electric multipole moments
*
            APART =   CLRX (KAP1,K,KAP2)
     :              * OVLFAC
*
         ENDIF
*
      ELSE
*
         APART = 0.0D 00
         GJPART = 0.0D 00
         DGJPART = 0.0D 00
*
      ENDIF
*
      RETURN
      END
