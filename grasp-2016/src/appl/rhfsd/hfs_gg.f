************************************************************************
*                                                                      *
      SUBROUTINE HFS(NAME)
*                                                                      *
*   This routine controls the main sequence of routine calls for the   *
*   calculation  of the  hyperfine structure parameters.               *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, GRACAH1, ISPAR, ITJPO, *
*                        ONEPARTICLEJJ.                                *
*               [HFS92]: MATELT, RINT, RINTHF.                         *
*                                                                      *
*   Written by Per Jonsson and Farid A. Parpia                         *
*                                                                      *
*   Modified by Per Jonsson to evaluate g_j factors                    * 
*                                                                      *
*   Modified by Per Jonsson to read angular data from file             *
*                                                     AUGGUST 2011     *
*                                                                      *
*   Modified by G. Gaigalas by includeing the new spin-angular         *
*   libraries.                                    Vilnius May 2012     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      PARAMETER (KEY = KEYORB)
      POINTER (PNTRIQ,RIQDUMMY)
      INTEGER FFMIN,FFMAX,FF
      CHARACTER*24 NAME
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,AVAIL
*
      DIMENSION TSHELL(NNNW),RINTME(2,NNNW,NNNW),AMELT(2,NNNW,NNNW)
      DIMENSION TSHELL_S(NNNW),IA_S(NNNW)
      DIMENSION RINTGJ(NNNW,NNNW),RINTDGJ(NNNW,NNNW)
      DIMENSION GJMELT(NNNW,NNNW),DGJMELT(NNNW,NNNW)
*
      POINTER (PNTHFC,HFC(2,*))
      POINTER (PNTGJC,GJC(*))
      POINTER (PNTDGJC,DGJC(*))
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF3/EMPAM,RBCM
     :      /DEF9/CVAC,PI
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI
     :      /DEF11/FMTOAU,B1
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FOPARM/ICCUT
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /NSMDAT/HFSI,HFSD,HFSQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
     :      /OPT6/NTC(10)
     :      /ORB2/NCF,NW,PNTRIQ
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-10)
*
*   Allocate storage for local arrays
*
      CALL ALLOC (PNTHFC,2*NVEC,8)
      CALL ALLOC (PNTGJC,NVEC,8)
      CALL ALLOC (PNTDGJC,NVEC,8)
*
*   Initialise
*
      DO I = 1,2
        DO J = 1,NVEC
           HFC(I,J) = 0.0D 00
        END DO 
      END DO
*
      DO J = 1,NVEC
        GJC(J) = 0.0D 00
        DGJC(J) = 0.0D 00
      END DO

*
*   Calculate and save the radial integrals and angular
*   matrix elements for the two multipolarities
*
      DO KT = 1,2
        DO I = 1,NW
          DO J = 1,NW
            IF (KT .EQ. 1) THEN
              RINTME(KT,I,J) = RINTHF (I,J,-2)
              RINTGJ(I,J) = RINTHF(I,J,1)
            ELSE
              RINTME(KT,I,J) = RINT (I,J,-3)
              RINTDGJ(I,J) = RINT(I,J,0)
            ENDIF
            CALL MATELT (I,KT,J,APART,GJPART,DGJPART)
            AMELT(KT,I,J) = APART
            IF (KT.EQ.1) THEN
              GJMELT(I,J) = GJPART
              DGJMELT(I,J) = DGJPART
            ENDIF
          END DO
        END DO 
      END DO
*
*   Set the parity of the one-body operators
*
      IPT = 1
*
*   Check if angular data is available and appropriate
*
      CALL ANGDATA(NAME,AVAIL)
     
      IF (.NOT.AVAIL) THEN
*
*   Sweep through the Hamiltonian matrix to determine the
*   diagonal and off-diagonal hyperfine constants
*
        DO 11 IC = 1,NCF
*
*   Output IC on the screen to show how far the calculation has preceede
*
          IF (MOD(IC,100) .EQ. 0) THEN
            CALL CONVRT (IC,CNUM,LCNUM)
            PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
          ENDIF
*
          DO 10 IR = 1,NCF
*
*   If LFORDR is .TRUE., a `first order' calculation is indicated;
*   only the CSFs with serial numbers exceeding IC are treated specially
*   only diagonal elements are evaluated for the `first order' CSFs
*
            IF (  LFORDR .AND.
     :         (IC .GT. ICCUT) .AND.
     :         (IC .NE. IR) ) GOTO 10
*
            ISPARC = ISPAR (IC)
            ITJPOC = ITJPO (IC)
            ITJPOR = ITJPO (IR)
            IDIFF = ITJPOC - ITJPOR
*
*   Loop over the multipolarities
*
            DO 9 KT = 1,2
*
*   Initialise the accumulator
*
              ELEMNT = 0.0D 00
              ELEMNTGJ = 0.0D 00
              ELEMNTDGJ = 0.0D 00
*
*   Consider 1 case
*               (k)
*   (1) < J || T   || J >    , k = 1,2  and IR >= IC
*
              IF ((IDIFF .EQ. 0) .AND. (IR. GE. IC)) THEN
*
                CALL ONEPARTICLEJJ(KT,IPT,IC,IR,IA,IB,TSHELL)
*
*   Accumulate the contribution from the one-body operators;
*
                IF (IA .NE. 0) THEN
                  IF (IA .EQ. IB) THEN
                    NCOUNT = 0
                    DO 6 IA = 1,NW
                      IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                        NCOUNT = NCOUNT + 1
                        TSHELL_S(NCOUNT) = TSHELL(IA)
                        IA_S(NCOUNT) = IA
                        ELEMNT =   ELEMNT
     :                           + AMELT(KT,IA,IA)
     :                           * RINTME(KT,IA,IA)
     :                           * TSHELL(IA)
                        IF (KT.EQ.1) THEN
                          ELEMNTGJ = ELEMNTGJ 
     :                           + GJMELT(IA,IA)
     :                           * RINTGJ(IA,IA)
     :                           * TSHELL(IA)
                          ELEMNTDGJ = ELEMNTDGJ
     :                           + DGJMELT(IA,IA)
     :                           * RINTDGJ(IA,IA)
     :                           * TSHELL(IA)
                        ENDIF
                      ENDIF
    6               CONTINUE
                    IF (KT.EQ.1) THEN
                      WRITE(129) IC,IR,-NCOUNT
                    ELSE
                      WRITE(129) IC,IR,NCOUNT
                    END IF
                    DO I = 1,NCOUNT
                      LAB = IA_S(I)*KEY + IA_S(I)
                      WRITE(129) TSHELL_S(I),LAB
                    END DO
                  ELSE
                    IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
                      IF (KT.EQ.1) THEN
                        WRITE(129) IC,IR,-1
                      ELSE
                        WRITE(129) IC,IR,1
                      END IF
                      LAB = IA*KEY + IB
                      WRITE(129) TSHELL(1),LAB
                      ELEMNT =   ELEMNT
     :                       + AMELT(KT,IA,IB)
     :                       * RINTME(KT,IA,IB)
     :                       * TSHELL(1)
                      IF (KT.EQ.1) THEN
                        ELEMNTGJ = ELEMNTGJ
     :                       + GJMELT(IA,IB)
     :                       * RINTGJ(IA,IB)
     :                       * TSHELL(1) 
                        ELEMNTDGJ = ELEMNTDGJ
     :                       + DGJMELT(IA,IB)
     :                       * RINTDGJ(IA,IB)
     :                       * TSHELL(1) 
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
*
*   Multiply with the configuration expansion coefficients and add the
*   contributions from the matrix elements to obtain total contributions
*
                DO 8 K = 1,NVEC
                  LOC1 = (K-1)*NCF
                  IF (IR .NE. IC) THEN
                    CONTR =   ELEMNT*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                    CONTRGJ = ELEMNTGJ*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                   CONTRDGJ = ELEMNTDGJ*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                  ELSE
                    CONTR =  ELEMNT*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                    CONTRGJ = ELEMNTGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                    CONTRDGJ = ELEMNTDGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                  ENDIF
*
*   Magnetic dipole and the two operators of the g_j factor
*
                  IF (KT .EQ. 1) THEN
                    HFC(1,K) =   HFC(1,K) + CONTR
                    GJC(K) = GJC(K) + CONTRGJ
                    DGJC(K) = DGJC(K) + CONTRDGJ
*
*   Electric quadrupole
*
                  ELSEIF (KT .EQ. 2) THEN
                    HFC(2,K) =   HFC(2,K) + CONTR
                  ENDIF
    8           CONTINUE
*
              ENDIF
    9       CONTINUE
*
   10     CONTINUE
   11   CONTINUE

      ELSE

*    Read angular data from file and compute matrix elements

        ICOLD = 0
        IROLD = 0
        KTOLD = 0
        DO 
          READ(129,END=999) IC,IR,NCOUNT
          IF (NCOUNT.LT.0) THEN
            KT = 1
            NCOUNT = -NCOUNT
          ELSE
            KT = 2
          END IF
*
*   Initialise the accumulator
*
          IF ((IC.NE.ICOLD).OR.(IR.NE.IROLD).OR.(KT.NE.KTOLD)) THEN
            ISPARC = ISPAR (IC)
            ITJPOC = ITJPO (IC)
            ITJPOR = ITJPO (IR)
            IDIFF = ITJPOC - ITJPOR
            ICOLD = IC
            IROLD = IR
            KTOLD = KT
            ELEMNT = 0.0D 00
            ELEMNTGJ = 0.0D 00
            ELEMNTDGJ = 0.0D 00
          END IF
*
*   Accumulate the contribution from the one-body operators;
*
          DO I = 1,NCOUNT
            READ(129) TSHELL_R,LAB
            IA = LAB/KEY
            IB = MOD(LAB,KEY)
            ELEMNT =   ELEMNT
     :          + AMELT(KT,IA,IB)
     :          * RINTME(KT,IA,IB)
     :          * TSHELL_R
            IF (KT.EQ.1) THEN
              ELEMNTGJ = ELEMNTGJ
     :             + GJMELT(IA,IB)
     :             * RINTGJ(IA,IB)
     :             * TSHELL_R 
              ELEMNTDGJ = ELEMNTDGJ
     :             + DGJMELT(IA,IB)
     :             * RINTDGJ(IA,IB)
     :             * TSHELL_R
            END IF
          END DO
*
*   Multiply with the configuration expansion coefficients and add the
*   contributions from the matrix elements to obtain total contributions
*
          DO K = 1,NVEC
            LOC1 = (K-1)*NCF
            IF (IR .NE. IC) THEN
               CONTR =   ELEMNT*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
               CONTRGJ = ELEMNTGJ*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
               CONTRDGJ = ELEMNTDGJ*2.D0*EVEC(IC+LOC1)*EVEC(IR+LOC1)
            ELSE
               CONTR =  ELEMNT*EVEC(IC+LOC1)*EVEC(IR+LOC1)
               CONTRGJ = ELEMNTGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
               CONTRDGJ = ELEMNTDGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
            ENDIF
*
*   Magnetic dipole and the two operators of the g_j factor
*
            IF (KT .EQ. 1) THEN
              HFC(1,K) =   HFC(1,K) + CONTR
              GJC(K) = GJC(K) + CONTRGJ
              DGJC(K) = DGJC(K) + CONTRDGJ
*
*   Electric quadrupole
*
            ELSEIF (KT .EQ. 2) THEN
              HFC(2,K) =   HFC(2,K) + CONTR
            ENDIF
          END DO 
        END DO
  999   CONTINUE
      END IF

*
*   These are the conversion factors to obtain the hyperfine
*   constants in MHz
*
      AUMHZ = AUCM*CCMS*1.0D-06
      BARNAU = 1.0D-24/RBCM**2
      DNMAU = B1/(2.0D 00*CVAC*EMPAM)
*
      GFAC = AUMHZ*DNMAU*HFSD/HFSI
      HFAC = AUMHZ*2.0D 00*HFSQ*BARNAU
*
*   Output the hyperfine interaction constants
*
      WRITE (29,402)
*
      DO 13 I = 1,NVEC
*
        JJ = IATJPO(I)
*
        IF (JJ.GT.1) THEN
*
          FJ = 0.5D 00*DBLE (JJ-1)
*
          GJA1 =  SQRT (
     :                       1.0D 00
     :                     / (FJ*(FJ+1.0D 00))
     :                    )

          AFA1 =   GFAC*GJA1
            
          BFA1 =   HFAC
     :             * SQRT (
     :                       (FJ*(2.0D 00*FJ-1.0D 00))
     :                     / (   (FJ+1.0D 00)
     :                         * (2.0D 00*FJ+3.0D 00) )
     :                    )
*
*   Output diagonal hfs and g_j factors to file <name>.h or <name>.ch
*
          GJ = CVAC*GJA1*GJC(I)
          DGJ = 0.001160D0*GJA1*DGJC(I)  
          WRITE (29,403) IVEC(I),LABJ(JJ),
     :           LABP((IASPAR(I)+3)/2),
     :           AFA1*HFC(1,I),
     :           BFA1*HFC(2,I),
     :           GJ+DGJ
        END IF
   13 CONTINUE
*
      CALL DALLOC (PNTHFC)
      CALL DALLOC (PNTGJC)
      CALL DALLOC (PNTDGJ)



      RETURN
*
  402 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity ',8X,'A (MHz)',13X,'B (MHz)',13X,'g_J'/)
  403 FORMAT (1X,1I3,5X,2A4,1P,3D20.10)
*
      END
