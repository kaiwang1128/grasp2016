************************************************************************
*                                                                      *
      SUBROUTINE HFS (host)
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
*   Modified by G. Gaigalas by includeing the new spin-angular         *
*   libraries.                                      Vilnius May 2012   *
*                                                                      *
*   Updated by Per Jonsson,    July 2015                               *         
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
c
cbieron include 'parameters.def'
c
      include 'parameters.def'

      POINTER (PNTRIQ,RIQDUMMY)
      INTEGER FFMIN,FFMAX,FF
      CHARACTER*11 CNUM
      CHARACTER*4 JLBL,LABJ,LABP
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
*
      DIMENSION TSHELL(NNNW),RINTME(2,NNNW,NNNW),AMELT(2,NNNW,NNNW)
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
      include 'mpif.h'
      integer  myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr 

      integer lenhost
c     character*80  host
c     CHARACTER host*(MPI_MAX_PROCESSOR_NAME)
      character*(*) host

      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond
      character msg*80

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

      DO J = 1,NVEC
        GJC(J) = 0.0D 00                                                                                                                                                                    
        DGJC(J) = 0.0D 00
      END DO

*
*   Calculate and save the radial integrals and angular
*   matrix elements for the two multipolarities
*
      DO 5 KT = 1,2
        DO 4 I = 1,NW
          DO 3 J = 1,NW
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
    3     CONTINUE
    4   CONTINUE
    5 CONTINUE

*
*   Set the parity of the one-body operators
*
      IPT = 1
*
*   Sweep through the Hamiltonian matrix to determine the
*   diagonal and off-diagonal hyperfine constants
*
cb MPI
      DO 11 IC = myid+1, NCF, nprocs
*
*   Output IC on the screen to show how far the calculation has preceede
*
C         IF (IC .LE. 3*nprocs .OR. IC .GE. NCF-3*nprocs
C     &       .OR. MOD (IC-1,100*nprocs) .EQ. myid) THEN
      if (mod(IC,100).eq.0) then
        CALL CONVRT (IC,CNUM,LCNUM)
        PRINT *, 'myid=',myid, ' Column '//CNUM(1:LCNUM)//' complete'
       endif
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
                  DO 6 IA = 1,NW
                    IF (ABS (TSHELL(IA)) .GT. CUTOFF) THEN
                      ELEMNT =   ELEMNT
     :                         + AMELT(KT,IA,IA)
     :                         * RINTME(KT,IA,IA)
     :                         * TSHELL(IA)
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
    6             CONTINUE
                ELSE
                  IF (ABS (TSHELL(1)) .GT. CUTOFF) THEN
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
                  CONTR = ELEMNT*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                  CONTRGJ = ELEMNTGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                  CONTRDGJ = ELEMNTDGJ*EVEC(IC+LOC1)*EVEC(IR+LOC1)
                ENDIF
*
*   Magnetic dipole
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
    8         CONTINUE
*
            ENDIF
    9     CONTINUE
*
   10   CONTINUE
   11 CONTINUE

      CALL MPI_ALLREDUCE (MPI_IN_PLACE,HFC,2*NVEC,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_SUM, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE (MPI_IN_PLACE,GJC,NVEC,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_SUM, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE (MPI_IN_PLACE,DGJC,NVEC,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_SUM, MPI_COMM_WORLD, ierr)



cb 
cb from here to END only node0 works
      if (myid .eq. 0) then !node0
*
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  start: ' // 
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime 
      msgLength = len_trim (msg)
      write(*,*)
      print *, 'Only node0 works:'
      print *, msg(1:msgLength)      

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
      WRITE (24,402)
*
      DO 13 I = 1,NVEC
*
        JJ = IATJPO(I)
*
        IF (JJ.GT.1) THEN
*
          FJ = 0.5D 00*DBLE (JJ-1)
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
          WRITE (24,403) IVEC(I),LABJ(JJ),
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
cb
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  stop : ' // 
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime 
      msgLength = len_trim (msg)
      print *, msg(1:msgLength)      

      endif ! node0

      RETURN
*
  402 FORMAT (//' Interaction constants:'//
     :' Level1  J Parity ',8X,'A (MHz)',13X,'B (MHz)',13X,'g_J'/)
  403 FORMAT (1X,1I3,5X,2A4,1P,3D20.10)       
*
      call MPI_Get_processor_name (host, lenhost, ierr)
      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      END
