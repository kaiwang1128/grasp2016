************************************************************************
*                                                                      *
      SUBROUTINE SETHAM (myid, nprocs, jblock, ELSTO,ICSTRT, nelmntt
     &                   , atwinv,slf_en)
*                                                                      *
*   Sets up the Hamiltonian matrix and determines the average energy.  *
*                                                                      *
*   Serial I/O moved out; able to run on single/multiple processors    *
*   For this purpose a new common /setham_to_genmat2/ is created       *
*                                                                      *
*   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
*               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
*   Block version by Xinghong He          Last revision: 15 Jun 1998   *
*   Shift diagonal elements by Per J                      March 2007   *
*                                                                      *
************************************************************************
*
      use symexpand_mod  
      IMPLICIT REAL*8          (A-H, O-Z)

      integer*8 nelmnt,nelmntt,nelmnttmp

      include 'parameters.def'
      POINTER (PINDT1,INDT1DUMMY)
      POINTER (PINDT2,INDT2DUMMY)
      POINTER (PINDT3,INDT3DUMMY)
      POINTER (PINDT4,INDT4DUMMY)
      POINTER (PINDT5,INDT5DUMMY)
      POINTER (PINDT6,INDT6DUMMY)
      POINTER (PVALT1,VALT1DUMMY)
      POINTER (PVALT2,VALT2DUMMY)
      POINTER (PVALT3,VALT3DUMMY)
      POINTER (PVALT4,VALT4DUMMY)
      POINTER (PVALT5,VALT5DUMMY)
      POINTER (PVALT6,VALT6DUMMY)
      POINTER (PCOEIL,COEILDUMMY)
      POINTER (PCOEVL,COEVLDUMMY)
      POINTER (PCTEIL,CTEILDUMMY)
      POINTER (PCTEVL,CTEVLDUMMY)
      POINTER (PNEVAL,EVALDUMMY)
      POINTER (PINDKE,INDKEDUMMY)
      POINTER (PVALKE,VALKEDUMMY)
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      POINTER (PNTJQS,JQSDUMMY)
      POINTER (PNJCUP,JCUPDUMMY)
      POINTER (PVINIL,VINILDUMMY)
      POINTER (PVINVL,VINVLDUMMY)
      POINTER (PINDVP,INDVPDUMMY)
      POINTER (PVALVP,VALVPDUMMY)
      POINTER (PNTRPF,RPFDUMMY)
      POINTER (PNTRQF,RQFDUMMY)
      POINTER (PNIVEC,NVECMXDUMMY)

      POINTER (PCTEVLRK,VALTEIRK(1))                                  
      POINTER (PCTEILRK, INDTEIRK(1))

      INTEGER, DIMENSION(:), ALLOCATABLE :: MAP
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NTYPE
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IROWSYM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: EMTSYM
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: EMTBLOCK

      LOGICAL FIRST,FRSTCO,FRSTCT,FRSTKI,FRSTVI,FRSTVP,
     :        LDBPA,lshift,
     :        LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES,GETYN
      CHARACTER*2 NH
*
      EXTERNAL BREID,CORD
*
      DIMENSION TSHELL(NNNW),SLF_EN(1)

      POINTER (PNTEMT,EMT(1))
      POINTER (PNIROW,IROW(1))
*
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
*
      COMMON/BCORE/ICORE(NNNW)
     :      /BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /COEILS/NDCOEA,NCOEI,PCOEIL,PCOEVL,FRSTCO
     :      /CTEILS/NDCTEA,NCTEI,PCTEIL,PCTEVL,FRSTCT
     :      /DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEBUGA/LDBPA(5)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /FOPARM/ICCUT
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI
     :      /NCDIST/ZDIST(NNNP)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /STAT/PNTJQS,PNJCUP
     :      /STOR/KEEP(2,2)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /VINLST/NDVIN,NVINTI,PVINIL,PVINVL,FRSTVI
     :      /VPILST/NDVPA,NVPI,PINDVP,PVALVP,FRSTVP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /CTEILSRK/PCTEILRK,PCTEVLRK
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF
     :      /DEF10/AUCM,AUEV,CCMS,FASI,FBSI

* Per common for shift

      COMMON/hamshiftj/nshiftj(100),nasfshift(100,100),
     :       asfenergy(100,100),lshift

* Per end

*     ...For pre-run
      POINTER (PNEVEC1,EVEC1(1))
      COMMON/EIGVEC1/PNEVEC1

*
      PARAMETER (KEY = KEYORB)
*
*   Matrix elements smaller than CUTOFF are not accumulated
*
      PARAMETER (CUTOFF = 1.0D-20)

      COMMON/setham_to_genmat2/CUTOFFtmp,
     &  NCOEItmp, NCOECtmp, NCTEItmp, NCTECtmp, NTPItmp(6), NMCBPtmp, 
     &  NCOREtmp, NVPItmp, NKEItmp, NVINTItmp, NELMNTtmp, ncftmp
*
!-----------------------------------------------------------------------
      PRINT *, 'Calling setham ...'


!      WRITE(*,*) 'For now we define the map by hand'
!      WRITE(*,*) 'Number of CSFs in added lis', NCF
!      WRITE(*,*) 'Give the number of CSFs in the original small list'
!      READ(*,*) NCFSMALL
       NCFSMALL=ncsfsymperblock(jblock)
      ALLOCATE(MAP(NCFSMALL))

      DO I = 1,NCFSMALL
!         WRITE(*,*) 'Give position of CSF',I,'in the added list'
!         READ(*,*) MAP(I)
         MAP(I)=map1(i,jblock)
      END DO

!  Call findtype in order to define the type of CSF and largest value
!  for the orbitals of the symbolic CSF. Note that we define the type
!  for all CSFs in the added list although we only will need the
!  information for the CSFs defined by the mapping.
!
      ALLOCATE(NTYPE(6,NCF))

      CALL FINDTYPE(NTYPE,MAP,NCFSMALL,NCFGEN,NCFTOT,MAXSPAN,NORBGEN)

      WRITE(*,*) 'In setham'
      DO K = 1,NCFSMALL
        J = MAP(K)
        WRITE(*,*) 'CSF in added list',J
        WRITE(*,*) NTYPE(:,J)                                                                                                                                                         
      END DO

      WRITE(*,*) 'NCFGEN',NCFGEN
      WRITE(*,*) 'NCFTOT',NCFTOT
      WRITE(*,*) 'MAXSPAN',MAXSPAN
      WRITE(*,*) 'NORBGEN',NORBGEN

!  Allocate three arrays
!  ESYM(MAXSPAN,MAXSPAN)   matrix elements between CSFs spanned by
!                          two symbolic CSFs
!  EMTSYM(NCFTOT,MAXSPAN)  holds in sparse form all matrix elements
!                          involving a symbolic CSF
!  IROWSYM(NCFTOT,MAXSPAN) holds corresponding rows     

      ALLOCATE(IROWSYM(NCFTOT,MAXSPAN))
      ALLOCATE(EMTSYM(NCFTOT,MAXSPAN))
      ALLOCATE(EMTBLOCK(MAXSPAN,MAXSPAN))
!     
**************************
*
*     CSF(i) = symbolic CSF that spans N(i) CSFs      
*     <CSF   |H|CSF(i)> generates an 1 x N(i) array of elements
*     <CSF(i)|H|CSF(i)> generates an upper triangular N(i) x N(i) array
*     <CSF(i)|H|CSF(j)> generates an N(i) x N(j) array of elements      
*
*     Allocate an array of some largest dimension Nmax x Nmax as computed by
*     FINDTYPE. Store the elements in this array
*
*     Aoocate a double precision and an integer array with the largest 
*     dimension NCSFtot x Nmax for storing the matrix elements and the
*     row positions.
*     NCSFtot is computed by FINDTYPE. In reality we can take NCSFtot
*     smaller as the matrix will be sparse.     
*     10 000 000 x 400 = 32 Gb is very much nore than we need.
*
*     As the computation of the interactions between the symbolic CSFs
*     proceeds fill the arrays above with non-zero elements. We can then 
*     write to rci.res in the same format as before
*
*      
      nelmnt = nelmntt
      
      ATWINV = 1.D0/EMN
      !Per Bug fix 30/4-2013

      IPRERUN = 0   ! Pre-run has been disabled
*
*   Allocate storage to arrays in COMMON/BUFFER/; these are
*   used for the Coulomb and transverse two-electron integrals
*
      CALL ALCBUF (1)

*     ...Locals
      CALL alloc (pntemt, ncf, 8) !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!
      CALL alloc (pnirow, ncf, 4) !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!
*
      INC1 = 1
      INC2 = 1
*
*   Initialisations for contributions from the Dirac-Coulomb
*   operator
*
      KT  = 0
      IPT = 1
*
      INCOR = 1

      NCOEC = 0
*
      NCTEC   = 0

      IF (LTRANS) THEN

*        ...Initialisations for transverse interaction correction
        DO 2 I = 1, NW
          ICORE(I) = 0
          DO J = 1, NCF
            IF (ICHOP (I,J) .LE. 0) GOTO 2
          ENDDO
          ICORE(I) = 1
    2   CONTINUE

        NMCBP = 0
        NCORE = 0
      ENDIF

! Loop over rows of the Hamiltonian matrix - distributed

      ICTOT = 0
      DO 10 ICC = icstrt, NCFSMALL, nprocs

        IC = MAP(ICC) 

        NELC = 0    ! counter - Number of non-zeros of this row

! Loop over columns of the current row

        IRTOT = 0
        DO 85 IRR = 1, ICC

          IR = MAP(IRR) 

          IF (LFORDR .AND. (IR .GT. ICCUT)) THEN
            IF (IR.NE.IC) CYCLE
          END IF             

          EMTBLOCK = 0.D0     ! accumulates various contributions to H 

          CALL ONESCALAR11(IC,IR,NCOEC,EMTBLOCK(1,1))
          CALL RKCO_GG11(IC,IR,INCOR,NCTEC,INC2,NMCBP,NCORE,
     :                   ELSTO,EMTBLOCK(1,1))
          IF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.1)) THEN
            CALL MATRIXBLOCK12(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :                         NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN)
!          ELSEIF ((NTYPE(1,IC).EQ.2).AND.(NTYPE(1,IR).EQ.2)) THEN
!             IF (IC.EQ.IR) THEN   ! Diagonal block of type 2
!              CALL MATRIXBLOCK2(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!             ELSE IF (IC.NE.IR) THEN   
!              CALL MATRIXBLOCK22(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!            END IF
!          ELSEIF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.1)) THEN
!              CALL MATRIXBLOCK13(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!          ELSEIF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.2)) THEN
!              CALL MATRIXBLOCK23(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!          ELSEIF ((NTYPE(1,IC).EQ.3).AND.(NTYPE(1,IR).EQ.3)) THEN
!             IF (IC.EQ.IR) THEN   ! Diagonal block of type 3
!              CALL MATRIXBLOCK3(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!             END IF
!          ELSEIF ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.1)) THEN
!              CALL MATRIXBLOCK14(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!          ELSEIF ((NTYPE(1,IC).EQ.4).AND.(NTYPE(1,IR).EQ.2)) THEN
!              CALL MATRIXBLOCK24(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
!          ELSEIF ((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.1)) THEN
!              CALL MATRIXBLOCK15(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
!     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)
          ELSEIF ((NTYPE(1,IC).EQ.5).AND.(NTYPE(1,IR).EQ.2)) THEN
              CALL MATRIXBLOCK25(IC,IR,NCOEC,INCOR,NCTEC,INC2,NMCBP,
     :                NCORE,ELSTO,NTYPE,EMTBLOCK,MAXSPAN,NORBGEN)

          END IF
          IF ((IR .EQ. IC) .OR. (DABS (EMTBLOCK(1,1)) .GT. CUTOFF)) THEN
            NELC       = NELC + 1
            EMT(NELC)  = EMTBLOCK(1,1)
            IROW(NELC) = IR
          ENDIF

   85   CONTINUE            

c zou
        IF(LSE) EMT(NELC) = EMT(NELC) + SLF_EN(IC)
c zou
*
*   This column is done; write it to disk
*
        WRITE (imcdf) NELC, ELSTO, (EMT(IR), IR = 1, NELC),
     :                             (IROW(IR), IR = 1, NELC)
! This EAV (and the above EMT) does not have ELSTO.
        EAV = EAV + EMT(NELC)
*
        IF (MOD (IC, 100) .EQ. 0 .OR.
     &    IC .LT. nprocs*2 .OR. IC .GT. (NCF-nprocs*2)) THEN
          PRINT *, 'Row ', IC, ': ', NELC, ' nonzero elements;'
     &             , '  block = ', jblock
        ENDIF
*
*   Update the counter for the total number of elements
*
        NELMNT = NELMNT + NELC
*
   10 CONTINUE
*
*   Deallocate storage for the arrays in /BUFFER/
*
      CALL ALCBUF (3)

*     ...Locals
      CALL DALLOC (PNTEMT)
      CALL DALLOC (PNIROW)

*  Deallocate arrays related to the symbolic part

      DEALLOCATE(NTYPE)
      DEALLOCATE(IROWSYM)
      DEALLOCATE(EMTSYM)
      DEALLOCATE(EMTBLOCK)
      DEALLOCATE(MAP)

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      CUTOFFtmp = CUTOFF
      NCOEItmp = NCOEI
      NCOECtmp = NCOEC
      NCTEItmp = NCTEI
      NCTECtmp = NCTEC
      NTPItmp = NTPI
      NMCBPtmp = NMCBP
      NCOREtmp = NCORE
      NVPItmp = NVPI
      NKEItmp = NKEI
      NVINTItmp = NVINTI
      NELMNTtmp = NELMNT
      NCFtmp = NCF

      RETURN
      END
