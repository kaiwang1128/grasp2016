************************************************************************
*                                                                      *
      SUBROUTINE lodres
*                                                                      *
*   Loads the data from the  .res  file. A number of checks are made   *
*   to ensure correctness and consistency.                             *
*   This subroutine is called by every node and the only differences
*   from the non-mpi version are:
*     1) STOP and stopmpi (not necessary)
*     2) LSE needs coomunication - can be moved out of this routine. 
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*   MPI version by Xinghong He            Last revision:  1 Jun 1998   *
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016            *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL GETYN,LFORDR,LTRANS,LVP,LSE,LNMS,LSMS,YES
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /FOPARM/ICCUT
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WFAC/WFACT
     :      /WHERE/IMCDF

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
      L_MPI = MPI_LOGICAL
*
*   Read the basic parameters of the electron cloud; check these
*   against those deduced from the  .csl file
*
      READ (imcdf) NELECR, NCFRES, NWRES, nblockres

      IF (NELECR .NE. NELEC .OR. NCFRES .NE. NCF .OR.
     &    NWRES .NE. NW .OR. nblockres .NE. nblock) THEN
         CALL stopmpi ('lodres: NELEC/NCF/NW does not match', myid)
      ENDIF
*
*   Read the nuclear parameters
*
      READ (imcdf) Z,EMN
      READ (imcdf) NPARM,(PARM(I),I = 1,NPARM)
      READ (imcdf) N,(ZZ(I),I = 1,N),NNUC

      IF (N .GT. NNNP) THEN
         CALL stopmpi ('lodres: N greater than NNNP', myid)
      ENDIF
*
*   Read the physical effects specifications
*   iccutblk() is now an array of length nblock.
*
      READ (imcdf) C, LFORDR, (ICCUTblk(i), i = 1, nblock),
     &             LTRANS, WFACT, LVP, LNMS, LSMS
*
*   Read the remaining parameters controlling the radial grid and the
*   grid arrays
*
      NP10 = N+10
      READ (imcdf) RNT,H,HP,
     :          (R(I),I = 1,NP10),
     :          (RP(I),I = 1,NP10),
     :          (RPOR(I),I = 1,NP10)
*
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Allocate storage for the radial wavefunction arrays
*
      CALL ALLOC (PNTRPF,NNNP*NW,8)
      CALL ALLOC (PNTRQF,NNNP*NW,8)
*
*   Read the orbital wavefunctions and the associated arrays
*
      DO 1 J = 1, NW
         READ (imcdf) E(J), GAMA(J), PZ(J), MF(J)
         READ (imcdf) (PF(I,J), I = 1, MF(J)), (QF(I,J), I = 1, MF(J))
    1 CONTINUE
*
*   Determine if the self-energy contribution is to be estimated
*
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) 'Estimate contributions from the self-energy?'
         LSE = GETYN ()
      ENDIF

      CALL MPI_Bcast (LSE,1,L_MPI,0,MPI_COMM_WORLD,ierr)

      RETURN
      END
