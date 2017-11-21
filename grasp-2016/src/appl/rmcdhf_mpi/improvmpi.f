************************************************************************
      SUBROUTINE improvmpi (EOL, J, lsort, DAMPMX)
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INCLUDE 'mpif.h'
      LOGICAL EOL, lsort
      INTEGER J

*   The difference from the serial version is that it calls MPI        *
*   version subroutines (setlagmpi, cofpotmpi, matrixmpi, newcompi).   *
*                                                                      *
*   Improve the orbital J.                                             *
*                                                                      *
*   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, matrixmpi,    *
*                         newcompi, ORTHOR, ROTATE, SETCOF, setlagmpi, *
*                         SOLVE, XPOT, YPOT.                           *
*               [LIB92]: ORTHSC, QUAD.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*   Modified by Xinghong He                 Last update: 05 Aug 1988   *
*   Modified for ifort -i8 by A. Kramida (AK) Last update 22 Mar 2016  *
*                                                                      *
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PCDAMP,CDAMPDUM)
      POINTER (PNTRIQ,RIQDUM)
      LOGICAL FAIL,FIRST,ORTHST
      CHARACTER*2 NH

      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      POINTER (PNTNDA,NDA(1)), (PNTNXA,NXADUMMY),(PNTNYA,NYADUMMY)
      POINTER (PNTRDA,DA(NDDIM)), (PNTRXA,XADUMMY),(PNTRYA,YADUMMY)
cbq MPI_ALLREDUCE NDCOF
c      pointer (pda_buffer,da_buffer)
c      pointer (pnda_buffer,nda_buffer)
c      dimension da_buffer(NDDIM)
c      dimension nda_buffer(NDDIM)
c      integer ndcof_buffer, indcof,nda_buffer
c
cAK Dynamic arrays for consolidation of DA and NDA from all nodes
      pointer (pda_buffer,da_buffer)
      pointer (pnda_buffer,nda_buffer)
      pointer (pndcof_buffer,ndcof_buffer)
      integer ndcof_buffer(*),nda_buffer(*)
      DOUBLE PRECISIOn da_buffer(*)
c      PARAMETER (NBUFSIZE = 10000)
c      integer nda_buffer(NBUFSIZE),ndcof_buffer(NBUFSIZE)
c      DOUBLE PRECISIOn da_buffer(NBUFSIZE)
c      integer, dimension(:), allocatable :: nda_buffer,
c     & ndcof_buffer,nda_buf_local,ndcof_buf_local
c      double precision, dimension(:), allocatable :: da_buffer,
c     & da_buffer_loc 

      PARAMETER (P2    = 2.0D-01,
     :           P005  = 5.0D-03,
     :           P0001 = 1.0D-04)
*
*   C Froese Fischer's IPR and ED1 parameter
*
      DATA IPR /0/
      DATA ED1 /0.D0/
      DATA FIRST /.FALSE./

      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/

      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
      ISIZE = sizeof(NCF)
      I_MPI = MPI_INTEGER
      if (ISIZE.EQ.8) I_MPI = MPI_INTEGER8
c
      GAMAJ = GAMA(J)
*
*   C Froese Fischer's parameters IPR, ED1, ED2 are set and
*   used in this routine and in DAMPCK
*
    1 ED2 = E(J)
*
*   Set up the exchange potential and arrays XU, XV as appropriate
*
*   Set coefficients for YPOT, XPOT, DACON
*   Compute direct potential, exchange potential
*   Add in Lagrange-multiplier contribution
*   Add in derivative-terms contribution
*
      npts = N
!      print *, "improvmpi0:NDCOF",NDCOF,myid
      CALL cofpotmpi (EOL, J, npts)
!      print *, "improvmpi1:DA(1),NDCOF",DA(1),NDCOF,myid
*
*   Calculate deferred corrections
*
      CALL DEFCOR (J)
*
*   Solve the Dirac equation

cAK     NDCOF, NDA, and DA cannot be reduced by finding the max of NDCOF
c          and summing up DA from all nodes. The correct replacement is below
cbq MPI_ALLREDUCE NDCOF
c      call MPI_ALLREDUCE(ndcof,ndcof_buffer,1,
c     :          I_MPI,MPI_MAX,MPI_COMM_WORLD,ierror)
c      if (ndcof_buffer .gt. ndcof) then
!        print *, 'improvmpi: ndcof, ndcof_buffer, myid=',
!     &                     ndcof, ndcof_buffer, myid
c        write (*,'(A35,4I10)') 'improv11: myid,NDCOF,ndcbfr=',
c     &    myid,NDCOF,ndcof_buffer
c        do indcof = ndcof+1, ndcof_buffer
c           da(indcof) = 0.0
c           nda(indcof) = 0
c        enddo
c        ndcof = ndcof_buffer
c      endif
*
      INV = 0
*
cbq MPI_ALLREDUCE DA
c      if (ndcof.gt.0) then
c         call alloc(pda_buffer, ndcof, 8)
c         call alloc(pnda_buffer, ndcof, ISIZE)
c         call MPI_ALLREDUCE(da,da_buffer,ndcof,
c     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
c         da(1:ndcof) = da_buffer(1:ndcof)
c      nda(1:ndcof) = nda_buffer(1:ndcof)
c         call dalloc(pda_buffer)
c         call dalloc(pnda_buffer)
c      ENDIF
cAK      Consolidate NDA and DA from all nodes
      if (nprocs.gt.1) then
         call alloc(pndcof_buffer, nprocs, ISIZE)
         call MPI_GATHER(ndcof, 1, I_MPI, ndcof_buffer, 
     &        1, I_MPI, 0, MPI_COMM_WORLD, ierr)
         CALL MPI_Bcast (ndcof_buffer, nprocs, I_MPI, 0, 
     &                         MPI_COMM_WORLD, ierr)
         ndcof_max = 0
         do i = 1, nprocs
            if (ndcof_buffer(i).gt.ndcof_max)
     &          ndcof_max = ndcof_buffer(i)
         ENDDO
         if (ndcof_max.gt.0) then
            ntot = ndcof_max*nprocs
            call alloc(pnda_buffer, ntot, ISIZE)
            call MPI_GATHER(nda, ndcof_max, I_MPI, nda_buffer, 
     &           ndcof_max, I_MPI, 0, MPI_COMM_WORLD, ierr)
            CALL MPI_Bcast (nda_buffer, ntot, I_MPI, 0, 
     &                         MPI_COMM_WORLD, ierr)
            call alloc(pda_buffer, ntot, 8)
            call MPI_GATHER(da, ndcof_max, MPI_DOUBLE_PRECISION, 
     &        da_buffer, ndcof_max, MPI_DOUBLE_PRECISION, 0, 
     &        MPI_COMM_WORLD, ierr)
            CALL MPI_Bcast (da_buffer, ntot, MPI_DOUBLE_PRECISION, 
     &                         0, MPI_COMM_WORLD, ierr)
            i_last = 0
            do iproc = 1, nprocs
               ndcip = ndcof_buffer(iproc)
               do jproc = 1, ndcip
                  ifound = 0
                  ind_buf = (iproc - 1)*ndcof_max + jproc
                  do k = 1, i_last
                     if (nda_buffer(ind_buf) .eq. NDA(k)) THEN
                        DA(k) = DA(k) + da_buffer(ind_buf)
                        ifound = k
                        exit
                     endif
                  enddo
                  if (ifound.eq.0) then
                     i_last = i_last + 1
                     do while (i_last .GT. NDDIM)
                        IF (NDDIM .GT. 0) THEN
                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 2)
                        ELSE
                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 1)
                        ENDIF
                     ENDdo
                     NDA(i_last) = nda_buffer(ind_buf)
                     DA(i_last) = da_buffer(ind_buf)
                  endif
               enddo
            enddo
            NDCOF = ndcof_max
            call dalloc(pndcof_buffer)
            call dalloc(pnda_buffer)
            call dalloc(pda_buffer)
         endif
      endif

!     print *,"improvmpi before: p0=",p0, "myid=",myid
!     call flush(6)
      CALL SOLVE (J, FAIL, INV, JP, NNP)
!      print *,"improvmpi after: p0=",p0, "myid=",myid
!     call flush(6)
*
*   Upon failure issue message; take corrective action if possible
*
      IF (FAIL) THEN
         IF (myid .EQ. 0) WRITE (*,300) NP(J),NH(J),METHOD(J)
         IF (METHOD(J) .NE. 2) THEN
            METHOD(J) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            CALL ORTHSC
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN
!               CALL matrixmpi
!               CALL newcompi(WTAEV)
!            ENDIF
            CALL setlagmpi (EOL)
            GOTO 1
         ELSE
            IF (myid .EQ. 0) WRITE (*,301)
            !CALL TIMER (0)
            STOP
         ENDIF
      ENDIF
*
*   Compute norm of radial function
*
      TA(1) = 0.D0
      DO I = 2, MTP0
         TA(I) = (P(I)**2 + Q(I)**2) * RP(I)
      ENDDO
      MTP = MTP0

      CALL QUAD (DNORM)
!      pRINT *, dnorm,'dnorm',myid,'?'

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      CALL CONSIS (J)
*
*   Normalize
*
      DNFAC = 1.D0 / DSQRT (DNORM)
      P0 = P0 * DNFAC
!      print*,'dnfac=', dnfac,'myid=', myid
      DO I = 1, MTP0
         P(I) = P(I) * DNFAC
         Q(I) = Q(I) * DNFAC
      ENDDO
*
*   Check if different method should be used or if improvement
*   count should be reduced
*
      DEL1 = DABS (1.D0 - ED2 / E(J))
      IF (METHOD(J) .EQ. 1) THEN
         DEL2 = DMAX1 (DABS (1.D0 - DSQRT (DNORM)),
     :               DABS (DNFAC - 1.D0))
         IF ((DEL1 .LT. P005) .AND. (DEL2 .GT. P2)) THEN
            METHOD(J) = 2
!           print*,method(j),'=method(j)','?',myid
            GOTO 1
         ENDIF
      ELSE
         IF ((DEL1 .LT. P0001) .AND. (NSIC .GT. 1)) NSIC = NSIC-1
      ENDIF
*
*   Damp the orbital --- if not converged
*
       IF (scnsty(J) .GT. ACCY) THEN
          CALL DAMPCK (IPR, J, ED1, ED2)
          odampj = DABS (odamp(j))
       ELSE
          odampj = 0.D0    ! take the whole new orbital
       ENDIF
      CALL DAMPOR (J, INV, odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      IF (ORTHST) THEN
         !CALL orthor (J, inv)
         nwww = nw
         CALL orthy (nwww, J, lsort)
      ENDIF
*
*   Print details of iteration
*
      IF (myid .EQ. 0)
     &   WRITE (*,302) NP(J),NH(J),E(J),METHOD(J),PZ(J),SCNSTY(J),
     &              DNORM-1,ODAMPJ,JP,MF(J),INV,NNP
      DAMPMX = DMAX1(DAMPMX,DABS(ODAMPJ))

  300 FORMAT (/' Failure; equation for orbital ',1I2,1A2,
     :         ' could not be solved using method ',1I1)
  301 FORMAT (//' ****** Error in SUBROUTINE IMPROV ******'
     :          /' Convergence not obtained'/)
cb
cb 80-col
cb
c 302 FORMAT (1X,1I2,1A2,1P,1D18.9,1x,1I2,D11.3,1D10.2,1D10.2,
  302 FORMAT (1X,1I2,1A2,1P,1D16.7,1x,1I2,D11.3,1D10.2,1D10.2,
c    :        0P,F7.4,1x,1I3,1x,1I3,1x,1I2,2x,1I2)
c    :        0P,F6.3,1x,1I4,1x,1I4,1x,1I2,1x,1I2)
     :        0P,F6.3,1x,1I5,1x,1I5,1x,1I2,1x,1I2)

      RETURN
      END

