      SUBROUTINE MANEIG (iatjpo, iaspar, NELMNT_a )
*   eval(), evec(), iatjpo, iaspar
*                                                                      *
*   This module  manages the  operation of the  eigensolvers and the   *
*   storage of the eigenpairs.  There are two principal branches:      *
*                                                                      *
*      (1) Matrix of order 1: the trivial case                         *
*      (2) Matrix of order greater than 1: there are two branches      *
*             (i) Matrices of order less than or equal to IOLPCK:      *
*                 eigenpairs are found using LAPACK SUBROUTINEs        *
*            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
*                 are found using DVDSON; this involves up to three    *
*                 steps:                                               *
*                    (a) The matrix is analysed to determine its       *
*                        block structure (only irreducibe matrices     *
*                        are correctly treated by DVDSON)              *
*                    (b) Eigenpairs are extracted for each block       *
*                    (c) The appropriate eigenpairs are selected and   *
*                        stored                                        *
*                 Different methods of storage and different           *
*                 versions of the matrix-vector multiply are used      *
*                 depending upon the order and density of the matrix   *
*   
*   iatjpo - Output. (2j+1) of the dominant coefficients
*   iaspar - Output, Parity (1 or -1) of the dominant coefficients
*   nelmnt - Number of non-zero matrix elements on the current process *
*            stored in the common block                                *
*   nelmnt_a - The maximum number of non-zero matrix elements for all  *
*             nodes                                                    *
*                                                                      *
*     value  meaning  of  method used IV                               *
*       1    LAPACK                  =>  ncf < IOLPCK                  * 
*       2    Davidson, Dense-Memory  =>  dense requires less memory    *
*       3    Davidson, Sparse-Memory =>  sparse requires less memory   *
*       4    Davidson, Sparse-Disk   =>  memory requirement too large  *
*                                        for at least one  process     *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, posfile,         *
*                        RALLOC.                                       *
*               [RCI92]: DNICMV, SPICMVmpi, SPODMV.                    *
*               [DVDSON]: DVDSON.                                      *
*               [AUXBLAS]: DINIT/SINIT.                                *
*               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
*               [LAPACK]: DSPEVX/SSPEVX.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*   Modified Misha Saparov                                  Feb 1997   *
*            Charlotte F. Fischer                           May 1997   *
*                 Except for the disk version, all matrices have       *
*                 diagonals  shifted by EAV                            *
*  All arrays allocated here are de-allocated except pneval and pnevec *
*  which will be de-allocated in matrix.                               *
*   MPI Version  By Xinghong He           Last revision: 29 Jul 1998   *
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016                   *
*   Modified for automatic method use:    CFF, April 2016              *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      EXTERNAL DNICMV,SPICMVmpi,SPODMV

      POINTER (PNEVAL,EVAL(1))
      COMMON/EIGVAL/EAV,PNEVAL

      POINTER (PNEVEC,EVEC(1))
      COMMON/EIGVEC/PNEVEC

!  nposition+1 is the current position of the .res file
!  It is set in matrix and used in maneig, spodmv

      COMMON/fposition/nposition

      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/ncf,NW,PNTRIQ

      POINTER (PNIVEC,IVEC(1))
      COMMON/PRNT/NVEC,PNIVEC,NVECMX
cAK    INTEGER*4 variable required for compatibility with ifort -i8
      INTEGER*4 IVEC

      COMMON/WCHBLK/JBLOCK
      COMMON/WHERE/IMCDF
      COMMON/iounit/istdi,istdo,istde

      LOGICAL HIEND,LDISC,SPARSE
      CHARACTER*8 CNUM

      POINTER (PONTRW,W(1))
      POINTER (PONTRZ,Z(1))
      POINTER (PNWORK,WORK(1))
      POINTER (PIWORK,IWORK(1))
      POINTER (PIFAIL,IFAIL(1))
      POINTER (PNDIAG,DIAG(1))
      POINTER (PJWORK,JWORK(1))
!     length of work arrays (real, integer) needed for Davidson
      INTEGER*8   LWORK, LIWORK

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER msg*80

!  ABSTOL is the absolute error tolerance for the eigenvalues
!  used by lapack routine
      PARAMETER (ABSTOL = 1.0D-10)

!  IOLPCK  is the largest matrix to be solved using LAPACK routines
      INTEGER  IOLPCK 
      PARAMETER (IOLPCK = 2000)

!  NINCOR  is the maximum amount of memory that can be allocated to
!     to representation of the distributed Hamiltonian, for the 
!     current process
      INTEGER*8  NINCOR
      PARAMETER (NINCOR = 134217728) ! = 1 GB or  memory

!cff   ... variable for deciding between sparse or dense
!     Memory need:  
!     sparse : nstore_s = nelmnt + (isize*(nelmnt + ncf+1))/8
!     dense  : nstore_d = (ncf*ncf+1))/(2*nprocs)    (on average)
!     Use Davidson (disk)  if nstore_s > NINCOR
!     Use Davidson (memory, sparse) if nstore_s < nstore_d
!     Use Davidson (memory, dense) if nstore_d < nstore_s
      INTEGER*8 NSTORE_s, NSTORE_d
      INTEGER*8 NELMNT_a, NELMNT

!cff  ... add NDENSE, the memory needed if LAPACK is used.
!     ... This is the memory needed for a trianglular matrix of
!     ... size < IOLPCK
      INTEGER  NDENSE
      integer*8 ncft

!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
CFF   Set default integer length to 4
      ISIZE = 4
c
       IF (myid .EQ. 0) PRINT *, 'Calling maneig...'
       ncft = ncf
!  (nrows+1) is the number of records of the present block's .res file

      nrows = (ncf - myid - 1 + nprocs) / nprocs
      IF (ncf .LT. nprocs) nrows = ncf / (myid+1)
      !CALL posfile (1, imcdf, nrows+1)
      CALL posfile (0, imcdf, nposition)

      IF (ncf .EQ. 1) THEN
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         IF (myid .EQ. 0) WRITE (24,*) 'Trivial eigenvalue problem.'

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         CALL ALLOC (PNEVAL, 1, 8)
         CALL ALLOC (PNEVEC, 1, 8)
         EVAL(1) = 0.D0
         EVEC(1) = 1.D0

         ! Still read through the .res file
CGG
CGG Gediminas NIST 2005.11.03
CGG         READ (imcdf)
         DO i = 1, nrows + 1
            READ (imcdf)
         ENDDO

      ELSE	!if-2
!-----------------------------------------------------------------------
!
! (2) - Non trivial
!
!-------------------------------------------------------
         if(myid==0) write(*,*) "ncf=",ncf, " iolpck=", iolpck
         IF (ncf .LE. IOLPCK) THEN
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory, 
!
!-------------------------------------------------------
            IF (myid .EQ. 0) THEN
               PRINT *, 'LAPACK routine DSPEVX'
     :                 //' selected for eigenvalue problem.'
               WRITE (24,*) 'LAPACK routine DSPEVX'
     :                 //' selected for eigenvalue problem.'
            ENDIF

!   Allocate storage for the dense representation of the matrix
!   and initialize emt
            NDENSE = (ncf*(ncf+1))/2
            CALL ALLOC (PNTEMT,NDENSE,8)
            CALL DINIT (NDENSE,0.0D 00,EMT,1)
!
!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            CALL ALLOC (PNWORK,ncf,8)
            CALL ALLOC (PNIROW,ncf,ISIZE)
            READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
            IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs) 
     &          STOP 'maneig:1'

            DO 2 I = myid + 1, ncf, nprocs
               iofset = (i*(i-1))/2
               READ (IMCDF) NELC,ELSTO,(WORK(IR),IR = 1,NELC),
     :                                 (IROW(IR),IR = 1,NELC)
               ! In the row-mode of the lower triangle, 
               ! diagonal is the last one
               DO 1 IR = 1, NELC - 1
                  EMT(iofset+IROW(IR)) = WORK(IR)
    1          CONTINUE
               EMT(iofset+IROW(NELC)) = WORK(NELC)-EAV

    2       CONTINUE

            CALL DALLOC (PNWORK)
            CALL DALLOC (PNIROW)

! Let each node have a complete copy of EMT

            CALL gdsummpi (EMT, NDENSE)
!
!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    nvecmx - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1
!
            NVECMN = ncf
            DO 3 I = 1,NVEC
               NVECMN = MIN (NVECMN,IVEC(I))
    3       CONTINUE
            NVEX = NVECMX-NVECMN+1
            CALL ALLOC (PONTRW,    NVEX,8)
            CALL ALLOC (PONTRZ,ncf*NVEX,8)
            CALL ALLOC (PNWORK,ncf*8   ,8)
            CALL ALLOC (PIWORK,ncf*5   ,ISIZE)
            CALL ALLOC (PIFAIL,    NVEX,ISIZE)
            CALL DSPEVX ('V','I','U',ncf,EMT,DUMMY,DUMMY,
     :                           NVECMN,NVECMX,ABSTOL,M,W,Z,ncf,
     :                           WORK,IWORK,IFAIL,INFO)
            IF (INFO .NE. 0) THEN
               CALL stopmpi ('maneig: Failure in DSPEVX [LAPACK]', myid)
            ENDIF
            CALL DALLOC (PNWORK)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PIFAIL)
            CALL DALLOC (PNTEMT)

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            CALL ALLOC (PNEVAL,    NVEC,8)
            CALL ALLOC (PNEVEC,ncf*NVEC,8)
            DO 4 I = 1,NVEC
               LOC = IVEC(I)
               EVAL(I) = W(LOC-NVECMN+1)
               iofset = ncf*(I-1)
               LOC = ncf*(LOC-NVECMN)
               CALL DCOPY (ncf,Z(LOC+1),1,EVEC(iofset+1),1)
    4       CONTINUE
            CALL DALLOC (PONTRW)
            CALL DALLOC (PONTRZ)

         ELSE
!-----------------------------------------------------------------------
!
! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
            IF (myid .EQ. 0) WRITE (24,*) 'DVDSON routine selected'
     :                 //' for eigenvalue problem;'
!--------------------------------------------------------------
            IF (myid .EQ. 0) print *,'zou: nelmnt',NELMNT    
            if (myid .eq. 0) 
     &         print *, 'myid, nprocs,nincor', myid, nprocs, nincor
        
c           If on any process more memory is needed than NINCOR, use disc
!           If not, should the matrix be sparse or dense (distributed)
!           sparse requires = nelmnt + (isize*(nelmnt+ncf+1))/8
!           dense  storage = (ncf*(ncf+1)/2/nproc - on average
!           (because of integer arithmetic nproc-1 is added to numerator
!           
            NSTORE_s = NELMNT_a + ((NELMNT_a + ncf+1)*ISIZE )/8
            NSTORE_d = ((ncf+nprocs-1)/nprocs)
            NSTORE_d = NSTORE_d*((ncf+1)/2)
            IF (myid .EQ. 0) write(*,*) "nstore_s", nstore_s,
     :                                     "nstore_d", nstore_d
            LDISC = .FALSE.
            !IF ( NSTORE_s > NINCOR) THEN
!              .. we need to use disk
            !   LDISC =.TRUE.
               SPARSE=.TRUE.
            !ELSE
            
            !   IF (NSTORE_s .LT. NSTORE_d ) THEN
            !      SPARSE = .TRUE.
            !   ELSE
            !      SPARSE = .FALSE.
            !   ENDIF
            !ENDIF

            CALL ALLOC (PNDIAG,ncf,8)
            CALL dinit (ncf, 0.d0, diag, 1)
            IF (LDISC) THEN
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal   ( sparse)
!
!-------------------------------------------------------
               IF (myid .EQ. 0)  THEN
                   WRITE (*,*) ' matrix stored on disc;'
                   WRITE (24,*) ' matrix stored on disc;'
               END IF
!
!   Disk storage; necessarily sparse; one column of the matrix in
!   memory
               IMV = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block 
              
               print *, myid, 'Reading diagonal, unit imcdf =', imcdf
               READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
               IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs) 
     &            STOP 'maneig:2'

               DO I = myid + 1, ncf, nprocs
                  READ (IMCDF) NELC,ELSTO,(dummy,IR=2,NELC), diatmp
     &                                 , (idummy, ir=1,nelc)
                  DIAG(I) = diatmp - EAV
               ENDDO
               print *, ' finished'

            ELSE
!-----------------------------------------------------------------------
!
! (2.2.2) - DVDSON --- Memory, load all
!
!-------------------------------------------------------
!
!   Core storage; load matrix into memory

               LDISC = .FALSE.
               IF (SPARSE) THEN
!-----------------------------------------------------------------------
!
! (2.2.2.1) - DVDSON --- Memory, load all, sparse
!
!-------------------------------------------------------
                  IF (myid .EQ. 0) THEN
                     WRITE (24,*) ' sparse matrix stored in memory'
                     WRITE ( *,*) ' sparse matrix stored in memory'
                  END IF

                  IMV = 2

                  CALL ALLOC2 (PNTEMT,NELMNT,8)
                  !IF (myid .EQ. 0) THEN
                  !   WRITE (24,*) ' Exit from alloc pntemt'
                  !   WRITE ( *,*) ' Exit from alloc pntemt'
                  !END IF
                  CALL ALLOC2 (PNIROW,NELMNT,ISIZE)
                  !IF (myid .EQ. 0) THEN
                  !   WRITE (24,*) ' Exit from alloc pnirow'
                  !   WRITE ( *,*) ' Exit from alloc pnirow'
                  !END IF
                  CALL ALLOC (PIENDC,ncf+1,ISIZE)
                  !IF (myid .EQ. 0) THEN
                  !   WRITE (24,*) ' Exit from alloc piendc'
                  !   WRITE ( *,*) ' Exit from alloc piendc'
                  !END IF
                  iofset = 0
                  IENDC(0) = 0
                  READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
                  IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &                  .OR. nprocsdum .NE. nprocs) 
     &               STOP 'maneig:3'
                  DO I = myid + 1, ncf, nprocs
                     READ (IMCDF) NELC,ELSTO,
     :                            (EMT(IR+iofset),IR = 1,NELC),
     :                           (IROW(IR+iofset),IR = 1,NELC)
                     EMT(NELC+iofset) = EMT(NELC+iofset) - EAV
                     DIAG(I) = EMT(NELC+iofset)
                     iofset = iofset + NELC
                     if (iofset.lt.0) then
                        write(24,*)'Error, iofset.lt.0 1!'
                        write( *,*)'Error, iofset.lt.0 1!'
                        stop 
                     endif
                     IENDC(I) = iofset
                  ENDDO
               ELSE
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  IF (myid .EQ. 0) THEN
                     WRITE (24,*) ' full matrix stored in in memory'
                     WRITE ( *,*) ' full matrix stored in in memory'
                  END IF

                  IMV = 3

! Find NDENSE_L, the number of elements on the node (dense form)

                  NDENSE_L = 0
                  DO i = myid + 1, ncf, nprocs
                     NDENSE_L = NDENSE_L + i
                  ENDDO

                  CALL ALLOC2 (PNTEMT,NDENSE_L,8)
                  CALL DINIT (NDENSE_L,0.0D 00,EMT,1)
                  CALL ALLOC (PNWORK,ncf,8)
                  CALL ALLOC (PNIROW,ncf,ISIZE)

                  READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
                  IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &                  .OR. nprocsdum .NE. nprocs) 
     &               STOP 'maneig:4'

                  iofset = 0
                  DO 8 I = myid + 1, ncf, nprocs
                     READ (IMCDF) NELC,ELSTO,
     :                            (WORK(IR),IR = 1,NELC),
     :                            (IROW(IR),IR = 1,NELC)
                     work(nelc) = work(nelc) - eav
                     DIAG(I) = WORK(NELC)
                     DO IR = 1, NELC
                        EMT(iofset+IROW(IR)) = WORK(IR)
                     ENDDO
                     iofset = iofset + i
    8             CONTINUE
                  CALL DALLOC (PNWORK)
                  CALL DALLOC (PNIROW)

               ENDIF
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            ENDIF
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------

! Make diagonals global, no matter it is disk or memory mode
 
            CALL gdsummpi (DIAG, ncf)
            if (myid.eq.0) print *, 'Got diag, imv=',imv

!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained

CPJ            LIM = MIN (ncf,NVECMX+20)
            LIM = MIN (ncf,2*NVECMX+60)
c           LWORK = 2*ncf*LIM+LIM*LIM+(NVECMX+10)*LIM+NVECMX
            LWORK = 2*ncft*LIM+LIM*LIM*2+11*LIM+NVECMX
            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' LWORK=', LWORK
               WRITE (24,*) '   NCF=', NCF   
               WRITE (24,*) '   LIM=', LIM   
               WRITE (24,*) 'NVECMX=', NVECMX
               WRITE ( *,*) ' LWORK=', LWORK
               WRITE ( *,*) '   NCF=', NCF   
               WRITE ( *,*) '   LIM=', LIM   
               WRITE ( *,*) 'NVECMX=', NVECMX
            END IF
            if (myid.eq.0.and.LWORK.lt.0) stop "ERROR, LWORK.lt.0"
            CALL ALLOC2 (PNWORK,LWORK,8)
            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' Exit from alloc pnwork'
               WRITE ( *,*) ' Exit from alloc pnwork'
            END IF
            work(1:lwork) = 0.0d0 
            LIWORK = 6*LIM+NVECMX
            CALL ALLOC2 (PIWORK,LIWORK,ISIZE)
            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' Exit from alloc piwork'
               WRITE ( *,*) ' Exit from alloc piwork'
            END IF
!*changed by Misha 02/12/97
            CRITE = 1.0D-17
            CRITC = 1.0D-09
            CRITR = 1.0D-09
!            ORTHO = max(1D-8,CRITR)
            ORTHO = dmax1(1D-9,CRITR)
! end of changes

            MAXITR = MAX (NVECMX*100,200)
            CALL ALLOC (PJWORK,LIM,ISIZE)

            CALL ALLOC (PNEVAL,    NVECMX,8)
            CALL ALLOC2 (PNEVEC,ncft*NVECMX,8)

            DMUNGO = 10.D99
            CALL DINIT (NVECMX,DMUNGO,EVAL,1)

            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' Exit from alloc PJWORK, PNEVAL, PNEVEC'
               WRITE ( *,*) ' Exit from alloc PJWORK, PNEVAL, PNEVEC'
            END IF
            
!   Compute the eigenpairs in each block

            NVEX = nvecmx
            IF (LDISC) THEN
               MBLOCK = NVEX
            ELSE
               MBLOCK = 1
            ENDIF
            NEND = ncf*NVEX

            ILOW = 1
            IHIGH = NVEX
            NIV = NVEX
!***********************************************************************
!
!   Call Davidson eigensolver
!
            IF (IMV .EQ. 1) THEN
!******************** sparse and matrix on disk **********************
              IF (myid .EQ. 0) print *, ' Sparse - Disk, iniestsd'
              CALL posfile (0, imcdf, nposition)  ! was within iniestsd before
!             print *, 'Calling iniestsd '
              CALL INIESTSD (IOLPCK, ncf, myid, nprocs, NIV, work, 
     :            IMCDF, EAV)
              print *, 'Returned from iniestsd '
              if (ncf.gt. IOLPCK) then
                print *, 'Calling GDVD'
                CALL GDVD (SPODMV,ncf,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
               print *, 'Returned from GDVD'
              end if
            ELSEIF (IMV .EQ. 2) THEN
                !print*, ' imv = 2'
!******************** sparse and matrix in memory ********************
               IF (myid .EQ. 0) print *, ' Sparse - Memory, iniestmpi'
               CALL iniestmpi (IOLPCK, ncf,NIV,WORK,EMT,IENDC,IROW)
               IF (myid .EQ. 0) THEN
                  WRITE (24,*) ' Exit from iniestmpi'
                  WRITE ( *,*) ' Exit from iniestmpi'
               END IF

               if(ncf.gt. IOLPCK) then

                 IF (myid .EQ. 0) THEN
                     WRITE (24,*) ' Calling GDVD'
                     WRITE ( *,*) ' Calling GDVD'
                 END IF
                 CALL GDVD (SPICMVmpi,ncf,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
                 IF (myid .EQ. 0) THEN
                    WRITE (24,*) ' Exit from GDVD'
                    WRITE ( *,*) ' Exit from GDVD'
                 END IF
               end if

               CALL dalloc (pntemt)
               CALL dalloc (pnirow)
               CALL dalloc (piendc)

            ELSEIF (IMV .EQ. 3) THEN
              !print*, ' imv = 3'
!*************************** dense and in memory **********************
c                  IF (myid .EQ. 0) print *, ' Dense - Memory, iniestdm'
!                  CALL INIESTDM (2000,ncf,NIV,WORK,EMT)
              CALL INIESTDM (IOLPCK,ncf,NIV,WORK,EMT)
              if (ncf.gt. IOLPCK) then
                CALL GDVD (DNICMV,ncf,LIM,DIAG,ILOW,IHIGH,
     :              JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :              WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :              NMV,IERR)
              end if
              CALL dalloc (pntemt)
            ENDIF
!***********************************************************************
            CALL DALLOC (PNDIAG)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PJWORK)

            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' ', nloops, ' iterations;'
               WRITE (24,*) ' ', nmv,' matrix-vector multiplies.'
            ENDIF

            IF (IERR .NE. 0) THEN
               WRITE (istde,*) 'MANEIG: Returned from DVDSON with'
               WRITE (istde,*) ' IERR = ',IERR,'.'
               CALL stopmpi ('maneig: DVDSON wrong', myid)
            ENDIF

!   Put the eigenpairs in order, overwriting as necessary
            if (NEND.lt.0) stop "Error, NEND.lt.0"
            CALL DCOPY (NVEX    , WORK(NEND+1), 1, EVAL, 1)
            CALL DCOPYi8 (ncft*NVEX, WORK(1)     , 1, EVEC, 1) !dcopyi8 in dvdson.f
            CALL DALLOC (PNWORK)

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

            IF (NVEC .LT. NVECMX) THEN
               CALL ALLOC (PIWORK,NVECMX,ISIZE)
               DO I = 1,NVECMX
                  IWORK(I) = I
               ENDDO
               DO I = 1,NVEC
                  iofset = IVEC(I)
                  LOC = IWORK(I)
                  IF (iofset .NE. LOC) THEN
                     CALL DSWAP (1,EVAL(iofset),1,
     :                             EVAL(I     ),1)
                     IWORK(I) = IWORK(iofset)
                     IWORK(iofset) = LOC
                     iofset = ncf*(iofset-1)
                     if (iofset.lt.0) stop "Error, iofset.lt.0 2!"
                     LOC = ncf*(I-1)
                     CALL DSWAP (ncf,EVEC(iofset+1),1,
     :                               EVEC(LOC   +1),1)
                  ENDIF
               ENDDO
               if (ncf*NVECMX.lt.0) stop "ncf*NVECMX.lt.0"
               CALL DALLOC (PIWORK)
               CALL RALLOC (PNEVAL,    NVECMX,    NVEC,8)
               CALL RALLOC (PNEVEC,ncf*NVECMX,ncf*NVEC,8)
            ENDIF

         ENDIF
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      ENDIF
! (3e)               *** E N D   O F    A L L

!--------------------------------------------------------------------
! Only the following quantities are needed after this routine is
! finished:
!   eval(), evec(), iatjpo, iaspar
!--------------------------------------------------------------------
!
!   Clean up eigenvectors; determine their J/P values
!
      DO 23 J = 1, NVEC

!         Find the dominant component of each eigenvector

         iofset = (J-1)*ncf
         if (iofset.lt.0) stop "Error, iofset.lt.0 3!"

         AMAX = 0.d0
         DO  I = 1, ncf
            WA = DABS (EVEC(I+iofset))
            IF (WA .GT. AMAX) THEN
               AMAX = WA
               IA = I
            ENDIF
         ENDDO

!          Find the angular momentum and parity of the dominant component

         iatjpo = ITJPO (IA)
         iaspar = ISPAR (IA)

!          Change sign of eigenvactor if dominant component is negative

         IF (EVEC(IA+iofset) .LT. 0.d0) THEN
           DO I = 1, ncf
              EVEC(I+iofset) = -EVEC(I+iofset)
           ENDDO
         ENDIF
   23 CONTINUE

      RETURN
      END
