************************************************************************
*
      SUBROUTINE INIESTSD (nmax, ncf, myid, nprocs,
     :    NIV, BASIS, IMCDF, EAV)
      IMPLICIT REAL*8          (A-H, O-Z)
        include 'mpif.h'
*
*    Routine for providing initial estimates from upper left corner
*    of the matrix. It is exact (not estimates) if ncf <= nmax which
*    is set in the calling routine to the IOLPCK parameter .
*
*    Matrix is sparse and on the disk
*
*   MPI version by Xinghong He            Last revision: 14 Dec 1998
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016            *
*
************************************************************************
      DIMENSION basis(*)
      ! Locals
      POINTER (iqap,ap(1)),(iqeig,eigval(1)),(iqvec,vec(1))
      POINTER (iqwork,work(1)),(iqiwork,iwork(1)),(iqif,IFAIL(1))
      POINTER (phmx, hmx(1)), (pirow, irow(1))
!      pointer(iqap_buffer,ap_buffer(1))
!-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
cff   Set the default integer length to 4
      ISIZE = 4
c
      NS = min (nmax, ncf)
      nnn = (NS*(NS+1))/2
      CALL alloc (iqap, (NS*(NS+1))/2, 8)
      call alloc(iqap_buffer,nnn, 8)

      CALL dinit ((NS*(NS+1))/2, 0.d0, ap, 1)

***** separate upper left block of size NS*NS

      CALL alloc (phmx, ncf, 8)
      CALL alloc (pirow, ncf, ISIZE)
      READ (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
c      IF (myid .EQ. 0) PRINT *, 'iniestsd ...........'
      IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs)
     &   STOP 'iniestsd: ncf read wrong'

*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.

      DO j = myid + 1, ns, nprocs
         joff = (j*(j-1))/2
         READ (IMCDF) NELC,ELSTO,(HMX(IR),IR=1,NELC),
     :                          (IROW(IR),IR=1,NELC)
         HMX(NELC) = HMX(NELC) - EAV ! Shift the diagonal
         DO ir = 1, nelc
            ap(irow(ir) + joff) = hmx(ir)
         ENDDO
      ENDDO
      nnn = (NS*(NS+1))/2
! Let each node have a complete copy of ap
!      call MPI_ALLREDUCE(ap,ap_buffer,nnn,
!     :       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
!      ap(1:nnn)=ap_buffer(1:nnn);

      CALL gdsummpi (ap, (NS*(NS+1))/2)

! To be in step with other cases, go through the whole block.
!
! This is not necessary since currently the file pointer is moved
! to the absolute position of the .res files which is always counted
! from the begining of the .res files of each node. Besides, the
! following segment seems not working properly for the last block.
! Xinghong He 98-12-14

      !mylast = j - nprocs
      !DO j = mylast, ncf, nprocs
      !   READ (imcdf)
      !ENDDO

      CALL dalloc (phmx)
      CALL dalloc (pirow)

      CALL alloc (iqeig,NS,8)
      CALL alloc (iqvec,NS*NIV,8)
      CALL alloc (iqwork,8*NS,8)
      CALL alloc (iqiwork,5*NS,ISIZE)
      CALL alloc (iqif,NS,ISIZE)

!      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
!     :          NS,AP,-1.,-1.,1,NIV,0.d0,
!     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      ABSTOL = 2*SLAMCH('S')
      CALL DSPEVX ('V','I','U',
     :          NS,AP,-1.,-1.,1,NIV,ABSTOL,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
       !print*,' In iniestsdi'
       !print*, 'myid:',myid, ' INFO=',info
      IERR = -ABS (INFO)

*******************************************************************

*       ..Build the Basis.

      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)
*       ...scatter the vectors
      DO J = 1, NIV
         CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO
      CALL dcopy (NIV, EIGVAL,1,BASIS(NIV*ncf+1),1)

      CALL dalloc (iqap)
      CALL dalloc (iqeig)
      CALL dalloc (iqvec)
      CALL dalloc (iqwork)
      CALL dalloc (iqiwork)
      CALL dalloc (iqif)

      RETURN
      END
