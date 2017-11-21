*************************************************************************
*
      SUBROUTINE INIESTdm (nmax, ncf, NIV, BASIS, hmx)

*  Note that hmx is not global for the mpi version. i.e., each
*  node only has its elements and the size of hmx is also local.
*
*   MPI version by Xinghong He            Last revision: 18 Jun 1998
*   Modified for ifort -i8 by A. Kramida 22 Mar 2016            *
*
*************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
      DIMENSION basis(*)
      POINTER (iqap,ap(1)),(iqeig,eigval(1)),(iqvec,vec(1))
      POINTER (iqwork,work(1)),(iqiwork,iwork(1)),(iqif,IFAIL(1))
      DIMENSION HMX(*)
      POINTER (PNEVAL, EVALDUM)
      COMMON/EIGVAL/EAV,PNEVAL

      !INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
*-----------------------------------------------------------------------
cAK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
cff    Set default integer length to 4
      ISIZE = 4
c
      NS = MIN (nmax, ncf)

      CALL alloc (iqap, NS*(NS+1)/2, 8)
      CALL dinit ((NS*(NS+1))/2, 0.d0, ap, 1)
	
*  Get the upper left sub-matrix
*  ap is global, hmx is local !!!

      ihmx = 0
      DO i = myid + 1, ns, nprocs
         iap = i*(i-1)/2
         DO j = 1, i
            ap(j+iap) = hmx(j+ihmx)
         ENDDO
         ihmx = ihmx + i
      ENDDO

*  Merge ap from all nodes and then send to all nodes

      CALL gdsummpi (ap, (NS*(NS+1))/2)

      CALL alloc (iqeig,   NS,     8)
      CALL alloc (iqvec,   NS*NIV, 8)
      CALL alloc (iqwork,  8*NS,   8)
      CALL alloc (iqiwork, 5*NS,   ISIZE)
      CALL alloc (iqif,     NS,    ISIZE)

!      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
!     :          NS,AP,-1.,-1.,1,NIV,0.d0,
!     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      ABSTOL = 2*SLAMCH('S')
      CALL DSPEVX ('V','I','U',
     :          NS,AP,-1.,-1.,1,NIV,ABSTOL,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      !print*,' in iniestdm '
      !print*, 'myid:',myid, ' INFO=',info
           IERR = -ABS (INFO)


*  Build the Basis. 

      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)

*  scatter the vectors

      DO J = 1, NIV
        CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO

      CALL dcopy (NIV, EIGVAL, 1, BASIS(NIV*ncf+1), 1)

      CALL dalloc (iqap)
      CALL dalloc (iqeig)
      CALL dalloc (iqvec)
      CALL dalloc (iqwork)
      CALL dalloc (iqiwork)
      CALL dalloc (iqif)

      RETURN
      END
