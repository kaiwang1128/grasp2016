************************************************************************
*                                                                      *
      SUBROUTINE SPODMV (N,M,B,C)
*                                                                      *
*   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
*   lower triangle of the  (NxN)  matrix  A  is read from the disk     *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
*               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI.                  *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 13 Oct 1992   *
*   MPI Version by Xinghong He            Last revision: 17 Aug 1998   *
*   Modified for ifort -i8 by A. Kramida,                 8 Apr 2016   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

      POINTER (PNEVAL,EVALDUMMY)
      COMMON/EIGVAL/EAV,PNEVAL
     :      /WHERE/IMCDF

      !...nposition+1 is the current position of the .res file
      !...It is set in matrix and used in maneig, spodmv
      COMMON/fposition/nposition

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

      POINTER (PNTEMT,EMT(1))     ! AK: Restored, as it is necessary for ifort
      POINTER (PNIROW,IROW(1))   ! AK: Restored, as it is necessary for ifort
cAK      DIMENSION EMT(N), IROW(N)
      DIMENSION B(N,M), C(N,M)
!-----------------------------------------------------------------------
      ISIZE = sizeof(N)  ! AK: handling the ifort -i8 compiler option

      ncf = n

*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON

      CALL DINIT (N*M, 0.D0, C, 1)

      !...moved from maneig before "CALL GDVD (SPODMV..."
      CALL posfile (0, imcdf, nposition)

      READ (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
      IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &      .OR. nprocsdum .NE. nprocs) 
     &   STOP 'spodmv: ncf read wrong'

      CALL alloc (pntemt, ncf, 8)    ! AK: Restored, as it is necessary for ifort
      CALL alloc (pnirow, ncf, ISIZE)  ! AK: Restored, as it is necessary for ifort

      DO ICOL = myid + 1, N, nprocs
         READ (IMCDF) NELC,ELSTO,(EMT(IR),IR = 1,NELC),
     :                          (IROW(IR),IR = 1,NELC)
         DO IV = 1, M
            DIAG = C(ICOL,IV) + (EMT(nelc)-EAV)*B(ICOL,IV)
            CALL DMERGE (NELC-1,B(1,IV),C(1,IV),
     :                  IROW(1),EMT(1),B(ICOL,IV),DL)
            C(ICOL,IV) = DIAG + DL
         ENDDO
      ENDDO

      CALL dalloc (pntemt)   ! AK: Restored, as it is necessary for ifort
      CALL dalloc (pnirow)   ! AK: Restored, as it is necessary for ifort

      CALL gdsummpi (C, N*M)

      RETURN
      END
