************************************************************************
*                                                                      *
      SUBROUTINE ALLOC (PTR,NMLOCS,LENGTH)
*                                                                      *
*   This  routine allocates  NMLOCS memory locations, each of length   *
*   LENGTH bytes to the pointee of the pointer PTR.   LENGTH is  not   *
*   required for the  CFT77  version --- the  8-byte word  length is   *
*   assumed common to all data types; it has been retained to ensure   *
*   compatibility with the code that calls this subprogram.            *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This version by Farid A. Parpia.      Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      * 
*  Modified by C.F. Fischer
*  This version of the code has the default integer size (I4 or I8)    *
*  A new version has NMLOCS  as I8 even though the default may be I4   *
************************************************************************
*
cbieron DEC pointer
c
c     INTEGER PTR                 ! DEC and CRAY require ?
      pointer (PTR,ptrdummy)
      INTEGER NMLOCS, LENGTH
      INTEGER*8 NBYTES
      INTEGER*8 MALLOC,MMM
c
*
*   Compute the size of the memory region in bytes; this is retained for
*   calling compatibility and debugging purposes in the Cray version
*   Also, all allocations are multiples of 8 bytes for alignment
*   purposes.
*
cgd   NBYTES = NMLOCS*LENGTH
      NBYTES = NMLOCS
	  NBYTES = NBYTES*LENGTH
      IF (MOD(NBYTES,8) .NE. 0) NBYTES = 8*(NBYTES/8 + 1)
*
      IF (NBYTES. LE. 0) THEN
         PRINT *, 'ALLOC-i4: Invalid memory request:'
         PRINT *, ' PTR = ',PTR,', NMLOCS = ',NMLOCS,
     :               ', LENGTH = ',LENGTH,'.'
         STOP
*
      ELSE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         PTR = MMM (NBYTES)                ! SUN & DEC
!         PTR = MALLOC(%VAL(NBYTES))        ! IBM
!         IF (PTR .EQ. 0) THEN
!
!         IABORT = 0
!         CALL HPALLOC(PTR,NMLOCS,IERR,IABORT) ! CRAY 
!         IF (ierr.ne.0) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Choose one from above

         PTR = MMM (NBYTES)                ! SUN & DEC
         IF (PTR .EQ. 0) THEN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            PRINT *, 'ALLOC-i4: Unable to allocate memory:'
            PRINT *, ' PTR = ',PTR,', NMLOCS =',NMLOCS,
     :                  ', LENGTH = ',LENGTH,
     :                  '.'
            STOP
         ENDIF
      ENDIF
      RETURN
      END

************************************************************************
*                                                                      *
      SUBROUTINE ALLOC2 (PTR,NMLOCS,LENGTH)
*                                                                      *
*   This  routine allocates  NMLOCS memory locations, each of length   *
*   LENGTH bytes to the pointee of the pointer PTR.                    *
*   In this version it is assumed the argument is an I8 integer.       *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This version by Farid A. Parpia.      Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*   Modified by C. F. Fischer                            29 Apr 2016   * 
************************************************************************
*
      pointer (PTR,ptrdummy)
      INTEGER LENGTH
      INTEGER*8 NBYTES,  NMLOCS
      INTEGER*8 MALLOC,MMM
c
*
*   Compute the size of the memory region in bytes;

CFF   change the following to avoid integer expressions of diff. types
*     NBYTES = NMLOCS
*     NBYTES = NBYTES*LENGTH
      NBYTES = LENGTH
      NBYTES = NBYTES*NMLOCS
      IF (MOD(NBYTES,8) .NE. 0) NBYTES = 8*(NBYTES/8 + 1)
*
      IF (NBYTES. LE. 0) THEN
         PRINT *, 'ALLOC-i8: Invalid memory request:'
         PRINT *, ' PTR = ',PTR,', NMLOCS = ',NMLOCS,
     :               ', LENGTH = ',LENGTH,'.',
     :               ', NBYTES = ',NBYTES,'.'
         STOP
*
      ELSE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         PTR = MMM (NBYTES)                ! SUN & DEC
!         PTR = MALLOC(%VAL(NBYTES))        ! IBM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Choose one from above

         PTR = MMM (NBYTES)                ! SUN & DEC
         IF (PTR .EQ. 0) THEN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            PRINT *, 'ALLOC-i8: Unable to allocate memory:'
            PRINT *, ' PTR = ',PTR,', NMLOCS =',NMLOCS,
     :                  ', LENGTH = ',LENGTH, 
     :                  'NBYTES=', NBYTES, '.'
            STOP
         ENDIF
      ENDIF
      RETURN
      END

