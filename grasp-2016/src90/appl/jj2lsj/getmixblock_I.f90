      MODULE getmixblock_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  18:32:57   1/ 6/07  
      SUBROUTINE getmixblock (NAME, NCI) 
      CHARACTER (LEN = 24), INTENT(IN) :: NAME 
      INTEGER, INTENT(IN) :: NCI 
!VAST.../DEF1/ NELEC(INOUT)
!VAST.../EIGVAL/ EAV(OUT), PNEVAL(INOUT)
!VAST.../EIGVEC/ PNEVEC(INOUT)
!VAST.../ORB2/ NCF(OUT), NW(INOUT)
!VAST.../PRNT/ NVEC(OUT), PNIVEC(INOUT)
!VAST.../SYMA/ PIATJP(INOUT), PIASPA(INOUT)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST...Calls: OPENFL, ALLOC, IVEC
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
