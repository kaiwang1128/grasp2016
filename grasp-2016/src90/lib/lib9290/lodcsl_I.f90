      MODULE lodcsl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  13:07:22   2/14/04  
      SUBROUTINE lodcsl (NCORE) 
      INTEGER NNNW 
      PARAMETER(NNNW=214) 
      INTEGER, INTENT(OUT) :: NCORE 
!VAST.../DEBUGA/ LDBPA(IN)
!VAST.../DEF1/ NELEC(OUT)
!VAST.../ORB2/ NCF(OUT), NW(OUT), PNTRIQ(INOUT)
!VAST.../ORB4/ NP(INOUT), NAK(INOUT)
!VAST.../ORB5/ NKL(INOUT), NKJ(INOUT)
!VAST.../ORB10/ NH(INOUT)
!VAST.../STAT/ PNTJQS(INOUT), PNJCUP(INOUT)
!VAST.../TERMS/ JTAB(IN), NTAB(IN)
!VAST.../IOUNIT/ ISTDE(IN)
!VAST.../BLK/ NBLOCK(OUT), NCFBLK(OUT)
!VAST...Calls: PRSRSL, CONVRT, ALLOC, PRSRCN, PARSJL, RALC2D
!VAST...Calls: PACK, IQ, JQS, JCUP, ITJPO, ISPAR
!...This routine performs I/O.
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
