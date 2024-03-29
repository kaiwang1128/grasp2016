************************************************************************
*                                                                      *
      SUBROUTINE CANONICAL
*                                                                      *
*     Sorts the orbitals in canonical order                            *
*                                                                      *
*   Written by Per Jonsson, March 2014                                 *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
*
      CHARACTER*2 NH,NHTMP(NNNW)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
      INTEGER NPTMP(NNNW),NAKTMP(NNNW),NKLTMP(NNNW),NKJTMP(NNNW),
     :        MAP(NNNW),NPC(1000),NAKC(1000)
*     :        MAP(NNNW),NPC(NNNW),NAKC(NNNW)

*
*  Orbitals in canonical order
*      1s   
*      2s   2p-  2p   
*      3s   3p-  3p   3d-  3d   
*      4s   4p-  4p   4d-  4d   4f-  4f   
*      5s   5p-  5p   5d-  5d   5f-  5f   5g-  5g 
*      6s   6p-  6p   6d-  6d   6f-  6f   6g-  6g   6h-  6h  
*      7s   7p-  7p   7d-  7d   7f-  7f   7g-  7g   7h-  7h   7i-  7i
*      8s   8p-  8p   8d-  8d   8f-  8f   8g-  8g   8h-  8h   8i-  8i   8k-  8k
*      9s   9p-  9p   9d-  9d   9f-  9f   9g-  9g   9h-  9h   9i-  9i   9k-  9k   9l-  9l
*     10s  10p- 10p  10d- 10d  10f- 10f  10g- 10g  10h- 10h  10i- 10i  10k- 10k  10l- 10l  10m- 10m
*     11s  11p- 11p  11d- 11d  11f- 11f  11g- 11g  11h- 11h  11i- 11i  11k- 11k  11l- 11l  11m- 11m 
*     12s  12p- 12p  12d- 12d  12f- 12f  12g- 12g  12h- 12h  12i- 12i  12k- 12k  12l- 12l  12m- 12m
*     13s  13p- 13p  13d- 13d  13f- 13f  13g- 13g  13h- 13h  13i- 13i  13k- 13k  13l- 13l  13m- 13m
*     14s  14p- 14p  14d- 14d  14f- 14f  14g- 14g  14h- 14h  14i- 14i  14k- 14k  14l- 14l  14m- 14m
*     15s  15p- 15p  15d- 15d  15f- 15f  15g- 15g  15h- 15h  15i- 15i  15k- 15k  15l- 15l  15m- 15m
*      
*  Kappa
*      -1    1   -2    2   -3    3   -4    4   -5    5   -6    6   -7    7   -8    8   -9    9  -10

*  We arrange the canonical orbitals in order

      N = 0
*      DO I = 1,15
      DO I = 1,30
         DO J = 1,MIN(2*I-1,19)
            N = N + 1
            NPC(N) = I
            IF (MOD(J,2).EQ.0) THEN
               NAKC(N) = J/2
            ELSE
               NAKC(N) = -(J+1)/2
            END IF
         END DO
      END DO

*      DO I = 1,N
*         WRITE(*,*) I,NPC(I),NAKC(I)
*      END DO

*      WRITE(*,*) N

*      WRITE(*,*) 'NP,NAK,NKL,NKJ,NH'
*      DO I = 1,NW
*         WRITE(*,*) NP(I),NAK(I),NKL(I),NKJ(I),NH(I)
*      END DO

*      PAUSE

*  Map the orbitals on the canonical orbitals

      M = 0
      DO I = 1,N
         DO J = 1,NW
            IF ((NP(J).EQ.NPC(I)).AND.(NAK(J).EQ.NAKC(I))) THEN
               M = M + 1
               MAP(M) = J
            END IF
         END DO
      END DO

      DO I = 1,NW
         NPTMP(I)  = NP(MAP(I))
         NAKTMP(I) = NAK(MAP(I))
         NKLTMP(I) = NKL(MAP(I))
         NKJTMP(I) = NKJ(MAP(I))
         NHTMP(I)  = NH(MAP(I))
      END DO

      DO I = 1,NW
         NP(I)  = NPTMP(I)
         NAK(I) = NAKTMP(I)
         NKL(I) = NKLTMP(I)
         NKJ(I) = NKJTMP(I)
         NH(I)  = NHTMP(I)
      END DO


*      WRITE(*,*) 'Final order'

*      DO I = 1,NW
*         WRITE(*,*) NP(I),NAK(I),NKL(I),NKJ(I),NH(I)
*      END DO

*      pause

      END
