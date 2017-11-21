************************************************************************
*                                                                      *
        SUBROUTINE genintbreit1wrap (myid, nprocs, j2max)

*   Written by Per Jonsson                Last revision: October 2014  *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      POINTER (PINDT1,INDTP1(1))
      POINTER (PVALT1,VALTP1(1))
      POINTER (PINDT2,INDTP2(1))
      POINTER (PVALT2,VALTP2(1))

      POINTER (PINDT3,INDTP3DUMMY(1))
      POINTER (PVALT3,VALTP3DUMMY(1))
      POINTER (PINDT4,INDTP4DUMMY(1))
      POINTER (PVALT4,VALTP4DUMMY(1))
      POINTER (PINDT5,INDTP5DUMMY(1))
      POINTER (PVALT5,VALTP5DUMMY(1))
      POINTER (PINDT6,INDTP6DUMMY(1))
      POINTER (PVALT6,VALTP6DUMMY(1))

      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)

!-----------------------------------------------------------------------

      CALL genintbreit1 ((myid), (nprocs), N, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      CALL gisummpi (INDTP1, N)
      CALL gdsummpi (VALTP1, N)

      RETURN
      END
