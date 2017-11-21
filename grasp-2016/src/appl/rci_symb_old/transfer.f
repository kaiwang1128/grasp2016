***********************************************************************
*                                                                     *      
      SUBROUTINE TRANSFER(A,IAROW,IACOUNT,B,M1,M2,N,IROWACCUM)
*                                                                     *
*  This subroutine transfers the non-zero data in B(M2,N) along       *
*  with the corresponding ACCUMULATED row position to the arrays      *
*  A(M1,N) and IROW(M1,N).                                            *
*  IACOUNT(N) is a counter array that keeps track of                  *
*  the number of non-zero elements in A(M1,N) for each column         *
*                                                                     *
*  For a documentation see xxxx                                       *         
*                                                                     *
*  Per Jönsson, Malmö University,                April 2017           *
*                                                                     *
***********************************************************************         

      IMPLICIT NONE   
      REAL*8 A(M1,N), B(M2,N)
      INTEGER I,J,M1,M2,N,IROWACCUM
      INTEGER IAROW(M1,N),IACOUNT(N)

      DO I = 1,N 
        DO J = 1,M2
          IF (B(J,I).NE.0.D0) THEN
            IACOUNT(I) = IACOUNT(I) + 1
            A(IACOUNT(I),I) = B(J,I)
            IAROW(IACOUNT(I),I) = J + IROWACCUM
          END IF
        END DO
      END DO

      IROWACCUM = IROWACCUM + M2

      RETURN
      END




