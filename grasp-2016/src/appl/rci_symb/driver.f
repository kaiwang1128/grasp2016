***********************************************************************
*                                                                     *      
      PROGRAM DRIVER
         
*  Per Jönsson, Malmö University,                April 2017           *
*                                                                     *
***********************************************************************         

      REAL*8 EA(20,20), EB(2,10)
      INTEGER I,J,IROWACCUM
      INTEGER IEAROW(20,20),IEACOUNT(20)

      EA = 0.D0
      EB = 0.D0
      IEAROW = 0
      IACOUNT = 0
      IROWACCUM = 0

      EB(1,1) = 1.D0
      EB(1,2) = 2.D0
      EB(1,3) = 0.D0
      EB(2,1) = 0.D0
      EB(2,2) = 4.D0
      EB(2,3) = 0.D0

*     EB = [1 2 0
*           0 4 0]      

      CALL TRANSFER(EA(1:20,1:3),IEAROW(1:20,1:3),
     :        IEACOUNT(1:3),EB(1:2,1:3),20,2,3,IROWACCUM)

      WRITE(*,*) 'First result'

      DO I = 1,3
        DO J = 1,IEACOUNT(I)
           WRITE(*,*) 'J,I,EA(J,I),IEAROW(J,1)',J,I,EA(J,I),
     :                                           IEAROW(J,I)
        END DO
      END DO

      DO I = 1,3
         WRITE(*,*) 'I,IEACOUNT(I)',IEACOUNT(I)
      END DO

      WRITE(*,*) 

      EB(1,1) = -1.D0
      EB(1,2) = 0.D0
      EB(1,3) = 5.D0

*     EB = [-1 0 5 ]

      WRITE(*,*) 'Second result'
      CALL TRANSFER(EA(1:20,1:3),IEAROW(1:20,1:3),
     :        IEACOUNT(1:3),EB(1:1,1:3),20,1,3,IROWACCUM)

      DO I = 1,3
         WRITE(*,*) 'I,IEACOUNT(I)',IEACOUNT(I)
      END DO

      WRITE(*,*) 

      DO I = 1,3
        DO J = 1,IEACOUNT(I)
           WRITE(*,*) 'J,I,EA(J,I),IEAROW(J,1)',J,I,EA(J,I),
     :                                           IEAROW(J,I)
        END DO
      END DO

      END




