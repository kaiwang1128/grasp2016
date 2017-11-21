************************************************************************
*                                                                      *
      SUBROUTINE FINDTYPE(NTYPE,MAP,NCFSMALL,NCFGEN,NCFTOT,MAXSPAN,
     :               NORBGEN)
*                                                                      *
*   Reads the array with occupation numbers and deterimines the type   *
*   of CSF and the ordernumber of the orbitals in ORB_SD-MR            *
*   In the examplification of the span we assume that ORB_GEN          *
*   Is 1s,2s,2p-,2p                                                    *
*                                                                      *
*   Type 1: Ordinary CSF, with no orbital in ORB_SD-MR, e.g. 2s2p      *         
*   Type 2: Symbolic CSF with one orbital in ORB_SD-MR, e.g. 2s7p      *
*           spanning 2s3p, 2s4p, ..., 2s7p                             *
*   Type 3: Symbolic CSF with two orbitals of different symmetry       *
*           in ORB_SD-MR, e.g. 7s7p                                    *
*           spanning 3s3p, 3s4p,..., 3s7p, 4s3p, 4s4p, ...,4s7p,       *
*           7s3p,7s4p,...,7s7p                                         *
*   Type 4: Symbolic CSF with two orbitals of the same symmetry        *
*           in ORB_SD-MR, e.g. 6s7s                                    *
*           spanning 3s4s,3s5s,..,3s7s,4s5s,4s6s,4s7s,5s6s,5s7s,6s7s   *
*   Type 5: SymbolicSF with a doubly occupied orbital in ORB_SD_MR     *
*           e.g. 7s(2)                                                 *
*           spanning 3s(2),4s(2),...,7s(2)                             *         
*                                                                      *
*   NCF    = number of CSFs in the symbolic list                       *         
*   NCFGEN = number of CSFs of type 1 in the symbolic list             *         
*   NCFTOT = total number of CSFs spanned by the CSFs in the           *
*            symbolic list                                             *         
*                                                                      *
*   NTYPE(6,NCF)                                                       *
*                                                                      *
*   NTYPE(1,J) = type of CSF J                                         *
*   NTYPE(2,J) = number of spanned CSFs from symbolic CSF J            *       
*   NTYPE(3,J) = lower position of first  symbolic orb in CSF J        *       
*   NTYPE(4,J) = upper position of first  symbolic orb in CSF J        *         
*   NTYPE(5,J) = lower position of second symbolic orb in CSF J        *       
*   NTYPE(6,J) = upper position of second symbolic orb in CSF J        *         
*                                                                      *
*   Written by Per Jonsson, Malmo University April 2017                *
*                                                                      *
************************************************************************
*
      use symexpand_mod 
      include 'parameters.def'

      INTEGER*4 IQA
      POINTER (PNTRIQ,IQA(NNNWP,1))

      CHARACTER*2 NH

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)      
     :      /ORB10/NH(NNNW)

      INTEGER NTYPE(6,NCF),MAP(NCFSMALL)
    
      NTYPE = 0

* Define the highest principal quantum number in ORB_GEN. 

      WRITE(*,*) 'In findtype, give highest n in ORB_GEN'
!      READ(*,*) NMAX

* Work out the number of orbitals in ORB_GEN assuming that all
* orbitals with n <= nmax are included

* NMAX = 2 --> 1s,2s,2p-2p                  --> NORBGEN = 4      
* NMAX = 3 --> 1s,2s,2p-2p,3s,3p-,3p,3d-,3d --> NORBGEN = 9      
* etc      
      if(nonsymb.eq.4)then
        NMAX = 2
      else if(nonsymb.eq.9)then
        NMAX = 3
      else
       write(*,*)'stop in findtype.f'
       stop
      end if
      SELECT CASE(NMAX)
        CASE(1)
          NORBGEN = 1
        CASE(2)
          NORBGEN = 4
        CASE(3)
          NORBGEN = 9
        CASE(4)
          NORBGEN = 16
        CASE(5)
          NORBGEN = 25
        CASE(6)
          NORBGEN = 36
        CASE(7)
          NORBGEN = 49
        CASE(8)
          NORBGEN = 64
        CASE(9)
          NORBGEN = 81
      END SELECT    

* Loop over CSFs (including added ones) for current block 

      DO J = 1,NCF
        N = 0

* Scan occupation numbers. Put upper orbital positions in element 4 and 6
* We only have to loop over the orbitals not in ORB_GEN
 
        DO I = NORBGEN+1,NW
          IF (IQ(I,J).EQ.2) THEN  ! Doubly occupied, We know that we have type 5
            N = 5                 ! Set N = 5 to flag this
            NTYPE(4,J) = I
            EXIT                  ! We are done with this CSF       
          ELSEIF (IQ(I,J).EQ.1) THEN
            N = N + 1
            NTYPE(2*N+2,J) = I    ! Element 4 for n = 1 and 6 for n = 2
            IF (N.EQ.2) EXIT      ! We are done with this CSF
          END IF
        END DO

* Now compute the type and lower positions. The latter will be in element 3 and 5         

        IF (N.EQ.0) THEN                
          NTYPE(1,J) = 1
        ELSEIF (N.EQ.1) THEN
          NTYPE(1,J) = 2
          NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1
        ELSEIF (N.EQ.2) THEN      ! Two cases: equal or unequal symmetry
          IF (NH(NTYPE(4,J)).EQ.NH(NTYPE(6,J))) THEN  ! Equal
            NTYPE(1,J) = 4
            NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1 
            NTYPE(5,J) = NTYPE(6,J) - NP(NTYPE(6,J)) + NMAX + 2
          ELSE                                        ! Unequal
            NTYPE(1,J) = 3
            NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1 
            NTYPE(5,J) = NTYPE(6,J) - NP(NTYPE(6,J)) + NMAX + 1

          END IF                                        
        ELSEIF (N.EQ.5) THEN
          NTYPE(1,J) = 5
          NTYPE(3,J) = NTYPE(4,J) - NP(NTYPE(4,J)) + NMAX + 1
        END IF
      END DO        

* Count number of CSFs of type 1 as well as the total number of
* CSFs spanned by the symbolic CSFs      

      NCFTOT = 0
      NCFGEN = 0
      MAXSPAN = 0
      DO K = 1,NCFSMALL
        J = MAP(K)
        SELECT CASE (NTYPE(1,J))
          CASE(1)
            NDIFF = 1 
            NCFGEN = NCFGEN + NDIFF
          CASE(2)
            NDIFF = NTYPE(4,J) - NTYPE(3,J) + 1 
          CASE(3)  
            NDIFF = (NTYPE(4,J) - NTYPE(3,J) + 1)*
     :              (NTYPE(6,J) - NTYPE(5,J) + 1)       
          CASE(4) 
            NDIFF = ((NTYPE(6,J) - NTYPE(5,J) + 2)*
     :                         (NTYPE(6,J) - NTYPE(5,J)+1 ))/2
          CASE(5)
            NDIFF = NTYPE(4,J) - NTYPE(3,J) + 1
        END SELECT
        NTYPE(2,J) = NDIFF
        NCFTOT = NCFTOT + NDIFF
        IF (NDIFF.GT.MAXSPAN) THEN
          MAXSPAN = NDIFF
        END IF

      END DO
!      IF(NTYPE(1,J) .EQ. 4)Then
!        NTYPE(5,J)=NTYPE(5,J)+1
!      END IF
      RETURN
      END 
