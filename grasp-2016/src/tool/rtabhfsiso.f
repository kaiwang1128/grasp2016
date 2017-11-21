************************************************************************
*                                                                      *
      PROGRAM RTABHFSISO

*                                                                      *
*   This program reads the output from rhfs_lsj and ris_lsj            *
*   and produces LaTeX tables of hfs and isdata                        * 
*                                                                      *
*   Written by Per Jonsson, Malmo University, May 2015                 *
*                                                                      *
************************************************************************


      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,NCOUNT1,NCOUNT2,IEND,JPOS,JPOS1,JPOSNEW,PPOS
      INTEGER APOS,BPOS,NUNIT,INDX(10000),MATCH,NSKIP
      DOUBLE PRECISION A,B,C,ENERGY1(10000),ENERGY2(10000)
      CHARACTER*1 PSTRING1(10000),PSTRING2(10000)
      CHARACTER*4 JSTRING1(10000),JSTRING2(10000)
      CHARACTER*14 ASTRING1(10000),BSTRING1(10000),ASTRINGDUMMY
      CHARACTER*14 ASTRING2(10000),BSTRING2(10000)
      CHARACTER*15 CSTRING1(10000),CSTRING2(10000),CSTRINGDUMMY
      CHARACTER*80 HFSFILE(20),ISFILE(20)
      CHARACTER*500 STRING1(10000),STRING2(10000),STRINGANALYZE
      CHARACTER*64 LATEXSTRING(10000),INPUTSTRING,BLANKSTRING

      WRITE(*,*) ' RTABHFSISO    '
      WRITE(*,*) ' This program reads the output from rhfs_lsj and '
      WRITE(*,*) ' ris_lsj and produces LaTeX tables of hfs and is '
      WRITE(*,*) ' data.'
      WRITE(*,*) ' Input files: name.(c)hlsj, name.(c)ilsj         '
      WRITE(*,*) '    that are files produced by rhfs_lsj and ris_lsj '
      WRITE(*,*) ' Output file: hfsiso.tex '

      OPEN(12,FILE='hfsiso.tex',FORM='FORMATTED',STATUS='UNKNOWN')

*---- Define blank strings for later use ----------------------------

      DO I = 1,14
         ASTRINGDUMMY(I:I) = ' '
      END DO
      DO I = 1,15
         CSTRINGDUMMY(I:I) = ' '
      END DO
      DO I = 1,64
         BLANKSTRING(I:I) = ' '
      END DO

*---- Define input files ---------------------------------------------      

      write(*,*) 
      WRITE(*,'(A)') '  HFS data (1), IS data (2), HFS and IS data (3)'
      READ(*,*) N
      IF (N.EQ.1) THEN
         WRITE(*,'(A)') '  How many HFS files ?'
         READ(*,*) M
         DO I = 1,M
            WRITE(*,*) ' Full name of HFS file',I
            READ(*,'(A)') HFSFILE(I)
            OPEN(40+I,FILE=TRIM(HFSFILE(I)),FORM='FORMATTED',
     :          STATUS = 'OLD')
         END DO
      ELSEIF (N.EQ.2) THEN
         WRITE(*,'(A)') '  How many IS files ?'
         READ(*,*) M
         DO I = 1,M
            WRITE(*,*) ' Full name of IS file',I
            READ(*,'(A)') ISFILE(I)
            OPEN(80+I,FILE=TRIM(ISFILE(I)),FORM='FORMATTED',
     :          STATUS = 'OLD')
         END DO
      ELSE
         WRITE(*,'(A)') '  How many HFS and IS files ?'
         READ(*,*) M
         DO I = 1,M
            WRITE(*,*) '  Full name of HFS file',I
            READ(*,'(A)') HFSFILE(I)
            WRITE(*,*) '  Full name of corresponding IS file',I
            READ(*,'(A)') ISFILE(I)
            OPEN(40+I,FILE=TRIM(HFSFILE(I)),FORM='FORMATTED',
     :          STATUS = 'OLD')
            OPEN(80+I,FILE=TRIM(ISFILE(I)),FORM='FORMATTED',
     :          STATUS = 'OLD')
         END DO
      END IF

      write(*,*) 
      write(*,*) ' Inspect the name.(c)hlsj or name.(c)ilsj files and '
      write(*,*) ' determine how many positions should be skipped in '
      write(*,*) ' the string that determines the label. For example'
      write(*,*) ' if the string is 1s(2).2s_2S.2p(2)3P2_4P and 1s(2) '
      write(*,*) ' is a core then you want to skip 1s(2). i.e. 6'
      write(*,*) ' positions'
      write(*,*) 
      write(*,*) ' How many positions should be skipped?'
      read(*,*) NSKIP

*----- See write up for RIS for units ------------------------------      


*---- Write header to LaTeX table ----------------------------------

      write(12,'(A)') '\documentclass[10pt]{article}'
      write(12,'(A)') '\usepackage{longtable}'
      write(12,'(A)') '\begin{document}'          
      IF (N.EQ.1) THEN
         WRITE(12,'(A)') '\begin{longtable}{lrrrr} \hline'
         WRITE(12,'(A)') 'State & $E$(a.u.) & $A$(MHz) & $B$(MHz) & '
         WRITE(12,'(A)') '$g_J$ \\ \hline'
      ELSEIF (N.EQ.2) THEN
         WRITE(12,'(A)') '\begin{longtable}{lrrrr} \hline'
         WRITE(12,'(A)') 'State & $E$(a.u.) & $K_{NMS}$(a.u.) & '
         WRITE(12,'(A)') '$K_{SMS}$(a.u.) & DENS(a.u.) \\ \hline'
      ELSE
         WRITE(12,'(A)') '\begin{longtable}{lrrrrrrr} \hline'
         WRITE(12,'(A)') 'State & $E$(a.u.) & $K_{NMS}$(a.u.) & '
         WRITE(12,'(A)') '$K_{SMS}$(a.u.) & DENS(a.u.) & $A$(MHz) & '
         WRITE(12,'(A)') '$B$(MHz) & $g_J$ \\ \hline'
      END IF

*---- Start to read data from the files -----------------------------

      IF ((N.EQ.1).OR.(N.EQ.2)) THEN
         J = 1
         IF (N.EQ.1) THEN
            NUNIT = 40
         ELSE
            NUNIT = 80
         END IF

         DO I = 1,M
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            IF (N.EQ.1) THEN
               READ(NUNIT+I,*)
               READ(NUNIT+I,*)
               READ(NUNIT+I,*)
            END IF
            READ(NUNIT+I,'(A)') STRINGANALYZE 
C            WRITE(*,*) TRIM(STRINGANALYZE)
            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'J') THEN
                  JPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'JPOS',JPOS
            JPOS = JPOS - 17

C            JPOS = JPOS - 18 - 4

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'P') THEN
                  PPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'PPOS',PPOS
            PPOS = PPOS - 17

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  APOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'APOS',APOS
            APOS = APOS - 17

            K = K + 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  BPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'BPOS',BPOS
            BPOS = BPOS - 17

            DO 
               READ(NUNIT+I,405,IOSTAT=IEND) ENERGY1(J),STRING1(J)
               PSTRING1(J) = STRING1(J)(PPOS:PPOS)
               JSTRING1(J) = STRING1(J)(JPOS-3:JPOS)
               ASTRING1(J) = STRING1(J)(APOS-13:APOS)
               BSTRING1(J) = STRING1(J)(BPOS-13:BPOS)
               CSTRING1(J) = STRING1(J)(BPOS+2:BPOS+16)
               IF (IEND<0) GOTO 99
C               WRITE(*,*) J,ENERGY1(J),TRIM(STRING1(J))
               J = J + 1
            END DO
99          CONTINUE
         END DO
      
         NCOUNT1 = J - 1

C         WRITE(*,*) NCOUNT1

*---- Now, sort all the energies -----------------------------------

         CALL INDEXS(NCOUNT1,ENERGY1,.FALSE.,INDX)

         JPOS = JPOS - 5

*--- Write energy ordered -----------------------------------------

         DO J = 1,NCOUNT1
            I = INDX(J)
            INPUTSTRING = BLANKSTRING
            JPOSNEW = LEN_TRIM(STRING1(I)(1:JPOS))
            INPUTSTRING(1:JPOSNEW-NSKIP) = STRING1(I)(1+NSKIP:JPOSNEW)
C            INPUTSTRING(1:JPOSNEW) = STRING1(I)(1:JPOSNEW)
            CALL LATEXCONVERT(INPUTSTRING,LATEXSTRING(I))
            LATEXSTRING(I) = '$'//TRIM(LATEXSTRING(I))//
     :             '_{'//JSTRING1(I)
            IF (PSTRING1(I).EQ.'-') THEN
               LATEXSTRING(I) = TRIM(LATEXSTRING(I))//'}^o$'
            ELSE
               LATEXSTRING(I) = TRIM(LATEXSTRING(I))//'}$'
            END IF

            WRITE(12,600) TRIM(LATEXSTRING(I)),' &', ENERGY1(I),' &',
     :        ASTRING1(I),' &',BSTRING1(I),' &', CSTRING1(I),' \\'


         END DO

*----------------------

      ELSE


*---Start to read the iso files -----------------------------------

         J = 1
         NUNIT = 80
         DO I = 1,M
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            IF (N.EQ.1) THEN
               READ(NUNIT+I,*)
               READ(NUNIT+I,*)
               READ(NUNIT+I,*)
            END IF
            READ(NUNIT+I,'(A)') STRINGANALYZE 
C            WRITE(*,*) TRIM(STRINGANALYZE)
            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'J') THEN
                  JPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'JPOS',JPOS
            JPOS = JPOS - 17

C            JPOS = JPOS - 18 - 4

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'P') THEN
                  PPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'PPOS',PPOS
            PPOS = PPOS - 17

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  APOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'APOS',APOS
            APOS = APOS - 17

            K = K + 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  BPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'BPOS',BPOS
            BPOS = BPOS - 17

            DO 
               READ(NUNIT+I,405,IOSTAT=IEND) ENERGY1(J),STRING1(J)
               PSTRING1(J) = STRING1(J)(PPOS:PPOS)
               JSTRING1(J) = STRING1(J)(JPOS-3:JPOS)
               ASTRING1(J) = STRING1(J)(APOS-13:APOS)
               BSTRING1(J) = STRING1(J)(BPOS-13:BPOS)
               CSTRING1(J) = STRING1(J)(BPOS+2:BPOS+16)
               IF (IEND<0) GOTO 990
C               WRITE(*,*) J,ENERGY1(J),TRIM(STRING1(J))
               J = J + 1
            END DO
990         CONTINUE
         END DO

         NCOUNT1 = J - 1

C         WRITE(*,*) NCOUNT1

         JPOS1 = JPOS - 5

*---- Read hfs file ------------------------------------------

         J = 1
         NUNIT = 40
         DO I = 1,M
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            READ(NUNIT+I,*)
            READ(NUNIT+I,'(A)') STRINGANALYZE 
C            WRITE(*,*) TRIM(STRINGANALYZE)
            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'J') THEN
                  JPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'JPOS',JPOS
            JPOS = JPOS - 17

C            JPOS = JPOS - 18 - 4

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.'P') THEN
                  PPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'PPOS',PPOS
            PPOS = PPOS - 17

            K = 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  APOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'APOS',APOS
            APOS = APOS - 17

            K = K + 1
            DO
               IF (STRINGANALYZE(K:K).EQ.')') THEN
                  BPOS = K
                  EXIT
               END IF
               K = K + 1
            END DO
C            WRITE(*,*) 'BPOS',BPOS
            BPOS = BPOS - 17

            DO 
               READ(NUNIT+I,405,IOSTAT=IEND) ENERGY2(J),STRING2(J)
               PSTRING2(J) = STRING2(J)(PPOS:PPOS)
               JSTRING2(J) = STRING2(J)(JPOS-3:JPOS)
               ASTRING2(J) = STRING2(J)(APOS-13:APOS)
               BSTRING2(J) = STRING2(J)(BPOS-13:BPOS)
               CSTRING2(J) = STRING2(J)(BPOS+2:BPOS+16)
               IF (IEND<0) GOTO 991
C               WRITE(*,*) J,ENERGY2(J),TRIM(STRING2(J))
               J = J + 1
            END DO
991         CONTINUE
         END DO

         NCOUNT2 = J - 1

C         WRITE(*,*) NCOUNT2
         

*---- Now, sort all the energies from the iso file ----------------

         CALL INDEXS(NCOUNT1,ENERGY1,.FALSE.,INDX)

*--- Write energy ordered and add hfs data if energies match

         DO J = 1,NCOUNT1
            I = INDX(J)
            INPUTSTRING = BLANKSTRING
            JPOSNEW = LEN_TRIM(STRING1(I)(1:JPOS1))
            INPUTSTRING(1:JPOSNEW-NSKIP) = 
     :                     STRING1(I)(1+NSKIP:JPOSNEW)
C            INPUTSTRING(1:JPOSNEW) = STRING1(I)(1:JPOSNEW)
            CALL LATEXCONVERT(INPUTSTRING,LATEXSTRING(I))
            LATEXSTRING(I) = '$'//TRIM(LATEXSTRING(I))//
     :        '_{'//JSTRING1(I)
            IF (PSTRING1(I).EQ.'-') THEN
               LATEXSTRING(I) = TRIM(LATEXSTRING(I))//'}^o$'
            ELSE
               LATEXSTRING(I) = TRIM(LATEXSTRING(I))//'}$'
            END IF

            MATCH = 0
            DO K = 1,NCOUNT2
               IF (ENERGY1(I).EQ.ENERGY2(K)) THEN
                  MATCH = K
                  EXIT
               END IF
            END DO

            IF (MATCH.GT.0) THEN
               WRITE(12,700) TRIM(LATEXSTRING(I)),' &',ENERGY1(I),' &',
     :          ASTRING1(I),' &',BSTRING1(I),' &',CSTRING1(I),' &', 
     :          ASTRING2(K),' &',BSTRING2(K),' &',CSTRING2(K),' \\'
            ELSE
               WRITE(12,700) TRIM(LATEXSTRING(I)),' &',ENERGY1(I),' &',
     :          ASTRING1(I),' &',BSTRING1(I),' &',CSTRING1(I),' &', 
     :          ASTRINGDUMMY,' &',ASTRINGDUMMY,' &',CSTRINGDUMMY,' \\'
            END IF
         END DO
      END IF

      WRITE(12,'(A)') '\hline\\'
      IF (N.EQ.1) THEN         
         WRITE(12,'(A)') '\caption{Hyperfine interaction constants}'
      ELSEIF (N.EQ.2) THEN
         WRITE(12,'(A)') '\caption{Isotope shift parameters}'
      ELSE
         WRITE(12,'(A)') '\caption{Isotope shift and hfs parameters}'
      END IF
      WRITE(12,'(A)') '\end{longtable}'
      WRITE(12,'(A)') '\end{document}'

  405 FORMAT (1X,F14.7,2X,A) 
  406 FORMAT (1X,F14.7,A) 
  600 FORMAT (1X,2A,1X,F14.7,7A)
  700 FORMAT (1X,2A,1X,F14.7,13A)
C     405 FORMAT (1X,F14.7,2X,A,2A4,1P,3D16.7) 

      END PROGRAM RTABHFSISO

************************************************************************
*                                                                      *
      SUBROUTINE indexS(n,a,ldown,indx)
*                                                                      *
*     Sort out the order of array a and store the index in indx        *
*                                                        (a pointer)   *
*     The input array a is unchanged  written in the bases of UpDown   *
*                                                                      *
*     !$Id: rlevels.f,v 1.2 2003/10/02 07:56:22 per Exp $              *
*     $Log: rlevels.f,v $                                              *
*     Revision 1.2  2003/10/02 07:56:22  per                           *
*     *** empty log message ***                                        *
*                                                                      *
*     Revision 1.1.1.1  2003/01/04 21:45:39  georgio                   *
*                                                                      *
************************************************************************
      IMPLICIT NONE
      CHARACTER*(*)    RCSID
      PARAMETER        ( RCSID
     & ='$Id: rlevels.f,v 1.2 2003/10/02 07:56:22 per Exp $'
     & )
      LOGICAL          ldown              ! .TRUE. then Big ---> Small
      INTEGER          n, indx(n)
      DOUBLE PRECISION a(n)
      INTEGER          i, j, ipos, jpos, jhere
      DOUBLE PRECISION aimx
!
!     Initialize the index array
      DO i = 1, n
         indx(i) = i
      ENDDO
      IF (ldown) THEN
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .GT. aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ELSE
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .LT. aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ENDIF
      RETURN
      END

******************************************************************

      subroutine latexconvert(labelstring,latexstring)

! This subroutine converts a label string to latex
! It is basically the same routine as in renergytable.f90      
! Per Jonssson, Malmo University, November 2014

      implicit none
      integer :: i,j,k,l,ncase
      character(len=64) :: labelstring, dummystring, latexstring
      character(len=1) :: char1, char2, char3

      do i = 1,61

!  Replace (n) with ^n

         if ((labelstring(i:i).eq.'(').and.
     &      (labelstring(i+2:i+2).eq.')')) then
            labelstring(i:i) = '^'
            labelstring(i+2:i+2) = ' '
         end if
      end do

      do i = 1,61

!  Replace . with \,

         if (labelstring(i:i).eq.'.') then
            dummystring = labelstring
            labelstring(1:i-1) = dummystring(1:i-1)
            labelstring(i:i) = '\'
            labelstring(i+1:i+1) = ','
            labelstring(i+2:64) = dummystring(i+1:62) 
         end if
      end do

      do i = 1,61

!  Replace _ with ~

         if (labelstring(i:i).eq.'_') labelstring(i:i) = '~'
      end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and integer2 replace with (^integer1_integer2S), (^integer1_integer2P), etc

      do l = 1,15
         ncase = 0
         do i = 1,61
            do j = 48,57
               do k = 48,57 
                  char1 = labelstring(i:i)
                  char2 = labelstring(i+1:i+1)
                  char3 = labelstring(i+2:i+2)
                  if ((ichar(char1).eq.j).and.(ichar(char3).eq.k).and.
     &               ((char2.ne.'~').and.(char2.ne.' ').and.
     &               (char2.ne.'_'))) then
                     dummystring = labelstring
                     labelstring(1:i-1) = dummystring(1:i-1)
                     labelstring(i:i+6) = 
     &                    '(^'//char1//'_'//char3//char2//')'
                     labelstring(i+7:64) = dummystring(i+3:60)
                     ncase = ncase + 1
                  end if
               end do
            end do
            if (ncase.eq.1) exit
         end do

      end do

!  If integer1 and S, P, D, F, G, H, I, K, L, M, N and not integer2 replace with ^integer1S, ^integer1P, etc 

      do i = 1,61
!  
         if (labelstring(i:i).eq.'~') then
            dummystring = labelstring
            labelstring(1:i) = dummystring(1:i)
            labelstring(i+1:i+1) = '^'
            labelstring(i+2:64) = dummystring(i+1:62)
         end if
      end do

      latexstring = trim(labelstring)

      return
      end      


