!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     rcsfsymexpandfull
!
!     This program reads the symbolic list rcsf.inp and adds
!     CSFs to obtain a full normal list for checking and validation
!      
!     WRitten by Per Jonsson, Malmo University, Sweden
!     May 2017
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      program rcsfsymexpandfull
      use symexpand_mod
      implicit none
      character*1500 :: string(5)
      character*300 :: line1,line2,line3,lineadd
      character*2 :: string1,string2,string3,string4,sym(21)
      integer :: i,j,nlength,n,keep,nc,m,iunit,norb(21)
      integer :: npos1,npos2,npos3,npos4,ntime,n1,n2
      integer :: maxcsf,ncount1,ncount2

      write(*,*)
      write(*,*) 'rcsfsymexpandfull'
      write(*,*)

! Intitialize symmetry and number

      sym(1)  = 's '                                                        
      sym(2)  = 'p-'
      sym(3)  = 'p '
      sym(4)  = 'd-'
      sym(5)  = 'd '
      sym(6)  = 'f-'
      sym(7)  = 'f '
      sym(8)  = 'g-'
      sym(9)  = 'g '
      sym(10) = 'h-'
      sym(11) = 'h '
      sym(12) = 'i-'
      sym(13) = 'i '
      sym(14) = 'k-'
      sym(15) = 'k '
      sym(16) = 'l-'
      sym(17) = 'l '
      sym(18) = 'm-'
      sym(19) = 'm '
      sym(20) = 'n-'      
      sym(21) = 'n '      

      norb = 0

! Open files    

      open (19,file = 'rcsf.inp', status = 'old', form = 'formatted')
      open (20,file = 'rcsf.exp' ,status = 'unknown',form = 'formatted')

! Start by looping through the file to count the number of blocks and
! CSFs in each block

      do i = 1, 5
         read(19,'(a)') string(i)
      end do

      ncsfsymperblock = 0
      nblock = 1
      do 
        read(19,'(a)',end=9) line1
        if (line1(1:2).eq.' *') then
          nblock = nblock + 1 
          read(19,'(a)') line1
        end if
        read(19,'(a)') line2  
        read(19,'(a)') line3  
        ncsfsymperblock(nblock) = ncsfsymperblock(nblock) + 1
      end do

9     continue

! Allaocate space for map

      do i = 1,nblock
         write(*,*) 'Block',i,'number of CSFs',ncsfsymperblock(i)
      end do
      maxcsf = maxval(ncsfsymperblock)
      write(*,*) 'maxcsf',maxcsf
      allocate(map(maxcsf,nblock))

! Rewind

      rewind(19)

! Read and analyze header

      do i = 1, 5
         read(19,'(a)') string(i)
         write(20,'(a)') trim(string(i))
      end do

      nlength = len_trim(string(4)) + 1

      write(*,*) trim(string(4))

! Determine the number of canonical orbitals by looking for 
! adjacent orbitals with the same symmetry (excluding 1s)

      string1 = '  '
      do i = 2,nlength/5
        string2 = string(4)(5*(i-1)+4:5*i)   ! symmmetry string
        if (string1.eq.string2) then         ! adjacent symmetry
          nonsymb = i - 2 
          exit
        end if
        string1 = string2
      end do

      write(*,*) ' Number of non-symbolic orbitals'
      write(*,*) nonsymb

! Loop through the symbolic orbitals and deterine the number of orbitals for
! each symmetry

      do i = nonsymb+1,nlength/5
        do j = 1,21
          if (string(4)(5*(i-1)+4:5*i).eq.sym(j)) then
            norb(j) = norb(j) + 1
          end if
        end do
      end do

      do i = 1,21
        write(*,*) norb(i),sym(i) 
      end do

! Loop through the CSFs and determine the type

      nblock = 1
      ncount1 = 0
      ncount2 = 0
      do 
        read(19,'(a)',end=99) line1
        if (line1(1:2).eq.' *') then
          write(20,'(a)') line1
          read(19,'(a)') line1
          nblock = nblock + 1
          ncount1 = 0
          ncount2 = 0
        end if
        read(19,'(a)') line2  
        read(19,'(a)') line3  

! Loop through the symbolic orbitals

        npos1 = 0
        npos2 = 0
        npos3 = 0
        npos4 = 0
        ntime = 0
        do i = 1,len_trim(line1)/9
          do j = nonsymb+1,nlength/5
            if (line1(9*(i-1)+1:9*i-4).eq.string(4)(5*(j-1)+1:5*j)) then
              if (ntime.eq.0) then
                npos1 = i
                npos2 = j 
                ntime = ntime + 1
              else if (ntime.eq.1) then  
                npos3 = i
                npos4 = j 
                ntime = ntime + 1
              end if
            end if
          end do
        end do 
C        write(*,*) npos1,npos2,npos3,npos4,ntime

!  Write CSFs of type 1

        if (ntime.eq.0) then
          write(20,'(a)') trim(line1)
          write(20,'(a)') trim(line2)
          write(20,'(a)') trim(line3)
        end if
      
        if (ntime.eq.1) then  ! CSF of type 2 or 5

!  Check the number of symbolic orbitals with this symmetry

          do j = 1,21
            if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
              n = norb(j)
              exit
            end if
          end do

!  Replacements 7s    --> 4s,5s,6s,7s           
!               7s(2) --> 4s(2),5s(2),6s(2),7s(2)            

!  Singe loop over symbolic orbitals to do the replacement

          lineadd = line1 
          do j = 1,n
            lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :        string(4)(5*(npos2-(n-j)-1)+1:5*(npos2-(n-j))) 
C     :        string(4)(5*(npos2-2)+1:5*(npos2-1)) 
            write(20,'(a)') trim(lineadd)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)
            ncount2 = ncount2 + 1
          end do
        end if

        if (ntime.eq.2) then ! Type 3,4
          lineadd = line1 
          if (npos4-1.ne.npos2) then ! Two symmetries e.g. 7s7p

!  Check the number of symbolic orbitals with these two symmetries

            do j = 1,21
              if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
                n1 = norb(j)
                exit
              end if
            end do

            do j = 1,21
              if (string(4)(5*(npos4-1)+4:5*npos4).eq.sym(j)) then
                n2 = norb(j)
                exit
              end if
            end do

!  Replacements 7s7p --> 4s4p,4s5p,4s6p,4s7p,5s4p,5s5p,5s6p,5s7,          
!                        ..... 7s6p,7s7p

!  Double loop over symbolic orbitals to do the replacement

            lineadd = line1 
            do i = 1,n1
              lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :          string(4)(5*(npos2-(n1-i)-1)+1:5*(npos2-(n1-i))) 
              do j = 1,n2 
                lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :            string(4)(5*(npos4-(n2-j)-1)+1:5*(npos4-(n2-j))) 
                write(20,'(a)') trim(lineadd)
                write(20,'(a)') trim(line2)
                write(20,'(a)') trim(line3)
                ncount2 = ncount2 + 1
              end do
            end do
 
          else        ! Same symmetry e.g. 6s7s

!  Check the number of symbolic orbitals with this symmetry

            do j = 1,21
              if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
                n = norb(j)
            write(713,*) 'Hello',sym(j),n
                exit
              end if
            end do

!  Replacements 6s7s --> 4s5s,4s6s,4s7s,5s6s,5s7s,6s7s

!  Double loop over symbolic orbitals to do the replacement

            lineadd = line1 

!  Please note that npos4 is the postion of 7s and that npos2 that
!  of 6s. We do the distribution in relation to npos4           

            do i = 1,n-1
              lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :          string(4)(5*(npos4-(n-i)-1)+1:5*(npos4-(n-i))) 
              do j = i+1,n 
                lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :            string(4)(5*(npos4-(n-j)-1)+1:5*(npos4-(n-j))) 
                write(20,'(a)') trim(lineadd)
                write(20,'(a)') trim(line2)
                write(20,'(a)') trim(line3)
                ncount2 = ncount2 + 1
              end do
            end do

          end if 
        end if
 
        ncount1 = ncount1 + 1

        map(ncount1,nblock) = ncount2

      end do


99    continue

! This is for testing the module      
C      call testmodule


      end
    
