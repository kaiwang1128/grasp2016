!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     rcsfsymexpand
!
!     This program reads the symbolic list rcsf.inp and adds
!     CSFs needed to complete the calculation of matrix elements
!     between CSFs in symboic blocks
!      
!     WRitten by Per Jonsson, Malmo University, Sweden
!     May 2017
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      subroutine rcsfsymexpand(name)
      implicit none
      character*1500 :: string(5)
      character*300 :: line1,line2,line3,lineadd
      CHARACTER*128 NAME
      character*2 :: string1,string2,string3,string4,sym(21)
      integer :: i,j,nlength,nonsymb,n,keep,nc,m,iunit,norb(21)
      integer :: npos1,npos2,npos3,npos4,ntime

      write(*,*)
      write(*,*) 'rcsfsymexpand'
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

      open (19,file = trim(name)//'.c', status = 'old', 
     :      form = 'formatted')
      open (20,file = trim(name)//'.exp' ,status = 'unknown',
     :      form = 'formatted')

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

      do 
        read(19,'(a)',end=99) line1
        if (line1(1:2).eq.' *') then
          write(20,'(a)') line1
          read(19,'(a)') line1
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
        if (ntime.eq.1) then
          lineadd = line1 
          lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :        string(4)(5*(npos2-2)+1:5*(npos2-1)) 
          write(20,'(a)') trim(lineadd)
          write(20,'(a)') trim(line2)
          write(20,'(a)') trim(line3)
        end if
        if (ntime.eq.2) then
          lineadd = line1 
          if (npos4-1.ne.npos2) then ! Two symmetries e.g. 7s7p

!           Both orbitals one step down e.g. 7s7p --> 6s6p            

            lineadd = line1 
            lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :        string(4)(5*(npos2-2)+1:5*(npos2-1)) 
            lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :        string(4)(5*(npos4-2)+1:5*(npos4-1)) 
            write(20,'(a)') trim(lineadd)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)
 
!           Left orbital one step down e.g. 7s7p --> 6s7p            

            lineadd = line1 
            lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :        string(4)(5*(npos2-2)+1:5*(npos2-1)) 
            write(20,'(a)') trim(lineadd)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)
  
!           Right orbital one step down e.g. 7s7p --> 7s6p            

            lineadd = line1 
            lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :        string(4)(5*(npos4-2)+1:5*(npos4-1)) 
            write(20,'(a)') trim(lineadd)
            write(20,'(a)') trim(line2)
            write(20,'(a)') trim(line3)

          else        ! Same symmetry e.g. 6s7s

!  Check the number of symbolic orbitals with this symmetry

            do j = 1,21
              if (string(4)(5*(npos2-1)+4:5*npos2).eq.sym(j)) then
                n = norb(j)
                exit
              end if
            end do
CDEBUG            write(713,*) trim(line1),npos1,npos2,npos3,npos4,ntime,n

            if (n.gt.3) then

!           Both orbitals two steps down e.g. 6s7s --> 4s5s        

              lineadd = line1 
              lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :          string(4)(5*(npos2-3)+1:5*(npos2-2)) 
              lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :          string(4)(5*(npos4-3)+1:5*(npos4-2)) 
              write(20,'(a)') trim(lineadd)
              write(20,'(a)') trim(line2)
              write(20,'(a)') trim(line3)
            end if

            if (n.gt.2) then

!           Both orbitals one step down e.g. 6s7s --> 5s6s        

              lineadd = line1 
              lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :          string(4)(5*(npos2-2)+1:5*(npos2-1)) 
              lineadd(9*(npos3-1)+1:9*npos3-4) = 
     :          string(4)(5*(npos4-2)+1:5*(npos4-1)) 
              write(20,'(a)') trim(lineadd)
              write(20,'(a)') trim(line2)
              write(20,'(a)') trim(line3)

!           Left orbital one step down e.g. 6s7s --> 5s7s           

              lineadd = line1 
              lineadd(9*(npos1-1)+1:9*npos1-4) = 
     :          string(4)(5*(npos2-2)+1:5*(npos2-1)) 
              write(20,'(a)') trim(lineadd)
              write(20,'(a)') trim(line2)
              write(20,'(a)') trim(line3)
  
            end if  
          end if 
        end if
 
        write(20,'(a)') trim(line1)
        write(20,'(a)') trim(line2)
        write(20,'(a)') trim(line3)
      end do


99    continue

      close(19)
      close(20)

      end
    
