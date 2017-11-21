!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     RCSFSYM
!
!     This subroutine removes non-symbolic CSFs
!      
!     Let n = 7 be the max principal quantum number for the symbolic 
!     orbitals then 1s(2)6s7s should be kept wheras 1s(2)6s7d, 1s(2)7s6d 
!     should be removed and only 1s(2)7s7d kept      
!     
!     WRitten by Per Jonsson, Malmo University, Sweden
!     May 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      subroutine rcsfsym
      implicit none
      character*1500 :: string(5),stringappend
      character*300 :: line1,line2,line3
      character*2 :: string1,string2,string3,string4,nmax(21),sym(21) 
      character*2 :: nmaxnon(21),nstring
      integer :: i,j,nlength,nonsymb,n,keep,nc,m,iunit
      integer :: nmaxi(21),nmaxnoni(21),n0(21),nmaxcan

      common/nmaxval/nmaxcan

      write(*,*)
      write(*,*) 'Weed out CSFs not needed'
      write(*,*)

      stringappend = repeat(' ',1500)  ! Initialize stringappend
      nmax = '  '
      nmaxnon = '  '

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

      n0(1)  = 1
      n0(2)  = 2
      n0(3)  = 2
      n0(4)  = 3
      n0(5)  = 3
      n0(6)  = 4
      n0(7)  = 4
      n0(8)  = 5
      n0(9)  = 5
      n0(10) = 6
      n0(11) = 6
      n0(12) = 7
      n0(13) = 7
      n0(14) = 8
      n0(15) = 8
      n0(16) = 9
      n0(17) = 9
      n0(18) = 10
      n0(19) = 10
      n0(20) = 11
      n0(21) = 11

      open (19, file = 'rcsf.out', status = 'old', form = 'formatted',
     :          action='readwrite')
      open (1312, status = 'scratch', form = 'formatted')
      open (1313, status = 'scratch', form = 'formatted')


! Read and analyze header to determine the highest principal quantum
! number for each symmetry. Also add orbitals to the header

      do i = 1, 5
         read (19,'(a)') string(i)
         write(*,*) trim(string(i))
      end do

      nlength = len_trim(string(4)) + 1

CDEBUG      write(*,*) trim(string(5))

C      write(*,*) ' Give the number of non-symbolic orbitals'

      select case(nmaxcan)
        case(1)
          nonsymb = 1
        case(2)
          nonsymb = 4
        case(3)
          nonsymb = 9
        case(4)
          nonsymb = 16
        case(5)
          nonsymb = 25
        case(6)
          nonsymb = 36
        case(7)
          nonsymb = 49
        case(8)
          nonsymb = 64
        case(9)
          nonsymb = 81
      end select

! Analyze the non-symbolic orbitals

      do i = 1,nonsymb
         string1 = string(4)(5*(i-1)+4:5*i)   ! symmmetry string
         string2 = string(4)(5*(i-1)+2:5*i-2) ! n string
         do j = 1,21
            if (string1.eq.sym(j)) then
               nmaxnon(j) = string2
            end if
         end do
      end do

CDEBUG      do i = 1,21
CDEBUG        if (nmaxnon(i).ne.'  ') then
CDEBUG          write(*,*) nmaxnon(i),sym(i)
CDEBUG        end if
CDEBUG      end do

! Analyze the orbtitals and determine the maximum n for each symmetry

      do i = nonsymb+1,nlength/5
         string1 = string(4)(5*(i-1)+4:5*i)   ! symmmetry string
         string2 = string(4)(5*(i-1)+2:5*i-2) ! n string
         do j = 1,21
            if (string1.eq.sym(j)) then
               nmax(j) = string2
            end if
         end do
      end do

CDEBGUG      do i = 1,21
CDEBUG        if (nmax(i).ne.'  ') then
CDEBUG          write(*,*) nmax(i),sym(i)
CDEBUG        end if
CDEBUG      end do

   10 continue
      read (19, '(a)', end = 99) line1
CDEBUG      write(*,*) trim(line1)

* Count the number of symbolic orbitals in the CSF 

      n = 0
      do i = nonsymb+1,nlength/5
         if (index(line1,string(4)(5*(i-1)+1:5*i)).gt.0) n = n + 1
      end do

      keep = 1
      if (n.eq.1) then
        nc = 0
        m = len_trim(line1)
        string1 = line1(m-5:m-4)         ! sym of symbolic orbital 1
        string2 = line1(m-7:m-6)         ! n of symbolic orbital 1
        do i = 1,21
          if ((string1.eq.sym(i)).and.(string2.eq.nmax(i))) then
            nc = nc + 1
          end if   
        end do   
        if (nc.lt.1) keep = 0            ! symbolic orbital should have
      elseif (n.eq.2) then               ! max n for keeping
        nc = 0 
        m = len_trim(line1)
        string1 = line1(m-5:m-4)         ! sym of symbolic orbital 1
        string2 = line1(m-7:m-6)         ! n of symbolic orbital 1
        string3 = line1(m-14:m-13)       ! sym of symbolic orbital 2
        string4 = line1(m-16:m-15)       ! n of symbolic orbital 2
        if (string1.ne.string3) then
          do i = 1,21
             if ((string1.eq.sym(i)).and.(string2.eq.nmax(i))) then
                nc = nc + 1
             end if   
             if ((string3.eq.sym(i)).and.(string4.eq.nmax(i))) then
                nc = nc + 1
             end if   
          end do   
          if (nc.lt.2) keep = 0     ! Both symbolic orbitals should
        end if                      ! have max n for keeping
      end if                     
            
      read (19, '(a)') line2
      read (19, '(a)') line3
CDEBUG      write(*,*) trim(line2)
CDEBUG      write(*,*) trim(line3)

      if (keep.eq.1) then
         if (n.eq.0) then 
            iunit = 1312
         else
            iunit = 1313
         end if
         write(iunit,'(a)') trim(line1)
         write(iunit,'(a)') trim(line2)
         write(iunit,'(a)') trim(line3)
      end if

      goto 10

  99  continue

      rewind(19)
      rewind(1312)
      rewind(1313)


      do i = 1,3
         write(19,'(a)') trim(string(i))
      end do

* Do some manipulations to write out all orbitals on line 4

      do i = 1,21
         read(nmax(i),'(i2)') nmaxi(i) ! Convert nmax string to integer
         read(nmaxnon(i),'(i2)') nmaxnoni(i)
CDEBUG         write(*,*) sym(i),nmaxi(i),nmaxnoni(i)
      end do

      do i = 1,21
        if (nmaxi(i).gt.0) then
          do j = max0(n0(i),nmaxnoni(i)+1),nmaxi(i)
CDEBUG            write(*,*) j,sym(i)
            write(nstring,'(i2)') j  ! Convert integer j to string
            if (sym(i)(2:2).eq.'') sym(i)(2:2) = '+'
            stringappend = trim(stringappend)//' '//nstring//sym(i)
          end do
        end if
      end do

! Remove + from the string

      do i = 1,1500
        if (stringappend(i:i).eq.'+') stringappend(i:i) = ' '
      end do
      
CDEBUG      write(*,'(a)') trim(stringappend)

! Now write line 4  and 5

      write(19,'(a)') string(4)(1:5*nonsymb)//trim(stringappend)
      write(19,'(a)') string(5)

! Start by writing all the non-symbolic CSFs

   11 continue
      read (1312,'(a)', end = 998) line1
      read (1312,'(a)') line2
      read (1312,'(a)') line3
      write(19,'(a)') trim(line1)
      write(19,'(a)') trim(line2)
      write(19,'(a)') trim(line3)
      go to 11  
  998 continue     

! Here comes the symbolic CSFs

   12 continue
      read (1313,'(a)', end = 999) line1
      read (1313,'(a)') line2
      read (1313,'(a)') line3
      write(19,'(a)') trim(line1)
      write(19,'(a)') trim(line2)
      write(19,'(a)') trim(line3)
      go to 12  
  999 continue     


      close(19)

      end
    
