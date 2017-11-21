*-----------------------------------------------------------------------
*
*     RPCFGENERATE -- A PROGRAM FOR GENERATING PAIR CORRELATION
*                     FUNCTION EXPANSIONS      
*
*     Per Jönsson
*     Malmö University
*     S-20506 Malmö, Sweden
*
*     This program takes a SD-MR list and partitions it into a 
*     sequence of PCFs that each includes the CSFs in the MR. 
*     In this version of the code all valence electrons are treated
*     together. The program generates a sequence of PCFs based on
*     relativistic orbital pairs but in the end these are merged to
*     PCFs that adheres to anon-relativistic notation. 
*     For example, the program internally generates PCFs 2p-2p, 2p-2p, 2p2p 
*     etc but they are merged to a common 2p2p that covers all three
*     relativistic PCFs at the end.
*
*     September 2016
*
*-----------------------------------------------------------------------
*
      program rpcfgenerate 
      implicit double precision (a-h,o-z)
      integer     qval(20,10000,120),qcore(40),qorb(120),qorbcsf(120)
      integer     q(16),nval(20),ncsfmr(20),nbounds(20,2)
      integer     mrpos(10000),nfile,nduplicate(1600),nblock
      integer     qvalorb(20,98),nfilepointer(200)
      character*1 ans
      character*4  val(20,120),core(40),orb(120),elc(16),valorb(98)
      character*3  occupation(4)
      character*24 name,label,namenonrel(1600)
      character*300 line1,line2,line3,line4,line5,corestring
      character*300 line1save,line2save,line3save,line4save,line5save
      character*300 mrline1(10000),mrline2(10000),mrline3(10000)
      data valorb/
     :  ' 1s ',' 2s ',' 2p-',' 2p ',' 3s ',' 3p-',' 3p ',' 3d-',' 3d ',
     :  ' 4s ',' 4p-',' 4p ',' 4d-',' 4d ',' 4f-',' 4f ',' 5s ',' 5p-', 
     :  ' 5p ',' 5d-',' 5d ',' 5f-',' 5f ',' 5g-',' 5g ',' 6s ',' 6p-',
     :  ' 6p ',' 6d-',' 6d ',' 6f-',' 6f ',' 6g-',' 6g ',' 6h-',' 6h ', 
     :  ' 7s ',' 7p-',' 7p ',' 7d-',' 7d ',' 7f-',' 7f ',' 7g-',' 7g ',
     :  ' 7h-',' 7h ',' 7i-',' 7i ',' 8s ',' 8p-',' 8p ',' 8d-',' 8d ',
     :  ' 8f-',' 8f ',' 8g-',' 8g ',' 8h-',' 8h ',' 8i-',' 8i ',' 8k-',
     :  ' 8k ',' 9s ',' 9p-',' 9p ',' 9d-',' 9d ',' 9f-',' 9f ',' 9g-',
     :  ' 9g ',' 9h-',' 9h ',' 9i-',' 9i ',' 9k-',' 9k ',' 9l-',' 9l ',
     :  '10s ','10p-','10p ','10d-','10d ','10f-','10f ','10g-','10g ', 
     :  '10h-','10h ','10i-','10i ','10k-','10k ','10l-','10l '/
      data occupation/'(2)','(4)','(6)','(8)'/


      write(*,*)
      write(*,*) ' Welcome to the RPCFGENERATE program'
      write(*,*) 
      write(*,*) ' A configuration list rcsf.inp that contains CSFs  '
      write(*,*) ' from an rcsfgenerate SD-MR run combined with      '
      write(*,*) ' angular reduction using rcsfinteract should be    '
      write(*,*) ' available. The CSFs in rcsf.inp are partitioned   '
      write(*,*) ' into PCFs based on the available MR given in      '
      write(*,*) ' rcsfmr.inp. In the current version all outer      '
      write(*,*) ' valence orbitals are treated together.            '
      write(*,*) ' A selected number of PCFs can be concatenated to  '
      write(*,*) ' a file concat.c                                   '
      write(*,*) '                                                   '
      write(*,*) ' Current restrictions: 120  orbitals               '
      write(*,*) '                     10000  reference CSFs         '
      write(*,*) '                        16  subshells              '
      write(*,*) '                        20  blocks                 '
      write(*,*) '                        10l highest orb in MR      '
      write(*,*) '                                                   '
      write(*,*) ' Input files: rcsfmr.inp, rcsf.inp                 '
      write(*,*) ' Output files: a number of PCF files with          '
      write(*,*) ' typical names 1s1snonrel.c, 1s2snonrel.c etc      '
      write(*,*) ' along with the file concat.c                      '
      write(*,*)

* --- Open configuration list file rcsl.inp  

      open(unit=7,file='rcsf.inp',status='old')

* --- Open multireference list file mrlist

      open(unit=8,file='rcsfmr.inp',status='old')

* --- Check that rcsf.inp has been reduced

      write(*,*) ' Has rcsf.inp been reduced using rcsfinteract (y/n) ?'      
      read(*,'(a)') ans
      if ((ans.eq.'N').or.(ans.eq.'n')) then
        write(*,*) ' rcsf.inp must be reduced using rcsfinteract else'
        write(*,*) ' the partition in PCFs can not be done'
        STOP
      end if

* --- Read initial lines of configuration list

      read(7,'(a)') line1save
      read(7,'(a)') line2save
      read(7,'(a)') line3save
      read(7,'(a)') line4save
      read(7,'(a)') line5save

* --- Define valence occupation for multireference

      read(8,'(a)') line1
      read(8,'(a)') line2
      read(8,'(a)') line3
      read(8,'(a)') line4
      read(8,'(a)') line5

* --- Make sure we have the same core

      if (line2save.ne.line2) then 
        write(*,*) ' The two lists should have the same core'
        stop
      end if

* --- Read label for output files

      write(*,*) ' Specify label for output files'
      read(5,'(a)') label

* --- Define core orbitals

      write(*,*) ' Specify the active core in relativistic notation '
      write(*,*) ' e.g. 2s(2)2p-(2)2p(4)                            '
      read(*,'(a)') corestring

* --- Analyze corestring, each orbital should be saved in an a4 string 
*     and centered in such a way that the l-quantum number is in
*     position 3, e.g. ' 2p-',' 2p '

      corestring = adjustl(corestring)
      nlength = len_trim(corestring)
      ncore = 0

* --- The loop over orbitals needs to be divided in two parts dependent
*     on how many digits the principal quantum number has

      do i = 1,81
        nfind = index(corestring,trim(valorb(i)(2:4))//'(')  
        if (nfind.gt.0) then
          ncore = ncore + 1
          core(ncore) = valorb(i)
        end if
      end do

      do i = 82,98
        nfind = index(corestring,trim(valorb(i)(1:4))//'(')  
        if (nfind.gt.0) then
          ncore = ncore + 1
          core(ncore) = valorb(i)
        end if
      end do

* --- Now find the occupation of the core orbitals      

      ncore = 0
      do i = 1,nlength-2
        do j = 1,4
          if (corestring(i:i+2).eq.occupation(j)) then
            ncore = ncore + 1
            qcore(ncore) = 2*j
          end if
        end do
      end do

* --- Write out for debug

      write(*,*)
      write(*,*) ' Label and occupation for core orbitals' 
      do i = 1,ncore
        write(*,'(2x,a,i3)') core(i),qcore(i)
      end do

* --- Define valence orbitals
   
      qvalorb = 0 
      nblock = 1
      do 
        q = 0
        elc = '   '

* --- Handle block separators

        read(8,'(a)',end=7) line1
        if (line1(2:2).eq.'*') then
          read(8,'(16(1x,a4,1x,i2,1x))',end=7)  (elc(k),q(k),k=1,16)
          nblock = nblock + 1
        else
          backspace 8
          read(8,'(16(1x,a4,1x,i2,1x))',end=7)  (elc(k),q(k),k=1,16)
        end if
        backspace 8
        read(8,'(a)') line1
        read(8,'(a)') line2
        read(8,'(a)',end=7) line3
        do j = 1,98
          do k = 1,16
            if (valorb(j).eq.elc(k)) then
              qvalorb(nblock,j) = 1
            end if
          end do
        end do
      end do

    7 continue

      do i = 1,nblock
        nval(i) = 0
        do j = 1,98
          if ((qvalorb(i,j).eq.1)) then
            nfound = 0
            do l = 1,ncore 
              if (valorb(j).eq.core(l)) then
                nfound = 1
              end if
            end do 
            if (nfound.eq.0) then
              nval(i) = nval(i) + 1
              val(i,nval(i)) = valorb(j)
            end if
          end if
        end do
      end do 

* --- Debug write out 

      write(*,*) 
      write(*,*) ' Valence orbitals for the different blocks'
      do i = 1,nblock
        do j = 1,nval(i)
          write(*,'(a,i3,a,a)') '  Block',i, ' valence orbital ',
     :       val(i,j)
C          write(*,*) ' Block, valence orbital',i,val(i,j)
        end do
      end do    

* --- Define valence occupation for multireference      

* --- Loop over CSFs in MR

      rewind(8)

      read(8,'(a)') line1
      read(8,'(a)') line2
      read(8,'(a)') line3
      read(8,'(a)') line4
      read(8,'(a)') line5

      qval = 0
      mr = 0
      nblock = 1
      ncsfmr = 0
      do 
        q = 0
        elc = '   '

* --- Handle block separators

        read(8,'(a)',end=5) line1
        if (line1(2:2).eq.'*') then
          read(8,'(16(1x,a4,1x,i2,1x))',end=5)  (elc(k),q(k),k=1,16)
          nblock = nblock + 1
        else
          backspace 8
          read(8,'(16(1x,a4,1x,i2,1x))',end=5)  (elc(k),q(k),k=1,16)
        end if
        backspace 8
        read(8,'(a)') line1
        read(8,'(a)') line2
        read(8,'(a)',end=5) line3
        ncsfmr(nblock) = ncsfmr(nblock) + 1
        mr = mr + 1
        mrline1(mr) = line1
        mrline2(mr) = line2
        mrline3(mr) = line3
        do j = ncore+1,16
          do k = 1,nval(nblock)
            if (elc(j).eq.val(nblock,k)) qval(nblock,mr,k) = q(j)
          end do
        end do
      end do

    5 continue

* --- Derive bounds for the CSFs in the different MR blocks

      nbounds(1,1) = 1
      nbounds(1,2) = ncsfmr(1)
      do k = 2,nblock
         nbounds(k,1) = nbounds(k-1,2)+1
         nbounds(k,2) = nbounds(k,1) + ncsfmr(k) - 1
      end do

* --- Debug print

C      write(*,*)
C      do k = 1,nblock
C         write(*,*) ' Block bounds in MR',k,nbounds(k,1),nbounds(k,2)
C      end do      
C      write(*,*)

* --- Write out for debug

C      do k = 1,nblock
C        write(*,*) ' Block ',k,'number of CSFs in MR',ncsfmr(k) 
C        write(*,*)
C        do i = nbounds(k,1),nbounds(k,2)
C          write(*,*) ' Label and occupation for valence orbitals in CSF'
C     :       ,i 
C          do j = 1,nval(k)
C            write(*,*) val(k,j),qval(k,i,j)
C          end do
C        end do 
C        write(*,*)
C      end do

* --- Open relativistic PCFs involving core core

      do i = 1,ncore
        do j = i,ncore
          name = trim(adjustl(core(i)))//trim(adjustl(core(j)))//
     :           trim(label)//'.c'
          open(unit=10+i+40*j,file=name,status='unknown')
        end do
      end do

* --- Open relativistic PCFs involving core valence and valence valence

      do i = 1,ncore
        name = trim(adjustl(core(i)))//'v'//trim(label)//'.c'
        open(unit=10000+i,file=name,status='unknown')
      end do    

      name = 'vv'//trim(label)//'.c'
      open(unit=20000,file=name,status='unknown')

* --- Loop over CSFs from configuration list and classify to which pair

      ncsf = 0
      nblock = 1
      do 
        q = 0
        elc = '   '

* --- Handle block separators

        read(7,'(a)',end=10) line1
        if (line1(2:2).eq.'*') then
          read(7,'(16(1x,a4,1x,i2,1x))',end=10)  (elc(k),q(k),k=1,16)
          nblock = nblock + 1
        else
          backspace 7
          read(7,'(16(1x,a4,1x,i2,1x))',end=10)  (elc(k),q(k),k=1,16)
        end if
        backspace 7
        read(7,'(a)') line1
        read(7,'(a)') line2
        read(7,'(a)',end=10) line3
        ncsf = ncsf + 1

* --- Check if CSFs in multireference. If so skip and we will write this
*     later

        inmr = 0
        do i = 1,mr
          if ((line1.eq.mrline1(i)).and.(line2.eq.mrline2(i)).and.
     :        (line3.eq.mrline3(i))) then
            inmr = 1
            goto 100
          end if
        end do

* --- Loop over CSFs in multireference in the same block

        do i = nbounds(nblock,1),nbounds(nblock,2)

* --- Define orbtals and occupation for current reference CSF

          norb = ncore + nval(nblock)
          do j = 1,norb
            if (j.le.ncore) then
              orb(j) = core(j)
              qorb(j) = qcore(j)
            else
              orb(j) = val(nblock,j-ncore)
              qorb(j) = qval(nblock,i,j-ncore)
            end if
          end do

* --- Match orbitals and set ccupation for current CSF

          qorbcsf = 0
          do j = 1,16
            do k = 1,norb
              if (elc(j).eq.orb(k)) then
                qorbcsf(k) = q(j)
              end if
            end do
          end do


* --- Loop over orbital pairs from which the excitations are done

          do j = 1,norb
            do k = j,norb

* --- Determine which case
*      case 1: excitations from the same core shell
*      case 2: excitations from two different core shells
*      case 3: excitations from a core shell and from the valence
*      case 4: one excitation from the valence shell
*      case 5: two excitations from the valence shell

              ncase = 0
              if ((k.eq.j).and.(j.le.ncore)) ncase=1
              if ((k.gt.j).and.((k.le.ncore).and.(j.le.ncore))) ncase=2
              if ((j.le.ncore).and.(k.gt.ncore)) ncase=3
              if ((k.eq.j).and.(j.gt.ncore)) ncase=4
              if ((k.gt.j).and.(j.gt.ncore)) ncase=5

* --- Process the different cases

              if (ncase.eq.1) then

* --- Occupation of j in CSF one or two less than in reference

                if (qorbcsf(j).eq.(qorb(j)-1)) then

                  n1 = 0
                  do l = 1,ncore
                    if (l.ne.j) then
                      if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                    end if
                  end do

                  n2 = 0
                  n3 = 0
                  do l = ncore+1,norb
                    if (qorbcsf(l).ne.qorb(l)) then 
                      n2 = n2 + 1
                      if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99

                else if (qorbcsf(j).eq.(qorb(j)-2)) then

                  n1 = 0
                  do l = 1,ncore
                    if (l.ne.j) then
                      if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                    end if
                  end do

                  n2 = 0
                  n3 = 0
                  n4 = 0
                  do l = ncore+1,norb
                    if (qorbcsf(l).ne.qorb(l)) then 
                      n2 = n2 + 1
                      if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                      if (qorbcsf(l).eq.(qorb(l)+2)) n4 = n4 + 1
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n4.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.2).and.(n3.eq.2))) goto 99

                end if

              else if (ncase.eq.2) then

                if ((qorbcsf(j).eq.(qorb(j)-1)).and.
     :              (qorbcsf(k).eq.(qorb(k)-1))) then

                  n1 = 0
                  do l = 1,ncore
                    if ((l.ne.j).and.(l.ne.k)) then
                      if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                    end if
                  end do

                  n2 = 0
                  n3 = 0
                  n4 = 0
                  do l = ncore+1,norb
                    if (qorbcsf(l).ne.qorb(l)) then 
                      n2 = n2 + 1
                      if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                      if (qorbcsf(l).eq.(qorb(l)+2)) n4 = n4 + 1
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n4.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.2).and.(n3.eq.2))) goto 99

                end if

              else if (ncase.eq.3) then

                if ((qorbcsf(j).eq.(qorb(j)-1)).and.
     :              (qorbcsf(k).eq.(qorb(k)-1))) then

                  n1 = 0
                  do l = 1,ncore
                    if (l.ne.j) then
                      if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                    end if
                  end do

                  n2 = 0
                  n3 = 0
                  n4 = 0
                  do l = ncore+1,norb
                    if (l.ne.k) then  
                      if (qorbcsf(l).ne.qorb(l)) then 
                        n2 = n2 + 1
                        if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                        if (qorbcsf(l).eq.(qorb(l)+2)) n4 = n4 + 1
                      end if
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n4.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.2).and.(n3.eq.2))) goto 99

                end if

              else if (ncase.eq.4) then

                if (qorbcsf(j).eq.(qorb(j)-1)) then

                  n1 = 0
                  do l = 1,ncore
                    if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                  end do

                  n2 = 0
                  n3 = 0
                  do l = ncore+1,norb
                    if (l.ne.j) then
                      if (qorbcsf(l).ne.qorb(l)) then 
                        n2 = n2 + 1
                        if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                      end if
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99

                else if (qorbcsf(j).eq.(qorb(j)-2)) then

                  n1 = 0
                  do l = 1,ncore
                    if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                  end do

                  n2 = 0
                  n3 = 0
                  n4 = 0
                  do l = ncore+1,norb
                    if (l.ne.j) then
                      if (qorbcsf(l).ne.qorb(l)) then 
                        n2 = n2 + 1
                        if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                        if (qorbcsf(l).eq.(qorb(l)+2)) n4 = n4 + 1
                      end if
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n4.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.2).and.(n3.eq.2))) goto 99

                end if

              else if (ncase.eq.5) then

                if ((qorbcsf(j).eq.(qorb(j)-1)).and.
     :              (qorbcsf(k).eq.(qorb(k)-1))) then

                  n1 = 0
                  do l = 1,ncore
                    if (qorbcsf(l).ne.qorb(l)) n1 = n1 + 1
                  end do

                  n2 = 0
                  n3 = 0
                  n4 = 0
                  do l = ncore+1,norb
                    if ((l.ne.j).and.(l.ne.k)) then
                      if (qorbcsf(l).ne.qorb(l)) then 
                        n2 = n2 + 1
                        if (qorbcsf(l).eq.(qorb(l)+1)) n3 = n3 + 1
                        if (qorbcsf(l).eq.(qorb(l)+2)) n4 = n4 + 1
                      end if
                    end if
                  end do

                  if ((n1.eq.0).and.(n2.eq.0)) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n3.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.1).and.(n4.eq.1))) goto 99
                  if ((n1.eq.0).and.((n2.eq.2).and.(n3.eq.2))) goto 99

                end if

              end if 
            
            end do
          end do
        end do

   99   continue

* --- Write CSF to correct LCF file

        if ((j.le.ncore).and.(k.le.ncore)) then
          write(10+j+40*k,'(a)') trim(line1)
          write(10+j+40*k,'(a)') trim(line2)
          write(10+j+40*k,'(a)') trim(line3)
        end if

        if ((j.le.ncore).and.(k.gt.ncore)) then
          write(10000+j,'(a)') trim(line1)
          write(10000+j,'(a)') trim(line2)
          write(10000+j,'(a)') trim(line3)
        end if

        if ((j.gt.ncore).and.(k.gt.ncore)) then
          write(20000,'(a)') trim(line1)
          write(20000,'(a)') trim(line2)
          write(20000,'(a)') trim(line3)
        end if

  100   continue

* --- Write CSFs in MR to all PCF files

C        if (inmr.eq.1) then 
C          do i = 1,ncore
C            do j = i,ncore
C              write(10+i+40*j,'(a)') trim(line1) 
C              write(10+i+40*j,'(a)') trim(line2) 
C              write(10+i+40*j,'(a)') trim(line3) 
C            end do
C          end do

C          do i = 1,ncore
C            write(10000+i,'(a)') trim(line1) 
C            write(10000+i,'(a)') trim(line2) 
C            write(10000+i,'(a)') trim(line3) 
C          end do    

C          write(20000,'(a)') trim(line1) 
C          write(20000,'(a)') trim(line2) 
C          write(20000,'(a)') trim(line3) 
C        end if

      end do

   10 continue


* ---- Start by clearing any old non-relativistic PCF files that might exist

      nfile = 0
      do i = 1,ncore
        do j = i,ncore
          nfile = nfile + 1
          name = trim(adjustl(core(i)(1:3)))
     :      //trim(adjustl(core(j)(1:3)))//'nonrel'//trim(label)//'.c'
          namenonrel(nfile) = name
          open(unit=99999,file=name,status='unknown')
          close(99999,status='delete')
        end do
      end do

      do i = 1,ncore
        nfile = nfile + 1 
        name = trim(adjustl(core(i)(1:3)))//'v'
     :    //'nonrel'//trim(label)//'.c'
        namenonrel(nfile) = name
        open(unit=99999,file=name,status='unknown')
        close(99999,status='delete')
      end do    

      nfile = nfile + 1
      name = 'vvnonrel'//trim(label)//'.c'
      namenonrel(nfile) = name
      open(unit=99999,file=name,status='unknown')
      close(99999,status='delete')

C --- Add the first five lines and the MR to each PCF file      

      nduplicate = 0
      do i = 1,nfile
        do j = i+1,nfile
          if (trim(namenonrel(i)).eq.trim(namenonrel(j))) then
            nduplicate(j) = 1
          end if
        end do
      end do

      write(*,*)
      write(*,*) ' Generated PCF files'

      do i = 1,nfile
        if (nduplicate(i).eq.0) then
          write(*,'(a,i3,a,a)') '  (',i,')  ',namenonrel(i) 
          open(unit=99999,file=trim(namenonrel(i)),status='unknown',
     :     access='append')
          write(99999,'(a)') line1save
          write(99999,'(a)') line2save
          write(99999,'(a)') line3save
          write(99999,'(a)') line4save
          write(99999,'(a)') line5save
          rewind(8) 
C --- Read MR header
          read(8,'(a)') line1
          read(8,'(a)') line2
          read(8,'(a)') line3
          read(8,'(a)') line4
          read(8,'(a)') line5
          do 
             read(8,'(a)',end=81) line1
             if (line1(2:2).eq.'*') then
               read(8,'(a)') line1
             end if
             read(8,'(a)') line2
             read(8,'(a)') line3
             write(99999,'(a)') trim(line1)
             write(99999,'(a)') trim(line2)
             write(99999,'(a)') trim(line3)
          end do
81        continue          

          close(99999)
        end if

      end do

C --- Count how many lines the header and MR occupies (will be used
C     later)      

      rewind(8) 
      read(8,'(a)') line1
      read(8,'(a)') line2
      read(8,'(a)') line3
      read(8,'(a)') line4
      read(8,'(a)') line5
      nheader = 5
      do 
        read(8,'(a)',end=83) line1
        if (line1(2:2).eq.'*') then
          read(8,'(a)') line1
        end if
        read(8,'(a)') line2
        read(8,'(a)') line3
        nheader = nheader + 3
      end do
83    continue          

C      write(*,*) 'Number of lines in header and MR',nheader

      write(*,*) 
      write(*,*) ' Number of PCF files you want to concatenate'
      read(*,*) nconcat
      if (nconcat.gt.0) then
        write(*,'(a,i4,a)') ' Order numbers of the ',nconcat,
     :    ' PCF files to concatenate' 
        read(*,*) (nfilepointer(i),i=1,nconcat)
      end if
C      do i = 1,nconcat
C        write(*,*) 'Number of file',i
C        read(*,*) nfilepointer(i)
C      end do

* --- Now open non-relativistic PCFs involving core in append mode and
*     write content to non-relativistic PCF

      do i = 1,ncore
        do j = i,ncore
          name = trim(adjustl(core(i)(1:3)))
     :      //trim(adjustl(core(j)(1:3)))//'nonrel'//trim(label)//'.c'
          open(unit=99999,file=name,status='unknown',access='append')
          rewind(10+i+40*j)
          do 
            read(10+i+40*j,'(a)',end=59) line1
            write(99999,'(a)') trim(line1)
          end do
59        continue 
          close(99999)
        end do
      end do

C --- Non-relativistic PCFs involving core valence and valence valence in append mode

      do i = 1,ncore
        name = trim(adjustl(core(i)(1:3)))//'v'
     :    //'nonrel'//trim(label)//'.c'
        open(unit=99999,file=name,status='unknown',access='append')
        rewind(unit=10000+i)
        do 
          read(10000+i,'(a)',end=69) line1
          write(99999,'(a)') trim(line1)
        end do
69      continue 
        close(99999)
      end do    

      name = 'vvnonrel'//trim(label)//'.c'
      open(unit=99999,file=name,status='unknown',access='append')
      rewind(unit=20000)
      do 
        read(20000,'(a)',end=79) line1
        write(99999,'(a)') trim(line1)
      end do
79    continue 
      close(99999)

C --- Open concatenated file in append mode

      if (nconcat.gt.0) then
        open(unit=99999,file='concat.c',status='unknown')
        close(99999,status='delete')

       open(unit=99999,file='concat.c',status='unknown',access='append')
        do i = 1,nconcat
          open(unit=99999+i,file=namenonrel(nfilepointer(i)),
     :      status='unknown')

C --- Make sure that we write the header and the MR only once in the
C    concatenated list

          if (i.ne.1) then
            do j = 1,nheader
              read(99999+i,*)
            end do
          end if
          do 
            read(99999+i,'(a)',end=73) line1
            write(99999,'(a)') trim(line1)
          end do
73        continue 
          close(99999+i)
        end do
        close(99999)
      end if

C --- Finally put the lists on block form

      write(*,*) 
      write(*,*) 'Shown are the accumulated number of CSFs as the'
      write(*,*) 'number of PCFs increase. Note that the CSFs in the'
      write(*,*) 'MR will be counted for every PCF'
      write(*,*) 

      do i = nfile,1,-1
        if (nduplicate(i).eq.0) then
          write(*,*) 'Calling block for     ',namenonrel(i)
          call rcsfblock(namenonrel(i))
          write(*,*)
          write(*,*)
        end if
      end do

C --- If concatenated list put it on block form

      if (nconcat.gt.0) then
        write(*,*) 'Calling block for concat.c'
        call rcsfblock('concat.c')
      end if

* --- Close and delete relativistic PCFs

      do i = 1,ncore
        close(30000+i,status='delete')
      end do

      do i = 1,ncore
        do j = i,ncore
          close(10+i+40*j,status='delete')
        end do
      end do

      do i = 1,ncore
        close(10000+i,status='delete')
      end do    

      close(20000,status='delete')

      stop
      end
