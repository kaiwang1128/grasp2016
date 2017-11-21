      subroutine testmodule 
      use symexpand_mod
      implicit none
      integer :: i,j

      write(*,*)
      write(*,*) 'in testmodule'
      write(*,*)

      write(*,*) ' Number of non-symbolic orbitals',nonsymb
      write(*,*) 'nblocks',nblock
      write(*,*) 'ncsfsymperblock'
      do i = 1,nblock
         write(*,*) ncsfsymperblock(i)
      end do

      do j = 1,nblock
        do i = 1,ncsfsymperblock(j)
           write(*,*) i,j,map(i,j)
        end do
      end do


      end
    
