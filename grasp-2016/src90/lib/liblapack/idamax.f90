      integer function idamax (n, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:40   2/12/04  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
!jbq idamax dx(1) -> dx(*) 
      real(double) , intent(in) :: dx(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix 
      real(double) :: dmax 
!-----------------------------------------------
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.
!
!
      idamax = 0 
      if (n < 1) return  
      idamax = 1 
      if (n == 1) return  
      if (incx /= 1) then 
!
!        code for increment not equal to 1
!
         ix = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         dmax = dabs(dx(ix)) 
         ix = ix + incx 
         do i = 2, n 
            if (dabs(dx(ix)) > dmax) then 
               idamax = i 
               dmax = dabs(dx(ix)) 
            endif 
            ix = ix + incx 
         end do 
         return  
      endif 
!
!        code for increment equal to 1
!
      dmax = dabs(dx(1)) 
      do i = 2, n 
         if (dabs(dx(i)) <= dmax) cycle  
         idamax = i 
         dmax = dabs(dx(i)) 
      end do 
      return  
      end function idamax 
