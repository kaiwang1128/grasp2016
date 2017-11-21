      subroutine dscal(n, da, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:35   2/12/04  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      real(double) , intent(in) :: da 
!jbq dscal dx(1) -> dx(*) 
      real(double) , intent(inout) :: dx(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, m, mp1 
!-----------------------------------------------
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.
!
!
      if (n <= 0) return  
      if (incx /= 1) then 
!
!        code for increment not equal to 1
!
         ix = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         dx(ix:(n-1)*incx+ix:incx) = da*dx(ix:(n-1)*incx+ix:incx) 
         return  
      endif 
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
      m = mod(n,5) 
      if (m /= 0) then 
         dx(:m) = da*dx(:m) 
         if (n < 5) return  
      endif 
      mp1 = m + 1 
      dx(mp1:((n-mp1+5)/5)*5-1+mp1) = da*dx(mp1:((n-mp1+5)/5)*5-1+mp1) 
      return  
      end subroutine dscal 
