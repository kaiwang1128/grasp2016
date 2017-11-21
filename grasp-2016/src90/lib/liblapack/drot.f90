      subroutine drot(n, dx, incx, dy, incy, c, s) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:34   2/12/04  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(in) :: c 
      real(double) , intent(in) :: s 
!jbq drot dx(1) dy(1) -> dx(*) dy(*)
      real(double) , intent(inout) :: dx(*) 
      real(double) , intent(inout) :: dy(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1 
         iy = 1 
         if (incx < 0) ix = ((-n) + 1)*incx + 1 
         if (incy < 0) iy = ((-n) + 1)*incy + 1 
         do i = 1, n 
            dtemp = c*dx(ix) + s*dy(iy) 
            dy(iy) = c*dy(iy) - s*dx(ix) 
            dx(ix) = dtemp 
            ix = ix + incx 
            iy = iy + incy 
         end do 
         return  
      endif 
!
!       code for both increments equal to 1
!
      do i = 1, n 
         dtemp = c*dx(i) + s*dy(i) 
         dy(i) = c*dy(i) - s*dx(i) 
         dx(i) = dtemp 
      end do 
      return  
      end subroutine drot 
