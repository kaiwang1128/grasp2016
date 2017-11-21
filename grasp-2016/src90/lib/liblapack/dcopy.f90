      subroutine  dcopy(nSIZE,dx,incx,dy,incy)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:28   2/12/04
!...Switches:
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nSIZE
      integer , intent(in) :: incx
      integer , intent(in) :: incy
!jbq dcopy dx(1) dy(1) -> dx(*) dy(*)
!     real(double) , intent(in) :: dx(1)
!     real(double) , intent(out) :: dy(1)
      real(double) , intent(in) :: dx(*)
      real(double) , intent(out) :: dy(*)
!     double precision dx(1),dy(1)
!      integer i,incx,incy,ix,iy,m,mp1,n
!jbq
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1
!-----------------------------------------------
!
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      if(nSIZE.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-nSIZE+1)*incx + 1
      if(incy.lt.0)iy = (-nSIZE+1)*incy + 1
      do 10 i = 1,nSIZE
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(nSIZE,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( nSIZE .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,nSIZE,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
