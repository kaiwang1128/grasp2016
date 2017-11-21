      real(kind(0.0d0)) function dmach (job) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:31:33   2/12/04  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: job 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: eps, tiny, huge, s 
!-----------------------------------------------
!
!     smach computes machine parameters of floating point
!     arithmetic for use in testing only.  not required by
!     linpack proper.
!
!     if trouble with automatic computation of these quantities,
!     they can be set by direct assignment statements.
!     assume the computer has
!
!        b = base of arithmetic
!        t = number of base  b  digits
!        l = smallest possible exponent
!        u = largest possible exponent
!
!     then
!
!        eps = b**(1-t)
!        tiny = 100.0*b**(-l+t)
!        huge = 0.01*b**(u-t)
!
!     dmach same as smach except t, l, u apply to
!     double precision.
!
!     cmach same as smach except if complex division
!     is done by
!
!        1/(x+i*y) = (x-i*y)/(x**2+y**2)
!
!     then
!
!        tiny = sqrt(tiny)
!        huge = sqrt(huge)
!
!
!     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
!
!
      eps = 1.0D0 
      eps = eps/2.0D0 
      s = 1.0D0 + eps 
      do while(s > 1.0D0) 
         eps = eps/2.0D0 
         s = 1.0D0 + eps 
      end do 
      eps = 2.0D0*eps 
!
      s = 1.0D0 
      tiny = s 
      s = s/16.0D0 
      do while(s*1.0 /= 0.0D0) 
         tiny = s 
         s = s/16.0D0 
      end do 
      tiny = (tiny/eps)*100.0 
      huge = 1.0D0/tiny 
!
      if (job == 1) dmach = eps 
      if (job == 2) dmach = tiny 
      if (job == 3) dmach = huge 
      return  
      end function dmach 
