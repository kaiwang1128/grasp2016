      subroutine drotg(da, db, c, s) 
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
      real(double) , intent(inout) :: da 
      real(double) , intent(inout) :: db 
      real(double) , intent(out) :: c 
      real(double) , intent(out) :: s 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: roe, scale, r, z 
!-----------------------------------------------
!
!     construct givens plane rotation.
!     jack dongarra, linpack, 3/11/78.
!                    modified 9/27/86.
!
!
      roe = db 
      if (dabs(da) > dabs(db)) roe = da 
      scale = dabs(da) + dabs(db) 
      if (scale == 0.0D0) then 
         c = 1.0D0 
         s = 0.0D0 
         r = 0.0D0 
      else 
         r = scale*dsqrt((da/scale)**2 + (db/scale)**2) 
         r = dsign(1.0D0,roe)*r 
         c = da/r 
         s = db/r 
      endif 
      z = s 
      if (dabs(c)>0.0D0 .and. dabs(c)<=s) z = 1.0D0/c 
      da = r 
      db = z 
      return  
      end subroutine drotg 
