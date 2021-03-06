! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.h"

subroutine selfs(i,j,r2,iovpar,em,ne)
   !implicit _REAL_ (a-h,o-z)
   implicit none

   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,1)
   _REAL_ :: parm(13)
   integer :: id
   integer, external :: iovpar  
!   dimension parm(13),em(ne,1)
  ! external iovpar

   id=iovpar(i,j,r2,parm)
   em(1,1)=parm(13)
   return
end subroutine selfs


subroutine selfp(i,j,r2,iovpar,em,ne)
   implicit _REAL_ (a-h,o-z)
   dimension parm(13),em(ne,3)
   external iovpar

   id=iovpar(i,j,r2,parm)
   do 11 l=1,3
      do 10 m=1,3
         em(l,m)=0.0
10    enddo
      em(l,l)=parm(12)
11 enddo
   return
end subroutine selfp


subroutine selfd(i,j,r2,iovpar,em,ne)
   implicit _REAL_ (a-h,o-z)
   dimension parm(13),em(ne,5)
   external iovpar

   id=iovpar(i,j,r2,parm)
   do 11 l=1,5
      do 10 m=1,5
         em(l,m)=0.0
10    enddo
      em(l,l)=parm(11)
11 enddo
   return
end subroutine selfd
