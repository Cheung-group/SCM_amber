! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdwrit here]
subroutine mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb,x,v,xc, box,tt,temp0)

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : ncsu_on_mdwrit => on_mdwrit
#endif /* DISABLE_NCSU */

   implicit none
   integer nstep,nrp,nr,nres,ntxo,ntr,ntb
   _REAL_ x,v,xc,box,tt,temp0
   character(len=89) restrt2
   character(len=12) num
   integer istart,iend
   logical first
   save first
   data first/.true./
#  include "files.h"
   
   !     -- open/write/close the restrt file:

#ifndef DISABLE_NCSU
   call ncsu_on_mdwrit()
#endif /* DISABLE_NCSU */
   
   if( first ) then
      if (ntxo == 0) then
         call amopen(16,restrt,owrite,'U','W')
      else
         call amopen(16,restrt,owrite,'F','W')
      end if
      first = .false.
   else
      if (ntxo == 0) then
         call amopen(16,restrt,'O','U','W')
      else
         call amopen(16,restrt,'O','F','W')
      end if
   end if
   call mdwri2(16,nrp,nr,nres,ntxo,ntr,ntb,x,v,xc, &
         box,tt,temp0)
   close(16)
   
   !     -- consider whether to save secondary restrt;
   
   if (ntwr >= 0) return
   
   do iend=1,80
      if (restrt(iend:iend) <= ' ') goto 1
   end do
   1 continue
   iend = iend - 1
   write(num,'(i12)') nstep
   do istart=1,12
      if (num(istart:istart) /= ' ') goto 2
   end do
   2 continue
   write(restrt2, '(a,a,a)') restrt(1:iend), '_', num(istart:12)
   write(6,'(a,a)') ' writing ', restrt2
   if (ntxo == 0) then
      call amopen(17,restrt2,owrite,'U','W')
   else
      call amopen(17,restrt2,owrite,'F','W')
   end if
   call mdwri2(17,nrp,nr,nres,ntxo,ntr,ntb,x,v,xc, &
         box,tt,temp0)
   close(17)
   return
end subroutine mdwrit 
!------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdwri2 here]
subroutine mdwri2(nf,nrp,nr,nres,ntxo,ntr,ntb,x,v,xc, box,tt,temp0)
   use nblist, only: a,b,c,alpha,beta,gamma
#ifdef MPI
   use remd, only : rem
#endif
   implicit none
   integer nf,nrp,nr,nres,ntxo,ntr,ntb
   
   !     ----- ROUTINE TO WRITE FINAL COORDINATES AND VELOCITIES -----
   
   _REAL_ x(*),v(*),xc(*),box(3),tt,temp0
   integer nr3,i,i3,j,j3
#  include "files.h"
   
   nr3 = 3*nr
   if(ntxo /= 0) then
      
      !     ----- FORMATTED WRITING -----
      
      write(nf,9008) title

#ifdef MPI
      if(rem > 0) then
         if( nr > 99999 ) then
           write(nf,9019) nr,tt,temp0
         else
           write(nf,9018) nr,tt,temp0
         end if
      else
#endif
         if( nr > 99999 ) then
           write(nf,9019) nr,tt
         else
           write(nf,9018) nr,tt
         end if
#ifdef MPI
      end if
#endif
      write(nf,9028) (x(i),i=1,nr3)
      write(nf,9028) (v(i),i=1,nr3)

      if ( ntb /= 0 ) write(nf,9028) a,b,c,alpha,beta,gamma
      
   else
      
      !     ----- BINARY WRITING -----
      write(nf) title
      write(nf) nr,tt
      write(nf) (x(i),i = 1,nr3)
      write(nf) (v(i),i = 1,nr3)
      if ( ntb /= 0 ) write(nf) a,b,c,alpha,beta,gamma
   end if
   
   9008 format(a80)
   9018 format(i5,2e15.7)
   9019 format(i6,2e15.7)
   9028 format(6f12.7)
   return
end subroutine mdwri2 
