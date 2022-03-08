! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mix the two states with a lambda-based mixing
subroutine mix_frcti(f,ener,fcopy,ecopy,nr3,clambda,klambda)

   use softcore, only : ifsc, mix_frc_sc,log_dvdl, dvdl_norest

   implicit none
   integer nr3,klambda
   _REAL_ f(*),ener(*),fcopy(*),ecopy(*)
   _REAL_ clambda,dvdl,w0,w1

   !     ---we are in final state:   

   w0 = (1.d0 - clambda)**klambda
   w1 = 1.d0 - w0
   if( klambda == 1 )then
     dvdl = ener(23)-ecopy(23)
   else
     dvdl = klambda*(ener(23)-ecopy(23))*(1.d0 - clambda)**(klambda-1)
   end if
      
   if (dvdl_norest /=0) then
      dvdl = dvdl - ecopy(32)
   end if

   ener(1:51) = w1*ener(1:51) + w0*ecopy(1:51)

   if (ifsc == 1) then
      ! This subroutine also adds the softcore contribution to dvdl
      call mix_frc_sc(dvdl,w0,w1,f,fcopy)
   else
      f(1:nr3) = w1*f(1:nr3) + w0*fcopy(1:nr3)
   end if

   ener(39) = dvdl

   ecopy(1:51)  = ener(1:51)

   if (ifsc == 0) then
      fcopy(1:nr3) = f(1:nr3)
   end if

   call log_dvdl(dvdl)      

   return
end subroutine mix_frcti 

