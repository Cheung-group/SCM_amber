! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bond here]
subroutine bond(nbin,ib,jb,icb,x,xx,ix,f,eb)
   

   !     ----- ROUTINE TO GET BOND ENERGY AND FORCES FOR THE POTENTIAL
   !           OF CB*(B-B0)**2
   !------------------------------------------------------------------------

   use decomp, only : decpr, decpair
   use parms , only: req, rk
   use poisson_boltzmann, only : outflagorig

!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  EVB modules                                                     ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#ifdef LES
   use pimd_vars, only: nrg_all,nbead,nbead_inv,ipimd,ineb
#  ifdef MPI
      use remd, only : rem
#  endif
#endif
   use pimd_vars, only: bnd_vir
#ifdef MPI /* SOFT CORE */
   use softcore, only: ifsc, nsc, sc_ener, oneweight
#endif
   implicit none
#ifdef MPI
#  include "parallel.h"
#endif
#  include "nmr.h"
#  include "md.h"
#  include "files.h"
   !-------------passed-in variables  ---------------------------
   integer nbin
   integer ib(*),jb(*),icb(*),ix(*)
   _REAL_ x(*),xx(*),f(*),eb
  
   
#  include "box.h"
#  include "flocntrl.h"
#ifdef LES
#  include "les.h"
#endif
   
   !------------- local variables  ---------------------------
   integer max190
   parameter (max190=190)
   integer nb,ist,maxlen,ksum,i3,j3,ic,jn,ii,jj
   _REAL_ xij(max190),yij(max190),zij(max190),eaw(max190), &
         rij(max190),dfw(max190)
   _REAL_ da,df,xa,ya,za,rij0,ebl
   logical skip
   integer piece,start,end,newnb

!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  EVB variables                                                   ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   integer :: ndx

#ifndef MPI
   integer numtasks,mytaskid
   
   numtasks=1
   mytaskid=0
#endif

   if ( do_bond == 0 .or. nbin == 0 ) return
   
   if(numtasks > 1)then
      piece = nbin/numtasks
      start = mytaskid*piece+1
      end   = mytaskid*piece+piece
      if(mytaskid == (numtasks-1)) end = nbin
   else
      start=1
      end=nbin
   end if
   nb = end
   ist = start -1
   
   ebl = 0.0d+00
#ifdef LES
   elesb = 0.0d+00
#endif
   
   !     ----- GRAND LOOP FOR THE bond STUFF -----
   
   4200 continue
   maxlen = max190
   skip = (ist+maxlen) > nb
   if (skip) maxlen = nb-ist
   if (maxlen > 0) then
      
      do jn = 1,maxlen
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)
         
         !           ----- CALCULATION OF THE bond vector -----
         
         xij(jn) = x(i3+1)-x(j3+1)
         yij(jn) = x(i3+2)-x(j3+2)
         zij(jn) = x(i3+3)-x(j3+3)
      end do
      
      do jn = 1,maxlen
         rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
         rij(jn) = sqrt(rij0)
      end do
      
      !         ----- CALCULATION OF THE ENERGY AND DER -----
      
      do jn = 1,maxlen
         ii = (ib(jn+ist) + 3)/3
         jj = (jb(jn+ist) + 3)/3

         ic = icb(jn+ist)
         rij0 = rij(jn)
         da = rij0-req(ic)
         if((ifcap == 2 .or. ifcap == 5))then
            if( (outflagorig(ii) == 1 .or. outflagorig(jj) == 1)) then
               da = 0.0d0
            end if
         endif
         !  for rms deviation from ideal bonds:
         ebdev = ebdev + da*da
         df = rk(ic)*da
         eaw(jn) = df*da
#ifdef MPI /* SOFT CORE */
         ! For dual-topology softcore runs, bonds involving sc atoms are modified here
         if (ifsc /= 0) then
            ! Check if a softcore atom is involved in this bond
            if (nsc(ii) == 1 .or. nsc(jj) == 1) then
               ! This bond needs to 
               ! a) get its energy removed from Ebond
               ! b) get its force scaled up by 1/weight
               sc_ener(1) = sc_ener(1) + eaw(jn)
               eaw(jn) = 0.0d0
               df = df * oneweight
            end if
         end if
         if(decpr .and. (idecomp == 1 .or. idecomp == 2)) then
            if(ifsc == 0 ) then
               call decpair(4,ii,jj,eaw(jn)/(nstlim/ntpr))
            else if (nsc(ii) /= 1 .and. nsc(jj) /= 1) then
               call decpair(4,ii,jj,eaw(jn)/(nstlim/ntpr))
            end if
         end if
#else
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(4,ii,jj,eaw(jn))
         end if
#endif
#ifdef MPI
#  ifdef LES
         if(rem == 2) then
           if(cnum(ii) > 0 .or. cnum(jj) > 0) then
              elesb = elesb + eaw(jn)
           endif
         end if
#  endif
#endif

#ifdef LES
!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!!
      if(ipimd>0.or.ineb>0) then
         if( cnum(ii) == 0 .and. cnum(jj) == 0 ) then
             nrg_all(1:nbead) = nrg_all(1:nbead) + eaw(jn) * nbead_inv 
         else
             if( cnum(ii) == 0 ) then 
                ndx = cnum(jj)
             else
                ndx = cnum(ii)
             end if

             nrg_all(ndx) = nrg_all(ndx) + eaw(jn)
         end if
      end if
!!
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#endif /* LES */

!! Amber harmonic bond interaction for comparison with morsify ''''''''''''''''''''''
!     if( jn == 1 ) then                                                    ! DEBUG ;
!        write(6,'(A)') '|'                                                 ! DEBUG ;
!        write(6,'(A)') '| ^^^ Amber harmonic energy ///////////////////'   ! DEBUG ;
!        write(6,'(A)') '|'                                                 ! DEBUG ;
!     endif                                                                 ! DEBUG ;
!     write(6,'(A, I8,F14.8)') '| NRG_BND =  ', jn, eaw(jn)                 ! DEBUG ;
!! Amber harmonic bond interaction for comparison with morsify ''''''''''''''''''''''

         dfw(jn) = (df+df)/rij0
      end do
      
      !         ----- CALCULATION OF THE FORCE -----
     
      do jn = 1,maxlen
         i3 = ib(jn+ist)
         j3 = jb(jn+ist)
         df = dfw(jn)
         xa = df*xij(jn)
         ya = df*yij(jn)
         za = df*zij(jn)
         f(i3+1) = f(i3+1)-xa
         f(i3+2) = f(i3+2)-ya
         f(i3+3) = f(i3+3)-za
         f(j3+1) = f(j3+1)+xa
         f(j3+2) = f(j3+2)+ya
         f(j3+3) = f(j3+3)+za
         bnd_vir(1,1)=bnd_vir(1,1)+xa*xij(jn)
         bnd_vir(1,2)=bnd_vir(1,2)+xa*yij(jn)
         bnd_vir(1,3)=bnd_vir(1,3)+xa*zij(jn)
         bnd_vir(2,1)=bnd_vir(2,1)+ya*xij(jn)
         bnd_vir(2,2)=bnd_vir(2,2)+ya*yij(jn)
         bnd_vir(2,3)=bnd_vir(2,3)+ya*zij(jn)
         bnd_vir(3,1)=bnd_vir(3,1)+za*xij(jn)
         bnd_vir(3,2)=bnd_vir(3,2)+za*yij(jn)
         bnd_vir(3,3)=bnd_vir(3,3)+za*zij(jn)

!! Amber harmonic bond interaction for comparison with morsify ''''''''''''''''''''''
!     if( jn == 1 ) then                                                    ! DEBUG ;
!        write(6,'(A)') '|'                                                 ! DEBUG ;
!        write(6,'(A)') '| ^^^ Amber harmonic force ////////////////////'   ! DEBUG ;
!        write(6,'(A)') '|'                                                 ! DEBUG ;
!     endif                                                                 ! DEBUG ;
!     write(6,'(A,2I8,F14.8)') '| FX = ', i3+1, j3+1, xa                    ! DEBUG ;
!     write(6,'(A,2I8,F14.8)') '| FY = ', i3+2, j3+2, ya                    ! DEBUG ;
!     write(6,'(A,2I8,F14.8)') '| FZ = ', i3+3, j3+3, za                    ! DEBUG ;
!     write(6,'(A)') '|'                                                    ! DEBUG ;
!! Amber harmonic bond interaction for comparison with morsify ''''''''''''''''''''''

      end do

      do ksum = 1, maxlen
         ebl = ebl + eaw(ksum)
      end do
      ist = ist+maxlen
   end if  ! (maxlen > 0)
   if(.not.skip) goto 4200
   
   !     ----- ALL DONE -----
   
   eb = ebl
   return
end subroutine bond 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                     ANGL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine angl here]
subroutine angl(nbain,it,jt,kt,ict,x,xx,ix,f,eba)

   use decomp, only : decpr, decpair, decangle
   use parms, only: teq, tk, rub,rkub, charmm
   use constants, only: one, third

#ifdef LES
   use pimd_vars, only: ipimd,ineb,nrg_all, nbead, nbead_inv
#  ifdef MPI
      use remd, only : rem
#  endif
#endif
#ifdef MPI /* SOFT CORE */
   use softcore, only: ifsc, nsc, sc_ener, oneweight
#endif
   implicit none

#ifdef MPI
#  include "parallel.h"
#endif

   !-------------passed-in variables  ---------------------------
   _REAL_ x(*),xx(*),f(*),eba
   integer ix(*),it(*),jt(*),kt(*),ict(*),nbain
   !------------- local variables  ---------------------------
   logical skip

   !     ----- ROUTINE TO GET THE ANGLE ENERGIES AND FORCES FOR THE
   !           POTENTIAL OF THE TYPE CT*(T-T0)**2

#  include "box.h"
#  include "nmr.h"
#  include "md.h"
#  include "files.h"
#  include "flocntrl.h"
#ifdef LES
#  include "les.h" 
#endif
   integer,parameter :: maxangles=190
   integer nba,maxlen,lenc,ksum,i3,j3,k3,ic,istc,jn,ist
   _REAL_ dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9
   _REAL_ sth,cik,cii,ckk,da,df,ant0,ct0,ct1,ct2,st
   _REAL_ rij0,rik0,rkj0,ebal,ecnl,eub,ecn
   _REAL_ xij(maxangles),yij(maxangles),zij(maxangles),xkj(maxangles), &
          ykj(maxangles), zkj(maxangles),cst(maxangles),eaw(maxangles), &
          rij(maxangles),rkj(maxangles),rik(maxangles), dfw(maxangles), &
          ant(maxangles),rij1(maxangles)

   _REAL_ ebw(maxangles),xa,ya,za !For charmm

   ! previously declared, see passed-in variables.  srb
   !     DIMENSION IT(*),JT(*),KT(*),ICT(*),X(*),XX(*),IX(*),F(*)
   integer piece,start,end,nbtmp,newnb
   integer ii,jj,kk
   _REAL_ pt999
   data pt999 /0.9990d0/

!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  EVB variables                                                   ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   integer :: ndx

   if ( do_angle == 0 .or. nbain == 0 )return

   nba=nbain
#ifdef MPI
   piece = nba/numtasks
   start = mytaskid*piece+1
   end   = mytaskid*piece+piece
   if(mytaskid == (numtasks-1)) end = nba
   nbtmp = nba
   !     JV Use actual count this PE will do, reset at end of routine
   nba = end
   ist = start -1
#else
   ! FLOW CONTROL FLAG (debug)
   ist = 0
#endif

   ebal = 0.0d0
   ecnl = 0.0d0
   eub  = 0.0d0
#ifdef LES
   elesa = 0.0d0
#endif

   !     ----- GRAND LOOP FOR THE angle STUFF -----

   4200 continue
   maxlen = maxangles
   skip = (ist+maxlen) > nba
   if (skip) maxlen = nba-ist
   if (maxlen <= 0) goto 220

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = jt(jn+ist)
      k3 = kt(jn+ist)

      !           ----- CALCULATION OF THE angle -----

      xij(jn) = x(i3+1)-x(j3+1)
      yij(jn) = x(i3+2)-x(j3+2)
      zij(jn) = x(i3+3)-x(j3+3)
      xkj(jn) = x(k3+1)-x(j3+1)
      ykj(jn) = x(k3+2)-x(j3+2)
      zkj(jn) = x(k3+3)-x(j3+3)
   end do

   do jn = 1,maxlen
      rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
      rkj0 = xkj(jn)*xkj(jn)+ykj(jn)*ykj(jn)+zkj(jn)*zkj(jn)
      rik0 = sqrt(rij0*rkj0)
      ct0 = (xij(jn)*xkj(jn)+yij(jn)*ykj(jn)+zij(jn)*zkj(jn))/rik0
      ct1 = max(-pt999,ct0)
      ct2 = min(pt999,ct1)
      cst(jn) = ct2
      ant(jn) = acos(ct2)
      rij(jn) = rij0
      rkj(jn) = rkj0
      rik(jn) = rik0
   end do

   !         ----- CALCULATION OF THE ENERGY AND DER -----

   do jn = 1,maxlen
      ic = ict(jn+ist)
      ant0 = ant(jn)
      da = ant0-teq(ic)
      !                                 for rms deviation from ideal angles:
      eadev = eadev + da*da
      df = tk(ic)*da
      eaw(jn) = df*da
#ifdef MPI /* SOFT CORE */
      ! For dual-topology softcore runs, angles involving sc atoms are modified here
      if (ifsc /= 0) then
         ii = (it(jn+ist) + 3)/3
         jj = (jt(jn+ist) + 3)/3
         kk = (kt(jn+ist) + 3)/3
         ! Check if a softcore atom is involved in this angle
         if (nsc(ii) == 1 .or. nsc(jj) == 1 .or. nsc(kk) == 1) then
            ! This angle needs to 
            ! a) get its energy removed from Eangle
            ! b) get its force scaled up by 1/weight
            sc_ener(2) = sc_ener(2) + eaw(jn)
            eaw(jn) = 0.0d0
            df = df * oneweight
         end if
      end if
      if(decpr .and. (idecomp == 1 .or. idecomp == 2)) then
         ii = (it(jn+ist) + 3)/3
         jj = (jt(jn+ist) + 3)/3
         kk = (kt(jn+ist) + 3)/3
         if(ifsc == 0) then
            call decangle(ii,jj,kk,eaw(jn)/(nstlim/ntpr))
         else if(nsc(ii) /= 1 .and. nsc(jj) /= 1 .and. nsc(kk) /= 1) then
            call decangle(ii,jj,kk,eaw(jn)/(nstlim/ntpr))
         end if
      end if
#else
      if(idecomp == 1 .or. idecomp == 2) then
         ii = (it(jn+ist) + 3)/3
         jj = (jt(jn+ist) + 3)/3
         kk = (kt(jn+ist) + 3)/3
         call decangle(ii,jj,kk,eaw(jn))
      end if
#endif
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ii = (it(jn+ist) + 3)/3                                                                    
         jj = (jt(jn+ist) + 3)/3                                                                    
         kk = (kt(jn+ist) + 3)/3
         if(cnum(ii) > 0 .or. cnum(jj) > 0 .or. cnum(kk) > 0) then
           elesa = elesa + eaw(jn)  
         endif
      end if
#  endif
#endif

#if defined(LES)
!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      if( ipimd>0.or.ineb>0 ) then
         ii = (it(jn+ist) + 3)/3
         jj = (jt(jn+ist) + 3)/3
         kk = (kt(jn+ist) + 3)/3
  
         if( cnum(ii) == 0 .and. cnum(jj) == 0 .and. cnum(kk) == 0 ) then
            nrg_all(1:nbead) = nrg_all(1:nbead) + eaw(jn) * nbead_inv
         else
            if( cnum(ii) /= 0 ) then
               ndx = cnum(ii)
            else if( cnum(jj) /= 0 ) then
               ndx = cnum(jj)
            else 
               ndx = cnum(kk)
            endif     
            nrg_all(ndx) = nrg_all(ndx) + eaw(jn)
         endif
      endif
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#endif /* LES */

      dfw(jn) = -(df+df)/sin(ant0)
   end do

   !         ----- CALCULATION OF THE FORCE -----

   do jn = 1,maxlen
      i3 = it(jn+ist)
      j3 = jt(jn+ist)
      k3 = kt(jn+ist)

      st = dfw(jn)
      sth = st*cst(jn)
      cik = st/rik(jn)
      cii = sth/rij(jn)
      ckk = sth/rkj(jn)
      dt1 = cik*xkj(jn)-cii*xij(jn)
      dt2 = cik*ykj(jn)-cii*yij(jn)
      dt3 = cik*zkj(jn)-cii*zij(jn)
      dt7 = cik*xij(jn)-ckk*xkj(jn)
      dt8 = cik*yij(jn)-ckk*ykj(jn)
      dt9 = cik*zij(jn)-ckk*zkj(jn)
      dt4 = -dt1-dt7
      dt5 = -dt2-dt8
      dt6 = -dt3-dt9

      f(i3+1) = f(i3+1)-dt1
      f(i3+2) = f(i3+2)-dt2
      f(i3+3) = f(i3+3)-dt3
      f(j3+1) = f(j3+1)-dt4
      f(j3+2) = f(j3+2)-dt5
      f(j3+3) = f(j3+3)-dt6
      f(k3+1) = f(k3+1)-dt7
      f(k3+2) = f(k3+2)-dt8
      f(k3+3) = f(k3+3)-dt9
   end do

   if (charmm) then

   !         --- include Urey-Bradley terms:

     do jn = 1,maxlen
        i3 = it(jn+ist)
        j3 = kt(jn+ist)
        xij(jn) = x(i3+1)-x(j3+1)
        yij(jn) = x(i3+2)-x(j3+2)
        zij(jn) = x(i3+3)-x(j3+3)
     end do

     do jn = 1,maxlen
        rij0 = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
        rij1(jn) = one/sqrt(rij0)
     end do

     !         ----- CALCULATION OF THE ENERGY AND DER -----

     do jn = 1,maxlen
        ic = ict(jn+ist)
        rij0 = one/rij1(jn)
        da = rij0-rub(ic)
        df = rkub(ic)*da
        ebw(jn) = df*da
#       ifdef MPI
        if(decpr .and. (idecomp == 1 .or. idecomp == 2)) then
           ii = (it(jn+ist) + 3) *third
           kk = (kt(jn+ist) + 3) *third
           call decpair(4,ii,kk,ebw(jn)/(nstlim/ntpr))
        end if
#       else
        if(idecomp == 1 .or. idecomp == 2) then
           ii = (it(jn+ist) + 3) *third
           kk = (kt(jn+ist) + 3) *third
           call decpair(4,ii,kk,ebw(jn))
        end if
#       endif
        dfw(jn) = (df+df)*rij1(jn)
     end do

     !         ----- CALCULATION OF THE FORCE -----

     do jn = 1,maxlen
        i3 = it(jn+ist)
        j3 = kt(jn+ist)
        df = dfw(jn)
        xa = df*xij(jn)
        ya = df*yij(jn)
        za = df*zij(jn)
        f(i3+1) = f(i3+1)-xa
        f(i3+2) = f(i3+2)-ya
        f(i3+3) = f(i3+3)-za
        f(j3+1) = f(j3+1)+xa
        f(j3+2) = f(j3+2)+ya
        f(j3+3) = f(j3+3)+za
     end do
     do ksum = 1, maxlen
        eub = eub + ebw(ksum)
     end do
     !-------------- end urey bradley terms ----------------------------
   end if !charmm

   do ksum=1, maxlen
      ebal = ebal + eaw(ksum)
   end do
   ist = ist+maxlen

   220 if(.not.skip) goto 4200
   eba = ebal
   if (charmm) then
        write(6,*) 'Urey-Bradley energy: ',eub
     eba = eba + eub
   end if

#ifdef MPI
   nba = nbtmp
#endif

   return
end subroutine angl
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   EPHI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ephi here]
subroutine ephi(nphiin,ip,jp,kp,lp,icp,cg,iac,x,xx,ix,f,dvdl, &
      ep,enbp,eelp,dcharge)

   use decomp, only : decpr, decpair, decphi
   use parms, only: ipn,pn,pk,gamc,gams,fmn,cn1,cn2,cn114,cn214, charmm
   use constants, only : zero, one, two, six, twelve, PI

!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  EVB modules                                                     ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#ifdef LES
   use pimd_vars, only: ipimd,ineb,nrg_all, nbead, nbead_inv
#  ifdef MPI
      use remd, only : rem
#  endif
#endif
#ifdef MPI /* SOFT CORE */
   use softcore, only: ifsc, nsc, sc_ener, oneweight
#endif
   implicit none

#ifdef MPI
#  include "parallel.h"
#else
   integer numtasks,mytaskid
   parameter (numtasks=1,mytaskid=0)
#endif
#if defined(LES)
#  include "les.h"
#  include "memory.h"
#else

   _REAL_ lfac
#endif

#  include "ew_frc.h"
#  include "box.h"
! Antonios added
#  include "HB.1.h"
! Antonios end
#  include "md.h"
#  include "files.h"
#  include "flocntrl.h"
#  include "ew_cntrl.h"
#  include "ew_erfc_spline.h"
#  include "extra_pts.h"

   !-------------passed-in variables  ---------------------------
   integer ix(*),ip(*),jp(*),kp(*),lp(*),icp(*),iac(*)
   _REAL_   ep,enbp,eelp,ecn
   _REAL_ cg(*),dcharge(*),dvdl
   _REAL_ x(*),xx(*),f(*),eba
   integer nphiin,mba,mphi

   !------------- local variables  ---------------------------
   logical skip

   _REAL_ intdieli
   integer max190,maxlen
   parameter (max190=190)
   _REAL_ xij(max190),yij(max190),zij(max190),xkj(max190), &
         ykj(max190),zkj(max190),xkl(max190),ykl(max190), &
         zkl(max190),dx(max190),dy(max190),dz(max190),gx(max190), &
         gy(max190),gz(max190),ct(max190),cphi(max190),sphi(max190), &
         z1(max190),z2(max190),fxi(max190),fyi(max190),fzi(max190), &
         fxj(max190),fyj(max190),fzj(max190),fxk(max190),fyk(max190), &
         fzk(max190),fxl(max190),fyl(max190),fzl(max190),df(max190)

   _REAL_ xa,ya,za,f1,f2,r1,r2,r6,r12,dfn,g,fmuln
   _REAL_ dr1,dr2,dr3,dr4,dr5,dr6,drx,dry,drz
   _REAL_ dc1,dc2,dc3,dc4,dc5,dc6
   _REAL_ df0,df1,dflim,dums,cosnp,sinnp
   _REAL_ z10,z20,z12,z11,z22,ap0,ap1,ct0,ct1,s,ftem
   _REAL_ epl,ecnl,enbpl,eelpl,scnb0,scee0
   integer ksum,lenc,ia1,ia2,ibig,isml,ii,jj,kk,ll
   integer ic,inc,ic0,iduml,idumi,kdiv
   integer i3,j3,k3,l3,k3t,l3t
   integer jn,istc,ist,nphi


   _REAL_ epw(max190)
   _REAL_ gmul(10)


   data gmul/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00, &
         0.0d+00,8.0d+00,0.0d+00,10.0d+00/
   _REAL_ tm24,tm06,tenm3
   data tm24,tm06,tenm3/1.0d-18,1.0d-06,1.0d-03/

   !     ---- ARRAYS GAMC = PK*COS(PHASE) AND GAMS = PK*SIN(PHASE) ----

   integer piece,start,end,nbtmp,newnb

!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  EVB variables                                                   ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   integer :: ndx

   if ( do_ephi == 0 .or. nphiin == 0 )return
   nbtmp=nphiin
   if(numtasks > 1)then
      piece = nphiin/numtasks
      start = mytaskid*piece+1
      end= mytaskid*piece+piece
      if(mytaskid == (numtasks-1)) end = nphiin
   else
      start=1
      end=nphiin
   end if
   nphi = end
   ist = start -1

   intdieli = one/intdiel
   epl = zero
   ecnl = zero
   enbpl = zero
   eelpl = zero
#ifdef LES
   elesd = zero
#endif
   scnb0 = one/scnb
   scee0 = one/scee

   !     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----

   4200 continue
   maxlen = 190
   skip = (ist+maxlen) > nphi
   if(skip) maxlen = nphi-ist
   if(maxlen <= 0) goto 820

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      j3 = jp(jn+ist)
      k3t = kp(jn+ist)
      l3t = lp(jn+ist)
      k3 = iabs(k3t)
      l3 = iabs(l3t)

      !           ----- CALCULATION OF ij, kj, kl VECTORS -----

      xij(jn) = x(i3+1)-x(j3+1)
      yij(jn) = x(i3+2)-x(j3+2)
      zij(jn) = x(i3+3)-x(j3+3)
      xkj(jn) = x(k3+1)-x(j3+1)
      ykj(jn) = x(k3+2)-x(j3+2)
      zkj(jn) = x(k3+3)-x(j3+3)
      xkl(jn) = x(k3+1)-x(l3+1)
      ykl(jn) = x(k3+2)-x(l3+2)
      zkl(jn) = x(k3+3)-x(l3+3)
   end do

   !         ----- GET THE NORMAL VECTOR -----

   do jn = 1,maxlen
      dx(jn) = yij(jn)*zkj(jn)-zij(jn)*ykj(jn)
      dy(jn) = zij(jn)*xkj(jn)-xij(jn)*zkj(jn)
      dz(jn) = xij(jn)*ykj(jn)-yij(jn)*xkj(jn)
      gx(jn) = zkj(jn)*ykl(jn)-ykj(jn)*zkl(jn)
      gy(jn) = xkj(jn)*zkl(jn)-zkj(jn)*xkl(jn)
      gz(jn) = ykj(jn)*xkl(jn)-xkj(jn)*ykl(jn)
   end do

   do jn = 1,maxlen
      fxi(jn) = sqrt(dx(jn)*dx(jn) &
            +dy(jn)*dy(jn) &
            +dz(jn)*dz(jn)+tm24)
      fyi(jn) = sqrt(gx(jn)*gx(jn) &
            +gy(jn)*gy(jn) &
            +gz(jn)*gz(jn)+tm24)
      ct(jn) = dx(jn)*gx(jn)+dy(jn)*gy(jn)+dz(jn)*gz(jn)
   end do

   !         ----- BRANCH IF LINEAR DIHEDRAL -----

   do jn = 1,maxlen
#ifdef CRAY_PVP
      bit = one/fxi(jn)
      bik = one/fyi(jn)
      z10 = cvmgt(zero,bit,tenm3 > fxi(jn))
      z20 = cvmgt(zero,bik,tenm3 > fyi(jn))
#else
      z10 = one/fxi(jn)
      z20 = one/fyi(jn)
      if (tenm3 > fxi(jn)) z10 = zero
      if (tenm3 > fyi(jn)) z20 = zero
#endif
      z12 = z10*z20
      z1(jn) = z10
      z2(jn) = z20
#ifdef CRAY_PVP
      ftem = cvmgz(zero,one,z12)
#else
      ftem = zero
      if (z12 /= zero) ftem = one
#endif
      fzi(jn) = ftem
      ct0 = min(one,ct(jn)*z12)
      ct1 = max(-one,ct0)
      s = xkj(jn)*(dz(jn)*gy(jn)-dy(jn)*gz(jn))+ &
            ykj(jn)*(dx(jn)*gz(jn)-dz(jn)*gx(jn))+ &
            zkj(jn)*(dy(jn)*gx(jn)-dx(jn)*gy(jn))
      ap0 = acos(ct1)
      ap1 = pi-sign(ap0,s)
      ct(jn) = ap1
      cphi(jn) = cos(ap1)
      sphi(jn) = sin(ap1)
   end do

   !         ----- CALCULATE THE ENERGY AND THE DERIVATIVES WITH RESPECT TO
   !               COSPHI -----

   do jn = 1,maxlen
      ic = icp(jn+ist)
      inc = ipn(ic)
      ct0 = pn(ic)*ct(jn)
      cosnp = cos(ct0)
      sinnp = sin(ct0)
      epw(jn) = (pk(ic)+cosnp*gamc(ic)+sinnp*gams(ic))*fzi(jn)
#ifdef MPI
      if(decpr .and. (idecomp == 1 .or. idecomp == 2)) then
         ii = (ip(jn+ist) + 3)/3
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3
         if(ifsc == 0) then
            call decphi(ii,jj,kk,ll,epw(jn)/(nstlim/ntpr))
         else if( nsc(ii) /= 1 .and. nsc(jj) /= 1 .and. nsc(kk) /= 1 &
                               .and. nsc(ll) /= 1) then
            call decphi(ii,jj,kk,ll,epw(jn)/(nstlim/ntpr))
         end if
      end if
#else
      if(idecomp == 1 .or. idecomp == 2) then
         ii = (ip(jn+ist) + 3)/3
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3
         call decphi(ii,jj,kk,ll,epw(jn))
      end if
#endif
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ii = (ip(jn+ist) + 3)/3 
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3 
         if(cnum(ii) > 0 .or. cnum(jj) > 0 .or. cnum(kk) > 0 .or. &
         cnum(ll) > 0) then
           elesd = elesd + epw(jn)
         endif
      end if
#  endif
#endif

#if defined(LES)
!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      if(ipimd>0.or.ineb>0) then
         ii = (ip(jn+ist) + 3)/3
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3
 
         if( cnum(ii) == 0 .and. cnum(jj) == 0 .and. &
             cnum(kk) == 0 .and. cnum(ll) == 0       ) then
            nrg_all(1:nbead) = nrg_all(1:nbead) + epw(jn) * nbead_inv
         else
            if( cnum(ii) /= 0 ) then
               ndx = cnum(ii)
            else if( cnum(jj) /= 0 ) then
               ndx = cnum(jj)
            else if( cnum(kk) /= 0 ) then
               ndx = cnum(kk)
            else
               ndx = cnum(ll)
            endif
            nrg_all(ndx) = nrg_all(ndx) + epw(jn)
         endif
      endif
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#endif /* LES */

      df0 = pn(ic)*(gamc(ic)*sinnp-gams(ic)*cosnp)
      dums = sphi(jn)+sign(tm24,sphi(jn))
      dflim = gamc(ic)*(pn(ic)-gmul(inc)+gmul(inc)*cphi(jn))
#ifdef CRAY_PVP
      df1 = cvmgt(dflim,df0/dums,tm06 > abs(dums))
#else
      df1 = df0/dums
      if(tm06 > abs(dums)) df1 = dflim
#endif
#ifdef MPI /* SOFT CORE */
      ! For dual-topology softcore runs, dihedrals involving sc atoms are modified here
      if (ifsc /= 0) then
         ii = (ip(jn+ist) + 3)/3
         jj = (jp(jn+ist) + 3)/3
         kk = (iabs(kp(jn+ist)) + 3)/3
         ll = (iabs(lp(jn+ist)) + 3)/3
         ! Check if a softcore atom is involved in this dihedral
         if (nsc(ii) == 1 .or. nsc(jj) == 1 .or. nsc(kk) == 1 .or. nsc(ll) == 1) then
            ! This dihedral needs to 
            ! a) get its energy removed from Edihed
            ! b) get its force scaled up by 1/weight
            sc_ener(3) = sc_ener(3) + epw(jn)
            epw(jn) = 0.0d0
            df1 = df1 * oneweight
         end if
      end if
#endif
      df(jn) = df1*fzi(jn)
   end do

   !         ----- NOW DO TORSIONAL FIRST DERIVATIVES -----

   do jn = 1,maxlen

      !           ----- NOW, SET UP ARRAY DC = FIRST DER. OF COSPHI W/RESPECT
      !                 TO THE CARTESIAN DIFFERENCES T -----

      z11 = z1(jn)*z1(jn)
      z12 = z1(jn)*z2(jn)
      z22 = z2(jn)*z2(jn)
      dc1 = -gx(jn)*z12-cphi(jn)*dx(jn)*z11
      dc2 = -gy(jn)*z12-cphi(jn)*dy(jn)*z11
      dc3 = -gz(jn)*z12-cphi(jn)*dz(jn)*z11
      dc4 =  dx(jn)*z12+cphi(jn)*gx(jn)*z22
      dc5 =  dy(jn)*z12+cphi(jn)*gy(jn)*z22
      dc6 =  dz(jn)*z12+cphi(jn)*gz(jn)*z22

      !           ----- UPDATE THE FIRST DERIVATIVE ARRAY -----

      dr1 = df(jn)*(dc3*ykj(jn) - dc2*zkj(jn))
      dr2 = df(jn)*(dc1*zkj(jn) - dc3*xkj(jn))
      dr3 = df(jn)*(dc2*xkj(jn) - dc1*ykj(jn))
      dr4 = df(jn)*(dc6*ykj(jn) - dc5*zkj(jn))
      dr5 = df(jn)*(dc4*zkj(jn) - dc6*xkj(jn))
      dr6 = df(jn)*(dc5*xkj(jn) - dc4*ykj(jn))
      drx = df(jn)*(-dc2*zij(jn) + dc3*yij(jn) + &
            dc5*zkl(jn) - dc6*ykl(jn))
      dry = df(jn)*( dc1*zij(jn) - dc3*xij(jn) - &
            dc4*zkl(jn) + dc6*xkl(jn))
      drz = df(jn)*(-dc1*yij(jn) + dc2*xij(jn) + &
            dc4*ykl(jn) - dc5*xkl(jn))
      fxi(jn) = - dr1
      fyi(jn) = - dr2
      fzi(jn) = - dr3
      fxj(jn) = - drx + dr1
      fyj(jn) = - dry + dr2
      fzj(jn) = - drz + dr3
      fxk(jn) = + drx + dr4
      fyk(jn) = + dry + dr5
      fzk(jn) = + drz + dr6
      fxl(jn) = - dr4
      fyl(jn) = - dr5
      fzl(jn) = - dr6
   end do  !  jn = 1,maxlen

   !         ----- END OF A STRIP OF DIHEDRALS -----

   !--------------- 1-4 NONBOND  INTERACTIONS -----------------------------
   !         for ewald calcs, 1-4 interactions are done in get_14_cg()
   !                                                    or get_14_dipole().

#if !defined(LES)
   igbNEzero: if( igb /= 0 ) then
#endif

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      l3t = lp(jn+ist)
      l3 = iabs(l3t)
      xij(jn) = x(i3+1)-x(l3+1)
      yij(jn) = x(i3+2)-x(l3+2)
      zij(jn) = x(i3+3)-x(l3+3)
   end do

   !             ----- NOW LOOP OVER ALL THE DIHEDRALS (CONSTANT DIEL) -----

   do jn = 1,maxlen
      ct(jn) = xij(jn)*xij(jn)+yij(jn)*yij(jn)+zij(jn)*zij(jn)
   end do

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      k3t = kp(jn+ist)
      l3t = lp(jn+ist)
      ic0 = icp(jn+ist)
      idumi = sign(1,k3t)
      iduml = sign(1,l3t)
      kdiv = (2+idumi+iduml)/4
      l3 = iabs(l3t)
      fmuln = kdiv*fmn(ic0)

      ii = (i3+3)/3
      jj = (l3+3)/3
      ia1 = iac(ii)
      ia2 = iac(jj)
      ibig = max0(ia1,ia2)
      isml = min0(ia1,ia2)
      ic = ibig*(ibig-1)/2+isml

      !             ----- CALCULATE THE 14-EEL ENERGY -----

      r2 = fmuln/ct(jn)
      r1 = sqrt(r2)

#ifdef LES
      lfac=lesfac(nlesty*(lestyp(ii)-1)+lestyp(jj))
#else 
      lfac = 1.d0
#endif

      if( eedmeth == 5 ) then
         g = cg(ii)*cg(jj)*r2*lfac
      else
         g = cg(ii)*cg(jj)*r1*lfac*intdieli
      end if

      if (icnstph /= 0 .and. mod(irespa,ntcnstph) == 0) then
         dvdl = dvdl + scee0*(dcharge(ii)* &
               dcharge(jj)*r1*lfac*intdieli-g)
      end if
      sphi(jn) = g
      if(idecomp > 0) then
         ii = (i3+3)/3
         jj = (l3+3)/3
#        ifdef MPI
         if(decpr .and. idecomp == 1) then
            call decpair(4,ii,jj,sphi(jn)*scee0/(nstlim/ntpr))
         else if(decpr .and. idecomp == 2) then
            call decpair(2,ii,jj,sphi(jn)*scee0/(nstlim/ntpr))
         endif
#        else
         if(idecomp == 1) then
            call decpair(4,ii,jj,sphi(jn)*scee0)
         else if(idecomp == 2) then
            call decpair(2,ii,jj,sphi(jn)*scee0)
            !             else if(idecomp.eq.3) then
            !               --- not considered since
            !                     no pairwise decomp for internal energies
         else if(idecomp == 4) then
            call decpair(-2,ii,jj,sphi(jn)*scee0)
         end if
#        endif
      end if
      r6 = r2*r2*r2
      r12 = r6*r6
      if (charmm) then
        f1 = cn114(ic)*r12*lfac
        f2 = cn214(ic)*r6*lfac
      else
        f1 = cn1(ic)*r12*lfac
        f2 = cn2(ic)*r6*lfac
      end if
      cphi(jn) = f1-f2
!Antonios added
                  if(maphb(II,JJ).lt.1) then
                  EHBV=EHBV+F1-F2
                  else
                  EHBA=EHBA+F1-F2
                  endif
!Antonios end

      if(idecomp > 0) then
         ii = (i3+3)/3
         jj = (l3+3)/3
#        ifdef MPI
         if(decpr .and. idecomp == 1) then
            call decpair(4,ii,jj,cphi(jn)*scnb0/(nstlim/ntpr))
         else if(decpr .and. idecomp == 2) then
            call decpair(3,ii,jj,cphi(jn)*scnb0/(nstlim/ntpr))
         end if
#        else
         if(idecomp == 1) then
            call decpair(4,ii,jj,cphi(jn)*scnb0)
         else if(idecomp == 2) then
            call decpair(3,ii,jj,cphi(jn)*scnb0)
            !             else if(idecomp.eq.3) then
            !               --- not considered since
            !                     no pairwise decomp for internal energies
         else if(idecomp == 4) then
            call decpair(-3,ii,jj,cphi(jn)*scnb0)
         end if
#        endif
      end if

#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ii = (i3+3)/3
         jj = (l3+3)/3
       if(cnum(ii) > 0 .or. cnum(jj) > 0) then
         elesd = elesd + cphi(jn)*scnb0 + sphi(jn)*scee0
       endif
      endif
#  endif
#endif

#if defined(LES)
      !! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      if(ipimd>0.or.ineb>0) then
         if( cnum(ii) == 0 .and. cnum(jj) == 0 ) then
            nrg_all(1:nbead) = nrg_all(1:nbead) &
               + ( cphi(jn) * scnb0 + sphi(jn) * scee0 ) * nbead_inv
         else
            if( cnum(ii) /= 0 ) then
               ndx = cnum(ii)
            else
               ndx = cnum(jj)
            endif
            nrg_all(ndx) = nrg_all(ndx) + cphi(jn) * scnb0 & 
                                        + sphi(jn) * scee0 
         endif
      end if
      !! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#endif /* LES */

      if( eedmeth == 5 ) then
         dfn =((-twelve*f1+six*f2)*scnb0-two*g*scee0)*r2
      else
         dfn =((-twelve*f1+six*f2)*scnb0-g*scee0)*r2
      end if

      xa = xij(jn)*dfn
      ya = yij(jn)*dfn
      za = zij(jn)*dfn
      fxi(jn) = fxi(jn)-xa
      fyi(jn) = fyi(jn)-ya
      fzi(jn) = fzi(jn)-za
      fxl(jn) = fxl(jn)+xa
      fyl(jn) = fyl(jn)+ya
      fzl(jn) = fzl(jn)+za
      e14vir(1,1)=e14vir(1,1)+xa*xij(jn)
      e14vir(2,1)=e14vir(2,1)+ya*xij(jn)
      e14vir(3,1)=e14vir(3,1)+za*xij(jn)
      e14vir(1,2)=e14vir(1,2)+xa*yij(jn)
      e14vir(2,2)=e14vir(2,2)+ya*yij(jn)
      e14vir(3,2)=e14vir(3,2)+za*yij(jn)
      e14vir(1,3)=e14vir(1,3)+xa*zij(jn)
      e14vir(2,3)=e14vir(2,3)+ya*zij(jn)
      e14vir(3,3)=e14vir(3,3)+za*zij(jn)
   end do  !  jn = 1,maxlen
   !
   !------------------- End of 1-4 interactions -----------------------------
   !
#if !defined(LES)
   endif igbNEzero
#endif
   !         ----- SUMUP ALL THE GRADIENTS -----

   do jn = 1,maxlen
      i3 = ip(jn+ist)
      j3 = jp(jn+ist)
      k3 = iabs(kp(jn+ist))
      l3 = iabs(lp(jn+ist))

      f(i3+1) = f(i3+1) + fxi(jn)
      f(i3+2) = f(i3+2) + fyi(jn)
      f(i3+3) = f(i3+3) + fzi(jn)
      f(j3+1) = f(j3+1) + fxj(jn)
      f(j3+2) = f(j3+2) + fyj(jn)
      f(j3+3) = f(j3+3) + fzj(jn)
      f(k3+1) = f(k3+1) + fxk(jn)
      f(k3+2) = f(k3+2) + fyk(jn)
      f(k3+3) = f(k3+3) + fzk(jn)
      f(l3+1) = f(l3+1) + fxl(jn)
      f(l3+2) = f(l3+2) + fyl(jn)
      f(l3+3) = f(l3+3) + fzl(jn)
   end do

   do ksum = 1,maxlen
      enbpl = enbpl+cphi(ksum)
      eelpl = eelpl+sphi(ksum)
      epl   = epl  +epw(ksum)
   end do
   ist = ist+maxlen
   820 if(.not.skip) goto 4200

   !     ---- ALL DONE -----

   enbp = enbpl*scnb0
   eelp = eelpl*scee0
   ep   = epl
   ecn = ecnl
#ifdef MPI
   nphiin = nbtmp
#endif
   return
end subroutine ephi
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine capwat here]
subroutine capwat(nat,x,f,ecap)

   use poisson_boltzmann, only : outflagorig

   implicit _REAL_ (a-h,o-z)
#ifdef MPI
#  include "parallel.h"
#endif

   !     ----- ROUTINE TO CALCULATE THE CAP FORCE -----

#  include "box.h"
#  include "flocntrl.h"
   dimension x(3,*),f(3,*)
   data tm34,zero/1.0d-34,0.0d0/

   ecap = 0.d0

   ! FLOW CONTROL FLAG (debug)
   if ( do_cap == 0 )return
#ifdef MPI
   do i = natcap+1+mytaskid,nat,numtasks
#else
   do i = natcap+1,nat
#endif

      xa = xcap-x(1,i)
      ya = ycap-x(2,i)
      za = zcap-x(3,i)
      da = sqrt(xa*xa+ya*ya+za*za+tm34)
      delta = max(zero,da-cutcap)
      if( ifcap == 2 .or. ifcap == 5) then
         if(outflagorig(i) == 1) then
            delta = 0.0d0
         end if
      end if
      ecap = ecap + 0.5*fcap*delta*delta
      df = fcap*delta/da
      f(1,i) = f(1,i)+df*xa
      f(2,i) = f(2,i)+df*ya
      f(3,i) = f(3,i)+df*za
   end do
   return
end subroutine capwat
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ 
subroutine xconst(natc,econ,igroup,x,f,xc,weit & 
     ,natom)
   implicit _REAL_ (a-h,o-z)
#ifdef MPI
#  include "parallel.h"
#endif
#  include "flocntrl.h"
   integer,intent(in) :: natom
   logical myatm(natom)

   !     ----- ROUTINE TO PUT HARMONIC CONSTRAINTS FOR POSITION -----

   dimension igroup(*),x(*),f(*),xc(*),weit(*)

   ! FLOW CONTROL FLAG (debug)
   if ( doxconst == 0 )return
   econ = 0.0d+00

# ifdef MPI
      do ii = 1+mytaskid,natc,numtasks
# else
      do ii = 1,natc
# endif
         i = igroup(ii)
         wt = weit(ii)
         i3 = 3*i-3
         ax = x(i3+1)-xc(i3+1)
         ay = x(i3+2)-xc(i3+2)
         az = x(i3+3)-xc(i3+3)
         wx = wt*ax
         wy = wt*ay
         wz = wt*az
         eadd = wx*ax+wy*ay+wz*az
         econ = econ+eadd
         f(i3+1) = f(i3+1)-(wx+wx)
         f(i3+2) = f(i3+2)-(wy+wy)
         f(i3+3) = f(i3+3)-(wz+wz)
      end do
   return
end subroutine xconst
!-----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ targeted MD simulations with restraints based on RMSD 
subroutine xtgtmd(econ,rmsgroup,x,f,xc,mass)
! uses reference coords just like positional restraints
!  energy = 0.5 * tgtmdfrc * nattgtrms * (rmsd(t)-tgtrmsd)**2
!  d(energy)/dxi = k * (1 - tgtrmsd/rmsdvalue) * mass(i)/totmass * (xi-xiref)
! Note: this routine modeled on xconst() and is very similar

   implicit none

   ! nattgtrms from "memory.h"
   ! ntr from "md.h"
#  include "memory.h"     
#  include "md.h"
#  include "tgtmd.h"
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#endif
#  include "flocntrl.h"

   ! variables in "tgtmd.h":
   ! integer itgtmd; logical dotgtmd,rmsok; 
   ! real tgtrmsd,tgtmdfrc,rmsdvalue

   ! routine's formal variables:
   integer rmsgroup(*)
   _REAL_  econ,x(*),f(*),xc(*),mass(*)

   ! local variables:
   _REAL_ factor,totmass,massfac, ax,ay,az
   integer ii, i, i3
   
   ! flow control flag (debug)
   if ( do_tgt == 0 ) return

   ! ENERGY DOES NOT DEPEND ON ATOMIC TERMS
   ! since 1 term, do not do in parallel
   ! others will get it from energy reduction in MPI code

   if (ntr == 0) then  ! If ntr=1, ADD to restraint energy from xconst()
      econ = 0.d0      ! otherwise initialize econ to zero.
   end if

   if (master) then
      econ = econ + 0.5d0 * nattgtrms * tgtmdfrc * (rmsdvalue-tgtrmsd)**2
   end if

   ! We need to test for small rmsdvalue since we divide by it
   ! and user is likely to use the initial structure as refc
   ! therefore giving rmsd=0. If so, skip it on this step and wait
   ! for the difference vectors to become large enough to define forces.
   ! *** This means that energy will not be conserved for such steps. ***

   if (rmsdvalue < 0.001d0) return

   factor = nattgtrms * tgtmdfrc * (1.d0 - tgtrmsd/rmsdvalue)

   totmass = 0.d0
   do ii=1,nattgtrms
      i = rmsgroup(ii)
      totmass = totmass + mass(i)
   end do

#ifdef MPI
   do ii = 1+mytaskid,nattgtrms,numtasks
#else
   do ii = 1,nattgtrms
#endif

      i = rmsgroup(ii)   ! i is atom number
      i3 = 3*i-3         ! i3+1 is start of coordinate location

      ax = x(i3+1) - xc(i3+1)
      ay = x(i3+2) - xc(i3+2)
      az = x(i3+3) - xc(i3+3)

      massfac = mass(i)/totmass

      ! ax, ay and az are distances to reference coordinate. Force is 
      ! negative derivative, positive derivs are calculated above.
      f(i3+1) = f(i3+1) - factor*ax*massfac
      f(i3+2) = f(i3+2) - factor*ay*massfac
      f(i3+3) = f(i3+3) - factor*az*massfac

   end do

   return
end subroutine xtgtmd
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ overlap and calculate RMSD of current coords to reference coords
!
subroutine rmsfit(xc,x,mass,fitgroup,rmsgroup,rmsdvalue,rmsok)
!
! Carlos added code for targeted MD (modified by VH).
! Fit XC='refc' to X='inpcrd': do not change X='inpcrd' but 
! translate/rotate XC='refc' each step. The fit is done on 'nattgtfit'
! atoms but rmsd value is calculated for 'nattgtrms' atoms
   use pimd_vars, only: ineb
   use neb_vars,  only: neb_rotat
   implicit none
   ! natom and nattgtfit,nattgtrms are from "memory.h"
   ! ntr is from "md.h"
#  include "memory.h"     
#  include "md.h"
#  include "flocntrl.h"
#ifdef MPI
#  include "parallel.h"
#endif

   ! routine's formal variables:
   _REAL_  xc(*), x(*), mass(*), rotat(3,3), rmsdvalue, rotat_in(3,3)
   logical rmsok
   integer fitgroup(*), rmsgroup(*)

   ! local variables:
   _REAL_ kabs(3,3), e(3), b(3,3)
   _REAL_ xcm2, ycm2, zcm2, det, norm
   _REAL_ xcm1, ycm1, zcm1, tmpx, tmpy, tmpz
   _REAL_ small, totmass
   _REAL_ tmprms, tmp
   integer i, ii, i3, j, k, ismall

   ! following vars needed for lapack diagonalization
   integer lwork, ldum, info
   parameter (lwork=48)
   _REAL_  kabs2(3,3), work(lwork), dum
   _REAL_  wr(3), wi(3), a(3,3)


   rmsok = .true.
   small = 1.d20

   ! zero some variables
   do i=1,3
      e(i) = 0.d0
      do j=1,3
         a(i,j)     = 0.d0
         b(i,j)     = 0.d0
         kabs(i,j)  = 0.d0
         kabs2(i,j) = 0.d0
         rotat(i,j) = 0.d0
      end do
   end do

   ! do the fitting (=overlap) for ntr=0 (no restraints)
   if (ntr == 0) then 
      ! STEPS:
      ! 1) center both coord sets on respective CM of fit regions
      ! 2) determine overlap (fit) matrix
      ! 3) rotate reference coords and get rmsdvalue
      ! 4) restore CM of both coords

      ! calculate center of mass of fit region, and center both
      ! molecules (whole molecules not just fit regions)
      xcm1 = 0.d0; ycm1 = 0.d0; zcm1 = 0.d0
      xcm2 = 0.d0; ycm2 = 0.d0; zcm2 = 0.d0
      
      totmass = 0.d0   ! totmass for fit region

      do ii=1,nattgtfit
         i = fitgroup(ii)   ! 'i' is actual atom number
         i3 = 3*i-3
         totmass = totmass + mass(i)
         xcm1 = mass(i) * x(i3+1)  + xcm1
         ycm1 = mass(i) * x(i3+2)  + ycm1
         zcm1 = mass(i) * x(i3+3)  + zcm1
         xcm2 = mass(i) * xc(i3+1)  + xcm2
         ycm2 = mass(i) * xc(i3+2)  + ycm2
         zcm2 = mass(i) * xc(i3+3)  + zcm2
      end do
      xcm1 = xcm1 / totmass
      ycm1 = ycm1 / totmass
      zcm1 = zcm1 / totmass
      xcm2 = xcm2 / totmass
      ycm2 = ycm2 / totmass
      zcm2 = zcm2 / totmass

      ! Move both molecules (all atoms) with respect to CM of fit regions.
      ! Don't modify xc (refc) coordinates, instead do all transformation
      ! that modify xc in a temporary array xctmp().
      do i=1,natom
         i3 = 3*i-3
         x(i3+1) = x(i3+1) - xcm1
         x(i3+2) = x(i3+2) - ycm1
         x(i3+3) = x(i3+3) - zcm1
         xc(i3+1) = xc(i3+1) - xcm2
         xc(i3+2) = xc(i3+2) - ycm2
         xc(i3+3) = xc(i3+3) - zcm2
      end do
      ! calculate the Kabsch matrix
      do ii=1,nattgtfit
         i = fitgroup(ii)
         i3 = 3*i-3
         kabs(1,1) = kabs(1,1) + mass(i) * x(i3+1)*xc(i3+1)
         kabs(1,2) = kabs(1,2) + mass(i) * x(i3+1)*xc(i3+2)
         kabs(1,3) = kabs(1,3) + mass(i) * x(i3+1)*xc(i3+3)
         kabs(2,1) = kabs(2,1) + mass(i) * x(i3+2)*xc(i3+1)
         kabs(2,2) = kabs(2,2) + mass(i) * x(i3+2)*xc(i3+2)
         kabs(2,3) = kabs(2,3) + mass(i) * x(i3+2)*xc(i3+3)
         kabs(3,1) = kabs(3,1) + mass(i) * x(i3+3)*xc(i3+1)
         kabs(3,2) = kabs(3,2) + mass(i) * x(i3+3)*xc(i3+2)
         kabs(3,3) = kabs(3,3) + mass(i) * x(i3+3)*xc(i3+3)
      end do
   
      ! check that the determinant is not zero
      det = kabs(1,1)*kabs(2,2)*kabs(3,3) - kabs(1,1)*kabs(2,3)*kabs(3,2) -  &
            kabs(1,2)*kabs(2,1)*kabs(3,3) + kabs(1,2)*kabs(2,3)*kabs(3,1) +  &
            kabs(1,3)*kabs(2,1)*kabs(3,2) - kabs(1,3)*kabs(2,2)*kabs(3,1)
   
      if (abs(det) < 1.d-5) then
         rmsok = .false.
         write (6,*) "small determinant in rmsfit():", abs(det)
         goto 99
      end if
   
      ! construct a positive definite matrix by multiplying kabs
      ! matrix by its transpose
      do i=1,3
         do j=1,3
            do k=1,3
               kabs2(i,j) = kabs2(i,j) + kabs(k,i)*kabs(k,j)
            end do
         end do
      end do
   
      ! Diagonalize kabs2 (matrix RtR in the above reference). On output
      ! kabs2 is destroyed and the eigenvectors are provided in 'a': 
      ! a(i=1,3 ; j) is the j-th eigenvector, 'wr' has the eigenvalues.
   
      ldum = 1
      call D_OR_S()geev('N','V',3,kabs2,3,wr,wi,dum,ldum,a,3,work,lwork,info)
   
      if (info /= 0) then
         write(6,'(" Error in diagonalization routine dgeev")')
         rmsok = .false.
         goto 99
      end if

      ! find the smallest eigenvalue
      do i=1,3
         if (wr(i) < small) then
            ismall = i
            small  = wr(i)
         end if
      end do
   
      ! generate the b vectors
      do j=1,3
         ! 'wr' very small and negative sometimes -> abs()
         norm = 1.d0/sqrt(abs(wr(j)))
         do i=1,3
            do k=1,3
               b(i,j) = b(i,j) + kabs(i,k)*a(k,j)*norm
            end do
         end do
      end do
   
      ! calculate the rotation matrix rotat
   
200   continue
   
      rotat(:, :) = 0.d0  ! array assignment
   
      do i=1,3
         do j=1,3
            do k=1,3
               rotat(i,j) = rotat(i,j) + b(i,k)*a(j,k)
            end do
         end do
      end do

      ! check if the determinant is negative
      det = rotat(1,1)*rotat(2,2)*rotat(3,3) -  &
            rotat(1,1)*rotat(2,3)*rotat(3,2) -  &
            rotat(1,2)*rotat(2,1)*rotat(3,3) +  &
            rotat(1,2)*rotat(2,3)*rotat(3,1) +  &
            rotat(1,3)*rotat(2,1)*rotat(3,2) -  &
            rotat(1,3)*rotat(2,2)*rotat(3,1)
   
      if (abs(det) < 1.d-10) then
         rmsok = .false.
         goto 99
      end if
   
      ! If the determinant is negative, invert all elements to
      ! obtain rotation transformation (excluding inversion)
      if (det < 0) then
         do i=1,3
            b(i,ismall) = -b(i,ismall)
         end do
         goto 200
      end if

      ! Rotate reference coordinates: modify xc().
      do i=1,natom
         i3 = 3*i-3
         tmpx = rotat(1,1)*xc(i3+1) + rotat(1,2)*xc(i3+2) + rotat(1,3)*xc(i3+3)
         tmpy = rotat(2,1)*xc(i3+1) + rotat(2,2)*xc(i3+2) + rotat(2,3)*xc(i3+3)
         tmpz = rotat(3,1)*xc(i3+1) + rotat(3,2)*xc(i3+2) + rotat(3,3)*xc(i3+3)
         xc(i3+1) = tmpx
         xc(i3+2) = tmpy
         xc(i3+3) = tmpz
      end do


      ! calculate current rmsd for rms region (itgtrmsgp)
      totmass = 0.d0      ! totmass for rms region
      rmsdvalue = 0.d0


      do ii=1,nattgtrms
         i = rmsgroup(ii)
         i3 = 3*i-3

         totmass = totmass + mass(i)

         tmpx = xc(i3+1) - x(i3+1)
         tmpy = xc(i3+2) - x(i3+2)
         tmpz = xc(i3+3) - x(i3+3)

         rmsdvalue = rmsdvalue + mass(i)*(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
      end do


      if( ineb>0 ) then
         neb_rotat = rotat
      else
         ! Now we need to translate both simulation coordinates X and reference
         ! coordinates XC back: remove the centering done previously. We need to
         ! do this to reference coordinates so they are still fit to simulation
         ! coordinates for energy/force calculation (in routines xconst() and 
         ! xtgtmd()) - that means ADD the *simulation* coords CM.
         do i=1,natom
            i3 = 3*i-3
            x(i3+1) = x(i3+1) + xcm1
            x(i3+2) = x(i3+2) + ycm1
            x(i3+3) = x(i3+3) + zcm1
            xc(i3+1) = xc(i3+1) + xcm1  ! yes, it's CM1 (not CM2)
            xc(i3+2) = xc(i3+2) + ycm1
            xc(i3+3) = xc(i3+3) + zcm1
         end do
      end if
  
   else 
      ! ntr=1 -> don't fit just calculate rmsdvalue. This does not
      ! change X() or XC().
   
      ! calculate current rmsd for rms region (itgtrmsgp)
      rmsdvalue = 0.d0
      totmass = 0.d0
      do ii=1,nattgtrms
         i = rmsgroup(ii)
         i3 = 3*i-3
         totmass = totmass + mass(i)
         tmpx = xc(i3+1) - x(i3+1)
         tmpy = xc(i3+2) - x(i3+2)
         tmpz = xc(i3+3) - x(i3+3)
         rmsdvalue = rmsdvalue + mass(i)*(tmpx*tmpx + tmpy*tmpy + tmpz*tmpz)
      end do
      
   end if
   
   rmsdvalue = sqrt(rmsdvalue/totmass)
   
   return

99 rmsdvalue=-99.
   rmsok = .false.
   return

end subroutine rmsfit
!-----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine bellyf here]
subroutine bellyf(nat,igrp,f)

   implicit _REAL_ (a-h,o-z)
   dimension igrp(*),f(*)
   data zero/0.0d0/
   do i = 1,nat
      if(igrp(i) == 0) then
         i3 = 3*i-3
         f(i3+1) = zero
         f(i3+2) = zero
         f(i3+3) = zero
      end if
   end do
   return
end subroutine bellyf
