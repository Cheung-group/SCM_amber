!<compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"

!----------------------------------------------------------------------------
!     --- GET_NB_ENERGY --
!----------------------------------------------------------------------------
!     ...the main routine for non bond energy (vdw and hbond)
!     as well as direct part of ewald sum. It is structured for parallelism.
!----------------------------------------------------------------------------
subroutine get_nb_energy(iac,ico,ntypes,charge, &
      cn1,cn2,asol,bsol,force,numatoms, &
      nucgrd, &
      ipairs, &
      ewaldcof,eedtbdns,eed_cub,eed_lin, &
      maxnblst,eelt,evdw,ehb,dir_vir,eedvir, &
      filter_cut,ee_type,eedmeth,dxdr, &
      epol,dipole,field,mpoltype)

   use nblist, only: imagcrds,bckptr,nlogrid,nhigrid,numvdw,numhbnd, &
                      myindexlo,myindexhi,numimg, numsc
   use constants, only : zero
   use stack
#ifdef MPI /* SOFT CORE */
   use softcore, only : sc_ener
#endif
   implicit none
   character(kind=1,len=13) :: routine="get_nb_energy"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
#  include "flocntrl.h"
#  include "def_time.h"
#  include "sgld.h"
! Antonios added
#  include "HB.h"
#  include "CROWD.h"
#  include "CHI.h"
! Antonios end   
   
   integer l_real_df,l_real_x,l_real_y,l_real_z,l_real_r2,l_int
   integer numatoms,maxnblst,mpoltype
   integer iac(*),ico(*),ntypes,nucgrd,ee_type,eedmeth
   _REAL_ charge(*),cn1(*),cn2(*),asol(*),bsol(*)
   _REAL_ ewaldcof,eedtbdns,dxdr, &
         eed_cub(4,*),eed_lin(2,*),dir_vir(3,3)
   integer ipairs(maxnblst)
   _REAL_ force(3,numatoms),eelt,epol,evdw,ehb
   _REAL_ eedvir,filter_cut, &
         dipole(3,*),field(3,*)
   
   integer index,numpack,i,k,l,ncell_lo,ncell_hi,ntot,nvdw,nhbnd
   _REAL_ xk,yk,zk
   integer nn,ncache
   
   !     ---FLOW CONTROL FLAG (debug and future respa)
   
   if ( do_dir == 0 )return
   
   eelt = zero
   epol = zero
   evdw = zero
! Antonios added
   EHBV = zero
   EHBA = zero
   ECHI = zero
   ENPC = zero
   ENCC = zero
! Antonios end

   ehb = zero
   eedvir = zero
#ifdef MPI /* SOFT CORE */
   sc_ener(7)=0.0d0
   sc_ener(8)=0.0d0
#endif
   dir_vir(1:3,1:3) = zero
   numpack = 1
   call timer_start(TIME_SHORT_ENE)
   do index = myindexlo,myindexhi
      if ( numimg(index) > 0 )then
         ncell_lo = nlogrid(index)
         ncell_hi = nhigrid(index)
         do k = ncell_lo,ncell_hi
            i = bckptr(k)
            xk = imagcrds(1,k)
            yk = imagcrds(2,k)
            zk = imagcrds(3,k)
#ifdef MPI /* SOFT CORE */
            ! SOFT CORE contribution in numsc
            ntot = numvdw(i) + numhbnd(i) + numsc(i)
#else
            ntot = numvdw(i) + numhbnd(i)
#endif
            nvdw = numvdw(i)
            nhbnd = numhbnd(i)
! Pengzhi
!            if ( ntot > 0 )then
               if ( mpoltype == 0 )then
                  ! allocate 6 temporary caches for performance optimizations
#ifdef MPI /* SOFT CORE */
                  if (nhbnd > numsc(i)) then
#endif
                     ncache = max( nvdw, numhbnd(i) )
#ifdef MPI /* SOFT CORE */
                  else
                     ncache = max( nvdw, numsc(i) )
                  end if
#endif
                  call get_stack(l_real_df,ncache,routine)
                  call get_stack(l_real_x,ncache,routine)
                  call get_stack(l_real_y,ncache,routine)
                  call get_stack(l_real_z,ncache,routine)
                  call get_stack(l_real_r2,ncache,routine)
                  call get_istack(l_int,ncache,routine)
                  if(.not. rstack_ok)then
                     deallocate(r_stack)
                     allocate(r_stack(1:lastrst),stat=alloc_ier)
                     call reassign_rstack(routine)
                  endif
                  if(.not. istack_ok)then
                     deallocate(i_stack)
                     allocate(i_stack(1:lastist),stat=alloc_ier)
                     call reassign_istack(routine)
                  endif
                  REQUIRE(rstack_ok)
                  REQUIRE(istack_ok)
                  call short_ene(i,xk,yk,zk,ipairs(numpack),ntot,nvdw,nhbnd, &
                        ewaldcof,eedtbdns, &
                        eed_cub,eed_lin,charge, &
                        ntypes,iac,ico,cn1,cn2,asol,bsol,filter_cut, &
                        eelt,evdw,ehb,force,dir_vir,ee_type,eedmeth,dxdr, &
                        eedvir,r_stack(l_real_df),r_stack(l_real_x), &
                        r_stack(l_real_y),r_stack(l_real_z), &
                        r_stack(l_real_r2),i_stack(l_int) )
                  call free_stack(l_real_r2,routine)
                  call free_stack(l_real_z,routine)
                  call free_stack(l_real_y,routine)
                  call free_stack(l_real_x,routine)
                  call free_stack(l_real_df,routine)
                  call free_istack(l_int,routine)
               else if ( mpoltype == 1 )then
                  call short_ene_dip(i,xk,yk,zk,ipairs(numpack),ntot,nvdw, &
                        ewaldcof,eedtbdns, &
                        eed_cub,eed_lin,charge,dipole, &
                        ntypes,iac,ico,cn1,cn2,asol,bsol,filter_cut, &
                        eelt,epol,evdw,ehb,force,field,dir_vir, &
                        ee_type,eedmeth,dxdr,eedvir)
               end if
               numpack = numpack + ntot
!            end if  ! ( ntot > 0 )
! Pengzhi end
         end do  !  k = ncell_lo,ncell_hi
      end if  ! ( numimg(k) > 0 )
   end do  !  index = myindexlo,myindexhi

!9.2.2010
!Antonios added the calculation of the triple scalar product that
!         accounts for the correct chirality in amino acid pairs.

! Margaret
! add chiral term
! U(chiral)=0.5*xkchi*(triple-triple_o)^2
! triple is the triple scaler product
! triple=(AxB)*C, where A=CB-CA,B=NC-CA, C=CC-CA
! CA=xi,CB=x(i+1),NC=x(i+2),CC=x(i+3)
! Ax=x(i+1)-xi
! Ay=y(i+1)-yi
! Az=z(i+1)-zi
! Bx=x(i+2)-xi
! By=y(i+2)-yi
! Bz=z(i+2)-zi
! Cx=x(i+3)-xi
! Cy=y(i+3)-yi
! Cz=z(i+3)-zi
! F_xi=-dU(chiral)/dxi
!     =-dU/dAx * dAx/dxi - dU/dBx* dBx/dxi - dU/dCx* dCx/dxi
!     = dU/dAx+dU/dBx+dU/dCx
! F_xi1=-dU(chiral)/dxi1
!       = -dU/dAx * dAx/dxi1
!       = -dU/dAx 
! F_xi+2=-dU(chiral)/dxi2
!       = -dU/dBx * dBx/dxi2
!       = -dU/dBx 
! F_xi+3=-dU(chiral)/dxi3
!       = -dU/dCx * dCx/dxi3
!       = -dU/dCx 

! Antonios end
   call timer_stop(TIME_SHORT_ENE)
   return
!AS mpi-debug write(72,*)evdw Output is partial
end subroutine get_nb_energy 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate the direct Ewald component of the potentials.
!
!-----------------------------------------------------------------------
!     --- SHORT_ENE ---
!-----------------------------------------------------------------------
! short_ene subroutine modified by Nathalie Godbout, sgi 04/26/00.
! additional optimizations aimed at IA32 SSE2, Scott Brozell, TSRI Oct 2002.
! more tweaks, Ross Walker, TSRI 2005.

subroutine short_ene(i,xk,yk,zk,ipairs,ntot,nvdw,nhbnd, &
      ewaldcof,eedtbdns, &
      eed_cub,eed_lin,charge, &
      ntypes,iac,ico,cn1,cn2,asol,bsol,filter_cut, &
      eelt,evdw,ehb,force,dir_vir, &
      ee_type,eedmeth,dxdr,eedvir, &
      cache_df,cache_x,cache_y,cache_z,cache_r2,cache_bckptr )
   use nblist, only: imagcrds,bckptr,tranvec,cutoffnb,volume
   use constants, only : zero, one, two, half, third, TWOPI
   use pimd_vars, only : ipimd,ineb, nrg_all,nbead,nbead_inv
   use decomp, only: decpr, decpair
#ifdef MPI /* SOFT CORE */
   use softcore, only : scalpha, sigma6, foureps, sc_dvdl, isProcessV1, &
        sc_ener, nsc, oneweight, weight0, weight1
#endif
   implicit none
   _REAL_ xk,yk,zk
   integer i,nvdw,nhbnd,ntot
   integer ipairs(*),ee_type,eedmeth
   _REAL_ ewaldcof,eed_cub(*),eed_lin(2,*), &
         charge(*),dir_vir(3,3),eedvir
   _REAL_ eedtbdns,filter_cut,dxdr
   integer ntypes,iac(*),ico(*)
   _REAL_ cn1(*),cn2(*),asol(*),bsol(*), &
         eelt,evdw,ehb,force(3,*)
   integer ic,j,m,n,ind,iaci
   _REAL_ del,delrinv
   _REAL_ switch,d_switch_dx
   _REAL_ b0,b1
   _REAL_ ee_vir_iso
   _REAL_ filter_cut2,xx
   _REAL_ comm1
   _REAL_ xktran(3,18)
   _REAL_ e3dx,e4dx
#ifdef TVDW
   _REAL_ r4,r6pinv
#endif
   ! SOFT CORE
   _REAL_ rfour,denom, denom2, denom3
   _REAL_ ecur
   integer, parameter :: mask27 = 2**27 - 1
   integer nprefetch
#ifdef LES
#  include "les.h"
#endif
   _REAL_ delx,dely,delz,delr, delr2,    f10,   r10,  cgi, &
        delr2inv, r6,    f6,   f12, df,       dfee,  dx,   x, &
        dfx,      vxx,   vxy,  vxz, dfy,      vyy,   vyz,  dumy, &
        dfz,      vzz,   dumz, dumx
#ifdef LES
#  include "ew_cntrl.h"
#endif
! Antonios added
#  include "HB.h"
#  include "CROWD.h"
#  include "CHI.h"
! Antonios end   
#include "md.h"
#include "files.h"
#include "sgld.h"
   _REAL_ uips,uips2,uips6,twou,twou2,twou6,twou12
   _REAL_ pipse,dpipse,eipse,dipse,pvc,pvcu,dvcu,pva,pvau,dvau

   integer itran

   ! variables for conditionally cached data.
   _REAL_ cache_df(*),cache_x(*),cache_y(*)
   _REAL_ cache_z(*),cache_r2(*)
   integer im_new,icount
   integer cache_bckptr(*)


   vxx = zero
   vxy = zero
   vxz = zero
   vyy = zero
   vyz = zero
   vzz = zero
   ee_vir_iso = zero
   del = one / eedtbdns
   dumx = zero
   dumy = zero
   dumz = zero
   filter_cut2 = filter_cut*filter_cut
   cgi = charge(i)
   iaci = ntypes * (iac(i) - 1)
! Qian change,dextran has 2 typies beads for crowders
!   NPP  = ntypes - 1
   NPP  = ntypes - 2
! end   
   do m=1,18
      xktran(1,m) = tranvec(1,m) - xk
      xktran(2,m) = tranvec(2,m) - yk
      xktran(3,m) = tranvec(3,m) - zk
   end do

#ifdef LES
   lestmp=nlesty*(lestyp(i)-1)
#endif
   
   ! The "eedmeth" decision is unrolled here, since this provides
   ! significant efficiency improvements, in spite of greatly
   ! increasing the code length.
   ! Each eedmeth case calculates the 12-6 LJ terms and 12-10 LJ terms
   ! in separate sets of fissioned loops.
   ! Each loop set consists of three loop fission fragments:
   ! prologue, electrostatic, and epilogue.
   ! The prologue loop computes delta r**2 and conditionally caches it,
   ! delta x, delta y, delta z, and the bckptr index.
   ! The electrostatic loop calculates and caches the direct Ewald sum.
   ! The epilogue loop computes the van der Waals interaction and updates
   ! the forces using the previously cached data.
   ! This substantial code tuning is designed to enable software
   ! pipelining of loops for superscalar and similar architectures.
   
   ! For eedmeth = 1 12-6 LJ, the most common case, the electrostatic loop
   ! is fissioned into two fragments. The first loop calculates and
   ! caches principally the reciprocal square root of delta r**2.
   ! The second loop completes the electrostatic work.
   ! This code tuning is designed to enable SIMD vectorization of the
   ! reciprocal square root on IA32 SSE2 compatible platforms.
   

   if ( eedmeth == 1 )then
      
      !-------------------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 1
      !-------------------------------------------------------------
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      !
      !  calculation starts: loop over the data gathered in the temporary
      !  array caches.
      call vdinvsqrt( icount, cache_r2, cache_df )
      cache_r2(1:icount) = dxdr*cache_df(1:icount)*cache_r2(1:icount)
      !
      ! the df cache now stores the reciprocal of delta r, variable delrinv.
      ! the r2 cache contains delta r times the derivative;
      ! this product, dxdr*delr, is denoted as variable x below.
      !
      do im_new = 1,icount
         j = cache_bckptr(im_new)
         delrinv = cache_df(im_new)
         x = cache_r2(im_new)
         !
         ! -- cubic spline on switch:
         !
         ind = eedtbdns*x
         dx = x - ind*del
         ind = 4*ind

         e3dx = dx*eed_cub(3+ind)
         e4dx = dx*dx*eed_cub(4+ind)
         switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
               (e3dx + e4dx*third)*half)
         d_switch_dx = eed_cub(2+ind) + e3dx+ e4dx*half
         
         !---TD Got the idea for B_l from Walter Smith's CCP5 article 1982
         !   Ewald for point multipoles
         
         b0 = switch*delrinv
         b1 = b0 - d_switch_dx*dxdr
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0 ) then
            ! ---if we are using PME, then the correction for lfac will
            !    be done after the reciprocal space calculation is done,
            !    so no need for it here
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1

#ifdef LES
# include "ene_decomp.h"
#endif

#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee
#endif
         delr2inv = delrinv*delrinv
         dfee = dfee*delr2inv

         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee
      end do  !  im_new = 1,icount

      if( tvips )then
         ! Use IPS for L-J energy:
#        include "ips_lj.h"
      else
         ! regular epilogue:
#        include "ew_directe.h"
      end if
      
      !--------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 1
      !--------------------------------------------------------
      
      icount = 0
      do m = nvdw+1,nvdw+nhbnd
#     include "ew_directp.h"
      end do
      
      call vdinvsqrt( icount, cache_r2, cache_df )
      cache_r2(1:icount) = dxdr*cache_df(1:icount)*cache_r2(1:icount)

      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delrinv = cache_df(im_new)
         x = cache_r2(im_new)
         
         !           -- cubic spline on switch:
         
         ind = eedtbdns*x
         dx = x - ind*del
         ind = 4*ind

         e3dx = dx*eed_cub(3+ind)
         e4dx = dx*dx*eed_cub(4+ind)
         switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
               (e3dx + e4dx*third)*half)
         d_switch_dx = eed_cub(2+ind) + e3dx+ e4dx*half
         
         !---TD Got the idea for B_l from Walter Smith's CCP5 article 1982
         !   Ewald for point multipoles
         
         b0 = switch*delrinv
         b1 = b0 - d_switch_dx*dxdr
#ifdef LES
         if (use_pme /= 0 .and.ipimd==0.and.ineb==0 ) then
            
            !---if we are using PME, then the correction for lfac will
            !   be done after the reciprocal space calculation is done,
            !   so no need for it here
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0  
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1

#ifdef LES
# include "ene_decomp.h"
#endif


#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee
#endif
         delr2inv = delrinv*delrinv
         dfee = dfee*delr2inv
         
#ifdef HAS_10_12
         cache_r2(im_new)=delr2inv
#endif
         cache_df(im_new)=dfee

      end do  !  im_new = 1, icount

#     include "ew_directe2.h"

#ifdef MPI /* SOFT CORE */
      !--------------------------------------------------------
      !     Now loop over the softcore (modified 12-6) LJ terms for eedmeth = 1
      !--------------------------------------------------------
      
      icount = 0
      ! run over the 3rd subset of the pairlist
      do m = nvdw+nhbnd+1,ntot      
#     include "ew_directp.h"
      end do
      
      call vdinvsqrt( icount, cache_r2, cache_df )
      cache_r2(1:icount) = dxdr*cache_df(1:icount)*cache_r2(1:icount) ! cache_r2 is saved for the vdW interactions below

      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delrinv = cache_df(im_new)
         x = cache_r2(im_new)
         
         !           -- cubic spline on switch:
         
         ind = eedtbdns*x
         dx = x - ind*del
         ind = 4*ind

         e3dx = dx*eed_cub(3+ind)
         e4dx = dx*dx*eed_cub(4+ind)
         switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
               (e3dx + e4dx*third)*half)
         d_switch_dx = eed_cub(2+ind) + e3dx+ e4dx*half
         
         !---TD Got the idea for B_l from Walter Smith's CCP5 article 1982
         !   Ewald for point multipoles
         
         b0 = switch*delrinv
         b1 = b0 - d_switch_dx*dxdr
#if defined(LES) 
         if (use_pme /= 0 .and. ipimd.eq.0 ) then
            
            !---if we are using PME, then the correction for lfac will
            !   be done after the reciprocal space calculation is done,
            !   so no need for it here
            
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
!#   ifdef PIMD 
!           if ( ievb /= 0 .and. ipimd.ne.0 ) then
!            if(cnum(i).eq.0.and.cnum(j).eq.0) then
!               nrg_ele(1:nbead)=nrg_ele(1:nbead) + comm1*b0
!            else 
!               if(cnum(i).ne.0) then
!                  nrg_ele(cnum(i)) = nrg_ele(cnum(i)) + comm1*b0
!               else
!                  nrg_ele(cnum(j)) = nrg_ele(cnum(j)) + comm1*b0
!               end if
!            end if
!
!#   endif
         end if
#else
         comm1 = cgi*charge(j)
#endif
         
         eelt = eelt + comm1*b0
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,comm1*b0/(nstlim/ntpr))
         dfee = comm1*b1
#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee
#endif
         delr2inv = delrinv*delrinv
         dfee = dfee*delr2inv
         
         cache_r2(im_new)=cache_x(im_new)*cache_x(im_new)+cache_y(im_new)*cache_y(im_new)+ &
              cache_z(im_new)*cache_z(im_new)     ! inefficient, change this to use another cache later. 
                                                  ! The following ew_directe3.h or directe4.h needs r^2 in cache_r2
         cache_df(im_new)=dfee

      end do  !  im_new = 1, icount

      ! V1 uses ew_directe3.h, in which softcore atoms are treated as appearing,
      ! i.e. fully interacting at lambda=1 and 'soft' at small lambda
      ! V0 uses ew_directe4.h, in which softcore atoms are treated as vanishing, 
      ! i.e. fully interacting at lambda=0 and 'soft' at large lambda

      if ( isProcessV1 ) then
#include "ew_directe3.h"
      else
#include "ew_directe4.h"
      end if

#endif /* MPI for SOFT CORE  */

   else if ( eedmeth == 2 )then
      !--------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 2
      !--------------------------------------------------
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      !  calculation starts: loop over the data gathered in the temporary
      !  array caches.
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1,icount
         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         !---- linear lookup on switch:
         delrinv = one/sqrt(delr2)
         delr = delr2*delrinv
         delr2inv = delrinv*delrinv
         x = dxdr*delr
         xx = eedtbdns*x + 1
         ind = xx
         dx = xx - ind
         switch = (one - dx)*eed_lin(1,ind) + &
               dx*eed_lin(1,ind+1)
         d_switch_dx = (one - dx)*eed_lin(2,ind) + &
               dx*eed_lin(2,ind+1)
         b0 = switch*delrinv
         b1 = (b0 - d_switch_dx*dxdr)*delr2inv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            !---if we are using PME, then the correction for lfac will
            !   be done after the reciprocal space calculation is done,
            !   so no need for it here
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0       
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1

#ifdef LES
# include "ene_decomp.h"
#endif


#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee*delr2
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do  !  im_new = 1,icount

      if( tvips )then
         ! Use IPS for L-J energy:
#        include "ips_lj.h"
      else
         ! regular epilogue:
#        include "ew_directe.h"
      end if

      !---------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 2
      !---------------------------------------------------------
      icount = 0
      do m = nvdw+1,ntot
#     include "ew_directp.h"
      end do
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !            -- linear lookup on switch:
         
         delrinv = one/sqrt(delr2)
         delr = delr2*delrinv
         delr2inv = delrinv*delrinv
         x = dxdr*delr
         xx = eedtbdns*x + 1
         ind = xx
         dx = xx - ind
         switch = (one - dx)*eed_lin(1,ind) + &
               dx*eed_lin(1,ind+1)
         d_switch_dx = (one - dx)*eed_lin(2,ind) + &
               dx*eed_lin(2,ind+1)
         b0 = switch*delrinv
         b1 = (b0 - d_switch_dx*dxdr)*delr2inv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            !---if we are using PME, then the correction for lfac will
            !   be done after the reciprocal space calculation is done,
            !   so no need for it here
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0        
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1
#ifdef LES
# include "ene_decomp.h"
#endif

#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee*delr2
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do  !  im_new = 1, icount

#     include "ew_directe2.h"

   else if ( eedmeth == 3 )then

      !---------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 3
      !---------------------------------------------------
      
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      
      !  calculation starts: loop over the data gathered in the temporary
      !  array caches.
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1,icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !            -- explicit function call:
         
         delrinv = one/sqrt(delr2)
         delr = delr2*delrinv
         delr2inv = delrinv*delrinv
         x = dxdr*delr
         call get_ee_func(x,switch,d_switch_dx,ee_type)
         
         b0 = switch*delrinv
         b1 = (b0 - d_switch_dx*dxdr)*delr2inv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            
            !             ---if we are using PME, then the correction for lfac will
            !                be done after the reciprocal space calculation is done,
            !                so no need for it here
            
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0        
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1
#ifdef LES
# include "ene_decomp.h"
#endif

#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee*delr2
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee
      end do  !  im_new = 1,icount

      if( tvips )then
         ! Use IPS for L-J energy:
#        include "ips_lj.h"
      else
         ! regular epilogue:
#        include "ew_directe.h"
      end if
      
      !---------------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 3
      !---------------------------------------------------------------
      
      icount = 0
      do m = nvdw+1,ntot
#     include "ew_directp.h"
      end do
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !            -- explicit function call:
         
         delrinv = one/sqrt(delr2)
         delr = delr2*delrinv
         delr2inv = delrinv*delrinv
         x = dxdr*delr
         call get_ee_func(x,switch,d_switch_dx,ee_type)
         
         b0 = switch*delrinv
         b1 = (b0 - d_switch_dx*dxdr)*delr2inv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            
            !             ---if we are using PME, then the correction for lfac will
            !                be done after the reciprocal space calculation is done,
            !                so no need for it here
            
            comm1 = cgi*charge(j)
         else
            lfac=lesfac(lestmp+lestyp(j))
            comm1 = cgi*charge(j)*lfac
         end if
#else
         comm1 = cgi*charge(j)
#endif
         ecur = comm1*b0        
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = comm1*b1
#ifdef LES
# include "ene_decomp.h"
#endif

#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfee*delr2
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do  !  im_new = 1, icount

#     include "ew_directe2.h"

   else if ( eedmeth == 4 )then

      !-------------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 4
      !-------------------------------------------------------
      
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      
      !  calculation starts: loop over the data gathered in the temporary
      !  array caches.
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1,icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !             -- don't use a switch:
         !             -- straight Coulomb
         
         delrinv = one/sqrt(delr2)
         delr2inv = delrinv*delrinv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            b0 = cgi*charge(j)*delrinv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delrinv
         end if
#else
         b0 = cgi*charge(j)*delrinv
#endif
         ecur = b0
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = b0*delr2inv
        
#ifdef LES
# include "ene_decomp.h"
#endif 
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee
      end do

      if( tvips )then
         ! Use IPS for L-J energy:
#        include "ips_lj.h"
      else
         ! regular epilogue:
#        include "ew_directe.h"
      end if
      
      !--------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 4
      !--------------------------------------------------------
      
      icount = 0
      do m = nvdw+1,nvdw+nhbnd
#     include "ew_directp.h"
      end do
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !             -- don't use a switch:
         !             -- straight Coulomb
         
         delrinv = one/sqrt(delr2)
         delr2inv = delrinv*delrinv
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            b0 = cgi*charge(j)*delrinv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delrinv
         end if
#else
         b0 = cgi*charge(j)*delrinv
#endif
         ecur = b0
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = b0*delr2inv
#ifdef LES
# include "ene_decomp.h"
#endif         
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do

#     include "ew_directe2.h"

#ifdef MPI /* SOFT CORE */
      !--------------------------------------------------------
      !     Now loop over the softcore (modified 12-6) LJ terms for eedmeth = 4
      !--------------------------------------------------------
      
      icount = 0
      do m = nvdw+nhbnd+1,ntot
#     include "ew_directp.h"
      end do
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)

         !             -- don't use a switch:
         !             -- straight Coulomb
         
         delrinv = one/sqrt(delr2)
         delr2inv = delrinv*delrinv
#if defined(LES) 
         if (use_pme /= 0.and.ipimd.eq.0) then
            b0 = cgi*charge(j)*delrinv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delrinv
!#   ifdef PIMD
!            if(cnum(i).eq.0.and.cnum(j).eq.0) then
!               nrg_ele(1:nbead)=nrg_ele(1:nbead) + comm1*b0
!            else 
!               if(cnum(i).ne.0) then
!                  nrg_ele(cnum(i)) = nrg_ele(cnum(i)) + comm1*b0
!               else
!                  nrg_ele(cnum(j)) = nrg_ele(cnum(j)) + comm1*b0
!               end if
!            end if
!
!#   endif
         end if
#else
         b0 = cgi*charge(j)*delrinv
#endif
         eelt = eelt + b0
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,b0/(nstlim/ntpr))
         dfee = b0*delr2inv

         cache_r2(im_new)=delr2                                               
         ! Contrary to the 12-6 and 12-10 cases above, cache_r2 contains r^2 here
         cache_df(im_new)=dfee

      end do

      ! V1 uses ew_directe3.h, in which softcore atoms are treated as appearing,
      ! i.e. fully interacting at lambda=1 and 'soft' at small lambda
      ! V0 uses ew_directe4.h, in which softcore atoms are treated as vanishing, 
      ! i.e. fully interacting at lambda=0 and 'soft' at large lambda

      if ( isProcessV1 ) then
#include "ew_directe3.h"
      else
#include "ew_directe4.h"
      end if

#endif /* MPI for SOFT CORE */

   else if ( eedmeth == 5 )then

      !---------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 5
      !---------------------------------------------------
      
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      
      !  calculation starts: loop over the data gathered in the temporary
      !  array caches.
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1,icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !             -- use dielectric of 1/r:
         
         delr2inv = one/delr2
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            b0 = cgi*charge(j)*delr2inv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delr2inv
         end if
#else
         b0 = cgi*charge(j)*delr2inv
#endif 
         ecur = b0
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = two*b0*delr2inv

#ifdef LES
# include "ene_decomp.h"
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee
      end do

      if( tvips )then
         ! Use IPS for L-J energy:
#        include "ips_lj.h"
      else
         ! regular epilogue:
#        include "ew_directe.h"
      end if
      
      !---------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 5
      !---------------------------------------------------------
      
      icount = 0
      do m = nvdw+1,ntot
#     include "ew_directp.h"
      end do
      
      ! SGI compiler directive to prevent compiler loop fusioning.
      !*$* NO FUSION
      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         
         !             -- use dielectric of 1/r:
         
         delr2inv = one/delr2
#ifdef LES
         if (use_pme /= 0 .and. ipimd==0) then
            b0 = cgi*charge(j)*delr2inv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delr2inv
         end if
#else
         b0 = cgi*charge(j)*delr2inv
#endif
         ecur = b0
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = two*b0*delr2inv
#ifdef LES
# include "ene_decomp.h"
#endif
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do

#     include "ew_directe2.h"

   else if ( eedmeth == 6 )then

      !-------------------------------------------------------
      !     Loop over the 12-6 LJ terms for eedmeth = 6
      !-------------------------------------------------------
      
      icount = 0
      do m = 1,nvdw
#     include "ew_directp.h"
      end do
      
      call vdinvsqrt( icount, cache_r2, cache_df )
      ! the df cache now stores the reciprocal of delta r, variable delrinv.

      if( teips .or. tvips ) nnbips=nnbips+icount*2

      do im_new = 1,icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         delrinv = cache_df(im_new)
         delr2inv = delrinv*delrinv

         !  -- use a ipse long range potential:
         !          eipse=e0+(b0+b1r^2+b2r^4+b3r^6)/sqrt(2-r^2)
         ! compare: ele=1/r  fele=-1/r^2

         uips=ripsinv*delr2*delrinv
         uips2=delr2*rips2inv
         twou2=1.d0/(2.0d0-uips2)
         twou=sqrt(twou2)
#ifdef LES
         if (use_pme /= 0.and.ipimd==0) then
            b0 = cgi*charge(j)*delrinv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delrinv
         end if
#else
         b0 = cgi*charge(j)*delrinv
#endif
         pipse = bipse0 + uips2*(bipse1 + uips2*(bipse2 + uips2*bipse3))
         dpipse = 2.0d0*bipse1 + uips2*(4.0d0*bipse2 + 6.0d0*uips2*bipse3)
         dipse = uips*(dpipse + pipse*twou2)*twou*rips2inv
         eipse = b0*uips*(pipse*twou - pipsec)
         ecur = b0 + eipse
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = b0*(delr2inv - dipse)
#ifdef LES
# include "ene_decomp.h"
#endif         
         cache_r2(im_new)=delr2inv
         cache_df(im_new)=dfee

      end do

      if(TVIPS)then
         ! Use IPS for L-J energy
#        include "ips_lj.h"
      else
         ! epilogue
#        include "ew_directe.h"
      endif

      !--------------------------------------------------------
      !     Now loop over the 12-10 LJ terms for eedmeth = 6
      !--------------------------------------------------------
      
      icount = 0
      do m = nvdw+1,ntot
#     include "ew_directp.h"
      end do
      call vdinvsqrt( icount, cache_r2, cache_df )
      
      if( teips .or. tvips ) nnbips=nnbips+icount*2

      do im_new = 1, icount

         j = cache_bckptr(im_new)
         delr2 = cache_r2(im_new)
         delrinv = cache_df(im_new)
         delr2inv = delrinv*delrinv
#ifdef LES
         if (use_pme /= 0.and.ipimd==0) then
            b0 = cgi*charge(j)*delrinv
         else
            lfac=lesfac(lestmp+lestyp(j))
            b0 = cgi*charge(j)*lfac*delrinv
         end if
#else
         b0 = cgi*charge(j)*delrinv
#endif

         !  -- use a ipse long range potential:
         !           eipse=e0+(b0+b1r^2+b2r^4+b3r^6)/sqrt(2-r^2)
         !  compare: ele=1/r  fele=-1/r^2

         uips = ripsinv*delr2*delrinv
         uips2 = delr2*rips2inv
         twou2 = 1.d0/(2.0d0 - uips2)
         twou = sqrt(twou2)
         pipse = bipse0 + uips2*(bipse1 + uips2*(bipse2 + uips2*bipse3))
         dpipse = 2.0d0*bipse1 + uips2*(4.0d0*bipse2 + 6.0d0*uips2*bipse3)
         dipse = uips*(dpipse + pipse*twou2)*twou*rips2inv
         eipse = b0*uips*(pipse*twou - pipsec)
         ecur = b0 + eipse
         eelt = eelt + ecur
         ! -- ti decomp
         if(decpr .and. idecomp > 0) call decpair(2,i,j,ecur/(nstlim/ntpr))
         dfee = b0*(delr2inv - dipse)
#ifdef LES
# include "ene_decomp.h"
#endif

#ifdef HAS_10_12
         cache_r2(im_new)=delr2inv
#endif
         cache_df(im_new)=dfee

      end do

#     include "ew_directe2.h"

   end if  ! ( eedmeth )

   force(1,i) = force(1,i) - dumx
   force(2,i) = force(2,i) - dumy
   force(3,i) = force(3,i) - dumz
#ifndef noVIRIAL
   dir_vir(1,1) = dir_vir(1,1) + vxx
   dir_vir(1,2) = dir_vir(1,2) + vxy
   dir_vir(2,1) = dir_vir(2,1) + vxy
   dir_vir(1,3) = dir_vir(1,3) + vxz
   dir_vir(3,1) = dir_vir(3,1) + vxz
   dir_vir(2,2) = dir_vir(2,2) + vyy
   dir_vir(2,3) = dir_vir(2,3) + vyz
   dir_vir(3,2) = dir_vir(3,2) + vyz
   dir_vir(3,3) = dir_vir(3,3) + vzz
   eedvir = eedvir + ee_vir_iso
#endif
   return
end subroutine short_ene 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate the direct Ewald component of the potentials 
!+     with polarizabilities.
!
!-----------------------------------------------------------------------
!     --- SHORT_ENE_DIP ---
!-----------------------------------------------------------------------

subroutine short_ene_dip(i,xk,yk,zk,ipairs,numtot,numvdw, &
      ewaldcof,eedtbdns, &
      eed_cub,eed_lin,charge,dipole, &
      ntypes,iac,ico,cn1,cn2,asol,bsol,filter_cut, &
      eelt,epol,evdw,ehb,frc,field,dir_vir, &
      ee_type,eedmeth,dxdr,eedvir)

   use nblist, only: bckptr,imagcrds,tranvec
   use constants, only : zero, one, two, three, five, six, twelve, third, half
   implicit none
   _REAL_ xk,yk,zk
   integer i,numvdw,numtot
   integer ipairs(*),ee_type,eedmeth
   _REAL_ ewaldcof,eed_cub(4,*),eed_lin(2,*), &
         charge(*),dipole(3,*),dir_vir(3,3),eedvir
   _REAL_ eedtbdns,filter_cut,dxdr
   integer ntypes,iac(*),ico(*)
   _REAL_ cn1(*),cn2(*),asol(*),bsol(*), &
         eelt,epol,evdw,ehb,frc(3,*),field(3,*)
   integer ic,j,m,n,ind,iaci
   _REAL_ del
   _REAL_ switch,d_switch_dx
   _REAL_ ee_vir_iso
   _REAL_ edx,edy,edz
   _REAL_ b0,b1,b2,b3,fac,dotir,dotjr,dotij,fact
   _REAL_ dphii_dx,dphii_dy,dphii_dz, &
         dphij_dx,dphij_dy,dphij_dz
   _REAL_  term,term0,term1,termi,termj,cgj

   _REAL_ filter_cut2,xx
   _REAL_ xktran(3,18)
   integer, parameter :: mask27 = 2**27 - 1
#ifdef LES
#  include "les.h"
#endif
   _REAL_ delx,dely,delz,delr, delr2,    f10,   r10,  cgi, &
        delr2inv, r6,    f6,   f12, df,       dfee,  dx,   x, &
        dfx,      vxx,   vxy,  vxz, dfy,      vyy,   vyz,  dumy, &
        dfz,      vzz,   dumz, dumx
   integer itran
   
   fac = two*ewaldcof*ewaldcof
   ee_vir_iso = zero
   del = one / eedtbdns
   dumx = zero
   dumy = zero
   dumz = zero
   edx = zero
   edy = zero
   edz = zero
   filter_cut2 = filter_cut*filter_cut
   cgi = charge(i)
   iaci = ntypes * (iac(i) - 1)

   do m=1,18
      xktran(1,m) = tranvec(1,m) - xk
      xktran(2,m) = tranvec(2,m) - yk
      xktran(3,m) = tranvec(3,m) - zk
   end do
#ifdef LES
   lestmp=nlesty*(lestyp(i)-1)
#endif
   do m = 1,numvdw
      n = ipairs(m)
      itran=ishft(n,-27)
      n = iand(n,mask27)
      j = bckptr(n)
      delx = imagcrds(1,n) + xktran(1,itran)
      dely = imagcrds(2,n) + xktran(2,itran)
      delz = imagcrds(3,n) + xktran(3,itran)
      delr2 = delx*delx + dely*dely+delz*delz
      if ( delr2 < filter_cut2 )then
         delr = sqrt(delr2)
         delr2inv = one/delr2
         x = dxdr*delr
         cgj = charge(j)
         if ( eedmeth == 1 )then
            
            !           -- cubic spline on switch
            
            ind = eedtbdns*x + 1
            dx = x - (ind-one)*del
            switch = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                  dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
            d_switch_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                  dx*eed_cub(4,ind)*half)
            
         else if ( eedmeth == 2 )then
            
            !           ---linear lookup on switch, deriv
            
            xx = eedtbdns*x + 1
            ind = xx
            dx = xx - ind
            switch = (one - dx)*eed_lin(1,ind) + &
                  dx*eed_lin(1,ind+1)
            d_switch_dx = (one - dx)*eed_lin(2,ind) + &
                  dx*eed_lin(2,ind+1)
            
         else if ( eedmeth == 3 )then
            
            !           ---direct function call:
            
            call get_ee_func(x,switch,d_switch_dx,ee_type)
            
         else if ( eedmeth == 4 ) then
            
            !            ---use un-modified Coulomb interaction, no switch
            
            switch = one
            d_switch_dx = zero
            
         else
            
            write(6,*) 'bad eedmeth in ew_short_dip: ',eedmeth
            call mexit( 6,1 )
            
         end if  ! ( eedmeth == 1 )
         
         ! TD Got the idea for B_l from Walter Smith's CCP5 article 1982
         ! Ewald for point multipoles
         ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
         ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
         
         b0 = switch*delr*delr2inv
         fact = d_switch_dx*dxdr
         b1 = (b0 - fact)*delr2inv
         fact = fac*fact
         b2 = (three*b1 - fact)*delr2inv
         fact = fac*fact
         b3 = (five*b2 - fact)*delr2inv
         
         ! B1 = (B0 - d_switch_dx*dxdr)*delr2inv
         ! B2 = (Three*B1 - fac*ewaldcof*d_switch_dx)*delr2inv
         ! B3 = (Five*B2 - fac*fac*ewaldcof*d_switch_dx)*delr2inv
         
         ! epol = dip_i dot grad_i of phii
         ! phii is direct sum electrostatic potential at i due to j
         ! so phii = cgj*B0 + dipj dot gradj of B0 = cgj*B0 - dotjr*B1
         ! dphii_dx etc are derivatives with respect to r_i
         ! phij is direct sum electrostatic potential at j due to i
         ! dphij_dx etc are derivatives with respect to r_j
         
         dotjr = dipole(1,j)*delx+dipole(2,j)*dely+dipole(3,j)*delz
         dotir = dipole(1,i)*delx+dipole(2,i)*dely+dipole(3,i)*delz
         dotij = dipole(1,i)*dipole(1,j)+dipole(2,i)*dipole(2,j)+ &
               dipole(3,i)*dipole(3,j)
         
         ! gradi phii = cgj*rij*B1 + dipj*B1 - dotjr*rij*B2
         ! so epol = -cgi*dotjr*B1 + (cgj*B1 - dotjr*B2)*dotir + dotij*B1
         
         eelt = eelt + cgi*cgj*b0
         term = cgj*dotir-cgi*dotjr+dotij
         epol = epol + term*b1 - dotir*dotjr*b2
         term0 = cgi*cgj*b0 + term*b1 - dotir*dotjr*b2
         
         ! so ene = ene + term0; dfx = dterm0_dx etc
         ! grad term0 = term1*rij + B1*grad_i term - grad_i dotir*dotjr*B2
         ! grad_i term = -cgj*dip_i + cgi*dip_j
         ! grad_i dotir = -dip_i; similar for dotjr
         ! grad_i term0 = term1*rij + (-cgj*B1+dotjr*B2)*dip_i +
         !                (cgi*B1+dotir*B2)*dip_j
         
         term1 = cgi*cgj*b1 + term*b2 - dotir*dotjr*b3
         termi = cgi*b1+dotir*b2
         termj = cgj*b1-dotjr*b2
         dfx = term1*delx + termi*dipole(1,j) - termj*dipole(1,i)
         dfy = term1*dely + termi*dipole(2,j) - termj*dipole(2,i)
         dfz = term1*delz + termi*dipole(3,j) - termj*dipole(3,i)
#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfx*delx - dfy*dely - dfz*delz
#endif
         ic = ico(iaci+iac(j))
         r6 = delr2inv*delr2inv*delr2inv
#ifdef LES
         lfac=lesfac(lestmp+lestyp(j))
         f6 = cn2(ic)*r6*lfac
         f12 = cn1(ic)*(r6*r6)*lfac
#else
         f6 = cn2(ic)*r6
         f12 = cn1(ic)*(r6*r6)
#endif
         evdw = evdw + f12 - f6
         
         !         ---force related quantities
         
         df = (twelve*f12 - six*f6)*delr2inv
         dfx = dfx + df*delx
         dfy = dfy + df*dely
         dfz = dfz + df*delz
#ifndef noVIRIAL
         dir_vir(1,1) = dir_vir(1,1) - dfx*delx
         dir_vir(1,2) = dir_vir(1,2) - dfx*dely
         dir_vir(1,3) = dir_vir(1,3) - dfx*delz
         dir_vir(2,1) = dir_vir(2,1) - dfy*delx
         dir_vir(2,2) = dir_vir(2,2) - dfy*dely
         dir_vir(2,3) = dir_vir(2,3) - dfy*delz
         dir_vir(3,1) = dir_vir(3,1) - dfz*delx
         dir_vir(3,2) = dir_vir(3,2) - dfz*dely
         dir_vir(3,3) = dir_vir(3,3) - dfz*delz
#endif
         frc(1,j) = frc(1,j) + dfx
         frc(2,j) = frc(2,j) + dfy
         frc(3,j) = frc(3,j) + dfz
         dumx = dumx + dfx
         dumy = dumy + dfy
         dumz = dumz + dfz
         
         !         ---field related quantities
         
         dphii_dx = termj*delx + b1*dipole(1,j)
         dphii_dy = termj*dely + b1*dipole(2,j)
         dphii_dz = termj*delz + b1*dipole(3,j)
         dphij_dx = -termi*delx + b1*dipole(1,i)
         dphij_dy = -termi*dely + b1*dipole(2,i)
         dphij_dz = -termi*delz + b1*dipole(3,i)
         edx = edx + dphii_dx
         edy = edy + dphii_dy
         edz = edz + dphii_dz
         field(1,j) = field(1,j) - dphij_dx
         field(2,j) = field(2,j) - dphij_dy
         field(3,j) = field(3,j) - dphij_dz
      end if  ! ( delr2 < filter_cut2 )
   end do  !  m = 1,numvdw
   
   do m = numvdw+1,numtot
      n = ipairs(m)
      itran=ishft(n,-27)
      n = iand(n,mask27)
      j = bckptr(n)
      delx = imagcrds(1,n) + xktran(1,itran)
      dely = imagcrds(2,n) + xktran(2,itran)
      delz = imagcrds(3,n) + xktran(3,itran)
      delr2 = delx*delx + dely*dely+delz*delz
      if ( delr2 < filter_cut2 )then
         delr = sqrt(delr2)
         delr2inv = one/delr2
         x = dxdr*delr
         cgj = charge(j)

         if ( eedmeth == 1 )then
            
            !           -- cubic spline on switch
            
            ind = eedtbdns*x + 1
            dx = x - (ind-one)*del
            switch = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                  dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*third)*half)
            d_switch_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                  dx*eed_cub(4,ind)*half)
            
         else if ( eedmeth == 2 )then
            
            !           ---linear lookup on switch, deriv
            
            xx = eedtbdns*x + 1
            ind = xx
            dx = xx - ind
            switch = (one - dx)*eed_lin(1,ind) + &
                  dx*eed_lin(1,ind+1)
            d_switch_dx = (one - dx)*eed_lin(2,ind) + &
                  dx*eed_lin(2,ind+1)
            
         else if ( eedmeth == 3 )then
            
            !           ---direct function call:
            
            call get_ee_func(x,switch,d_switch_dx,ee_type)
            
         else if ( eedmeth == 4 ) then
            
            !            ---use un-modified Coulomb interaction, no switch
            
            switch = one
            d_switch_dx = zero
            
         else
            
            write(6,*) 'bad eedmeth in ew_short_dip: ',eedmeth
            call mexit( 6,1 )
            
         end if  ! ( eedmeth == 1 )
         
         b0 = switch*delr*delr2inv
         fact = d_switch_dx*dxdr
         b1 = (b0 - fact)*delr2inv
         fact = fac*fact
         b2 = (three*b1 - fact)*delr2inv
         fact = fac*fact
         b3 = (five*b2 - fact)*delr2inv
         !         B0 = switch*delr*delr2inv
         !         B1 = (B0 - d_switch_dx*dxdr)*delr2inv
         !         B2 = (Three*B1 - fac*ewaldcof*d_switch_dx)*delr2inv
         !         B3 = (Five*B2 - fac*fac*ewaldcof*d_switch_dx)*delr2inv
         
         dotjr = dipole(1,j)*delx+dipole(2,j)*dely+dipole(3,j)*delz
         dotir = dipole(1,i)*delx+dipole(2,i)*dely+dipole(3,i)*delz
         dotij = dipole(1,i)*dipole(1,j)+dipole(2,i)*dipole(2,j)+ &
               dipole(3,i)*dipole(3,j)
         eelt = eelt + cgi*cgj*b0
         term = cgj*dotir-cgi*dotjr+dotij
         epol = epol + term*b1 - dotir*dotjr*b2
         term0 = cgi*cgj*b0 + term*b1 - dotir*dotjr*b2
         term1 = cgi*cgj*b1 + term*b2 - dotir*dotjr*b3
         termi = cgi*b1+dotir*b2
         termj = cgj*b1-dotjr*b2
         dfx = (term1)*delx + termi*dipole(1,j) - termj*dipole(1,i)
         dfy = (term1)*dely + termi*dipole(2,j) - termj*dipole(2,i)
         dfz = (term1)*delz + termi*dipole(3,j) - termj*dipole(3,i)
#ifndef noVIRIAL
         ee_vir_iso = ee_vir_iso - dfx*delx - dfy*dely - dfz*delz
#endif
#ifdef HAS_10_12
         
         ! --- this code allows 10-12 terms; in many (most?) (all?) cases, the
         !     only "nominal" 10-12 terms are on waters, where the asol and bsol
         !     parameters are always zero; hence we can skip the L-J part; note
         !     that we still have to compute the electrostatic interactions
         
         ic = -ico(iaci+iac(j))
         r10 = delr2inv*delr2inv*delr2inv*delr2inv*delr2inv
         f10 = bsol(ic)*r10
         f12 = asol(ic)*(r10*delr2inv)
         ehb = ehb + f12 - f10
         df = (twelve*f12 - ten*f10)*delr2inv
#else
         df = zero
#endif
         
         !         ---force related quantities
         
         dfx = dfx + df*delx
         dfy = dfy + df*dely
         dfz = dfz + df*delz
#ifndef noVIRIAL
         dir_vir(1,1) = dir_vir(1,1) - dfx*delx
         dir_vir(1,2) = dir_vir(1,2) - dfx*dely
         dir_vir(1,3) = dir_vir(1,3) - dfx*delz
         dir_vir(2,1) = dir_vir(2,1) - dfy*delx
         dir_vir(2,2) = dir_vir(2,2) - dfy*dely
         dir_vir(2,3) = dir_vir(2,3) - dfy*delz
         dir_vir(3,1) = dir_vir(3,1) - dfz*delx
         dir_vir(3,2) = dir_vir(3,2) - dfz*dely
         dir_vir(3,3) = dir_vir(3,3) - dfz*delz
#endif
         frc(1,j) = frc(1,j) + dfx
         frc(2,j) = frc(2,j) + dfy
         frc(3,j) = frc(3,j) + dfz
         dumx = dumx + dfx
         dumy = dumy + dfy
         dumz = dumz + dfz
         
         !         ---field related quantities
         
         dphii_dx = termj*delx + b1*dipole(1,j)
         dphii_dy = termj*dely + b1*dipole(2,j)
         dphii_dz = termj*delz + b1*dipole(3,j)
         dphij_dx = -termi*delx + b1*dipole(1,i)
         dphij_dy = -termi*dely + b1*dipole(2,i)
         dphij_dz = -termi*delz + b1*dipole(3,i)
         edx = edx + dphii_dx
         edy = edy + dphii_dy
         edz = edz + dphii_dz
         field(1,j) = field(1,j) - dphij_dx
         field(2,j) = field(2,j) - dphij_dy
         field(3,j) = field(3,j) - dphij_dz
      end if  ! ( delr2 < filter_cut2 )
   end do  !  m = numvdw+1,numtot
   frc(1,i) = frc(1,i) - dumx
   frc(2,i) = frc(2,i) - dumy
   frc(3,i) = frc(3,i) - dumz
   field(1,i) = field(1,i) - edx
   field(2,i) = field(2,i) - edy
   field(3,i) = field(3,i) - edz
#ifndef noVIRIAL
   eedvir = eedvir + ee_vir_iso
#endif
   return
end subroutine short_ene_dip 

