#include "copyright.h"
#include "dprec.h"
#include "assert.h"
!*********************************************************************
!               SUBROUTINE TRAJENE
!*********************************************************************
! carlos add trajene routine for processing trajectory energies

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trajene here]
subroutine trajene(x,ix,ih,ipairs,ene,ok,qsetup)

   use nblist,only: fill_tranvec,a,b,c!,alpha,beta,gamma,nbflag,cutoffnb
   use constants,only: one,half

   implicit none
   integer ipairs(*)
   _REAL_ carrms

   _REAL_ x(*),ene(*),oldbox(3),rmu(3)
   integer ix(*)
   character(len=4) ih(*)
   logical ok,qsetup,loutfm
   integer member,j

#  include "memory.h"
#  include "tgtmd.h"
! DAN ROE: Is extra.h still needed?
#  include "extra.h"
! needed for ntb, box
#  include "box.h" 
!#  include "ew_cntrl.h"
!#  include "ew_erfc_spline.h"
! needed for nfft
#  include "ew_pme_recip.h"
! Needed for INPTRAJ_UNIT, inptraj
#  include "files.h"

   loutfm = ioutfm <= 0

   ! Open the inptraj file for coordinate reading
   call amopen(INPTRAJ_UNIT,inptraj,'O','F','R')

   read(INPTRAJ_UNIT,100) title
   write (6,100) title
   100 format(a80)

   member=0

   !     loop over trajectory file, exiting only on error or end of file

   do while ( .true. )

      !       --- read next coordinate set from trajectory

      read(INPTRAJ_UNIT,110,end=1000,err=1010) (x(j),j=lcrd,lcrd+natom*3-1)

      ! DAN ROE:
      ! If box coords are present in trajectory (ntb>0), read them in and
      !  update the unit cell and grid size.
      if (ntb>0) then
         ! 1- Save old box coords.
         oldbox(1) = box(1)
         oldbox(2) = box(2)
         oldbox(3) = box(3)
         !write(6,*) "DEBUG: OLDBOX: ",oldbox(1),oldbox(2),oldbox(3)
         !write(6,*) "DEBUG: OLDabc: ",a,b,c

         ! 2- Read in current box coords.
         read(INPTRAJ_UNIT,'(3f8.3)',end=1000,err=1010) box(1), box(2), box(3)
         !write(6,*) "DEBUG: BOX: ",box(1),box(2),box(3)

         ! 3- If the box size has changed, some parameters have to be recalc.
         if (box(1)/=oldbox(1).or.box(2)/=oldbox(2).or.box(3)/=oldbox(3)) then

            ! 3a- Calculate box scaling factors
            rmu(1) = box(1) / oldbox(1)
            rmu(2) = box(2) / oldbox(2)
            rmu(3) = box(3) / oldbox(3)
            !write(6,*) "DEBUG: RMU: ",rmu(1),rmu(2),rmu(3)

            ! 3b- Update the box
            ! DAN ROE: Make sure fill_tranvec is needed.
            call redo_ucell(rmu)
            call fill_tranvec()
            !write(6,*) "DEBUG: NEWabc: ",a,b,c

            ! 3c- Recompute grid sizes
            !   RCFFT needs an even dimension for x direction
            call compute_nfft((a + one)*half   ,nfft1)
            nfft1=nfft1*2
            call compute_nfft(b,nfft2)
            call compute_nfft(c,nfft3)
         endif ! box sizes have changed

         ! DAN ROE: More Debug
         !write(6,'(/a)') 'Ewald parameters:'
         !write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
         !      ', ew_type =',ew_type,', nbflag  =',nbflag, &
         !      ', use_pme =',use_pme
         !write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
         !      ', eedmeth =',eedmeth,', netfrc  =',netfrc
         !write(6, 9002) a, b, c
         !write(6, 9003) alpha, beta, gamma
         !write(6, 9004) nfft1, nfft2, nfft3
         !write(6, 9006) cutoffnb, dsum_tol
         !write(6, 9007) ew_coeff
         !write(6, 9005) order
         !9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
         !9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
         !9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
         !9005 format (5x,'Interpolation order =',i5)
         !9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
         !9007 format (5x,'Ewald Coefficient =',f9.5)
      endif ! ntb>0

      110 format(10f8.3)

      member=member+1

      write (6,'(a,i6)') 'minimizing coord set #',member

      call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
            ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(ibellygp), &
            x(l95),ene,carrms,qsetup)

      write (6,364) ene(23),carrms
      364 format ('minimization completed, ENE=',1x,e14.8, &
            1x,'RMS=',1x,e12.6)

      if (master .and. itgtmd == 1) then
         write (6,'(a,f8.3)') "Final RMSD from reference: ",rmsdvalue
      end if

      ! write the frame to the mdcrd file if the user has set ntwx=1
      ! don't worry about imaging etc since we write the same way it came in

      if (master .and. ntwx >= 1) then
         call corpac(x(lcrd),1,natom*3,MDCRD_UNIT,loutfm)
         if(ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
      !elseif (master) then
      !   write (6,*) "Not writing coordinates to mdcrd due to NTWX value"
      endif

      !       ---loop for next coordinate set

   end do

   !     ---end of trajectory file

   1000 write (6,*) "Trajectory file ended"
   ok=.true.
   return

   1010 write (6,*) "Error in trajectory file"
   ok=.false.
   return

end subroutine trajene

