#include "copyright.h"
#include "dprec.h"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit the final minimization report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_RESULTS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! optionally do pseudocontact shift constraints,
! optionally decompose energies for mm_pbsa,
! and emit nmr related information (nmrptx).

subroutine report_min_results( nstep, gradient_rms, coordinates, &
      forces, energies, igraph, xx, ix, ih )

   use decomp, only : checkdec, printdec, printpdec
   implicit none

   integer, intent(in)             :: nstep
   _REAL_,  intent(in)             :: gradient_rms
   _REAL_,  intent(in)             :: coordinates(*)
   _REAL_,  intent(in)             :: forces(*)
   _REAL_,  intent(in)             :: energies(51)
   character(len=4), intent(in)    :: igraph(*)    ! atom name map
   _REAL_,  intent(inout)          :: xx(*)        ! real dynamic memory
   integer, intent(inout)          :: ix(*)        ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*)        ! hollerith dynamic memory

#  include "box.h"
#  include "extra.h"
#  include "files.h"
#  include "md.h"
#  include "memory.h"
#  include "nmr.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max

   if (master) then
      write(6, '(/ /20x,a,/)' ) 'FINAL RESULTS'
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (idecomp > 0) call checkdec(idecomp)
      if (idecomp == 1 .or. idecomp == 2) call printdec(ix)
      if (idecomp == 3 .or. idecomp == 4) call printpdec(ix)
      if (nmropt > 0) then
         if (iredir(7) /= 0) call pcshift(-1,coordinates,forces)
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
         call ndvptx(coordinates,forces,ih(m04),ih(m02),ix(i02),nres, &
               xx(l95),natom,xx(lwinv),ntb,xx(lnmr01),ix(inmr02),6)
      end if
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_results


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit a minimization progress report 
!-----------------------------------------------------------------------
!     --- REPORT_MIN_PROGRESS ---
!-----------------------------------------------------------------------
! Find the maximum component of the gradient (grdmax),
! emit this and other gradient details (printe),
! and emit nmr related information (nmrptx).

subroutine report_min_progress( nstep, gradient_rms, forces, energies, igraph )

   implicit none

   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: forces(*)
   _REAL_,  intent(in)    :: energies(51)
   character(len=4), intent(in)    :: igraph(*)    ! atom name map

#  include "extra.h"
#  include "files.h"
#  include "md.h"
#  include "nmr.h"

   character(len=4) :: atom_name_of_gmax
   integer          :: atom_number_of_gmax
   _REAL_           :: gradient_max

   if (master) then
      call grdmax( forces, gradient_max, atom_number_of_gmax )
      atom_name_of_gmax = igraph(atom_number_of_gmax)
      if (imin /= 5) rewind(MDINFO_UNIT)
      call printe( nstep, gradient_rms, gradient_max, energies, &
             atom_number_of_gmax, atom_name_of_gmax )
      if (nmropt > 0) then
         call nmrptx(6)
         call nmrptx(MDINFO_UNIT)
      end if
      if (imin /= 5) call amflsh(MDINFO_UNIT)
   end if

   return
end subroutine report_min_progress


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the maximum gradient component and the corresponding atom
!-----------------------------------------------------------------------
!     --- GRDMAX ---
!-----------------------------------------------------------------------

subroutine grdmax( gradient, gradient_max, atom_number_of_gmax )

   use constants, only : zero
   implicit none
   _REAL_,  intent(in)  :: gradient(*)
   !     This is actually two-dimensional (3,natoms), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   _REAL_,  intent(out) :: gradient_max
   integer, intent(out) :: atom_number_of_gmax

#  include "memory.h"

   integer :: i
   _REAL_  :: gi

   gradient_max = ZERO
   atom_number_of_gmax = 1
   do i = 1,3*natom
      gi = abs(gradient(i))
      if (gi > gradient_max) then
         gradient_max = gi
         atom_number_of_gmax = i
      end if
   end do
   atom_number_of_gmax = (atom_number_of_gmax - 1)/3 + 1

   return
end subroutine grdmax


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print status of a minimization calculation.
!-----------------------------------------------------------------------
!     --- PRINTE ---
!-----------------------------------------------------------------------
! Output the step number, the gradient rms, the max gradient component,
! and the atom label and atom number of the max gradient component.

subroutine printe( nstep, gradient_rms, gradient_max, energies, &
      atom_number_of_gmax, atom_name_of_gmax )
   
   use qmmm_module, only : qmmm_nml,qmmm_struct, PM3, AM1, MNDO, PDDGPM3, &
                           PDDGMNDO, PM3CARB1, DFTB, RM1
   use cns_xref

   implicit none
   integer, intent(in)          :: nstep
   _REAL_,  intent(in)          :: gradient_rms
   _REAL_,  intent(in)          :: gradient_max
   _REAL_,  intent(in)          :: energies(51)
   integer, intent(in)          :: atom_number_of_gmax
   character(len=4), intent(in) :: atom_name_of_gmax

#  include "md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"
#ifdef APBS
#  include "files.h"
#endif /* APBS */

   _REAL_  epot,enonb,enele,ehbond,ebond,eangle,edihed,enb14,eel14
   _REAL_  econst,epolar,aveper,aveind,avetot,esurf,edisp,diprms,dipiter
   _REAL_  dipole_temp,escf,dvdl

   !     ----- Extract Energies. -----

   epot   = energies(1)
   enonb  = energies(2)
   enele  = energies(3)
   ehbond = energies(4)
   ebond  = energies(5)
   eangle = energies(6)
   edihed = energies(7)
   enb14  = energies(8)
   eel14  = energies(9)
   econst = energies(10)
   epolar = energies(11)
   aveper = energies(12)
   aveind = energies(13)
   avetot = energies(14)
   esurf  = energies(15)
   dvdl   = energies(17)
   diprms = energies(22)
   dipiter = energies(23)
   dipole_temp = energies(24)
   escf = energies(25)
   edisp  = energies(28)

   write(6,9018)
   write(6,9028) nstep, epot, gradient_rms, gradient_max, &
         atom_name_of_gmax, atom_number_of_gmax
   write(6,9038) ebond,eangle,edihed
   if( igb == 0 ) then
      write(6,9048) enonb,enele,ehbond
   else if ( igb == 10 ) then
      write(6,9050) enonb,enele,ehbond
#ifdef APBS
   else if ( igb == 6 .and. mdin_apbs ) then
      write(6,9050) enonb,enele,ehbond
#endif /* APBS */
   else
      write(6,9049) enonb,enele,ehbond
   end if
   write(6,9058) enb14,eel14,econst
   if (qmmm_nml%ifqnt) then
     !write the SCF energy
     if (qmmm_nml%qmtheory == PM3) then
        write(6,9080) escf
     else if (qmmm_nml%qmtheory == AM1) then
        write(6,9081) escf
     else if (qmmm_nml%qmtheory == MNDO) then
        write(6,9082) escf
     else if (qmmm_nml%qmtheory == PDDGPM3) then
        write(6,9083) escf
     else if (qmmm_nml%qmtheory == PDDGMNDO) then
        write(6,9084) escf
     else if (qmmm_nml%qmtheory == PM3CARB1) then
        write(6,9085) escf
     else if (qmmm_nml%qmtheory == DFTB) then
        write(6,9086) escf
     else if (qmmm_nml%qmtheory == RM1) then
        write(6,9087) escf
     else
        write(6,'(" ERROR - UNKNOWN QM THEORY")')
     end if
   end if
#ifdef PUPIL_SUPPORT
   write(6,9900) escf
#endif
   if( gbsa > 0 ) write(6,9077) esurf
   if (igb == 10) write(6,9074) esurf,edisp
   call write_cns_xref_min_energies( energies(47), energies(49), &
      energies(50), epot-econst-energies(47) )
#ifdef APBS
   if (igb == 6 .and. mdin_apbs ) write(6,9069) esurf
#endif /* APBS */
   if (epolar /= 0.0) write(6,9068) epolar
   if (econst /= 0.0) write(6,9078) epot-econst
   if ( dvdl /= 0.d0) write(6,9089) dvdl
   if (induced > 0.and.indmeth < 3) write(6,9090)diprms,dipiter
   if (induced > 0.and.indmeth == 3) write(6,9091)diprms, &
         dipole_temp
   if (itgtmd == 1) then
      write (6,'(a,f8.3)') "Current RMSD from reference: ",rmsdvalue
      write (6,'(a,f8.3)') "Current target RMSD:         ",tgtrmsd
   end if

   call amflsh(6)

   !     ----- SEND IT TO THE INFO FILE -----

   if (imin /= 5) then
      write(7,9018)
      write(7,9028) nstep, epot, gradient_rms, gradient_max, &
            atom_name_of_gmax, atom_number_of_gmax
      write(7,9038) ebond,eangle,edihed
      if( igb == 0 ) then
         write(7,9048) enonb,enele,ehbond
      else if ( igb == 10 ) then
         write(7,9050) enonb,enele,ehbond
#ifdef APBS
      else if ( igb == 6 .and. mdin_apbs ) then
         write(7,9050) enonb,enele,ehbond
#endif /* APBS */
      else
         write(7,9049) enonb,enele,ehbond
      end if
      write(7,9058) enb14,eel14,econst
      if (qmmm_nml%ifqnt) then
        !write the SCF energy
        if (qmmm_nml%qmtheory == PM3) then
           write(7,9080) escf
        else if (qmmm_nml%qmtheory == AM1) then
           write(7,9081) escf
        else if (qmmm_nml%qmtheory == MNDO) then
           write(7,9082) escf
        else if (qmmm_nml%qmtheory == PDDGPM3) then
           write(7,9083) escf
        else if (qmmm_nml%qmtheory == PDDGMNDO) then
           write(7,9084) escf
        else if (qmmm_nml%qmtheory == PM3CARB1) then
           write(7,9085) escf
        else if (qmmm_nml%qmtheory == DFTB) then
           write(7,9086) escf
        else if (qmmm_nml%qmtheory == RM1) then
           write(7,9087) escf
        else
           write(7,'(" ERROR - UNKNOWN QM THEORY")')
        end if
      end if
#ifdef PUPIL_SUPPORT
      write(7,9900) escf
#endif
      if( gbsa > 0 ) write(7,9077) esurf
      if ( igb == 10 ) write(7,9074) esurf,edisp
#ifdef APBS
      if (igb == 6 .and. mdin_apbs ) write(7,9069) esurf
#endif /* APBS */
      if (epolar /= 0.0) write(7,9068) epolar
      if (econst /= 0.0) write(7,9078) epot-econst
      if ( dvdl /= 0.d0) write(7,9089) dvdl
      if (induced > 0.and.indmeth < 3) write(7,9090)diprms,dipiter
      if (induced > 0.and.indmeth == 3) write(7,9091)diprms, &
            dipole_temp
   end if

   9018 format (/ /,3x,'NSTEP',7x,'ENERGY',10x,'RMS',12x,'GMAX',9x, &
         'NAME',4x,'NUMBER')
   9028 format(1x,i6,2x,3(2x,1pe13.4),5x,a4,2x,i7,/)
   9038 format (1x,'BOND    = ',f13.4,2x,'ANGLE   = ',f13.4,2x, &
         'DIHED      = ',f13.4)
   9048 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'HBOND      = ',f13.4)
   9049 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EGB        = ',f13.4)
   9050 format (1x,'VDWAALS = ',f13.4,2x,'EEL     = ',f13.4,2x, &
         'EPB        = ',f13.4)
   9058 format (1x,'1-4 VDW = ',f13.4,2x,'1-4 EEL = ',f13.4,2x, &
         'RESTRAINT  = ',f13.4)
   9068 format (1x,'EPOLAR  = ',f13.4)
#ifdef APBS
   9069 format (1x,'ENPOLAR = ',f13.4)
#endif /* APBS */
   9074 format (1x,'ECAVITY = ',f13.4,2x,'EDISPER = ',f13.4)
   9077 format (1x,'ESURF   = ',f13.4)
   9078 format (1x,'EAMBER  = ',f13.4)
   9080 format (1x,'PM3ESCF =',f14.4)
   9081 format (1x,'AM1ESCF =',f14.4)
   9082 format (1x,'MNDOESCF=',f14.4)
   9083 format (1x,'PDDGPM3-ESCF=',f14.4)
   9084 format (1x,'PDDGMNDO-ESCF=',f14.4)
   9085 format (1x,'PM3CARB1-ESCF=',f14.4)
   9086 format (1x,'DFTBESCF=',f14.4)
   9087 format (1x,'RM1ESCF =',f14.4)
   9089 format (1x,'DV/DL  = ',f14.4)
   9090 format(1x,'Dipole convergence: rms = ',e10.3,' iters = ',f6.2)
   9091 format(1x,'Dipole convergence: rms = ',e10.3, &
         ' temperature = ',f6.2)
#ifdef PUPIL_SUPPORT
   9900 format (1x,'PUPESCF =',f14.4)
#endif
   return
end subroutine printe
