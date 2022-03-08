! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Identify a titrating proton in each residue to use as center for updatepairs
subroutine cnstphinit(stateinf,trescnt,chrgdat,x,hindex,tpair)
   implicit none
#  include "dynph.h"
   type(const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent(in) :: trescnt
   _REAL_, intent(in) ::  chrgdat(0:*), x(*)
   integer, intent(out) :: hindex(0:*), tpair(0:*)
   integer :: res, atom, state
   ! Locate the (last) titrating proton of each residue
   do res = 0, trescnt - 1
      hindex(res) = stateinf(res)%first_atom !Guard against poorly defined titrating groups
      do atom = 0, stateinf(res)%num_atoms - 1
         do state = 0, stateinf(res)%num_states - 1
            if (0.d0 == chrgdat(stateinf(res)%first_charge &
                             + state * stateinf(res)%num_atoms &
                             + atom)) then
               hindex(res) = stateinf(res)%first_atom + atom
            end if
         end do
      end do
   end do
   call cnstphupdatepairs(x,trescnt,hindex,tpair)
end subroutine cnstphinit
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Locate neighboring titrating residues for multi-site MC moves
!+ Called from cntsphinit and on steps where restarts are written
subroutine cnstphupdatepairs(x,trescnt,hindex,tpair)
   implicit none
   _REAL_, intent(in) :: x(*)
   integer, intent(in) :: trescnt, hindex(0:*)
   integer, intent(out) :: tpair(0:*)
   ! tpair is the array holding information about interacting
   ! (neighboring) pairs of titrating residues. It has two sections,
   ! a header of length trescnt+1 and a dynamic list section.
   ! The header gives indices into the list. Each element of the
   ! header gives the index of the first element in the pair list for
   ! that residue
   ! e.g. tpair(3) gives the index of the first "neighbor" of res 3
   ! The index of the last "neighbor" can be determined by looking at
   ! tpair(i+1) (see cnstphbeginstep for example. tpair(trescnt)
   ! contains a tail index, so this trick works for the last residue
   ! as well.
   ! The elements of the list section hold the (titrating) residue
   ! numbers of the neighbors
   integer :: i, j, pi,  ti, tj
   _REAL_ :: xij, yij, zij

   pi = trescnt + 1             !Start dynamic list after header

   do i = 0, trescnt - 1
      tpair(i) = pi             !Start index into list for current residue
      ti = 3 * hindex(i)        !coordinate index for titrating proton of residue i
      do j = 0, trescnt - 1      !Optimization possible, but probably not worth it
         if (i == j) cycle
         tj = 3 * hindex(j)
         xij = x(ti-2) - x(tj-2)
         yij = x(ti-1) - x(tj-1)
         zij = x(ti)   - x(tj)
         ! If titrating hydrogens are within 5.0 angstrom cut-off
         if (5.d0 > sqrt(xij*xij + yij*yij + zij*zij)) then
            if (pi < TITR_STATES_C * 4) then 
               tpair(pi) = j
               pi = pi + 1
            else
               write(6,*) "Non-fatal constant pH pair list overflow. Increase TITR_RES_C in dynph.h; recompile"
            end if
         end if
      end do
   end do
   tpair(trescnt) = pi
 end subroutine cnstphupdatepairs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnstphbeginstep(stateinf,resstate,trescnt, &
      statene,chrgdat,dcharge, iselres, iselstat, gen, tpair)
   implicit none
#  include "dynph.h"
#  include "random.h"
#  include "md.h"   
   type(const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent(in) :: resstate(0:*), trescnt, tpair(0:*)
   integer, intent(out) :: iselres(0:*), iselstat(0:*)
   _REAL_, intent(in) ::  statene(0:*), chrgdat(0:*)
   _REAL_, intent(out) :: dcharge(1:*)
   type (rand_gen_state), intent(inout) :: gen
   integer :: i,j, neighborcount
   _REAL_ randval

   call amrand_gen(gen, randval)
   if (randval < 0.25d0) then   !25% chance of multi-site move
      iselres(0) = 2
   else
      iselres(0) = 1
   end if
   
   call amrand_gen(gen,randval)         !Select residue
   iselres(1) = int((randval*0.9999999d0)*trescnt)
   call amrand_gen(gen,randval)         !Select new state (always different from current)
   iselstat(1) = int((randval*0.9999999d0) &
         *(stateinf(iselres(1))%num_states - 1))
   if (iselstat(1) >= resstate(iselres(1))) iselstat(1) = iselstat(1) + 1
   call cnstphupdatechrg(dcharge, iselres(1), iselstat(1), chrgdat, stateinf)
   
   if (iselres(0) == 2) then    !Do two-site MC trial
      neighborcount = tpair(iselres(1)+1)-tpair(iselres(1))
      if (neighborcount > 0) then
         call amrand_gen(gen,randval)
         iselres(2) = tpair(tpair(iselres(1)) & !get neighbor res number from base index
               + int(randval*0.9999999d0*neighborcount)) !plus neighbor index
         call amrand_gen(gen,randval)         !Select new state (always different from current)
         iselstat(2) = int((randval*0.9999999d0) &
               *(stateinf(iselres(2))%num_states - 1))
         if (iselstat(2) >= resstate(iselres(2))) iselstat(2) = iselstat(2) + 1
         call cnstphupdatechrg(dcharge, iselres(2), iselstat(2), chrgdat, stateinf)
      else
         iselres(0) = 1         !No neighbors, so no multi-site MC
      end if
   end if
end subroutine cnstphbeginstep


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Evaluate transition energy for proposed protonation state, accept or revert
subroutine cnstphendstep (stateinf,resstate,protcnt, &
      trescnt,statene,chrgdat,dcharge,charge, dvdl, &
      iselres,iselstat, gen)
   implicit none
#  include "dynph.h"
#  include "random.h"
#  include "md.h"   
   type(const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent(inout) :: resstate(0:*)
   _REAL_, intent(inout) :: dvdl
   integer, intent(in) :: protcnt(0:*)
   integer, intent(in) :: trescnt, iselres(0:*),iselstat(0:*)
   _REAL_, intent(in) :: statene(0:*), chrgdat(0:*)
   _REAL_, intent(out) :: charge(1:*),dcharge(1:*)
   type (rand_gen_state), intent(inout) :: gen
   _REAL_ :: randval, deltae
   integer ::  i, statebase
   
   deltae = 0d0
   do i=1,iselres(0)
      statebase = stateinf(iselres(i))%first_state
      !     delta energy = E(proposed state) - E(current state)
      deltae = deltae + statene(iselstat(i)+statebase) &
            -statene(resstate(iselres(i))+statebase)
      !     Correct for pH (delta protons * pH * 2.303RT)
      deltae = deltae -(protcnt(iselstat(i)+statebase) &
            -protcnt(resstate(iselres(i))+statebase)) &
            *solvph*2.303d0*2.d-3*temp0
   end do
   call amrand_gen(gen,randval)

   dvdl = dvdl - deltae         !Adjust transition energy for non-elec factors

   if ((dvdl < 0) .or. (randval <= exp(-dvdl/(2.d-3*temp0)))) then
      !     Make new charges permanent
      do i=1,iselres(0)         !Generalized to n-site MC move
         call cnstphupdatechrg(charge, iselres(i), iselstat(i), chrgdat, stateinf)
         resstate(iselres(i)) = iselstat(i)
      end do
   else
      !     Revert dcharge to old charges
      do i=1,iselres(0)
         call cnstphupdatechrg(dcharge, iselres(i), resstate(iselres(i)), chrgdat, stateinf)
      end do
   end if
end subroutine cnstphendstep 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Update charges to reflect new state
subroutine cnstphupdatechrg(charge,iselres,inewstate,chrgdat,stateinf)
   implicit none
#  include "dynph.h"
   _REAL_, intent(out) :: charge(1:*)
   integer, intent(in) :: iselres,inewstate
   _REAL_, intent(in) :: chrgdat(0:*)
   type(const_ph_info), intent(in) :: stateinf(0:*)

   integer i
   do i = 0, stateinf(iselres)%num_atoms-1
      charge(i+stateinf(iselres)%first_atom) &
            =chrgdat(stateinf(iselres)%first_charge &
            + inewstate * stateinf(iselres)%num_atoms &
            + i)
   end do
end subroutine cnstphupdatechrg

subroutine cnstphwriterestart(stateinf, resstate, protcnt, trescnt, statene, inchrgdat)
   use constants, only : AMBER_ELECTROSTATIC
   implicit none
#  include "dynph.h"
#  include "files.h"
   
   type (const_ph_info), intent (in) :: stateinf(0:TITR_RES_C-1)
   integer, intent(in) :: resstate(0:TITR_RES_C-1), protcnt(0:TITR_STATES_C-1)
   integer, intent(in) :: trescnt
   _REAL_, intent(in) ::  statene(0:TITR_STATES_C-1), inchrgdat(0:ATOM_CHRG_C-1)
   _REAL_ :: chrgdat(0:ATOM_CHRG_C-1)
   integer :: iatom
   logical, save :: first = .true.
   character(len=40) :: resname(0:TITR_RES_C-1)
   character(len=7) :: stat
   
   common /cnstphresname/ resname
  
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene,trescnt, resname

   do iatom=0, ATOM_CHRG_C-1
      chrgdat(iatom) = inchrgdat(iatom) / AMBER_ELECTROSTATIC
   end do
   
   if (first) then              ! Need DELIM, so can't use amopen
      if (owrite == 'N') then
         stat = 'NEW'
      else if (owrite == 'O') then
         stat = 'OLD'
      else if (owrite == 'R') then
         stat = 'REPLACE'
      else if (owrite == 'U') then
         stat = 'UNKNOWN'
      end if
      open(unit=CNSTPH_UNIT, file=cprestrt, status=stat, form='FORMATTED', DELIM='APOSTROPHE')
      first = .false.
   else
      open(unit=CNSTPH_UNIT, file=cprestrt, status='OLD', form='FORMATTED', DELIM='APOSTROPHE')
   end if
   write(CNSTPH_UNIT, nml=cnstph)
   close(CNSTPH_UNIT)

   call amflsh(CPOUT_UNIT)       ! Make sure all cpout data up to this restart point is on disk
end subroutine cnstphwriterestart

subroutine cnstphwrite(resstate,iselres,trescnt)
#  include "dynph.h"
#  include "md.h"
#  include "extra.h"
#  include "files.h"
  
   integer, intent(in) :: resstate(0:TITR_RES_C-1)
   integer, intent(in) :: trescnt, iselres(0:*)
   
   logical :: full
   logical, save :: first = .true.
   integer :: i
   
   if (.not. master) then 
      return
   end if

   full = (first .or. irespa == nstlim .or. mod(irespa,ntwr)==0)
   
   if (full) then
       write(CPOUT_UNIT, '(a,f8.5)') 'Solvent pH: ',solvph
      write(CPOUT_UNIT, '(a,i8)') 'Monte Carlo step size: ',ntcnstph
      write(CPOUT_UNIT, '(a,i8)') 'Time step: ',irespa
      write(CPOUT_UNIT, '(a,f10.3)') 'Time: ',t
      do i=0,trescnt-1
         write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ',i,' State: ',resstate(i)
      end do
   else
      do i=1,iselres(0)
         write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ',iselres(i),' State: ',resstate(iselres(i))
      end do
   end if
   write(CPOUT_UNIT, '()')
   first = .false.
end subroutine cnstphwrite

  
