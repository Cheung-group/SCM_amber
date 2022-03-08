! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "assert.h"
!QMMM Charge adjustment routines: Written by Ross Walker (TSRI 2005)
!Zeros out charges on QM atoms
subroutine qm_zero_charges(charges)

!Zero's the Q charges in array charges - before it does this though
!it copies the current QM atom charges to the array qmmm_struct%qm_resp_charges

  use qmmm_module, only : qmmm_nml,qmmm_struct
  use constants, only : zero
  implicit none

!Passed in
  _REAL_ charges(*)

!local
  integer i
  qmmm_struct%qm_resp_charge_sum = zero
  do i=1,qmmm_struct%nquant    !loop over all quantum atoms
    qmmm_struct%qm_resp_charges(i) = charges(qmmm_nml%iqmatoms(i))
    qmmm_struct%qm_resp_charge_sum = qmmm_struct%qm_resp_charge_sum + &
                                     charges(qmmm_nml%iqmatoms(i))
    charges(qmmm_nml%iqmatoms(i)) = zero
  enddo

end subroutine qm_zero_charges

subroutine qm_zero_mm_link_pair_main_chg(nlink,link_pairs,charges)
!This subroutine is used to zero the charge on any MM atom that is part
!of a QM-MM link pair. The charge is zeroed in the main charge array. This has the
!effect of removing the interaction between the MM link pair atom and all the Real
!QM atoms AND MM-MMlinkpair interactions.

  use qmmm_module, only : qmmm_struct 
  implicit none

!Passed in
  integer, intent(in) :: nlink
  integer, intent(in) :: link_pairs(2,nlink)
  _REAL_, intent(inout) :: charges(*)

!Local
  integer i,y,ier

  if (qmmm_struct%zero_link_charges_first_call) then
! Allocate and store the original charges
    allocate(qmmm_struct%mm_link_pair_resp_charges(nlink),stat=ier)
    REQUIRE(ier == 0)
    qmmm_struct%zero_link_charges_first_call = .false.
  end if

  do i = 1, nlink
    y = link_pairs(1,i) !The  id of the MM link pair atom
    qmmm_struct%mm_link_pair_resp_charges(i)=charges(y)
    charges(y) = 0.0d0 !Zero the charge on this MM atom.
  end do

  return

end subroutine qm_zero_mm_link_pair_main_chg

subroutine qmmm_adjust_q(adjust_q, natom, nquant, nquant_nlink, nlink, charge, iqmatoms, &
                         qmcharge, atom_mask, mm_link_mask, master, coords)
!
!This routine must be called BEFORE zeroing the QM and MML charges.
!
!This subroutine is responsible for adjusting the resp charge of certain MM atoms in order
!to ensure that the total system charge is conserved. This is particularly important in QMMM
!GB simulations.
!
!When a QM region is selected RESP charges that are the sum of all QM + all MM link pair atoms
!are deleted from the system and replaced with the value of QM charge. Thus the missing charge
!is the difference between these values. E.g. if QM+MM link summed to +0.35 electrons and the
!QM charge was zero then +0.35 would need to be added to the remaining MM region somehow.

  use constants, only : zero, INV_AMBER_ELECTROSTATIC, AMBER_ELECTROSTATIC
  implicit none

!Passed in
  integer, intent(in) :: adjust_q,natom, nquant, nquant_nlink, nlink
  integer, intent(in) :: iqmatoms(nquant_nlink)
  integer, intent(in) :: qmcharge
  _REAL_, intent(inout) :: charge(natom), coords(3,natom)
  logical, intent(in) :: atom_mask(natom), mm_link_mask(natom)
  logical, intent(in) :: master

!Local
  integer :: i,j, closest_id, link_no
   _REAL_ :: q_correction_sum, q_correction, correction, final_q_sum
   _REAL_ :: smallest_dist, vec(3), link_coord(3), r2

  !If this is a pure QM calculation then there is nothing to do here.
  if (nquant == natom) return

  !q_correction = Sum(QM+MML) - QM_Charge

  q_correction_sum = zero

  !iqmatoms contains the MM number of the QM atom or the MML number of the link atom.
  do i = 1, nquant_nlink
    q_correction_sum = q_correction_sum + charge(iqmatoms(i))
  end do
  q_correction = q_correction_sum - dble(qmcharge)*AMBER_ELECTROSTATIC

  if (adjust_q == 1) then
    !if adjust_q == 1 then we divide this up between the atom nearest to each of the nlink
    !MM link pair atoms.
    !Sanity check - need natom-nquant_nlink >= nlink and nlink>0
    if ((natom - nquant_nlink)<nlink .or. nlink<1) &
        call sander_bomb('qmmm_adjust_q','QMMM: Error adjust_q=1 requires natom-(nquant+nlink)>=nlink','and nlink>0 abort')

    if (master) then
       write(6,'("QMMM: ADJUSTING CHARGES")')
       write(6,'("QMMM: ----------------------------------------------------------------------")')
       write(6,'("QMMM: adjust_q = 1")')
       write(6,'("QMMM: Adjusting the charge of closest nlink MM atoms to MM link pairs")')
       write(6,'("QMMM: to conserve total charge.")')
       write(6,'("QMMM: Atoms being adjusted = ")',ADVANCE='NO')
    end if
    correction = q_correction / dble(nlink)
    do i = 1, nlink
      !Go through each link atom in turn and find the MM atom that is closest to it.
      !This is an expensive procedure but since this is part of the initial setup it
      !doesn't matter too much.
      closest_id = -1
      smallest_dist = 1.0d30
      link_no = iqmatoms(nquant+i)
      link_coord(1:3) = coords(1:3,link_no)
      do j = 1, natom
        if ((.not. atom_mask(j)) .and. (.not. mm_link_mask(j))) then
          !It is a REAL (non link) MM atom.
          !Calculate the distance^2 and see if it is less than our smallest_dist
          vec(1:3) = coords(1:3,j) - link_coord(1:3)
          r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
          if (r2 <= smallest_dist) then
            smallest_dist=r2
            closest_id = j
          end if
        end if
      end do
      !Sanity check - if closest_id is zero we have a problem.
      REQUIRE(closest_id>0)
      if (master) write(6,'(i6)',ADVANCE='NO') closest_id
      charge(closest_id) = charge(closest_id) + correction
    end do
    if (master) then
      final_q_sum = qmcharge * AMBER_ELECTROSTATIC
      do i = 1, natom
        if ((.not. atom_mask(i)) .and. (.not. mm_link_mask(i))) final_q_sum = final_q_sum + charge(i)
      end do
      write(6,'(1x)')
      write(6,'("QMMM:                                  qm_charge = ",i4)') qmcharge
      write(6,'("QMMM:      QM atom RESP charge sum (inc MM link) = ",f8.3)') q_correction_sum*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM: Adjusting selected MM atom resp charges by = ",f8.3)') correction*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM:               Sum of MM + QM region is now = ",f8.3)') final_q_sum*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM: ----------------------------------------------------------------------")')
    end if
  else if (adjust_q == 2) then
    !if adjust_q == 2 then we simply divide the correction sum by the number of remaining
    !atoms and then add that value to each atom.
    !Sanity check - need natom-nquant_nlink >= 1
    if ((natom - nquant_nlink)<nlink) &
        call sander_bomb('qmmm_adjust_q','QMMM: Error adjust_q=2 requires natom-(nquant+nlink)>=1','abort')
    correction = q_correction / dble(natom-nquant_nlink) !Should never divide by zero here as code
                                                         !returns above if nquant==natom. - i.e. pure QM run.

    final_q_sum = qmcharge * AMBER_ELECTROSTATIC
    do i = 1, natom
      !Go through each atom adding correction to it's charge if it
      !is not a QM or MM Link atom.
      if ((.not. atom_mask(i)) .and. (.not. mm_link_mask(i))) then
        charge(i) = charge(i) + correction
        final_q_sum = final_q_sum+charge(i)
      end if
    end do

    if (master) then
      write(6,'("QMMM: ADJUSTING CHARGES")')
      write(6,'("QMMM: ----------------------------------------------------------------------")')
      write(6,'("QMMM: adjust_q = 2")')
      write(6,'("QMMM: Uniformally adjusting the charge of MM atoms to conserve total charge.")')
      write(6,'("QMMM:                             qm_charge = ",i4)') qmcharge
      write(6,'("QMMM: QM atom RESP charge sum (inc MM link) = ",f8.3)' ) q_correction_sum*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM: Adjusting each MM atom resp charge by = ",f8.3)') correction*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM:          Sum of MM + QM region is now = ",f8.3)') final_q_sum*INV_AMBER_ELECTROSTATIC
      write(6,'("QMMM: ----------------------------------------------------------------------")')
    end if

  end if


  return

end subroutine qmmm_adjust_q


