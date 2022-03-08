! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "assert.h"
subroutine qm2_calc_rij_and_eqns(coords, nquant_nlink, crdsmm, natom, npairs)

!-----------------------------------------------------------------
!Written by Ross Walker (TSRI, 2005)
!
!This routine should be called on each call to QM_MM. It's purpose
!is to calculate RIJ for each QM-QM pair. It stores this is in
!the qm2_rij_eqns structure. It also calculated a number of RIJ
!related equations that involve exponentials, sqrts etc. In this
!way they only need to be calculated once rather than several times
!during the energy calculation and then in the derivative code.
!
!-----------------------------------------------------------------

  use qmmm_module, only : qmmm_nml, qmmm_struct, qm2_rij_eqns, qm2_params, &
                          qmmm_mpi,alph_mm, PDDGPM3, PDDGMNDO
  use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, one, two
  implicit none

! Passed in
   integer, intent(in) :: nquant_nlink, natom, npairs
   _REAL_, intent(in) :: coords(3,nquant_nlink), crdsmm(*)
                         !crdsmm array is laid out as x,y,z,chg,x,y,z,chg...
! Local
   integer :: i,j, loop_count, four_npairs
   logical :: heavy_atom
   _REAL_ :: r2, rr2, rr, vec(3), onerij, rij, qmi_alpa, qmj_alpa
   _REAL_ :: qmi_oneBDD1, qmi_oneBDD2, qmi_oneBDD3
   _REAL_ :: qmj_oneBDD1, pddge1i, pddge2i
   _REAL_ :: ijBDD1, qmi_DD, qmi_QQ, qmi_QQ2
   _REAL_ :: RRADD, RRMDD, RRAQQ, RRMQQ
   _REAL_ :: RIJ_temp, qmx(3)
   integer :: ier=0
   integer :: num_per_thread, jstart, jend

#include "qm2_array_locations.h"
 
!For ew_coeff
#include "ew_pme_recip.h"

   if (qmmm_nml%qmqmrij_incore) then
     !Memory allocation if needed. This array is a static size since the 
     !number of QM-QM interactions does not change during a simulation.
     !Thus only need to do the allocation on the first call.
     if (qmmm_struct%qm2_calc_rij_eqns_first_call) then
       call qm2_allocate_qm2_qmqm_rij_eqns(qmmm_struct%PDDG_IN_USE)
     end if
     loop_count = 0
#ifdef MPI
     do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
        jstart =  qmmm_mpi%nquant_nlink_jrange(1,i)
        jend = qmmm_mpi%nquant_nlink_jrange(2,i)
#else
     do I=2,qmmm_struct%nquant_nlink
        jstart = 1
        jend = i-1
#endif
        qmi_alpa=qm2_params%cc_exp_params(i)
        qmi_oneBDD1=qm2_params%multip_2c_elec_params(3,i)
        qmx(1:3) = coords(1:3,i)
        do j=jstart,jend
           qmj_alpa=qm2_params%cc_exp_params(j)
           qmj_oneBDD1=qm2_params%multip_2c_elec_params(3,j)
           ijBDD1=qmi_oneBDD1+qmj_oneBDD1
           ijBDD1=ijBDD1*ijBDD1
           loop_count=loop_count+1
           vec(1:3)=qmx(1:3)-coords(1:3,j)
           r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
           rr2=r2*A2_TO_BOHRS2
           qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS2,loop_count) = rr2
           onerij = one/sqrt(r2)
           qm2_rij_eqns%qmqmrijdata(QMQMONERIJ,loop_count)=onerij
           rij = r2*onerij !one/onerij
           qm2_rij_eqns%qmqmrijdata(QMQMRIJ,loop_count)=rij
           qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS,loop_count) = rij*A_TO_BOHRS
           qm2_rij_eqns%qmqmrijdata(QMQMEXP1I,loop_count)=EXP(-qmi_alpa*rij)
           qm2_rij_eqns%qmqmrijdata(QMQMEXP1J,loop_count)=EXP(-qmj_alpa*rij)
           qm2_rij_eqns%qmqmrijdata(QMQMSQRTAEE,loop_count)=one/sqrt(RR2+ijBDD1)
        end do
     end do
     !Fill the PDDG data if PDDG-PM3 or PDDG-MNDO is requested.
     if (qmmm_struct%PDDG_IN_USE) then
        loop_count = 0
#ifdef MPI
        do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
           jstart =  qmmm_mpi%nquant_nlink_jrange(1,i)
           jend = qmmm_mpi%nquant_nlink_jrange(2,i)
#else
        do I=2,qmmm_struct%nquant_nlink
           jstart = 1
           jend = i-1
#endif
          pddge1i = qm2_params%pddge1(i)
          pddge2i = qm2_params%pddge2(i)
          do j=jstart,jend
            loop_count=loop_count+1
            RIJ_temp=qm2_rij_eqns%qmqmrijdata(QMQMRIJ,loop_count)
            qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP1,loop_count) = &
                         EXP(-10.0D0 * (RIJ_temp - pddge1i - qm2_params%pddge1(j))**2)
            qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP2,loop_count) = &
                         EXP(-10.0D0 * (RIJ_temp - pddge1i - qm2_params%pddge2(j))**2)
            qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP3,loop_count) = &
                         EXP(-10.0D0 * (RIJ_temp - pddge2i - qm2_params%pddge1(j))**2)
            qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP4,loop_count) = &
                         EXP(-10.0D0 * (RIJ_temp - pddge2i - qm2_params%pddge2(j))**2)
          end do
       end do
     end if
   end if
   !Only full QM-MM interaction (multipole) benefits from qmmmrij_incore
   if (qmmm_nml%qmmmrij_incore) then
     !Allocate memory here for the qmmmrijdata array
     !If this is the first call it will be allocated as either
     !min(npairs*qmmm_struct%nquant+npairs,qmmm_struct%nquant * (natom - qmmm_struct%nquant))
     !Since the number of pairs can change on each call we should check each time to see
     !if we need to re-allocate.
     num_per_thread = qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1
     if (qmmm_struct%qm2_calc_rij_eqns_first_call) then
       call qm2_allocate_qm2_qmmm_rij_eqns(natom,npairs)
       qmmm_struct%qm2_calc_rij_eqns_first_call = .false.
     else if (min(npairs*(num_per_thread)+npairs, &
              (num_per_thread) * (natom - qmmm_struct%nquant)) > qm2_rij_eqns%qmmmrij_allocated) then
       !We need to de-allocate the array and then call the allocation routine again so it gets
       !allocated larger.
       deallocate(qm2_rij_eqns%qmmmrijdata,stat=ier)
       REQUIRE(ier==0)
       call qm2_allocate_qm2_qmmm_rij_eqns(natom,npairs)
     end if

     loop_count = 0
     four_npairs=ishft(npairs,2) !*4

     do i=qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end 
       heavy_atom = (qm2_params%natomic_orbs(i) > 1)
       qmi_alpa=qm2_params%cc_exp_params(i)
       qmi_DD=qm2_params%multip_2c_elec_params(1,i)
       qmi_QQ=qm2_params%multip_2c_elec_params(2,i)*two
       qmi_oneBDD1=qm2_params%multip_2c_elec_params(3,i)**2
       qmi_oneBDD2=qm2_params%multip_2c_elec_params(4,i)**2
       qmi_oneBDD3=qm2_params%multip_2c_elec_params(5,i)**2
       qmi_QQ2=qmi_QQ*qmi_QQ+qmi_oneBDD3
       qmx(1:3)=coords(1:3,i)
       !RCW: We have duplication of code below but avoids us having to have the
       !if statement inside the inner loop
       if (heavy_atom) then
         do j=1,four_npairs,4
           loop_count=loop_count+1
           vec(1)=qmx(1)-crdsmm(j)
           vec(2)=qmx(2)-crdsmm(j+1)
           vec(3)=qmx(3)-crdsmm(j+2)
           r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
           rr2 = r2*A2_TO_BOHRS2
           onerij=one/sqrt(r2)
           qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)=onerij
           !rij=one/onerij
           rij=r2*onerij
           qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)=rij
           rr=rij*A_TO_BOHRS
           qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count) = exp(-qmi_alpa*rij)
           qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count) = exp(-ALPH_MM*rij)
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count) = one/sqrt(RR2+qmi_oneBDD1)
           !Heavy atom specific stuff
             RRADD = RR+qmi_DD
             RRADD = RRADD*RRADD+qmi_oneBDD2
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count) = one/sqrt(RRADD)
             RRMDD = RR-qmi_DD
             RRMDD = RRMDD*RRMDD+qmi_oneBDD2
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count) = one/sqrt(RRMDD)
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count) = one/SQRT(RR2+qmi_oneBDD3)
             RRAQQ=RR+qmi_QQ
             RRAQQ=RRAQQ*RRAQQ+qmi_oneBDD3
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count) = one/SQRT(RRAQQ)
             RRMQQ=RR-qmi_QQ
             RRMQQ=RRMQQ*RRMQQ+qmi_oneBDD3
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE,loop_count) = one/SQRT(RRMQQ)
             qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count) = one/SQRT(RR2+qmi_QQ2)
           !End Heavy atom specific stuff
         end do
       else !(heavy_atom)
         do j=1,four_npairs,4
           loop_count=loop_count+1
           vec(1)=qmx(1)-crdsmm(j)
           vec(2)=qmx(2)-crdsmm(j+1)
           vec(3)=qmx(3)-crdsmm(j+2)
           r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
           rr2 = r2*A2_TO_BOHRS2
           onerij=one/sqrt(r2)
           qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)=onerij
           rij=r2*onerij !one/onerij
           qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)=rij
           qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count) = exp(-qmi_alpa*rij)
           qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count) = exp(-ALPH_MM*rij)
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count) = one/sqrt(RR2+qmi_oneBDD1)
         end do
       end if !(heavy_atom)
     end do
   end if
   return
end subroutine qm2_calc_rij_and_eqns

!-----------------------------------------------------------------
!Written by Ross Walker (TSRI, 2005)
!Allocates memory for qm2_rij_eqns structure.
!-----------------------------------------------------------------

subroutine qm2_allocate_qm2_qmqm_rij_eqns(PDDG_IN_USE)

  use qmmm_module, only : qmmm_nml, qm2_rij_eqns, qmmm_struct, qmmm_mpi, &
                          PDDGPM3, PDDGMNDO
  implicit none

!Passed in
  logical, intent(in) :: PDDG_IN_USE

!Local
  integer :: ier=0
  integer :: array_size

  if (qmmm_nml%qmqmrij_incore) then
     !In parallel it only needs to be
     !qmmm_mpi%nquant_nlink_loop_extent_end-qmmm_mpi%nquant_nlink_loop_extent_begin+1
     !long.
#ifdef MPI
    array_size=qmmm_mpi%nquant_nlink_loop_extent_end-qmmm_mpi%nquant_nlink_loop_extent_begin+1
#else
    array_size = qmmm_struct%nquant_nlink * (qmmm_struct%nquant_nlink-1)/2
#endif
    if (PDDG_IN_USE) then
      allocate(qm2_rij_eqns%qmqmrijdata(QMQMNORIJ+QMQMNOPDDG,array_size),stat=ier)
      REQUIRE(ier==0)
    else
      allocate(qm2_rij_eqns%qmqmrijdata(QMQMNORIJ,array_size),stat=ier)
      REQUIRE(ier==0)
    end if
  end if
  return
end subroutine qm2_allocate_qm2_qmqm_rij_eqns

subroutine qm2_allocate_qm2_qmmm_rij_eqns(natom,npairs)

   use qmmm_module, only : qmmm_nml,qm2_rij_eqns, qmmm_struct, qmmm_mpi
   implicit none
!Passed in
  integer, intent(in) :: natom, npairs

!Local
  integer :: ier=0
  integer array_size, num_per_thread

  !We will allocate the qmmmrijdata array enough to hold the npairs+a bit
  !or qmmm_struct%nquant * (natom - qmmm_struct%nquant)
  !whichever is smaller. We then store the amount we allocated
  !in qm2_rij_eqns%qmmmrij_allocated. In this way when the number
  !of pairs changes this value can be checked and the array reallocated
  !large if necessary.

  !In parallel it only needs to be (qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)*qmmm_struct%qm_mm_pairs
  !in size

  if (qmmm_nml%qmmmrij_incore) then
     num_per_thread = qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1
     array_size = npairs*(num_per_thread)+npairs
     array_size = min(array_size,(num_per_thread) &
                               * (natom - qmmm_struct%nquant))
     allocate(qm2_rij_eqns%qmmmrijdata(QMMMNORIJ,array_size),stat=ier)
     REQUIRE(ier==0)
     qm2_rij_eqns%qmmmrij_allocated=array_size
     if (qmmm_nml%verbosity>1 .and. qmmm_mpi%commqmmm_master) then
       write(6,*) 'QMMM: Allocating qmmmrijdata array as ',QMMMNORIJ,' x ',array_size
       write(6,*) 'QMMM: to hold min npairs of ',npairs,' and nquant per thread of ',num_per_thread
     end if
  end if
  return
end subroutine qm2_allocate_qm2_qmmm_rij_eqns

