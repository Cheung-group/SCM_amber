! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "def_time.h"
!-----------------------------------
!Principal code for calculating
!QM potential for QMMM simulations.
!Principal Authors of current code:
!           Ross Walker
!           Mike Crowley
!
! Please send all comments or
! queries to: ross@rosswalker.co.uk
!-----------------------------------

subroutine qm_mm(coords,natom,scaled_mm_charges,f,escf, &
                 periodic,born_radii,one_born_radii, &
                 intdiel,extdiel,Arad,mmcut2, scf_mchg,atom_type) 
!
!     Argument list variables:
!
!     coords(natom*3)                 - Cartesian coordinates for all atoms.
!                                       Amber array
!     natom                           - Total number of REAL atoms.
!     qmmm_struct%nquant              - Number of REAL quantum atoms as specified in mdin.
!     iqmatoms(qmmm_struct%nquant)
!                                     - Atom numbers for quantum atoms link atoms given values of -1
!     scaled_mm_charges(natom)          - Atomic charges for mm atoms. (Scaled to elec units)
!     qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant) - Atomic numbers for qm atoms.
!     qmmm_struct%nlink               - Number of link atoms.
!     f((natom)*3)                    - Atomic forces.
!     qmmm_nml%qmcut               - cutoff in angstroms.
!     qmmm_nml%qmcut2              - cutoff^2 in angstroms^2.
!     escf                            - Heat of formation from QM.
!     qm_mm_pair_list(*)              - Atom-based list for quantum atoms:
!                                       i1,i2,i3. Each QM atom shares same list.
!     qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink)  
!                                     - Cartesian coordinates of quantum atoms.
!                                       (Extracted from coords by qm_extract_coords)

!     Locally defined arrays:
!     dxyzqm(3,qmmm_struct%nquant+qmmm_struct%nlink)     
!                                - Quantum mechanical derivatives from qm-mm
!                                                        interactions.
!     dxyzcl(3,natom)          - Classical derivatives from qm-mm interaction.
!     qm_xcrd(4,natom)         - imaged coords array from amber's array - note is 
!                                ordered in the same order as the pair_list. This
!                                is a DIFFERENT ORDER to the coords array.
!                                1->3 = coordinates, 4 = scaled_mm_charge
!    born_radii(1->natom)      - Effective GB radii - only used when doing qm with gb (and qm_gb==2)
!                                Calculated via an initial call to egb.
!    one_born_radii(1->natom)  - 1.0d0/born_radii(i)
!    scf_mchg                  - nquant long, gets filled with the mulliken charges during scf.
!    mmcut2 - cut^2 in angstroms^2 - cut from cntrl namelist

   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_rij_eqns, qmewald,  &
                           allocate_qmmm_pair_list, element_sym, qm_gb, DFTB, qmmm_mpi, qmmm_scratch
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
   use qm2_dftb_module, only: cm3
   use nblist,only: alpha,beta,gamma
  
   implicit none

#include "assert.h"

#ifdef MPI
#include "mpif.h"
#endif

! Passed in
   integer, intent(in) :: natom,periodic
   _REAL_ , intent(inout)  :: coords(natom*3) !Amber array - adjusted for link atoms
   _REAL_ , intent(in)  :: scaled_mm_charges(natom)
   _REAL_ , intent(out) :: f(natom*3)
   _REAL_ , intent(out) :: escf
   _REAL_ , intent(in) :: born_radii(natom), one_born_radii(natom)
   _REAL_ , intent(in) :: intdiel, extdiel, Arad, mmcut2
   _REAL_ , intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)
   character(len=4), intent(in) :: atom_type(*)

!Locals
   _REAL_ :: total_energy
   _REAL_ :: alpb_beta

   integer :: ier=0
   integer i, j, m, offset, qm_no

!! (GMS) For DFTB / CM3 charges
   _REAL_ :: cm3_chg, total_cm3_chg

!Locals for link atoms
   _REAL_ :: forcemod(3)
   integer :: lnk_no, mm_no

!=============================================================================
!                   START OF QMMM SETUP: allocate list memory
!=============================================================================

   call timer_start(TIME_QMMMSETUP)

   !Increment the counter of how many time qm_mm routine has been called.
   qmmm_struct%num_qmmm_calls = qmmm_struct%num_qmmm_calls + 1

!  If this is the first call to the routine, do some initial allocation
!  that has not been done elsewhere.
   if (qmmm_struct%qm_mm_first_call) then

     call allocate_qmmm_pair_list( natom )  !Also allocated npair scratch arrays
     allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink), stat=ier )
                !Stores the REAL and link atom qm coordinates
     REQUIRE(ier == 0)

     !Do some initial setup for qm_ewald if in use
     if (qmmm_nml%qm_ewald>0) then
        call timer_start(TIME_QMMMEWALDSETUP)
        qmewald%ewald_startup = .true. !Specify that we haven't done any qm_ewald stuff before this point.
                                        !Essentially that this is the first MD step.
        qmewald%natom = natom !QM ewald needs access to natom in some deep QM routines where it would
                               !not normally be available
        call qm_ewald_setup(qmewald%totkq,qmmm_nml%kmaxqx,qmmm_nml%kmaxqy,qmmm_nml%kmaxqz, &
                            qmmm_nml%ksqmaxq,natom,qmmm_struct%nquant,qmmm_struct%nlink) 
                                               !Also allocates kvec memory (totkq reals)
                                               !,KVec table which is 6 lots of (natom,totkq) reals.
                                               !This is also responsible for dividing up kvectors between
                                               !cpus.
        !If diagnostics are on we can print some info here.
        if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 2) then
           write(6,*) ''
           write(6,*) 'QMMM: Ewald - kmaxq(x,y,z) = ', qmmm_nml%kmaxqx, ',', qmmm_nml%kmaxqy, ',', &
                      qmmm_nml%kmaxqz
           write(6,*) 'QMMM: Ewald -      ksqmaxq = ', qmmm_nml%ksqmaxq
           write(6,*) 'QMMM: Ewald - Total number of k vectors = ', qmewald%totkq
           write(6,*) ''
        end if

        !If we are NOT updating the qm image atom charges on every SCF step for qm_ewald
        !then we initially need to zero the scf_mchg array. Since in this situation on the first step
        !we allow the QM image charges to vary with the SCF. It is only on steps 2 -> N that we
        !keep them fixed during the SCF. This avoids the need for an explicit test of first
        !call in the qm_ewald_calc_mm_pot code.
        if (qmmm_nml%qm_ewald==2) scf_mchg(1:qmmm_struct%nquant_nlink) = zero

        call timer_stop(TIME_QMMMEWALDSETUP)
     else
       qmewald%ewald_startup = .false.
     end if
     
     !Allocation for QM_GB (qmgb==2)
     if (qmmm_nml%qmgb == 2) then
       !Calculate dielectric factor
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = mmcut2
     end if
   end if ! ---- first call endif ----------
   call timer_stop(TIME_QMMMSETUP)

!=============================================================================
!                   BUILD NONBOND LIST for QM region
!=============================================================================
   if(periodic ==1) then
      !---- check periodic status and run some checks on the system
      call timer_start(TIME_QMMMCOORDSX)
      !Note also builds list and starts / stops relevant timers.
      call qm_fill_qm_xcrd_periodic(coords,natom, &
           qmmm_nml%iqmatoms,scaled_mm_charges,qmmm_scratch%qm_real_scratch)

      call timer_stop(TIME_QMMMLISTBUILD) !Since QMMMCOORDSX was stopped in qm_fill_qm_xcrd_periodic
                                          !and QMMMLISTBUILD was started.
   else
      !---------------Not Periodic ------------------------------------
      call timer_start(TIME_QMMMLISTBUILD)
      !---- cutoff based on distance from whole QM region
      !nb_list also fills qm_xcrd array and extracts qm atoms.
      call qm_fill_qm_xcrd( coords, natom, scaled_mm_charges) 
      call timer_stop(TIME_QMMMLISTBUILD)
   endif

   call timer_start(TIME_QMMMCOORDSX)

   !If verbosity is on print the number of QMMM pairs
   if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 1) then
      write(6,*) 'QMMM: No. QMMM Pairs per QM atom: ', qmmm_struct%qm_mm_pairs
   end if

   !Finally we need to position the link atoms for this step
   !We base this on the imaged coordinates so we don't need
   !to worry about periodic boundaries. Writes the link atoms
   !to the end of the qmmm_struct%qm_coords array.
   if ( qmmm_struct%nlink > 0 ) then
     !Note, this is essentially doing work that has
     !already been done (in force) and could be replaced with a simple 
     !array copy with a simple shift - but for now it doesn't hurt to
     !do the work twice.
     call position_link_atoms(coords)
   end if

   if (qmmm_mpi%commqmmm_master .AND. (qmmm_nml%verbosity>1 .OR. qmmm_struct%qm_mm_first_call)) then
      call print_link_atom_info( qmmm_struct%qm_coords,atom_type )
   end if
       
!=============================================================================
!                   START OF REST OF QMMM SETUP
!=============================================================================
   call timer_stop_start(TIME_QMMMCOORDSX,TIME_QMMMSETUP)
   if(qmmm_struct%qm_mm_first_call) then 
       if (qmmm_mpi%commqmmm_master) write(6,'(/80(1H-)/''  3.1 QM CALCULATION INFO'',/80(1H-)/)')
       call qm2_load_params_and_allocate() !Load the parameters
                                !Also does a lot of memory allocation and pre-calculates all
                                !the STO-6G orbital expansions.
       !Print a summary about memory usage
       !WARNING - FOR THE NUMBERS PRODUCED BY THE PRINT DYN MEM ROUTINE TO BE ACCURATE ALL
       !MEMORY ALLOCATION MUST HAVE BEEN DONE BY THIS STAGE.
       if (qmmm_mpi%commqmmm_master) then
          call qm_print_dyn_mem(natom,qmmm_struct%qm_mm_pairs)
          !Print the initial QM region coordinates
          call qm_print_coords(qmmm_struct%qm_coords)
          !Finally print the result header that was skipped in sander.
           write(6,'(/80(1H-)/''   4.  RESULTS'',/80(1H-)/)')
        end if
   end if !if (qmmm_struct%qm_mm_first_call)

   !======================
   !  Setup for QM EWALD
   !======================
   !If we are doing QM Ewald then we need to calculate the Kvectors. We do this on
   !every step since the box dimensions may change and these values depend on the
   !box dimensions.
   if (qmmm_nml%qm_ewald>0) then
     call timer_stop_start(TIME_QMMMSETUP,TIME_QMMMEWALDKTABLE)
!Parallel
     call qm_ewald_calc_kvec(qmewald%kvec, qmewald%dkvec, qmewald%dmkv, qmmm_nml%kmaxqx, &
                             qmmm_nml%kmaxqy, qmmm_nml%kmaxqz,qmmm_nml%ksqmaxq)
     !Next we calculate the KTABLE which is an array of exponentials (in complex
     !Sin, cos notation) in the form of 1->6 by 1->NKvectors wide by 1->Natom (all atoms) long.
     !Note, this routine is very memory intensive and very time consuming. Although if
     !we are using Walker and Crowley PME implementation then we only do the QM-QM
     !table which is fast.
!Parallel
     call qm_ewald_calc_ktable(natom, qmmm_struct%nquant, qmmm_struct%nlink, coords, qmewald%dmkv)
                                 !Will fill the kvector tables with the exponential values
                                 !memory for ktable_x_cos... should already have been allocated.

     call timer_stop_start(TIME_QMMMEWALDKTABLE,TIME_QMMMSETUP)

   end if
   !==========================
   !  End Setup for QM EWALD
   !==========================

   call timer_stop(TIME_QMMMSETUP)
!======================END OF QMMM SETUP ======================================

   !Calculate RIJ and many related equations here. Necessary memory allocation
   !is done inside the routine.
   call timer_start(TIME_QMMMRIJEQNS)
!Parallel
   call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, qmmm_struct%nquant_nlink,qmmm_struct%qm_xcrd, &
                              natom, qmmm_struct%qm_mm_pairs)
                                !and store them in memory to save time later.
   call timer_stop_start(TIME_QMMMRIJEQNS,TIME_QMMMENERGY)

   !============================
   ! Start Calculate SCF Energy
   !============================
!Parallel
   call qm2_energy(escf, scf_mchg, natom,born_radii, one_born_radii, coords, scaled_mm_charges)
   !============================
   ! End Calculate SCF Energy
   !============================

   call timer_stop(TIME_QMMMENERGY)

   !=============================
   ! Start Calculation of Forces
   !=============================

   !CALCULATE FORCES USING DENSITY MATRIX FROM SCF
   call timer_start(TIME_QMMMFQM)
   qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory==DFTB) then
     !We are doing DFTB
!Partially Parallel
     call qm2_dftb_get_qm_forces(qmmm_struct%dxyzqm)
   else
     !standard semi-empirical
!Parallel
     call qm2_get_qm_forces(qmmm_struct%dxyzqm)
   end if

   call timer_stop(TIME_QMMMFQM)

   call timer_start(TIME_QMMMFQMMM)

   if (qmmm_nml%qmmm_int>0) then
      qmmm_struct%dxyzcl=zero
      if (qmmm_nml%qmtheory == DFTB) then
        !We are doing DFTB
!Parallel
        call qm2_dftb_get_qmmm_forces(qmmm_struct%dxyzcl,qmmm_struct%dxyzqm, qmmm_scratch%qm_real_scratch(1), &
                                  qmmm_scratch%qm_real_scratch(natom+1), &
                                  qmmm_scratch%qm_real_scratch(2*natom+1), &
                                  qmmm_scratch%qm_real_scratch(3*natom+1))
      else
!Parallel
        call qm2_get_qmmm_forces(qmmm_struct%dxyzqm,qmmm_struct%qm_xcrd,qmmm_struct%dxyzcl)
      end if
   end if

   call timer_stop(TIME_QMMMFQMMM)

!If we are doing qm_ewald then we need to calculate the gradients due to the ewald potential here.
!Only available analytically, no numerical gradients are available.
!With qm_pme the forces are calculated in ew_recip using the mulliken charges
   if (qmmm_nml%qm_ewald>0) then
     call timer_start(TIME_QMMMFQMEWALD)
!Parallel
     call qm_ewald_get_forces(qmmm_struct%qm_xcrd, qmmm_struct%qm_coords,&
                              natom, scf_mchg, &
                              qmmm_nml%qmmmrij_incore, &
                              qmmm_struct%dxyzqm, qmmm_struct%dxyzcl, &
                              qmewald%dkvec, scaled_mm_charges)
     call timer_stop(TIME_QMMMFQMEWALD)
   end if

   !=============================
   ! End Calculation of Forces
   !=============================

!NOW WE NEED TO PUT THE CALCULATED FORCES INTO THE SANDER FORCE ARRAY
   call timer_start(TIME_QMMMCOLLATEF)
   do i=1,qmmm_struct%nquant
     m = qmmm_nml%iqmatoms(i)
     m = (m-1)*3
     f(m+1) = f(m+1) - qmmm_struct%dxyzqm(1,i)
     f(m+2) = f(m+2) - qmmm_struct%dxyzqm(2,i)
     f(m+3) = f(m+3) - qmmm_struct%dxyzqm(3,i)
   enddo
!Only need to do MM atoms that are in the list. Only need to do this if the QMMM interaction
!is being calculated using the full orbital interaction.
   if (qmmm_nml%qmmm_int >0) then
     do i = 1,qmmm_struct%qm_mm_pairs
        m = (qmmm_struct%qm_mm_pair_list(i)-1)*3
        f(m+1) = f(m+1) - qmmm_struct%dxyzcl(1,i)
        f(m+2) = f(m+2) - qmmm_struct%dxyzcl(2,i)
        f(m+3) = f(m+3) - qmmm_struct%dxyzcl(3,i)
     end do
   end if


   if (qmmm_nml%qm_ewald>0 .and. .not. qmmm_nml%qm_pme) then
!If we are doing QM ewald then we have a set of forces on all MM atoms that we need to put
!into the main force array. For qm_pme this is done in ew_recip later on in force.
!Parallel division in the ewald code is over kvectors so all threads have some forces on all
!atoms.

     do i = 1, natom
       offset = (3*i) - 2
       f(offset)   = f(offset)   - qmewald%d_ewald_mm(1,i)
       f(offset+1) = f(offset+1) - qmewald%d_ewald_mm(2,i)
       f(offset+2) = f(offset+2) - qmewald%d_ewald_mm(3,i)
     end do
   end if

   !Return the mm link pair coordinates back to the original mm atoms 
   !Needed before distributing link atom forces since it uses the unimaged array.
   call rst_mm_link_pair_crd(coords)
   
   !We need to divide the force on the link atom up between the
   !QM and MM atom.
   do i=1,qmmm_struct%nlink
     mm_no = 3*qmmm_struct%link_pairs(1,i)-2  !location of atom in x array
     lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom bound to link atom
     qm_no = 3*qmmm_nml%iqmatoms(lnk_no)-2
     !Note this routine uses the flink in the form -flink. 
     call distribute_lnk_f(forcemod,qmmm_struct%dxyzqm(1,qmmm_struct%nquant+i),coords(mm_no), &
                           coords(qm_no),qmmm_nml%lnk_dis)

     !NOTE: forces are reversed in QM calc with respect to amber force array
     !so we subtract forcemod from MM atom and add it to QM atom.
     !MM atom's new force = FMM(x,y,z) - FORCEMOD(x,y,z)
     f(mm_no) = f(mm_no) - forcemod(1)
     f(mm_no+1) = f(mm_no+1) - forcemod(2)
     f(mm_no+2) = f(mm_no+2) - forcemod(3)

     !QM atom's new force = FQM(x,y,z) - Flink(x,y,z) + FORCEMOD(x,y,z)
     !Note QM forces should be subtracted from sander F array to leave total force.
     f(qm_no) = f(qm_no) - qmmm_struct%dxyzqm(1,qmmm_struct%nquant+i) + forcemod(1)
     f(qm_no+1) = f(qm_no+1) - qmmm_struct%dxyzqm(2,qmmm_struct%nquant+i) + forcemod(2)
     f(qm_no+2) = f(qm_no+2) - qmmm_struct%dxyzqm(3,qmmm_struct%nquant+i) + forcemod(3)

   end do

   call timer_stop(TIME_QMMMCOLLATEF)

   !=============================================================
   !   Calculate Mulliken Charges for Later Printing if Needed
   !   Note: at present we calculate the mulliken charges even
   !         if we don't need them for printing since one might
   !         want to use them for dipole calculations etc. It is
   !         not very expensive to calculate them so might as well
   !         do it on every MD step.
   !=============================================================
   if (qmmm_mpi%commqmmm_master) then
      if ( qmmm_nml%qmtheory == DFTB .or. qm2_struct%calc_mchg_scf &
         .or. qmmm_nml%qm_ewald == 2) then
        !Mulliken charges have already been calculated and stored.
      else
        do i=1,qmmm_struct%nquant_nlink
          !Need to calculate Mulliken charges here.
          call qm2_calc_mulliken(i,scf_mchg(i))
        end do
      end if
   end if
   !=============================================================
   ! End Calculate Mulliken Charges for Later Printing if Needed
   !=============================================================

!Print some extra information if verbosity level is > 0
   if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 0) then
      !Verbosity level of 1 or more = print more accurate SCF energy
      write (6,'("QMMM:")')
      write (6,'("QMMM: SCF Energy =",f22.14," KCal/mol, ",f22.14," KJ/mol")') escf, escf*4.184d0
      !If verbosity level is greater than 1 we also print the nuclear and electronic energies.
      if (qmmm_nml%verbosity > 1) then
         write (6,'("QMMM:")')
         write (6,'("QMMM:        Electronic energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
              qmmm_struct%elec_eng, qmmm_struct%elec_eng*EV_TO_KCAL
         if  ( qmmm_nml%qmtheory == DFTB ) then
            write (6,'("QMMM:         Repulsive energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                 qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
!!            write (6,'("QMMM:        Careful: Dispersion Energy Already Included.")')
!!            write (6,'("QMMM:        Dispersion Energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
!!                 dftb_edisp*AU_TO_KCAL
            total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmqm
         else
           write (6,'("QMMM: QM core - QM core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
           write (6,'("QMMM: QM core - MM atom energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmmm,qmmm_struct%enuclr_qmmm*EV_TO_KCAL
           write (6,'("QMMM: Total core - core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                  qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm, &
                  (qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm)*EV_TO_KCAL
           total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmmm + qmmm_struct%enuclr_qmqm
         end if
         write (6,'("QMMM:             Total energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                total_energy, total_energy*EV_TO_KCAL
         !If verbosity level is greater than 3 we also print the force array on the QM atoms
         if (qmmm_nml%verbosity > 3) then
            write (6,'("QMMM:")')
            write (6,'("QMMM: Forces on QM atoms from SCF calculation")')
            write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j), qmmm_struct%dxyzqm(2,j), &
                   qmmm_struct%dxyzqm(3,j), j=1,qmmm_struct%nquant_nlink)
            if (qmmm_nml%verbosity > 4) then
               !Also print info in KJ/mol
               write (6,'("QMMM:")')
               write (6,'("QMMM: Forces on QM atoms from SCF calculation (KJ/mol)")')
               write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j)*4.184d0, &
                      qmmm_struct%dxyzqm(2,j)*4.184d0, qmmm_struct%dxyzqm(3,j)*4.184d0, &
                      j=1,qmmm_struct%nquant_nlink)
            end if   
         end if
      end if
   end if

   qmmm_struct%qm_mm_first_call = .false.
   qmewald%ewald_startup = .false.

   return

end subroutine qm_mm

!======================END OF QM_MM ======================================

!============== NON PERIODIC QM_XCRD, QM_COORDS and LIST =================
!Build non-periodic list - fill qm_xcrd and extract qm_coords.
subroutine qm_fill_qm_xcrd( x, natom, scaled_mm_charges ) 
   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_scratch
   implicit none
!#  include "memory.h" 
   !This routine calculates a qm_mm_pair_list which is of the form
   !listtype atom atom atom ... 
   !for qmmm_struct%nquant atoms (I.e each QM atom excluding any link atoms

   !The list of MM atoms that interact with each QM atom is based
   !on the cut off distance. It then fills qm_xcrd.

   !It also extracts the qm coordinates from the amber x array.

   !Each QM atom will use the identical list since to be included an MM atom
   !only needs to be within cut of any QM atom

   !qm_mm_pair_list comes from qmmm_module.

!Passed in
   integer natom
   _REAL_ , intent(in) ,dimension(3,natom) :: x
   _REAL_, intent(in), dimension(natom) :: scaled_mm_charges

   !Local Variables!
   integer j,m,i,n1
   _REAL_ , dimension(6) :: bxbnd
   _REAL_ x_qm, y_qm, z_qm, dx2, xbnd0, xbnd1, ybnd0, ybnd1, zbnd0, zbnd1
   logical include_atom

!     Find the bounding box limits of the QM region, then find all atoms
!     inside the box + cutoff before calculating or testing distances.
   m=qmmm_nml%iqmatoms(1)
   xbnd0=x(1,m)
   xbnd1=x(1,m)
   ybnd0=x(2,m)
   ybnd1=x(2,m)
   zbnd0=x(3,m)
   zbnd1=x(3,m)
   do j=2,qmmm_struct%nquant
      m = qmmm_nml%iqmatoms(j)
      xbnd0=min(xbnd0,x(1,m))
      xbnd1=max(xbnd1,x(1,m))
      ybnd0=min(ybnd0,x(2,m))
      ybnd1=max(ybnd1,x(2,m))
      zbnd0=min(zbnd0,x(3,m))
      zbnd1=max(zbnd1,x(3,m))
   enddo

   bxbnd(1)=xbnd0-qmmm_nml%qmcut
   bxbnd(2)=xbnd1+qmmm_nml%qmcut
   bxbnd(3)=ybnd0-qmmm_nml%qmcut
   bxbnd(4)=ybnd1+qmmm_nml%qmcut
   bxbnd(5)=zbnd0-qmmm_nml%qmcut
   bxbnd(6)=zbnd1+qmmm_nml%qmcut

   include_atom = .false.

   n1 = 0  ! index of qm_mm_pair_list to which we will put the current MM atom
 !---------- FIRST PASS ----------------------------------------------
   qmmm_scratch%qm_int_scratch(1:natom)=0 !Used for a mask
   do m=1,natom 
      !No short circuit evaluation in fortran so having a series of
      !sperate if statements here should be faster.
      if ( x(1,m) <= bxbnd(1) ) cycle
      if ( x(1,m) >= bxbnd(2) ) cycle
      if ( x(2,m) <= bxbnd(3) ) cycle
      if ( x(2,m) >= bxbnd(4) ) cycle
      if ( x(3,m) <= bxbnd(5) ) cycle
      if ( x(3,m) >= bxbnd(6) ) cycle
          
      qmmm_scratch%qm_int_scratch(m)=1

   enddo
   !Set QM atoms in qm_int_scratch to zero so they don't get counted as pairs.
   qmmm_scratch%qm_int_scratch(qmmm_nml%iqmatoms(1:qmmm_struct%nquant)) = 0

 !---------- Second PASS ----------------------------------------------
   do m=1,natom          ! we loop over all atoms - excluding QM atoms
       if ( qmmm_scratch%qm_int_scratch(m) == 1 ) then
          check_cut: do j=1,qmmm_struct%nquant
             i = qmmm_nml%iqmatoms(j) ! the atom number of the first QM atom
                              ! Find all MM atoms that are within CUT of this QM atom
             x_qm = x(1,i)    ! Get coordinates of QM atom
             y_qm = x(2,i)
             z_qm = x(3,i)
                              ! calculate the distance and see if it is within the cut off
             dx2 = ( ( x_qm - x(1,m) ) * ( x_qm - x(1,m) ) &
                   + ( y_qm - x(2,m) ) * ( y_qm - x(2,m) ) &
                   + ( z_qm - x(3,m) ) * ( z_qm - x(3,m) ) &
                   ) 
             if ( dx2 < qmmm_nml%qmcut2 ) then
                !We include the atom.
                !however, if this is a mm link pair
                !atom then we don't include it.
                if (qmmm_struct%mm_link_mask(m)) then
                  !We don't include it since it is a link atom
                  exit check_cut
                end if
                include_atom = .true.
                exit check_cut
             end if
          end do check_cut

          if ( include_atom ) then
                              !include this mm atom in the list
             n1 = n1+1
             qmmm_struct%qm_mm_pair_list( n1 ) = m
             qmmm_struct%qm_xcrd(1:3,n1) = x(1:3,m)
             qmmm_struct%qm_xcrd(4,n1) = scaled_mm_charges(m) 
             include_atom=.false.
          end if
       end if
    end do
    qmmm_struct%qm_mm_pairs = n1

    !Extract QM atoms from x into qmcoords array.
    do m=1,qmmm_struct%nquant
        i=qmmm_nml%iqmatoms(m)
        qmmm_struct%qm_coords(1:3,m) = x(1:3,i)
    enddo

    return

end subroutine qm_fill_qm_xcrd

!=============================================================================
!             QM_FILL_QM_XCRD PERIODIC
!=============================================================================
subroutine qm_fill_qm_xcrd_periodic(x,natom, &
                                    iqmatoms,scaled_mm_chrgs,real_scratch)

   use qmmm_module, only: qmmm_nml,qmmm_struct, qmmm_scratch
   use nblist, only: a,b,c,alpha,beta,gamma,ucell,recip,sphere
   use constants, only : zero, one, half, two

   implicit none
   integer , intent(in) :: natom
   _REAL_ , intent(in), dimension(3,natom) :: x
   integer, intent(in), dimension(qmmm_struct%nquant) :: iqmatoms
   _REAL_, intent(in), dimension(natom) :: scaled_mm_chrgs
   _REAL_, intent(out), dimension(3,natom) :: real_scratch 

   integer j,iqm_one, jqmatom, ier
   _REAL_ :: xbnd0,xbnd1,ybnd0,ybnd1,zbnd0,zbnd1
   _REAL_ :: offset(3), frac(3), xx,yy,zz
   _REAL_ , dimension(6) :: bxbnd

   integer :: m,n, n1
   _REAL_ :: fbndx0, fbndx1, fbndy0, fbndy1, fbndz0, fbndz1, dx2
   logical :: include_atom

   real_scratch(1:3,1) = zero

!Move the first QM atom to be at the origin.
   iqm_one = iqmatoms(1)
   frac(1:3) = x(1,iqm_one)*recip(1,1:3)+x(2,iqm_one)*recip(2,1:3)+ &
               x(3,iqm_one)*recip(3,1:3)
   frac(1:3) = frac(1:3) - anint(frac(1:3))
   offset(1) = frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
   offset(2) = frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
   offset(3) = frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)

   do j = 1, qmmm_struct%nquant
      jqmatom = iqmatoms(j)
      xx = x(1,jqmatom) - offset(1)
      yy = x(2,jqmatom) - offset(2)
      zz = x(3,jqmatom) - offset(3)
      frac(1:3) = xx*recip(1,1:3) + yy*recip(2,1:3) + zz*recip(3,1:3)
      frac(1:3) = frac(1:3) - anint(frac(1:3))
      real_scratch(1,j)=frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
      real_scratch(2,j)=frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
      real_scratch(3,j)=frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)
   end do

   xbnd0=zero; xbnd1=zero; ybnd0=zero; ybnd1=zero; zbnd0=zero; zbnd1=zero
   
   do j=1,qmmm_struct%nquant
      xbnd0=min(xbnd0,real_scratch(1,j))
      xbnd1=max(xbnd1,real_scratch(1,j))
      ybnd0=min(ybnd0,real_scratch(2,j))
      ybnd1=max(ybnd1,real_scratch(2,j))
      zbnd0=min(zbnd0,real_scratch(3,j))
      zbnd1=max(zbnd1,real_scratch(3,j))
   enddo
   xbnd0=xbnd0-qmmm_nml%qmcut
   ybnd0=ybnd0-qmmm_nml%qmcut
   zbnd0=zbnd0-qmmm_nml%qmcut
   xbnd1=xbnd1+qmmm_nml%qmcut
   ybnd1=ybnd1+qmmm_nml%qmcut
   zbnd1=zbnd1+qmmm_nml%qmcut

  !------ Check if QM region plus cutoff around it is too large for this box
  !---       The sphere method is the simplest check, but we can get more
  !---       sophisticated using the distance between parallel faces later if we need to...mfc
  !---       sphere is calculated in ew_box.f and used here.
   if( (xbnd1-xbnd0 > two*sphere) .or. (ybnd1-ybnd0 > two*sphere) .or. (zbnd1-zbnd0 > two*sphere) )then
      write(6,*) " ****************************************************"
      write(6,*) " ERROR: QM region + cutoff larger than box dimension:"
      write(6,'(2X,"QM-MM Cutoff = ",f8.4)') qmmm_nml%qmcut
      write(6,*) "  Coord   Lower     Upper    Size    Radius of largest sphere inside unit cell"
      write(6,'(5X,"X",4(2X,f8.3))') xbnd0, xbnd1, xbnd1-xbnd0, sphere
      write(6,'(5X,"Y",4(2X,f8.3))') ybnd0, ybnd1, ybnd1-ybnd0, sphere
      write(6,'(5X,"Z",4(2X,f8.3))') zbnd0, zbnd1, zbnd1-zbnd0, sphere
      write(6,*) " ****************************************************"
      call sander_bomb("QM_CHECK_PERIODIC<qm_mm.f>", &
        "QM region + cutoff larger than box", &
        "cannot continue, need larger box.")
   endif

   bxbnd(1)=xbnd0+offset(1)
   bxbnd(2)=xbnd1+offset(1)
   bxbnd(3)=ybnd0+offset(2)
   bxbnd(4)=ybnd1+offset(2)
   bxbnd(5)=zbnd0+offset(3)
   bxbnd(6)=zbnd1+offset(3)

!FILL QM_XCRD Periodic
  !---- First move center of QM region to origin: find offset
   offset(1)=(bxbnd(2)+bxbnd(1))*half
   offset(2)=(bxbnd(4)+bxbnd(3))*half
   offset(3)=(bxbnd(6)+bxbnd(5))*half

!  !---- Create new fractional coordinates with new origin
!    !--- find bounds and cut in fracs
!   fbndx0 = (bxbnd(1)-offset(1))*bxinv(1)
!   fbndx1 = (bxbnd(2)-offset(1))*bxinv(1)
!   fbndy0 = (bxbnd(3)-offset(2))*bxinv(2)
!   fbndy1 = (bxbnd(4)-offset(2))*bxinv(2)
!   fbndz0 = (bxbnd(5)-offset(3))*bxinv(3)
!   fbndz1 = (bxbnd(6)-offset(3))*bxinv(3)

    !---- run through coords to get fracs, select those within
    !     bounding box + cutoff

   qmmm_scratch%qm_int_scratch(1:natom) = 0  !Used for a mask
   do m=1,natom
      xx=x(1,m)-offset(1)
      yy=x(2,m)-offset(2)
      zz=x(3,m)-offset(3)
      frac(1:3) = xx*recip(1,1:3) + yy*recip(2,1:3) + zz*recip(3,1:3)
      frac(1:3) = frac(1:3) - anint(frac(1:3))
      xx = frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
      if( (xx > bxbnd(1)-offset(1)) .and. (xx < bxbnd(2)-offset(1)) ) then
         yy=frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
         if( ( yy > bxbnd(3)-offset(2)) .and. (yy < bxbnd(4)-offset(2)) ) then
            zz=frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)
            if( (zz > bxbnd(5)-offset(3)) .and. (zz < bxbnd(6) -offset(3)) ) then
               !------- This one inside box ------------------
               qmmm_scratch%qm_int_scratch(m)=1
               real_scratch(1,m) = xx
               real_scratch(2,m) = yy
               real_scratch(3,m) = zz
            endif
         endif
      endif
   enddo

!Fill the qm coordinate array with the imaged QM atoms
   do n=1,qmmm_struct%nquant
      m=iqmatoms(n)
      qmmm_struct%qm_coords(1:3,n) = real_scratch(1:3,m)
      qmmm_scratch%qm_int_scratch(iqmatoms(n)) = 0
   enddo
   call timer_stop_start(TIME_QMMMCOORDSX,TIME_QMMMLISTBUILD)
!Now calculate the pair list and fill qm_xcrd
   include_atom = .false.
   n1 = 0
 !---------- Second PASS ----------------------------------------------
   do m=1,natom          ! we loop over all atoms - excluding QM atoms
      if ( qmmm_scratch%qm_int_scratch(m) > 0 ) then
         check_cut: do j=1,qmmm_struct%nquant
            ! Find all MM atoms that are within CUT of this QM atom
            xx = qmmm_struct%qm_coords(1,j) ! Get coordinates of QM atom
            yy = qmmm_struct%qm_coords(2,j)
            zz = qmmm_struct%qm_coords(3,j)
            ! calculate the distance and see if it is within the cut off
            dx2 = ( ( xx - real_scratch(1,m) ) * ( xx - real_scratch(1,m) ) &
               + ( yy - real_scratch(2,m) ) * ( yy - real_scratch(2,m) ) &
               + ( zz - real_scratch(3,m) ) * ( zz - real_scratch(3,m) ) &
               )
            if ( dx2 < qmmm_nml%qmcut2 ) then
              !We include the atom.
              !however, if this is a mm link pair
              !atom then we don't include it.
              if (qmmm_struct%mm_link_mask(m)) then
                !We don't include it since it is a link atom
                exit check_cut
              end if
              include_atom = .true.
              exit check_cut
            end if
         end do check_cut
         if ( include_atom ) then
            n1 = n1+1
            qmmm_struct%qm_mm_pair_list( n1 ) = m
            qmmm_struct%qm_xcrd(1,n1)=real_scratch(1,m)
            qmmm_struct%qm_xcrd(2,n1)=real_scratch(2,m)
            qmmm_struct%qm_xcrd(3,n1)=real_scratch(3,m)
            qmmm_struct%qm_xcrd(4,n1)=scaled_mm_chrgs(m)
            include_atom=.false.
         end if
      end if                    ! endif of skip_m
   end do
   qmmm_struct%qm_mm_pairs = n1


   return
end subroutine qm_fill_qm_xcrd_periodic

