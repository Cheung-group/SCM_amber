! <compile=optimized>
#include "copyright.h"
#include "dprec.h"

!This file contains several routines for printing information
!about QMMM simulations.

!------------------------------------------------------------------------------
subroutine qm_print_coords()
!This routine prints the qm region coordinates.
!Including link atoms
!
!Written by Ross Walker, TSRI, 2004
!
!================================================

  use qmmm_module, only : qmmm_struct, element_sym, qmmm_nml
  implicit none

!Local
      integer i,j

      write(6,'(/,''  QMMM: QM Region Cartesian Coordinates (*=link atom) '')')
      write(6,'(''  QMMM: QM_NO.   MM_NO.'',2X,''ATOM'',9X,''X'',9X,''Y'',9X,''Z'')')
      do I=1,qmmm_struct%nquant
         WRITE(6,'("  QMMM:",I6,2X,I7,6X,A2,3X,3F10.4)') I, &
               qmmm_nml%iqmatoms(i), &
               element_sym(qmmm_struct%iqm_atomic_numbers(i)), &
               (qmmm_struct%qm_coords(J,I),J=1,3)
      end do
      do i=qmmm_struct%nquant+1,qmmm_struct%nquant_nlink
         WRITE(6,'("  QMMM:",I6,2X,7X,5X,"*",A2,3X,3F10.4)') I, &
               element_sym(qmmm_struct%iqm_atomic_numbers(i)), &
               (qmmm_struct%qm_coords(J,I),J=1,3)
      end do

      if (qmmm_nml%writepdb) then
        !We write a pdb of the coordinates
        call qm_write_pdb()
      end if

end subroutine qm_print_coords

!------------------------------------------------------------------------------

subroutine qm_write_pdb()
! This routine will write a crude pdb of the coordinates representing
! the QM atoms and link atoms. This allows the user to visually check that
! the selected QM region is what they expected.

! Written by Ross Walker (TSRI, 2005)

  use qmmm_module, only : qmmm_struct, element_sym
  implicit none

!Local
  integer :: i,j

  call amopen(23,'qmmm_region.pdb','R','F','W')
  write(23,'(a)') 'REMARK'
  do i=1,qmmm_struct%nquant_nlink
    write(23,60) i,element_sym(qmmm_struct%iqm_atomic_numbers(i)),'QM',1, &
                 (qmmm_struct%qm_coords(J,I),J=1,3)
  end do
  write(23,'(a)') 'END'

  close(unit=23) 

  60 FORMAT('ATOM',2X,I5,1X,A4,1X,A4,1X, I4,4X,3F8.3)

end subroutine qm_write_pdb

!------------------------------------------------------------------------------

subroutine qm_print_dyn_mem(natom,npairs)
!This routine prints a summary of the dynamic
!memory allocated for use in the QM calculation.
!Written by Ross Walker, TSRI, 2004
!Assumes _REAL_ is double precision = 8 bytes.
!================================================

  use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, qm2_rij_eqns, &
                          qmewald, qm_gb, qmmm_mpi, qmmm_scratch, &
                          AM1, PM3, PDDGPM3, PDDGMNDO, PM3CARB1, DFTB, RM1
  implicit none

#include "qm2_array_locations.h"

!Passed in
  integer, intent(in) :: natom, npairs

!Local
  integer total_memory, element_memory
  _REAL_, parameter :: bytes_to_mb = 1.0D0/(1024D0*1024D0)
  integer bytes_per_int, bytes_per_real, bytes_per_logical
  total_memory = 0
  bytes_per_int = bit_size(element_memory)/8 
  bytes_per_real = 8     !Assume size of _REAL_ is 8 bytes
  bytes_per_logical = 1  !Assume size of logical is a single byte

  write(6,'(/"| QMMM: Estimated QM Dynamic Memory Usage (per thread)")')
  write(6,'("| QMMM: ---------------------------------------------------")')

  element_memory = size(qmmm_struct%qm_atom_type)*bytes_per_real &
                  +size(qmmm_struct%qm_type_id)*bytes_per_real
  write(6,'("| QMMM:              QM Atom Type Info : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%qm_resp_charges)*bytes_per_real
  write(6,'("| QMMM:         QM RESP Charge Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_nml%iqmatoms)*bytes_per_int
  write(6,'("| QMMM:            QM Atom Number List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = 0
  if (qmmm_struct%nlink > 0) then  
    element_memory = size(qmmm_struct%link_pairs)*bytes_per_int
    element_memory = element_memory + &
                                 size(qmmm_struct%mm_link_pair_resp_charges) * bytes_per_real
    element_memory = element_memory + &
                                 size(qmmm_struct%mm_link_pair_saved_coords) * bytes_per_real
  end if
  write(6,'("| QMMM:                Link Atom Pairs : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if(qmmm_nml%peptide_corr .and. qm2_struct%n_peptide_links>0) then
    element_memory = size(qm2_struct%peptide_links)*bytes_per_int
    write(6,'("| QMMM:       Peptide Linkage Identity : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if

  element_memory = size(qmmm_struct%iqm_atomic_numbers)*bytes_per_int
  write(6,'("| QMMM:          QM Atomic Number List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  !Pair list is always allocated as at least 1 even if natom = nquant.
  element_memory = size(qmmm_struct%qm_mm_pair_list)*bytes_per_int
  write(6,'("| QMMM:                QM-MM Pair List : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = (size(qmmm_struct%atom_mask)+size(qmmm_struct%mm_link_mask))*bytes_per_logical
  write(6,'("| QMMM:                   QM Atom Mask : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%qm_coords)*bytes_per_real
!include the qm_xcrd(4,natom) array here. This is declared
!in qm_mm.f
  element_memory = element_memory+4*natom*bytes_per_real
  write(6,'("| QMMM:           QM Coordinate Arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_struct%scaled_mm_charges)*bytes_per_real
  write(6,'("| QMMM:         Scaled MM Charge Array : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  !SCF Mulliken Charge array
  element_memory = size(qm2_struct%scf_mchg)*bytes_per_real
  write(6,'("| QMMM:    SCF Mulliken Charge Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory  
  
  !qm ewald arrays
  if (qmmm_nml%qm_ewald>0) then
    element_memory = size(qmewald%kvec)+size(qmewald%dkvec)+size(qmewald%dmkv)
    if (.not. qmmm_nml%qm_pme) element_memory= element_memory + size(qmewald%ktable) + &
                                               size(qmewald%d_ewald_mm)
    element_memory = element_memory + size(qmewald%qmktable) + &
                     size(qmewald%mmpot) + size(qmewald%qmpot)
    element_memory = element_memory * bytes_per_real
    write(6,'("| QMMM:                QM Ewald Arrays : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if

  !QM GB arrays
  if (qmmm_nml%qmgb==2) then
    element_memory = size(qm_gb%qmqm_onefij)+size(qm_gb%gb_mmpot)+size(qm_gb%gb_qmpot)
    if (qm_gb%saltcon_on) element_memory = element_memory + size(qm_gb%qmqm_kappafij)
    element_memory = element_memory*bytes_per_real
    element_memory = element_memory + size(qm_gb%qmqm_gb_list)*bytes_per_int
    write(6,'("| QMMM:                   QM GB Arrays : ",i12," bytes")') element_memory
    total_memory = total_memory + element_memory
  end if
   
!Include the local force matrices, dxyzqm and dxyzcl here as well
  element_memory = 3*natom*bytes_per_real
  if (qmmm_nml%qmmm_int /= 0 ) then
    element_memory =element_memory+3*qmmm_struct%nquant_nlink*bytes_per_real
  end if
  write(6,'("| QMMM:                QM Force Arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%den_matrix)*bytes_per_real
  if (qmmm_nml%density_predict == 1) then
    element_memory = element_memory + (( size(qm2_struct%md_den_mat_guess1) + &
                                        size(qm2_struct%md_den_mat_guess2) ) * bytes_per_real)
  end if
  write(6,'("| QMMM:                 Density Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = (size(qm2_struct%old_den_matrix)+size(qm2_struct%old2_density))* &
                   bytes_per_real
  write(6,'("| QMMM:          Density Matrix Copies : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = qmmm_struct%nquant_nlink*16*bytes_per_real
  write(6,'("| QMMM: Fock2 Density Matrix Workspace : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%fock_matrix)*bytes_per_real
  write(6,'("| QMMM:                    Fock Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmtheory == DFTB) then
    element_memory = 0
  else
    element_memory = size(qm2_struct%eigen_vectors)*bytes_per_real
  end if
  write(6,'("| QMMM:           Eigen Vector Storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmqm_erep_incore) then
     element_memory = size(qm2_struct%qm_qm_e_repul)*bytes_per_real
  else
     element_memory = 0
  end if
  write(6,'("| QMMM: QM-QM Elec Repulsion Integrals : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%qm_qm_2e_repul)*bytes_per_real
  write(6,'("| QMMM:  QM 2-Elec Repulsion Integrals : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_struct%hmatrix)*bytes_per_real
  write(6,'("| QMMM:              1-Electron Matrix : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_params%core_chg) + size(qm2_params%betasas) + &
                      size(qm2_params%orb_elec_ke) + &
                      size(qm2_params%onec2elec_params) + &
                      size(qm2_params%multip_2c_elec_params) + &
                      size(qm2_params%betasap) + size(qm2_params%betapap)
  if (qmmm_nml%qmtheory /= DFTB) then
    element_memory = element_memory + size(qm2_params%atom_orb_zz_sxs_over_sas) + &
                      size(qm2_params%atom_orb_zz_sxp_over_sap) + &
                      size(qm2_params%atom_orb_zz_pxp_over_pap) + &
                      size(qm2_params%atom_orb_ss_eqn) + &
                      size(qm2_params%atom_orb_sp_ovlp) + &
                      size(qm2_params%atom_orb_pp_ovlp_inj) + &
                      size(qm2_params%atom_orb_pp_ovlp_ieqj1) + &
                      size(qm2_params%atom_orb_pp_ovlp_ieqj2) + &
                      size(qm2_params%atom_orb_ss_eqn_adb) + &
                      size(qm2_params%atom_orb_sp_eqn_xy)+ &
                      size(qm2_params%atom_orb_sp_eqn_xx1) + &
                      size(qm2_params%atom_orb_sp_eqn_xx2) + &
                      size(qm2_params%atom_orb_pp_eqn_xxy1) + &
                      size(qm2_params%atom_orb_pp_eqn_xxy2)
  end if
  if (qmmm_nml%qmtheory==PDDGPM3 .OR. qmmm_nml%qmtheory==PDDGMNDO) then
    element_memory = element_memory + size(qm2_params%pddge1) + size(qm2_params%pddge2)+ &
                     size(qm2_params%pddg_term1) +size(qm2_params%pddg_term2) + &
                     size(qm2_params%pddg_term3)+ size(qm2_params%pddg_term4)
  end if
  if (qmmm_nml%qmtheory==AM1 .OR. qmmm_nml%qmtheory==PM3 .OR. qmmm_nml%qmtheory==PDDGPM3 .OR. &
      qmmm_nml%qmtheory==PM3CARB1 .OR. qmmm_nml%qmtheory==RM1 ) then
    element_memory = element_memory + size(qm2_params%FN1) + size(qm2_params%FN2) + &
                     size(qm2_params%FN3)
  end if
  element_memory = element_memory * bytes_per_real
  write(6,'("| QMMM:       _REAL_ parameter storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qm2_params%natomic_orbs)+size(qm2_params%orb_loc)+ &
                   size(qm2_params%pascal_tri1)+size(qm2_params%pascal_tri2) 
  element_memory = element_memory*bytes_per_int
  write(6,'("| QMMM:      integer parameter storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmqmrij_incore) then
#ifdef MPI
    element_memory = qmmm_mpi%nquant_nlink_loop_extent_end-qmmm_mpi%nquant_nlink_loop_extent_begin+1
#else
    element_memory = qmmm_struct%nquant_nlink * (qmmm_struct%nquant_nlink-1)/2
#endif
    if (qmmm_nml%qmtheory == PDDGPM3 .OR. qmmm_nml%qmtheory == PDDGMNDO) then
      element_memory = element_memory*(QMQMNORIJ+QMQMNOPDDG)*bytes_per_real
    else
      element_memory = element_memory*QMQMNORIJ*bytes_per_real
    end if
  else
    element_memory = 0
  end if
  write(6,'("| QMMM:         QM-QM RIJ Eqns storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  if (qmmm_nml%qmmmrij_incore) then
    element_memory = min(npairs*(qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)+npairs, &
                                 (qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)*(natom-qmmm_struct%nquant_nlink))
    element_memory = QMMMNORIJ*element_memory*bytes_per_real
  else
    element_memory = 0
  end if
  write(6,'("| QMMM:         QM-MM RIJ Eqns storage : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory
  
  if (qmmm_nml%qmtheory == DFTB) then
    element_memory = 0
  else 
    !qmmm_scratch%mat_diag_workspace
    element_memory=size(qmmm_scratch%mat_diag_workspace)*bytes_per_real
  end if

  element_memory = element_memory+size(qmmm_scratch%qm_real_scratch)*bytes_per_real
#ifdef MPI
# ifndef USE_MPI_IN_PLACE
  !Add in the reduction matsize temp array - only allocated if we can't do MPI_IN_PLACE
  element_memory = element_memory + size(qmmm_scratch%matsize_red_scratch)*bytes_per_real
# endif
#endif
  !Add in any scratch space for divide and conquer lapack diagonaliser
  element_memory = element_memory + qmmm_scratch%lapack_dc_real_scr_aloc*bytes_per_real
  !Add in scratch space for the vectmp by npairs arrays - these are always
  !allocated even if there are no mm atoms in which case they have size 1.
  if (qmmm_nml%allow_pseudo_diag) then
    element_memory = element_memory + size(qmmm_scratch%pdiag_scr_norbs_norbs)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_scr_noccupied_norbs)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp1)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp2)*bytes_per_real
    element_memory = element_memory + size(qmmm_scratch%pdiag_vectmp3)*bytes_per_real
  end if
  write(6,'("| QMMM:          _REAL_ Scratch arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  element_memory = size(qmmm_scratch%qm_int_scratch)*bytes_per_int
#ifdef MPI
  !Add in the jrange arrays to the scratch count.
  element_memory = element_memory + size(qmmm_mpi%nquant_nlink_jrange)*bytes_per_int
#endif
  !Add in any scratch space for divide and conquer lapack diagonaliser
  element_memory = element_memory + qmmm_scratch%lapack_dc_int_scr_aloc*bytes_per_int
  if (qmmm_nml%allow_pseudo_diag) &
    element_memory = element_memory + size(qmmm_scratch%pdiag_vecjs)*bytes_per_int
  write(6,'("| QMMM:         Integer Scratch arrays : ",i12," bytes")') element_memory
  total_memory = total_memory + element_memory

  write(6,'("| QMMM: ---------------------------------------------------")')
  write(6,'("| QMMM:        Total Dynamic Memory Usage: ",f10.3," Mb")') &
            total_memory * bytes_to_mb

end subroutine qm_print_dyn_mem

!------------------------------------------------------------------------------

subroutine qm_print_ref(ref_index,atomic_number)
!This routine prints the reference corresponding to ref_index
!The idea here is that in the parameter list a paper index is
!given to each parameter set and the relevant reference can
!be printed at startup for each element in the calculation.
!
!Written by Ross Walker, TSRI, 2005
!
!================================================

  use qmmm_module, only : element_sym
  implicit none

!Passed in
      integer, intent(in) :: ref_index, atomic_number

!Local
      !Reference 1 - Most MNDO parameters
      !M.J.S. DEWAR, W. THIEL, JACS., 99, 4899, (1977)
      if (ref_index == 1) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 99, 4899, (1977)")') &
                      element_sym(atomic_number)
      else if (ref_index == 2) then
        write(6,'("| QMMM: ",A2,": TAKEN FROM MNDOC BY W.THIEL, QCPE 438, 2, p63, (1982)")') &
                      element_sym(atomic_number)
      else if (ref_index == 3) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 100, 777, (1978)")') &
                      element_sym(atomic_number)
      else if (ref_index == 4) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 99, 5231, (1977)")') &
                      element_sym(atomic_number)
      else if (ref_index == 5) then
        write(6,'("| QMMM: ",A2,": L.P.DAVIS et al. JCC, 2, 433, (1981)")') &
                      element_sym(atomic_number)
      else if (ref_index == 6) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS 5, 375 (1986)")') &
                      element_sym(atomic_number)
      else if (ref_index == 7) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 100, 3607, (1978)")') &
                      element_sym(atomic_number)
      else if (ref_index == 8) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 7, 140, (1986)")') &
                      element_sym(atomic_number)
      else if (ref_index == 9) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 4, 158, (1983)")') &
                      element_sym(atomic_number)
      else if (ref_index == 10) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 5, 1494, (1986)")') &
                      element_sym(atomic_number)
      else if (ref_index == 11) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 6, 186, (1987)")') &
                      element_sym(atomic_number)
      else if (ref_index == 12) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 4, 542, (1983)")') &
                      element_sym(atomic_number)
      else if (ref_index == 13) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS,106, 6771, (1984)")') &
                      element_sym(atomic_number)
      else if (ref_index == 14) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JCC, 5,358, (1984)")') &
                      element_sym(atomic_number)
      else if (ref_index == 15) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 4, 1964, (1985)")') &
                      element_sym(atomic_number)
      else if (ref_index == 16) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR, et al. ORGANOMETALLICS, 4, 1973, (1985)")') &
                      element_sym(atomic_number)
      else if (ref_index == 17) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. JACS, 107, 3902, (1985)")') &
                      element_sym(atomic_number)
      else if (ref_index == 18) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. THEOCHEM, 180, 1, (1988)")') &
                      element_sym(atomic_number)
      else if (ref_index == 19) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 9, 508, (1990)")') &
                      element_sym(atomic_number)
      else if (ref_index == 20) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 6, 1486, (1987)")') &
                      element_sym(atomic_number)
      else if (ref_index == 21) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. THEOCHEM, 187,1, (1989)")') &
                      element_sym(atomic_number)
      else if (ref_index == 22) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. INORG. CHEM., 29, 3881, (1990)")') &
                      element_sym(atomic_number)
      else if (ref_index == 23) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 7, 522, (1988)")') &
                      element_sym(atomic_number)
      else if (ref_index == 24) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 8, 1544, (1989)")') &
                      element_sym(atomic_number)
      else if (ref_index == 25) then
        write(6,'("| QMMM: ",A2,": M.J.S.DEWAR et al. ORGANOMETALLICS, 8, 1547, (1989)")') &
                      element_sym(atomic_number)
      else if (ref_index == 26) then
        write(6,'("| QMMM: ",A2,": J.J.P.STEWART, JCC, 10, 209 (1989)")') &
                      element_sym(atomic_number)
      else if (ref_index == 27) then
        write(6,'("| QMMM: ",A2,": J.J.P.STEWART, JCC, 3, 320 (1991)")') &
                      element_sym(atomic_number)
      else if (ref_index == 28) then
        write(6,'("| QMMM: ",A2,": J.P.MCNAMARA et al. CHEM. PHYS. LETT., 394, 429, (2004)")') &
                      element_sym(atomic_number)
      else if (ref_index == 29) then
        write(6,'("| QMMM: ",A2,": REPASKY et al. JCC, 23, 1601, (2002)")') &
                      element_sym(atomic_number)
      else if (ref_index == 30) then
        write(6,'("| QMMM: ",A2,": TUBERT-BROHMAN et al. JCC, 25, 138, (2003)")') &
                      element_sym(atomic_number)
      else if (ref_index == 31) then
        write(6,'("| QMMM: ",A2,": TUBERT-BROHMAN et al. JCTC, 1, 817, (2005)")') &
                      element_sym(atomic_number)
      else if (ref_index == 32) then
        write(6,'("| QMMM: ",A2,": J.P.MCNAMARA et al. J. MOL. GRA. MOD., 24, 128, (2005)")') &
                      element_sym(atomic_number)
      else if (ref_index == 33) then
        write(6,'("| QMMM: ",A2,": G.B.ROCHA et al. J. COMP. CHEM., 27, 1101, (2006)")') &
                      element_sym(atomic_number)
      else
        !UNKNOWN REFERENCE
        write(6,'(''| QMMM: '',A2,'': REFERENCE UNKNOWN.'')') element_sym(atomic_number)
      endif

      return

end subroutine qm_print_ref

!------------------------------------------------------------------------------
#ifdef QMMM_OMP
subroutine qm_print_omp_info()

   implicit none

   integer omp_get_num_threads !omp_get_num_threads is a function

!if OpenMP is going to be in use (to do the diagonalization) then the master thread
!when it gets to doing the diagonalisation will spawn shared memory threads to match
!the number of threads specified with the environment variable OMP_NUM_THREADS.
!Here we will print an informational message.
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
   write(6,'(/a)') '| QMMM(OpenMP): SMP Multithreading is in use for matrix diagonalization routines.'
   write(6,'(a,i4,a/)') '| QMMM(OpenMP): ',omp_get_num_threads(),' threads will be spawned by the master MPI thread as needed.'
!$OMP END MASTER
!$OMP END PARALLEL

   return

end subroutine qm_print_omp_info
#endif
