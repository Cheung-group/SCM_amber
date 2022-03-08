! <compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains allocatable arrays
!and variables used for coupled potential
!qmmm calculations. If you need access to one
!of the public variables from module you should
!include it in your code with
!use qmmm_module
!
!Main Authors:
!     Ross Walker
!     Mike Crowley
!+++++++++++++++++++++++++++++++++++++++++++++

module qmmm_module

use constants, only : A2_TO_BOHRS2, one, zero

!--------- CUTOFF PARAMETERS ---------------

! Exponential decay factor for the MM atom in the core-core interaction
_REAL_, parameter ::  ALPH_MM = 5.0d0

_REAL_, parameter :: AXIS_TOL = 1.0d-8  !Tolerance at which to define a vector is along the axis.
_REAL_, parameter :: OVERLAP_CUTOFF = 100.0d0*A2_TO_BOHRS2 !Distance^2 in bohrs at which to assume
                                                           !Gaussian overlap is zero.
_REAL_, parameter :: EXPONENTIAL_CUTOFF = 30.0d0 !Value of x at which to assume Exp(-x) = zero.


! ---- Common constants and definitions ----

integer, parameter :: PM3 = 1
integer, parameter :: AM1 = 2
integer, parameter :: MNDO = 3
integer, parameter :: PDDGPM3 = 4
integer, parameter :: PDDGMNDO = 5
integer, parameter :: PM3CARB1 = 6
integer, parameter :: DFTB = 7
integer, parameter :: RM1 = 8

!-------------------------------------------

!---------------- ELEMENTS -----------------

!Total number of elements
INTEGER, PARAMETER :: nelements = 86  !see also parameter in parameters.h
                                                                                                                                                             
! . The element symbols.
CHARACTER ( LEN = 2 ), DIMENSION(1:nelements), PARAMETER :: ELEMENT_SYM = (/ &
           'H ', 'He', &
           'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',  &
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar',  &
           'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  &
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', &
           'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
           'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', &
           'Cs', 'Ba', 'La', &
           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', &
           'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  &
           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'/)

!-------------------------------------------


!--------- SIZE / ARRAY LIMITATIONS -----------------------

      ! Locks maximum valence orbitals at 4 = S,P (no D or F - code is not present)
integer, parameter :: MAX_VALENCE_ORBITALS =4 
integer, parameter :: MAX_VALENCE_DIMENSION = MAX_VALENCE_ORBITALS*(MAX_VALENCE_ORBITALS+1)/2

!----------------------------------------------------------

!GLOBAL VARIABLES / ARRAYS - AVAILABLE TO ANY CODE THAT USES THIS MODULE

!REMEMBER TO MODIFY qmmm_mpi_setup below if you modify this structure
type qmmm_namelist
   !Contains QMMM namelist variables - all cpus should have this info available.
 _REAL_ :: qmcut     !Cutoff in angstroms for QM-MM electrostatics - default = same as MM-MM cutoff
 _REAL_ :: qmcut2    !Cutoff^2
 _REAL_ :: lnk_dis   !Distance in angstroms between QM atom and link atom
 _REAL_ :: scfconv    !SCF Convergence criteria for SCF routine - Minimum (tightest conv) = 1.0D-12, max = 1.0D0
                      !Default = 1.0D-8
 _REAL_ :: density_conv !Convergence criteria for density matrix. If tight_p_conv = true then = scfconv
                        !else it is 0.5*sqrt(scf_conv)
 _REAL_  :: pseudo_diag_criteria           !Criteria for the maximum change in the density matrix between two 
                                           !scf iterations before allowing pseudo diagonalisations.
                                           !Default = 0.05.
 integer :: lnk_atomic_no !Atomic number of link atom.
 integer :: qmgb      !Method for GB in QMMM (when igb /= 0 and /= 6) - 0 (default) leave QM charges at 0
                      ! effectively doing gas phase QM (EGB(MM-MM only). 1 - use resp charges from prmtop
                      ! for qmmm charges. 2 - use mulliken charges for QM atoms consistent with GB field, modified
                      ! fock matrix. 3 - (FOR_DEBUG_ ONLY) use mulliken charges from gas phase QMMM calc, not consistent
                      ! with GB field and so gradients will not be accurate.
 integer :: qmtheory  !Level of theory to use for QM region (Hamiltonian) 1 = PM3, 2 = AM1, 3 = MNDO (Default = PM3)
                      !                                                   4 = PDDG/PM3, 5 = PDDG/MNDO, 6 = PM3CARB1
                      !                                                   7 = DFTB,     8 = RM1
 integer :: qmcharge  !Charge of the QM region in electron units - must be an integer charge. (Default = 0)
 integer :: spin      !Spin state - default = 1 (singlet). Current Valid values = 1 (alternate options not currently available).
 integer :: verbosity !Controls amount of info printed about qm part of calc - (Default = 0)
 integer :: itrmax    !Maximum number of SCF cycles to conduct before assuming convergence has failed. (Default = 1000).
 integer :: qmshake   !Whether to shake qm atoms if ntc>1 (default = 1 - shake QM atoms) 0 = do not shake.
 integer :: kmaxqx, kmaxqy, kmaxqz !Used for qmewald - Maximum K space vectors
 integer :: ksqmaxq                !Used for qmewald - Maximum K squared values for spherical cutoff in k space.
 integer :: qm_ewald          !Flag for doing an ewald sum for periodic QM-QM interactions - default = 1
 integer :: qmmm_int  !Flag controlling the way QM-MM interactions are handled:
                      !   0 = mechanical embedding. No electrostatic QM-MM is calculated. Only interactions
                      !       come from VDW and bonds, angles and dihedrals that cross the boundary.
                      !   1 = Full QM-MM interaction is calculated by placing a 1S gaussian orbital on the MM
                      !       atom and calculating the full multipole interaction (Default)
                      !   2 = As 1 but also include extra Gaussian core-core terms if AM1, PM3 or derivative
                      !       type hamiltonians are in used.
 integer :: adjust_q  !Flag for whether to adjust the charges of certain MM atoms so that charge is conserved.
                      !   0 = No charge adjustment
                      !   1 = Missing charge is distributed over the nearest atom to each MM link pair.
                      !   2 = Missing charge is distributed evenly over all MM atoms - excluding MM link pairs. (default)

 integer :: diag_routine  !Flag controlling which diagonalization routine to use when doing full diag in SCF.
                          !  0 = Automatically pick fastest routine.
                          !  1 = Use internal diagonalization routine. (default)
                          !  2 = Use lapack dspev.
                          !  3 = Use lapack dspevd.
                          !  4 = Use lapack dspevx.
                          !  5 = Use lapack dsyev.
                          !  6 = Use lapack dsyevd.
                          !  7 = Use lapack dsyevr. (not currently implemented)
                          !  8 = Use lapack dsyevx. (not currently implemented)
 integer :: density_predict !Flag controlling the way in which the density matrix for MD step t+dt is predicted based
                            !on previous MD steps.
                            ! 0 = Use final converged density matrix from previous MD step. (default)
                            ! 1 = Use predictive algorithm based on Phys. Rev. Lett., 2006, 97, 123001
                            !     Here Pguess(t) = 2Pconv(t-dt) - Pguess(t-2dt)

 integer :: idc       !sah added for divcon/MOPAC
 integer :: divpb     !added for DivPB

 integer, dimension(:), pointer :: iqmatoms !integer list of atom numbers of the qm atoms as read from mdin (nquant)
                                            !nquant_nlink long, link atoms are on the end and contain the number of the
                                            !mm pair making up the link atom.
 
! GMS: Charge scaling factor for QM-FEP
 _REAL_  :: chg_lambda      ! Charge scaling for Free Energy calculation
!!DFTB Options
 integer :: dftb_maxiter    ! Maximum number of SCC iterations before resetting Broyden (default: 70 )
 integer :: dftb_disper     ! Use dispersion?  (default: 0 = false)
 integer :: dftb_chg        ! DFTB Charge output: (default: 0 = Mulliken, 1 = CM3)
 _REAL_  :: dftb_telec      ! Electronic temperature, in Kelvins. Default: 0.0K
 _REAL_  :: dftb_telec_step ! Step size for automatic telec increse, for convergence. Default = 0.0K
 !!End DFT Options

 logical :: ifqnt     !Flag for whether QMMM is on (True) or off (False)
 logical :: qmqm_analyt    !Flag for analytical vs (pseudo) numerical qm-qm derivates in qm2 routines. (Default = true, analytical).
 logical :: tight_p_conv !Flag to control convergence criteria for density matrix. If 0 SCF routine will converge energy
                         !to SCFCRT (def = 1*10^-8) and density to 0.05 * Sqrt(SCFCONV)
                         !If 1 then both energy and density will be converged to SCFCONV.
 logical :: printcharges !Flag to control the printing of mulliken, cm1a and cm2a charges. If set to true (1) then
                         !charges are calculated and printed on every step. Otherwise no charges are printed.
 logical :: peptide_corr !Default = .false. - add correction to peptide linkages?
 logical :: qmqm_erep_incore               !Store qmqm 1-electron repulsion integrals in memory?
 logical :: allow_pseudo_diag              !Whether or not to allow pseudo diagonalisations in SCF when possible.
 logical :: qmqmrij_incore   !Flag to store qmqm rij and related equations in memory - default = true.
 logical :: qmmmrij_incore   !Flag to store qmmm rij and related equations in memory - default = true.
 logical :: writepdb  !if true then a crude pdb of the qm system will be written to file qmmm_region.pdb
 logical :: qm_pme    !Flag for using PME instead of regular ewald for QM-MM interactions when qm_ewald>0.
                      !Default = true


end type qmmm_namelist

type ( qmmm_namelist ) qmmm_nml



type qmmm_structure
 _REAL_ :: enuclr_qmqm, enuclr_qmmm !Calculated Core-Core energy for QM-QM nuclei-nuclei and QM nuclei - MM charge interactions.
                                    !in (eV)
 _REAL_ :: elec_eng                 !Electronic energy (in eV)
 _REAL_, dimension(:), pointer :: qm_resp_charges  !Nquant long - contains the original resp charges for the
                                                       !QM atoms as read from the prmtop file.
 _REAL_ :: qm_resp_charge_sum !Sum of the resp charges making up the quantum region. - In AMBERELECTROSTATIC UNITS
 _REAL_, dimension(:), pointer :: mm_link_pair_resp_charges
 _REAL_, dimension(:,:), pointer :: mm_link_pair_saved_coords !The coordinates of the mm link pair atoms as they would be for
                                                            !a given MD step in amber's unimaged coordinate array. Extracted
                                                            !by position link atoms and restored back
                                                            !into the x array by restore_mm_link_pair_coords.
 _REAL_, dimension(:,:), pointer :: qm_coords     !Cartesian coordinates of ALL (real+link) qm atoms [3*(nquant+nlink) long]
 _REAL_, dimension(:), pointer :: scaled_mm_charges !MM charges scaled by one scale to make them electron units
 _REAL_, dimension(:,:), pointer :: dxyzqm, dxyzcl  !Used to store the forces generated by qm_mm before adding them to the main f array.
 _REAL_, dimension(:,:), pointer :: qm_xcrd         !Contains imaged mm coordinates and scaled mm charges.
 integer :: nquant    !Total number of quantum atoms (excluding link atoms) as read from mdin
 integer :: nlink     !Total number of link atoms
 integer :: nquant_nlink !Total number of quantum atoms = nquant+nlink
 integer :: qm_ntypes                           !The number of atom types present. 
 integer, dimension(nelements) :: qm_type_id    !The id of each type, essentially the atomic number of that
                                                !type, e.g. type 1 may have atomic number = 8.
 integer, dimension(:), pointer :: qm_atom_type !The type of each qm atom, essentially a re-basing of the
                                                !atomic numbers to minimise memory usage.

 integer, dimension(:,:), pointer :: link_pairs ! list of MM-QM atoms for which link atoms were added
 integer, dimension(:), pointer :: iqm_atomic_numbers !integer list of atomic numbers for the qm atoms 
 integer                        :: qm_mm_pairs     !Number of pairs per QM atom. - length of pair_list. 
 integer, dimension(:), pointer :: qm_mm_pair_list !Non bond pair list for each QM atom
 integer, dimension(:), pointer :: all_atom_numbers !atomic numbers of ALL atoms (MM and QM)
 integer :: ntotatm
 integer, dimension(:), pointer :: jqatms !Array for dcqtp to keep track of QM atoms -- allocated in qm_div.f
 integer :: num_qmmm_calls       !Counter of the number of times qm_mm has been called - effectively nstep.
 logical, dimension(:), pointer :: atom_mask !True / false mask specifying if atom is a QM atom. True = QM atom (natom long)
 logical, dimension(:), pointer :: mm_link_mask !True / false mask specifying if atom is a MM link pair atom. True = MM link pair atom (natom long)
 logical :: qm_mm_first_call !Set to true at beginning of sander subroutine and then set to false at the end of qm_mm. Used for allocation purposes.
 logical :: fock_first_call !Set to true at beginning of sander subroutine and then set to false at the end of qm_mm. Used for allocation purposes.
 logical :: fock2_2atm_first_call 
 logical :: qm2_allocate_e_repul_first_call
 logical :: qm2_calc_rij_eqns_first_call
 logical :: qm2_scf_first_call
 logical :: zero_link_charges_first_call
 logical :: adj_mm_link_pair_crd_first_call
 logical :: AM1_OR_PM3 !Set to True if theory is AM1, PM3, PDDG/PM3, PM3CARB1, RM1
 logical :: PDDG_IN_USE !Set to True if theory is PDDG/PM3 or PDDG/MNDO
end type qmmm_structure

type ( qmmm_structure ) qmmm_struct

type qm2_structure  !Variables that are specific to qm_routine=2 (qm2)
 _REAL_, dimension(:), pointer :: den_matrix      !The total density matrix
 _REAL_, dimension(:), pointer :: old_den_matrix  !The old total density matrix from preious step
                                                  !Allocated in qm2_load_params on first call - deallocated
                                                  !by deallocate_qmmm
 _REAL_, dimension(:), pointer :: old2_density  !Used by qm2_cnvg as workspace, norbs

 _REAL_, dimension(:), pointer :: md_den_mat_guess1 !These two guesses are only used when density_predict=1
 _REAL_, dimension(:), pointer :: md_den_mat_guess2 !They contain Pguess(t-1) and Pguess(t-2).

 _REAL_, dimension(:), pointer :: fock_matrix     !Fock matrix
 _REAL_, dimension(:,:), pointer :: qm_mm_e_repul !Array containing the QM-MM electron repulsion integrals
 _REAL_, dimension(:), pointer :: qm_qm_2e_repul  !Array containing the QM-QM 2-electron repulsion integrals
                                                  !This is a big array, allocated by qm_mm and deallocated
                                                  !by deallocate_qmmm - it needs to be a total of n2el long.
 _REAL_, dimension(:), pointer :: hmatrix  !The 1-electron matrix - used by routines called from qm2_energy
                                           !allocated in qm_mm on first call - deallocated by  deallocate_qmmm
 _REAL_, dimension(:,:), pointer :: qm_qm_e_repul !Array containing the QM-QM electron repulsion integrals
                                                  !This was originally written as a file to disk in the energy
                                                  !routines and then re-read in the derivative routines. Now
                                                  !it is stored in memory. Allocated in qm_mm on first call
                                                  ! - deallocated by deallocate_qmmm
 _REAL_, dimension(:,:), pointer :: fock2_ptot2   !Used in qm2_fock2 routine - 16,nquant_nlink allocated in qm2_load_params
                                                  !deallocated in deallocate_qmmm
 _REAL_, dimension(:,:), pointer :: eigen_vectors !Holds the eigen vectors during the SCF - allocated in qm2_load_params
 _REAL_, dimension(:), pointer :: scf_mchg !Hold the mulliken charges at each scf step if calc_mchg_scf is true.
                                           !Otherwise only do it at the end of the scf.
 integer :: matsize                        !Size of the various packed symmetric matrices. (norbs(norbs+1)/2)
 integer :: n2el                           !Number of 2 electron repulsion integrals, calculated by
                                           !moldat = 
                                           !50*nheavy(nheavy-1)+10*nheavy*nlight+(nlight*(nlight-1))/2
 integer :: norbs                          !Total Number of atomic orbitals.
 integer :: nclosed                        !Number of doubly occupied orbitals
 integer :: nopenclosed                    !Number of doubly occupied and singly occupied orbitals
 integer :: qm_mm_e_repul_allocated        !This was originally written as a file to disk in the energy
                                           !routines and then re-read in the derivative routines. Now
                                           !it is stored in memory. Allocated in qm_mm on first call
                                           !or if pair list changes too much. See qm_mm_e_repul array above
                                           ! - deallocated by deallocate_qmmm
 integer                         :: n_peptide_links !Number of peptide linkages in QM region to apply MM correction to.
 integer, dimension(:,:), pointer :: peptide_links !Identity of peptide linkages 4,n_peptide_linkages. 1 to 4 = H-N-C-O atom
                                                   !numbers.
 logical calc_mchg_scf                     !If set to true the mulliken charges will be calculated on each SCF iteration.
end type qm2_structure

type ( qm2_structure ) qm2_struct

type qm2_params_structure !Parameters for each atom in the qm region
 _REAL_ :: tot_heat_form !Calculated in qm2_load_params - independent of structure so constant for a sander run.
 _REAL_, dimension(:), pointer :: core_chg !The core charge on each atom as seen by the electrons
                                           !allocated as nquant_nlink long by qm2_load_params - deallocated
                                           !by deallocate_qmmm
 _REAL_, dimension(:,:), pointer :: orb_elec_ke !Orbital electron kinetic energy integrals (2,nquant_nlink) - allocated in
                                              !qm2_load_params. Deallocated by deallocate qmmm (1=s, 2=p
 _REAL_, dimension(:,:), pointer :: betasas !betas(ni)+betas(nj) qm_ntypes x qm_ntypes
 _REAL_, dimension(:,:), pointer :: betasap !betas(ni)+betap(nj) qm_ntypes x qm_ntypes
 _REAL_, dimension(:,:), pointer :: betapap !betap(ni)+betap(nj) qm_ntypes x qm_ntypes
 _REAL_, dimension(:,:), pointer :: FN1     !PM3 / AM1 specific parameters for core-core repulsions.
 _REAL_, dimension(:,:), pointer :: FN2     !PM3 / AM1 specific parameters for core-core repulsions.
 _REAL_, dimension(:,:), pointer :: FN3     !PM3 / AM1 specific parameters for core-core repulsions.
 _REAL_, dimension(:,:), pointer :: onec2elec_params
                                            !Coulomb and exchange one centre-two electron integral params.
                                            !Allocated as nquant_nlink long in qm2_load_params. Deallocated
                                            !in deallocate QMMM.
 _REAL_, dimension(:,:), pointer :: multip_2c_elec_params !Parameters for the multipole expansion of the 2 centre 2
                                                   !electron integerals. 9, nquant_nlink in order
                                                   !DD,QQ,AM,AD,AQ,AM2,AD2,AQ2 allocated in qm2_load_params.
 _REAL_, dimension(:), pointer :: cc_exp_params    !Exponents for core core repulsions.

 _REAL_, dimension(:), pointer :: s_orb_exp_by_type !S orbital expansion coefficients for the Slater orbital expansion
 _REAL_, dimension(:), pointer :: p_orb_exp_by_type !P orbital expansion coefficients for the Slater orbital expansion
                                            !both are allocated as qm_ntypes long by qm2_load_params
                                            !deallocated by qm2_setup_orb_exp as it is not needed after this
                                            !is called.
!Arrays for PDDG Hamiltonians
 _REAL_, dimension(:), pointer :: pddge1, pddge2
!Arrays for pre-computed orbital interactions
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_sxs_over_sas !atom_orb_zz_s_x_s/atom_orb_zz_one_s_a_s
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_sxp_over_sap !atom_orb_zz_s_x_p/atom_orb_zz_one_s_a_p
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_zz_pxp_over_pap !atom_orb_zz_p_x_p/atom_orb_zz_one_p_a_p
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_ss_eqn !sqrt((two*sqrt(atom_orb_zz_s_x_s)*atom_orb_zz_one_s_a_s)**3
                                                        !*atom_orb_cc_s_x_s
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_ovlp!2.0D0 * qm2_params%atom_orb_zz_s_by_type(K,qmitype)*
                                                        !SQRT(qm2_params%atom_orb_zz_p_by_type(L,qmjtype))*
                                                        !atom_orb_zz_one_s_a_p(k,l,qmitype,qmjtype)
                                                        !*atom_orb_sp_eqn
                                                        !Used in gover for Si-Pj overlap energy and -Sj-Pi overlap.

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_inj !-4.0D0*sqrt(atom_orb_zz_p_x_p)*
             !atom_orb_zz_one_p_a_p*qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)*
             !atom_orb_pp_eqn
             !Used in gover for Pi-Pj overlap energy when i!=j

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_ieqj1 !-4.0D0*atom_orb_pp_eqn*
             !sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p*
             !qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)
             !Used in gover for Pi-Pj overlap energy when i==j

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_ovlp_ieqj2 !2.0D0*atom_orb_pp_eqn*
             !sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_ss_eqn_adb !-2.0D0*A2_TO_BOHRS2*qm2_params%atom_orb_ss_eqn_adb(i,j,qmitype,qmjtype)*
                                                            !qm2_params%atom_orb_zz_sxs_over_sas(i,j,qmitype,qmjtype)
                                                            !Used for S-S overlap in QM-QM derivatives

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xy  !-four*A3_TO_BOHRS3*qm2_params%atom_orb_zz_sxp_over_sap(i,j,qmitype,qmjtype)**2* &
              !(one/(SQRT(qm2_params%atom_orb_zz_p_by_type(J,qmjtype))))*atom_orb_sp_eqn
              !Used for S-P overlap in QM-QM derivatives where P... /= axis
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xx1 !Used for S-P overlap in QM-QM derivatives where P... == axis
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_sp_eqn_xx2 

 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_eqn_xxy1 !-four*A2_TO_BOHRS2*(ADB_array(inner_index)**2)*
                                                             !(one/(SQRT(atom_orb_zz_p_x_p)))*atom_orb_pp_eqn
                                                             !Used for P-P overlap in QM-QM derivatives where P... = P... /= axis
 _REAL_, dimension(:,:,:,:), pointer :: atom_orb_pp_eqn_xxy2 !eight*A2_TO_BOHRS2*A2_TO_BOHRS2*(ADB_array(inner_index)**3)*
                                                             !(one/(SQRT(atom_orb_zz_p_x_p)))*atom_orb_pp_eqn
 
!End arrays for pre-computed orbital interactions

!Pre-computed PDDG parameters if PDDG Hamiltonian is in use.
 _REAL_, dimension(:,:), pointer :: pddg_term1, pddg_term2, pddg_term3, pddg_term4
                                    !Ntypes*ntypes - allocated in qm2_load_params if PDDG is in use. Stores the pre-exponential
                                    !part of the PDDG equation.
 integer, dimension(:), pointer :: natomic_orbs !Number of atomic orbitals on atom.
 integer, dimension(:,:), pointer :: orb_loc    !Locations of orbitals. 2,nquant_nlink. 1,x gives beginning of
                                                !of orbitals on atom x. 2,x gives last orbital on atom x.
 integer, dimension(:), pointer :: pascal_tri1 !Lower half triangle indices (pascal's triangle)
 integer, dimension(:), pointer :: pascal_tri2 !Allocated in load params.
 integer, dimension(:), pointer :: NUM_FN   !Number of FNX terms (first dimension) that are not zero.
end type qm2_params_structure

type (qm2_params_structure) qm2_params

type  qm2_rij_eqns_structure !This structure is used to store RIJ info for each QM-QM pair and related equations
!QM-QM
  _REAL_, dimension(:,:), pointer :: qmqmrijdata !These two arrays store QMQM and QMMM RIJ values and related
!QM-MM                                           !equations. See array_locations.h for the first dimension
 _REAL_, dimension(:,:), pointer ::  qmmmrijdata !offsets.
 integer :: qmmmrij_allocated
end type qm2_rij_eqns_structure

type (qm2_rij_eqns_structure) qm2_rij_eqns

!QMEwald specific structure
type qm_ewald_structure
   _REAL_, dimension(:), pointer :: kvec !Array storing K vectors (totkq long
   _REAL_, dimension(:,:), pointer :: dkvec !Array storing derivative K vectors (3,totkq long
   _REAL_, dimension(:,:), pointer :: dmkv !used for calculating the ktable
   _REAL_, dimension(:,:,:), pointer :: ktable !Table for storing complex exp(ik,r[j])
                                               !dimensions = 6,natom,totkq
                                               !1,x,y = x_cos
                                               !2,x,y = x_sin
                                               !3,x,y = y_cos
                                               !4,x,y = y_sin
                                               !5,x,y = z_cos
                                               !6,x,y = z_sin
   _REAL_, dimension(:,:,:), pointer :: qmktable !As Ktable but stores the qmatom copies in a linear 1->nquant fashion.
   _REAL_, dimension(:), pointer :: mmpot !Nquant long, stores the potential at each QM atom due to the MM field.
   _REAL_, dimension(:), pointer :: qmpot  !Nquant long, stores the self energy of the QM atoms to avoid double counting
   _REAL_, dimension(:,:), pointer :: d_ewald_mm !3,natom long stores gradients on MM atoms due to QM-MM Ewald field.
                                                 !Reciprocal forces.
   _REAL_ :: ewald_core !Ewald Potential with QM CORE charges - energy in eV.
   _REAL_ :: mm_recip_e !Reciprocal energy from MM atoms - qm_pme.
   integer :: totkq  !Total number of kspace vectors
   integer :: natom  !Same as sander's natom, copied here by qm_mm for convenience.
   logical :: ewald_startup !True if this is the very first MD step and we are doing qmewald.
end type qm_ewald_structure

type (qm_ewald_structure) qmewald

type qm_gb_structure
   _REAL_, dimension(:), pointer :: qmqm_onefij !Stores the 1.0/fij equations for qm-qm pairs. Since these
                                             !values only depend on Rij and the effective radii 
                                             !they remain fixed during the SCF and so are calculated
                                             !once outside the SCF and then reused inside.
   _REAL_, dimension(:), pointer :: qmqm_kappafij !Stores exp(-kappa*fij) - These are calculated outside of the
                                                  !scf and reused inside. Only allocated and calculated if kappa/=0.0d0
                                                  !in other words saltcon /= 0.0d0.
   _REAL_, dimension(:), pointer :: gb_mmpot  !Nquant long, stores GB potential at each QM atom due to MM atoms.
   _REAL_, dimension(:), pointer :: gb_qmpot  !Nquant long, stores GB potential at each QM atom due to QM atoms.
   _REAL_ :: intdieli, extdieli !1.0d0/intdiel and extdiel respectively
   _REAL_ :: kappa    !Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
                      !T = 298.15, epsext=78.5, kappa = sqrt( 0.10806d0 * saltcon )
                      !scaled kappa by 0.73 to account(?) for lack of ion exlcusions. 
   _REAL_ :: mmcut2   !cut^2 in angstroms^2 from cntrl namelist
   _REAL_ :: one_Arad_beta !alpb_beta/Arad when alpb/=0
   integer, dimension(:,:), pointer :: qmqm_gb_list !1+nquant,nquant - list of qm atoms that interact with current qm atom.
                                                    !+1 because the first entry stores the number of interactions for QM atom y. 
   logical :: saltcon_on !True if saltcon /= 0.0d0
   logical :: alpb_on !True if alpb = 1
end type qm_gb_structure

type (qm_gb_structure) qm_gb

type qmmm_mpi_structure
  integer :: commqmmm !Communications within a given set of QMMM threads potentially a subset of processors of commsander.
  integer :: numthreads !Number of threads in commqmmm.
  integer :: mytaskid   !Task id of this thread in commqmmm.
  integer :: openmp_numthreads !When openmp is in use this will contain how many threads would be spawned by the master.
  integer :: natom_start !Where this thread should start for 1,natom loops
  integer :: natom_end   !Where this thread should end for 1,natom loops
  integer :: nquant_nlink_start !Where this thread should start for 1,nquant_nlink loops
  integer :: nquant_nlink_end !Where this thread should end for 1,nquant_nlink loops
  integer :: totkq_count
  integer :: kvec_start !Where this thread starts and ends in 1 to totkq loops.
  integer :: kvec_end
  integer :: two_e_offset !Offset into two electron matrix for this thread

!Below are for matrix type double loop load balancing
  integer ::                        nquant_nlink_istart !These three control load balancing
  integer ::                        nquant_nlink_iend
  integer :: nquant_nlink_loop_extent_begin !
  integer :: nquant_nlink_loop_extent_end !
  integer, dimension(:,:), pointer :: nquant_nlink_jrange !within a do i = 1, nquant_nlink
                                                        !           do j=1,i-1
                                                        !loop. Basically istart and iend define
                                                        !the extent of the outer loop and then
                                                        !(1,i) and (2,i) of jrange define the 
                                                        !limits of the inner loop for this processor.
  logical :: commqmmm_master !True if master thread of commqmmm

end type qmmm_mpi_structure

type (qmmm_mpi_structure) qmmm_mpi

type qmmm_scratch_structure
!Various scratch arrays used as part of QMMM - one should typically assume that upon leaving a routine
!the contents of these arrays can be assumed to be junk.
  _REAL_, dimension(:), pointer  :: matsize_red_scratch !Allocated as qm2_struct%matsize when doing MPI during the 
                                                        !load parameters routine. ONLY ALLOCATED IF WE CAN'T
                                                        !DO MPI_IN_PLACE.
                                                        !+1 in size so we can pack extra energies on the end etc.
  _REAL_, dimension(:), pointer :: qm_pme_scratch !natom long scratch for qm_pme - only allocated if doing qm_pme
  _REAL_, dimension(:,:), pointer :: mat_diag_workspace !Matrix diagonalisation workspace - allocated in qm2_load_params
                                                        !norbs,6 for internal diag
                                                        !norbs,1 if lapack diag.
!The Pseudo diagonalizer needs a total of 5 real scratch arrays and 1 integer scratch array.
!Passed in Scratch Arrays - these are all allocated in qm2_load_params and only if allow_pseudo_diag is true.
!ONLY ALLOCATED ON COMMQMMM MASTER THREAD
  _REAL_, dimension(:,:), pointer :: pdiag_scr_norbs_norbs     !(norbs,norbs)
  _REAL_, dimension(:,:), pointer :: pdiag_scr_noccupied_norbs !(noccupied,norbs)
  _REAL_, dimension(:), pointer :: pdiag_vectmp1               !(noccupied*(norbs-noccupied))
  _REAL_, dimension(:), pointer :: pdiag_vectmp2               !(noccupied*(norbs-noccupied))
  _REAL_, dimension(:), pointer :: pdiag_vectmp3               !(noccupied*(norbs-noccupied))
  integer, dimension(:,:), pointer :: pdiag_vecjs                !(2,noccupied*(norbs-noccupied))

 _REAL_, dimension(:), pointer :: lapack_dc_real_scr
 _REAL_, dimension(:), pointer :: lapack_dc_int_scr
!END ONLY ALLOCATED ON COMMQMMM MASTER THREAD
!Scratch Arrays - These arrays are available for any subroutine to use - they are allocated by allocate_qmmm
!Each routine that uses one of these arrays for temporary storage should assume that the contents of such array
!will be garbage once the routine is left or another routine is called.
 _REAL_, dimension(:), pointer :: qm_real_scratch !Real scratch array - 4*natom.
 integer, dimension(:), pointer :: qm_int_scratch !Integer scratch array - 3*natom.
 integer :: lapack_dc_real_scr_aloc !Number of reals allocated for lapack_dc_real_scr
 integer :: lapack_dc_int_scr_aloc !Number of ints allocated for lapack_dc_int_scr_aloc
 integer :: qm_mm_pairs_allocated !Size of expected qm_mm_pairs that scratch arrays and calc_rij_array was allocated to.
end type qmmm_scratch_structure

type (qmmm_scratch_structure) qmmm_scratch

!END GLOBALS

!BEGIN SUBROUTINES

contains


#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Does QMMM specific mpi setup and broadcasts
subroutine qmmm_mpi_setup( master )
  implicit none

#  include "parallel.h"
#  include "mpif.h"
#  include "memory.h"
!Passed in
   logical, intent(in) :: master
!Local
   integer :: mpi_division, i, istartend(2)
   integer :: loop_extent, loop_extent_begin, loop_extent_end
   integer :: j, jstart, jend
   integer :: ier=0, istatus=0
   integer :: natom_replicates, nquant_replicates, nquant_nlink_replicates
   integer :: nlink_replicates

   !Broadcast all of the info in qmmm_nml
   !I don't think we can assume the structure is linear in memory so we have to broadcast each seperately
   call mpi_bcast(qmmm_nml%qmcut, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%qmcut2, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%lnk_dis, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%scfconv, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%density_conv, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%pseudo_diag_criteria, 1, MPI_DOUBLE_PRECISION, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%lnk_atomic_no, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmgb, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmtheory, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmcharge, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%spin, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%verbosity, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%itrmax, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%adjust_q, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%diag_routine, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%density_predict, 1, mpi_integer, 0, commsander, ier) 

   call mpi_bcast(qmmm_nml%idc, 1, mpi_integer, 0, commsander, ier) 

   call mpi_bcast(qmmm_nml%qmshake, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qm_ewald, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%kmaxqx, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%kmaxqy, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%kmaxqz, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%ksqmaxq, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmmm_int, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmqm_analyt, 1, mpi_logical, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%tight_p_conv, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%printcharges, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%peptide_corr, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%allow_pseudo_diag, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmqm_erep_incore, 1,mpi_logical, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%qmqmrij_incore, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qmmmrij_incore, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%qm_pme, 1, mpi_logical, 0, commsander, ier) 

   call mpi_bcast(qmmm_nml%writepdb, 1, mpi_logical, 0, commsander, ier) 

   call mpi_bcast(qmmm_struct%nquant, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_struct%nlink, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_struct%nquant_nlink, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_struct%qm_ntypes, 1, mpi_integer, 0, commsander, ier) 
   call mpi_bcast(qmmm_struct%qm_type_id, nelements, mpi_integer, 0, commsander, ier) 

   call mpi_bcast(qmmm_struct%AM1_OR_PM3, 1, mpi_logical, 0, commsander, ier) 
   call mpi_bcast(qmmm_struct%PDDG_IN_USE, 1, mpi_logical, 0, commsander, ier) 

   call mpi_bcast(qm2_struct%calc_mchg_scf,1, mpi_logical, 0, commsander, ier)

   call mpi_bcast(qm_gb%alpb_on,1,mpi_logical, 0, commsander, ier)
   call mpi_bcast(qm_gb%saltcon_on,1,mpi_logical, 0, commsander, ier)
   call mpi_bcast(qm_gb%kappa,1,MPI_DOUBLE_PRECISION, 0, commsander, ier)

!! GMS
   call mpi_bcast(qmmm_nml%chg_lambda , 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
!! DFTB options
   call mpi_bcast(qmmm_nml%dftb_maxiter    , 1, mpi_integer         , 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%dftb_disper     , 1, mpi_integer         , 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%dftb_chg        , 1, mpi_integer         , 0, commsander, ier) 
   call mpi_bcast(qmmm_nml%dftb_telec      , 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%dftb_telec_step , 1, MPI_DOUBLE_PRECISION, 0, commsander, ier)
!!

   !Now we do the arrays, note we have to allocate on all the non-master nodes
   !Note: we don't need to allocate the pair list as this is done in qm_mm routine.
   !It is also filled here so we don't need to broadcast it at this point.

   if ( .not. master ) then
      call allocate_qmmm( natom )
      if (qmmm_struct%nlink > 0) then
        allocate ( qmmm_struct%link_pairs(2,qmmm_struct%nlink), stat=ier )   
        REQUIRE(ier == 0)
      end if
      allocate ( qmmm_struct%qm_atom_type(qmmm_struct%nquant_nlink),stat=ier )
      REQUIRE(ier == 0)
   end if

!scf_mchg has not been allocated yet on non-master threads so allocate it.
   allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0) !Deallocated in deallocate qmmm

   !Now we can broadcast each of the arrays
   call mpi_bcast(qmmm_struct%qm_atom_type, qmmm_struct%nquant_nlink, MPI_INTEGER, 0, commsander, ier)
!All nodes now zero the charges and copy them.
!   call mpi_bcast(qmmm_struct%qm_resp_charges, qmmm_struct%nquant, MPI_DOUBLE_PRECISION, 0, commsander, ier)
   call mpi_bcast(qmmm_nml%iqmatoms, qmmm_struct%nquant_nlink, mpi_integer, 0, commsander, ier)
   if (qmmm_struct%nlink > 0) then
      call mpi_bcast(qmmm_struct%link_pairs, 2*qmmm_struct%nlink, mpi_integer, 0, commsander, ier)
   end if
   call mpi_bcast(qmmm_struct%iqm_atomic_numbers, qmmm_struct%nquant_nlink, mpi_integer, 0, commsander, ier)
   call mpi_bcast(qmmm_struct%atom_mask, natom, mpi_logical, 0, commsander, ier)
   call mpi_bcast(qmmm_struct%mm_link_mask, natom, mpi_logical, 0, commsander, ier)
   call mpi_bcast(qmmm_struct%scaled_mm_charges, natom, MPI_DOUBLE_PRECISION, 0, commsander, ier)

!Set the master flag on all threads
   qmmm_mpi%commqmmm_master = master

   natom_replicates = natom
   nquant_replicates = qmmm_struct%nquant
   nquant_nlink_replicates = qmmm_struct%nquant_nlink
   nlink_replicates = qmmm_struct%nlink 
!Setup the commqmmm communicator. For regular QMMM MD simulations this will
!be the same size as commsander and have the same members. For QMMM LES
!simulations this will contain only the current processor and be of numtasks=1

!In theory we could use this later to make a smaller commqmmm than commsander for efficiency reasons.

   !A single commqmm that is the same on all threads that make up commsander
   qmmm_mpi%commqmmm = commsander
   qmmm_mpi%numthreads = sandersize
   qmmm_mpi%mytaskid = sanderrank
   qmmm_mpi%commqmmm_master = master

   !Divide up i=2,nquant_nlink
   !             j=1,i-1
   !loops

   !Matrix looks like this. E.g. for 4 cpus and 12 QM atoms - thread id's listed
   !Total elements = 66
   !cpus 0,1 and 2 do 17 elements each
   !cpu 3 does 15 elements
   ! i 1 2 3 4 5 6 7 8 9 10 11 12
   !j1   0 0 0 0 0 0 1 1 2  2  3
   ! 2     0 0 0 0 0 1 1 2  2  3
   ! 3       0 0 0 1 1 1 2  2  3
   ! 4         0 0 1 1 1 2  2  3
   ! 5           0 1 1 1 2  2  3
   ! 6             1 1 1 2  2  3
   ! 7               1 2 2  3  3
   ! 8                 2 2  3  3
   ! 9                   2  3  3
   !10                      3  3
   !11                         3
   !12 

   !first of all work out out the total extent of the loop.
   !allocate the memory for the jrange array

   allocate(qmmm_mpi%nquant_nlink_jrange(2,nquant_nlink_replicates),stat=ier)
   REQUIRE(ier==0)
   loop_extent = nquant_nlink_replicates*(nquant_nlink_replicates-1)/2
   mpi_division = (loop_extent+(qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
   loop_extent_end = min(mpi_division*(qmmm_mpi%mytaskid+1),loop_extent)
   loop_extent_begin = mpi_division*qmmm_mpi%mytaskid+1
   !loop_extent_begin = (istart-1)(istart-2)/2 + jstart
   !loop_extent_end = (iend-1)(iend-2)/2 + jend
   !s = 1+sqrt(1+8x)/2
   !i = int(s) - ROUNDED UP
   qmmm_mpi%nquant_nlink_istart = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_begin)))/2.0d0)
   qmmm_mpi%nquant_nlink_iend   = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_end)))/2.d0)
   qmmm_mpi%nquant_nlink_loop_extent_begin = loop_extent_begin
   qmmm_mpi%nquant_nlink_loop_extent_end = loop_extent_end

   !Now we need to work out what range of j values we do for each i we will be doing.
   !What value of j would, when coupled with our istart give us loop_extent_begin?
   ! j = loop_extent_begin -((-i-1)(i-2)/2)

   jstart = loop_extent_begin - ((qmmm_mpi%nquant_nlink_istart-1)*(qmmm_mpi%nquant_nlink_istart-2)/2)
   jend   = loop_extent_end - ((qmmm_mpi%nquant_nlink_iend-1)*(qmmm_mpi%nquant_nlink_iend-2)/2)

   do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend

     if (i == qmmm_mpi%nquant_nlink_istart) then
       qmmm_mpi%nquant_nlink_jrange(1,i) = jstart
     else
       qmmm_mpi%nquant_nlink_jrange(1,i) = 1
     end if

     if (i == qmmm_mpi%nquant_nlink_iend) then
       qmmm_mpi%nquant_nlink_jrange(2,i) = jend
     else
       qmmm_mpi%nquant_nlink_jrange(2,i) = i-1
     end if

   end do

   !------- End Matrix Type Calc Load Balancing ------------

   !Now divide up the atoms between threads
   !We will divide up evenly as best we can. E.g. with 1603 atoms on 4 threads
   !we would ideally do 1->401, 402->802, 803->1203, 1204->1603
   !But this is quite difficult to manage. For the moment I will just divide
   !up using a simple formula that ensures that most threads have an even spread.
   !Thus threads 1 to 3 do 401 atoms and thread 4 does 400 atoms. This approach
   !allows linear movement in memory.
   !Note the last thread is responsible for ALL link atoms.

   !NATOM division first.
   mpi_division = (natom_replicates + (qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
   qmmm_mpi%natom_end = min(mpi_division*(qmmm_mpi%mytaskid+1),natom_replicates)
   qmmm_mpi%natom_start = min(mpi_division*qmmm_mpi%mytaskid+1,natom_replicates+1)

   !write info about atom division - Get each thread to send the master it's values
   !This allows a sanity check.
   if (qmmm_mpi%commqmmm_master) then
      write (6,'(/a,i4,a)') '|QMMM: Running QMMM calculation in parallel mode on ',qmmm_mpi%numthreads,' threads.'
      write (6,'(a)') '|QMMM: All atom division among threads:'
      write (6,'(a)') '|QMMM:                  Start       End      Count'
      !Already know my own.
      write(6,'(a,i8,a,i8,a,i8,a)') &
            '|QMMM: Thread(   0): ',qmmm_mpi%natom_start,'->',qmmm_mpi%natom_end, &
                                  '  (',qmmm_mpi%natom_end-qmmm_mpi%natom_start+1,')'
      do i = 1, qmmm_mpi%numthreads-1
         call mpi_recv(istartend,2,mpi_integer,i,0,qmmm_mpi%commqmmm,istatus,ier)
         write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
             '|QMMM: Thread(',i,'): ',istartend(1),'->',istartend(2), &
                                  '  (',istartend(2)-istartend(1)+1,')'
     end do
   else
     !Send a message to the master with our counts in.
     istartend(1) = qmmm_mpi%natom_start
     istartend(2) = qmmm_mpi%natom_end
     call mpi_send(istartend,2,mpi_integer,0,0,qmmm_mpi%commqmmm,ier)
   end if

   !Nquant
   mpi_division = (nquant_replicates + (qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
   qmmm_mpi%nquant_nlink_end = min(mpi_division*(qmmm_mpi%mytaskid+1),nquant_replicates)
   qmmm_mpi%nquant_nlink_start = min(mpi_division*qmmm_mpi%mytaskid+1,nquant_replicates+1)

   !Nquant_nlink - Last thread gets all the link atoms
   if (qmmm_mpi%mytaskid == qmmm_mpi%numthreads - 1) then
     !Last thread
     qmmm_mpi%nquant_nlink_end = qmmm_mpi%nquant_nlink_end + nlink_replicates
   end if

   if (qmmm_mpi%commqmmm_master) then
     write (6,'(/a)') '|QMMM: Quantum atom + link atom division among threads:'
     write (6,'(a)') '|QMMM:                  Start       End      Count'
     !Already know my own
     write(6,'(a,i8,a,i8,a,i8,a)') &
           '|QMMM: Thread(   0): ',qmmm_mpi%nquant_nlink_start,'->',qmmm_mpi%nquant_nlink_end, &
           '  (',qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1,')'
     do i = 1, qmmm_mpi%numthreads-1
       call mpi_recv(istartend,2,mpi_integer,i,0,qmmm_mpi%commqmmm,istatus,ier)
       write(6,'(a,i4,a,i8,a,i8,a,i8,a)') &
           '|QMMM: Thread(',i,'): ',istartend(1),'->',istartend(2), &
                                  '  (',istartend(2)-istartend(1)+1,')'
     end do
   else
     !Send a message to the master with our counts in.
     istartend(1) = qmmm_mpi%nquant_nlink_start
     istartend(2) = qmmm_mpi%nquant_nlink_end
     call mpi_send(istartend,2,mpi_integer,0,0,qmmm_mpi%commqmmm,ier)
   end if

!Now we need to calculate the offset into the 2 electron matrix that
!this thread will have - this depends on the number of light and heavy
!atoms and so is calculated at the end of qm2_load_params_and_allocate.

   return

end subroutine qmmm_mpi_setup
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates space for qmmm variables and arrays that depend only on nquant or natom
subroutine allocate_qmmm( natm )
   implicit none
!passed in
   integer , intent(in) :: natm
!Local
#  include "memory.h"
   integer :: ier=0
   integer :: natom_replicate, nquant_replicate

   allocate ( qmmm_struct%qm_resp_charges(qmmm_struct%nquant), stat=ier )
   REQUIRE(ier == 0)
 
   !iqmatoms and iqm_atomic_numbers are allocated here as nquant+nlink long but 
   !initially nlink is zero so it only gets allocated as nquant long on the master
   !thread. It is then resized when nlink are known about. When all other mpi threads
   !call this allocation routine nquant+nlink should include both qm and link atoms.
   allocate ( qmmm_nml%iqmatoms(qmmm_struct%nquant+qmmm_struct%nlink), stat=ier )
   REQUIRE(ier == 0)
   allocate ( qmmm_struct%iqm_atomic_numbers((qmmm_struct%nquant+qmmm_struct%nlink)), stat=ier )
   REQUIRE(ier == 0)
   allocate ( qmmm_struct%atom_mask(natm), stat=ier )
   REQUIRE(ier == 0)
   allocate ( qmmm_struct%mm_link_mask(natm), stat=ier )
   REQUIRE(ier == 0)
   allocate (qmmm_struct%scaled_mm_charges(natm), stat=ier )
   REQUIRE(ier == 0) 
   allocate (qmmm_scratch%qm_real_scratch(4*natm), stat=ier )
   REQUIRE(ier == 0)
   allocate (qmmm_scratch%qm_int_scratch(3*natm), stat=ier )
   REQUIRE(ier == 0)

   nquant_replicate = qmmm_struct%nquant
   natom_replicate = natm 

   !dxyzcl array only actually needs to be 3,qm_mm_pairs long..
   if (qmmm_nml%qmmm_int /= 0) then
     allocate ( qmmm_struct%dxyzcl(3,natom_replicate), stat=ier )
     REQUIRE(ier == 0) !Deallocated in deallocate qmmm
   end if
   !qm_xcrd only actually needs to be 4,qm_mm_pairs long...
   allocate (qmmm_struct%qm_xcrd(4,natom_replicate), stat=ier )
   REQUIRE(ier == 0) !Dealocated in deallocate qmmm

   return
end subroutine allocate_qmmm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates space for qmmm pair list
subroutine allocate_qmmm_pair_list( natom )
                                                                                                                                                
   implicit none
!Passed in
   integer, intent(in) :: natom  !Total number of REAL atoms in sim - MM + QM

!Local
   integer :: ier=0
   integer :: npairmax
   !NOTE we are overestimating the memory needed here.

!   allocate ( qmmm_struct%qm_mm_pair_list( (natom - qmmm_struct%nquant + 1) * qmmm_struct%nquant ),    &
!                                              stat=ier )
!All QM atoms share the same pair list:
   npairmax = natom - qmmm_struct%nquant+1
   allocate ( qmmm_struct%qm_mm_pair_list( npairmax ), stat=ier ) 
   REQUIRE(ier == 0)
   return
end subroutine allocate_qmmm_pair_list
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates space for qmmm variables and arrays
subroutine deallocate_qmmm()

   implicit none

!Local
   integer :: ier=0

   !If this is a parallel run the non master threads will only have
   !allocated this memory if LES is on since otherwise QM calc is
   !currently only done on master thread.
   !Deallocate pointers 
   if (qmmm_nml%qmqm_erep_incore) then
     deallocate ( qm2_struct%qm_qm_e_repul, stat = ier )
     REQUIRE(ier == 0)
   end if
   if (qmmm_nml%idc==0) then
     if (qm2_struct%n_peptide_links>0) then
       deallocate ( qm2_struct%peptide_links, stat=ier )
       REQUIRE(ier == 0)
     end if
     if (qmmm_nml%qmtheory /= DFTB) then
       deallocate ( qm2_struct%eigen_vectors, stat=ier )
       REQUIRE(ier == 0)
       if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
         deallocate ( qmmm_scratch%mat_diag_workspace, stat=ier )
         REQUIRE(ier == 0)
       end if
     end if
     if (qmmm_scratch%lapack_dc_real_scr_aloc /= 0) then
       deallocate ( qmmm_scratch%lapack_dc_real_scr, stat=ier)
       REQUIRE(ier == 0)
     end if
     if (qmmm_scratch%lapack_dc_int_scr_aloc /= 0) then
       deallocate ( qmmm_scratch%lapack_dc_int_scr, stat=ier)
       REQUIRE(ier == 0)
     end if
     if (qmmm_nml%allow_pseudo_diag .and. qmmm_mpi%commqmmm_master) then
       deallocate(qmmm_scratch%pdiag_scr_norbs_norbs, &
                  qmmm_scratch%pdiag_scr_noccupied_norbs, &
                  qmmm_scratch%pdiag_vectmp1, &
                  qmmm_scratch%pdiag_vectmp2, &
                  qmmm_scratch%pdiag_vectmp3, &
                  qmmm_scratch%pdiag_vecjs, &
                  stat=ier)
       REQUIRE(ier == 0)
     end if
     if (qmmm_mpi%commqmmm_master .and. qmmm_nml%diag_routine==7 .and. (.not. qmmm_nml%allow_pseudo_diag)) then
       !we allocated pdiag_scr_norbs_norbs to store the unpacked matrix so make sure
       !we deallocate it.
       deallocate(qmmm_scratch%pdiag_scr_norbs_norbs, stat=ier)
       REQUIRE(ier == 0)
     end if
     deallocate ( qm2_struct%fock2_ptot2, stat=ier )
     REQUIRE(ier == 0)
     deallocate ( qm2_struct%hmatrix, stat = ier )
     REQUIRE(ier == 0)
     deallocate (qm2_struct%qm_qm_2e_repul, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm2_struct%fock_matrix, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm2_struct%old2_density, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm2_struct%old_den_matrix, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm2_struct%den_matrix, stat = ier )
     REQUIRE(ier == 0)

     if (qmmm_nml%density_predict == 1) then
       !We are using Niklasson et al density matrix prediction algorithm.
       deallocate ( qm2_struct%md_den_mat_guess1, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qm2_struct%md_den_mat_guess2, stat = ier )
       REQUIRE(ier == 0)
     end if

     !parameter array
     if (qmmm_nml%qmtheory /= DFTB) then
       deallocate (qm2_params%atom_orb_pp_eqn_xxy2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_pp_eqn_xxy1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_sp_eqn_xx2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_sp_eqn_xx1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_sp_eqn_xy, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_ss_eqn_adb, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_pp_ovlp_ieqj2, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_pp_ovlp_ieqj1, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_pp_ovlp_inj, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_sp_ovlp, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_ss_eqn, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_zz_pxp_over_pap, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_zz_sxp_over_sap, stat = ier)
       REQUIRE(ier == 0)
       deallocate (qm2_params%atom_orb_zz_sxs_over_sas, stat = ier)
       REQUIRE(ier == 0)
       if (qmmm_nml%qmtheory==AM1 .OR. qmmm_nml%qmtheory==PM3 .OR. qmmm_nml%qmtheory==PDDGPM3 .OR. &
           qmmm_nml%qmtheory == PM3CARB1 .OR. qmmm_nml%qmtheory==RM1) then
          deallocate (qm2_params%NUM_FN, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%FN3, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%FN2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%FN1, stat = ier)
          REQUIRE(ier == 0)
       end if
       if (qmmm_nml%qmtheory==PDDGPM3 .OR. qmmm_nml%qmtheory==PDDGMNDO) then
          deallocate (qm2_params%pddge1, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%pddge2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%pddg_term1, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%pddg_term2, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%pddg_term3, stat = ier)
          REQUIRE(ier == 0)
          deallocate (qm2_params%pddg_term4, stat = ier)
          REQUIRE(ier == 0)
       end if
     end if
     deallocate (qm2_params%cc_exp_params, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%multip_2c_elec_params, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%onec2elec_params, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%orb_elec_ke, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%betasas, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%betasap, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%betapap, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%orb_loc, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%natomic_orbs, stat = ier)
     REQUIRE(ier == 0)
     deallocate (qm2_params%core_chg, stat = ier )
     REQUIRE(ier == 0)
     deallocate (qm2_params%pascal_tri1, stat = ier )
     REQUIRE(ier == 0)
     deallocate (qm2_params%pascal_tri2, stat = ier )
     REQUIRE(ier == 0)
   end if !if idc==0

   if (qmmm_nml%qmqmrij_incore) then
     if (qmmm_mpi%commqmmm_master) then
       deallocate (qm2_rij_eqns%qmqmrijdata, stat = ier )
       REQUIRE(ier == 0)
     end if
   end if
   if (qmmm_nml%qmmmrij_incore) then
     if (qmmm_mpi%commqmmm_master) then
       deallocate (qm2_rij_eqns%qmmmrijdata, stat = ier )
       REQUIRE(ier == 0)
     end if
   end if

   deallocate ( qmmm_struct%qm_coords, stat = ier)
   REQUIRE(ier == 0)   

   deallocate ( qmmm_struct%qm_xcrd, stat = ier )
   REQUIRE(ier == 0)

   if (qmmm_nml%qmmm_int /= 0) then
     deallocate ( qmmm_struct%dxyzcl, stat = ier )
     REQUIRE(ier == 0)
   end if

   deallocate ( qmmm_scratch%qm_int_scratch, stat = ier )
   REQUIRE(ier == 0)

   deallocate ( qmmm_scratch%qm_real_scratch, stat = ier )
   REQUIRE(ier == 0)

   deallocate ( qmmm_struct%dxyzqm, stat = ier )
   REQUIRE(ier == 0)

   deallocate ( qmmm_struct%scaled_mm_charges, stat = ier )
   REQUIRE(ier == 0)

   deallocate ( qmmm_struct%mm_link_mask, stat = ier)
   REQUIRE(ier == 0)   

   deallocate ( qmmm_struct%atom_mask, stat = ier)
   REQUIRE(ier == 0)   

   deallocate ( qmmm_struct%qm_resp_charges, stat = ier)
   REQUIRE(ier == 0)   

   if (qmmm_struct%nlink > 0 .and. qmmm_nml%idc==0) then
     deallocate ( qmmm_struct%link_pairs, stat = ier)
     REQUIRE(ier == 0)
     deallocate ( qmmm_struct%mm_link_pair_resp_charges, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qmmm_struct%mm_link_pair_saved_coords, stat = ier )
     REQUIRE(ier == 0)
   end if

   deallocate ( qmmm_nml%iqmatoms, stat = ier)
   REQUIRE(ier == 0)

   deallocate ( qmmm_struct%iqm_atomic_numbers, stat = ier)
   REQUIRE(ier == 0)

   deallocate (qmmm_struct%qm_atom_type, stat = ier)
   REQUIRE(ier == 0)

   deallocate ( qmmm_struct%qm_mm_pair_list, stat = ier)
   REQUIRE(ier == 0)

   !Deallocate qm_ewald memory if it is on.
   if ( qmmm_nml%qm_ewald>0 ) then
     deallocate ( qmewald%dmkv, stat = ier)
     REQUIRE(ier == 0)
     deallocate ( qmewald%dkvec, stat = ier)
     REQUIRE(ier == 0)
     deallocate ( qmewald%kvec, stat = ier)
     REQUIRE(ier == 0)
     if (.not. qmmm_nml%qm_pme) then
       deallocate ( qmewald%ktable, stat = ier )
       REQUIRE(ier == 0)
       deallocate ( qmewald%d_ewald_mm, stat = ier )
       REQUIRE(ier == 0)
     else
       deallocate ( qmmm_scratch%qm_pme_scratch, stat = ier )
       REQUIRE(ier == 0)
     end if
     deallocate ( qmewald%qmktable, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qmewald%mmpot, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qmewald%qmpot, stat = ier )
     REQUIRE(ier == 0)
   end if
   !Deallocate the scf Mulliken charge array.
   !Was allocated on all cpus.
   deallocate ( qm2_struct%scf_mchg, stat = ier )
   REQUIRE(ier == 0)

   !Deallocate QM-GB arrays - only used with qmgb=2.
   if ( qmmm_nml%qmgb == 2 ) then
     deallocate ( qm_gb%qmqm_onefij, stat = ier )
     REQUIRE(ier == 0)
     if ( qm_gb%saltcon_on) then
       deallocate ( qm_gb%qmqm_kappafij, stat = ier )
       REQUIRE(ier == 0)
     end if
     deallocate ( qm_gb%gb_mmpot, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm_gb%gb_qmpot, stat = ier )
     REQUIRE(ier == 0)
     deallocate ( qm_gb%qmqm_gb_list, stat = ier )
     REQUIRE(ier == 0)
   end if

!MPI Specific deallocations
#ifdef MPI
   deallocate(qmmm_mpi%nquant_nlink_jrange, stat=ier )
   REQUIRE(ier == 0)
# ifndef USE_MPI_IN_PLACE
   deallocate(qmmm_scratch%matsize_red_scratch, stat=ier )
   REQUIRE(ier == 0)
# endif
#endif
 
   ! GMS: Dellocate DFTB arrays
   if ((qmmm_nml%qmtheory == DFTB).and.(qmmm_mpi%commqmmm_master)) call qm2_dftb_deallocate

   return

end subroutine deallocate_qmmm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Reads the qmmm namelist and calls the qmmm memory allocation routines
#ifdef SQM
subroutine read_qmmm_namelist_and_allocate( natom_in, igb, atnam )
#else
subroutine read_qmmm_namelist_and_allocate( igb, ih, ix, x, cut, use_pme, ntb )
#endif

   use findmask
   use constants, only : RETIRED_INPUT_OPTION
   implicit none

!STATIC MEMORY
integer :: max_quantum_atoms  !Needed in read qmmm namelist since namelists cannot contain pointers
parameter ( max_quantum_atoms = 10000 )
!END STATIC MEMORY

!Passed in
   integer :: igb        !Value of igb from cntrl namelist
#ifdef SQM
   integer use_pme, ntb
   integer :: natom_in 
   character(len=8) atnam(*)
#else
   character(len=4) ih(*)
   integer, intent(in) :: ix(*)
   _REAL_, intent(in) :: x(*)
   _REAL_, intent(in) :: cut !MM-MM cutoff in angstroms
   integer, intent(in) :: use_pme, ntb
#endif

!local  
   _REAL_ :: qmcut      ! local copied to qmmm_nml%qmcut - specified cutoff to use for QM-MM electrostatics.
                         ! Default = same as regular MM cutoff.
   _REAL_ :: lnk_dis     ! Distance from the QM atom to place link atom.
   _REAL_ :: scfconv     ! local copied to qmmm_nml%scfconv - Convergence criteria for SCF routine. Default = 1.0D-8.
                         ! Minimum (tightest criteria) = 1.0D-16
   _REAL_ :: pseudo_diag_criteria !Criteria - maximum change in density matrix between successive SCF iterations
                                  !in which to allow pseudo diagonalisation. Default = 0.05.
   integer :: lnk_atomic_no !Atomic number of link atom
   integer :: qmgb       ! local copied to qmmm_nml%qmgb - flag for type of GB do with QM region
   integer :: qmtheory   ! local copied to qmmm_nml%qmtheory - flag for level of theory to use for QM region
   integer :: qmcharge   ! local copied to qmmm_nml%qmcharge - value of charge on QM system
   integer :: spin       ! local copied to qmmm_nml%spin - spin state of system
   integer :: i,j        ! temporary counter for local loops
   integer :: ifind
   integer :: qmqmdx     ! local copied to qmmm_nml%qmqm_analyt - 1 = analytical, 2 = numerical QM-QM derivatives in qm2
   integer :: verbosity  ! local copied to qmmm_nml%verbosity - Controls amount of info about QM part of calc that is printed (0=Def)
   integer :: tight_p_conv ! local copied to qmmm_nml%tight_p_conv - Controls convergence of density matrix. 1 (Def) = 0.05*sqrt(SCFCRT)
                           ! 2 = SCFCRT
   integer :: printcharges !Local copied to qmmm_nml%printcharges as a logical. 1 = true - print mulliken and cm1a and cm2a charges
                           !on every step. 0 = false = don't print charges. Default = 0 (.false.)
   integer :: peptide_corr !Local copied to the logical qmmm_nml%peptide_corr
                           !Add MM correction to peptide linkages 0 = No (Default), 1 = Yes.
   integer :: itrmax       !Local copied to qmmm_nml%itrmax - Maximum number of scf cycles to run
                           !before assuming convergence has failed (default = 1000)
   integer :: qmshake      !Local copied to qmmm_nml%qmshake - shake QM atoms if ntc>1?
   integer :: qmqmrij_incore !Flag to store rij between qm-qm pairs and related equations in memory.
                             !1 (default) = store in memory. 0 = calc on fly.
   integer :: qmmmrij_incore !Flag to store rij between qm-mm pairs and related equations in memory.
                             !1 (default) = store in memory. 0 = calc on fly.
   integer :: qmqm_erep_incore !Flag to store QM-QM 1 electron repulsion integrals in memory or to calculate
                               !them on the fly. Only available with QM-QM analytical derivatives.
                               !1 (default) = store in memory. 0 = calc on fly.
   integer :: pseudo_diag      !Whether to allow pseudo diagonalisations to be done when possible in SCF.
                               !0 (default) = Always do full diagonalisations.
                               !1 = do pseudo diagonalisations when possible.
   integer :: qm_ewald          !0 (default) do only regular QM-MM interaction in periodic calculations.
                                !1           do ewald based periodic QM-MM interactions.
                                !2           do ewald based periodic QM-MM but with the QM image charges
                                !            fixed at the previous steps mulliken charges during the SCF.
   integer :: qm_pme            !0 use regular Ewald for doing QM-MM interactions when qm_ewald>0.
                                !1 (default) use PME to do the reciprocal sum.
   integer :: kmaxqx, kmaxqy, kmaxqz !Maximum K space vectors
   integer :: ksqmaxq !Maximum K squared values for spherical cutoff in k space.
   integer :: writepdb
   integer :: qmmm_int !QM-MM interaction method
   integer :: adjust_q
   integer :: diag_routine !Controls diagonalization routine to use in SCF.
   integer :: density_predict !Controls prediction of density matrix for next SCF step.
   logical :: mdin_qmmm=.false.

   integer :: idc
   integer :: divpb
!! (GMS)
   _REAL_  :: chg_lambda       ! Charge scaling factor for free energy calculation
!! DFTB options
   integer :: dftb_maxiter     ! Max # of iterations before resetting Broyden (default: 70 ) ==> qmmm_nml%dftb_maxiter
   integer :: dftb_disper      ! Use dispersion?  (default: 0 = false) ==> qmmm_nml%dftb_disper
   integer :: dftb_chg         ! DFTB CM3 charges (default: 0 = Mulliken, 1 = CM3) ==> qmmm_nml%dftb_chg
   _REAL_  :: dftb_telec       ! Electronic temperature, in Kelvins. (Default = 0.0K) ==> qmmm_nml%dftb_telec 
   _REAL_  :: dftb_telec_step  ! Telec step size for convergence accelerator (Default = 0.0K) ==> qmmm_nml%dftb_telec_step
#include "memory.h"

!Apparently you can't use a pointer in a namelist :-( Therefore
!we need a local scratch array that will be big enough that 
!the iqmatoms list never exceeds it
   integer :: iqmatoms( max_quantum_atoms )
   character(len=1024) :: qmmask
   character(len=12) :: qm_theory !Options=PM3,AM1,MNDO,PDDG-PM3,PM3PDDG,PDDG-MNDO,PDDGMNDO,PM3-CARB1,PM3CARB1
                                !        DFTB,SCC-DFTB, RM1
   integer, dimension(:), pointer :: isqm
   integer :: ier=0

   namelist /qmmm/ qmcut, iqmatoms,qmmask,qmgb,qm_theory, qmtheory, &
                   qmcharge, qmqmdx, verbosity, tight_p_conv, scfconv, &
                   printcharges, peptide_corr, itrmax, qmshake, qmqmrij_incore, &
                   qmqm_erep_incore, qmmmrij_incore, &
                   lnk_dis, lnk_atomic_no, spin, pseudo_diag,   &
                   pseudo_diag_criteria, &
                   qm_ewald, qm_pme, kmaxqx, kmaxqy, kmaxqz, ksqmaxq, &
                   writepdb, qmmm_int, adjust_q, diag_routine, density_predict, idc, divpb, &
                   dftb_maxiter, dftb_disper, dftb_chg, dftb_telec, dftb_telec_step, &
                   chg_lambda

!Setup defaults
#ifdef SQM
   qmcut = 9999.d0
   use_pme = 0
   ntb = 0
#else
   qmcut = cut
#endif
   lnk_dis=1.09d0  !Methyl C-H distance
   lnk_atomic_no=1 !Hydrogen
   qmgb = 2 !Gets set to zero if igb==6 or igb==0.
   qm_theory = ''
   qmtheory = RETIRED_INPUT_OPTION 
   qmcharge = 0
   spin = 1
   qmqmdx = 1
   verbosity = 0
   tight_p_conv = 0
   scfconv = 1.0D-8
   printcharges = 0
   peptide_corr = 0
   itrmax = 1000
   qmshake = 1
   qmmask=''
   iqmatoms(1:max_quantum_atoms) = 0
   qmqmrij_incore = 1
   qmmmrij_incore = 1
   qmqm_erep_incore = 1
   pseudo_diag = 1
   pseudo_diag_criteria = 0.05
   qm_ewald=1 !Default is to do QMEwald, with varying charges, if ntb=0 or use_pme=0 then this will get turned off
   qm_pme = 1 !use pme for QM-MM
   kmaxqx=5; kmaxqy=5; kmaxqz=5    !Maximum K space vectors
   ksqmaxq=27 !Maximum K squared values for spherical cutoff in k space.
   writepdb = 0 !Set to 1 to write a pdb on the first step with just the QM region in it.
   qmmm_int = 1 !Default, do full interaction without extra Gaussian terms for PM3 / AM1 etc.
   adjust_q = 2 !Default adjust q over all atoms.
   diag_routine = 1 !Use default internal diagonalizer.
   density_predict = 0 !Use density matrix from previous MD step.
   idc = 0
   divpb = 0

   !DFTB
   dftb_maxiter     = 70   
   dftb_disper      = 0
   dftb_chg         = 0
   dftb_telec       = 0.0d0
   dftb_telec_step  = 0.0d0
   chg_lambda  = 1.0d0

   !Read qmmm namelist
   rewind 5

   call nmlsrc('qmmm',5,ifind)
   if (ifind /= 0) mdin_qmmm=.true.

   !Read qmmm namelist
   rewind 5
   if ( mdin_qmmm ) then
     read(5,nml=qmmm)
   else
     write(6, '(1x,a,/)') 'Could not find qmmm namelist'
     call mexit(6,1)
   endif

   !Parse qm_theory
   !1) if qm_theory = '' then the user didn't set anything - check qmtheory value for
   !   backwards compatibility.
   if ( qm_theory == '' ) then
     if ( qmtheory == RETIRED_INPUT_OPTION ) then
       !Neither was specified so we use the default of PM3.
       qmtheory = PM3
     end if
   else
     !Check if qmtheory was also specified along with qm_theory.
     if ( qmtheory /= RETIRED_INPUT_OPTION ) then
       call sander_bomb('read_qmmm_namelist','Both qm_theory and qmtheory are specified in the qmmm namelist.', &
                         'Specify ONE of either qm_theory or qmtheory. Note use of qmtheory is now deprecated.')
     end if
  
     !No qmtheory specified in input file so pass qm_theory text string.
     if ( qm_theory == 'PM3' ) then
       qmtheory = PM3
     else if ( qm_theory == 'AM1' ) then
       qmtheory = AM1
     else if ( qm_theory == 'MNDO' ) then
       qmtheory = MNDO
     else if ( qm_theory == 'PM3-PDDG' .OR. qm_theory == 'PM3PDDG' .OR. qm_theory == 'PM3_PDDG' &
          .OR. qm_theory == 'PDDG-PM3' .OR. qm_theory == 'PDDGPM3' .OR. qm_theory == 'PDDG_PM3' ) then
       qmtheory = PDDGPM3
     else if ( qm_theory == 'MNDO-PDDG' .OR. qm_theory == 'MNDOPDDG' .OR. qm_theory == 'MNDO_PDDG' &
          .OR. qm_theory == 'PDDG-MNDO' .OR. qm_theory == 'PDDGMNDO' .OR. qm_theory == 'PDDG_MNDO' ) then
       qmtheory = PDDGMNDO
     else if ( qm_theory == 'PM3-CARB1' .OR. qm_theory == 'PM3CARB1' .OR. qm_theory == 'PM3_CARB1' ) then
       qmtheory = PM3CARB1
     else if ( qm_theory == 'DFTB' .OR. qm_theory == 'SCCDFTB' .OR. qm_theory == 'SCC-DFTB' .OR. &
               qm_theory == 'SCC_DFTB' ) then
       qmtheory = DFTB
     else if ( qm_theory == 'RM1' ) then
       qmtheory = RM1
     else
        call sander_bomb('read_qmmm_namelist','Unknown method specified for qm_theory', &
                         'Valid options are: PM3, AM1, RM1, MNDO, PM3-PDDG, MNDO-PDDG, PM3-CARB1 and DFTB')
     endif
   end if !( qm_theory == '' ) 
   
#ifdef SQM
      natom = natom_in
      qmmm_struct%nquant = natom
      qmmm_struct%nlink = 0
      qmmm_struct%nquant_nlink = natom
      do i=1,natom
         iqmatoms(i) = i
      end do
#else
   if( qmmask /= '' ) then  !  get the quantum atoms from the mask
      write(6,'(a)') ''
      write(6,'(a)') 'LOADING THE QUANTUM ATOMS AS GROUPS'

      allocate(isqm( natom ), stat=ier)
      REQUIRE(ier==0)

      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), qmmask, isqm )
      qmmm_struct%nquant = sum(isqm(1:natom))
      write(6,'(a,a,a,i5,a)') '     Mask ', qmmask(1:len_trim(qmmask)), &
         ' matches ',qmmm_struct%nquant,' atoms'

      j = 0
      do i=1,natom
         if( isqm(i)>0 ) then
            j = j+1
            iqmatoms(j) = i
         end if
      end do

      deallocate(isqm, stat=ier)
      REQUIRE(ier==0)

   else  !  get the count from the input iqmatoms array

      do i = 1,max_quantum_atoms
         if( iqmatoms(i) == 0 ) exit
      end do
      qmmm_struct%nquant = i-1

   end if

!Initialize nlink to 0
   qmmm_struct%nlink = 0
   qmmm_struct%nquant_nlink = qmmm_struct%nquant

#endif

! Test to see if QM atom selection is legal.
   call validate_qm_atoms(iqmatoms,qmmm_struct%nquant,natom)

!  check we don't bust our statically allocated max_quantum_atoms
   call int_legal_range('QMMM: (number of quantum atoms) ', &
      qmmm_struct%nquant, 1, max_quantum_atoms )

   call qmsort(iqmatoms) !ensure the list of qm atoms is sorted numerically

   call float_legal_range('QMMM: (QM-MM Cutoff) ', qmcut,0.0D0,1.0D30)
   call int_legal_range('QMMM: (QM GB Method) ', qmgb,0,3)
   call int_legal_range('QMMM: (QM Theory) ', qmtheory,1,8)
   call int_legal_range('QMMM: (QM-QM Derivatives) ', qmqmdx,1,2)
   call int_legal_range('QMMM: (Verbosity) ', verbosity,0,5)
   call int_legal_range('QMMM: (Max SCF Iterations) ', itrmax,1,10000000)
   call int_legal_range('QMMM: (Shake on QM atoms) ', qmshake,0,1)
   call int_legal_range('QMMM: (Density Matrix Convergence) ', tight_p_conv,0,1)
   call float_legal_range('QMMM: (SCF Convergence) ', scfconv,1.0D-16,1.0D0)
   call int_legal_range('QMMM: (PRINT CHARGES) ', printcharges,0,1)
   call int_legal_range('QMMM: (Spin State) ', spin,1,1)
!RCW: Currently limit spin state to singlets only since the code for spin>1 does not exist / work at present.
!     WARNING - IF WE LATER ALLOW SPIN>1 qm2_densit will need updating.
   call int_legal_range('QMMM: (Peptide Correction) ',peptide_corr,0,1)
   call int_legal_range('QMMM: (QM-QM RIJ in Core) ',qmqmrij_incore,0,1)
   call int_legal_range('QMMM: (QM-MM RIJ in Core) ',qmmmrij_incore,0,1)
   call int_legal_range('QMMM: (QM-QM E-Rep in Core) ',qmqm_erep_incore,0,1)
   call int_legal_range('QMMM: (Link Atomic Number) ',lnk_atomic_no,1,nelements)
   call int_legal_range('QMMM: (Pseudo Diag) ',pseudo_diag,0,1)
   call int_legal_range('QMMM: (QM Ewald) ',qm_ewald,0,2)
   call int_legal_range('QMMM: (QM PME) ',qm_pme,0,1)
   call int_legal_range('QMMM: (QM Ewald kmaxqx) ',kmaxqx,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqy) ',kmaxqy,1,99999999)
   call int_legal_range('QMMM: (QM Ewald kmaxqz) ',kmaxqz,1,99999999)
   call int_legal_range('QMMM: (QM Ewald ksqmaxq) ',ksqmaxq,1,kmaxqx*kmaxqy*kmaxqz)
   call int_legal_range('QMMM: (QM-MM qmmm_int) ',qmmm_int,0,2)
   call int_legal_range('QMMM: (QM-MM adjust_q) ',adjust_q,0,2)
   call int_legal_range('QMMM: (QM-MM diag_routine) ',diag_routine,0,7)
   call int_legal_range('QMMM: (QM-MM density_predict) ',density_predict,0,1)
   call float_legal_range('QMMM: (Pseudo Diag Criteria) ',pseudo_diag_criteria,1.0D-12,1.0D0)
   call float_legal_range('QMMM: (Link Atom Distance) ',lnk_dis,0.7D0,4.0D0)

!! GMS
   call float_legal_range('QMMM: (QM-MM chg_lambda)'     , chg_lambda     , 0.0D0 , 1.0D0  )
   qmmm_nml%chg_lambda  = chg_lambda
!! DFTB
   call int_legal_range(  'QMMM: (QM-MM dftb_maxiter ) ' , dftb_maxiter   , 1     , 10000  )
   call int_legal_range(  'QMMM: (QM-MM dftb_disper) '   , dftb_disper    , 0     , 1      )
   call int_legal_range(  'QMMM: (QM-MM dftb_chg   ) '   , dftb_chg       , 0     , 1      )
   call float_legal_range('QMMM: (QM-MM dftb_telec)'     , dftb_telec     , 0.0D0 , 1.0D4  )

   qmmm_nml%dftb_maxiter   = dftb_maxiter
   qmmm_nml%dftb_disper      = dftb_disper
   qmmm_nml%dftb_chg         = dftb_chg
   qmmm_nml%dftb_telec       = dftb_telec
   qmmm_nml%dftb_telec_step  = dftb_telec_step

   if (dftb_chg > 0) printcharges=1
!!

   qmmm_nml%qmcut = qmcut
   qmmm_nml%qmcut2 = qmcut*qmcut
   qmmm_nml%lnk_dis = lnk_dis
   qmmm_nml%lnk_atomic_no = lnk_atomic_no

!DFTB Limitations - current things not supported in DFTB
!These are silent limitations that are non fatal - just to
!avoid problems with the default. Fatal errors are handled
!later in this routine.
   if (qmtheory == DFTB) then
     qmqmrij_incore=0
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
     qmqmdx=1
   end if

!Divcon limitations - current things not supported in sander.DIVCON
!These are silent changes to remove defaults. Fatal errors are handled
!later in this routine.
   if (idc /= 0) then
     qmqmrij_incore=0
     qmmmrij_incore=0
     qmqm_erep_incore=0
     pseudo_diag=0
   end if

   !qmgb values:
   ! 0 - do GB but leave QM charges as zero and add nothing to Fock matrix. This is like a vacuum
   !     QM molecule in a solvated MM system.
   ! 1 - do GB using the prmtop fixed resp charges for the GB calculation.
   ! 2 - do GB using Mulliken charges that are consistent with the GB field by modifying the fock
   !     matrix at every SCF step. (default)
   ! 3 - do GB using QM gas phase Mulliken charges - This is really a debugging option since the charges
   !     will not be consistent with the GB field since the fock matrix is not modified. This similarly
   !     means that the gradients will not be accurate. A warning will be printed at every QM call if this
   !     option is selected.

   !Make sure igb in &cntrl namelist is compatible with qmgb setting.
   if (igb==0 .or. igb==6) then
      !no qmgb available
      qmgb = 0
   end if
   !Print warning about qmgb being for debugging only.
   if (qmgb==3) then
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
     write(6,*) "QMMM: qmgb = 3 is designed for debugging purposes only. It gives GB"
     write(6,*) "QMMM:          energies based on gas phase QM Mulliken charges. These charges"
     write(6,*) "QMMM:          are NOT consistent with the GB field felt by the QM region and"
     write(6,*) "QMMM:          so any gradients calculated using this approach will"
     write(6,*) "QMMM:          NOT BE ACCURATE."
     write(6,*) "QMMM:          This option is really designed for:"
     write(6,*) "QMMM:                SINGLE POINT ENERGY EVALUATIONS ONLY"
     write(6,*) "QMMM: ------------------------------ WARNING --------------------------------"
   end if
   qmmm_nml%qmgb = qmgb

   qmmm_nml%qmtheory = qmtheory
   qmmm_struct%AM1_OR_PM3 = (qmmm_nml%qmtheory == AM1 .OR. qmmm_nml%qmtheory == PM3 &
                             .OR. qmmm_nml%qmtheory == PDDGPM3 .OR. qmmm_nml%qmtheory == PM3CARB1 .OR. &
                             qmmm_nml%qmtheory == RM1)
   qmmm_struct%PDDG_IN_USE = (qmmm_nml%qmtheory == PDDGPM3 .OR. qmmm_nml%qmtheory == PDDGMNDO )
   qmmm_nml%qmcharge = qmcharge
   qmmm_nml%spin = spin
   qmmm_nml%verbosity = verbosity
   qmmm_nml%itrmax = itrmax
   qmmm_nml%qmshake = qmshake
   qmmm_nml%pseudo_diag_criteria = pseudo_diag_criteria
   if (qmqmdx /= 1) then
      qmmm_nml%qmqm_analyt = .false. !Do numerical QM-QM derivatives in qm2
   else
      qmmm_nml%qmqm_analyt = .true.  !Do analytical QM-QM dericatives in qm2
   end if
   if (tight_p_conv /= 1) then
      qmmm_nml%tight_p_conv = .false. !Loose density matrix convergence (0.05*sqrt(SCFCRT))
   else
      qmmm_nml%tight_p_conv = .true.  !Tight density matrix convergence (SCFCRT)
   end if
   !Write a warning about excessively tight convergence requests.
   if ( scfconv < 1.0D-12 ) then
     write(6,'(" QMMM: WARNING - SCF Conv = ",G8.2)') scfconv
     write(6,*) "QMMM:           There is a risk of convergence problems when the"
     write(6,*) "QMMM:           requested convergence is less that 1.0D-12 kcal/mol."
   end if
   qmmm_nml%scfconv = scfconv
   !How tight do we want the density convergence?
   if (qmmm_nml%tight_p_conv) then
      qmmm_nml%density_conv = qmmm_nml%scfconv
   else
      qmmm_nml%density_conv = 0.05D0 * sqrt(qmmm_nml%scfconv)
   end if

   if ( printcharges /= 1) then
      qmmm_nml%printcharges=.false.
   else
      qmmm_nml%printcharges=.true.
   end if
   if ( peptide_corr == 0) then
      qmmm_nml%peptide_corr = .false.
   else
      qmmm_nml%peptide_corr =  .true.
   end if 
   if ( qmqmrij_incore == 0 .or. qmqmdx == 2 ) then
      qmmm_nml%qmqmrij_incore = .false.
   else
      !Only available with analytical derivatives
      qmmm_nml%qmqmrij_incore = .true.
   end if
   if ( qmmmrij_incore == 0 .or. qmmm_int==0 ) then
      qmmm_nml%qmmmrij_incore = .false.
   else
      qmmm_nml%qmmmrij_incore = .true. !Only available with qmmm_int=1 or qmmm_int=2
   end if
   if ( qmqm_erep_incore == 0 .or. qmqmdx == 2 ) then
      qmmm_nml%qmqm_erep_incore = .false.
   else
      !Only available with analytical derivatives.
      qmmm_nml%qmqm_erep_incore = .true.
   end if
   if ( pseudo_diag == 1 ) then
     qmmm_nml%allow_pseudo_diag = .true.
   else
     qmmm_nml%allow_pseudo_diag = .false.
   end if
 
   qmmm_nml%qm_ewald = qm_ewald
   qmmm_nml%ksqmaxq = ksqmaxq
   qmmm_nml%kmaxqx = kmaxqx
   qmmm_nml%kmaxqy = kmaxqy
   qmmm_nml%kmaxqz = kmaxqz
   !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
   !have put in the namelist and set the value to false.
   if (ntb==0 .or. use_pme==0) then
     qmmm_nml%qm_ewald = 0
     qmmm_nml%qm_pme = .false.
   end if
   if (qmmm_nml%qm_ewald>0 .and. qm_pme>0) then
     qmmm_nml%qm_pme=.true.
   else
     qmmm_nml%qm_pme=.false.
   end if

   if ( writepdb == 0 ) then
     qmmm_nml%writepdb=.false.
   else
     qmmm_nml%writepdb=.true.
   end if
   qmmm_nml%qmmm_int = qmmm_int

   qmmm_nml%idc = idc
   qmmm_nml%divpb = divpb

!Setup some specific calculation flags that depend on namelist variables.
!Need to make sure these get copied to other threads in an MPI run.

   !Will we be calculating the Mulliken charges on every SCF iteration?
   !Default is no. Will be set to true in a bit if certain options, such as qm_ewald
   !require it.
   if (qmmm_nml%qm_ewald==1 .or. qmmm_nml%qmgb>1) then
     !We will be needing the mulliken charges on every SCF iteration
     qm2_struct%calc_mchg_scf = .true.
   else
     qm2_struct%calc_mchg_scf = .false.
   end if

   !DFTB Calculates Mulliken charges anyway so we might as well store them in the correct place.
   if (qmmm_nml%qmtheory == DFTB) qm2_struct%calc_mchg_scf = .true.

!At this point we know nquant and natom so we can allocate our arrays that depend on nquant or natom
!Note if this is a LES run qmmm_struct%nquant is nqaunt
!Note non master mpi threads need to call this allocation routine manually themselves.
   call allocate_qmmm( natom )

#ifndef SQM
!Copy the iqmatoms from our static local array to our allocated array.
   if(qmmm_nml%divpb == 1)then
      !the +100 is for link atoms, needs some way to figure out actual number of link atoms
      allocate(qmmm_struct%all_atom_numbers(natom+100), stat=ier)
      REQUIRE(ier == 0)

      do i=1,natom
         call get_atomic_number(ih(m04+i-1),x(lmass+i-1),qmmm_struct%all_atom_numbers(i))
      enddo
   endif
#endif

   do i = 1, qmmm_struct%nquant
      qmmm_nml%iqmatoms(i) = iqmatoms(i)
      !Get the atomic numbers (used to be done in rdparm2...)
#ifdef SQM
      call get_atomic_number( atnam(i), x(lmass+i-1), qmmm_struct%iqm_atomic_numbers(i) )
#else
      j = iqmatoms(i)
      call get_atomic_number( ih(m04+j-1), x(lmass+j-1), qmmm_struct%iqm_atomic_numbers(i) )
#endif
   end do

!Now we have a list of atom numbers for QM atoms we can build a true false (natom long) list
!specifying what the quantum atoms are. Useful for doing quick .OR. operations against other
!lists.
   qmmm_struct%atom_mask = .false. !Note, sets entire natom long array to false

   do i = 1, qmmm_struct%nquant
     qmmm_struct%atom_mask(qmmm_nml%iqmatoms(i)) = .true.
   end do

   qmmm_nml%adjust_q = adjust_q

   qmmm_nml%diag_routine = diag_routine
   qmmm_nml%density_predict = density_predict

! --- CHECK FOR LIMITATIONS ---


  !--- DIVCON LIMITATIONS ---
  !Divcon currently only works with gas phase simulations.
  if (qmmm_nml%idc>0) then
    if (ntb /= 0) then
      !This covers qmewald, qm_pme as well.
      call sander_bomb('read_qmmm_namelist','idc /= 0 (Divcon on) but periodic boundaries are in use.', &
                       'Periodic boundaries are currently only supported when idc == 0')
    end if
    if (qmmm_nml%qmgb/=0) then
      call sander_bomb('read_qmmm_namelist','idc /= 0 (Divcon on) but qmgb/=0. QMMM GB is currently', &
                       'only supported when idc == 0')
    end if
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_namelist','idc /= 0 (Divcon on) but peptide_corr /= 0.', &
                       'Peptide correction for divcon is handled by the divcon.in file. Set peptide_corr = 0 to proceed.')
    end if
    if (qmmm_nml%printcharges) then
      call sander_bomb('read_qmmm_namelist','idc /= 0 (Divcon on) but printcharges /= 0.', &
                       'Printcharges is not available with idc > 0.')
    end if
#ifdef MPI
    !No support for parallel Divcon
    write(6,*) 'Divcon capability (idc>0) can only run in serial mode for now'
    call mexit(6,1)
#endif
    !Write a warning about qmtheory being ignored with Divcon.
    write(6,'("|QMMM: WARNING DIVCON IN USE")')
    write(6,'("|QMMM: qm_theory IS IGNORED WHEN USING DIVCON - QM HAMILTONIAN MUST BE SELECTED")')
    write(6,'("|QMMM: IN DIVCON.IN FILE.")')
  end if
  !--- END DIVCON LIMITATIONS ---

  !--- DFTB LIMITATIONS ---
  if (qmmm_nml%qmtheory == DFTB ) then
    !if (qmmm_nml%qmgb > 1) then
    !  call sander_bomb('read_qmmm_namelist', 'qmtheory = DFTB but qmgb > 1. DFTB does not currently', &
    !                   'support full QM GB. Only qmgb values of 0 or 1 are supported.')
    !end if
    if (qmmm_nml%peptide_corr) then
      call sander_bomb('read_qmmm_namelist','qm_theory=DFTB but peptide_corr /= 0.', &
                       'Peptide correction is not available, or required for  DFTB. Set peptide_corr = 0 to proceed.')
    end if
  end if
  !--- END DFTB LIMITATIONS ---

! --- END CHECK FOR LIMITATIONS ---

  return

end subroutine read_qmmm_namelist_and_allocate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+Sorts the array iqmatoms into numerical order
subroutine qmsort( iqmatoms)

   implicit none
   integer, intent(inout) :: iqmatoms(*)

! Local
   integer i,j,lcurrent
!
!   sort array iqmatoms in ascending order
!
   do i = 1, qmmm_struct%nquant
      lcurrent = iqmatoms(i)
      do j = i+1,qmmm_struct%nquant
         if (lcurrent.gt.iqmatoms(j)) then
            iqmatoms(i) = iqmatoms(j)
            iqmatoms(j) = lcurrent
            lcurrent = iqmatoms(i)
         endif
      end do
   end do
end subroutine qmsort
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+reallocation of integer array
subroutine reallocate_iqm_atno_array(orig_size,new_size)
      !This routine reallocates iqm_atomic_numbers which was
      !originally allocated as orig_size as new_size
      !array ends up being padded with zeros
                                                                                                                                                             
      implicit none

      !Passed in
      integer, intent(in) :: orig_size, new_size
                                                                                                                                                             
      !Local
      integer, dimension(:), pointer :: temp_array
      integer i
      integer :: ier=0

      allocate(temp_array(new_size),stat=ier)
      REQUIRE(ier==0)

      temp_array = zero 
                                                                                                                                                             
      do i=1, orig_size   !Copy out the data to a temporary array
         temp_array(i) = qmmm_struct%iqm_atomic_numbers(i)
      end do
                                                                                                                                                             
      !Deallocate and allocate it as bigger
      deallocate ( qmmm_struct%iqm_atomic_numbers, stat = ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_struct%iqm_atomic_numbers(new_size), stat=ier )
      REQUIRE(ier == 0)

      !copy in the data
      do i=1, new_size
         qmmm_struct%iqm_atomic_numbers(i) = temp_array(i) 
      end do

      deallocate(temp_array,stat=ier)
      REQUIRE(ier==0)

      return
                                                                                                                                                             
end subroutine reallocate_iqm_atno_array

subroutine reallocate_iqmatoms_array(orig_size,new_size)
      !This routine reallocates iqmatoms which was
      !originally allocated as orig_size as new_size
      !array ends up being padded with zeros
                                                                                                                                                             
      implicit none
                                                                                                                                                             
      !Passed in
      integer, intent(in) :: orig_size, new_size
                                                                                                                                                             
      !Local
      integer, dimension(:), pointer :: temp_array
      integer i
      integer :: ier=0
                                                                                                                                                             
      allocate(temp_array(new_size),stat=ier)
      REQUIRE(ier==0)
                                                                                                                                                             
      temp_array = zero
                                                                                                                                                             
      do i=1, orig_size   !Copy out the data to a temporary array
         temp_array(i) = qmmm_nml%iqmatoms(i)
      end do
                                                                                                                                                             
      !Deallocate and allocate it as bigger
      deallocate ( qmmm_nml%iqmatoms, stat = ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_nml%iqmatoms(new_size), stat=ier )
      REQUIRE(ier == 0)
                                                                                                                                                             
      !copy in the data
      do i=1, new_size
         qmmm_nml%iqmatoms(i) = temp_array(i)
      end do
                                                                                                                                                             
      deallocate(temp_array,stat=ier)
      REQUIRE(ier==0)
                                                                                                                                                             
      return
                                                                                                                                                             
end subroutine reallocate_iqmatoms_array

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+Identify qm atom elements.
      subroutine get_atomic_number(atom_name,atom_mass,atomic_number)
!
!
!     This subroutine assigns atomic number based upon the first
!     letter of the atom symbol read in from the topology file and the mass.
!     The assumption here is that an atom name matches the first letter of the
!     element. If it doesn't then this routine will need to be modified.

!     atom_name ==>   character containing the atomic name
!     iqm_atomic_numbers ==>   integer array of atomic numbers assigned to atoms
!
      implicit none

!passed in
      integer :: atomic_number
      _REAL_ :: atom_mass
      character(len=4) :: atom_name

!local
      integer :: i, j

!Nobel Gasses are not supported.
!Lanthanides are not supported.
!Actinides are not supported.

      if(atom_name(1:1) .eq. 'a' .or. atom_name(1:1) .eq. 'A') then
        if(atom_mass > 24.0d0 .and. atom_mass <= 28.0d0) then
          atomic_number =  13 !Aluminium
        elseif(atom_mass > 73.0d0 .and. atom_mass <= 77.0d0) then
          atomic_number =  33 !Arsenic
        elseif(atom_mass > 106.0d0 .and. atom_mass <= 109.0d0) then
          atomic_number =  47 !Silver
        elseif(atom_mass > 195.0d0 .and. atom_mass <= 199.0d0) then
          atomic_number =  79 !Gold
        elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number =  85 !Astatine
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'b' .or. atom_name(1:1) .eq. 'B') then
        if(atom_mass > 8.0d0 .and. atom_mass <= 10.0d0) then
          atomic_number =  4 !Beryllium
        elseif(atom_mass > 10.0d0 .and. atom_mass <= 12.0d0) then
          atomic_number =  5 !Boron
        elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number =  35 !Bromine
        elseif(atom_mass > 135.0d0 .and. atom_mass <= 139.0d0) then
          atomic_number =  56 !Barium
        elseif(atom_mass > 207.0d0 .and. atom_mass <= 211.0d0) then
          atomic_number =  83 !Bismuth
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'c' .or. atom_name(1:1) .eq. 'C') then
        if(atom_mass > 10.0d0 .and. atom_mass <= 14.0d0) then
          atomic_number =  6 !Carbon
        elseif(atom_mass > 33.0d0 .and. atom_mass <= 37.0d0) then
          atomic_number =  17 !Chlorine
        elseif(atom_mass > 38.0d0 .and. atom_mass <= 42.0d0) then
          atomic_number =  20 !Calcium
        elseif(atom_mass > 50.0d0 .and. atom_mass <= 54.0d0) then
          atomic_number =  24 !Chromium
        elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number =  27 !Cobalt
        elseif(atom_mass > 61.0d0 .and. atom_mass <= 65.0d0) then
          atomic_number =  29 !Copper
        elseif(atom_mass > 110.0d0 .and. atom_mass <= 114.0d0) then
          atomic_number =  48 !Cadmium
        elseif(atom_mass > 131.0d0 .and. atom_mass <= 135.0d0) then
          atomic_number =  55 !Cesium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'd' .or. atom_name(1:1) .eq. 'D') then
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'e' .or. atom_name(1:1) .eq. 'E') then
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'f' .or. atom_name(1:1) .eq. 'F') then
        if(atom_mass > 17.0d0 .and. atom_mass <= 21.0d0) then
          atomic_number =  9 !Fluorine
        elseif(atom_mass > 54.0d0 .and. atom_mass <= 58.0d0) then
          atomic_number =  26 !Iron
        elseif(atom_mass > 218.0d0 .and. atom_mass <= 228.0d0) then
          atomic_number =  87 !Francium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'g' .or. atom_name(1:1) .eq. 'G') then
        if(atom_mass > 67.0d0 .and. atom_mass <= 71.0d0) then
          atomic_number =  31 !Gallium
        elseif(atom_mass > 71.0d0 .and. atom_mass <= 75.0d0) then
          atomic_number =  32 !Germanium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'h' .or. atom_name(1:1) .eq. 'H') then
        if(atom_mass > 0.0d0 .and. atom_mass <= 2.0d0) then
          atomic_number =  1 !Hydrogen
        elseif(atom_mass > 3.0d0 .and. atom_mass <= 5.0d0) then
          atomic_number =  2 !Helium
        elseif(atom_mass > 176.0d0 .and. atom_mass <= 180.0d0) then
          atomic_number =  72 !Hafnium
        elseif(atom_mass > 198.0d0 .and. atom_mass <= 202.0d0) then
          atomic_number =  80 !Mercury
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'i' .or. atom_name(1:1) .eq. 'I') then
        if(atom_mass > 112.0d0 .and. atom_mass <= 116.0d0) then
          atomic_number = 49 !Indium
        elseif(atom_mass > 125.0d0 .and. atom_mass <= 129.0d0) then
          atomic_number =  53 !Iodine
        elseif(atom_mass > 190.0d0 .and. atom_mass <= 194.0d0) then
          atomic_number =  77 !Iridium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'j' .or. atom_name(1:1) .eq. 'J') then
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'k' .or. atom_name(1:1) .eq. 'K') then
        if(atom_mass > 37.0d0 .and. atom_mass <= 41.0d0) then
          atomic_number = 19 !Potassium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'l' .or. atom_name(1:1) .eq. 'L') then
        if(atom_mass > 6.0d0 .and. atom_mass <= 8.0d0) then
          atomic_number = 3 !Lithium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'm' .or. atom_name(1:1) .eq. 'M') then
        if(atom_mass > 22.0d0 .and. atom_mass <= 26.0d0) then
          atomic_number = 12 !Magnesium
        elseif(atom_mass > 53.0d0 .and. atom_mass <= 57.0d0) then
          atomic_number = 25 !Manganese
        elseif(atom_mass > 94.0d0 .and. atom_mass <= 98.0d0) then
          atomic_number = 42 !Molybdenem
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'n' .or. atom_name(1:1) .eq. 'N') then
        if(atom_mass > 13.0d0 .and. atom_mass <= 15.0d0) then
          atomic_number = 7 !Nitrogen
        elseif(atom_mass > 21.0d0 .and. atom_mass <= 25.0d0) then
          atomic_number = 11 !Sodium
        elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number = 28 !Nickel
        elseif(atom_mass > 95.0d0 .and. atom_mass <= 99.0d0) then
          atomic_number = 41 !Niobium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'o' .or. atom_name(1:1) .eq. 'O') then
        if(atom_mass > 14.0d0 .and. atom_mass <= 18.0d0) then
          atomic_number = 8 !Oxygen
        elseif(atom_mass > 188.0d0 .and. atom_mass <= 192.0d0) then
          atomic_number = 76 !Osmium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'p' .or. atom_name(1:1) .eq. 'P') then
        if(atom_mass > 29.0d0 .and. atom_mass <= 33.0d0) then
          atomic_number = 15 !Phosphorus
        elseif(atom_mass > 104.0d0 .and. atom_mass <= 108.0d0) then
          atomic_number = 46 !Palladium
        elseif(atom_mass > 193.0d0 .and. atom_mass <= 197.0d0) then
          atomic_number = 78 !Platinum
        elseif(atom_mass > 205.0d0 .and. atom_mass <= 208.0d0) then
          atomic_number = 82 !Lead
        elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number = 84 !Polonium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'q' .or. atom_name(1:1) .eq. 'Q') then
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'r' .or. atom_name(1:1) .eq. 'R') then
        if(atom_mass > 84.0d0 .and. atom_mass <= 88.0d0) then
          atomic_number = 37 !Rubidium
        elseif(atom_mass > 99.0d0 .and. atom_mass <= 102.0d0) then
          atomic_number = 44 !Ruthenium
        elseif(atom_mass > 102.0d0 .and. atom_mass <= 105.0d0) then
          atomic_number = 45 !Rhodium
        elseif(atom_mass > 184.0d0 .and. atom_mass <= 188.0d0) then
          atomic_number = 75 !Rhenium
        elseif(atom_mass > 220.0d0 .and. atom_mass <= 230.0d0) then
          atomic_number = 88 !Radium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        end if
      elseif(atom_name(1:1) .eq. 's' .or. atom_name(1:1) .eq. 'S') then
        if(atom_mass > 26.0d0 .and. atom_mass <= 30.0d0) then
          atomic_number = 14 !Silicon
        elseif(atom_mass > 30.0d0 .and. atom_mass <= 34.0d0) then
          atomic_number = 16 !Sulphur
        elseif(atom_mass > 43.0d0 .and. atom_mass <= 47.0d0) then
          atomic_number = 21 !Scandium
        elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number = 34 !Selenium
        elseif(atom_mass > 86.0d0 .and. atom_mass <= 89.0d0) then
          atomic_number = 38 !Strontium
        elseif(atom_mass > 116.0d0 .and. atom_mass <= 120.0d0) then
          atomic_number = 50 !Tin
        elseif(atom_mass > 120.0d0 .and. atom_mass <= 124.0d0) then
          atomic_number = 51 !Antimony
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        end if
      elseif(atom_name(1:1) .eq. 't' .or. atom_name(1:1) .eq. 'T') then
        if(atom_mass > 46.0d0 .and. atom_mass <= 50.0d0) then
          atomic_number = 22 !Titanium
        elseif(atom_mass > 96.0d0 .and. atom_mass <= 100.0d0) then
          atomic_number = 43 !Technetium
        elseif(atom_mass > 125.0d0 .and. atom_mass <= 130.0d0) then
          atomic_number = 52 !Tellurium
        elseif(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 73 !Tantalum
        elseif(atom_mass > 201.0d0 .and. atom_mass <= 206.0d0) then
          atomic_number = 81 !Thallium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'u' .or. atom_name(1:1) .eq. 'U') then
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'v' .or. atom_name(1:1) .eq. 'V') then
        if(atom_mass > 49.0d0 .and. atom_mass <= 53.0d0) then
          atomic_number = 23 !Vanadium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'w' .or. atom_name(1:1) .eq. 'W') then
        if(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 74 !Tungsten
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        endif
      elseif(atom_name(1:1) .eq. 'x' .or. atom_name(1:1) .eq. 'X') then
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
      elseif(atom_name(1:1) .eq. 'y' .or. atom_name(1:1) .eq. 'Y') then
        if(atom_mass > 86.0d0 .and. atom_mass <= 90.0d0) then
          atomic_number = 39 !Yttrium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        end if
      elseif(atom_name(1:1) .eq. 'z' .or. atom_name(1:1) .eq. 'Z') then
        if(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
          atomic_number = 40 !Zirconium
        else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
        end if
      else
        write(6,*) 'Unable to correctly identify element ', atom_name
        call mexit(6,1)
      endif

      return

end subroutine get_atomic_number

subroutine validate_qm_atoms(iqmatoms, nquant, natom)
! This routine will check the list of atoms numbers stored in iqmatoms
! and check the following:
!
! 1) All are >= 1 .and. <= natoms
! 2) All are unique integer numbers
!
! Written by Ross Walker, TSRI, 2004
!
!=======================================================
     implicit none

!Passed in
   integer, intent(in) :: nquant, natom
   integer, intent(in) :: iqmatoms(nquant)

!Local
   integer icount1, icount2, iatom

   ! Sanity check 1, ensure nquant isn't bigger than natom (it can't be)
   if ((nquant < 1) .OR. (nquant > natom)) then
      write (6,'(" QM ATOM VALIDATION: nquant has a value of ",i8)') nquant
      write (6,'(" which is bigger than natom of ",i8,". Need 0 < nquant <= natom.")') natom
      call sander_bomb('validate_qm_atoms','nquant illegal', 'Need 0 < nquant <= natom')
   end if

   ! Check 2 - loop over nquant atoms and check it is > 1 and <= natom and it is unique
   do icount1=1,nquant
      !check atom number is legal
      iatom = iqmatoms(icount1)
      if ( (iatom > 0) .AND. (iatom <= natom)  ) then
        !QM atom ID is valid - check it is unique
        do icount2=(icount1+1),nquant
          if ( iatom == iqmatoms(icount2) ) then
            write (6,'(" QM ATOM VALIDATION: qm atom ",i8," is not unique.")') iatom
            call sander_bomb('validate_qm_atoms',&
                 'QM atoms specified with iqmatoms do not form a unique set.', &
                 'Require a unique list of qm atom numbers')
          end if
        end do
      else
        !it is not legal
        write (6,'(" QM ATOM VALIDATION: iqmatom ID number of ",i8," is not valid.")') iatom
        call sander_bomb('validate_qm_atoms','invalid QM atom ID', 'Need 0 < QM atom ID <= natom')
      end if
   end do

end subroutine validate_qm_atoms

!END SUBROUTINES
end module qmmm_module

