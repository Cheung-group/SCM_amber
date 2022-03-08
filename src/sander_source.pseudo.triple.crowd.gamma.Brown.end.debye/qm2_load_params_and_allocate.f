! <compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"
subroutine qm2_load_params_and_allocate()

! Written by: Ross Walker (TSRI, 2005)

! This routine should be called before running any qm2 routine calculations.
! It is responsible for filling the parameter arrays with the designated 
! parameters for the method chosen.

! All parameters are loaded into qm2_params structure.

! Note, allocation is done in this routine for all pointers in the qm2_params structure.
! Deallocation is done by deallocate qmmm routine.

! Parameter definitions:
!     Common params for all methods:
!          core(nquant_nlink) - The core charge on each atom
!  natomic_orbs(nquant_nlink) - Number of atomic orbitals on each atom.
!  orb_loc(2,nquant_nlink) - (1,x) = position of first orbital on atom x. 2,x = last orbital on x.
!  heat_of_form(nquant_nlink) - Gas Phase heat of formation for atom.

      use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, AM1, PM3, MNDO, PDDGPM3, PDDGMNDO, &
                              PM3CARB1, DFTB, RM1, nelements, MAX_VALENCE_ORBITALS, qmmm_mpi, qmmm_scratch
      use constants, only : EV_TO_KCAL, AU_TO_EV
      implicit none

!include the file containing all of the parameter constants.

!Locals
    _REAL_ :: pdiag_guess1, pdiag_guess2, pddg_zaf, pddg_zbf
    _REAL_ :: gssc, gspc, gppc, gp2c, hspc, elec_eng
    _REAL_ :: exponent_temp1, base_temp1, exponent_temp2, base_temp2
    _REAL_ :: HSP1_temp, HSP2_temp, DD1_temp, DD2_temp, DD_diff, DD3_temp, hpp
    _REAL_ :: aloc_real_scr
    integer :: ilaenv !is a function
    integer :: ilaenv_blocksize
    integer :: iostmp, ioptmp, int_dummy
    integer :: i, j, iqm_atomic, n_atomic_orb, first_orb, last_orb
    integer :: nheavy_atoms, nlight_atoms, nelectrons, nopen
    integer :: ier=0
    integer, allocatable :: reference_index(:) !Local, deallocated at end of routine

#ifdef MPI
#include "mpif.h"
    integer :: jstart, jend, ia, ib, ja, jb, status
    integer :: loop_extent_begin, loop_extent_end, loop_extent, loop_count
    integer, dimension(:), allocatable :: gather_array !Allocated and deallocated in this routine.
#endif

#include "qm2_parameters.h"

!Memory allocation - Needs cleaning up for DFTB at some point.
      allocate ( reference_index(qmmm_struct%qm_ntypes), stat = ier ) !Local
      REQUIRE(ier == 0)

      allocate ( qm2_params%core_chg(qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0) 
      allocate ( qm2_params%natomic_orbs(qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0) 
      allocate ( qm2_params%orb_loc(2,qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0) 
      allocate ( qm2_params%onec2elec_params(5,qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0)
      allocate ( qm2_params%multip_2c_elec_params(5,qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0)
      allocate ( qm2_params%cc_exp_params(qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0)
      allocate (qm2_params%orb_elec_ke(2,qmmm_struct%nquant_nlink), stat=ier )
      REQUIRE(ier == 0)
      allocate ( qm2_params%betasas(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
      REQUIRE(ier == 0) 
      allocate ( qm2_params%betasap(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
      REQUIRE(ier == 0) 
      allocate ( qm2_params%betapap(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
      REQUIRE(ier == 0) 
      if (qmmm_nml%qmtheory==AM1 .OR. qmmm_nml%qmtheory==PM3 .OR. qmmm_nml%qmtheory==PDDGPM3 .OR. &
          qmmm_nml%qmtheory==PM3CARB1 .OR. qmmm_nml%qmtheory==RM1 ) then
         allocate ( qm2_params%FN1(4,qmmm_struct%qm_ntypes), stat=ier ) 
         REQUIRE(ier == 0) 
         allocate ( qm2_params%FN2(4,qmmm_struct%qm_ntypes), stat=ier ) 
         REQUIRE(ier == 0) 
         allocate ( qm2_params%FN3(4,qmmm_struct%qm_ntypes), stat=ier ) 
         REQUIRE(ier == 0) 
         allocate ( qm2_params%NUM_FN(qmmm_struct%qm_ntypes), stat=ier ) 
         REQUIRE(ier == 0) 
      end if
      if (qmmm_nml%qmtheory==PDDGPM3 .OR. qmmm_nml%qmtheory==PDDGMNDO) then
         allocate ( qm2_params%pddge1(qmmm_struct%nquant_nlink), stat=ier )
         REQUIRE(ier == 0)
         allocate ( qm2_params%pddge2(qmmm_struct%nquant_nlink), stat=ier )
         REQUIRE(ier == 0)
         allocate ( qm2_params%pddg_term1(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
         REQUIRE(ier == 0)
         allocate ( qm2_params%pddg_term2(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
         REQUIRE(ier == 0)
         allocate ( qm2_params%pddg_term3(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
         REQUIRE(ier == 0)
         allocate ( qm2_params%pddg_term4(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier )
         REQUIRE(ier == 0)
      end if

!Note: s_orb_exp and p_orb_exp are only needed on the first call to fill 
!qm2_setup_orb_exp so after they are used in qm2_setup_orb_exp they are deallocated.
      allocate (qm2_params%s_orb_exp_by_type(nelements), stat=ier )
      REQUIRE(ier == 0) 
      allocate (qm2_params%p_orb_exp_by_type(nelements), stat=ier )
      REQUIRE(ier == 0) 
!Zero the total heat of formation before calculating it.
      qm2_params%tot_heat_form = 0.0D0
! We start by loading in the parameters that are common to all the semi-empirical methods.
!---------- COMMON PARAMS ---------------
      qm2_struct%norbs=0
      nheavy_atoms=0
      nelectrons=-qmmm_nml%qmcharge
      do i=1,qmmm_struct%nquant_nlink
         iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
         qm2_params%core_chg(i)=dble(core_chg(iqm_atomic))
         nelectrons=nelectrons+core_chg(iqm_atomic)
         n_atomic_orb=natomic_orbs(iqm_atomic)
         if (n_atomic_orb>1) nheavy_atoms=nheavy_atoms+1
         ! Check we don't bust any static arrays
         if (n_atomic_orb > MAX_VALENCE_ORBITALS) then
            write (6,*) 'n_atomic_orb of ',n_atomic_orb,' exceeds max_valence_orbitals of MAX_VALENCE_ORBITALS'
            call sander_bomb('qm2_load_params.f','exceeded max','Check qmmm_module.f and parameters.h')
         end if
         qm2_params%natomic_orbs(i)=n_atomic_orb
         qm2_params%orb_loc(1,i)=qm2_struct%norbs+1
         qm2_params%orb_loc(2,i)=qm2_struct%norbs+n_atomic_orb
         qm2_struct%norbs = qm2_struct%norbs+n_atomic_orb
         qm2_params%tot_heat_form=qm2_params%tot_heat_form+heat_of_form(iqm_atomic)
      end do !i=1,qmmm_struct%nquant_nlink

      !Work out how many 2 electron integrals there will be
      nlight_atoms=qmmm_struct%nquant_nlink-nheavy_atoms
      qm2_struct%n2el=50*nheavy_atoms*(nheavy_atoms-1)+10*nheavy_atoms*nlight_atoms+ &
                      ishft((nlight_atoms*(nlight_atoms-1)),-1)
      !QMMM e-repul memory depends on QM-MM pair list size so is
      !allocated later on and checked on every call.
      call qm2_allocate_qmqm_e_repul(qm2_struct%n2el)

      !Protect DUMB users from STUPID errors
      if (nelectrons > 2*qm2_struct%norbs) then
        if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
          write(6,'(''QMMM: ERROR-number of electrons: '',i5,'' is more'')') nelectrons
          write(6,'(''QMMM: than 2xnorbs of: '',i5)') qm2_struct%norbs
          write(6,'(''QMMM: Check qmcharge in qmmm namelist and rerun'')')                                            
          write(6,'(''QMMM: the calculation.'')')                                                                     
        end if                                                                                                    
        call mexit(6,1)
      end if  
      !Now we know the number of electrons work out how many closed and open shells there are.
      IF(qmmm_nml%spin==1 .OR. qmmm_nml%spin==3 .OR.qmmm_nml%spin==5)THEN
!        Make sure we have an even number of electrons
         if((nelectrons/2)*2 /= nelectrons) THEN
           if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
             write(6,'(''QMMM: System specified with odd number of electrons ('',i5,'')'')') nelectrons
             write(6,'(''QMMM: but odd spin ('',i3,''). Correct error and re-run calculation.'')') qmmm_nml%spin
           end if
           call mexit(6,1)
         end if
      else if (qmmm_nml%spin==2 .OR. qmmm_nml%spin==4 .OR. qmmm_nml%spin==6) then
!        Make sure we have an odd number of electrons.`
         if((nelectrons/2)*2 == nelectrons) then
           if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
             write(6,'(''QMMM: System specified with even number of electrons ('',i5,'')'')') nelectrons
             write(6,'(''QMMM: but even spin ('',i3,''). Correct error and re-run calculation.'')') qmmm_nml%spin
           end if
           call mexit(6,1)
         end if
      end if

      if (qmmm_nml%spin==1) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) write (6,'(''QMMM: SINGLET STATE CALCULATION'')')
         NOPEN=0
      else if (qmmm_nml%spin==2) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) write (6,'(''QMMM: DOUBLET STATE CALCULATION'')')
         NOPEN=1
      else if (qmmm_nml%spin==3) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) WRITE(6,'(''QMMM: TRIPLET STATE CALCULATION'')')
         NOPEN=2
      else if (qmmm_nml%spin==4) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) WRITE(6,'(''QMMM: QUARTET STATE CALCULATION'')')
         NOPEN=3
      else if (qmmm_nml%spin==5) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) WRITE(6,'(''QMMM: QUINTET STATE CALCULATION'')')
         NOPEN=4
      else if (qmmm_nml%spin==6) then
         if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) WRITE(6,'(''QMMM: SEXTET STATE CALCULATION'')')
         NOPEN=5
      ENDIF                                                                                                       
      qm2_struct%nclosed=nelectrons/2          
      if( NOPEN > 0 ) then
         qm2_struct%nclosed=qm2_struct%nclosed-nopen/2
         if(qm2_struct%nclosed+nopen>qm2_struct%norbs) then
            if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
              write(6,'(''QMMM: Number of doubly filled ('',i3,'') plus'')') qm2_struct%nclosed
              write(6,'(''QMMM: number of partly filled ('',i3,'') levels'')') nopen
              write(6,'(''QMMM: is greater than the total number of orbitals ('',i5,'').'')') qm2_struct%norbs
              write(6,'(''QMMM: Fix problem and re-run calculation.'')')
            end if
            stop
         end if
         if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
            write(6,'(''QMMM: ROHF CALCULATION'')')
            write(6,'(''QMMM: THERE ARE'',I3,'' DOUBLY FILLED LEVELS'')') qm2_struct%nclosed
            write(6,'(''QMMM: AND'')')
            write(6,'(''QMMM:          '',i3,'' SINGLY OCCUPIED LEVELS'')') NOPEN
        end if
      else
        if(qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call)  then
            write(6,'(''QMMM: RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS ='',I3)') qm2_struct%nclosed
        end if
      end if
      qm2_struct%nopenclosed=nopen+qm2_struct%nclosed

      !Allocate things that depend on norbs:

      allocate (qm2_params%pascal_tri1(qm2_struct%norbs), stat=ier )
      REQUIRE(ier == 0)
      allocate (qm2_params%pascal_tri2(qm2_struct%norbs), stat=ier )
      REQUIRE(ier == 0)

      qm2_struct%matsize = ishft(qm2_struct%norbs*(qm2_struct%norbs+1),-1) !ishift(x,-1) = integer divide by 2
      allocate ( qm2_struct%den_matrix(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
      allocate ( qm2_struct%old_den_matrix(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
      !zero the entire density matrix on the first call
      qm2_struct%den_matrix = 0.0d0; qm2_struct%old_den_matrix = 0.0d0;

      allocate ( qm2_struct%old2_density(qm2_struct%norbs), stat=ier )
      REQUIRE ( ier == 0 ) !Used by qm2_cnvg as workspace.

      if (qmmm_nml%density_predict == 1) then
        !We are using Niklasson et al. density matrix prediction.
        allocate ( qm2_struct%md_den_mat_guess1(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
        allocate ( qm2_struct%md_den_mat_guess2(qm2_struct%matsize), stat=ier )
        REQUIRE ( ier == 0 )
      end if

      allocate ( qm2_struct%fock_matrix(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
      allocate ( qm2_struct%hmatrix(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
#ifdef MPI
# ifndef USE_MPI_IN_PLACE
      !Allocate a temporary array for doing the reduce of P and F etc. Only needed if
      !we can't do things using MPI_IN_PLACE
      allocate ( qmmm_scratch%matsize_red_scratch(qm2_struct%matsize), stat=ier )
      REQUIRE ( ier == 0 )
# endif
#endif
      allocate (qm2_struct%fock2_ptot2(16,qmmm_struct%nquant_nlink),stat=ier)
      REQUIRE(ier==0)

      !set up array of lower half triangle indicies (Pascal's triangle)
      do I=1,qm2_struct%norbs
         qm2_params%pascal_tri1(I)=ishft((I*(I-1)),-1)
         qm2_params%pascal_tri2(I)=qm2_params%pascal_tri1(I)+I
      end do

      !Fill the diagonal of the density matrix with the first guess:

      pdiag_guess1=dble(qmmm_nml%qmcharge)/(qm2_struct%norbs+1.D-10)
      do i=1,qmmm_struct%nquant_nlink
         first_orb=qm2_params%orb_loc(1,i)
         last_orb=qm2_params%orb_loc(2,i)
         pdiag_guess2=(qm2_params%core_chg(i)/ &
            (qm2_params%natomic_orbs(i)+1.D-10))-pdiag_guess1
         do j=first_orb,last_orb
           qm2_struct%den_matrix(qm2_params%pascal_tri2(j))=pdiag_guess2
           qm2_struct%old_den_matrix(qm2_params%pascal_tri2(j))=pdiag_guess2
         end do
      end do

      if ( qmmm_nml%density_predict == 1) then
        !We are using Pguess(t) = 2Pconv(t-1) - Pguess(t-2)
        !in this case for an initial start we set
        !den_matrix = 0.5 initial guess
        !md_den_mat_guess1 = initial guess
        !md_den_mat_guess2 = 0.0d0
        ! then
        ! on step 1 we get: den_matrix = 2.0d0 * den_matrix - md_den_mat_guess2 (0,0d0)
        !                              = initial guess
        ! on step 2 we get: md_den_mat_guess2 = md_den_mat_guess1 = initial_guess
        ! den_matrix = 2.0d0 * den_matrix - md_den_mat_guess2 (initial_guess)

        qm2_struct%md_den_mat_guess2(1:qm2_struct%matsize) = 0.0d0
        qm2_struct%md_den_mat_guess1(1:qm2_struct%matsize) = 0.0d0
        qm2_struct%den_matrix(1:qm2_struct%matsize) = qm2_struct%den_matrix(1:qm2_struct%matsize) * 0.5d0

      end if

!----------------------------------------
!
! Now we fill up the data depending on the method we are using
!
!--------------------------------------
!           MNDO PARAMS               *
!--------------------------------------
      if (qmmm_nml%qmtheory == MNDO) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in MNDO
          if (.NOT. element_supported_mndo(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no MNDO parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM MNDO NOT AVAILABLE FOR THIS ATOM')
          end if

          
          !----------------------------------------
          ! Calculate parameters that are actually 
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
          !            + HSP*HSPC
          iostmp = ios(iqm_atomic)
          ioptmp = iop(iqm_atomic)
          gssc = dble(max(iostmp-1,0))
          gspc = dble(iostmp * ioptmp)
          gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                 + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          hspc = dble(-ioptmp)
          elec_eng = uss_mndo(iqm_atomic)*iostmp + upp_mndo(iqm_atomic)*ioptmp &
                   + gss_mndo(iqm_atomic)*gssc + gsp_mndo(iqm_atomic)*gspc &
                   + gpp_mndo(iqm_atomic)*gppc + gp2_mndo(iqm_atomic)*gp2c &
                   + hsp_mndo(iqm_atomic)*hspc

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_mndo(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_mndo(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_mndo(iqm_atomic)*p_orb_exp_mndo(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_mndo(iqm_atomic) + p_orb_exp_mndo(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_mndo(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_mndo(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_mndo(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            !AD
            dd1_temp = (HSP_mndo(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_mndo(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_mndo(iqm_atomic)-gp2_mndo(iqm_atomic)) 
            hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                 !match mopac results for Chlorine and potentially other elements.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if
      
          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

          qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_mndo(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_mndo(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_mndo(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_mndo(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_mndo(iqm_atomic)
          qm2_params%cc_exp_params(i) = alp_mndo(iqm_atomic)
          qm2_params%orb_elec_ke(1,i) = uss_mndo(iqm_atomic)
          qm2_params%orb_elec_ke(2,i) = upp_mndo(iqm_atomic)
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_mndo(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_mndo(i)
        end do

! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = betas_mndo(qmmm_struct%qm_type_id(i))+betas_mndo(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = betas_mndo(qmmm_struct%qm_type_id(i))+betap_mndo(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = betap_mndo(qmmm_struct%qm_type_id(i))+betap_mndo(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = mndo_ref_index(qmmm_struct%qm_type_id(i))
        end do

!--------------------------------------
!       END  MNDO PARAMS              *
!--------------------------------------

!--------------------------------------
!           AM1 PARAMS                *
!--------------------------------------
      else if (qmmm_nml%qmtheory == AM1) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in AM1
          if (.NOT. element_supported_am1(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no AM1 parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM AM1 NOT AVAILABLE FOR THIS ATOM')
          end if

          !----------------------------------------
          ! Calculate parameters that are actually
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
          !            + HSP*HSPC
          iostmp = ios(iqm_atomic)
          ioptmp = iop(iqm_atomic)
          gssc = dble(max(iostmp-1,0))
          gspc = dble(iostmp * ioptmp)
          gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                 + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          hspc = dble(-ioptmp)
          elec_eng = uss_am1(iqm_atomic)*iostmp + upp_am1(iqm_atomic)*ioptmp &
                   + gss_am1(iqm_atomic)*gssc + gsp_am1(iqm_atomic)*gspc &
                   + gpp_am1(iqm_atomic)*gppc + gp2_am1(iqm_atomic)*gp2c &
                   + hsp_am1(iqm_atomic)*hspc

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_am1(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_am1(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_am1(iqm_atomic)*p_orb_exp_am1(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_am1(iqm_atomic) + p_orb_exp_am1(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_am1(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_am1(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_am1(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            dd1_temp = (HSP_am1(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_am1(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_am1(iqm_atomic)-gp2_am1(iqm_atomic)) 
            hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                 !match mopac results for Chlorine and potentially other elements.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if

          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

! Next do the electronic energy - add it to the total heat of formation energy.
          qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_am1(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_am1(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_am1(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_am1(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_am1(iqm_atomic)
          qm2_params%cc_exp_params(i) = alp_am1(iqm_atomic)
          qm2_params%orb_elec_ke(1,i) = uss_am1(iqm_atomic)
          qm2_params%orb_elec_ke(2,i) = upp_am1(iqm_atomic)
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_am1(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_am1(i)
        end do

! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          qm2_params%NUM_FN(i) = NUM_FN_am1(qmmm_struct%qm_type_id(i))
          do j=1,4
             qm2_params%FN1(j,i) = FN1_am1(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN2(j,i) = FN2_am1(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN3(j,i) = FN3_am1(j,qmmm_struct%qm_type_id(i))
          end do
        end do

        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = betas_am1(qmmm_struct%qm_type_id(i))+betas_am1(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = betas_am1(qmmm_struct%qm_type_id(i))+betap_am1(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = betap_am1(qmmm_struct%qm_type_id(i))+betap_am1(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = am1_ref_index(qmmm_struct%qm_type_id(i))
        end do
!--------------------------------------
!        END  AM1 PARAMS              *
!--------------------------------------

!--------------------------------------
!      PM3 AND PM3CARB1 PARAMS        *
!--------------------------------------
      else if (qmmm_nml%qmtheory == PM3 .OR. qmmm_nml%qmtheory == PM3CARB1 ) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in PM3
          if (.NOT. element_supported_pm3(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no PM3 parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM PM3 NOT AVAILABLE FOR THIS ATOM')
          end if
          !----------------------------------------
          ! Calculate parameters that are actually
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
          !            + HSP*HSPC
          iostmp = ios(iqm_atomic)
          ioptmp = iop(iqm_atomic)
          gssc = dble(max(iostmp-1,0))
          gspc = dble(iostmp * ioptmp)
          gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                 + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          hspc = dble(-ioptmp)
          elec_eng = uss_pm3(iqm_atomic)*iostmp + upp_pm3(iqm_atomic)*ioptmp &
                   + gss_pm3(iqm_atomic)*gssc + gsp_pm3(iqm_atomic)*gspc &
                   + gpp_pm3(iqm_atomic)*gppc + gp2_pm3(iqm_atomic)*gp2c &
                   + hsp_pm3(iqm_atomic)*hspc

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_pm3(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_pm3(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_pm3(iqm_atomic)*p_orb_exp_pm3(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_pm3(iqm_atomic) + p_orb_exp_pm3(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pm3(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_pm3(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pm3(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            dd1_temp = (HSP_pm3(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_pm3(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_pm3(iqm_atomic)-gp2_pm3(iqm_atomic)) 
            hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                 !match mopac results for Chlorine and potentially other elements.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if

          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pm3(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_pm3(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pm3(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pm3(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pm3(iqm_atomic)
          if (qmmm_nml%qmtheory == PM3CARB1 .AND. (iqm_atomic==1 .OR. iqm_atomic==8)) then
            !Load the PM3CARB1 versions of O and H params in place of the default PM3 params
            elec_eng = uss_pm3carb1(iqm_atomic)*iostmp + upp_pm3carb1(iqm_atomic)*ioptmp &
                     + gss_pm3(iqm_atomic)*gssc + gsp_pm3(iqm_atomic)*gspc &
                     + gpp_pm3(iqm_atomic)*gppc + gp2_pm3(iqm_atomic)*gp2c &
                     + hsp_pm3(iqm_atomic)*hspc
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%cc_exp_params(i) = alp_pm3carb1(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pm3carb1(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pm3carb1(iqm_atomic)
          else
            qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
            qm2_params%cc_exp_params(i) = alp_pm3(iqm_atomic)
            qm2_params%orb_elec_ke(1,i) = uss_pm3(iqm_atomic)
            qm2_params%orb_elec_ke(2,i) = upp_pm3(iqm_atomic)
          end if
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pm3(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pm3(i)
        end do

! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          qm2_params%NUM_FN(i) = NUM_FN_pm3(qmmm_struct%qm_type_id(i))
          do j=1,4
             qm2_params%FN1(j,i) = FN1_pm3(j,qmmm_struct%qm_type_id(i))   
             qm2_params%FN2(j,i) = FN2_pm3(j,qmmm_struct%qm_type_id(i))   
             qm2_params%FN3(j,i) = FN3_pm3(j,qmmm_struct%qm_type_id(i))   
          end do
        end do

!RCW: PM3CARB1 update - we need to make sure we use the correct parameters here for hydrogen and oxygen atoms
        !Simple solution is to replace the data in the betas_pm3 and betap_pm3 arrays with the PM3CARB1 values
        !This avoids having to use complex if statements in the loop below.
        if (qmmm_nml%qmtheory == PM3CARB1) then
           !Replace PM3 betas and betap params for O and H with PM3CARB1 values BEFORE they are copied
           !into the working array.
           betas_pm3(1) = betas_pm3carb1(1)
           betas_pm3(8) = betas_pm3carb1(8)
           betap_pm3(1) = betap_pm3carb1(1)
           betap_pm3(8) = betap_pm3carb1(8)
        end if
        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = betas_pm3(qmmm_struct%qm_type_id(i))+betas_pm3(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = betas_pm3(qmmm_struct%qm_type_id(i))+betap_pm3(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = betap_pm3(qmmm_struct%qm_type_id(i))+betap_pm3(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = pm3_ref_index(qmmm_struct%qm_type_id(i))
        end do

        !Replace with PM3CARB1 reference if necessary
        if (qmmm_nml%qmtheory == PM3CARB1) then
          do i = 1, qmmm_struct%qm_ntypes
            if (qmmm_struct%qm_type_id(i) == 1 .OR. qmmm_struct%qm_type_id(i) == 8) then
              reference_index(i) = pm3carb1_ref_index(qmmm_struct%qm_type_id(i))
            end if
          end do
        end if

!--------------------------------------
!    END PM3 AND PM3 CARB1 PARAMS     *
!--------------------------------------

!--------------------------------------
!         PDDG/PM3 PARAMS             *
!--------------------------------------
      else if (qmmm_nml%qmtheory == PDDGPM3) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in PM3
          if (.NOT. element_supported_pddgpm3(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no PDDG-PM3 parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM PDDG-PM3 NOT AVAILABLE FOR THIS ATOM')
          end if
          !----------------------------------------
          ! Calculate parameters that are actually
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          ! In PM3/PDDG the electronic energy is treated as an independent parameter
          ! and so we have to have it explicitly defined in the parameter file and
          ! cannot just calculate it.

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_pddgpm3(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_pddgpm3(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_pddgpm3(iqm_atomic)*p_orb_exp_pddgpm3(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_pddgpm3(iqm_atomic) + p_orb_exp_pddgpm3(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pddgpm3(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_pddgpm3(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pddgpm3(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            dd1_temp = (HSP_pddgpm3(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_pddgpm3(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_pddgpm3(iqm_atomic)-gp2_pddgpm3(iqm_atomic)) 
!            hpp = max(0.1d0,hpp) - It seems that PM3/PDDG does not have this maximum
!                                   adding it changes Cl results.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if

          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

          qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng_pddgpm3(iqm_atomic)*EV_TO_KCAL)
          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pddgpm3(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_pddgpm3(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pddgpm3(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pddgpm3(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pddgpm3(iqm_atomic)
          qm2_params%cc_exp_params(i) = alp_pddgpm3(iqm_atomic)
          qm2_params%orb_elec_ke(1,i) = uss_pddgpm3(iqm_atomic)
          qm2_params%orb_elec_ke(2,i) = upp_pddgpm3(iqm_atomic)
          qm2_params%pddge1(i) = pddge1_pm3(iqm_atomic)
          qm2_params%pddge2(i) = pddge2_pm3(iqm_atomic)
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pddgpm3(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pddgpm3(i)
        end do

! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          qm2_params%NUM_FN(i) = NUM_FN_pddgpm3(qmmm_struct%qm_type_id(i))
          do j=1,4
             qm2_params%FN1(j,i) = FN1_pddgpm3(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN2(j,i) = FN2_pddgpm3(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN3(j,i) = FN3_pddgpm3(j,qmmm_struct%qm_type_id(i))
          end do
        end do

        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = betas_pddgpm3(qmmm_struct%qm_type_id(i)) &
               +betas_pddgpm3(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = betas_pddgpm3(qmmm_struct%qm_type_id(i)) &
               +betap_pddgpm3(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = betap_pddgpm3(qmmm_struct%qm_type_id(i)) &
               +betap_pddgpm3(qmmm_struct%qm_type_id(j))
            pddg_zaf = dble(core_chg(qmmm_struct%qm_type_id(i))) &
               /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
               dble(core_chg(qmmm_struct%qm_type_id(j))))
            pddg_zbf = dble(core_chg(qmmm_struct%qm_type_id(j))) &
               /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
               dble(core_chg(qmmm_struct%qm_type_id(j))))
            qm2_params%pddg_term1(i,j) = &
               pddg_zaf*pddgc1_pm3(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc1_pm3(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term2(i,j) = &
               pddg_zaf*pddgc1_pm3(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc2_pm3(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term3(i,j) = &
               pddg_zaf*pddgc2_pm3(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc1_pm3(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term4(i,j) = &
               pddg_zaf*pddgc2_pm3(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc2_pm3(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = pddgpm3_ref_index(qmmm_struct%qm_type_id(i))
        end do
!--------------------------------------
!          END PDDG/PM3 PARAMS        *
!--------------------------------------
!--------------------------------------
!           PDDG/MNDO PARAMS          *
!--------------------------------------
      elseif (qmmm_nml%qmtheory == PDDGMNDO) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in MNDO
          if (.NOT. element_supported_pddgmndo(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no PDDG-MNDO parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM PDDG-MNDO NOT AVAILABLE FOR THIS ATOM')
          end if

          !----------------------------------------
          ! Calculate parameters that are actually
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          ! In PM3/PDDG the electronic energy is treated as an independent parameter
          ! and so we have to have it explicitly defined in the parameter file and
          ! cannot just calculate it.

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_pddgmndo(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_pddgmndo(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_pddgmndo(iqm_atomic)*p_orb_exp_pddgmndo(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_pddgmndo(iqm_atomic) + p_orb_exp_pddgmndo(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_pddgmndo(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_pddgmndo(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_pddgmndo(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            dd1_temp = (HSP_pddgmndo(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_pddgmndo(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_pddgmndo(iqm_atomic)-gp2_pddgmndo(iqm_atomic)) 
!            hpp = max(0.1d0,hpp) - It seems that like PM3/PDDG, MNDO/PDDG does not have this maximum.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if

          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

! Next do the electronic energy - add it to the total heat of formation energy.
          qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng_pddgmndo(iqm_atomic)*EV_TO_KCAL)
          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_pddgmndo(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_pddgmndo(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_pddgmndo(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_pddgmndo(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_pddgmndo(iqm_atomic)
          qm2_params%cc_exp_params(i) = alp_pddgmndo(iqm_atomic)
          qm2_params%orb_elec_ke(1,i) = uss_pddgmndo(iqm_atomic)
          qm2_params%orb_elec_ke(2,i) = upp_pddgmndo(iqm_atomic)
          qm2_params%pddge1(i) = pddge1_mndo(iqm_atomic)
          qm2_params%pddge2(i) = pddge2_mndo(iqm_atomic)
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_pddgmndo(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_pddgmndo(i)
        end do
                                                                                                                                                             
! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = &
               betas_pddgmndo(qmmm_struct%qm_type_id(i)) &
               +betas_pddgmndo(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = &
               betas_pddgmndo(qmmm_struct%qm_type_id(i)) &
               +betap_pddgmndo(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = &
               betap_pddgmndo(qmmm_struct%qm_type_id(i)) &
               +betap_pddgmndo(qmmm_struct%qm_type_id(j))
            pddg_zaf = dble(core_chg(qmmm_struct%qm_type_id(i))) &
               /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
               dble(core_chg(qmmm_struct%qm_type_id(j))))
            pddg_zbf = dble(core_chg(qmmm_struct%qm_type_id(j))) &
               /(dble(core_chg(qmmm_struct%qm_type_id(i)))+ &
               dble(core_chg(qmmm_struct%qm_type_id(j))))
            qm2_params%pddg_term1(i,j) = &
               pddg_zaf*pddgc1_mndo(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc1_mndo(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term2(i,j) = &
               pddg_zaf*pddgc1_mndo(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc2_mndo(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term3(i,j) = &
               pddg_zaf*pddgc2_mndo(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc1_mndo(qmmm_struct%qm_type_id(j))
            qm2_params%pddg_term4(i,j) = &
               pddg_zaf*pddgc2_mndo(qmmm_struct%qm_type_id(i)) &
               +pddg_zbf*pddgc2_mndo(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = pddgmndo_ref_index(qmmm_struct%qm_type_id(i))
        end do
!--------------------------------------
!     END PDDG/MNDO PARAMS            *
!--------------------------------------

!--------------------------------------
!           RM1 PARAMS                *
!--------------------------------------
      else if (qmmm_nml%qmtheory == RM1) then
        do i = 1,qmmm_struct%nquant_nlink
          iqm_atomic=qmmm_struct%iqm_atomic_numbers(i)
! Check that parameters exist for this element in RM1
          if (.NOT. element_supported_rm1(iqm_atomic)) then
             write(6,'("QMMM: Atom number: ",i6," has atomic number ",i4,".")') i, iqm_atomic
             write(6,'("QMMM: There are no RM1 parameters for this element. Sorry.")')
             call sander_bomb('qm2_load_params','UNSUPPORTED ELEMENT','QM RM1 NOT AVAILABLE FOR THIS ATOM')
          end if

          !----------------------------------------
          ! Calculate parameters that are actually
          ! derived from other parameters.
          !----------------------------------------

          !1) Electronic Energy (EISOL)
          !   elec_eng = USS*IOS + UPP*IOP + UDD*IOD + GSS*GSSC + GPP*GPPC + GSP*GSPC + GP2*GP2C
          !            + HSP*HSPC
          iostmp = ios(iqm_atomic)
          ioptmp = iop(iqm_atomic)
          gssc = dble(max(iostmp-1,0))
          gspc = dble(iostmp * ioptmp)
          gp2c = dble((ioptmp * (ioptmp - 1))/2) &
                 + 0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          gppc = -0.5d0*dble(min(ioptmp,6-ioptmp)*(min(ioptmp,6-ioptmp)-1)/2)
          hspc = dble(-ioptmp)
          elec_eng = uss_rm1(iqm_atomic)*iostmp + upp_rm1(iqm_atomic)*ioptmp &
                   + gss_rm1(iqm_atomic)*gssc + gsp_rm1(iqm_atomic)*gspc &
                   + gpp_rm1(iqm_atomic)*gppc + gp2_rm1(iqm_atomic)*gp2c &
                   + hsp_rm1(iqm_atomic)*hspc

          !2) multip_2c_elec_params(1-5,i) (DD,QQ,AM,AD,AQ)
          !   DD = (( (4.0d0*s_orb_exp*p_orb_exp)**(nsshell+0.5d0) ) * (2.0d0*nsshell + 1)) &
          !       / (( (s_orb_exp + p_orb_exp)**(2.0d0*nsshell + 2.0d0) ) * sqrt(3.0d0))
          !
          !   QQ = sqrt((4.0d0*nsshell**2+6.0d0*nsshell+2.0d0)/20.0d0)/p_orb_exp
          !   AM = GSS/AU_TO_EV
          if (p_orb_exp_rm1(iqm_atomic) .ne. 0.0d0 .or. &
              s_orb_exp_rm1(iqm_atomic) .ne. 0.0d0) then
              exponent_temp1 = nsshell(iqm_atomic)+0.5d0
              base_temp1 = 4.0d0*s_orb_exp_rm1(iqm_atomic)*p_orb_exp_rm1(iqm_atomic)
              exponent_temp2 = 2.0d0*nsshell(iqm_atomic) + 2.0d0
              base_temp2 = s_orb_exp_rm1(iqm_atomic) + p_orb_exp_rm1(iqm_atomic)
              qm2_params%multip_2c_elec_params(1,i) = ((base_temp1**exponent_temp1)*(2.0d0*nsshell(iqm_atomic) + 1.0d0)) &
                                                      / ((base_temp2**exponent_temp2) * sqrt(3.0d0))
              qm2_params%multip_2c_elec_params(2,i) = &
                     sqrt((4.0d0*nsshell(iqm_atomic)**2+6.0d0*nsshell(iqm_atomic) &
                    +2.0d0)/20.0d0)/p_orb_exp_rm1(iqm_atomic)
          else
            qm2_params%multip_2c_elec_params(1,i)= 0.0d0
            qm2_params%multip_2c_elec_params(2,i)= 0.0d0
          end if
          if (GSS_rm1(iqm_atomic) .ne. 0.0d0 ) then
            qm2_params%multip_2c_elec_params(3,i) = (0.5d0*AU_TO_EV)/GSS_rm1(iqm_atomic) !AM
          else
            qm2_params%multip_2c_elec_params(3,i) = 0.0d0
          end if
          ! Calculation of AD and AQ
          if (iqm_atomic == 1) then
            qm2_params%multip_2c_elec_params(4,i) = qm2_params%multip_2c_elec_params(3,i) !AD for H
            qm2_params%multip_2c_elec_params(5,i) = qm2_params%multip_2c_elec_params(3,i) !AQ for H
          else
            dd1_temp = (HSP_rm1(iqm_atomic) &
                      /(AU_TO_EV*qm2_params%multip_2c_elec_params(1,i)**2))**(1.D0/3.D0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.5D0*dd1_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd1_temp**2)
              hsp2_temp = 0.5D0*dd2_temp &
                        - 0.5D0/sqrt(4.D0*qm2_params%multip_2c_elec_params(1,i)**2+1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(HSP_rm1(iqm_atomic)/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(4,i) = 0.5d0/dd2_temp
            !END AD
            !AQ
            hpp = 0.5D0*(gpp_rm1(iqm_atomic)-gp2_rm1(iqm_atomic)) 
            hpp = max(0.1d0,hpp) !I have no idea where this max comes from but it is required to
                                 !match mopac results for Chlorine and potentially other elements.
            dd1_temp = (16.0d0*hpp &
                       /(AU_TO_EV*48.0d0*qm2_params%multip_2c_elec_params(2,i)**4))**(1.0d0/5.0d0)
            dd2_temp = dd1_temp + 0.04d0
            do j = 1, 5
              dd_diff = dd2_temp - dd1_temp
              hsp1_temp = 0.25d0*dd1_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd1_temp**2)
              hsp2_temp = 0.25d0*dd2_temp - 0.5d0/sqrt(4.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                          + 1.0d0/dd2_temp**2) + 0.25d0/sqrt(8.0d0*qm2_params%multip_2c_elec_params(2,i)**2 &
                         + 1.0d0/dd2_temp**2)
              if (abs(hsp2_temp - hsp1_temp) < 1.0d-25) exit
              dd3_temp = dd1_temp + dd_diff*(hpp/AU_TO_EV-hsp1_temp) &
                       / (hsp2_temp - hsp1_temp)
              dd1_temp = dd2_temp
              dd2_temp = dd3_temp
            end do
            qm2_params%multip_2c_elec_params(5,i) = 0.5d0/dd2_temp
            !END AQ
          end if

          !----------------------------------------
          ! End calculation of derived parameters.
          !----------------------------------------

! Next do the electronic energy - add it to the total heat of formation energy.
          qm2_params%tot_heat_form = qm2_params%tot_heat_form-(elec_eng*EV_TO_KCAL)
          qm2_params%onec2elec_params(1,i) = 0.5D0*GSS_rm1(iqm_atomic)
          qm2_params%onec2elec_params(2,i) = GSP_rm1(iqm_atomic)
          qm2_params%onec2elec_params(3,i) = 0.5d0*GPP_rm1(iqm_atomic)
          qm2_params%onec2elec_params(4,i) = 1.25d0*GP2_rm1(iqm_atomic)
          qm2_params%onec2elec_params(5,i) = 0.5d0*HSP_rm1(iqm_atomic)
          qm2_params%cc_exp_params(i) = alp_rm1(iqm_atomic)
          qm2_params%orb_elec_ke(1,i) = uss_rm1(iqm_atomic)
          qm2_params%orb_elec_ke(2,i) = upp_rm1(iqm_atomic)
        end do
        !Loop over elements
        do i = 1,nelements
! Next get the Slater orbital expansion coefficients
           qm2_params%s_orb_exp_by_type(i) = s_orb_exp_rm1(i)
           qm2_params%p_orb_exp_by_type(i) = p_orb_exp_rm1(i)
        end do

! Precompute some parameters to save time later
! RCW: By rebasing the atomic numbers of each atoms as types we reduce the overall
! size of the arrays. While this does not save us much memory it greatly increases
! the chance of cache hits.
        do i=1,qmmm_struct%qm_ntypes
          qm2_params%NUM_FN(i) = NUM_FN_rm1(qmmm_struct%qm_type_id(i))
          do j=1,4
             qm2_params%FN1(j,i) = FN1_rm1(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN2(j,i) = FN2_rm1(j,qmmm_struct%qm_type_id(i))
             qm2_params%FN3(j,i) = FN3_rm1(j,qmmm_struct%qm_type_id(i))
          end do
        end do

        do i=1,qmmm_struct%qm_ntypes
          do j = 1,qmmm_struct%qm_ntypes
            qm2_params%betasas(i,j) = betas_rm1(qmmm_struct%qm_type_id(i))+betas_rm1(qmmm_struct%qm_type_id(j))
            qm2_params%betasap(i,j) = betas_rm1(qmmm_struct%qm_type_id(i))+betap_rm1(qmmm_struct%qm_type_id(j))
            qm2_params%betapap(i,j) = betap_rm1(qmmm_struct%qm_type_id(i))+betap_rm1(qmmm_struct%qm_type_id(j))
          end do
        end do

        do i = 1, qmmm_struct%qm_ntypes
          !Store local reference info here
          reference_index(i) = rm1_ref_index(qmmm_struct%qm_type_id(i))
        end do
!--------------------------------------
!        END  RM1 PARAMS              *
!--------------------------------------
      else if (qmmm_nml%qmtheory == DFTB) then
         call qm2_dftb_load_params
      else
        !UNKNOWN method - should never actually get this far but might as well call
        !sander bomb just in case.
        write (6,'("QMMM ERROR: Method ID: ",i5," is not supported.")') qmmm_nml%qmtheory
        call sander_bomb('qm2_load_params','UNSUPPORTED METHOD', &
                         'SELECTED LEVEL OF THEORY IS NOT AVAILABLE - PLEASE CHECK YOUR INPUT FILE')
      end if


      if (qmmm_nml%qmtheory /= DFTB) then
        !Now see if user wants an MM peptide torsion correction
        qm2_struct%n_peptide_links = 0
        if(qmmm_nml%peptide_corr) THEN
           if(qmmm_mpi%commqmmm_master) write (6,'(''QMMM: MOLECULAR MECHANICS CORRECTION APPLIED TO PEPTIDE LINKAGES'')')
           call qm2_identify_peptide_links(qm2_struct%n_peptide_links,qmmm_struct%qm_coords)
           if (qmmm_mpi%commqmmm_master) then
              write (6,'(''QMMM: '',i5,'' PEPTIDE LINKAGES HAVE BEEN FOUND:'')') qm2_struct%n_peptide_links
              do i=1,qm2_struct%n_peptide_links
                write(6,'(''QMMM:    '',i4,'' - '',i4,'' - '',i4,'' - '',i4)') qm2_struct%peptide_links(1,i), &
                      qm2_struct%peptide_links(2,i), qm2_struct%peptide_links(3,i), qm2_struct%peptide_links(4,i)
              end do
           end if
        end if
        call qm2_setup_orb_exp
        !Finally setup the STO-6G orbital expansions and allocate the memory required.
        !Setup the STO-6G orbital expansions and pre-calculate as many overlaps by type
        !as we can and store these in memory. This will help a lot with speed in the
        !energy and derivative code.

        !The last stage is to print out the references regarding which parameters are in use
        if (qmmm_mpi%commqmmm_master.and.qmmm_struct%qm_mm_first_call) then
          write(6,'(/,"| QMMM: Parameter sets in use:")')
          do i=1,qmmm_struct%qm_ntypes
            call qm_print_ref(reference_index(i), qmmm_struct%qm_type_id(i))
          end do
        end if
      end if !qmmm_nml%qmtheory /= DFTB

      deallocate(reference_index,stat=ier)
      REQUIRE(ier == 0)

!In Parallel calculate the offset into the two electron array for each thread.
!This depends on the number of light-light, light-heavy and heavy-heavy interactions
!that this thread will do.

!Simulate what my loop would be and work out what my ending offset would be.
#ifdef MPI
      loop_count = 0
      do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
        jstart = qmmm_mpi%nquant_nlink_jrange(1,i)
        jend = qmmm_mpi%nquant_nlink_jrange(2,i)
        ia = qm2_params%orb_loc(1,i)
        ib = qm2_params%orb_loc(2,i)
        do j = jstart, jend
          ja = qm2_params%orb_loc(1,j)
          jb = qm2_params%orb_loc(2,j)
          if (ib /= ia .and. ja/=jb) then
            !Heavy - Heavy
            loop_count = loop_count+100
          elseif (ia /= ib) then
            !Light - Heavy
            loop_count = loop_count+10
          elseif (ja /= jb) then
            !Heavy - Light
            loop_count = loop_count+10
          else
            !Light - Light
            loop_count = loop_count+1
          endif
        end do
      end do

      !At the end of this loop loop_count should be the starting value for the next cpu.
      !However, we need to add the offset from all the other cpus to this.
      !Each thread in turn passes it's total so far to the next thread.
      !This total becomes that threads offset
      allocate(gather_array(qmmm_mpi%numthreads),stat=ier)
      REQUIRE(ier==0)
      call mpi_allgather(loop_count,1,mpi_integer,gather_array,1,mpi_integer,qmmm_mpi%commqmmm,ier)

      !Now, our starting offset is the sum of all the previous thread's loop count
      qmmm_mpi%two_e_offset = 0
      do i = 1, qmmm_mpi%mytaskid
        qmmm_mpi%two_e_offset = qmmm_mpi%two_e_offset + gather_array(i)
      end do
      deallocate(gather_array,stat=ier)
      REQUIRE(ier==0)
#else
      qmmm_mpi%two_e_offset = 0
#endif

!#ifdef MPI
!     !NOT CURRENTLY USED
!     !Now we know the number of orbitals divide them up between cpus.
!     !allocate the memory for the jrange array
!
!     allocate(qmmm_mpi%norb_jrange(2,qm2_struct%norbs),stat=ier)
!     REQUIRE(ier==0)
!     loop_extent = qm2_struct%norbs*(qm2_struct%norbs-1)/2
!     mpi_division = (loop_extent+(qmmm_mpi%numthreads-1))/qmmm_mpi%numthreads
!     loop_extent_end = min(mpi_division*(qmmm_mpi%mytaskid+1),loop_extent)
!     loop_extent_begin = mpi_division*qmmm_mpi%mytaskid+1
!!loop_extent_begin = (istart-1)(istart-2)/2 + jstart
!!loop_extent_end = (iend-1)(iend-2)/2 + jend
!!s = 1+sqrt(1+8x)/2
!!i = int(s) - ROUNDED UP
!     qmmm_mpi%norb_istart = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_begin)))/2.0d0)
!     qmmm_mpi%norb_iend   = ceiling((1.0d0+sqrt(1.0d0+8.0d0*dble(loop_extent_end)))/2.d0)
!     qmmm_mpi%norb_loop_extent_begin = loop_extent_begin
!     qmmm_mpi%norb_loop_extent_end = loop_extent_end
!
!!Now we need to work out what range of j values we do for each i we will be doing.
!!What value of j would, when coupled with our istart give us loop_extent_begin?
!! j = loop_extent_begin -((-i-1)(i-2)/2)
!
!     jstart = loop_extent_begin - ((qmmm_mpi%norb_istart-1)*(qmmm_mpi%norb_istart-2)/2)
!     jend   = loop_extent_end - ((qmmm_mpi%norb_iend-1)*(qmmm_mpi%norb_iend-2)/2)
!
!     do i = qmmm_mpi%norb_istart, qmmm_mpi%norb_iend
!
!       if (i == qmmm_mpi%norb_istart) then
!         qmmm_mpi%norb_jrange(1,i) = jstart
!       else
!         qmmm_mpi%norb_jrange(1,i) = 1
!       end if
!
!       if (i == qmmm_mpi%norb_iend) then
!         qmmm_mpi%norb_jrange(2,i) = jend
!       else
!         qmmm_mpi%norb_jrange(2,i) = i-1
!       end if
!
!      end do
!
!#endif

      !****************************************
      !* DIAGONALIZATION ALLOCATION AND SETUP *
      !****************************************
      !****************************************
      !*** DOES NOT INCLUDE DFTB AT PRESENT ***
      !****************************************

      !Here we determine which diagonalization routine will be used and 
      !allocate the necessary scratch arrays.

      if (qmmm_nml%qmtheory /= DFTB) then
      !*** EIGEN VECTORS ***
        allocate (qm2_struct%eigen_vectors(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
        REQUIRE(ier==0)


        !*** FULL DIAGONALIZATIONS and PSEUDO DIAGONALIZATIONS ***
        if (qmmm_mpi%commqmmm_master) then
          write(6,*)
          write(6,'(''| QMMM: *** Diagonalization Routine Information ***'')')
          if (qmmm_nml%allow_pseudo_diag) &
                write(6,'(''| QMMM: Pseudo diagonalizations are allowed.'')')
        end if

        !*** PSEUDO DIAGONALIZATIONS ***

        if (qmmm_mpi%commqmmm_master .and. qmmm_nml%allow_pseudo_diag) then
          !only master needs to do the matrix diagonalizations
          allocate( qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs), &
                    qmmm_scratch%pdiag_scr_noccupied_norbs(qm2_struct%nopenclosed,qm2_struct%norbs), &
                    qmmm_scratch%pdiag_vectmp1(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
                    qmmm_scratch%pdiag_vectmp2(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
                    qmmm_scratch%pdiag_vectmp3(qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
                    qmmm_scratch%pdiag_vecjs(2,qm2_struct%nopenclosed*(qm2_struct%norbs-qm2_struct%nopenclosed)), &
                    stat=ier )
          REQUIRE(ier == 0)
        end if

        !*** FULL DIAGONALIZATIONS ***

        !  0 = Automatically pick fastest routine.
        !  1 = Use internal diagonalization routine. (default)
        !  2 = Use lapack dspev.
        !  3 = Use lapack dspevd.
        !  4 = Use lapack dspevx.
        !  5 = Use lapack dsyev.
        !  6 = Use lapack dsyevd.
        !  7 = Use lapack dsyevr.
        !  8 = Use lapack dsyevx. (not currently implemented)

        if (qmmm_nml%diag_routine == 0) then
          !automatic selection of diagonalization routine based on timings.
#ifdef NO_DETAILED_TIMINGS
          !not available with NO_DETAILED_TIMINGS - since we need a functioning wallclock call.
          call sander_bomb('qm2_load_params.f',&
                           'diag_routine=0 is not available when compiled with NO_DETAILED_TIMINGS', &
                           'Please manually select a diagonalization routine.')
#endif
          if (qmmm_mpi%commqmmm_master) then
            write(6,'(''| QMMM: Auto diagonalization routine selection is in use.'')')
            write(6,'(''| QMMM:'')')
!only master does diagonalization at present.
            call qm2_time_diag_routines(qmmm_nml%diag_routine, qmmm_nml%verbosity)
          end if
#ifdef MPI
!have a barrier here to make all threads wait for timing to have been done. We will also
!broadcast the final chosen routine. Not strictly necessary at present but may help with
!debugging later on.
          call mpi_bcast(qmmm_nml%diag_routine,1,MPI_INTEGER,0,qmmm_mpi%commqmmm,ier)
!Note we also need to broadcast the value of allow_pseudo_diag in case it changed.
          call mpi_bcast(qmmm_nml%allow_pseudo_diag, 1, mpi_logical, 0, qmmm_mpi%commqmmm, ier) 
#endif
        else
          if (qmmm_mpi%commqmmm_master) &
            write(6,'(''| QMMM: Auto diagonalization routine selection is disabled.'')')
        end if

        !Now we either have diag_routine set in qmmm_nml or was set by the auto selection routine
        !above we are now ready to do all our allocation etc.

        if (qmmm_nml%diag_routine == 1) then
          !Internal routine.
          qmmm_scratch%lapack_dc_int_scr_aloc=0
          qmmm_scratch%lapack_dc_real_scr_aloc=0
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            write(6,'(''| QMMM: Using internal diagonalization routine (diag_routine=1).'')')
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ier)
            REQUIRE(ier==0)
          end if
        else if (qmmm_nml%diag_routine == 2) then
          !DSPEV
          !dspev needs a real workspace array of size (1,3n) where n = norbs
          !We simply use the dimensions 2 to 4 of mat_diag_workspace for this
          qmmm_scratch%lapack_dc_int_scr_aloc=0
          qmmm_scratch%lapack_dc_real_scr_aloc=0
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            write(6,'(''| QMMM: Using dspev routine (diag_routine=2).'')')
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,4),stat=ier)
            REQUIRE(ier==0)
          end if
        else if (qmmm_nml%diag_routine == 3) then
          !DSPEVD
          !We will use a divide and conquer algorithm for the matrix diagonalisation.
          !now allocate the scratch arrays.
          !Divide and conquer algorithm needs a minimum of:
          !  3+5*norbs ints and 1+6*norbs+norbs**2 reals of scratch space
          !But it also supports the option of calling it to request how much space it needs.
          !We shall do this by default.
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            write(6,'(''| QMMM: Using dspevd routine (diag_routine=3).'')')
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: Calling dspevd to query required scratch array sizes.'')')
            end if
            call dspevd('V','U',qm2_struct%norbs,qm2_struct%fock_matrix, &
            qmmm_scratch%mat_diag_workspace(1,1), &
            qm2_struct%eigen_vectors, qm2_struct%norbs, &
            aloc_real_scr, -1, &
            qmmm_scratch%lapack_dc_int_scr_aloc , -1, &
            ier)
            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: dspevd required REAL scratch = '',i12,'' REALS.'')') &
                                                 qmmm_scratch%lapack_dc_real_scr_aloc
              write(6,'(''QMMM: dspevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                                 qmmm_scratch%lapack_dc_int_scr_aloc
            end if
!             qmmm_scratch%lapack_dc_real_scr_aloc=1+6*qm2_struct%norbs+(qm2_struct%norbs*qm2_struct%norbs)
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
!            qmmm_scratch%lapack_dc_int_scr_aloc=3+5*qm2_struct%norbs
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
            REQUIRE(ier==0)
          else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if
        else if (qmmm_nml%diag_routine == 4) then
          !DSPEVX
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            write(6,'(''| QMMM: Using dspevx routine (diag_routine=4).'')')

            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            qmmm_scratch%lapack_dc_real_scr_aloc = qm2_struct%norbs*(qm2_struct%norbs+1)/2
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            !Note for dspevx we pad lapack_dc_int_scr by an extra norbs so we can use the 
            !first norbs for IFAIL.
            qmmm_scratch%lapack_dc_int_scr_aloc=6*qm2_struct%norbs
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
            REQUIRE(ier==0)
          else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if
        else if (qmmm_nml%diag_routine == 5) then
          !DSYEV  - unpacked diagonalizer
          !DSYEV has complex requirements of scratch arrays. 
          !But it also supports the option of calling it to request how much space it needs.
          !We shall do this by default.
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            write(6,'(''| QMMM: Using dsyev routine (diag_routine=5).'')')
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: Calling dsyev to query required scratch array sizes.'')')
            end if
            call dsyev('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, &
                       qm2_struct%norbs, qmmm_scratch%mat_diag_workspace(1,1), &
                       aloc_real_scr, -1, ier)
            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: dsyev required REAL scratch = '',i12,'' REALS.'')') &
                                                 qmmm_scratch%lapack_dc_real_scr_aloc
            end if
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
          
          else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if
        else if (qmmm_nml%diag_routine == 6) then
          !DSYEVD - unpacked divide and conquor diagonalizer
          !DSYEVD has complex requirements of scratch arrays. 
          !But it also supports the option of calling it to request how much space it needs.
          !We shall do this by default.
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            write(6,'(''| QMMM: Using dsyevd routine (diag_routine=6).'')')
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: Calling dsyevd to query required scratch array sizes.'')')
            end if
            call dsyevd('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, qm2_struct%norbs, &
                        qmmm_scratch%mat_diag_workspace(1,1), aloc_real_scr, -1, &
                        qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)
            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: dsyevd required REAL scratch = '',i12,'' REALS.'')') &
                                                 qmmm_scratch%lapack_dc_real_scr_aloc
              write(6,'(''QMMM: dsyevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                                 qmmm_scratch%lapack_dc_int_scr_aloc
            end if     
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
            REQUIRE(ier==0)
          else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if
        else if (qmmm_nml%diag_routine == 7) then
          !DSYEVR - Reduction to tridiagonal form before diagonalization.
          !DSYEVR has complex requirements of scratch arrays. 
          !NOTE: for dsyevr we can't just put the matrix into the eigenvectors array as
          !      it uses the matrix for scratch space. We will assume we can use the
          !      pdiag_scr_norbs_norbs scratch - array which will get allocated if we
          !      are doing pseudo diags and/or we are choosing diag_routine == 7.
          !NOTE2: The lapack_dc_int array is actually made 2xnorbs bigger than what is
          !       requested since the first 2xnorbs are used for ISUPPZ.
          !But it also supports the option of calling it to request how much space it needs.
          !We shall do this by default.
          if (qmmm_mpi%commqmmm_master) then !Only master does diagonalisation
            allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
            REQUIRE(ier==0)
            write(6,'(''| QMMM: Using dsyevr routine (diag_routine=7).'')')
!            if (qmmm_nml%verbosity >= 2) then
!              write(6,'(''QMMM: Calling dsyevr to query required scratch array sizes.'')')
!            end if

! RCW 2008/01/22 - weird segfaults in qm2_hcore_qmqm when we call dsyevr to calculate the memory requirements
!                  so for the moment we will just use the default amounts.
            !The dimension of the array WORK.  LWORK >= max(1,26*N).
            !For optimal efficiency, LWORK >= (NB+6)*N,
            !where NB is the max of the blocksize for DSYTRD and DORMTR
            !returned by ILAENV.
            ilaenv_blocksize = ILAENV( 1, 'DSYTRD', 'U', qm2_struct%norbs, -1, -1, -1 )
            ilaenv_blocksize = max(ilaenv_blocksize, ILAENV( 1, 'DORMTR', 'U', qm2_struct%norbs, -1, -1, -1 ))

            qmmm_scratch%lapack_dc_real_scr_aloc=(ilaenv_blocksize+6)*qm2_struct%norbs

            !The dimension of the array IWORK.  LIWORK >= max(1,10*N).
            qmmm_scratch%lapack_dc_int_scr_aloc = 10*qm2_struct%norbs

!            abstol = 2.0d0 * dlamch('S')
!            call dsyevr('V','A','U',qm2_struct%norbs,qmmm_scratch%pdiag_scr_norbs_norbs,qm2_struct%norbs, &
!                 0.0d0, 0.0d0, 0, 0, abstol, int_dummy, qmmm_scratch%mat_diag_workspace(1,1), &
!                 qm2_struct%eigen_vectors, qm2_struct%norbs, qmmm_scratch%lapack_dc_int_scr_aloc, &
!                 aloc_real_scr, -1, qmmm_scratch%lapack_dc_int_scr_aloc, -1, ier)
!            qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)

            !for dsyevr we pad lapack_dc_int_scr with 2*norbs to allow use in ISUPPZ.
            qmmm_scratch%lapack_dc_int_scr_aloc = qmmm_scratch%lapack_dc_int_scr_aloc+2*qm2_struct%norbs
            if (qmmm_nml%verbosity >= 2) then
              write(6,'(''QMMM: dsyevr required REAL scratch = '',i12,'' REALS.'')') &
                                                 qmmm_scratch%lapack_dc_real_scr_aloc
              write(6,'(''QMMM: dsyevr required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                                 qmmm_scratch%lapack_dc_int_scr_aloc
            end if     
            allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
            REQUIRE(ier==0)
            if (.not. qmmm_nml%allow_pseudo_diag) then
              !qmmm_scratch%pdiag_scr_norbs_norbs will not have been allocated so make sure
              !we allocate it.
              allocate(qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
              REQUIRE(ier==0)
            end if
          else
            qmmm_scratch%lapack_dc_int_scr_aloc=0
            qmmm_scratch%lapack_dc_real_scr_aloc=0
          end if

        else
          call sander_bomb('qm2_load_params','METHOD NOT CURRENTLY IMPLEMENTED', &
                           'Selected Diagonalization routine is not currently implemented.')
        end if !diag_routine

      end if !not DFTB

      !********************************************
      !* END DIAGONALIZATION ALLOCATION AND SETUP *
      !********************************************

      return
end subroutine qm2_load_params_and_allocate

subroutine qm2_time_diag_routines(diag_routine,verbosity)

  !This routine tries a series of different diagonalizers and returns the value
  !of diag_method that is fastest as the first argument.

  use qmmm_module, only : qmmm_scratch, qm2_struct, qmmm_nml
  use constants, only : zero

  implicit none

  integer, intent(out) :: diag_routine
  integer, intent(in) :: verbosity

!locals
  _REAL_ :: abstol, smallsum, small
  _REAL_ :: dlamch  !is a function
  _REAL_ :: start_time, end_time, current_time, fastest_time
  _REAL_ :: aloc_real_scr
  integer :: ilaenv !is a function
  integer :: ilaenv_blocksize
  integer :: diag_iterations
  integer :: i, ier
  integer :: diag_routine_to_test

  abstol = 2.0d0 * dlamch('S') !tolerance for dspevr

  !We need to do enough calls to generate decent statistics. The Fortran timers do not appear to
  !be incredibly good so we will call the diagonalization a number of times based on the size of
  !the matrix.

  if (qm2_struct%norbs <= 50) then
    diag_iterations = 1000
  elseif (qm2_struct%norbs <= 100 ) then
    diag_iterations = 250
  elseif (qm2_struct%norbs <= 150 ) then
    diag_iterations = 100
  elseif (qm2_struct%norbs <= 200 ) then
    diag_iterations = 50
  elseif (qm2_struct%norbs <= 250 ) then
    diag_iterations = 40
  elseif (qm2_struct%norbs <= 300 ) then
    diag_iterations = 25
  elseif (qm2_struct%norbs <= 350 ) then
    diag_iterations = 15
  elseif (qm2_struct%norbs <= 400 ) then
    diag_iterations = 10
  elseif (qm2_struct%norbs <= 450 ) then
    diag_iterations = 5
  elseif (qm2_struct%norbs <= 500 ) then
    diag_iterations = 4
  elseif (qm2_struct%norbs <= 550 ) then
    diag_iterations = 3
  elseif (qm2_struct%norbs <= 600 ) then
    diag_iterations = 2
  else
    diag_iterations = 1
  endif

  write(6,'(''| QMMM: Timing diagonalization routines:'')')
  write(6,'(''| QMMM:                              norbs = '',i8)') qm2_struct%norbs
  write(6,'(''| QMMM:    diag iterations used for timing = '',i8)') diag_iterations
  write(6,'(''| QMMM:'')')

  !--------------------
  !1) INTERNAL ROUTINE
  !--------------------

  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ier)
  REQUIRE(ier==0)

  diag_routine_to_test = 1
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:              Internal diag routine = '',F8.2,'' seconds'')') current_time

  !initially just assume that the internal routine is the fastest.
  fastest_time = current_time
  diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.

  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)

  !--------------------
  !END INTERNAL ROUTINE
  !--------------------

  !--------------------
  !2) DSPEV ROUTINE
  !--------------------
  !DSPEV
  !dspev needs a real workspace array of size (1,3n) where n = norbs
  !We simply use the dimensions 2 to 4 of mat_diag_workspace for this
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,4),stat=ier)
  REQUIRE(ier==0)
  
  diag_routine_to_test = 2
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                 Dspev diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  
  deallocate (qmmm_scratch%mat_diag_workspace, stat=ier)
  REQUIRE(ier==0)
  !--------------------
  !END DSPEV ROUTINE
  !--------------------

  !--------------------
  !3) DSPEVD ROUTINE
  !--------------------
  if (verbosity >= 2) then
     write(6,'(''QMMM: Calling dspevd to query required scratch array sizes.'')')
  end if
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
  REQUIRE(ier==0)
  call dspevd('V','U',qm2_struct%norbs,qm2_struct%fock_matrix, &
  qmmm_scratch%mat_diag_workspace(1,1), &
  qm2_struct%eigen_vectors, qm2_struct%norbs, aloc_real_scr, -1, &
  qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)

  qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)

  if (verbosity >= 2) then
    write(6,'(''QMMM: dspevd required REAL scratch = '',i12,'' REALS.'')') &
                                       qmmm_scratch%lapack_dc_real_scr_aloc
    write(6,'(''QMMM: dspevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                       qmmm_scratch%lapack_dc_int_scr_aloc
  end if
  allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
  REQUIRE(ier==0)
  allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
  REQUIRE(ier==0)

  diag_routine_to_test = 3
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                Dspevd diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_real_scr_aloc = 0
  qmmm_scratch%lapack_dc_int_scr_aloc = 0

  !--------------------
  !END DSPEVD ROUTINE
  !--------------------

  !--------------------
  !4) DSPEVX ROUTINE
  !--------------------
  !Allocate mat diag workspace as 2 here since (x,2) is used for error reporting by dspevx
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,2),stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_real_scr_aloc = qm2_struct%norbs*(qm2_struct%norbs+1)/2
  allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_int_scr_aloc=5*qm2_struct%norbs
  allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc), stat=ier)
  REQUIRE(ier==0)

  diag_routine_to_test = 4
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                Dspevx diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_real_scr_aloc = 0
  qmmm_scratch%lapack_dc_int_scr_aloc = 0

  !--------------------
  !END DSPEVX ROUTINE
  !--------------------

  !--------------------
  !5) DSYEV ROUTINE
  !--------------------
  if (verbosity >= 2) then
     write(6,'(''QMMM: Calling dsyev to query required scratch array sizes.'')')
  end if
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
  REQUIRE(ier==0)
  call dsyev('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, &
             qm2_struct%norbs, qmmm_scratch%mat_diag_workspace(1,1), &
             aloc_real_scr, -1, ier)
  qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
  if (verbosity >= 2) then
      write(6,'(''QMMM: dsyev required REAL scratch = '',i12,'' REALS.'')') &
                                        qmmm_scratch%lapack_dc_real_scr_aloc
  end if
  allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
  REQUIRE(ier==0)

  diag_routine_to_test = 5
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                 Dsyev diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_real_scr_aloc = 0
  qmmm_scratch%lapack_dc_int_scr_aloc = 0

  !--------------------
  !END DSYEV ROUTINE
  !--------------------

  !--------------------
  !6) DSYEVD ROUTINE
  !--------------------
  if (verbosity >= 2) then
     write(6,'(''QMMM: Calling dsyevd to query required scratch array sizes.'')')
  end if
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
  REQUIRE(ier==0)
  call dsyevd('V','U',qm2_struct%norbs,qm2_struct%eigen_vectors, qm2_struct%norbs, &
  qmmm_scratch%mat_diag_workspace(1,1), aloc_real_scr, -1, &
  qmmm_scratch%lapack_dc_int_scr_aloc , -1, ier)
  qmmm_scratch%lapack_dc_real_scr_aloc = int(aloc_real_scr)
  if (verbosity >= 2) then
    write(6,'(''QMMM: dsyevd required REAL scratch = '',i12,'' REALS.'')') &
                                       qmmm_scratch%lapack_dc_real_scr_aloc
    write(6,'(''QMMM: dsyevd required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                       qmmm_scratch%lapack_dc_int_scr_aloc
  end if
  allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
  REQUIRE(ier==0)
  allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
  REQUIRE(ier==0)

  diag_routine_to_test = 6
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                Dsyevd diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)
  qmmm_scratch%lapack_dc_real_scr_aloc = 0
  qmmm_scratch%lapack_dc_int_scr_aloc = 0

  !--------------------
  !END DSYEVD ROUTINE
  !--------------------

  !--------------------
  !7) DSYEVR ROUTINE
  !--------------------
  allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
  REQUIRE(ier==0)
! RCW 2008/01/22 - weird segfaults in qm2_hcore_qmqm when we call dsyevr to calculate the memory requirements
!                  so for the moment we will just use the default amounts.
  !The dimension of the array WORK.  LWORK >= max(1,26*N).
  !For optimal efficiency, LWORK >= (NB+6)*N,
  !where NB is the max of the blocksize for DSYTRD and DORMTR
  !returned by ILAENV.
  ilaenv_blocksize = ILAENV( 1, 'DSYTRD', 'U', qm2_struct%norbs, -1, -1, -1 )
  ilaenv_blocksize = max(ilaenv_blocksize, ILAENV( 1, 'DORMTR', 'U', qm2_struct%norbs, -1, -1, -1 ))
  qmmm_scratch%lapack_dc_real_scr_aloc=(ilaenv_blocksize+6)*qm2_struct%norbs
  !The dimension of the array IWORK.  LIWORK >= max(1,10*N).
  qmmm_scratch%lapack_dc_int_scr_aloc = 10*qm2_struct%norbs
  !for dsyevr we pad lapack_dc_int_scr with 2*norbs to allow use in ISUPPZ.
  qmmm_scratch%lapack_dc_int_scr_aloc = qmmm_scratch%lapack_dc_int_scr_aloc+2*qm2_struct%norbs
  if (verbosity >= 2) then
     write(6,'(''QMMM: dsyevr required REAL scratch = '',i12,'' REALS.'')') &
                                      qmmm_scratch%lapack_dc_real_scr_aloc
     write(6,'(''QMMM: dsyevr required INTEGER scratch = '',i12,'' INTEGERS.'')') &
                                      qmmm_scratch%lapack_dc_int_scr_aloc
  end if
  allocate (qmmm_scratch%lapack_dc_real_scr(qmmm_scratch%lapack_dc_real_scr_aloc),stat=ier)
  REQUIRE(ier==0)
  allocate (qmmm_scratch%lapack_dc_int_scr(qmmm_scratch%lapack_dc_int_scr_aloc),stat=ier)
  REQUIRE(ier==0)
  if (.not. qmmm_nml%allow_pseudo_diag) then
     !qmmm_scratch%pdiag_scr_norbs_norbs will not have been allocated so make sure
     !we allocate it.
     allocate(qmmm_scratch%pdiag_scr_norbs_norbs(qm2_struct%norbs,qm2_struct%norbs),stat=ier)
     REQUIRE(ier==0)
  end if

  diag_routine_to_test = 7
  current_time=0.0d0
  do i = 1,diag_iterations
    !Need to rebuild the fock matrix each time since some diag routines can destroy it.
    call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
    call wallclock(start_time)
    call qm2_full_diagonalize(diag_routine_to_test,qm2_struct%fock_matrix,qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
    call wallclock(end_time)
    current_time = current_time + (end_time-start_time)
  end do

  write(6,'(''| QMMM:                Dsyevr diag routine = '',F8.2,'' seconds'')') current_time

  !check if this routine is faster than the ones tried so far.
  if (fastest_time >= current_time) then
    diag_routine = diag_routine_to_test !Note reason for this here is if the routine actually failed to converge then
                                        !diag_routine_to_test will have been set back to the internal diagonalizer
                                        !by the diag routine.
    fastest_time = current_time
  end if
  deallocate (qmmm_scratch%lapack_dc_int_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%lapack_dc_real_scr,stat=ier)
  REQUIRE(ier==0)
  deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
  REQUIRE(ier==0)
  if (.not. qmmm_nml%allow_pseudo_diag) then
     deallocate(qmmm_scratch%pdiag_scr_norbs_norbs,stat=ier)
     REQUIRE(ier==0)
  end if

  qmmm_scratch%lapack_dc_real_scr_aloc = 0
  qmmm_scratch%lapack_dc_int_scr_aloc = 0
 
  !--------------------
  !END DSYEVR ROUTINE
  !--------------------

  write(6,'(''| QMMM:'')')

  !-----------------------
  !PSEUDO DIAGONALIZATIONS
  !-----------------------
  !Finally we will test how long pseudo diagonalizations
  !take - this is not an exact science at all since the ratio
  !of pseudo diags to real diags does not remain constant. 
  !We will simply check how it compares like for like with
  !full diagonalizations. If the ful diagonalization is faster
  !than the pseudo diagonalization then we will turn off the
  !pseudo diagonalization. This will give faster execution but
  !may miss cases where it is still beneficial to turn off
  !the pseudo diag.
  if (qmmm_nml%allow_pseudo_diag) then
     allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,1),stat=ier)
     REQUIRE(ier==0)

      call qm2_smallest_number(small,smallsum)
      smallsum = max(10.0D0 * sqrt(smallsum),1.4000D-7)

      current_time=0.0d0
      do i = 1,diag_iterations
        !Need to rebuild the fock matrix each time since some diag routines can destroy it.
        call qm2_time_diag_routines_random(qm2_struct%matsize,qm2_struct%fock_matrix)
        call wallclock(start_time)
        call qm2_pseudo_diag(qm2_struct%fock_matrix,qm2_struct%eigen_vectors,qm2_struct%nopenclosed, &
                   qmmm_scratch%mat_diag_workspace(1,1),qm2_struct%norbs,smallsum, &
                   qmmm_scratch%pdiag_scr_norbs_norbs,qmmm_scratch%pdiag_scr_noccupied_norbs, &
                   qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2,qmmm_scratch%pdiag_vectmp3, &
                   qmmm_scratch%pdiag_vecjs)
        call wallclock(end_time)
        current_time = current_time + (end_time-start_time)
      end do

     write(6,'(''| QMMM:                Pseudo diag routine = '',F8.2,'' seconds'')') current_time

    !check if this routine is slower than the ones tried so far.
    if (fastest_time < current_time) then
      !The pseudo diag routine is slower than regular diags.
      write(6,'(''| QMMM: Pseudo diagonalization appears to be slower than regular'')')
      write(6,'(''| QMMM: diagonalization. Setting pseudo_diag=0 for optimum performance.'')')
      qmmm_nml%allow_pseudo_diag = .false.

      deallocate( qmmm_scratch%pdiag_scr_norbs_norbs, &
                  qmmm_scratch%pdiag_scr_noccupied_norbs, &
                  qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2, &
                  qmmm_scratch%pdiag_vectmp3, qmmm_scratch%pdiag_vecjs, &
                  stat=ier )
      REQUIRE(ier == 0)
    end if
    deallocate (qmmm_scratch%mat_diag_workspace,stat=ier)
    REQUIRE(ier==0)
  end if
  !---------------------------
  !END PSEUDO DIAGONALIZATIONS
  !---------------------------
  write(6,'(''| QMMM:'')')

  return

end subroutine qm2_time_diag_routines

subroutine qm2_time_diag_routines_random(matsize, fock_matrix)

  implicit none

  integer, intent(in) :: matsize
  _REAL_, intent(out) :: fock_matrix(matsize)

  integer :: i
  _REAL_ :: random_nmbr

  !Build a fock matrix to diagonalize
  !Since we don't have enough to build a real fock matrix here what we will
  !do is just fill the matrix elements with a random distribution.
  !the matrix is normally dense with elements ranging from -100 to +100
  !based on the 1NLN test case.
  do i = 1, matsize
    call amrand(random_nmbr)
    random_nmbr = 100.0d0 - ( 200.0d0 * random_nmbr )
    fock_matrix(i) = random_nmbr
  end do

  return

end subroutine qm2_time_diag_routines_random

