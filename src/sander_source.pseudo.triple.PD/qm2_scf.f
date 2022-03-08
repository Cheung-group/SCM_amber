! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "def_time.h"
#include "assert.h"
subroutine qm2_scf(fock_matrix, H, W, escf, den_matrix, scf_mchg)
!---------------------------------------------------------------------
! This is the main SCF routine.
! Written by Ross Walker (TSRI, 2005)
!
! This routine generates a self consistent field
! and returns the energy in the variable ESCF.
!
! The main arrays used are:
!     fock_matrix - Starts off containing the one-electron matrix
!                   and is used to hold the Fock matrix. Units = Ev/Bohr
!     H - Only ever contains the one electron matrix.
!     W - Only ever contains the two electron matrix.
!     P - Only ever contains the total density matrix of the 
!         current SCF step.
!
! qm2_struct%matsize = size of packed triangle (NORBS*(NORBS+1)/2)
! qm2_struct%n2el = Number of 2 electron integrals
!                    50*nheavy(nheavy-1)+10*nheavy*nlight+(nlight*(nlight-1))/2
!     
! References:
! PSEUDO DIAGONALISATION : "FAST SEMIEMPIRICAL CALCULATIONS",
!                           STEWART. J.J.P., CSASZAR, P., PULAY, P.,
!                           J. COMP. CHEM.,3, 227, (1982)
!---------------------------------------------------------------------

    use qmmm_module, only : qmmm_struct, qm2_struct, qm2_params, qmewald, qmmm_nml, &
                            qm_gb, qmmm_mpi, qmmm_scratch
    use constants, only : EV_TO_KCAL, zero, two
    implicit none


#ifdef MPI
#include "mpif.h"
integer :: ier
#endif

!Passed in
    _REAL_, intent(inout) :: fock_matrix(qm2_struct%matsize)
    _REAL_, intent(in) :: H(qm2_struct%matsize)
    _REAL_, intent(in) :: W(qm2_struct%n2el)
    _REAL_, intent(inout) :: den_matrix(qm2_struct%matsize)
    _REAL_, intent(inout) :: escf
    _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

!Local
    _REAL_ eold !SCF energy on previous step
    _REAL_ energy_diff, density_diff !Difference in energy and density from previous step.a
    _REAL_ small, smallsum, abstol !Precision limits of this machine.
    _REAL_ smallest_energy_diff(2) !Smallest energy diff found so far (1) and density diff for this step(2)
    _REAL_ scf_energy !Computed in parts on different cpus and then all reduced. Only master returns
                      !this value in escf.
    integer sm_energy_diff_step_number
    integer scf_iteration !Number of scf iteration that have been done.
    integer i
    logical converged, first_iteration
    logical doing_pseudo_diag !Initially false, set to true if we are doing a pseudo diagonalisation.
    logical allow_pseudo_diag !Initially set to false. Set to true after at least 2 SCF iterations have been done 
                              !and allowed by the namelist option. Set back to full on last step to force a full
                              !diagonalisation.
    logical pseudo_converged !Flag used to indicate when doing pseudo diagonalisations has converged the SCF.
                                  !Once this is set to try the value of allow_pseudo_diag will not be set to true
                                  !again. This is used to force the last few SCF steps to be full diagonalisations.

!qm2_Helect is a function
    _REAL_ qm2_HELECT

!dlamch is a function
    _REAL_ dlamch

#ifdef MPI
    _REAL_ tmp_recv(2)
# ifndef USE_MPI_IN_PLACE
    _REAL_ tmp_send(2)
# endif
#endif

    integer lapack_info

!Saves
    save smallsum
    save abstol !Underflow limit for dspevr

!Initialisation on first call
    if (qmmm_struct%qm2_scf_first_call) then
      !Find precision limits of this machine
      call qm2_smallest_number(small,smallsum)
      ! smallsum should be the smallest number for which 1.0D0 + smallsum /= 1.0D0
      ! or 1.0D-17 - whichever is larger.
      ! We will increase it in order to allow for roundoff in our calculations.
      ! we use max here to avoid problems which occur when smallsum is actually too small. 
      smallsum = max(10.0D0 * sqrt(smallsum),1.4000D-7)
!      smallsum = 10.0D0 * sqrt(smallsum)
      abstol = 2.0d0 * dlamch('S') !tolerance for dspevr
      qmmm_struct%qm2_scf_first_call=.false.
    end if
!Initialisation on every call
    converged = .false.
    first_iteration = .true.
    doing_pseudo_diag = .false.
    allow_pseudo_diag = .false.
    pseudo_converged = .false.
    eold = zero
    energy_diff = huge(energy_diff)
    density_diff = huge(density_diff)
    smallest_energy_diff(1) = huge(smallest_energy_diff(1))
    smallest_energy_diff(2) = huge(smallest_energy_diff(2))
    sm_energy_diff_step_number = 0

    if (qmmm_nml%verbosity > 2 .and. qmmm_mpi%commqmmm_master) then
      write(6,'("QMMM: ")')
      write(6,'("QMMM: SCF Convergence Information")')
      if (qmmm_nml%allow_pseudo_diag) write(6,'("QMMM: (*) = Pseudo Diagonalisation")')
      write(6,'("QMMM: Cycle         Energy       ")',ADVANCE='NO')
      write(6,'("           dE                    dP")')
    end if

! MAIN SCF LOOP
    do_scf: do scf_iteration=1, qmmm_nml%itrmax
      if (.NOT. first_iteration) then

        ! Step 1 - We haven't converged so we need to get a better Fock matrix.
        !         Diagonalise the RHF secular determinant
        ! Diagonalise the RHF secular determinant
        ! We have two options here. If we are allowed, and the density matrix 
        ! fluctuations are small enough we should do a pseudo diagonalisation 
        ! instead of a full one.
#ifdef MPI
        !Reduce the fock matrix to master thread to do diagonalisation
        call timer_start(TIME_QMMMENERGYSCFFOCKRED)
# ifdef USE_MPI_IN_PLACE
        if (qmmm_mpi%commqmmm_master) then
          call mpi_reduce(MPI_IN_PLACE,qm2_struct%fock_matrix,qm2_struct%matsize, &
                        MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        else
          call mpi_reduce(qm2_struct%fock_matrix,0,qm2_struct%matsize, &
                        MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        end if
# else
        call mpi_reduce(qm2_struct%fock_matrix,qmmm_scratch%matsize_red_scratch,qm2_struct%matsize, &
                        MPI_DOUBLE_PRECISION,mpi_sum,0,qmmm_mpi%commqmmm,ier)
        if (qmmm_mpi%commqmmm_master) &
           qm2_struct%fock_matrix(1:qm2_struct%matsize)=qmmm_scratch%matsize_red_scratch(1:qm2_struct%matsize)
# endif
        call timer_stop(TIME_QMMMENERGYSCFFOCKRED)
#endif

        if (allow_pseudo_diag .AND. density_diff <= qmmm_nml%pseudo_diag_criteria) then

          !We can do a pseudo diagonalisation.
          doing_pseudo_diag = .true. !Marker used to tell SCF routine to do a full diagonalisation before quitting.

          !Dimension 1 of mat_diag_workspace contains the eignvalues
          !Dimension 2-6 is used as scratch space.

#ifdef MPI
          !only master does the diagonalisation for time being.
          if (qmmm_mpi%commqmmm_master) then
#endif
          call timer_start(TIME_QMMMENERGYSCFPSEUDO)
!OPENMP PARALLEL
          call qm2_pseudo_diag(fock_matrix,qm2_struct%eigen_vectors,qm2_struct%nopenclosed, &
               qmmm_scratch%mat_diag_workspace(1,1),qm2_struct%norbs,smallsum, &
               qmmm_scratch%pdiag_scr_norbs_norbs,qmmm_scratch%pdiag_scr_noccupied_norbs, &
               qmmm_scratch%pdiag_vectmp1,qmmm_scratch%pdiag_vectmp2,qmmm_scratch%pdiag_vectmp3, &
               qmmm_scratch%pdiag_vecjs)
          call timer_stop(TIME_QMMMENERGYSCFPSEUDO)
#ifdef MPI
          end if
#endif
        else
#ifdef MPI
          !only master does the diagonalisation for time being.
          if (qmmm_mpi%commqmmm_master) then
#endif
          !Do a full diagonalisation
          call timer_start(TIME_QMMMENERGYSCFDIAG)
          call qm2_full_diagonalize(qmmm_nml%diag_routine,fock_matrix, &
                                    qm2_struct%norbs,qm2_struct%eigen_vectors,abstol)
          call timer_stop(TIME_QMMMENERGYSCFDIAG)
#ifdef MPI
          end if
#endif
        end if !pseudo diag.
        !End of step 1

        !Step 2 - Calculate the density matrix:
#ifdef MPI
        !only master calculates density matrix for time being.
        if (qmmm_mpi%commqmmm_master) then
#endif
        call timer_start(TIME_QMMMENERGYSCFDEN)
!OPENMP PARALLEL
        call qm2_densit(qm2_struct%eigen_vectors,qm2_struct%norbs,qm2_struct%nclosed,den_matrix,qm2_struct%matsize)
!OPENMP PARALLEL
        call qm2_cnvg(den_matrix, qm2_struct%old_den_matrix, qm2_struct%old2_density, &
             qm2_struct%norbs, scf_iteration, density_diff)
        call timer_stop(TIME_QMMMENERGYSCFDEN)
#ifdef MPI
        else
          density_diff=zero
        end if
        call timer_start(TIME_QMMMENERGYSCFDENBCAST)
        !only master calculated the density matrix
        call mpi_bcast(qm2_struct%den_matrix, qm2_struct%matsize, MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
        call timer_stop(TIME_QMMMENERGYSCFDENBCAST)
#endif
        !End of step 2
      end if ! if (.NOT. first_iteration)

      !Calculate the Mulliken charges for the current density matrix if we require
      !them on every SCF step. Save the results in the scf_mchg array
      !We also need to do this if we are doing QMewald with the image charges fixed but
      !it is our first MD step.
      !For the moment let all threads calculate mulliken charges as this is
      !probably quicker than doing a reduce.
      if (qm2_struct%calc_mchg_scf .or. qmewald%ewald_startup) then
        do i=1,qmmm_struct%nquant_nlink 
          call qm2_calc_mulliken(i,scf_mchg(i))
        end do
      end if

      !Step 3 - Build the fock matrix !Step 1 if this is our first SCF iteration.
      call timer_start(TIME_QMMMENERGYSCFFOCK)

      !Copy the one electron matrix into the FOCK matrix
      fock_matrix(1:qm2_struct%matsize)=H(1:qm2_struct%matsize)
!Parallel
      call qm2_fock2(fock_matrix,den_matrix,W,qm2_params%orb_loc)
!Parallel
      call qm2_fock1(fock_matrix,den_matrix)

      if ( qmmm_nml%qm_ewald>0 ) then
        call timer_start(TIME_QMMMENERGYSCFFOCKEWALD)
        !Add the qm_ewald contributions to the diagonal elements of
        !the fock matrix. This requires the Mulliken Charges
        !Note we skip this if we are keeping the image charges fixed during the SCF but
        !only after the first MD step has been done.
        if (qmmm_nml%qm_ewald==1 .OR. qmewald%ewald_startup) then
!Parallel
          call qm_ewald_qm_pot(qmmm_struct%nquant, qmmm_struct%nlink, scf_mchg, &
                                qmmm_struct%qm_coords,qmewald%kvec)
        end if
        call qm_ewald_add_fock(fock_matrix, qmewald%qmpot, qmewald%mmpot)
        call timer_stop(TIME_QMMMENERGYSCFFOCKEWALD)
      end if

      call timer_stop(TIME_QMMMENERGYSCFFOCK)
      !End step 3

      !Step 4 - Calculate the energy in KCal/mol
      call timer_start(TIME_QMMMENERGYSCFELEC)
      qmmm_struct%elec_eng = qm2_HELECT(qm2_struct%NORBS-1,den_matrix,H,fock_matrix)
      if (qmmm_nml%qm_ewald>0) call  qm_ewald_correct_ee(qmmm_struct%elec_eng, &
                                                       qmewald%mmpot, den_matrix)
#ifdef MPI
      !Open MP Diagonalisation - only master has density diff - reduce it along with energy
      !Reduce the elec_eng and density diff to all threads so they can check for convergence
# ifdef USE_MPI_IN_PLACE
      tmp_recv(1) = qmmm_struct%elec_eng
      tmp_recv(2) = density_diff
      call mpi_allreduce(MPI_IN_PLACE,tmp_recv,2, &
                        MPI_DOUBLE_PRECISION,mpi_sum,qmmm_mpi%commqmmm,ier)
      qmmm_struct%elec_eng=tmp_recv(1)
      density_diff=tmp_recv(2)
# else
      tmp_send(1) = qmmm_struct%elec_eng
      tmp_send(2) = density_diff
      call mpi_allreduce(tmp_send,tmp_recv,2, &
                        MPI_DOUBLE_PRECISION,mpi_sum,qmmm_mpi%commqmmm,ier)
      qmmm_struct%elec_eng=tmp_recv(1)
      density_diff=tmp_recv(2)
# endif
#endif
      scf_energy = qmmm_struct%elec_eng*EV_TO_KCAL
      energy_diff = scf_energy - eold
      if (qmmm_mpi%commqmmm_master) then
        if (abs(energy_diff) < abs(smallest_energy_diff(1))) then
          !Keep track of the smallest difference we have found. Useful for
          !when we fail to converge - we can tell the user how close we ever got.
          !Ignore changes that are zero.
          if (energy_diff/=zero)  smallest_energy_diff(1) = energy_diff
          if (density_diff/=zero) smallest_energy_diff(2) = density_diff
          if (density_diff/=zero .and. energy_diff/=zero) sm_energy_diff_step_number = scf_iteration
        end if
        !If verbosity is >2 then print some info about this SCF step.
        if (qmmm_nml%verbosity > 2 .and. qmmm_mpi%commqmmm_master) then
           if (doing_pseudo_diag) then
             write(6,'("QMMM: ",i5,"*",3G22.14)') scf_iteration, scf_energy, energy_diff, density_diff
           else
             write(6,'("QMMM: ",i5," ",3G22.14)') scf_iteration, scf_energy, energy_diff, density_diff
           end if
           if (qmmm_nml%verbosity > 4) then
              !also print info in KJ/mol
              write(6,'("QMMM: KJ/mol  ",3G22.14)') scf_energy*4.184d0,energy_diff*4.184d0,density_diff*4.184d0
           end if
        end if !(qmmm_nml%verbosity > 2)
        !End step 4
      end if
      call timer_stop(TIME_QMMMENERGYSCFELEC)
      ! Step 5 - Check if we have converged.
      !    check energy difference is less than qmmm_nml%scfconv
      !    check density difference is either zero or less than density_conv
      !    Make sure we have done at least 2 SCF iterations to avoid quiting
      !         due to some fluke.
      if (scf_iteration>2) then

        ! Since we have now done more than 2 iterations we can allow the
        ! pseudo diagonalisation as long as it is allowed from the 
        ! namelist pseudo_diag option.
        ! Also, if doing pseudo diagonalisations has allowed us to converge 
        ! and we are just doing a few more scf steps with full pseudo 
        ! diagonalisation then make sure we don't turn pseudo diagonalisation 
        ! back on by accident.

        if (qmmm_nml%allow_pseudo_diag .AND. (.not. pseudo_converged)) allow_pseudo_diag = .true.
        if (abs(energy_diff) < qmmm_nml%scfconv .AND. &
            (density_diff < qmmm_nml%density_conv .OR. density_diff == zero)) then

           ! We have converged. However, if the last SCF step was a pseudo
           ! diagonalisation we need to do a few more loops with full 
           ! diagonalisation.  Otherwise the forces will not be accurate enough.

           if (doing_pseudo_diag) then
             doing_pseudo_diag=.false.
             allow_pseudo_diag=.false.
             pseudo_converged=.true. !Stops any more pseudo diags being done for the rest of the SCF loop.
           else
             converged = .true.
             exit do_scf  !Break out of the SCF loop since we have converged.
           end if
        end if
      end if !scf_iteration>2
      !End of step 5
      eold = scf_energy !Copy the energy that this step got us for use next time.

      !Step 6 - add any extra terms to the SCF for which we only need density convergence.
      !Add in GB terms if qmgb=2
      call timer_start(TIME_QMMMENERGYSCFFOCK)
      if ( qmmm_nml%qmgb == 2 ) then
        call timer_start(TIME_QMMMENERGYSCFFOCKGB)
        !Step 1 - calculate the potential at QM atoms due to QM atoms with current Mulliken charges
        !Parallel
        qm_gb%gb_qmpot(1:qmmm_struct%nquant_nlink)=zero
        call qmgb_calc_qm_pot(qm_gb%gb_qmpot,qm_gb%qmqm_onefij,scf_mchg)
        !Step 2 - Add the mm potential and the qm potential to the fock matrix.
        call qmgb_add_fock(qmmm_struct%nquant_nlink, fock_matrix, qm_gb%gb_mmpot, qm_gb%gb_qmpot)
        call timer_stop(TIME_QMMMENERGYSCFFOCKGB)
      end if
      call timer_stop(TIME_QMMMENERGYSCFFOCK)
      !End step 6

      first_iteration = .false.
    end do do_scf ! do scf_iteration, qmmm_nml%itrmax 
! END MAIN SCF LOOP

!We get to this point in the code for 2 reasons.
! 1) because we exceeded the maximum number of scf iterations.
! 2) because we converged and so broke out of the above loop.
! Check which condition it is based on the value of the logical converged.
! Only the master thread returns the scf energy to force.
  if (qmmm_mpi%commqmmm_master) then 
    !Here this is just the electronic energy.
    escf = scf_energy
  else
    escf = zero
  end if

  if (.NOT. converged .and. qmmm_mpi%commqmmm_master) then
    !Convergence failed. Print a warning message and return with
    !the current unconverged density. This will mean the forces
    !are not accurate but in an MD simulation they may be good
    !enough to allow the user to get out of a potentially bad geometry.
    write(6,'(/,''QMMM: WARNING!'')')
    write(6,'(''QMMM: Unable to achieve self consistency to the tolerances specified'')')
    write(6,'(''QMMM: No convergence in SCF after '',i6,'' steps.'')') qmmm_nml%itrmax
    write(6,'(''QMMM: Job will continue with unconverged SCF. Warning energies'')')
    write(6,'(''QMMM: and forces for this step will not be accurate.'')')
    write(6,'(''QMMM: E = '',E12.4,'' DeltaE = '',E12.4,'' DeltaP = '',E12.4)') escf, energy_diff, density_diff
    write(6,'(''QMMM: Smallest DeltaE = '',E12.4,'' DeltaP = '',E12.4,'' Step = '',i6,/)') &
                      smallest_energy_diff(1), smallest_energy_diff(2), sm_energy_diff_step_number
  end if

  if (qmmm_nml%verbosity > 0 .and. converged .and. qmmm_mpi%commqmmm_master) then
     write(6,'("QMMM: SCF Converged to ",G10.4," in: ",i5," Cycles ")') qmmm_nml%scfconv,scf_iteration
  end if

  return
end subroutine qm2_scf

subroutine qm2_densit( eigen_vecs,norbs,ndubl, den_matrix,matsize)
!***********************************************************************        
!   
!   Computes the density matrix from the eigen vector matrix and the
!   MO occupancy.
!                                                                               
!   eigen_vecs = square eigen vector matrix, of size norbs BY norbs         
!                eigen vectors are stored in top left corner.
!        norbs = no. orbitals
!        ndubl = number of doubly occupied M.O.S.
!
!    Not currently used:
!        nsingl= no. double+single occupied M.O.S.
!                                                                               
!   returns den_matrix  = density matrix
!
! Routine and optimization / Open MP by Ross Walker (SDSC, 2007)                              
!                                                                               
!***********************************************************************        
!

      use qmmm_module, only : qm2_params
      implicit none

! Passed in
      integer, intent(in) :: norbs, ndubl, matsize
!      integer, intent(inout) :: nsingl
      _REAL_, intent(in) :: eigen_vecs(norbs,norbs)
      _REAL_, intent(out) :: den_matrix(matsize)

! Local
      integer unoc_start,i,j,k,l
      _REAL_ eigen_veci

!      integer :: sing_oc_start
!      !NSINGL=MAX(NDUBL,NSINGL)
!      unoc_start=NSINGL+1

      unoc_start = ndubl + 1

!if (ndubl == nsingl) then
!Currently only spin=1 is allowed.
!All  M.O.s doubly occupied.

#ifdef QMMM_OMP
!$OMP PARALLEL &
!$OMP DEFAULT(PRIVATE) &
!$OMP SHARED(den_matrix,matsize,unoc_start,norbs,qm2_params,eigen_vecs)
! den_matrix can be shared since no two threads will do the same value of L.
!$OMP DO SCHEDULE(static)
#endif
      do i = 1,matsize
        den_matrix(i)=0.0d0
      end do
#ifdef QMMM_OMP
!$OMP END DO
#endif

      do k=unoc_start,norbs
#ifdef QMMM_OMP
!$OMP DO SCHEDULE(guided)
        do i=1,norbs
          L=qm2_params%pascal_tri1(i)
#else
        L = 0
        do i=1,norbs
#endif
          eigen_veci = eigen_vecs(i,k)
          do j=1,i
            L=L+1
            den_matrix(L)=den_matrix(L)-2.0d0*eigen_veci*eigen_vecs(j,k)
          end do
          if (k==unoc_start) den_matrix(L)=2.0d0+den_matrix(L)
        end do
#ifdef QMMM_OMP
!$OMP END DO
#endif
      end do

#ifdef QMMM_OMP
!$OMP END PARALLEL
#endif

      return

end subroutine qm2_densit

!-----------------------------------------------------------------------

subroutine qm2_cnvg(den_matrix, old_den_matrix, old2_density,norbs, scf_iteration, density_diff)
!***********************************************************************        
!                                                                               
!  CNVG IS A TWO-POINT INTERPOLATION ROUTINE FOR SPEEDING CONVERGENCE           
!       OF THE DENSITY MATRIX.                                                  
!                                                                               
! ON OUTPUT den_matrix   = NEW DENSITY MATRIX                                    
!           old2_density  = DIAGONAL OF OLD DENSITY MATRIX                             
!           density_diff     = LARGEST DIFFERENCE BETWEEN OLD AND NEW DENSITY             
!                    MATRIX DIAGONAL ELEMENTS                                   
!***********************************************************************        
      use qmmm_module, only : qm2_params
      implicit none

!Passed in
      integer, intent(in) :: norbs, scf_iteration
      _REAL_, intent(inout) :: old2_density(norbs), old_den_matrix(*), den_matrix(*)
      _REAL_, intent(out) :: density_diff

!Local
      _REAL_ faca, damp, facb, facb_temp, fac, sum0, current_den_sum, sum2, sum3
      _REAL_ den_element_diff, den_matrix_element, old_den_matrix_element, old2_density_element
      integer iminus,k,i,j,ie

      if (scf_iteration > 4) then
        DAMP=0.05D0
      else
        DAMP=1.0d10
      endif

      density_diff=0.0D00
      FACA=0.0D00
      FACB=0.0D00
      FAC=0.0D00
      current_den_sum=0.0D0

#ifdef QMMM_OMP
!$OMP PARALLEL &
!$OMP DEFAULT(PRIVATE) &
!$OMP SHARED(norbs, qm2_params, scf_iteration, den_matrix, old_den_matrix, old2_density, current_den_sum, density_diff, FACA, FACB, FAC, DAMP, SUM2, SUM0, SUM3)
#endif
      if (MOD(scf_iteration-1,3) /= 0) then
#ifdef QMMM_OMP
!$OMP DO SCHEDULE(static) REDUCTION(+:current_den_sum) REDUCTION(max:density_diff)
        do i =1, norbs
           k=qm2_params%pascal_tri2(i)
#else
        K=0
        do I=1,norbs
           K=K+I !pascal_tri2
#endif
           den_matrix_element=den_matrix(K)
           old_den_matrix_element=old_den_matrix(K)
           old2_density(I)=old_den_matrix_element
           density_diff=max(ABS(den_matrix_element-old_den_matrix_element),density_diff)
           old_den_matrix(K)=den_matrix_element
           current_den_sum=current_den_sum+den_matrix_element
        end do
#ifdef QMMM_OMP
!$OMP END DO
#endif
      else
#ifdef QMMM_OMP
!$OMP DO SCHEDULE(static) REDUCTION(+:current_den_sum, faca, facb) REDUCTION(max:density_diff)
        do i =1, norbs
           k=qm2_params%pascal_tri2(i)
#else
        K=0
        do I=1,norbs
           K=K+I !pascal_tri2
#endif
           den_matrix_element=den_matrix(K)
           old_den_matrix_element=old_den_matrix(K)

           current_den_sum=current_den_sum+den_matrix_element
           den_element_diff=ABS(den_matrix_element-old_den_matrix_element)
           density_diff=max(den_element_diff,density_diff)
  
           FACA=FACA+den_element_diff*den_element_diff
           facb_temp = den_matrix_element-2.D00*old_den_matrix_element+old2_density(I)
           FACB=FACB+facb_temp*facb_temp

           old2_density(I)=old_den_matrix_element
           old_den_matrix(K)=den_matrix_element
        end do
!$OMP END DO
      end if
#ifdef QMMM_OMP
!$OMP SINGLE
#endif

      if (FACB > 0.0D00) then
        IF (FACA < (100.D00*FACB)) FAC=SQRT(FACA/FACB)
      end if

      !do i=1 case of loop below
      IF(ABS(old_den_matrix(1)-old2_density(1)) > DAMP) THEN
         old_den_matrix(1)=old2_density(1)+SIGN(DAMP,old_den_matrix(1)-old2_density(1))
      ELSE
         old_den_matrix(1)=old_den_matrix(1)+FAC*(old_den_matrix(1)-old2_density(1))
      ENDIF
      old_den_matrix(1)=MIN(2.0D0,MAX(old_den_matrix(1),0.D0))
      SUM2=old_den_matrix(1)
      den_matrix(1)=old_den_matrix(1)
#ifdef QMMM_OMP
!$OMP END SINGLE
!$OMP DO SCHEDULE(guided) REDUCTION(+:sum2)
      do i=2,norbs
        ie=qm2_params%pascal_tri1(i)
#else
      IE=1
      do I=2,norbs
#endif
         iminus=I-1
         do J=1,iminus
            IE=IE+1
            den_matrix_element=den_matrix(IE)
            old_den_matrix(IE)=den_matrix_element+FAC*(den_matrix_element-old_den_matrix(IE))
            den_matrix(IE)=old_den_matrix(IE)
         end do
         IE=IE+1
         old2_density_element = old2_density(I)
         old_den_matrix_element = old_den_matrix(IE)
         IF(ABS(old_den_matrix_element-old2_density_element) > DAMP) THEN
            old_den_matrix_element=old2_density_element+SIGN(DAMP,old_den_matrix_element-old2_density_element)
         ELSE
            old_den_matrix_element=old_den_matrix_element+FAC*(old_den_matrix_element-old2_density_element)
         ENDIF
         old_den_matrix_element=MIN(2.0D0,MAX(old_den_matrix_element,0.D0))
         old_den_matrix(IE)=old_den_matrix_element
         den_matrix(IE)=old_den_matrix_element
         SUM2=SUM2+old_den_matrix_element
      end do
#ifdef QMMM_OMP
!$OMP END DO
!$OMP END PARALLEL
#endif
!   RE-NORMALIZE IF ANY DENSITY MATRIX ELEMENTS HAVE BEEN TRUNCATED             

      SUM0=current_den_sum
      do
         IF(SUM2 > 1.D-3)THEN
           SUM3=current_den_sum/SUM2
         ELSE
           SUM3=0.D0
         ENDIF
         current_den_sum=SUM0
         IF(SUM2 <= 1.D-3 .OR. ABS(SUM3-1.D0) <= 1.D-5) exit
         SUM2=0.D0
         J=0
         do I=1,norbs
            J=J+I
            !J=qm2_params%pascal_tri2(i)
!  1.0d-20 is in case any occupancy is exactly zero
            old_den_matrix_element=max(old_den_matrix(J)*SUM3+1.D-20, 0.0d0)
                                                                               
!  SET UP RENORMALIZATION OVER PARTLY OCCUPIED M.O.'S ONLY.  FULL M.O.'S        
!  CAN'T BE FILLED ANY MORE                                                     
                                                                               
            IF(old_den_matrix_element > 2.0D0)THEN
               old_den_matrix_element=2.0D0
               current_den_sum=current_den_sum-2.0D0
            ELSE
               SUM2=SUM2+old_den_matrix_element
            ENDIF
            old_den_matrix(J)=old_den_matrix_element
            den_matrix(J)=old_den_matrix_element
         end do
      end do   

      RETURN
end subroutine qm2_cnvg

!---------------------------------------------------------------------------

subroutine qm2_mat_diag(a,n,m,v,e,w1,w2,w3,w4,w5)

      use qmmm_module, only : qm2_params
      use constants, only : zero, one, half
      implicit none

!Passed in
      integer, intent(in) :: m,n
      _REAL_, intent(inout) :: a(*) !Note original matrix will be corrupted by this routine.
      _REAL_, intent(out) ::  e(n), v(n,m)
      _REAL_, intent(out) ::  w1(n), w2(n), w3(n), w4(n), w5(n)

!**********************************************************************
!
! qm2_mat_diag is a diagonalisation routine, based on:
! Yoshitaka Beppu of Nagoya University, Japan.
!       For details see 'Computers & Chemistry' vol.6 1982. page 000.
!
! This version has been modernised and optimised by
! Ross Walker and David Case (TSRI 2005)
!
! on input    a       = matrix to be diagonalised (packed canonical)
!             n       = size of matrix to be diagonalised.
!             m       = number of eigenvectors needed.
!             e       = array of size at least n
!             v       = array of size at least nmax*m
!
! on output   e       = eigenvalues
!             v       = eigenvectors in array of size nmax*m
!
!***********************************************************************

!Local
      _REAL_ ssum, r, s, h, summ, ff, ww, hinv, rinv
      _REAL_ c, gersch, z, sinv, ee, ra, t, rn, vn, del, u
      integer nm1, nm2, irank, jrank, krank, i, j, k, l, kp1, ip1
      integer im1, ii, kk, kpiv, ig, itere, ll

      ! eps and eps3 are machine-precision dependent
      _REAL_, parameter :: eps=1.d-8
      _REAL_, parameter :: eps3=1.d-30

      ! Householder transformation

      nm1=n-1
      nm2=n-2
      krank=0
      do k=1,nm2
         kp1=k+1
         krank=krank+k
         w2(k)=a(krank)
         ssum=zero
         jrank=krank
         do j=kp1,n
            w2(j)=a(jrank+k)
            jrank=jrank+j
            ssum=w2(j)*w2(j)+ssum
         end do
         s=sign(sqrt(ssum),w2(kp1))
         w1(k)=-s
         w2(kp1)=w2(kp1)+s
         a(k+krank)=w2(kp1)
         h=w2(kp1)*s
         if(abs(h) < eps3) then
            a(krank) = h
         else
           hinv=one/h
           summ=zero
           irank=krank
           do i=kp1,nm1 !kp1->n-1
              ssum=zero
              do j=kp1,i
                 ssum=ssum+a(j+irank)*w2(j)
              end do
              ip1=i+1
!              jrank=ishft(i*(i+3),-1)
!              jrank=i*(i+3)/2
              jrank=qm2_params%pascal_tri2(i)+i
              do j=ip1,n
                 ssum=ssum+a(jrank)*w2(j)
                 jrank=jrank+j
              end do
              w1(i)=ssum*hinv
              irank=irank+i
              summ=w1(i)*w2(i)+summ
           end do
           !Now do the case for i=n
           ssum=zero
           do j=kp1,n
              ssum=ssum+a(j+irank)*w2(j)
           end do
           w1(n)=ssum*hinv
           summ=w1(n)*w2(n)+summ

           u=summ*half*hinv
           jrank=krank
           do j=kp1,n
              w1(j)=w2(j)*u-w1(j)
              do i=kp1,j
                 a(i+jrank)=w1(i)*w2(j)+w1(j)*w2(i)+a(i+jrank)
              end do
              jrank=jrank+j
           end do
           a(krank)=h
         end if !if (abs(h) < eps3)
      end do

      !w2(nm1)=a( ishft((nm1*(nm1+1)),-1) )
      w2(nm1)=a( qm2_params%pascal_tri2(nm1) )
      !w2(n)=a( ishft((n*(n+1)),-1) )
      w2(n)=a( qm2_params%pascal_tri2(n))
      !w1(nm1)=a( nm1+ishft((n*(n-1)),-1) )
      w1(nm1)=a( nm1+qm2_params%pascal_tri1(n) )
      w1(n)=zero
      gersch=abs(w2(1))+abs(w1(1))
      do i=1,nm1
         gersch=max(abs(w2(i+1))+abs(w1(i))+abs(w1(i+1)),gersch)
      end do
      del=eps*gersch
      w3(1:n)=w1(1:n)
      e(1:n)=w2(1:n)
      v(1:n,m)=e(1:n)
      if(abs(del) >= eps3) then

         ! QR-method with origin shift

         k=n
         do 
            l=k
            do
               if(abs(w3(l-1)) < del) exit
               l=l-1
               if(l == 1)  exit
            end do
            if(l /= k) then
               ww=(e(k-1)+e(k))*half
               r=e(k)-ww
               z=sign(sqrt(w3(k-1)**2+r*r),r)+ww
               ee=e(l)-z
               e(l)=ee
               ff=w3(l)
               rinv=one/sqrt(ee*ee+ff*ff)
               c=e(l)*rinv
               s=w3(l)*rinv
               ww=e(l+1)-z
               e(l)=(ff*c+ww*s)*s+ee+z
               e(l+1)=c*ww-s*ff
               do j=l+1,k-1
                 rinv=one/sqrt(e(j)*e(j)+w3(j)*w3(j))
                 r = one/rinv
                 w3(j-1)=s*r
                 ee=e(j)*c
                 ff=w3(j)*c
                 c=e(j)*rinv
                 s=w3(j)*rinv
                 ww=e(j+1)-z
                 e(j)=(ff*c+ww*s)*s+ee+z
                 e(j+1)=c*ww-s*ff
               end do
               w3(k-1)=e(k)*s
               e(k)=e(k)*c+z
               cycle
            end if
            k=k-1
            if(k == 1) exit
         end do

        ! At this point the array 'e' contains the un-ordered eigenvalues
        ! Straight selection sort of eigenvalues:

         j=n
         do
            l=1
            ii=1
            ll=1
            do i=2,j
              if( (e(i)-e(l)) >= zero ) then
                 l=i
              else
                ii=i
                ll=l
              end if
            end do
!           if(ii /= ll) then
            ww=e(ll)
            e(ll)=e(ii)
            e(ii)=ww
!           end if
            j=ii-1
            if(j < 2) exit
         end do
      end if

      ! Ordering of eigenvalues is complete.

      ! Inverse-iteration for eigenvectors:

      rn=zero
      ra=eps*0.6180339887485d0
                   ![0.618... is the fibonacci number (-1+sqrt(5))/2.]
      ig=1
      do i=1,m
         im1=i-1
         do j=1,n
            w3(j)=zero
            w4(j)=w1(j)
            w5(j)=v(j,m)-e(i)
            rn=rn+ra
            if(rn.ge.eps) rn=rn-eps
            v(j,i)=rn
         end do
         do j=1,nm1
            if(abs(w5(j)) < abs(w1(j))) then
               w2(j)=-w5(j)/w1(j)
               w5(j)=w1(j)
               t=w5(j+1)
               w5(j+1)=w4(j)
               w4(j)=t
               w3(j)=w4(j+1)
               if(abs(w3(j)).lt.eps3) w3(j)=del
               w4(j+1)=zero
            else
               if(abs(w5(j)).lt.eps3) w5(j)=del
               w2(j)=-w1(j)/w5(j)
            end if
            w4(j+1)=w3(j)*w2(j)+w4(j+1)
            w5(j+1)=w4(j)*w2(j)+w5(j+1)
         end do
         if(abs(w5(n)) < eps3) w5(n)=del
         !itere=1
         v(n,i)=v(n,i)/w5(n)
         v(nm1,i)=(v(nm1,i)-v(n,i)*w4(nm1))/w5(nm1)
         vn=max(abs(v(n,i)),abs(v(nm1,i)),1.0d-20)
         do k=nm2,1,-1 !n-2,1,-1
            v(k,i)=(v(k,i)-v(k+1,i)*w4(k)-v(k+2,i)*w3(k))/w5(k)
            vn=max(abs(v(k,i)),vn,1.d-20)
         end do
         s=1.0d-5/vn
         v(1:n,i)=v(1:n,i)*s

         do itere=2,5
            do j=1,nm1
               if(abs(w3(j)) >= eps3) then
                 t=v(j,i)
                 v(j,i)=v(j+1,i)
                 v(j+1,i)=t
                 v(j+1,i)=v(j,i)*w2(j)+v(j+1,i)
               end if
            end do
            v(n,i)=v(n,i)/w5(n)
            v(nm1,i)=(v(nm1,i)-v(n,i)*w4(nm1))/w5(nm1)
            vn=max(abs(v(n,i)),abs(v(nm1,i)),1.0d-20)
            do k=nm2,1,-1 !n-2,1,-1
               v(k,i)=(v(k,i)-v(k+1,i)*w4(k)-v(k+2,i)*w3(k))/w5(k)
               vn=max(abs(v(k,i)),vn,1.d-20)
            end do
            s=1.0d-5/vn
            v(1:n,i)=v(1:n,i)*s
            if(vn>1.0d0) exit
         end do

         ! Transformation of eigenvectors

         krank=ishft(nm2*(n+1),-1)
         kpiv=ishft(nm2*nm1,-1)
         do k=nm2,1,-1 !n-2,1,-1
            kp1=k+1
            if(abs(a(kpiv)) > eps3) then
               ssum=zero
               do kk=kp1,n
                  ssum=ssum+a(krank)*v(kk,i)
                  krank=krank+kk
               end do
               s=-ssum/a(kpiv)
               do kk=n,kp1,-1
                  krank=krank-kk
                  v(kk,i)=a(krank)*s+v(kk,i)
               end do
            end if
            kpiv=kpiv-k
            krank=krank-kp1
         end do
         do j=ig,im1 !i-1
            if( abs(e(j)-e(i)) < 0.05d0 ) exit
         end do
         ig=j
         if(ig /= i) then
            ! re-orthogonalisation:

            do k=ig,im1
               ssum=zero
               do j = 1,n
                 ssum=v(j,k)*v(j,i)+ssum
               end do
               s=-ssum
               v(1:n,i)=v(1:n,k)*s+v(1:n,i)
            end do
         end if

         ! Normalisation:

         ssum=1.d-24
         do j=1,n
            ssum=ssum+v(j,i)**2
         end do
         sinv=one/sqrt(ssum)
         v(1:n,i)=v(1:n,i)*sinv
      end do

      return
end subroutine qm2_mat_diag

!---------------------------------------------------------------------------

subroutine qm2_pseudo_diag(matrix,vectors,noccupied,eigen,norbs,smallsum, &
                           matrix_workspace, scratch_matrix, vectmp1, vectmp2, vectmp3, vecjs)

!------------------------------------------------
!  "FAST" but "APPROXIMATE" matrix diagonalisation.
!  See: Stewart, J.J.P., Csaszar, P., Pulay, P., J. Comp. Chem. 3:227,1982
!
!  Written by Ross Walker (TSRI, 2005)
!  Vector Version by Ross Walker (SDSC, 2006)
!  OpenMP Parallel Version by Ross Walker (SDSC, 2006)
!
!------------------------------------------------

!  qm2_pseudo_diag is a "FAST" but "approximate" matrix diagonalisation routine.
!  The vectors generated by this routine are more likely to be able to 
!  block-diagonalise the Fock matrix over the molecular orbitals than the 
!  starting vectors. It is based on the work by Stewart et al. and must be 
!  considered to be a pseudo diagonaliser because:
!
! 1) Eigenvectors are not generated. Only the occupied-virtual intersection is
!    diagonalised, not the secular determinant.
!
! 2) When rotation is used to eliminate elements of the secular determinant the
!    remainder is assumed to remain unchanged. Any elements that are created
!    are thus ignored.
!
! 3) The rotation required to eliminate those elements considered significant 
!    is approximated to using the eigenvalues of the exact diagonalisation 
!    throughout the rest of the iterative procedure. In other words the 
!    errors on each pseudo-step are propogated to the next SCF step. The errors
!    here are assumed to be smaller than the actual error in the SCF at this 
!    point in the SCF.
!

!Variables:
!
!  MATRIX - On input should contain the matrix to be diagonalised in the 
!           form of a packed lower half triangle.
!
  use qmmm_module, only : qm2_struct, qmmm_mpi
  implicit none

! Passed in
  _REAL_, intent(in) :: matrix(qm2_struct%matsize)
  integer, intent(in) :: noccupied, norbs
  _REAL_, intent(inout) :: vectors(norbs,norbs)
  _REAL_, intent(in) :: eigen(norbs)
  _REAL_, intent(in) :: smallsum

!Passed in Scratch Arrays
  _REAL_, intent(out) :: matrix_workspace(norbs,norbs)
  _REAL_, intent(out) :: scratch_matrix(noccupied,norbs)
  _REAL_, intent(out) :: vectmp1(noccupied*(norbs-noccupied))
  _REAL_, intent(out) :: vectmp2(noccupied*(norbs-noccupied))
  _REAL_, intent(out) :: vectmp3(noccupied*(norbs-noccupied))
  integer, intent(out) :: vecjs(2,noccupied*(norbs-noccupied))

! Local
  _REAL_ sum1,a,b,c,d,alpha,beta
  _REAL_ eigeni
  integer i,j,k,j1,k2,kk,m,lumo
  integer :: ii, jj, veccount

  lumo=noccupied+1
#ifdef QMMM_OMP
!$OMP PARALLEL &
!$OMP DEFAULT(PRIVATE) &
!$OMP SHARED(lumo, norbs, noccupied, matrix, vectors, scratch_matrix, matrix_workspace, veccount, eigeni, eigen, c,d, smallsum, vectmp1,vectmp2,vectmp3,vecjs)
!workspace can be shared for OMP since no two threads should do the same value of i.
!$OMP DO SCHEDULE(guided)
#endif
  do i = lumo,norbs
    kk=0
    do j=1,norbs
      sum1 = 0.0d0
      do k=1,j
        kk=kk+1
        sum1=sum1+matrix(kk)*vectors(k,i)
      end do
      j1 = j+1
      k2 = kk
      do k=j1,norbs
        k2=k2+k-1
        sum1=sum1+matrix(k2)*vectors(k,i)
      end do !k=j1,norbs
      matrix_workspace(j,i) = sum1
    end do !j=1,norbs
  end do
#ifdef QMMM_OMP
!$OMP END DO

!$OMP DO SCHEDULE(static)
#endif
  do i = lumo,norbs
    do j=1,noccupied
      scratch_matrix(j,i) = sum(matrix_workspace(1:norbs,i)*vectors(1:norbs,j))
    end do !j=1,noccupied
  end do !i=lumo,norbs
#ifdef QMMM_OMP
!$OMP END DO
!$OMP SINGLE
#endif
  !Vectored Version by Ross Walker
  veccount=0
  do i=lumo,norbs
    eigeni=eigen(i)
    do j=1,noccupied
      C=scratch_matrix(j,i) !Filling this array in a seperate loop above is faster in most cases
      D=eigen(j)-eigeni     !than doing everything here as a single loop because of dependencies.
      !Check the machine precision for whether to do a 2x2 rotation
      if (abs(C) >= (smallsum*abs(D))) then
        veccount = veccount+1
        vectmp1(veccount) = c
        vectmp2(veccount) = d
        vecjs(1,veccount) = i
        vecjs(2,veccount) = j
      end if
    end do
  end do
#ifdef QMMM_OMP
!$OMP END SINGLE
!$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) = 4.0d0*vectmp1(i)*vectmp1(i) + vectmp2(i)*vectmp2(i)
  end do 
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) = 1.0d0/sqrt(vectmp3(i))
  end do
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1,veccount
    vectmp3(i) =  0.5d0 + abs(0.5d0*vectmp2(i)*vectmp3(i))
  end do
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = 1.0d0-vectmp3(i)
  end do
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = sqrt(vectmp2(i))
  end do
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp2(i) = -sign(vectmp2(i),vectmp1(i))
  end do
!$OMP END DO
!$OMP DO SCHEDULE(static)
  do i = 1, veccount
    vectmp3(i) = sqrt(vectmp3(i))
  end do
!$OMP END DO
!$OMP END PARALLEL
#else
  vectmp3(1:veccount) = 4.0d0*vectmp1(1:veccount)*vectmp1(1:veccount)+ &
                        vectmp2(1:veccount)*vectmp2(1:veccount)
  call vdinvsqrt(veccount, vectmp3, vectmp3) !E=1.0D0/sqrt(4.0D0*c*c+d*d)
  vectmp3(1:veccount) = 0.5d0 + abs(0.5d0*vectmp2(1:veccount)*vectmp3(1:veccount)) !E=0.5D0+abs(0.5d0*D*E)
  vectmp2(1:veccount) = 1.0d0-vectmp3(1:veccount)
  call vdsqrt(veccount, vectmp2, vectmp2)
  vectmp2(1:veccount) = -sign(vectmp2(1:veccount),vectmp1(1:veccount)) !beta = -sign(sqrt(1.0d0-E),c)
  call vdsqrt(veccount, vectmp3, vectmp3) !alpha = sqrt(E)
  !vectmp2 = beta
  !vectmp3 = alpha
#endif

!-------------------------------------------------------------------
! Now we have done the squaring we can do a crude 2 by 2 rotation,
! which we will assume eliminates the significant elements.
!-------------------------------------------------------------------

  do i=1, veccount
    ii = vecjs(1,i)
    jj = vecjs(2,i)
    beta = vectmp2(i)
    alpha = vectmp3(i)
    do m=1, norbs
      a=vectors(m,jj)  !jj is not part of the set of ii.
      b=vectors(m,ii)
      vectors(m,jj)=alpha*a+beta*b
      vectors(m,ii)=alpha*b-beta*a
    end do !m=1,norb
  end do !i=1,veccount

  return
end subroutine qm2_pseudo_diag

function qm2_HELECT(nminus,den_matrix,H,F)
!***********************************************************************
!
!    SUBROUTINE CALCULATES THE ELECTRONIC ENERGY OF THE SYSTEM IN EV.
!
!    ON ENTRY Nminus = NUMBER OF ATOMIC ORBITALS - 1
!             den_matrix = DENSITY MATRIX, PACKED, LOWER TRIANGLE.
!             H = ONE-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
!             F = TWO-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
!    ON EXIT
!        qm2_HELECT = ELECTRONIC ENERGY.
!
!***********************************************************************
   implicit none

   _REAL_, intent(in) :: den_matrix(*), H(*), F(*)
   integer, intent(in) :: nminus

   _REAL_ ed, qm2_helect
   integer k,i,j

   ED=0.0D00
   qm2_helect=0.0D00
   K=0
   do I=1,nminus
      K=K+1
      ED=ED+den_matrix(K)*0.5D0*(H(K)+F(K))
      do J=1,I
         K=K+1
         qm2_helect=qm2_helect+den_matrix(K)*(H(K)+F(K))
      end do
   end do

   K=K+1
   ED=ED+den_matrix(K)*0.5D0*(H(K)+F(K))

   qm2_helect=qm2_helect+ED
   RETURN
end function qm2_helect

subroutine qm2_full_diagonalize(diag_routine,matrix,matrix_dimension,eigen_vectors,abstol)
!******************************
!This routine is a central driver routine for doing the diagonalization,
!essentially it can be called with the matrix to be diagonalized and the
!type of diagonalization to be done. Note it assumes the matrix and
!everything else is laid out correctly in memory and the relevant
!timers have been started or stopped.
!******************************

  use qmmm_module, only : qmmm_scratch, qm2_struct

  implicit none

!Passed in
  integer, intent(inout) :: diag_routine !Controls the diagonalization method to be used.
                                      !1=built in diagonalizer
                                      !2=dspev
                                      !3=dspevd
                                      !4=dspevx
                                      !5=dsyev
                                      !6=dsyevd

  _REAL_, intent(inout) :: matrix(*)  !The matrix to be diagonalized
  integer, intent(in) :: matrix_dimension !The dimension of the matrix (1D).
  _REAL_, intent(out) :: eigen_vectors(matrix_dimension,matrix_dimension) !Eigen vectors
  _REAL_, intent(in) :: abstol !Machine underflow limit for dspevr

!Local
  integer :: ierr, eigen_value_count

  if (diag_routine == 1) then
     ! Note: on our test machine - Intel P4 and Altix this qm2_mat_diag
     ! routine is generally faster than the LAPACK diagonaliser for up
     ! to around 90 QM atoms.

     ! For systems larger than 90 QM atoms the LAPACK diagonaliser
     ! may be quicker but thorough testing on your specific machine,
     ! compiler and LAPACK library combination is necessary
     ! to find the crossover point.

     !Dimension 1 of mat_diag_workspace contains the eignvalues
     call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )

  else if (diag_routine == 2) then
     !use dspev
     !expects mat diag workspace of (norbs,4) - norbs,1 used for eigenvalues, norbs,2->4
     !is used for scratch space.
     call dspev('V','U',matrix_dimension,matrix, &
                qmmm_scratch%mat_diag_workspace(1,1), &
                eigen_vectors, matrix_dimension, &
                qmmm_scratch%mat_diag_workspace(1,2), &
                ierr)
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dspev failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if

  else if (diag_routine == 3) then
     !use a divide and conquer lapack diagonaliser
     call dspevd('V','U',matrix_dimension,matrix, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          eigen_vectors, matrix_dimension, &
          qmmm_scratch%lapack_dc_real_scr, &
          qmmm_scratch%lapack_dc_real_scr_aloc, &
          qmmm_scratch%lapack_dc_int_scr,&
          qmmm_scratch%lapack_dc_int_scr_aloc, ierr)
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dspevd failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if

  else if (diag_routine == 4) then
     !Use dspevx
     call dspevx('V','A','U',matrix_dimension,matrix, &
                 0,0,0,0,abstol, &
                 eigen_value_count, &
                 qmmm_scratch%mat_diag_workspace(1,1), &
                 eigen_vectors, matrix_dimension, &
                 qmmm_scratch%lapack_dc_real_scr, &
                 qmmm_scratch%lapack_dc_int_scr(matrix_dimension+1), &
                 qmmm_scratch%lapack_dc_int_scr, ierr )
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dspevx failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if
  else if (diag_routine == 5) then
     !use dsyev - this takes a matrix that is not in packed storage so 
     !since everything is currently done in packed storage we need to unpack the matrix
     !first. We will unpack the matrix into the eigen vectors array since dsyev overwrites
     !the matrix with the eigen_vectors anyway.
     call qm2_unpack_matrix(matrix, eigen_vectors, matrix_dimension) 
     call dsyev('V','U',matrix_dimension,eigen_vectors,matrix_dimension, &
                        qmmm_scratch%mat_diag_workspace(1,1), &
                        qmmm_scratch%lapack_dc_real_scr, &
                        qmmm_scratch%lapack_dc_real_scr_aloc, &
                        ierr)
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dsyev failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if
  else if (diag_routine == 6) then
     !use dsyevd - this takes a matrix that is not in packed storage so 
     !since everything is currently done in packed storage we need to unpack the matrix
     !first. We will unpack the matrix into the eigen vectors array since dsyevd overwrites
     !the matrix with the eigen_vectors anyway.
     call qm2_unpack_matrix(matrix, eigen_vectors, matrix_dimension)
     call dsyevd('V','U',matrix_dimension,eigen_vectors, matrix_dimension, &
                        qmmm_scratch%mat_diag_workspace(1,1),& 
                        qmmm_scratch%lapack_dc_real_scr, &
                        qmmm_scratch%lapack_dc_real_scr_aloc, &
                        qmmm_scratch%lapack_dc_int_scr, &
                        qmmm_scratch%lapack_dc_int_scr_aloc, &
                        ierr)
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dsyevd failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if
  else if (diag_routine == 7) then
     !use dsyevr - this takes a matrix that is not in packed storage so we first need to unpack it.
     !             we can't use the eigenvectors array here since dsyevr expects that in a different
     !             memory space. Thus we will use pdiag_scr_norbs_norbs array.
     call qm2_unpack_matrix(matrix, qmmm_scratch%pdiag_scr_norbs_norbs, matrix_dimension)
     call dsyevr('V','A','U',matrix_dimension,qmmm_scratch%pdiag_scr_norbs_norbs,matrix_dimension, &
                 0.0d0, 0.0d0, 0, 0, abstol, eigen_value_count, qmmm_scratch%mat_diag_workspace(1,1), &
                 eigen_vectors, matrix_dimension, qmmm_scratch%lapack_dc_int_scr, &
                 qmmm_scratch%lapack_dc_real_scr, &
                 qmmm_scratch%lapack_dc_real_scr_aloc, &
                 qmmm_scratch%lapack_dc_int_scr(2*matrix_dimension+1), &
                 (qmmm_scratch%lapack_dc_int_scr_aloc-2*matrix_dimension), ierr)
     !Check ierr returned zero.
     if (ierr /= 0) then
       write(6,'("| QMMM: ERROR dsyevr failed to converge and returned ",i8)') ierr
       write(6,'("| QMMM: Falling back on internal diagonalizer for the remainder of this run.")')
       diag_routine = 1
       deallocate (qmmm_scratch%mat_diag_workspace,stat=ierr)
       REQUIRE(ierr==0)
       allocate (qmmm_scratch%mat_diag_workspace(qm2_struct%norbs,6),stat=ierr)
       REQUIRE(ierr==0)
       call qm2_mat_diag(matrix,matrix_dimension,matrix_dimension,eigen_vectors, &
          qmmm_scratch%mat_diag_workspace(1,1), &
          qmmm_scratch%mat_diag_workspace(1,2), &
          qmmm_scratch%mat_diag_workspace(1,3), &
          qmmm_scratch%mat_diag_workspace(1,4), &
          qmmm_scratch%mat_diag_workspace(1,5), &
          qmmm_scratch%mat_diag_workspace(1,6) )
     end if
  else
    !Unknown diag method
    call sander_bomb('qm2_full_diagonalize','UNKNOWN DIAGONALIZATION ROUTINE','RUN TERMINATED')
  end if

  return

end subroutine qm2_full_diagonalize

subroutine qm2_unpack_matrix(packed_matrix, unpacked_matrix, matrix_dimension)
  !takes a packed upper triangle in packed_matrix
  !and unpacks it into unpacked_matrix

  implicit none

  integer, intent(in) :: matrix_dimension !the dimension of the matrix
  _REAL_, intent(in) :: packed_matrix(*)
  _REAL_, intent(out) :: unpacked_matrix(matrix_dimension,matrix_dimension)

!local
  integer :: i, j, jcount

  jcount = 1

  do j = 1, matrix_dimension
    do i = 1, j
      unpacked_matrix(i,j) = packed_matrix(jcount+I-1)
    end do
    jcount = jcount + j
  end do

  return

end subroutine qm2_unpack_matrix


subroutine qm2_pack_matrix(packed_matrix, unpacked_matrix, matrix_dimension)
  !takes an unpacked matrix and packs it as an upper triangle

  implicit none

  integer, intent(in) :: matrix_dimension !the dimension of the matrix
  _REAL_, intent(out) :: packed_matrix(*)
  _REAL_, intent(in) :: unpacked_matrix(matrix_dimension,matrix_dimension)

!local
  integer :: i, j, jcount

  jcount = 1

  do j = 1, matrix_dimension
    do i = 1, j
      packed_matrix(jcount+I-1) = unpacked_matrix(i,j)
    end do
    jcount = jcount + j
  end do

  return

end subroutine qm2_pack_matrix
