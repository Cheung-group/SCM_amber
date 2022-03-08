! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "def_time.h"
subroutine qm2_dftb_scf(escf, elec_eng,enuclr_qmqm,scf_mchg)

   use qmmm_module, only : qmmm_nml, qmmm_mpi, qmmm_struct
   use qm2_dftb_module, only : mcharge, izp_str, fermi_str, mol, disper, espin
   use constants, only : AU_TO_EV, AU_TO_KCAL

   implicit none

!For irespa
#include "md.h" 

!Passed in
   _REAL_, intent(out) :: elec_eng, enuclr_qmqm, escf
   _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

!Local
   _REAL_ :: geseatom, total_e
   integer :: outer_scf_count, inner_scf_count, i
   logical :: scf_convgd

   total_e      = 0.0d0
   
   scf_convgd=.false.
   outer_scf_count = 0
   inner_scf_count = 0

   !---------------------
   !BEGIN OUTER SCF LOOP
   !---------------------
   do_scf_outer: do while ( (.not.scf_convgd) .and. (outer_scf_count < qmmm_nml%itrmax) )
     
      call eglcao(qmmm_struct%qm_coords,total_e,elec_eng,enuclr_qmqm,&
                  inner_scf_count, outer_scf_count, scf_convgd,mol%qmat,scf_mchg )

      if ( .not.scf_convgd ) then 
        ! Take steps to improve convergence

        ! Reset Broyden mixing by restarting the scf process in eglcao.
        if (qmmm_nml%verbosity > 0 .and. qmmm_mpi%commqmmm_master) then
           write(6,'(" QMMM SCC-DFTB: SCC-DFTB FOR STEP ",i5," DID NOT CONVERGE AFTER ",i3," cycles.")') &
                 irespa,outer_scf_count
           if (outer_scf_count < qmmm_nml%itrmax) write(6,'(" QMMM SCC-DFTB: Resetting Broyden mixing.")')
        end if

        ! In the first time, maybe the initial charges were way too wrong
        ! for this iteration, re-set qmat:
        if (outer_scf_count == qmmm_nml%dftb_maxiter .and. qmmm_mpi%commqmmm_master) then
           if (qmmm_nml%verbosity > 0) &
                 write(6,'(" QMMM SCC-DFTB: Resetting initial charges.")')
           mol%qmat(1:qmmm_struct%nquant_nlink) = mcharge%qzero( izp_str%izp(1:qmmm_struct%nquant_nlink) )
        end if

        ! Increase telec...
        ! WARNING: This is optional, and should be used only with great care. I
        !          have no idea what effects could happen.
        if (fermi_str%telec_step > 0.0d0 .and. qmmm_mpi%commqmmm_master) then
           ! Increase tlec and give a warning
           fermi_str%telec = fermi_str%telec + fermi_str%telec_step
           if (qmmm_nml%verbosity > 0) &
                 write(6,'(" QMMM SCC-DFTB: Increasing the electronic temperature of step ",i5," to ",f10.3," K")')&
                 irespa, fermi_str%telec
        end if
      end if

   end do do_scf_outer ! do while ( (.not.scf_convgd) .and. (outer_scf_count < qmmm_nml%itrmax) )

   !---------------------
   !END OF OUTER SCF LOOP
   !---------------------

   ! If we had to increase telec, now we restore it to the original value.
   if ( fermi_str%telec /= qmmm_nml%dftb_telec ) then
     if (qmmm_mpi%commqmmm_master) then
       write(6,'(" QMMM SCC-DFTB: **WARNING** : The energy for this step was calculated with a higher")')
       write(6,'(" QMMM SCC-DFTB:               elctronic temperature. You should check the results.")')
       write(6,'(" QMMM SCC-DFTB:               STEP=",i5,5X," TELEC=",f10.3)') irespa,fermi_str%telec
     end if
     ! NOT IMPLEMENTED YET:
     ! CALL PRINT_DFTB_CHARGES_AND_DIPOLE_FOR_CHECKING
     fermi_str%telec = qmmm_nml%dftb_telec
     if (qmmm_mpi%commqmmm_master) &
       write(6,'(" QMMM SCC-DFTB: Electronic temperature restored to ",f10.3," K.")') fermi_str%telec
   end if

   ! If it still didn't converge, something must be really wrong. Bomb the calculation.
   if ( outer_scf_count >= qmmm_nml%itrmax .or. .not.scf_convgd) call dftb_conv_failure("dylcao <qm2_dftb_main.f> : ", &
          "SCC Convergence failure - ITRMAX exceeded.", "Exiting")


   if (qmmm_nml%verbosity > 0 .and. qmmm_mpi%commqmmm_master) then
     if (scf_convgd) then
        write(6,'(" QMMM SCC-DFTB: SCC-DFTB for step ",i5," converged in ",i3," cycles.")') irespa,outer_scf_count
     else
        write(6,'(" QMMM SCC-DFTB: **WARNING** SCC-DFTB FOR STEP ",i5," DID NOT CONVERGE AFTER ",i3," cycles.")') &
              irespa,inner_scf_count
     end if
        write(6,*)
   end if


   if (qmmm_mpi%commqmmm_master) then

     !==============================
     !     Prints DFTB results
     !==============================
     if (qmmm_nml%verbosity > 1) then
        write(6,'(" QMMM SCC-DFTB:")')
        write(6,'(" QMMM SCC-DFTB:    Atomization Energy (a.u.)     = ",f20.12)') escf
        write(6,'(" QMMM SCC-DFTB:    Electronic Energy  (eV)       = ",f20.12)') elec_eng
        write(6,'(" QMMM SCC-DFTB:    Repulsive Energy   (eV)       = ",f20.12)') enuclr_qmqm
        if (qmmm_nml%dftb_disper ==1 )&
              write(6,'(" QMMM SCC-DFTB:    Dispersion Energy  (eV)       = ",f20.12)') disper%edis * AU_TO_EV
        write(6,'(" QMMM SCC-DFTB:    Total Energy       (eV)       = ",f20.12)') total_e*AU_TO_EV
        write(6,'(" QMMM SCC-DFTB:")')

        if ( qmmm_nml%verbosity > 3) then
           write(6,'(" QMMM SCC-DFTB:")')
           write(6,'(" QMMM SCC-DFTB:    Atomization Energy (kcal/mol) = ",f20.12)') escf
           write(6,'(" QMMM SCC-DFTB:    Electronic Energy  (a.u.)     = ",f20.12)') elec_eng
           write(6,'(" QMMM SCC-DFTB:    Repulsive Energy   (a.u.)     = ",f20.12)') enuclr_qmqm
           if (qmmm_nml%dftb_disper ==1 )&
                 write(6,'(" QMMM SCC-DFTB:    Dispersion Energy  (a.u.)     = ",f20.12)') disper%edis
           write(6,'(" QMMM SCC-DFTB:    Total Energy       (a.u.)     = ",f20.12)') total_e
           write(6,'(" QMMM SCC-DFTB:")')
        end if
     end if

     geseatom=0.0d0
     do i=1,qmmm_struct%nquant_nlink
        geseatom = geseatom+espin(izp_str%izp(i))
     end do

     ! Binding Energy. Will be put into Amber's escf energy.
     ! (This is done here mainly to shift the zero of energy,
     !  so that the DFTB energy scale is closer to the MM energy scale.)
     ! Convert results to units used by Amber
     ! (The results from dylcao come in a.u.)
     escf = (total_e-geseatom) * AU_TO_KCAL
     elec_eng  = elec_eng * AU_TO_EV
     enuclr_qmqm = enuclr_qmqm * AU_TO_EV

   end if ! (qmmm_mpi%commqmmm_master)

   return
end subroutine qm2_dftb_scf

subroutine eglcao(qm_coords,total_e,elec_eng,enuclr_qmqm, &
      inner_scf_count, outer_scf_count, scc_converged,qmat,scf_mchg)

! SUBROUTINE EGLCAO
! =================
!
! Copyright 1997 by Peter Blaudeck, Dirk Porezag, Michael Haugk,
! Joachim Elsner
!
! The original routine has been largely modified by Gustavo Seabra
! for inclusion in the Amber package. (2005)
!
! *********************************************************************
!
! PROGRAM CHARACTERISTICS
! -----------------------
!
! eglcao calculates energy and gradient for dylcao as shown by Seifert.
! The determination of the occupation numbers has been changed to be
! also valid for metallic systems.
!


!In parallel all threads enter here.

   use qm2_dftb_module, only: MDIM,LDIM,NDIM,disper, lmax, dacc, mcharge, &
         izp_str, ks_struct, fermi_str
   use qmmm_module, only: qmmm_nml, qmmm_struct, qm2_struct, qm_gb, qmewald, ELEMENT_SYM, qmmm_mpi,qm2_params
   use constants, only : BOHRS_TO_A, AU_TO_KCAL, AU_TO_EV

   implicit none


   ! Parameters passed in:
   ! =====================
   _REAL_ , intent(in ) :: qm_coords(3,qmmm_struct%nquant_nlink)       ! QM atoms coordinates
   _REAL_ , intent(out) :: total_e            ! Total energy
   _REAL_ , intent(out) :: elec_eng          ! Electronic energy
   _REAL_ , intent(out) :: enuclr_qmqm         ! Repulsive energy
   integer, intent(out) :: inner_scf_count        ! Number of SCC iterations performed in this trial
   integer, intent(out) :: outer_scf_count        ! Total number of SCC iterations performed
   logical, intent(out) :: scc_converged! SCC procedure converged?
   _REAL_ , intent(inout) :: qmat(*)    ! Electron population per atom
   _REAL_ , intent(out) :: scf_mchg(qmmm_struct%nquant_nlink) ! Mulliken charges per atom


   ! Locals
   ! ======
   integer :: indj1,indk1,lumo
   integer :: cai, liend, ljend
   integer :: j, izpj, k, li, lj, i
   integer :: n, m, norbs
   integer :: nstart, nend, mstart, mend
   integer :: indkn, indjm, indi, indj,indk, indili, indjlj
   _REAL_ :: shifti, shiftj

   !Ewald stuff
   _REAL_  :: ew_corr ! Ewald correction to energy in ev. Needed only if qm_ewald > 0.
   _REAL_  :: temp_pot
   integer :: IA, IB, i1, i2

   _REAL_  :: elec_eng_old, qtot, eext, egb
   _REAL_  :: ecoul, efermi, spro, dipabs

   ! SCC Charge convergency
   _REAL_  :: ediff
   _REAL_  :: chdiff, chdiff_tmp
   integer :: atdiff

   integer :: ier

   _REAL_ :: gb_escf_corr !GB correction for escf.

#ifdef MPI
#include "mpif.h"
#endif

   if (qmmm_mpi%commqmmm_master) then

     ! Setup of charge-independent part of H and S
     ! -------------------------------------------
     do j = 1,qmmm_struct%nquant_nlink
        do k = 1,j !--> Calculates only the lower triangle

           ! Gets the hamiltonian and overlap terms
           ! referring to this specific pair (j,k)

           call slkmatrices(j,k,qm_coords,ks_struct%au,ks_struct%bu,LDIM)

           ! Puts the calculated block inside the big
           ! hamiltonian and ovelap matrices
           nstart = 1
           nend = ks_struct%ind(k+1)-ks_struct%ind(k)

           do n = nstart, nend

              indkn = ks_struct%ind(k) + n
              mstart = 1
              mend = ks_struct%ind(j+1)-ks_struct%ind(j)

              do m = mstart, mend

                 indjm = ks_struct%ind(j) + m
                 ks_struct%hamil(indjm,indkn) = ks_struct%au(m,n)
                 ks_struct%overl(indjm,indkn) = ks_struct%bu(m,n)

                 ! The actual matrices are symmetric
                 ks_struct%hamil(indkn,indjm) = ks_struct%au(m,n)
                 ks_struct%overl(indkn,indjm) = ks_struct%bu(m,n)

              end do
           end do
        end do
     end do

     ! EXTERNAL FIELD OF MM POINT CHARGES
     ! ----------------------------------
     ! The effect of the MM atoms is taken as an extra term to 
     ! the Hamiltonian, i.e., as an external energy shift to H.
     !
     ! This is calculated only once, outside the SCC process
     ! because the charges don't move.
     !
     if ( qmmm_nml%qmmm_int > 0 ) then
        call externalshift(qm_coords,izp_str%izp,ks_struct%shiftE)
     else
        ks_struct%shiftE(1:qmmm_struct%nquant_nlink)=0.0d0
     endif

     elec_eng_old = 0.0d0
     ks_struct%shift  = 0.0d0
     ks_struct%qmold(1:qmmm_struct%nquant_nlink) = qmat(1:qmmm_struct%nquant_nlink)

   end if !(qmmm_mpi%commqmmm_master)

   call timer_start(TIME_QMMMENERGYSCF)

   scc_converged = .false.

!! =============================================================================
!!                           SCF Loop Starts here.
!! =============================================================================
   if ((qmmm_nml%verbosity > 2) .and. qmmm_mpi%commqmmm_master) then
      write (6,*)
      write (6,'(" QMMM SCC-DFTB:  SCC convergence criteria: ")')
      write (6,'(" QMMM SCC-DFTB:      Energy:", 1P,e10.1)') qmmm_nml%scfconv
      write (6,'(" QMMM SCC-DFTB:      Charge:", 1P,e10.1)') qmmm_nml%density_conv
      write (6,'(" QMMM SCC-DFTB:      Telec :", f10.3,"K")') fermi_str%telec
      write (6,'(" QMMM SCC-DFTB:  It#",20X,"Energy",19X,"Diff",14X,"MaxChDiff",3X,"At#",3X,"Index",3X,"Symb",3X,"MullikCh")')
   end if

   !For the moment the master thread does the full scf loop while the non-master
   !threads do a reduced version of it that basically just calls the qmewald code
   !in parallel and then reduces the qmpot array to the master.
   if (qmmm_mpi%commqmmm_master) then
     scc_loop: do inner_scf_count = 1,qmmm_nml%dftb_maxiter

        outer_scf_count = outer_scf_count + 1
        ks_struct%a = 0.0d0
        ks_struct%b = 0.0d0
  
        ! charge-independent part of H and S
        do j = 1,qmmm_struct%nquant_nlink
           indj  = ks_struct%ind(j)   ! Index for THIS atom (j)
           indj1 = ks_struct%ind(j+1) ! Beggining of NEXT atom on list
           do k = 1,j
              indk  = ks_struct%ind(k)   ! Index for THIS atom (k)
              indk1 = ks_struct%ind(k+1) ! Beggining of NEXT atom on list
  
              ! Copies the lower triangle 
              ! of:           into: 
              !   'H' (hamil) --> 'a'
              !   'S' (overl) --> 'b'
              ! This is needed because the matrices 'a' and 'b' are
              ! modified by the eigenvalue solver to contain, on exit:
              ! 'a' --> the eigenvectors; and 
              ! 'b' --> the "triangular factor for the Cholesky factorization"
              !         B = L*L**T. (Since DTB uses lower triangular S)
              nstart = 1
              nend   = indk1-indk

              do n = nstart,nend
                 indkn = indk + n
                 mstart = 1
                 mend   = indj1-indj

                 do m = mstart, mend
                    indjm = indj + m
                    ks_struct%a(indjm,indkn) = ks_struct%hamil(indjm,indkn)
                    ks_struct%b(indjm,indkn) = ks_struct%overl(indjm,indkn)
                 end do
              end do
           end do
        end do ! do j = 1, nquant_nlink

        ! Add charge dependent terms (Hubbard, etc. )

        ! Start by defining shift. 
        ! The 'shift' is the difference between the
        ! SCC and non-SCC energies, which is the H^1_{\mu\nu}
        ! in the secular equations.

        ! HAMILTONIAN SHIFT
        ! -----------------
        ! Calculate hamiltonian shift due to SCC - also calculate scf_mchg.

        ! zero the whole ks_struct%shift vector (qmmm_struct%nquant_nlink long)
        ks_struct%shift = 0.0d0

        call HAMILSHIFT(qm_coords, izp_str%izp,&
              mcharge%uhubb,inner_scf_count,ks_struct%gammamat, &
              ks_struct%shift, qm2_struct%scf_mchg)

        ! EXTERNAL CHARGES SHIFT
        ! ----------------------
        ! Add external charges (QM/MM)
        ! This is added as another energy shift in the Hamiltonian
        ! (Notice the different sign)

!           do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
         ks_struct%shift(1:qmmm_struct%nquant_nlink) =    &
              ks_struct%shift(1:qmmm_struct%nquant_nlink) &
            - ks_struct%shiftE(1:qmmm_struct%nquant_nlink)


        ! QM EWALD SHIFT OR GB SHIFT
        ! --------------------------
        if ( qmmm_nml%qm_ewald>0 .or. qmmm_nml%qmgb == 2) then
#ifdef MPI
          !At present only the master has up to date scf_mchg. Broadcast it
          !to all threads
          call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                         MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif

          if (qmmm_nml%qm_ewald>0) then
!Semi-Parallel
#ifndef SQM
            call qm2_dftb_ewald_shift(scf_mchg)
#endif
          else
            !Must be qmgb==2
!Semi-Parallel
            call qm2_dftb_gb_shift(scf_mchg)
          end if
        end if
  
        ! FINAL HAMILTONIAN
        ! -----------------
        ! Update hamiltonian matrix (here in "a"),
        ! with the charge-dependent part
        do i = 1,qmmm_struct%nquant_nlink
           indi = ks_struct%ind(i)
           liend = lmax( izp_str%izp(i) )**2
           shifti = ks_struct%shift(i)
           do li = 1,liend
              indili = indi + li
               do j = 1,i          
               ! --> Note: Only the lower triangle is used
                 indj = ks_struct%ind(j)
                 ljend = lmax( izp_str%izp(j) )**2
                 shiftj = ks_struct%shift(j)
                 do lj = 1,ljend
                    indjlj = indj + lj
                    ks_struct%a(indili,indjlj) = ks_struct%a(indili,indjlj) &
                          + 0.5*ks_struct%overl(indili,indjlj)*(shifti+shiftj)
                 end do
              end do
           end do
        end do

        call timer_start(TIME_QMMMENERGYSCFDIAG) 
        ! Now that H and S are built (a and b), 
        ! SOLVE the EIGENVALUE PROBLEM
        ! ----------------------------
        ier = 0
        call ewevge(MDIM,MDIM,ndim,ks_struct%a,ks_struct%b,ks_struct%ev, &
                    1,-1,ier)
        call timer_stop(TIME_QMMMENERGYSCFDIAG) 

        ! Convergence failure: 
        if (ier /= 0) then
           write(6,*)
           write(6,*)" QMMM SCC-DFTB: ***************************************************"
           write(6,*)" QMMM SCC-DFTB: ERROR ON EWEVGE (Eigenvalue solver). "
           write(6,*)" QMMM SCC-DFTB: ewevge: ier =",ier,"inner_scf_count=",inner_scf_count
           write(6,*)" QMMM SCC-DFTB: ***************************************************"
           call dftb_conv_failure("eglcao <qm2_dftb_eglcao.f> : ", &
                             "Convergence failure on EWEVGE (Eigenvalue solver).", "Exiting")
        endif
        ! Calculate occupation (occ) and fermi energy (efermi)
        ! using fermi-distribution
        call FERMI(ndim,izp_str%nel,ks_struct%ev,ks_struct%occ,efermi)

        ! Electronic Energy
        ! -----------------
        ! (sum of occupied eigenvalues)       
        elec_eng = 0.0d0
        qm2_struct%nclosed = 0
        qm2_struct%nopenclosed = 0

        do i = 1,ndim
           if (ks_struct%occ(i) < dacc) exit
           elec_eng = elec_eng + ks_struct%occ(i)*ks_struct%ev(i)
           if ( ks_struct%occ(i) == 2.0d0) &
                         qm2_struct%nclosed = qm2_struct%nclosed + 1
           qm2_struct%nopenclosed = qm2_struct%nopenclosed + 1
        end do

        ! Lowest unoccupied level
        lumo = i

        ! MULLIKEN CHARGES
        ! ----------------
        ! qmat will contain the electron populations per atom
        call MULLIKEN(qmmm_struct%nquant_nlink,NDIM,izp_str%izp, &
                      lmax,dacc,qmat,mcharge%qzero,scf_mchg) 

!!           ! OUTPUT EIGENVECTORS
!!           call outeigenvectors(a,ev,occ,ind,qmmm_struct%nquant_nlink)


        ! complete calculation of electronic energy:
        ! charge-dependent energy contribution
        ! warning: this will only lead to the right result if convergence
        ! has been reached
        ecoul = 0.0d0
        eext  = 0.0d0

        ! COULOMBIC ELECTRONIC REPULSION ENERGY
        ! -------------------------------------
        do i = 1, qmmm_struct%nquant_nlink
           ecoul = ecoul + &
                   ks_struct%shift(i) &
                     *(qmat(i)+mcharge%qzero( izp_str%izp(i)))
        end do

        ! EXTERNAL CHARGES (QM/MM) COULOMB ENERGY
        ! ---------------------------------------
        ! This term is a simple coulomb interaction between 
        ! the external charge (in ks_struct%shiftE) and the mulliken charge
        ! on the atom. 
!           do i= qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
        do i= 1, qmmm_struct%nquant_nlink
           eext = eext + ks_struct%shiftE(i)*(mcharge%qzero( izp_str%izp(i))-qmat(i))
        end do

        ! Electronic Energy
        ! -----------------
        ! remark: elec_eng contains ks_struct%shiftE aready via ev,
        ! ks_struct%shift also contains -ks_struct%shiftE, i.e. ecoul also
        ! contains contributions from EXT
        elec_eng = elec_eng-0.5d0*ecoul + 0.5d0*eext

        ! =================
        !  SCC CONVERGENCE
        ! =================

        ! Energy difference
        ediff = elec_eng-elec_eng_old

        ! Maximum charge difference
        chdiff = 0.0d0
        do i = 1, qmmm_struct%nquant_nlink
           chdiff_tmp = abs(ks_struct%qmold(i) - qmat(i))
           if (chdiff_tmp > chdiff) then
              chdiff = chdiff_tmp
              atdiff = i
           end if
        end do

        ! Print convergence progress
        if (qmmm_nml%verbosity > 2) then
           if (inner_scf_count > 1) then
              write (6,'(" QMMM SCC-DFTB: ",i4,3X,3(3X, f20.15),3X,I3,3X,I5,3X,A,3X,F10.5)') &
                    outer_scf_count, elec_eng, ediff, chdiff, atdiff, qmmm_nml%iqmatoms(atdiff),&
                    element_sym(qmmm_struct%iqm_atomic_numbers(atdiff)), scf_mchg(atdiff)
           else
              write (6,'(" QMMM SCC-DFTB: ", 1X,i3,6X, f20.15)') outer_scf_count, elec_eng
           end if

        end if

        ! Convergency Test
        if ( (abs(ediff) < qmmm_nml%scfconv) &
             .and. (chdiff < qmmm_nml%density_conv) ) then
           scc_converged = .true.
#ifdef MPI
           !Tell other threads that we have converged.
           call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, &
                          qmmm_mpi%commqmmm, ier)
#endif

           exit scc_loop
#ifdef MPI
        else
           !Tell other threads that we have not converged.
           call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, &
                          qmmm_mpi%commqmmm, ier)
#endif
        end if

!! ------------------------------
!! Here begins the next iteration
!! ------------------------------

        ! Save the last energy
        elec_eng_old = elec_eng

        ! BROYDEN MIXING
        ! --------------
        call broyden(inner_scf_count,qmmm_struct%nquant_nlink,ks_struct%qmold,qmat)
        qmat(1:qmmm_struct%nquant_nlink) &
            = ks_struct%qmold(1:qmmm_struct%nquant_nlink)

     end do scc_loop ! (end SCC)
   else !This is a slave thread, NOT the master
      scc_slave_loop: do inner_scf_count = 1,qmmm_nml%dftb_maxiter
        outer_scf_count = outer_scf_count + 1
        if ( qmmm_nml%qm_ewald>0 .or. qmmm_nml%qmgb == 2) then
#ifdef MPI
          !At present only the master has up to date scf_mchg. Broadcast it
          !to all threads
          call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                         MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif              
        
          if (qmmm_nml%qm_ewald>0) then
!Semi-Parallel
            call qm2_dftb_ewald_shift(scf_mchg)
          else
            !Must be qmgb==2
!Semi-Parallel
            call qm2_dftb_gb_shift(scf_mchg)
          end if             
        end if
        !Wait to be told by the master if we have converged - implicit barrier in the bcast.
#ifdef MPI
        call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, qmmm_mpi%commqmmm, ier)
#endif
        if (scc_converged) exit scc_slave_loop
      end do scc_slave_loop
   end if

   call timer_stop(TIME_QMMMENERGYSCF)

!! =============================================================================
!!                               End of SCF loop
!! =============================================================================

   ! Calculate the Mulliken charges from the final electron population.
   if (qmmm_mpi%commqmmm_master .and. .not. scc_converged) then
     !We need to calculate scf_mchg again because we called broyden on the way out
     !of scc_loop and this changed qmat.
     do j = 1, qmmm_struct%nquant_nlink
        scf_mchg(j) = mcharge%qzero( izp_str%izp(j) ) - qmat(j)
     end do
   endif 
#ifdef MPI
   !At present only the master has up to date scf_mchg. Broadcast it
   !to all threads
   if (.not. scc_converged .or. qmmm_nml%qm_ewald==0) &
        call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                       MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif

   if (qmmm_mpi%commqmmm_master) then

     ! REPULSIVE ENERGY
     ! ----------------
     ! Repulsive Energy is added after the SCC process is complete.
     !
     call dftb_repulsive(qmmm_struct%nquant_nlink,izp_str%izp,qm_coords,enuclr_qmqm)

     ! EWALD ENERGY CORRECTION
     ! -----------------------
     ! This correction is necessary because, up to this point,
     ! the energy for the Ewald sum only included half of the term from the MM atoms
     ! and half from the QM atoms. The QM atoms should contribute only half but the
     ! MM atoms should contribute full. The corrects qm_ewald_correct_ee routine 
     ! supposedly corrects for this.
     !
     ! That routine uses the density matrix, and expects it to be indexed according
     ! to the contents of qm2_params%orb_loc. 

#ifndef SQM
     if ( qmmm_nml%qm_ewald>0 ) then

        ! Calculate the correction in Hartree/Bohrs

        call qm2_dftb_ewald_corr(qmmm_struct%nquant_nlink, ew_corr, &
                                 qmewald%mmpot, scf_mchg)!qmat)
        ! Put correction into elec_eng.
        elec_eng = elec_eng + ew_corr

     end if

     if ( qmmm_nml%qmgb == 2) then
        !Calculate the correction for GB - removes double counting of
        !GB energy in both escf and egb.
        gb_escf_corr = 0.0d0
        do i = 1,qmmm_struct%nquant_nlink
           gb_escf_corr = gb_escf_corr + (qm_gb%gb_mmpot(i)+qm_gb%gb_qmpot(i))*scf_mchg(i)
        end do
        elec_eng = elec_eng + (0.5d0*gb_escf_corr*BOHRS_TO_A)
    end if
#endif

     ! total energy
     total_e = elec_eng + enuclr_qmqm

     !Calculate Dispersion Energy
     if (qmmm_nml%dftb_disper == 1) then
        call timer_start(TIME_QMMMDFTBDISPE)
        call  dispersion_energy(qmmm_struct%nquant_nlink,qm_coords)
        total_e = total_e + disper%Edis
        call timer_stop(TIME_QMMMDFTBDISPE)
     endif
   end if !(qmmm_mpi%commqmmm_master)

end subroutine eglcao


! Routine to give a message if convergence fails. It has the same syntax as
! sander_bomb.
subroutine dftb_conv_failure(string1,string2,string3)

   use qmmm_module, only:  qmmm_nml, element_sym, qmmm_struct, qm2_struct, qmmm_mpi

!! Passed in:
   ! The only arguments are the strings which are passed to
   ! sander_bomb in case the execution must be stopped.
   character(len=*) :: string1 ! Usually expected to be the name of the subroutine / file
   character(len=*) :: string2 ! Usually expected to be a message
   character(len=*) :: string3 ! Usually expected to be another message.

!! Locals
   logical :: continue_job
   
!! --
   continue_job = .true.
!    continue_job = .false.
   if (qmmm_mpi%commqmmm_master) then
      write(6,'(/," QMMM SCC-DFTB: !!!! ============= WARNING ============= !!!!")')
      write(6,'(  " QMMM SCC-DFTB: Convergence could not be achieved in this step.")')
      if (continue_job) then
         write(6,'(  " QMMM SCC-DFTB: The calculation will continue, but energies and ")')
         write(6,'(  " QMMM SCC-DFTB: forces for this step will not be accurate. ")')
      else
         write(6,'(  " QMMM SCC-DFTB: The calculation will stop.")')
         write(6,'(/," QMMM SCC-DFTB: Last QM Region Cartesian Coordinates ")')
         write(6,'(" QMMM SCC-DFTB: ","SYM",10X,"X",16X,"Y",16X,"Z",13X,"Charge")')
         do i = 1, qmmm_struct%nquant_nlink
            write(6,'(" QMMM Mullik: ",A2,2X,4F16.10)')  &
                  element_sym(qmmm_struct%iqm_atomic_numbers(i)), &
                  (qmmm_struct%qm_coords(j,i), j=1,3), qm2_struct%scf_mchg(i)
         end do
         
         if (qmmm_struct%qm_mm_pairs > 0) then
            write(6,*) ' QMMM SCC-DFTB: number of external charges', qmmm_struct%qm_mm_pairs
            write(6,*) ' QMMM SCC-DFTB: Coordinates of external charges (XYZ)'
            write(6,*)
            write(6,'(" QMMM SCC-DFTB: ",i3)') qmmm_struct%qm_mm_pairs
            do i=1,qmmm_struct%qm_mm_pairs
               write(6,'(" QMMM SCC-DFTB: ",4(2x,f10.6))') (qmmm_struct%qm_xcrd(j,i),j=1,3)
            end do
         end if
         write(6,*) "***************************************************"
         call sander_bomb(string1,string2,string3)
      end if ! if (continue_job) then ...
   end if

   
   
end subroutine dftb_conv_failure


!!=========================================
!! Currently, this subroutine is not used.
!! But I kept it here for the future.
!!=========================================
!subroutine outeigenvectors(a,ev,occ,ind,nn)
!
!   use qm2_dftb_module, only: izp_str
!   use qmmm_module, only: qm2_struct
!
!   implicit none
!   integer, parameter :: nndim = 100
!   integer, parameter :: mdim  = 400
!
!   integer :: ind(nndim+1), nn
!   integer :: i,j,k,l, kl
!   integer :: neig,filled,hfilled
!   _REAL_  :: ks_struct%a(mdim,mdim), ev(mdim), occ(mdim), nel
!   character*4 orbital(9)
!
!   orbital(1) = 'S   '
!   orbital(2) = 'Px  '
!   orbital(3) = 'Py  '
!   orbital(4) = 'Pz  '
!   orbital(5) = 'Dxx '
!   orbital(6) = 'Dxy '
!   orbital(7) = 'Dyy ' 
!   orbital(8) = 'Dyz '
!   orbital(9) = 'Dzz '
!
!   neig = ind(nn+1)
!   filled = 0
!   hfilled = 0
!   do i = 1,neig 
!      if(occ(i) .gt. 1)then
!         filled = filled+1
!      else
!         if(occ(i) .gt. 0)then
!            hfilled = hfilled+1
!         endif
!      endif
!   end do
!   open (1,file='EVC.DAT',status='unknown')
!   rewind 1
!   write (1,'(20x,a)') 'THE ONE-ELECTRON EIGENVALUES AND EIGENVECTORS'
!   write (1,*)
!   write (1,*)
!   write (1,*) 'NUMBER OF FILLED LEVELS:  ',filled
!   if(hfilled .gt. 0)then
!      write (1,*) 'NUMBER OF HALF OCCUPIED STATES: ',hfilled
!   endif
!   write (1,*)
!   write (1,*)
!
!   do i=1,neig
!      kl = 0
!      write(1,13) 'THE ',i,'. EIGENVALUE:',ev(i),' H',ev(i)*27.2116,' eV'
!      write(1,*) 'Atom No.   Atom Type'
!      do j=1,nn
!         write(1,*) '  ',j,'         ',izp_str%izp(j)
!         k=0
!         do l=ind(j)+1, ind(j+1)
!            k=k+1
!            kl=kl+1
!            write(1,'(a,f15.10)') orbital(k),ks_struct%a(l,i)
!            qm2_struct%eigen_vectors(kl,i) = ks_struct%a(l,i)
!         end do
!      end do
!   end do
!   write(1,*)
!   endfile 1
!   close (1)
!13 format(1x,a,i4,a,f15.8,a,4x,f12.5,a)
!end subroutine outeigenvectors



