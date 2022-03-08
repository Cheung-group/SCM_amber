#include "dprec.h"
#include "assert.h"
subroutine putforces(nt,qz,a,f,chg,e)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                          C
!C  PUPIL program ITR-MEDIUM PROJECT                                        C
!C                                                                          C
!C  Author:           J. Torras                                             C
!C  Version:          0.0.8.9                                               C
!C  Date:             09-08-2005                                            C
!C  Proposal Date:    08-08-2005                                            C
!C  Revising Date:    02-06-2006                                            C
!C                                                                          C
!C  torras@qtp.ufl.edu                                                      C
!C                                                                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                          C
!C   FILE WITH ALL THE SUBROUTINES THAT MUST BE IMPLEMENTED DEPENDING OF    C
!C   WHICH ROLE WILL BE PLAY THE BINARY THAT WILL BE LINKED WITH THE        C
!C   PUPIL LIBRARY   libPUPIL.so    (v0.0.8.4)                              C
!C                                                                          C
!C   MOLECULAR DYNAMICS    --> IMPLEMENT   PUTFORCES                        C
!C   DOMAIN IDENTIFICATION --> IMPLEMENT   PUTDATA                          C
!C   QUANTUM CALCULATION   --> IMPLEMENT   PUTCOORDINATES                   C
!C                                                                          C
!C                                                                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   ATTACHED SUBROUTINE THAT INCORPORATES                                  C
!C   EXTERNAL CALCULATION OF FORCE                                          C
!C                                                                          C
!C   NT  TOTAL NUMBER OF ATOMS ASSOC. WITH EXTERNAL CALCULATION OF FORCES   C
!C   QZ  IF QZ IS NON-ZERO, THE QUANTUM ZONE HAS BEEN CHANGED.              C
!C       THE NEW DEFINITION OF THE QUANTUM ZONE IS GIVEN IN AN...           C
!C   A   ...INTEGER VECTOR OF DIMENSION 2*NT WITH THE ATOM NUMBER AND ATOM  C
!C       ROLE                                                               C
!C       ATOM_NUM1,ATOM_ROLE1, ... , ATOM_NUMn,ATOM_ROLEn                   C
!C            Atom Role:   0 -> Quantum atom                                C
!C                         1 -> Pseudoatom                                  C
!C                         2 -> Link Atom                                   C
!C                         3 -> Point charge                                C
!C   F   DOUBLE PRECISION VECTOR WITH CARTESIAN FORCE COMPONENTS ORDERED    C
!C       THE SAME AS PARTICLES IN THE ATOM NUMBER VECTOR IN                 C
!C       eV*ANGSTROM^(-1)                                                   C
!C   CHG DOUBLE PRECISION VECTOR WITH THE CHARGE OF EACH QUANTUM PARTICLE   C
!C       IN UNITS OF ELECTRON CHARGE MAGNITUDE                              C
!C       CHG1, ... , CHGn                                                   C
!C   E   DOUBLE PRECISION VALUE WITH THE TOTAL ENERGY IN eV                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   use pupildata
   use constants  , only : EV_TO_KCAL,AMBER_ELECTROSTATIC,zero
   use qmmm_module, only : qmmm_nml,qmmm_struct,qm2_struct,qmewald,&
                                qmmm_scratch
                                  
#include "memory.h"
#include "debug.h"

   INTEGER   a  (nt*2),nt,qz
   _REAL_    f  (nt*3),e
   _REAL_    chg(nt)

   integer   i,k,totQM
   integer :: ier=0

!  Initial considerations under QZ Change......
   pupQZchange = qz

!  If we are debugging, assume the QM zone has changed.
   if ( do_debugf /= 0 ) pupQZchange = 1

   if (pupQZchange .ne.0) then ! QM-zone changed
      if (qmmm_nml%verbosity > 2) &
         write(6,"(a16,1x,i6,1x,a23)") 'Searching within', nt, 'particles for QM atoms.'
!     Restoring charges to the previous set of quantum particles...
      do i=1,pupqatoms
      !   write(6,"(a12,1x,i6,3x,a12,1x,f12.6,3x,a12,1x,f12.6)") &
      !  'Atom number:',pupqlist(i),'Old charge =',realStack(L15+pupqlist(i)-1),&
      !  'New charge =',pupchg(pupqlist(i))
         realStack(L15+pupqlist(i)-1) = pupchg(pupqlist(i))
      enddo

!         Initialize the pupil mask
      do i=1,natom 
         pupmask(i) = 0
      enddo

!         Counting and getting real quantum particles
      totQM = 0
      do i=1,nt
         pupqlist(i) = abs(a((i-1)*2+1))
         if(a((i-1)*2+2).eq.0) then
            pupmask(a((i-1)*2+1)) = 1
            if (qmmm_nml%verbosity > 2) &
               write(6,"(a7,i6,a1,1x,i6)") 'QM atom', a((i-1)*2+1), ':', i
            totQM = totQM + 1
         endif
      enddo

!         Deleting previous qmmm memory structures
      if(pupStep .gt. 1) then
        ier = 0
        deallocate (qmmm_scratch%qm_real_scratch, stat=ier)
        REQUIRE(ier == 0)
        deallocate ( qmmm_struct%qm_resp_charges,stat=ier)
        REQUIRE(ier == 0)
        deallocate ( qmmm_struct%iqm_atomic_numbers, stat=ier )
        REQUIRE(ier == 0)
        deallocate ( qmmm_nml%iqmatoms, stat=ier )
        REQUIRE(ier == 0)
        deallocate ( qmmm_struct%atom_mask, stat=ier )
        REQUIRE(ier == 0)
      endif

!     Initializing qmm_nml and qmmm_struct structures if these are present. 

      igb = 0
      ier = 0
      qmmm_struct%nquant       = totQM   ! Current number of QM atoms
      qmmm_struct%nlink        = 0       ! Number of link atoms (default: 0)
      qmmm_nml%qm_pme          = .false. ! Use PME to compute long-range QM-QM electrostatics? (default: No)
      qm2_struct%calc_mchg_scf = .false. ! Compute Mulliken charges? (default: No)
      qmmm_nml%qmshake         = 0       ! Use SHAKE? (default: No)
      qmmm_nml%qm_ewald        = 0       ! Use Ewald interpolation for long-range QM-QM electrostatics? (default: No)
      allocate ( qmmm_scratch%qm_real_scratch(3*natom),stat=ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_struct%iqm_atomic_numbers(totQM),stat=ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_struct%qm_resp_charges(totQM),stat=ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_nml%iqmatoms(totQM),stat=ier)
      REQUIRE(ier == 0)
      allocate ( qmmm_struct%atom_mask(natom),stat=ier)
      REQUIRE(ier == 0)
      qmmm_struct%atom_mask(1:natom)       = .false.
      qmmm_struct%qm_resp_charges(1:totQM) = zero

  endif

! From here on, do whether or not the QM zone has changed

! Determine the number of particles in the PUPIL QM region
  pupqatoms = nt

  totQM = 0
  do i=1,nt
     bs1 = (pupqlist(i)-1)*3
     bs2 = (i-1)*3
     do j=1,3
!      Calculate the quantum force
       qfpup(bs1+j) = f(bs2+j)*EV_TO_KCAL
     enddo
!          write(6,"(a22,2x,i4,3(2x,e16.10),2x,a8,2x,f10.7)") 'Quantum Particle		  ==>F:',pupqlist(i),(f(bs2+k)*EV_TO_KCAL,k=1,3),'Charge =',chg(i)			

    ! Manage the quantum list
    if (qmmm_nml%ifqnt .and. a((i-1)*2+2) .eq. 0) then
      if (pupQZchange .ne. 0) then
        totQM = totQM + 1
        qmmm_struct%iqm_atomic_numbers(totQM) = pupatm(pupqlist(i))
        qmmm_nml%iqmatoms(totQM)              = pupqlist(i) 
        qmmm_struct%atom_mask(pupqlist(i))    = .true.
        ! Putting quantum atom charges to zero
        ! BPR - try using Mulliken charges instead
        realStack(L15+pupqlist(i)-1)          = zero
      endif
      ! This instruction was originally inside the pupQZchange .ne. 0 loop.
      ! But, of course, Mulliken charges can change even if the atoms in the QM region do not.
      !write(6,*) 'Updating charges on QM atoms'
      !write(6,*) '----------------------------'
      !write(6,"(a12,2x,i6,2x,a4,f12.6,2x,a4,f12.6)") 'Charge Atom=',pupqlist(i),'Old=',realStack(L15+pupqlist(i)-1),'QM Charge=',chg(i)*AMBER_ELECTROSTATIC
      !realStack(L15+pupqlist(i)-1) = chg(i) * AMBER_ELECTROSTATIC
    endif
  enddo

  ! Print the total QM (SCF) energy to mdout
  qmEnergy = E * EV_TO_KCAL
      if (qmmm_nml%verbosity > 2) then
         write(6,"(a12,e18.10)") 'QM ENERGY = ', qmEnergy
         write(6,"")
         write(6,"")
      endif
  RETURN
END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
