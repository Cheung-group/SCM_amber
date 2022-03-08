#include "copyright.h"
#include "dprec.h"
#include "assert.h"
#include "ncsu-config.h"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ The Molecular Dynamics/NMR Refinement/Modeling Module of the AMBER
!-----------------------------------------------------------------------
!     --- SANDER ---
!-----------------------------------------------------------------------

subroutine sander()

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : &
      ncsu_on_sander_init => on_sander_init, &
      ncsu_on_sander_exit => on_sander_exit
#endif /* DISABLE_NCSU */

   use lmod_driver
! Antonios added use zero
   use constants, only : INV_AMBER_ELECTROSTATIC,zero
   ! The main qmmm_struct contains all the QMMM variables and arrays
   use qmmm_module, only : qmmm_nml,qmmm_struct, deallocate_qmmm, qmmm_mpi, &
         qm2_struct, qmewald, qm_gb, &
#ifdef MPI
         qmmm_mpi_setup, &
#endif
         read_qmmm_namelist_and_allocate

   use genborn
   use decomp, only : allocate_int_decomp, allocate_real_decomp, &
         deallocate_int_decomp, deallocate_real_decomp, &
#ifdef MPI
         synchronize_dec, build_dec_mask, decmask, &
#endif
         nat, nrs, jgroup
   use fastwt
   use relax_mat
   use nmr, only: nmrrad, impnum
   use ew_recip, only: deallocate_m1m2m3,first_pme
   use parms, only: rk,tk,pk,cn1,cn2,numbnd,numang,nimprp,nptra,asol,bsol, &
         nphb,req,charmm
   use nblist, only:cutoffnb,skinnb,nblist_allocate,nblist_deallocate, &
         nblist_allreal,nblist_allint, num_calls_nblist, first_list_flag
   use stack
   use amoeba_runmd, only : AM_RUNMD_get_coords,AM_RUNMD
   use amoeba_mdin, only : beeman_integrator,iamoeba,am_nbead
   use amoeba_interface, only:  &
         AMOEBA_deallocate,AMOEBA_readparm

#ifdef APBS
   use apbs
#endif /* APBS */

   use molecule

#ifdef MPI /* SOFT CORE */
   use softcore, only: setup_sc, cleanup_sc, ifsc, extra_atoms, sc_sync_x, &
        summarize_ti_changes, sc_check_perturbed_molecules
#endif

#if defined(MPI)
   use evb_parm, only: xch_type
#if defined(LES)
   use evb_pimd, only: evb_pimd_init, PE_slice, master_worldrank, jobs_per_node 
#endif
! REMD
   use remd, only : rem, mdloop, repnum, &
                    remd_setup, remd_exchange, remd_cleanup
#else
#  define rem 0
#endif /* MPI */

   use pimd_vars, only: ipimd,ineb
!!jtc ========================= PUPIL INTERFACE =========================
#ifdef PUPIL_SUPPORT
   !  Using data structure for PUPIL   
   use pupildata, x=>realStack , ix=>ixStack , ih=>ihStack
#endif /*PUPIL_SUPPORT*/
!!jtc ========================= PUPIL INTERFACE =========================

   implicit none

   logical belly, erstop
   integer ier,ifind,jn,ncalls,xmin_iter

   character(len=4) itest
   logical ok
   logical newstyle
#  include "files.h"
#  include "memory.h"
#  include "nmr.h"
#  include "box.h"
#  include "md.h"
#  include "extra.h"
#  include "tgtmd.h"
#  include "les.h"
#  include "sgld.h"

#  include "parallel.h"
!Antonios added
#  include "HB.h"
#  include "CHI.h"
integer k
!Antonios end
! Qian added
#  include "debye.h"
! Qian add end
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
#  include "mpif.h"
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#  ifdef MPI_BUFFER_SIZE
   integer*4 mpibuf(mpi_buffer_size)
#  endif
!  REMD: loop is the current exchange. runmd is called numexchg times. 
   integer loop

   integer nrank, istat
   _REAL_ ener(30),vir(4)
   integer ierr
   integer partner
   !     ========================= END AMBER/MPI =========================
#endif
#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "ew_erfc_spline.h"
#  include "ew_parallel.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "def_time.h"

   _REAL_ ene(51)
   integer native,nr3,nr

   ! nmrcal vars
   _REAL_ f,enmr,devdis,devang,devtor,devplpt,devpln,devgendis,ag,bg,cg
   ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints
   integer numphi,nttyp,nhb

   ! runmin/trajene var
   _REAL_ carrms

   ! dipole momemt stuff
   integer ngrp

   character(len=8) initial_date, setup_end_date, final_date
   character(len=10) initial_time, setup_end_time, final_time
   _REAL_ time0, time1

   integer idiff,i,j,istop,index,ierror,itemp

   !jtc ========================= PUPIL INTERFACE =========================
#ifndef PUPIL_SUPPORT
   !   To avoid conflict with the shared data in PUPIL data module
   _REAL_,  dimension(:), allocatable :: x
   integer, dimension(:), allocatable :: ix
   character(len=4), dimension(:), allocatable :: ih
#endif /*PUPIL_SUPPORT*/
   !jtc ========================= PUPIL INTERFACE =========================

   integer, dimension(:), allocatable :: ipairs
   logical qsetup
   logical :: do_list_update=.false.


   !     ---- HERE BEGIN THE EXECUTABLE STATEMENTS ----

   ! Initialize the cpu timer. Needed for machines where returned cpu times
   ! are relative.

   call date_and_time( initial_date, initial_time )
   call wallclock( time0 )
   call init_timers()

   !jtc ========================= PUPIL INTERFACE =========================
#ifdef PUPIL_SUPPORT
   !     captioning line commands and rise up the corba interface

   puperror = 0
   call fixport()
   call inicorbaintfcmd(puperror)
   if (puperror .ne. 0) then
      write(6,*) 'Error in the PUPIL interface initialization.'
      call mexit(6,1)
   endif
   pupactive = .true.
   write(6,*) 'PUPIL CORBA Interface initialized.'
#endif
   !jtc ========================= PUPIL INTERFACE =========================


   ! ==== Flag to tell list builder to print size of list on first call =======
   first_list_flag = .true.
   ! ==== Flag to tell recip space routines to allocate on first call =======
   first_pme = .true.


   ! ==== Initialise first_call flags for QMMM ====
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
   qmmm_struct%num_qmmm_calls = 0

#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     Parallel initialization (setup is now done in multisander).

   !     Make PE 0 the master
   master = mytaskid == 0

   if ( master .and. numtasks > MPI_MAX_PROCESSORS ) then
      write(0, '(a,i4,a,i4)') &
            'Error: the number of processors must not be greater than ', &
            MPI_MAX_PROCESSORS, ', but is ', numtasks
      call mexit(6,1)
   end if
#  ifdef MPI_BUFFER_SIZE
   call mpi_buffer_attach(mpibuf, mpi_buffer_size*4, ierr)
#  endif

   !     ========================= END AMBER/MPI =========================

#else   /* not MPI follows */

   !     in the single-threaded version, the one process is master
   master = .true.
#endif  /* MPI */

   erstop = .false.
   qsetup = .true.

   !     --- generic packing scheme ---

   nwdvar = 1
   native = 32
#ifdef ISTAR2

   !     --- Int*2 packing scheme ---

   nwdvar = 2
#endif  /*ISTAR2*/
   numpk = nwdvar
   nbit = native/numpk

   !     ----- Only the master node (only node when single-process)
   !           performs the initial setup and reading/writing -----

   call timer_start(TIME_TOTAL)
   masterwork: if (master) then

      !        ---- first, initial reads to determine memory sizes:

      call mdread1()
      call amopen(8,parm,'O','F','R')
      call rdparm1(8)
      !        --- now, we can allocate memory:
     
      call locmem()
      write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
      write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
      write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
      write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
      write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr

      !     --- dynamic memory allocation:

      allocate( x(lastr), ix(lasti), ipairs(lastpr), ih(lasth), stat = ier )
      REQUIRE( ier == 0 )
      ix(1:lasti) = 0

      if ((igb /= 0 .and. igb /= 10).or.hybridgb>0) &
         call allocate_gb( natom, ncopy )


      if( idecomp > 0 ) then
#        ifdef MPI
         if (ifsc > 0) then
            call synchronize_dec(natom, nres)
         else
            nat = natom
            nrs = nres
         end if
#        else
         nat = natom
         nrs = nres
#        endif
         call allocate_int_decomp(natom, nres)
      else
         call allocate_int_decomp(1, 1)
      endif

      write(6,'(a,5x,a,i14)') '|', 'nblistReal', nblist_allreal
      write(6,'(a,5x,a,i14)') '|', 'nblist Int', nblist_allint
      write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
            (8*(lastr+lastrst+nblist_allreal)  &
            + 4*(lasth+lasti+lastpr+lastist+nblist_allint))/1024, &
            ' kbytes'

      !        --- finish reading the prmtop file and other user input:
      call rdparm2(x,ix,ih,ipairs,8)

      call AMOEBA_readparm(8,ntf,ntc,natom)! ntf,ntc get reset if amoeba prmtop

      if (qmmm_nml%ifqnt) then
         call read_qmmm_namelist_and_allocate(igb, ih, ix, x, cut, use_pme, ntb)
      end if

      call mdread2(x,ix,ih,ipairs)

      !        --- alloc memory for decomp module that needs info from mdread2
      if( idecomp == 1 .or. idecomp == 2 ) then
         call allocate_real_decomp(nrs)
#        ifdef MPI
         ! -- ti decomp
         if(ifsc > 0) then
            partner = ieor(masterrank,1)
            if (nat == natom) then
               nrank = masterrank
            else
               nrank = partner
            end if
            call mpi_bcast(jgroup, nat, MPI_INTEGER, nrank, commmaster, ierr)
         end if
#        endif
      else if( idecomp == 3 .or. idecomp == 4 ) then
         call allocate_real_decomp(npdec*npdec)
      end if

      !        ----- EVALUATE SOME CONSTANTS FROM MDREAD SETTINGS -----

      nr = nrp
      nr3 = 3*nr
      belly = ibelly > 0

!! jtc ========================= PUPIL INTERFACE =========================
#ifdef PUPIL_SUPPORT
      !     Allocation of memory and initialization
      pupStep  = 0
      puperror = 0
      allocate (qcell   (12     ),stat=puperror)
      allocate (pupmask (natom  ),stat=puperror)
      allocate (pupqlist(natom  ),stat=puperror)
      allocate (pupatm  (natom  ),stat=puperror)
      allocate (pupchg  (natom  ),stat=puperror)
      allocate (qfpup   (natom*3),stat=puperror)
      allocate (qcdata  (natom*9),stat=puperror)
      allocate (keyMM   (natom  ),stat=puperror)
      allocate (pupres  (nres   ),stat=puperror)
      allocate (keyres  (nres   ),stat=puperror)

      if(puperror /= 0) then
         write(6,*) 'Error allocating PUPIL Interface memory.'
         call mexit(6,1)
      endif

      !     Initialization of the Atomic Numbers and quantum forces vector        
      pupqatoms = 0
      iresPup   = 1
      pupres(1) = 1
      do iPup=1,natom
         bs1  = (iPup-1)*3
         call get_atomic_number(ih(iPup+m06-1),pupatm(iPup))
         if(iresPup.lt.nres) then
            if(iPup.ge.ix(iresPup+i02)) then
              iresPup = iresPup + 1
              pupres(iresPup) = iPup
            end if
         endif
         write(strAux,"(A4,'.',A4)") trim(ih(iresPup+m02-1)),adjustl(ih(iPup+m04-1))
         keyres(iresPup) = trim(ih(iresPup+m02-1))
         keyMM(iPup)     = trim(strAux)

         !       Getting all the initial charges....
         pupchg(iPup) = x(L15+iPup-1)
         !   write(6,*) 'Atom num.',iPup,' Label ', keyMM(iPup), 'Charge', pupchg(iPup)

         do jPup=1,3
            qfpup(bs1+jPup) = 0.0d0
         enddo
      enddo
      write(6,*) ' Got all atomic numbers.'

      !     Initialization of the PUPIL cell
      do iPup=1,12
         qcell(iPup) = 0.0d0
      enddo

      !     Submitting the key MM particles and its atomic numbers to PUPIL
      puperror = 0
      call putatomtypes(natom,puperror,pupatm,keyMM)
      if (puperror .ne. 0) then
         write(6,*) 'Error sending MM atom types to PUPIL.'
         call mexit(6,1)
      endif

      !     Submitting the Residue Pointer vector to PUPIL
      write(6,*) 'Number of residues = ',nres,' Number of atoms= ',natom
      !do iPup=1,nres
      !  write(6,*) 'Residue  ',iPup,keyres(iPup),pupres(iPup)
      !enddo
      puperror = 0
      call putresiduetypes(nres,puperror,pupres,keyres)
      if (puperror .ne. 0) then
         write(6,*) 'Error sending MM residue types to PUPIL.'
         call mexit(6,1)
      endif

      write(6,*) 'PUPIL CORBA Interface initialized.'
      write(*,*) ' Initialization of PUPIL structure done.'
#endif
!!jtc ========================= PUPIL INTERFACE =========================


      !        --- seed the random number generator ---

#ifdef MPI
      if(rem <= 0) call amrset(ig)
#else
      call amrset(ig)
#endif

      if (nbit < 32 .and. nr > 32767) then
         write(6, *) '  Too many atoms for 16 bit pairlist -'
         write(6, *) '    Recompile without ISTAR2'
         call mexit(6, 1)
      end if

      if (ntp > 0.and.iabs(ntb) /= 2) then
         write(6,*) 'Input of NTP/NTB inconsistent'
         call mexit(6, 1)
      end if

      !        ----- READ COORDINATES AND VELOCITIES -----

      call timer_start(TIME_RDCRD)
#ifdef LES
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.TRUE.)
#else
      call AMOEBA_check_newstyle_inpcrd(inpcrd,newstyle)
      if ( newstyle )then
         call AM_RUNMD_get_coords(natom,t,irest,ntb,x(lcrd),x(lvel))
      else
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.TRUE.) 
      endif
#endif
     if( iamoeba>0 ) then
        natom = natom*am_nbead
        nrp   = nrp*am_nbead
        nr    = nr*am_nbead
        nr3   = nr3*am_nbead
        ncopy = am_nbead
     end if


      if( igb == 0 .and. induced == 1 ) call get_dips(x,nr)

#ifdef APBS
      ! APBS initialization
      IF ( mdin_apbs ) THEN
         ! in: natom, coords, charge and radii (from prmtop)
         ! out: pb charges and pb radii (via apbs_vars module)
         CALL apbs_init(natom, x(lcrd), x(l15), x(l97))
      END IF
#endif /* APBS */

      !        ----- SET THE INITIAL VELOCITIES -----

      if (ntx <= 3) then
         call setvel(nr,x(lvel),x(lwinv),tempi,init,iscale,scalm)
         ! random numbers may have been "used up" in setting the intial
         ! velocities; re-set the generator so that all nodes are back in
         ! sync

#ifdef MPI
         if(rem <= 0) call amrset(ig)
#else
         call amrset(ig)
#endif

      end if
      if (belly) call bellyf(natom,ix(ibellygp),x(lvel))
      call timer_stop(TIME_RDCRD)

      !        --- If we are reading NMR restraints/weight changes,
      !            read them now:

      if (nmropt >= 1) then
         call nmrcal(x(lcrd),f,ih(m04),ih(m02),ix(i02),x(lwinv),enmr, &
               devdis,devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb,x(lnmr01), &
               ix(inmr02),x(l95),5,6,rk,tk,pk,cn1,cn2, &
               ag,bg,cg,numbnd,numang,numphi,nimprp, &
               nttyp,nhb,natom,natom,ntypes,nres,rad,wel,radhb, &
               welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'READ')
         ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints

         !           --- Determine how many of the torsional parameters
         !               are impropers

         call impnum(ix(i46),ix(i56),ix(i48),ix(i58),nphih,nphia, &
               0,nptra,nimprp)
      end if

      !        --set up info related to weight changes for the non-bonds:

!Antonios commented out for large size to work
!      call nmrrad(rad,wel,cn1,cn2,ntypes,0,0.0d0)
! Antonios end
      call decnvh(asol,bsol,nphb,radhb,welhb)

      if (iredir(4) > 0) call noeread(x,ix,ih)
      if (iredir(8) > 0) call alignread(natom, x(lcrd))
      if (iredir(9) > 0) call csaread(natom, x(lcrd))

      !-----------------------------------------------------------------------
      !        --- Call FASTWAT, which will tag those bonds which are part
      !            of 3-point water molecules. Constraints will be effected
      !            for these waters using a fast analytic routine -- dap.


      call timer_start(TIME_FASTWT )

      call fastwat(ih(m04),nres,ix(i02),ih(m02), &
            nbonh,nbona,ix(iibh),ix(ijbh),ibelly,ix(ibellygp), &
            iwtnm,iowtnm,ihwtnm,jfastw,ix(iifstwt), &
            ix(iifstwr),ibgwat,ienwat,ibgion,ienion,iorwat, &
            6,natom)
      call timer_stop(TIME_FASTWT)

      call getwds(ih(m04)   ,nres      ,ix(i02)   ,ih(m02)   , &
            nbonh     ,nbona     ,0         ,ix(iibh)  ,ix(ijbh)  , &
            iwtnm     ,iowtnm    ,ihwtnm    ,jfastw    ,ix(iicbh) , &
            req       ,x(lwinv)  ,rbtarg    ,ibelly  ,ix(ibellygp), &
            6)

      ! assign link atoms between quantum mechanical and molecular mechanical
      ! atoms if quantum atoms are present
      !
      ! after assigning the link atoms delete all connectivity between the
      ! QM atoms

      if(qmmm_nml%ifqnt) then


         call identify_link_atoms(nbona,ix(iiba),ix(ijba))


         if(nbonh.gt.0) call setbon(nbonh,ix(iibh),ix(ijbh),ix(iicbh), &
               ix(ibellygp)) ! remove bonds between QM atoms from list

         if(nbona.gt.0) call setbon(nbona,ix(iiba),ix(ijba),ix(iicba), &
               ix(ibellygp)) ! remove bonds between QM atoms from list

         if(ntheth.gt.0) call setang(ntheth,ix(i24),ix(i26),ix(i28),ix(i30), &
               ix(ibellygp)) ! remove angles between QM atoms from list

         if(ntheta.gt.0) call setang(ntheta,ix(i32),ix(i34),ix(i36),ix(i38),&
               ix(ibellygp)) ! remove angles between QM atoms from list

         if(nphih.gt.0) call setdih(nphih,ix(i40),ix(i42),ix(i44),ix(i46), &
               ix(i48), ix(ibellygp)) ! remove dihedrals between QM atoms from list

         if(nphia.gt.0) call setdih(nphia,ix(i50),ix(i52),ix(i54),ix(i56), &
               ix(i58), ix(ibellygp)) ! remove dihedrals between QM atoms from list

         ! Now we should work out the type of each quantum atom present. 
         ! This is used for our arrays of pre-computed parameters. It is 
         ! essentially a re-basing of the atomic numbers and is done to save 
         ! memory. Note: qm_assign_atom_types will allocate the qm_atom_type 
         ! array for us. Only the master calls this routine. All other 
         ! threads get this allocated and broadcast to them by the mpi setup 
         ! routine.

         call qm_assign_atom_types

#ifndef LES
         ncopy = 1
#endif

         ! Set default QMMM MPI parameters - for single cpu operation.
         ! These will get overwritten by qmmm_mpi_setup if MPI is on.
         qmmm_mpi%commqmmm_master = master
         qmmm_mpi%numthreads = 1
         qmmm_mpi%mytaskid = 0
         qmmm_mpi%natom_start = 1
         qmmm_mpi%natom_end = natom
         qmmm_mpi%nquant_nlink_start = 1
         qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink
         !Now we know how many link atoms we can allocated the scf_mchg array...
         allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier )
         REQUIRE(ier == 0) !Deallocated in deallocate qmmm
         !We can also allocate ewald_memory
         if (qmmm_nml%qm_ewald > 0 ) then
            call allocate_qmewald(natom)
         end if
         call allocate_qmgb(qmmm_struct%nquant_nlink)

         allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
         REQUIRE(ier == 0) !Deallocated in deallocate qmmm

      endif !if (qmmm_nml%ifqnt)

      ! --- Open the data dumping files and position it depending
      !     on the type of run:

      call open_dump_files
      if (master) call amflsh(6)

      !        --- end of master process setup ---
   end if  masterwork ! (master)

! Pengzhi 9/20/14
! CHANGE                                        - STATUS
! Only master task should read input files      - done
! Commented out writing fo fortran unit 80      - done
! Added bcast of data after inputs are captured - done
! Added timers                                  - done
tacc_masterwork: if(master) then

! Qian added
! Reading debye.inp
  open(81,file='debye.inp')
  read(81,*) reladielec,strength
  close(81)
! Qian add end

! Antonios added 9/2/10
! Reading triple.inp for Chiral term

      open(78,file='triple.inp')
      read(78,*)nbeta,xxkchi

        XKCHI=xxkchi
        ICHI=nbeta
! Pengzhi 9/20/14
!      write(80,*)ICHI,XKCHI

      do j=1,nbeta
        read(78,*)IICA,IICB,IINC,IICC,xxchi
! Pengzhi 9/20/14
!        write(80,*)IICA,IICB,IINC,IICC,xxchi
       ICA(j) = iica
       ICB(j) = iicb
       INC(j) = iinc
       ICC(j) = iicc
       XCHI(j) = xxchi
      enddo
                                      
! Antonios added 2/9/10
! Reading pseudo.inp for Hydrogen Bonds

      open(79,file='pseudo.inp')
      read(79,*)nhb_pair,Amp
! Pengzhi 9/20/14
!      write(80,*)nhb_pair,Amp
      if(natom.gt.10000) then
      write(755,*) 'error:natom gt 10000,check HB.h'
      endif

      if(nhb_pair.gt.92000) then
      write(755,*) 'error:nhb_pair gt 92000,check HB.h'
      endif

      do j=1,natom
       do k=1,natom
       maphb(j,k) = zero
       enddo
      enddo



      do j=1,nhb_pair
       NXI(j) = zero
       NXI1(j) = zero
       NXI2(j) = zero
       NXI3(j) = zero
      enddo



! ixi1<ixi2
      do j=1,nhb_pair
        read(79,*)ixi,ixi1,ixi2,ixi3
! Pengzhi 9/20/14
!        write(80,*)ixi,ixi1,ixi2,ixi3

       maphb(ixi1,ixi2) = j
       maphb(ixi2,ixi1) = j

!        write(80,*)maphb(j,j)

       NXI(j) = ixi
       NXI1(j) = ixi1
       NXI2(j) = ixi2
       NXI3(j) = ixi3

      enddo

close(78)
close(79)
close(80)
! Antonios end

! Pengzhi 9/20/14
 end if tacc_masterwork

#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     NOTE: in the current AMBER/MPI implementation, two means of
   !     running in parallel within sander are supported. The value
   !     of mpi_orig determines which approach is used.
   !     This is turned on when minimization (imin .ne. 0) is requested,
   !     and is otherwise off.

   !     When running the mpi_orig case, a variable notdone is now
   !     set by the master and determines when to exit the force()
   !     loop.  When the master has finished calling force, the
   !     master changes notdone to 0 and broadcasts the data one more
   !     time to signal end of the loop.  force() is modified so that
   !     in the mpi_orig case, an initial broadcast is done to receive
   !     the value from the master to decide whether to do the work or
   !     simply exit.

   !     ...set up initial data and send all needed data to other nodes,
   !     now that the master has it

   !     First, broadcast parameters in memory.h, so that all processors
   !     will know how much memory to allocate:

! Broadcast data read from input files
! From debye.inp
 CALL MPI_Bcast( reladielec, 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr ) 
 CALL MPI_Bcast( strength,   1, MPI_DOUBLE_PRECISION, 0, commsander, ierr ) 
! From triple.inp
! Pengzhi 10/30/15
! added broadcasting of ICHI as well 
! nbeta = ICHI, but only ICHI was used later in extra_pts.h
 CALL MPI_Bcast( ICHI,    1, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( nbeta,    1, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( ICA,  nbeta, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( ICB,  nbeta, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( INC,  nbeta, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( ICC,  nbeta, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( XKCHI,    1, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
 CALL MPI_Bcast( XCHI, nbeta, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
! From pseudo.inp
! NOTE: I am skipping hte initialization of NXI, NXI1, NXI2, NXI3 to zero 
!       since they are read from file immediately afterwards
 CALL MPI_Bcast( natom,    1, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( nhb_pair, 1, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( NXI,  nhb_pair, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( NXI1, nhb_pair, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( NXI2, nhb_pair, MPI_INTEGER, 0, commsander, ierr )
 CALL MPI_Bcast( NXI3, nhb_pair, MPI_INTEGER, 0, commsander, ierr )
 DO j=1,natom
   DO k=1,natom
      maphb(j,k) = zero
   END DO
 END DO
 DO j=1,nhb_pair
       ixi1 = NXI1(j)
       ixi2 = NXI2(j)
       maphb(ixi1,ixi2) = j
       maphb(ixi2,ixi1) = j
 END DO
! Pengzhi end 9/20/14
   call mpi_bcast(natom,BC_MEMORY,mpi_integer,0,commsander,ierr)
   ! -- ti decomp
   call mpi_bcast(idecomp,1,mpi_integer,0,commsander,ierr)
   call mpi_bcast(nat,1,mpi_integer,0,commsander,ierr)
   call mpi_bcast(nrs,1,mpi_integer,0,commsander,ierr)

   !----- Set up integer stack initial size --------------    
   call mpi_bcast(lastist,1,mpi_integer,0,commsander,ierr)
   call mpi_bcast(lastrst,1,mpi_integer,0,commsander,ierr)

   call stack_setup()

   call mpi_bcast(charmm,1,MPI_LOGICAL,0,commsander,ierr)
   if(charmm)call mpi_bcast(nimphi,1,MPI_INTEGER,0,commsander,ierr)

   call mpi_barrier(commsander,ierr)
   !     ---allocate memory on the non-master nodes:

   if( .not.master ) then
      allocate( x(1:lastr), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ix(1:lasti), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ipairs(1:lastpr), stat = ier )
      REQUIRE( ier == 0 )

      allocate( ih(1:lasth), stat = ier )
      REQUIRE( ier == 0 )

      if( charmm ) then
         allocate(im(nimphi),jm(nimphi),km(nimphi),lm(nimphi),imp(nimphi), &
            stat=ier)
         REQUIRE( ier == 0 )
      end if
      ! -- ti decomp
      if( idecomp > 0 ) then
         call allocate_int_decomp(natom, nres)
         if( idecomp == 1 .or. idecomp == 2 ) then
            call allocate_real_decomp(nrs)
         else if( idecomp == 3 .or. idecomp == 4 ) then
            call allocate_real_decomp(npdec*npdec)
         end if
      endif
   end if  ! ( .not.master )

   if(idecomp == 1 .or. idecomp == 2) then
      call mpi_bcast(jgroup,nat,MPI_INTEGER,0,commsander,ierr)
   end if

   call startup_groups(ierr)
   call startup(x,ix,ih)

!  +---------------------------------------------------------------------------+
!  |  Broadcast EVB/PIMD inputs/parameters to all PEs                          |
!  +:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::+
!  |  Note: The masters have all required EVB/PIMD inputs/parameters via call  |
!  |  to mdread2 ( evb_input, evb_init, evb_pimd_init ).  For EVB/PIMD, all    |
!  |  PEs need the inputs/parameters ... so we perform this initialization     |
!  |  again for all PEs besides the masters.  The alternative is to use        |
!  |  MPI_BCAST.                                                               |
!  +---------------------------------------------------------------------------+

   call mpi_bcast ( ievb, 1, MPI_INTEGER, 0, commworld, ierr )
   call mpi_bcast ( ipimd, 1, MPI_INTEGER, 0, commworld, ierr )
!KFW
   if( ievb /= 0 ) then
      call mpi_bcast ( evbin, 256, MPI_CHARACTER, 0, commworld, ierr )
      if( .not. master ) then
         call evb_input
         call evb_init
      endif
   endif
!KFW
#if defined(LES)
   call mpi_bcast ( ncopy, 1, MPI_INTEGER, 0, commworld, ierr )
   call mpi_bcast ( cnum(1:natom), natom, MPI_INTEGER, 0, commworld, ierr )
   call mpi_bcast ( evbin, 256, MPI_CHARACTER, 0, commworld, ierr )
#endif /* LES */
   if( ievb /= 0 ) then
#if defined(LES)
!     call mpi_bcast ( mastersize, 1, MPI_INTEGER, 0, commworld, ierr )
!     call mpi_bcast ( jobs_per_node, 1, MPI_INTEGER, 0, commworld, ierr )
!     call mpi_bcast ( nsize, 1, MPI_INTEGER, 0, commworld, ierr )
      if( ipimd > 0 .and. .not. master ) then
         call evb_input
         call evb_init
!        call evb_pimd_init
      endif
         call evb_pimd_init
!     call mpi_bcast ( master_worldrank, mastersize, MPI_INTEGER, 0, commworld, ierr )
!     call mpi_bcast ( PE_slice, jobs_per_node*nsize, MPI_INTEGER, 0, commworld, ierr )
#endif /* LES */
      call evb_bcast
      call evb_alloc
   endif

!  +---------------------------------------------------------------------------+
!  |  Obtain B vector for Schlegel's distributed Gaussian method               |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( xch_type ) ) == "dist_gauss" ) call schlegel_dg

#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      ! multi-CPU minimization do not work with soft core !
      if (imin > 0 .and. numtasks > 1) then
         call sander_bomb('imin > 0 and numtasks > 1','TI minimizations cannot be performed with > 2 CPUs','')
      end if
      call setup_sc(natom, nres, ih(m04), ih(m06), &
          ix(i02), ih(m02), x(lcrd), ntypes, clambda, nstlim)
      if (ntp > 0 .and. master) then 
         ! check which molecules are perturbed in NPT runs 
         call sc_check_perturbed_molecules(nspm, ix(i70))
      end if
      ! -- ti decomp
      if (idecomp > 0) then
         if (sanderrank == 0) call build_dec_mask(natom)
         call mpi_bcast(decmask, natom, MPI_INTEGER, 0, commsander, ierr)
      end if
      ! Make sure all common atoms have the same v (that of V0) in TI runs
      if ( ifsc /= 2 ) then
         if (master) call sc_sync_x(x(lvel),nr3)
         !call mpi_barrier(commsander,ierr)
         if (numtasks > 1) then
            call mpi_bcast(nr3,1,MPI_INTEGER,0,commsander,ierr)
            call mpi_bcast(x(lvel),nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)               
         end if
      end if
      if (master) call setnoshake_sc(ix,ntc,num_noshake)
   else
      extra_atoms=0
   end if
#endif

  if(.not. master)   call nblist_allocate(natom,ntypes,num_direct,numtasks)

   !  ---allocate memory for GB on the non-master nodes:

   if( .not.master ) then
      if ((igb /= 0 .and. igb /= 10).or.hybridgb>0) &
         call allocate_gb( natom, ncopy )
   end if  ! ( .not.master )

   nr = nrp
   nr3 = 3*nr
   belly = ibelly > 0

   !  ---Do setup for QMMM in parallel if ifqnt is on - this is essentially 
   !     everything in the QMMM namelist and a few other bits and pieces.
   !  ---Note, currently only master node knows if qmmm_nml%ifqnt is 
   !     on so we need to broadcast this first and then make decisions based 
   !     on this.

   call mpi_bcast(qmmm_nml%ifqnt, 1, mpi_logical, 0, commsander, ierror)

   if (qmmm_nml%ifqnt) then
      ! Broadcast all of the stuff in qmmm_nml and allocate the relevant 
      ! arrays on all processors. This will also set up information for
      ! openmp on the master processor if it is in use.
      call qmmm_mpi_setup( master )
      if (qmmm_nml%qm_ewald > 0 .and. .not. master) then
         call allocate_qmewald(natom)
      end if
      if (qmmm_nml%qmgb==2 .and. .not. master) then
         call allocate_qmgb(qmmm_struct%nquant_nlink)
      end if

      allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
      REQUIRE(ier == 0) !Deallocated in deallocate qmmm

   end if

   !    --- END QMMM MPI SETUP ---

   ! REM: call amrset if not REM. there is no need to do this everytime
   ! DAN ROE: would only be done once now anyway
   if(rem <= 0) call amrset(ig+1) 

   if (nmropt >= 1) &
         call nmrcal(x(lcrd),f,ih(m04),ih(m02),ix(i02),x(lwinv),enmr, &
         devdis,devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb,x(lnmr01), &
         ix(inmr02),x(l95),5,6,rk,tk,pk,cn1,cn2, &
         ag,bg,cg,numbnd,numang,numphi,nimprp, &
         nttyp,nhb,natom,natom,ntypes,nres,rad,wel,radhb, &
         welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'MPI ')
     ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints

   !   ---------------- Check system is neutral and print warning message ------
   !   ---------------- adjust charges for roundoff error.                ------
   if( igb == 0 .and. iyammp == 0 ) call check_neutral(x(l15),natom)

   !    ---------------- Old parallel for minimization ----------------------

   if (imin /= 0) then
      mpi_orig = .true.
      notdone = 1
   else
      mpi_orig = .false.
   end if

   if (mpi_orig .and. .not.master) then

      ! ...all nodes only do the force calculations (JV)
      ! Minimisation so only only master gets past the loop below
      ! hence need to zero QM charges on non-master threads here.

      if (qmmm_nml%ifqnt) then
         !Apply charge correction if required.
         if (qmmm_nml%adjust_q>0) then
            call qmmm_adjust_q(qmmm_nml%adjust_q, natom, qmmm_struct%nquant, qmmm_struct%nquant_nlink, &
                  qmmm_struct%nlink, x(L15), &
                  qmmm_nml%iqmatoms, qmmm_nml%qmcharge, qmmm_struct%atom_mask, &
                  qmmm_struct%mm_link_mask, master,x(LCRD))
         end if
         ! zero out the charges on the quantum mechanical atoms
         call qm_zero_charges(x(L15))
         if (qmmm_struct%nlink > 0 ) then
            !We need to exclude all electrostatic
            !interactions with MM link pairs, both QM-MM and MM-MM. Do this by
            !zeroing the MM link pair charges in the main charge array.
            !These charges are stored in qmmm_struct%mm_link_pair_resp_charges in case
            !they are later needed.
            call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink,qmmm_struct%link_pairs,x(L15))
         end if

         ! At this point we can also fill the qmmm_struct%scaled_mm_charges
         ! array - we only need to do this once as the charges are constant
         ! during a run. Having a separate array of scaled charges saves us
         ! having to do it on every qmmm routine call.
         do i = 1, natom
            qmmm_struct%scaled_mm_charges(i) = x(L15+(i-1)) * INV_AMBER_ELECTROSTATIC &
                                                            * qmmm_nml%chg_lambda ! charge scaling factor for FEP
         end do
      end if

      if (igb == 7 ) call igb7_init(natom, ncopy, x(l97)) !x(l97) is rborn()

      do while( notdone == 1 )
         call force(x,ix,ih,ipairs,x(lcrd),x(lforce),ener,vir, &
               x(l96), x(l97), x(l98), x(l99), qsetup, &
               do_list_update)
      end do

      goto 999  ! deallocate and return
   end if
   !    ----------------------------------------------------------------------

   if (master) then
      write(6, '(a,i4,a,/)') '|  Running AMBER/MPI version on ',numtasks, ' nodes'
      !BTREE is selected by default if cpu is a power of two.
      !The number of processes is required to be a power of two for Btree
      !Print a warning about inefficiency with cpus not being a power of 2.
      if ( numtasks>1 .and. logtwo(numtasks) <= 0 ) then 
         write(6, '(a,i4,a,/)') '|  WARNING: The number of processors is not a power of 2'
         write(6, '(a,i4,a,/)') '|           this may be inefficient on some systems.'
      end if
   end if
   if (master .and. numgroup > 1) write(6, '(a,i4,a,i4,a,i4,a)') &
         '|  MULTISANDER: ', numgroup, ' groups. ', &
         numtasks, ' processors out of ', worldsize, ' total.'
   if(master)call amflsh(6)

   !     ========================= END AMBER/MPI =========================
#else
   !   debug needs to copy charges at start and they can't change later
   !   ---------------- Check system is neutral and print warning message ------
   !   ---------------- adjust charges for roundoff error.                ------
   if( igb == 0 .and. iyammp == 0 ) call check_neutral(x(l15),natom)

   call amrset(ig+1)
   call stack_setup()

#endif /* MPI */

#ifdef QMMM_OMP
      !If we are using openmp for matrix diagonalization print some information.
      if (qmmm_nml%ifqnt .and. master) call qm_print_omp_info()
#endif

   call date_and_time( setup_end_date, setup_end_time )

   ! ----------------------------------------------------------------------
   ! Now do the dynamics or minimization.
   ! ----------------------------------------------------------------------

   if (igb == 7 ) call igb7_init(natom, ncopy, x(l97)) !x(l97) is rborn()

   if (qmmm_nml%ifqnt) then
      !Apply charge correction if required.
      if (qmmm_nml%adjust_q>0) then
         call qmmm_adjust_q(qmmm_nml%adjust_q, natom, qmmm_struct%nquant, qmmm_struct%nquant_nlink, &
               qmmm_struct%nlink, x(L15), &
               qmmm_nml%iqmatoms, qmmm_nml%qmcharge, qmmm_struct%atom_mask, &
               qmmm_struct%mm_link_mask, master,x(LCRD))
      end if
      !Zeroing of QM charges MUST be done AFTER call to check_neutral.
      ! zero out the charges on the quantum mechanical atoms
      call qm_zero_charges(x(L15))
      if (qmmm_struct%nlink > 0) then
         !We need to exclude all electrostatic
         !interactions with MM link pairs, both QM-MM and MM-MM. Do this by
         !zeroing the MM link pair charges in the main charge array.
         !These charges are stored in qmmm_struct%mm_link_pair_resp_charges in case
         !they are later needed.
         call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink,qmmm_struct%link_pairs,x(L15))
      end if

      ! At this point we can also fill the qmmm_struct%scaled_mm_charges
      ! array - we only need to do this once as the charges are constant
      ! during a run. Having a separate array of scaled charges saves us
      ! having to do it on every qmmm routine call.
      do i = 1, natom
         qmmm_struct%scaled_mm_charges(i) = x(L15+(i-1)) * INV_AMBER_ELECTROSTATIC &
                                                         * qmmm_nml%chg_lambda ! charge scaling factor for FEP
      end do
   end if

   ! use the debugf namelist to activate
   call debug_frc(x,ix,ih,ipairs,x(lcrd),x(lforce), &
         cn1,cn2,qsetup)

   if( master .AND. (.NOT. qmmm_nml%ifqnt)) &
         write(6,'(/80(''-'')/,''   4.  RESULTS'',/80(''-'')/)')

   ! Input flag imin determines the type of calculation: MD, minimization, ...

   select case ( imin )
   case ( 0 )
      !        --- Dynamics:

      !  Prepare for Isotropic periodic sum of nonbonded interaction
      if( ips.gt.0 ) call ipssys(natom,ntypes,ntb,x(l15),  &
            cut,cn1,cn2,ix(i04),x(lcrd))
      !  Prepare for SGLD simulation
      if(tsgld)call psgld(natom,x(lmass),x(lvel),x(lvsg))

      call timer_start(TIME_RUNMD)

#ifdef MPI
      ! ----===== REMD =====----
      ! If this is not a REMD run, runmd is called only once. 
      ! If this is a REMD run, runmd is called 0 to numexchg times,
      !  where the 0th runmd is just for getting initial PEs (no dynamics).
      if (rem==0) then
         ! Not a REMD run. runmd will be called once.
         loop=0
      else
         ! This is a REMD run. runmd will be called numexchg times.
         loop=numexchg

         ! Set up temptable, open remlog, etc.
         call remd_setup(numexchg,hybridgb,numwatkeep,temp0,mxvar,natom)
      endif ! Replica run setup

      ! Loop over REMD exchanges
      do mdloop=0, loop

         ! ----===== REMD EXCHANGE HANDLING =====----
         ! Note: mdloop==0 is just used to get initial energies for the
         !        first exchange. 
         if (rem>0 .and. mdloop>0) then
            call remd_exchange(x(lcrd),x(lvel),nr3,natom,nr,temp0)
         endif ! rem>0 and mdloop>0 
         ! ----===== END REMD EXCHANGE HANDLING =====----

         if (rem>0.and.mdloop.eq.0.and.master) &
            write (6,'(a,i4)') "REMD: Getting initial energy for replica ",repnum
#endif  /* MPI */

      if ( ipimd > 0 ) then
         call pimd_init(natom,x(lmass),x(lwinv),x(lvel),ipimd)
      end if

      if(ineb>0) call neb_init()

#     ifndef DISABLE_NCSU
      call ncsu_on_sander_init(ih, x(lmass), x(lcrd), rem)
#     endif /* DISABLE_NCSU */

      if ( beeman_integrator == 1 )then
         call AM_RUNMD( ix,ih,ipairs, &
               x(lwinv),x(lmass),x, &
               x(lcrd),x(lvel),x(lforce),qsetup)
      else

         call runmd (x,ix,ih,ipairs, &
               x(lcrd),x(lwinv),x(lmass),x(lforce), &
               x(lvel),x(lvel2),x(l45),x(lcrdr), &
               x(l50),x(l95),ix(i70),x(l75), &
               erstop,qsetup)

      endif !beeman_integrator == 1

#     ifndef DISABLE_NCSU
      call ncsu_on_sander_exit()
#     endif /* DISABLE_NCSU */

#ifdef MPI
      enddo ! Loop over REMD exchanges (mdloop)

      ! Cleanup REMD files.
      if (rem>0) call remd_cleanup()

      ! ----===== END REMD =====----
#endif


      call timer_stop(TIME_RUNMD)

      if (master) call amflsh(6)

      if (erstop) then
         ! This error condition stems from subroutine shake;
         ! furthermore, it seems that erstop can never be true since shake
         ! can never return with its third last argument, niter, equal to 0.
         ! SRB, Sep 24, 2003
         if (master) then
            write(6, *) 'FATAL ERROR'
         end if
         call mexit(6,1)
      end if

   case ( 1 )

      !        --- Minimization:

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
               ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(ibellygp), &
               x(l95),ene, carrms, qsetup)
      case ( LMOD_NTMIN_XMIN )
         write(6,'(a,i4)') '  LMOD XMIN Minimization.'
         xmin_iter = 0
         call run_xmin( x, ix, ih, ipairs, &
               x(lcrd), x(lforce), ene, qsetup, xmin_iter )
      case ( LMOD_NTMIN_LMOD )
         write(6,'(a,i4)') '  LMOD LMOD Minimization.'
         call run_lmod( x, ix, ih, ipairs, &
               x(lcrd), x(lforce), ene, qsetup )
      case default
         ! invalid ntmin
         ! ntmin input validation occurs in mdread.f
         ASSERT( .false. )
      end select

      if (master) call minrit(x(lcrd))  ! Write the restart file

   case ( 5 )
      !       ---carlos modified for reading trajectories (trajene option)

      write (6,*) "POST-PROCESSING OF TRAJECTORY ENERGIES"

      !       ---read trajectories and calculate energies for each frame

      call trajene(x,ix,ih,ipairs, ene,ok,qsetup)

      if (.not.ok) then
         write (6,*) 'error in trajene()'
         call mexit(6,1)
      end if

   case default
      ! invalid imin
      ! imin input validation should be transferred to mdread.f
      write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (',imin,').'
      ASSERT( .false. )
   end select

#ifdef MPI /* SOFT CORE */
   if (master) then
      if (icfe /=0 .and. ifsc == 1) call summarize_ti_changes(natom,resat)
   end if
#endif

   !     -- calc time spent running vs setup


   call timer_stop(TIME_TOTAL)
   call wallclock( time1 )
   call date_and_time( final_date, final_time )
   call profile_time( time1 - time0 , num_calls_nblist,profile_mpi)

#ifdef MPI
   !     =========================== AMBER/MPI ===========================

   !     Set and broadcast notdone in mpi_orig case to inform
   !     other nodes that we are finished calling force(). (tec3)

   if ( mpi_orig ) then
      notdone = 0
      call mpi_bcast(notdone,1,mpi_integer,0, commsander,ierr)
   end if

   !     ========================= END AMBER/MPI =========================
#endif

!! jtc ========================= PUPIL INTERFACE =========================
#ifdef PUPIL_SUPPORT
!     Finalize Corba Interface
      puperror = 0
      call killcorbaintfc(puperror)
      if(puperror /= 0) then
         write(6,*) 'Error ending PUPIL CORBA Interface.'
      endif
      pupactive = .false.
      write(6, '(a)') 'PUPIL CORBA Interface finalized.'
#endif
!! jtc ========================= PUPIL INTERFACE =========================

   

   if( master ) then
      call close_dump_files

      !     --- write out final times

      write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
            ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
            initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
      write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
            ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
            setup_end_date(5:6), '/',setup_end_date(7:8),'/',setup_end_date(1:4)
      write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
            ':', final_time(3:4), ':', final_time(5:10), '  on ', &
            final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
      call nwallclock( ncalls )
      write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
      call amflsh(6)

      if(iesp > 0) then
         call esp(natom,x(lcrd),x(linddip))
      end if
   end if
999 continue  !     --- dynamic memory deallocation:

   if (qmmm_nml%ifqnt .and. .not. qmmm_struct%qm_mm_first_call) then   
      call deallocate_qmmm      !If first_call is still true then this thread
      !never actually called the QMMM routine. E.g. more
      !threads than PIMD  replicates
   end if

   if(ipimd>0) call pimd_finalize(ipimd)

   if(ineb>0) call neb_finalize()

   if( idecomp > 0 ) then
      call deallocate_real_decomp()
      call deallocate_int_decomp()
   endif
   if(master .and. idecomp == 0) call deallocate_int_decomp()

#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      call cleanup_sc()
   end if
#endif

   call nblist_deallocate()
   call deallocate_stacks()
   if (( igb /= 0 .and. igb /= 10 ).or.hybridgb>0) then
      call deallocate_gb( )
   end if
   if (master) then  
      if( igb == 10 ) then
         call pb_free( )
      end if
   end if
   deallocate( ih, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( ipairs, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( ix, stat = ier )
   REQUIRE( ier == 0 )
   deallocate( x, stat = ier )
   REQUIRE( ier == 0 )
   if( ntb>0 .and. ifbox==1 .and. ew_type==0 .and. mpoltype==0 ) &
         call deallocate_m1m2m3()
   call AMOEBA_deallocate

   if (master.and.mdout /= 'stdout') close(6)

   return

end subroutine sander

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine esp here]
subroutine esp(natom,x,mom_ind)

   ! routine to calculate the ESP due to the induced moments (only)
   ! at the same spatial points as the reference QM.
   use constants, only : zero, BOHRS_TO_A, INV_AMBER_ELECTROSTATIC
   implicit none
   integer natom
   _REAL_  x(3,*)
   _REAL_  mom_ind(3,*)

#  include "files.h"
#  include "ew_mpole.h"

   integer dat_unit, new_unit, minus_new_unit
   parameter(dat_unit=30, new_unit=31, minus_new_unit=33)


   integer inat, nesp, idum
   _REAL_  xin, yin, zin
   integer jn, kn
   _REAL_  esp_qm, xb_esp, yb_esp, zb_esp
   _REAL_  x_esp, y_esp, z_esp
   _REAL_  e_x, e_y, e_z, e_q, esp_new
   _REAL_  dist, dist3
   integer iptr

   call amopen(dat_unit,"esp.dat",'O','F','R')
   call amopen(new_unit,"esp.induced",owrite,'F','W')
   call amopen(minus_new_unit,"esp.qm-induced",owrite,'F','W')
   read (dat_unit,'(3i5)')inat,nesp,idum
   write(6,'(t2,''inat = '',i5)')inat
   write(6,'(t2,''nesp = '',i5)')nesp

   write(new_unit,'(2i5)')inat,nesp
   write(minus_new_unit,'(2i5)')inat,nesp

   if (inat /= natom) then
      write(6,'(t2,''natom mismatch with esp file'')')
      call mexit(6,1)
   end if
   do jn = 1,inat
      read (dat_unit,'(17x,3e16.0)')xin,yin,zin
      write(new_unit,'(17x,3e16.7)')xin,yin,zin
      write(minus_new_unit,'(17x,3e16.7)')xin,yin,zin
   end do
   do jn = 1,nesp
      e_x = zero
      e_y = zero
      e_z = zero
      e_q = zero
      read(dat_unit,'(1x,4e16.0)')esp_qm,xb_esp,yb_esp,zb_esp
      x_esp = xb_esp * BOHRS_TO_A
      y_esp = yb_esp * BOHRS_TO_A
      z_esp = zb_esp * BOHRS_TO_A
      do kn = 1,natom
         dist = (sqrt((x(1,kn)-x_esp)**2 + &
               (x(2,kn)-y_esp)**2 + &
               (x(3,kn)-z_esp)**2))
         dist3 = dist**3
         e_x = e_x - mom_ind(1,kn   )*(x(1,kn)-x_esp)/dist3
         e_y = e_y - mom_ind(2,kn   )*(x(2,kn)-y_esp)/dist3
         e_z = e_z - mom_ind(3,kn   )*(x(3,kn)-z_esp)/dist3
      end do
      e_x = e_x * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
      e_y = e_y * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
      e_z = e_z * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
      e_q = e_q * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
      esp_new = e_x + e_y + e_z

      write(new_unit,      '(1x,4e16.7)')esp_new &
            ,xb_esp,yb_esp,zb_esp
      write(minus_new_unit,'(1x,4e16.7)')esp_qm-esp_new &
            ,xb_esp,yb_esp,zb_esp
   end do

   close(dat_unit)
   close(new_unit)
   close(minus_new_unit)
   return
end subroutine esp

