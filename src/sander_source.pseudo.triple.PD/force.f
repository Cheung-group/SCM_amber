! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "assert.h"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine to compute energies and forces
subroutine force(xx,ix,ih,ipairs,x,f,ener,vir, &
      fs, rborn, reff, onereff, qsetup,do_list_update )

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only : ncsu_on_force => on_force
#  if !defined(LES) && defined(MPI)
   use remd, only : rem, mdloop
#  endif
#  ifdef MPI
   use ncsu_sander_proxy, only : ncsu_remember_initremd => remember_initremd
#  endif
#endif

   use genborn
   use poisson_boltzmann, only : pb_force
   use dispersion_cavity, only : npopt, np_force
#ifdef APBS
   use apbs
#endif /* APBS */
   use trace
   use stack
   use pimd_vars, only: ipimd,ineb,nbead,bnd_vir,Epot_spring,Epot_deriv,real_mass,equal_part,nrg_all,nebrms,itimass
   use full_pimd_vars, only: totener
   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, &
                           qmmm_mpi, qmewald
   use constants, only : zero
   use relax_mat
   use ew_recip
   use parms, only:cn1,cn2,asol,bsol,pk,rk,tk,numbnd,numang,nptra,nphb,nimprp
   use parms, only: cn114,cn214
#ifdef PUPIL_SUPPORT

   ! Using the ucell variable in nblist
   use nblist, only:nonbond_list,a,b,c,alpha,beta,gamma,ucell
#else
   use nblist, only:nonbond_list, &
                    a,b,c,alpha,beta,gamma !only used in qmmm call
#endif /*PUPIL_SUPPORT*/

#ifdef DSSP
   use dssp, only: fdssp, edssp, idssp
#endif


   use amoeba_interface,only: AM_VAL_eval,AM_NonBond_eval
   use amoeba_mdin, only : iamoeba,am_nbead


#if defined(MPI)
   use evb_data, only: nrg_frc
   use softcore, only: sc_ener
#  ifdef LES
      use miller, only: dlnQ_dl
      use remd, only : rem, mdloop
#  endif
#endif /* MPI */

#ifdef PUPIL_SUPPORT
   use pupildata
#endif /*PUPIL_SUPPORT*/

   use cns_xref


   implicit none

#ifdef PUPIL_SUPPORT
   character(kind=1,len=5) :: routine="force"
#endif
   integer   ipairs(*)
   _REAL_ xx(*)
   integer   ix(*)
   character(len=4) ih(*)
   _REAL_ fs(*),rborn(*),reff(*),dvdl
   _REAL_, intent(out) :: onereff(*)
#include "def_time.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "extra_pts.h"
#include "parallel.h"
#include "HB.h"

! Qian
#include "CROWD.h"
! end

#ifdef MPI
#include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
#include "mpif.h"
   integer gb_rad_mpistart
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

!  GMS: I moved that outside of the IFDEF, because it is still needed
!       by PUPIL.
   integer ierr
   _REAL_ etmp
   logical belly,nocrst
#include "md.h"
#include "pb_md.h"
#include "box.h"
#include "nmr.h"
#include "memory.h"
#include "files.h"
#include "extra.h"
#include "tgtmd.h"
#include "flocntrl.h"
#include "les.h"
#include "sgld.h"
   integer istart,iend
   _REAL_ evdwex, eelex, virex(3,3)
   _REAL_ escf

   logical qsetup,do_list_update

   _REAL_  enmr(6),devdis(4),devang(4),devtor(4),devpln(4),devplpt(4),devgendis(4),entr,ecap
   _REAL_  x(*),f(*),ene(30),vir(4)
   _REAL_  ener(*) ! offsets in this ener array = offsets in runmd ener - 22
   save ene

#ifdef LES
   _REAL_  :: nrg_bead(nbead)
#endif

   integer ibead,bead_per_node
   _REAL_  virtmp(4),am_Ebnd,am_Eang,am_Edih,am_Evdw,am_Eelt,am_Enb14,am_Ee14,am_Epolar, &
           am_Espring,am_Ederiv

   integer i,m,nttyp,i3
   _REAL_  virvsene,eelt,epol,esurf,edisp,enpol
   _REAL_  epolar,aveper,aveind,avetot,emtot,dipiter,dipole_temp
   integer l_r2x,l_rjx,l_tmp1,l_tmp2,l_tmp3,l_tmp4,l_tmp5
   integer l_tmp6,l_tmp7,l_tmp8,l_jj,l_skipv, l_kvls,l_jvls,l_psi
   integer l_da,l_sumd
   integer qmoffset
   integer newbalance
   save newbalance


#if defined(MPI)
   integer :: status( MPI_STATUS_SIZE )
   _REAL_ :: vel0_nrg_sum
#endif /* MPI */

   call trace_enter( 'force' )

   call timer_start(TIME_FORCE)
   ene(:) = ZERO

   belly = ibelly > 0
   nocrst = .false.
   nttyp = ishft(ntypes*(ntypes+1),-1) !Division by 2


#ifdef MPI
   if (mpi_orig) then

      !     Check to see if we are done yet in mpi_orig case (tec3).
      !     This is done by monitoring the status of an integer notdone.
      !     If notdone .eq. 1 then we keep going.  notdone is set to zero
      !     when we no longer want to call force().  This perhaps is not the
      !     most efficient means to implement this check...

      call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
      if (notdone /= 1) return

      !       Send copies of xyz coords, setbox common block, vir array
      !       and NTNB value to all nodes from master with a broadcast.

      if (numtasks > 1) then

         call mpi_bcast(box,BC_BOXR,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(ntb,BC_BOXI,mpi_integer,0,commsander,ierr)
         call mpi_bcast(vir,3,MPI_DOUBLE_PRECISION,0, &
               commsander,ierr)
         call mpi_bcast(xx(lcrd),3*natom,MPI_DOUBLE_PRECISION, &
               0,commsander,ierr)
         call mpi_bcast(ntnb,1,mpi_integer,0,commsander,ierr)
         if (iabs(ntb) >= 2) then
            call mpi_bcast(xx(l45),3*natom,MPI_DOUBLE_PRECISION, &
                  0,commsander,ierr)
         end if
      end if
   end if

   istart = iparpt(mytaskid) + 1
   iend = iparpt(mytaskid+1)
#else
   istart = 1
   iend = natom
#endif

   if(iamoeba.eq.1) then
      REQUIRE(am_nbead.eq.ncopy)
   end if

   !     ----- ZERO OUT THE ENERGIES AND FORCES -----

   aveper=0.d0
   aveind=0.d0
   avetot=0.d0
   dipiter=0.d0
   dvdl=0.d0
   dipole_temp=0.d0
   enmr(1:6) = 0.d0
   vir(1:4) = 0.d0
   virvsene = 0.d0
   f(1:3*natom+iscale) = 0.d0
#ifdef LES
   if( ipimd>0.or.ineb>0) nrg_all(1:nbead)=0.d0
#endif


! This shuld only happen if *not* using PUPIL,
! because PUPIL resets the nb list itself
#ifndef PUPIL_SUPPORT
   if( igb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_NONBON)
      call timer_start(TIME_LIST)

      call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
               ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
               ix(ibellygp),belly,newbalance,cn1, &
               xx(lvel),xx(lvel2),ntp,xx(l45), qsetup, &
               do_list_update)
      call timer_stop(TIME_LIST)
      call timer_stop(TIME_NONBON)
   end if
#endif

#ifndef DISABLE_NCSU
#  ifdef MPI
   call ncsu_remember_initremd(rem.gt.0.and.mdloop.eq.0)
#  endif
   call ncsu_on_force(x, f, vir)
#endif

   !-----------------------------------------------------------------
   ! QMMM Link atom positioning in main amber coordinate array
   !-----------------------------------------------------------------
   !We need to adjust the link atoms positions
   !in the main amber coordinate array, replacing the mm link pair
   !atom coordinates temporarily.
   if (qmmm_nml%ifqnt) then
     call timer_start(TIME_QMMM)
     call timer_start(TIME_QMMMCOORDSX)
     if( qmmm_nml%idc == 0 ) call adj_mm_link_pair_crd(x)
     call timer_stop(TIME_QMMMCOORDSX)
     call timer_stop(TIME_QMMM)
   end if

   ! ----------------------------------------------------------------
   ! Do weight changes, if requested
   ! ----------------------------------------------------------------

   if (nmropt > 0) &
         call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
         devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb,xx(lnmr01), &
         ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
         cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
         nimprp,nttyp,nphb,natom,natom,ntypes,nres, &
         rad,wel,radhb,welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'WEIT')
   ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints

   epolar = 0.d0

   ! -----------------------------------------------------------------
   ! EGB if igb>0 and /=6 then we need to calculate the GB radii for this
   ! structure
   ! -----------------------------------------------------------------

   ! If we are doing qm_gb=2 then we need to calculate the GB radii
   ! before calling qm_mm

   if( igb > 0 .and. igb /= 6 .and. igb /=10 .and. &
      ( irespa < 2 .or. mod(irespa,nrespai) == 0) ) then
#ifdef MPI
      gb_rad_mpistart = mytaskid+1
#endif
      call timer_start(TIME_EGB)
      call timer_start(TIME_GBRAD1)
      !If qmmm and then this will calculate us the radii for the
      !link atoms and not the mm link pair atoms. The initial radii used are those
      !of the mm link pair's atom type though.

      call egb_calc_radii(igb,natom,x,fs,reff, &
                     onereff,fsmax,rgbmax, rborn, offset,gbalpha, &
                     gbbeta,gbgamma,rbornstat,xx(l188),xx(l189), &
                     xx(l186),xx(l187), gbneckscale, ncopy, rdt &
#ifdef MPI
                     ,gb_rad_mpistart &
#endif
                       )
      call timer_stop(TIME_GBRAD1)
      call timer_stop(TIME_EGB)


   end if

   ! QM/MM Contributions are now calculated before the NON-Bond info.
   ! ----------------------------------------------------------------
   ! Calculate the qm/mm contributions
   ! ----------------------------------------------------------------

   if(qmmm_nml%ifqnt) then

      ! If we are doing periodic boundaries with QM/MM PME then we need to
      ! do the PME calculation twice. First here to get the potential at
      ! each QM atom due to the PME and then again after all the QM is done
      ! to get the MM-MM potential and all of the gradients.

      if(qmmm_nml%qm_pme) then
         ! Ewald force will put the potential into the qm_ewald%mmpot array.
         call timer_start(TIME_EWALD)
         call ewald_force(x,natom,ix(i04),ix(i06),ntypes, &
               xx(l15),cn1,cn2,asol,bsol,eelt,epolar, &
               f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),.true. &
               ,cn114,cn214 &
               )
         call timer_stop(TIME_EWALD)
      endif

      call timer_start(TIME_QMMM)

#ifndef LES
        !========================================================
        !                      REGULAR QMMM
        !========================================================
        if(qmmm_nml%idc>0)then
           call qm_div(x, ix, f, escf, ih(m06))
        else
           call qm_mm(x, natom,qmmm_struct%scaled_mm_charges, &
               f,escf,periodic,reff,onereff, &
               intdiel,extdiel,Arad, cut,qm2_struct%scf_mchg,ih(m06))
        endif
        ene(25) = escf
#endif

        !========================================================
        !                  END REGULAR QMMM
        !========================================================
      call timer_stop(TIME_QMMM)
   end if !if(qmmm_nml%ifqnt)

   !---------------------------------------------------------------
   !END qm/mm contributions
   !---------------------------------------------------------------

#ifdef PUPIL_SUPPORT

   !*****************************************************
   !     Getting the Quantum forces with PUPIL package
   !*****************************************************

!  Reconstructing the simulation cell if there is any change
!  call inipupcell(natms,qcell,cell,xxx,yyy,zzz)
   do iPup=1,3    !vector loop
     do jPup=1,3  !Component loop
       qcell((iPup-1)*3+jPup) = ucell(jPup,iPup)
     enddo
   enddo
!  minimum point of the box ..... we assume (0,0,0) ????
   qcell(10) = 0.0d0
   qcell(11) = 0.0d0
   qcell(12) = 0.0d0

!  temporary vector to wrap the real coordinates to pass through
!  PUPIL interface.
   call get_stack(l_puptmp,3*natom,routine)
   if(.not. rstack_ok)then
     deallocate(r_stack)
     allocate(r_stack(1:lastrst),stat=alloc_ier)
     call reassign_rstack(routine)
   endif
   REQUIRE(rstack_ok)
   do iPup=1,3*natom
     r_stack(l_puptmp + iPup - 1) = x(iPup)
   end do

   if(ntb > 0) then
     call wrap_molecules(nspm,ix(i70),r_stack(l_puptmp))
     if(ifbox == 2) call wrap_to(nspm,ix(i70),r_stack(l_puptmp),box)
   end if

!   write(6,*) ' Updating PUPIL data structures.'

!  Preparing the coordinates,velocity and classic forces
!  to get quantum  force
   do iPup=1,natom
     bs1 = (iPup-1)*9
     bs2 = (iPup-1)*3
     do jPup=1,3
       qcdata(bs1  +jPup) = r_stack(l_puptmp + bs2 +jPup - 1)
!       qcdata(bs1  +jPup) = x(bs2+jPup)
!	   write(6,*) 'Coordinate.',iPup,'==>',realStack(lcrd+bs2+jPup-1),x(bs2+jPup)
!	   write(6,*) 'Velocity...',iPup,'==>',realStack(lvel+bs2+jPup-1)
       qcdata(bs1+3+jPup) = realStack(lvel+bs2+jPup-1)
       qcdata(bs1+6+jPup) = f(bs2+jPup)
     enddo
   enddo

!  Deallocating temporary stack
   call free_stack(l_puptmp,routine)

!  We are going to use the qmmm_nml and qmmm_struct variables to skip quantum atoms
!  in the force calculation
   qmmm_nml%ifqnt = .true.
   if(pupStep .EQ. 0) then
!!!  To keep initial values from the MD step 1
     ierr = 0
     allocate ( pupnb14(numnb14*2),stat=ierr)
     REQUIRE(ierr == 0)

     allocate ( pupbonh(nbonh*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupbona(nbona*3),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheth(ntheth*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( puptheta(ntheta*4),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphih(nphih*5),stat=ierr)
     REQUIRE(ierr == 0)
     allocate ( pupphia(nphia*5),stat=ierr)
     REQUIRE(ierr == 0)
     pupnbonh   = nbonh
     pupnbona   = nbona
     pupntheth  = ntheth
     pupntheta  = ntheta
     pupnphih   = nphih
     pupnphia   = nphia
     pupnumnb14 = numnb14
     call copy_14nb(ix(inb_14),pupnb14,numnb14)
     do iPup = 1,nbonh
       bs1 = iPup-1
       bs2 = bs1*3
       pupbonh(bs2+1)  = ix(iibh +bs1)
       pupbonh(bs2+2)  = ix(ijbh +bs1)
       pupbonh(bs2+3)  = ix(iicbh+bs1)
     enddo
     do iPup = 1,nbona
       bs1 = iPup-1
       bs2 = bs1*3
       pupbona(bs2+1)  = ix(iiba +bs1)
       pupbona(bs2+2)  = ix(ijba +bs1)
       pupbona(bs2+3)  = ix(iicba+bs1)
     enddo
     do iPup = 1,ntheth
       bs1 = iPup-1
       bs2 = bs1*4
       puptheth(bs2+1) = ix(i24  +bs1)
       puptheth(bs2+2) = ix(i26  +bs1)
       puptheth(bs2+3) = ix(i28  +bs1)
       puptheth(bs2+4) = ix(i30  +bs1)
     enddo
     do iPup = 1,ntheta
       bs1 = iPup-1
       bs2 = bs1*4
       puptheta(bs2+1) = ix(i32  +bs1)
       puptheta(bs2+2) = ix(i34  +bs1)
       puptheta(bs2+3) = ix(i36  +bs1)
       puptheta(bs2+4) = ix(i38  +bs1)
     enddo
     do iPup = 1,nphih
       bs1 = iPup-1
       bs2 = bs1*5
       pupphih(bs2+1)  = ix(i40  +bs1)
       pupphih(bs2+2)  = ix(i42  +bs1)
       pupphih(bs2+3)  = ix(i44  +bs1)
       pupphih(bs2+4)  = ix(i46  +bs1)
       pupphih(bs2+5)  = ix(i48  +bs1)
     enddo
     do iPup = 1,nphia
       bs1 = iPup-1
       bs2 = bs1*5
       pupphia(bs2+1)  = ix(i50  +bs1)
       pupphia(bs2+2)  = ix(i52  +bs1)
       pupphia(bs2+3)  = ix(i54  +bs1)
       pupphia(bs2+4)  = ix(i56  +bs1)
       pupphia(bs2+5)  = ix(i58  +bs1)
     enddo
   endif

!  Getting the quantum forces for a specific quantumn domain
   pupStep  = pupStep + 1
   puperror = 0
   pupLevelData = 3
   call getquantumforces(natom,pupLevelData,pupStep,puperror,qcdata,qcell)
   if (puperror.ne.0) then
     write (6,*) "Fatal Error: Error getting quantum forces."
     call mexit(6,1)
   endif
 ! Quantum energy treatment....
   ene(25) = qmEnergy

!  Deleting interactions CL-QZ if a new list of quantum atoms is given
   if (pupQZchange .ne. 0) then

!    ********** Rebuilding the nonbonding 14 list ***************
!    deleting all connectivity between the QM atoms

!      reinitializing internal nb 14 list structures from the beginning
       numnb14= pupnumnb14
       nbonh  = pupnbonh
       nbona  = pupnbona
       ntheth = pupntheth
       ntheta = pupntheta
       nphih  = pupnphih
       nphia  = pupnphia
       call copy_14nb(pupnb14,ix(inb_14),pupnumnb14)
       do iPup = 1,nbonh
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iibh +bs1) = pupbonh(bs2+1)
         ix(ijbh +bs1) = pupbonh(bs2+2)
         ix(iicbh+bs1) = pupbonh(bs2+3)
       enddo
       do iPup = 1,nbona
         bs1 = iPup-1
         bs2 = bs1*3
         ix(iiba +bs1) = pupbona(bs2+1)
         ix(ijba +bs1) = pupbona(bs2+2)
         ix(iicba+bs1) = pupbona(bs2+3)
       enddo
       do iPup = 1,ntheth
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i24  +bs1) = puptheth(bs2+1)
         ix(i26  +bs1) = puptheth(bs2+2)
         ix(i28  +bs1) = puptheth(bs2+3)
         ix(i30  +bs1) = puptheth(bs2+4)
       enddo
       do iPup = 1,ntheta
         bs1 = iPup-1
         bs2 = bs1*4
         ix(i32  +bs1) = puptheta(bs2+1)
         ix(i34  +bs1) = puptheta(bs2+2)
         ix(i36  +bs1) = puptheta(bs2+3)
         ix(i38  +bs1) = puptheta(bs2+4)
       enddo
       do iPup = 1,nphih
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i40  +bs1) = pupphih(bs2+1)
         ix(i42  +bs1) = pupphih(bs2+2)
         ix(i44  +bs1) = pupphih(bs2+3)
         ix(i46  +bs1) = pupphih(bs2+4)
         ix(i48  +bs1) = pupphih(bs2+5)
       enddo
       do iPup = 1,nphia
         bs1 = iPup-1
         bs2 = bs1*5
         ix(i50  +bs1) = pupphia(bs2+1)
         ix(i52  +bs1) = pupphia(bs2+2)
         ix(i54  +bs1) = pupphia(bs2+3)
         ix(i56  +bs1) = pupphia(bs2+4)
         ix(i58  +bs1) = pupphia(bs2+5)
       enddo

       call deleting_qm_atoms()
       qsetup = .true.
!      Setting as current quantum zone
       pupQZchange = 0

   endif
   ! For PUPIL, rebuild the neighbour list 
   ! and zero the charges on QM atoms at every step,
   ! because the QM atoms list may have changed
   if ( igb == 0 .and. iyammp == 0 ) then
         
     ! (for GB: do all nonbondeds together below)
     call timer_start(TIME_NONBON)
     call timer_start(TIME_LIST)
     !do_list_update=.true.         
     call nonbond_list(x,ix(i04),ix(i06),ix(i08),ix(i10), &
                       ntypes,natom/am_nbead,xx,ix,ipairs,ntnb, &
                       ix(ibellygp),belly,newbalance,cn1, &
                       xx(lvel),xx(lvel2),ntp,xx(l45), qsetup, &
                       do_list_update)
     !call qm_zero_charges(x(L15))
     call timer_stop(TIME_LIST)
     call timer_stop(TIME_NONBON)
   end if

#endif /*PUPIL_SUPPORT*/


  if(iamoeba.eq.1) then
      vir(1:4)=0.0
   end if
   ! ----------------------------------------------------------------
   ! Calculate the non-bonded contributions
   ! ----------------------------------------------------------------
   call timer_start(TIME_NONBON)


   call timer_start(TIME_EEXIPS)
   if( ips > 0 ) then
      call eexips(evdwex,eelex,istart,iend, ntb, &
           ix(i04),ix(i08),ix(i10),xx(l15),cn1,cn2,f,x,virex)
      vir(1) = vir(1)+0.5D0*(virips+virex(1,1))
      vir(2) = vir(2)+0.5D0*(virips+virex(2,2))
      vir(3) = vir(3)+0.5D0*(virips+virex(3,3))
   endif
   call timer_stop(TIME_EEXIPS)

   if( igb == 0 .and. iyammp == 0 ) then

      ! (for GB: do all nonbondeds together below)

      call timer_start(TIME_EWALD)

      if ( iamoeba == 1 )then
         call AM_NonBond_eval(natom,x,f,vir,xx,ipairs, &
                               evdw,eelt,epolar,&
                               enb14,ee14,diprms,dipiter)
      else
         if ( induced == 1 )then
            call handle_induced(x,natom,ix(i04),ix(i06),ntypes, &
               xx(l15),cn1,cn2,asol,bsol, &
               eelt,epolar,f,xx,ix,ipairs,xx(lpol), &
               xx(l45),virvsene,ix(i02),ibgwat,nres, &
               aveper,aveind,avetot,emtot,diprms,dipiter,dipole_temp,dt, &
               scee,scnb,ntb)
         else
            call ewald_force(x,natom,ix(i04),ix(i06),ntypes, &
               xx(l15),cn1,cn2,asol,bsol,eelt,epolar, &
               f,xx,ix,ipairs,xx(l45),virvsene, xx(lpol),.false. &
               ,cn114,cn214 &
               )
         end if
      end if ! iamoeba == 1

      call timer_stop(TIME_EWALD)

#ifdef MPI
      if(mytaskid == 0)then
#endif
         ene(2) = evdw
         ene(3) = eelt
         ene(4) = ehb
#ifdef MPI
      else
         ! energies have already been reduced to the master
         ! node in ewald_force, so here we zero out elements
         ! on non-master nodes:
         ene(2) = 0.d0
         ene(3) = 0.d0
         ene(4) = 0.d0
      end if
#endif

   if( ips > 0 )then
         ene(2)=ene(2)+evdwex
         ene(3)=ene(3)+eelex
   endif
   end if  ! ( igb == 0 .and. iyammp == 0 )

   call timer_stop(TIME_NONBON)

   ! ----------------------------------------------------------------
   ! Calculate the other contributions
   ! ----------------------------------------------------------------

   !     -- when igb==10, all nonbonds are done in routine pb_force, and
   !                      all nonpolar interactions are done in np_force:
   !
   !     -- HG put this part here such that "outflag" is known from a call
   !        of pb_force; outflag is needed in the "bond" routine in the case
   !        of ifcap == 2,5 (i.e., ivcap == 1,5)

#ifdef MPI
   if(mytaskid == 0)then
#endif
      if( igb == 10 ) then

         call timer_start(TIME_PBFORCE)
         call pb_force(natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),ix(i10), &
                 cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
         if ( pbgrid ) pbgrid = .false.
         if ( pbinit ) pbinit = .false.
         ene(2) = evdw
         ene(3) = eelt
         ene(4) = epol
         call timer_stop(TIME_PBFORCE)

         call timer_start(TIME_NPFORCE)
         esurf = 0.0d0; edisp = 0.0d0
         if ( ifcap == 0  .and. npopt /= 0 ) &
            call np_force(natom,nres,ntypes,ix(i02),ix(i04),ix(i06),&
                 cn1,cn2,x,f,esurf,edisp)
         if ( pbprint ) pbprint = .false.
         ene(23) = esurf
         ene(26) = edisp
         call timer_stop(TIME_NPFORCE)

      end if  ! ( igb == 10 )

#ifdef MPI
   end if
#endif

!  +---------------------------------------------------------------+
!  |  Bonds with H                                                 |
!  +---------------------------------------------------------------+

   call timer_start(TIME_BOND)

   ! initialize bond virial
   if(ipimd>0) bnd_vir = zero

#ifdef MPI /* SOFT CORE */
   ! zero only once, sc bond energy is sum of H and non-H terms
   sc_ener(1) = 0.0d0
#endif

   if( ntf < 2 ) then

      ebdev = 0.d0
      call bond(nbonh,ix(iibh),ix(ijbh),ix(iicbh),x,xx,ix,f,ene(6))
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesb
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Bonds without H                                              |
!  +---------------------------------------------------------------+

   if( ntf < 3 ) then

      call bond(nbona+nbper,ix(iiba),ix(ijba),ix(iicba),x,xx,ix,f,ene(7))
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesb
      endif
#  endif
#endif
      if (nbonh+nbona > 0) ebdev = sqrt( ebdev/(nbonh+nbona) )
   end if
! Qian add
  ene(6)=ecwb
! end

!  +---------------------------------------------------------------+
!  |  Angles with H                                                |
!  +---------------------------------------------------------------+

   if( ntf < 4 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(2) = 0.0d0
#endif

      eadev = 0.d0
      call angl(ntheth,ix(i24),ix(i26),ix(i28),ix(i30),x,xx,ix,f,ene(8))
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesa
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Angles without H                                             |
!  +---------------------------------------------------------------+

   if( ntf < 5 ) then

      call angl(ntheta+ngper,ix(i32),ix(i34),ix(i36),ix(i38),x,xx,ix,f, &
           ene(9))
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesa
      endif
#  endif
#endif
      if (ntheth+ntheta > 0) eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
   end if
! Qian add
  ene(8)=ecwa
! end
!  +---------------------------------------------------------------+
!  |  Dihedrals with H                                             |
!  +---------------------------------------------------------------+
   ! initialize 14 nb energy virial
   if(ipimd>0) e14vir = zero

   if( ntf < 6 ) then

#ifdef MPI /* SOFT CORE */
      ! zero only once, sc bond energy is sum of H and non-H terms
      sc_ener(3) = 0.0d0
#endif

      call ephi(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
           xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(10),ene(11),ene(12),xx(l190))
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesd
      endif
#  endif
#endif
   end if

!  +---------------------------------------------------------------+
!  |  Dihedrals without H                                          |
!  +---------------------------------------------------------------+

   if( ntf < 7 ) then

      call ephi(nphia+ndper,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
           xx(l15),ix(i04),x,xx,ix,f,dvdl,ene(13),ene(14),ene(15),xx(l190))
! Qian add
  ene(10)=ecwd
! end
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
         ene(24) = ene(24) + elesd
      endif
#  endif
#endif
   end if


   if(iamoeba==1) then
      call AM_VAL_eval(x,f,vir,ene(6),ene(8),ene(10))
   end if

   call timer_stop(TIME_BOND)

   ! --- calculate the position constraint energy ---

   if(natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
       call xconst(natc,entr,ix(icnstrgp),x,f,xx(lcrdr),xx(l60),natom )
       ene(20) = entr
   end if

   if ( itgtmd==1 .and. (nattgtfit > 0 .or. nattgtrms > 0) ) then

      ! Calculate rmsd for targeted md (or minimization) if requested.
      ! All nodes do rms fit, could just be master then broadcast.
      ! All nodes need all coordinates for this.

      call rmsfit(xx(lcrdr),x,xx(lmass),ix(itgtfitgp),  &
                  ix(itgtrmsgp),rmsdvalue,rmsok)

      if (.not.rmsok) then
         if (master) write (6,*) 'Fatal Error calculating RMSD !'
         call mexit(6, 1)
      end if

      call xtgtmd(entr,ix(itgtrmsgp),x,f,xx(lcrdr),xx(lmass))
      ene(20) = entr

   end if

   if(ifcap == 1 .or. ifcap == 2) then
      call capwat(natom,x,f,ecap)
      ene(20) = ene(20) + ecap
   else if(ifcap == 3) then
      write(6,*) 'No energy expression for spherical boundary known yet'
      call mexit(6,1)
   else if(ifcap == 4) then
      write(6,*) 'No energy expression for orthorhombic boundary known yet'
      call mexit(6,1)
      !call orth(natom,ix(ibellygp),x,f,eorth)
      !ene(20) = ene(20) + eorth
   end if
   ! No energy expression for ifcap == 5 given because only
   !    one step of minimization is allowed with this.

   !  (this seems very weird: we have already done an allreduce on molvir
   !  in ewald_force(); this just collects it on processor 0 (with zeroes
   !  on all slave nodes), then later does an allreduce...)

   if( mytaskid == 0 .and. iamoeba == 0 ) then
      vir(1) = vir(1)+0.5d0*molvir(1,1)
      vir(2) = vir(2)+0.5d0*molvir(2,2)
      vir(3) = vir(3)+0.5d0*molvir(3,3)
   end if

   if( igb == 0 .and. iyammp == 0 ) then
      ener(21) = virvsene
      ener(22) = diprms
      ener(23) = dipiter
      ener(24) = dipole_temp
   end if

   !     ---- get the noesy volume penalty energy: ------

   enoe = 0.d0
   if( iredir(4) /= 0 ) then
      call timer_start(TIME_NOE)
      call noecalc(x,f,xx,ix)
      call timer_stop(TIME_NOE)
   end if
   ene(22) = enoe

   !     -- when igb!=0 and igb!=10, all nonbonds are done in routine egb:

   esurf = 0.d0
   if( igb /= 0 .and. igb /= 10 ) then
      call timer_start(TIME_EGB)
      call egb( x,f,rborn,fs,reff,onereff,xx(l15),ix(i04),ix(i06), &
            ix(i08),ix(i10),xx(l190), &
            cut,ntypes,natom,natbel,epol,eelt,evdw, &
            esurf,dvdl,xx(l165),ix(i82),xx(l170),xx(l175),xx(l180), &
            xx(l185), xx(l186),xx(l187),xx(l188),xx(l189),ncopy )

      ene(2) = evdw
      ene(3) = eelt
      ene(4) = epol
      ene(23) = esurf
      ene(21) = dvdl
      call timer_stop(TIME_EGB)
#ifdef MPI
#  ifdef LES
      if(rem == 2) then
        ene(24) = ene(24) + elesp
      endif
#  endif
#endif

   end if  ! ( igb /= 0 .and. (igb /= 10 .or. igb /= 11))

#ifdef APBS
! APBS forces
      if( mdin_apbs ) then
         if (igb /= 6) then
            write(6, '(a)') '&apbs keyword requires igb=6.'
            call mexit(6,1)
         end if
         call timer_start(TIME_PBFORCE)
! in: coords, radii, charges
! out: updated forces (via apbs_params) and solvation energy (pol + apolar)
         if (sp_apbs) then
            call apbs_spenergy(natom, x, f, eelt, enpol)
         else
            call apbs_force(natom, x, f, ene(2), eelt, enpol)
         end if
!         ene(2) =
!         ene(3) =
         ene(4) = eelt
         ene(23) = enpol
         call timer_stop(TIME_PBFORCE)

      end if  ! ( mdin_apbs )
#endif /* APBS */

   if( master ) then
      !  These parts of the NMR energies are not parallelized, so only
      !  are done on the master node:
      eshf = 0.d0
      epcshf = 0.d0
      ealign = 0.d0
      ecsa = 0.d0
      if (iredir(5) /= 0) call cshf(natom,x,f)
      if (iredir(7) /= 0) call pcshift(natom,x,f)
      if (iredir(9) /= 0) call csa1(natom,x,f)
      if (iredir(8) /= 0) call align1(natom,x,f,xx(lmass))
   end if

   !     if freezemol has been set, zero out all of the forces for the
   !     real atoms; (no longer necessary to set ibelly).
   if( ifreeze > 0 ) then
      do i=1,3*natom
         f(i) = 0.d0
      end do
   end if
#ifdef MPI

   call timer_barrier( commsander )
   call timer_start(TIME_COLLFRC)

   !     add force, ene, vir, copies from all nodes
   !            also add up newbalance for nonperiodic.
   call fdist(f,xx(lfrctmp),ene,vir,newbalance)
   call timer_stop(TIME_COLLFRC)

#endif
! Qian add
   ecwb=ene(6)
   ecwa=ene(8)
   ecwd=ene(10)
   ene(6)=0.0
   ene(8)=0.0
   ene(10)=0.0
!   write(84,*)ecwb,ecwa,ecwd ! sum all the node
! end

   ! ---- at this point, the parallel part of the force calculation is
   !      finished, and the forces have been distributed to their needed
   !      locations.  All forces below here are computed redundantly on
   !      all processors, and added into the force vector.  Hence, below
   !      is the place to put any component of the force calculation that
   !      has not (yet) been parallelized.

   ! Calculate the NMR restraint energy contributions, if requested.
   ! (Even though this is not parallelized, it needs to be run on all
   ! threads, since this code is needed for weight changes as well as
   ! for NMR restraint energy analysis.  The whole thing good stand a
   ! major re-write....)

   if (nmropt > 0) &
      call nmrcal(x,f,ih(m04),ih(m02),ix(i02),xx(lwinv),enmr,devdis, &
         devang,devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb, &
         xx(lnmr01),ix(inmr02),xx(l95),31,6,rk,tk,pk,cn1, &
         cn2,asol,bsol,xx(l15),numbnd,numang,nptra-nimprp, &
         nimprp,nttyp,nphb,natom,natom,ntypes,nres, &
         rad,wel,radhb,welhb,rwell,isftrp,tgtrmsd,temp0les,-1,'CALC')
   enoe = ene(22)  ! so all processors now have the full enoe value

#ifdef DSSP
   if( idssp > 0 ) then
      call fdssp( natom,x,f,edssp )
      write(6,*) 'edssp = ', edssp
   else
      edssp = 0.d0
   end if
#endif

   !     ----- CALCULATE TOTAL ENERGY AND GROUP THE COMPONENTS -----

#ifndef LES
   if( igb == 0 ) then
      ene(11) = enb14
      ene(14) = 0.d0
      ene(12) = ee14
      ene(15) = 0.d0
   end if
#endif

   do m = 2,15
      ene(1) = ene(1) + ene(m)
   end do

   ene(1) = ene(1) + epolar + ene(23) + ene(25) + ene(26)

   ene(5) = ene(6)+ene(7)
   ene(6) = ene(8)+ene(9)
   ene(7) = ene(10)+ene(13)
   ene(8) = ene(11)+ene(14)
   ene(9) = ene(12)+ene(15)
   ene(10) = ene(17)+ene(20)+eshf+epcshf+enoe+enmr(1)+enmr(2)+enmr(3)+ &
             enmr(4)+enmr(5)+enmr(6)+ealign+ecsa
#ifdef DSSP
   ene(10) = ene(10) + edssp
#endif
     ! Updated 9/2007 by Matthew Seetin to enable plane-point and plane-plane restraints
   ene(1) = ene(1)+ene(10)
   ene(16) = ene(25) + ene(26)

   !    Here is a summary of how the ene array is used.  For parallel runs,
   !    these values get summed then rebroadcast to all nodes (via
   !    mpi_allreduce).

   !    ene(1):    total energy
   !    ene(2):    van der Waals
   !    ene(3):    electrostatic energy
   !    ene(4):    10-12 (hb) energy, or GB energy when igb.gt.0
   !    ene(5):    bond energy
   !    ene(6):    angle energy
   !    ene(7):    torsion angle energy
   !    ene(8):    1-4 nonbonds
   !    ene(9):    1-4 electrostatics
   !    ene(10):   constraint energy
   !    ene(11-19):  used a scratch, but not needed further below
   !    ene(20):   position constraint energy + cap energy
   !    ene(21):   charging free energy result
   !    ene(22):   noe volume penalty
   !    ene(23):   surface-area dependent energy, or cavity energy
   !    ene(24):   potential energy for a subset of atoms
   !    ene(25):   SCF Energy when doing QMMM
   !    ene(26):   implicit solvation dispersion energy

   !     ----- TRANSFER THE ENERGIES TO THE ARRAY ENER, USED IN PRNTMD -----

   ener(1:10) = ene(1:10)
   ener(11) = epolar
   ener(12) = aveper
   ener(13) = aveind
   ener(14) = avetot
   ener(15) = ene(23)
   ener(17) = ene(21)
   ener(18) = ene(24)
   ener(25) = ene(25)
   ener(28) = ene(26)
   ener(16) = ene(25) + ene(26)   

#ifdef PUPIL_SUPPORT
   !*****************************************************
   !     Closing the qmmm structure consideration
   !*****************************************************
!  Adding the quantum forces from last QM calculation
   do iPup=1,pupqatoms
     bs1 = (abs(pupqlist(iPup))-1)*3
     !write(6,"(a10,2x,i4,3(2x,e16.10))") 'Classic F:',abs(pupqlist(iPup)), &
     !                                              (f(bs1+jPup),jPup=1,3)
     do jPup=1,3
       bs2    = bs1    + jPup
       f(bs2) = f(bs2) + qfpup(bs2)
     enddo
     !write(6,"(a10,2x,i4,3(2x,e16.10))") 'Quantum F:',abs(pupqlist(iPup)), &
     !                                          (qfpup(bs1+jPup),jPup=1,3)
   enddo

!Final forces
!   do iPup=1,natom
!     bs2 = (iPup-1)*3
!     write(6,"(a10,2x,i4,3(2x,e16.10))") 'Total F:',iPup,(f(bs2+jPup),jPup=1,3)
!   enddo

!  Disconnecting qmmmm interactions
   qmmm_nml%ifqnt = .false.

#endif

   ! ----ADD X-RAY TARGET FUNCTION AND GRADIENT
   call cns_xref_run(natom,ih(m04), x,f,ener)

   !     ----- IF BELLY IS ON THEN SET THE BELLY ATOM FORCES TO ZERO -----
   if (belly) call bellyf(natom,ix(ibellygp),f)

!  +---------------------------------------------------------------+
!  |  Interface to EVB                                             |
!  +---------------------------------------------------------------+

#if defined(MPI)
#ifdef LES
!KFW   if( ipimd>0.or.ineb>0) then
   if( nbead > 0 ) then
      call mpi_allreduce ( nrg_all, nrg_bead, nbead, MPI_DOUBLE_PRECISION &
                         , MPI_SUM, commsander, ierr )
      nrg_all(:) = nrg_bead(:)
   end if
#endif
   if( ievb /= 0 ) call evb_ntrfc ( x, f, ener, xx(lmass), ix, ipairs, vel0_nrg_sum )
#endif /* MPI */

   if( ipimd>0 ) then
#ifdef LES
      call pimd_part_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
#ifdef MPI
      if( ievb /= 0 ) then
         nrg_frc(3)= vel0_nrg_sum
         nrg_frc(2)= equal_part + Epot_deriv
         nrg_frc(1)= nrg_frc(3) + nrg_frc(2)
         dlnQ_dl = dvdl
      endif
#endif
#else
      ener(1:26) = ener(1:26)/nbead
      f(1:natom*3) = f(1:natom*3)/nbead
      vir(1:3) = vir(1:3) /nbead
      atvir = atvir/nbead
      e14vir = e14vir/nbead
      bnd_vir = bnd_vir/nbead
      call pimd_full_spring_force(x,f,real_mass,Epot_spring,Epot_deriv,dvdl)
# ifdef MPI
      if(master) call mpi_reduce(ener,totener(23),28,MPI_DOUBLE_PRECISION, &
                                 MPI_SUM,0,commmaster,ierr)
# endif
#endif
      ! Pass dvdl = dV/dl for TI w.r.t. mass.
      if (itimass > 0) ener(17) = dvdl
   end if

   if(ineb>0) then
#ifdef LES
      call part_neb_forces( x, f, xx(lvel), ener(27) )
#else
      call full_neb_forces( xx(lmass), x, f, xx(lvel), ener(1), ener(27) )
#ifdef MPI
      if(sanderrank.eq.0) then
         call mpi_reduce(nebrms,etmp,1,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,0,commmaster,ierr)
         call mpi_reduce(ener,totener(23),28,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,0,commmaster,ierr)
         nebrms = sqrt( etmp/(3*natom*nbead-6*natom) )
      end if
# endif
#endif
   end if

   call timer_stop(TIME_FORCE)
   call trace_exit( 'force' )
   return
end subroutine force


