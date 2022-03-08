! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
#include "assert.h"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ main driver routine for molecular dynamics
subroutine runmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop, qsetup)

   !  Runmd operates in kcal/mol units for energy, amu for masses,
   !     and angstoms for distances.  To convert the input time parameters
   !     from picoseconds to internal units, multiply by 20.455
   !     (which is 10.0*sqrt(4.184)).
   
#if !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD)
   use ncsu_sander_hooks, only : ncsu_on_mdstep => on_mdstep
#endif

   use cmd_vars, only: activate, file_pos_cmd, file_vel_cmd, file_nrg_cmd,  &
                       nstep_cmd, t_cmd, eq_cmd, restart_cmd,  &
                       etot_cmd, eke_cmd, adiab_param, temp_cmd

   use pimd_vars, only: ipimd, ineb, nbead, natomCL, &
                        bnd_vir, Eimp_virial, equal_part, Epot_deriv,  &
                        tau_vol, Epot_spring, NMPIMD, CMD, cartpos, cartvel, &
                        itimass, real_mass

   use lscivr_vars, only: ilscivr, ndof_lsc, natom_lsc, mass_lsc, v2_lsc, &
                          ilsc, x_lsc, f_lsc, dx_lsc

   use nose_hoover_module, only : thermo_lnv, x_lnv, x_lnv_old, v_lnv,  &
                                  f_lnv_p, f_lnv_v, c2_lnv, mass_lnv,  &
                                  Thermostat_init 
 
   use full_pimd_vars, only: totener,totenert,totenert2,mybeadid

   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_mpi, qm2_struct, &
                           element_sym

   use qm2_dftb_module, only: cm3
   use constants, only : third, ten_to_minus3
   use trace
   use stack
   use decomp, only : nat, nrs, decpr, jgroup, indx, irespw, &
#ifdef MPI
   ! -- ti decomp
                      collect_dec, &
#endif
                      checkdec, printdec
   use fastwt
   use bintraj, only: end_binary_frame
   use nblist,only: fill_tranvec,volume,oldrecip,ucell

   use nose_hoover_module, only: Thermostat_switch,  &
                                 Thermostat_integrate_1,  &
                                 Thermostat_integrate_2,  &
                                 Thermostat_hamiltonian
   use nose_hoover_vars, only: file_nhc, nchain, thermo, nthermo, Econserved

#ifdef MPI
   use evb_parm,  only: evb_dyn, nbias
   use evb_data,  only: evb_frc, evb_vel0, evb_bias, evb_nrg, evb_nrg_ave &
                      , evb_nrg_rms, evb_nrg_tmp, evb_nrg_old, evb_nrg_tmp2 &
                      , evb_nrg_old2
   use wigner,    only: rflux
   use remd, only : rem, mdloop, remd_ekmh, repnum, &
                    myeptot, mytemp, mytargettemp,  &
                    hybrid_remd_ene
#  ifdef LES
   use evb_pimd,  only: evb_pimd_dealloc
   use miller,    only: i_qi
#  endif
   use softcore, only: ifsc, sc_dvdl, sc_tot_dvdl, sc_tot_dvdl_partner, &
                       extra_atoms, mix_temp_scaling, sc_pscale, &
                       adj_dvdl_stat, sc_mix_velocities, &
                       sc_nomix_frc, sc_sync_x, sc_print_energies, &
                       calc_softcore_ekin, &
                       sc_ener, sc_ener_ave, sc_ener_rms, sc_lngdyn, &
                       sc_ener_tmp, sc_ener_tmp2, sc_ener_old, sc_ener_old2, &
                       sc_mix_position, sc_print_dvdl_values, &
                       sc_degrees_o_freedom, dynlmb, sc_change_clambda
#endif

   use amoeba_mdin, only: iamoeba
   use amoeba_runmd, only: AM_RUNMD_scale_cell

   implicit none
   character(kind=1,len=5) :: routine="runmd"
   integer   ipairs(*), ix(*)
   _REAL_ xx(*)
   character(len=4) ih(*)
   _REAL_ combination

#ifdef MPI
#  include "parallel.h"
#  include "mpif.h"
#  ifdef LES
   _REAL_  :: fbead(3,natomCL), xbead(3,natomCL)
   integer :: mm, n
#  endif   
   _REAL_ mpitmp(8) !Use for temporary packing of mpi messages.
   integer ist(MPI_STATUS_SIZE), partner, ierr
#endif

! The following variables are needed since nstep and nstlim
!  behave differently in a REMD run.
! In certain places where output is required, total_nstep and total_nstlim 
!  take the place of nstep and nstlim. This allows replica runs to output
!  less often than every exchange.
! They are the absolute step # of the REMD or MD simulation.
   integer total_nstep, total_nstlim

#include "files.h"
#include "md.h"
#include "box.h"
#include "nmr.h"
#include "memory.h"
#include "extra.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "ew_mpole.h"
#include "def_time.h"
#include "extra_pts.h"
#if defined(LES)
#  include "les.h"
#endif
#include "pb_md.h"
#include "random.h"   
#include "sgld.h"
!Antonios added
#include "HB.h"
#include "CHI.h"
#include "CROWD.h"
!Antonios end
! Qian added
#include "debye.h"
! Qian add end
! additional variables for PIMD output
   _REAL_  :: xcmd(3*natomCL),vcmd(3*natomCL)
   integer :: ncmd
   ! for const press PIMD
   _REAL_ tmpvir(3,3),atomvir

   _REAL_ rndfsgld,rndsg1,rndsg2
   _REAL_ sysx,sysy,sysz,sysrange(3,2)
   logical mv_flag

#ifdef MMTSB
#  include "mmtsb.h"
   logical is_done_mmtsb      ! MMTSB replica exchange calculation completed
   _REAL_  lambda_mmtsb       ! MMTSB replica exchange new lambda
   _REAL_  pert_pe_mmtsb      ! MMTSB lambda replica exchange perturbed PE
   _REAL_  temp_mmtsb         ! MMTSB replica exchange new temperature
   _REAL_  unpert_pe_mmtsb    ! MMTSB lambda replica exchange unperturbed PE
#endif

   _REAL_ , dimension(1) :: shkh
   integer, dimension(1) :: ifstwr2
   integer :: nshkh

   integer idx, iatom, iatomCL,m
   _REAL_  Ekin2_tot,tmp,f_lnv
   integer :: idim, ithermo
   _REAL_ :: E_nhc, exp1, exp2, v_sum
   
   logical ivscm
   logical qspatial
   character(len=6)fnam

   ! Constant pH
   integer :: icpselres(0:3), icpselstat(0:3) ! randomly selected residue and state
   type (rand_gen_state) :: cnstph_rand_gen
   
   logical resetvelo
   integer nshak
   _REAL_ ekgs,eold3,eold4,etot_save,ekpbs
   
   logical do_list_update
   logical skip(*),belly,lout,loutfm,erstop,vlim,onstep
   _REAL_ x(*),winv(*),amass(*),f(*),v(*),vold(*), &
         xr(*),xc(*),conp(*)
   _REAL_ enert(51),enert2(51),ener(51),vir(4),ekcmt(4)
   _REAL_ enert_old(51),enert2_old(51),ecopy(51),edvdl(51), &
         enert_tmp(51),enert2_tmp(51),etot_start,edvdl_r(51)
   _REAL_ pres(4),rmu(3),fac(3),onefac(3),clfac
   _REAL_ tma(*)

   _REAL_ tspan,atempdrop,fln,scaltp,scaltpo
   _REAL_ vel,vel2,vcmx,vcmy,vcmz,vmax,vx,vy,vz
   _REAL_ winf,aamass,rterm,ekmh,ekph,ekpht,wfac,rsd,ekav
   _REAL_ fit,fiti,fit2,vscalt
   _REAL_ gammai,c_implic,c_explic,c_ave,sdfac,ekins0
   _REAL_ dtx,dtxinv,dt5,factt,ekin0,ekinp0,dtcp,dttp
   _REAL_ rndf,rndfs,rndfp,boltz2,pconv,tempsu
   _REAL_ xcm(3),acm(3),ocm(3),vcm(3),ekcm,ekrot

   integer nsp(*)
   integer idumar(4)
   integer l_temp
   integer i,j,im,i3,nitp,nits
   integer nstep,nrep,nrek,nren,iend,istart3,iend3
   integer nrx,nr,nr3,ntcmt,izero,istart
   logical ixdump,ivdump,itdump
   logical qsetup
   
   equivalence (pres(1),ener(11)),(ekcmt(1),ener(15))
   equivalence (vir(1),ener(19))
   integer nvalid, nvalidi
   _REAL_ eke,eket
   _REAL_ extent

   _REAL_ xcen,ycen,zcen,extents(3,2),centertest
   _REAL_, allocatable, dimension(:) :: frcti
   integer ier

   _REAL_ small
   data small/1.0d-7/
   data nren/51/

   !--- Variables for Mulliken Charge Printing ---
   _REAL_ :: total_mulliken_charge, mulliken_charge, total_cm3_chg, cm3_chg
   
   !--- VARIABLES FOR DIPOLE PRINTING ---
   integer prndipngrp
   integer prndipfind
   character(len=4) prndiptest

   _REAL_,parameter :: pressure_constant = 6.85695d+4
   ! variables used in constant pressure PIMD
   _REAL_ :: Nkt,centvir,pressure, aa, arg2, poly, e2, e4, e6, e8 
   ! variable used in CMD
   real(8) :: tmp_eke_cmd !Use for temporary packing of mpi messages.

   !==========================================================================
   
   call trace_enter( 'runmd' )

   !     ----- INITIALIZE SOME VARIABLES -----
   
#ifdef MPI
   if( master ) then
      ! If remd, runmd will be called many times, so we dont want to open every
      !  time. For normal md, mdloop will just be 0.
      if (mdloop.eq.0) call amopen(7,mdinfo,'U','F',facc)
   endif
#else
   if( master ) call amopen(7,mdinfo,'U','F','W')
#endif
   vlim = vlimit > small
   ntcmt = 0
   izero = 0
   belly = ibelly > 0
   lout = .true.
   loutfm = ioutfm <= 0
   nr = nrp
   nr3 = 3*nr
   ekmh = 0.d0
#ifdef LES
   ekmhles = 0.d0
#endif
      do_list_update=.false.
#ifdef MPI
   if ( mpi_orig ) then
      istart = 1
      iend = natom
   else
      istart = iparpt(mytaskid) + 1
      iend = iparpt(mytaskid+1)
   end if
#else
   istart = 1
   iend = nr
#endif
   istart3 = 3*istart -2
   iend3 = 3*iend

#ifdef MPI
   if( icfe /= 0 ) then
      allocate( frcti( nr3+3*extra_atoms ), stat = ier )
      REQUIRE( ier == 0 )
   end if
#endif
      
   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx  = nr3
   if (ntwprt > 0) nrx = ntwprt*3
   
   ! Cleanup the velocity if belly run
   if(belly) call bellyf(nr,ix(ibellygp),v)
   
   !=======================================================================
   ! Determine system degrees of freedom (for T scaling, reporting)
   
   ! Call DEGCNT to get the actual number of degrees of freedom for the
   ! solute and solvent. This call returns the correct numbers for belly
   ! simulations and simulations with separate solute/solvent scaling -- dap
   ! "IDUMAR" is dummy array. Used since this routine was also used w/ GIBBS.
   
#ifdef LES
   ! return LES and non-LES degrees,
   ! since separate solvent coupling no longer used
   ! large changes to degcnt were made
   ! cnum is now passed (LES copy number of each atom)
   call degcnt(ibelly,nr,ix(ibellygp),nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,rndfp,rndfles,cnum,temp0les)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.
   ! RNDFLES = # degrees of freedom for LES groups
   
   ! temp0les was init to negative number to signify not to use a LES bath
   ! just do standard code (meaning use solute/solvent baths)
   ! any positive (or zero) means to use LES bath with that target
   
   ! degcnt returns rndfs or rndfles in the rndfles variable
   ! depending on whether a LES bath was specified
   ! do this instead of duplicating call with rndfs or rndfles
   
   if (temp0les < 0.d0) then
      rndfs=rndfles
      rndfles=0.d0
   else
      rndfs=0.d0
   end if

   if (master) then
      write (6,'(a,f8.0)') &
            "# degrees of freedom in non-LES region: ",rndfp
      write (6,'(a,f8.0)') &
            "# degrees of freedom in     LES region: ",rndfles
   end if
   
   !    modify RNDFP to reflect NDFMIN (set in mdread)
   
   rndfp = rndfp - ndfmin

   if (temp0les < 0.d0) then
      rndf = rndfp+rndfs
   else
      rndf = rndfp+rndfles
   end if
   
#else
   
   call degcnt(ibelly,nr,ix(ibellygp),nsolut,nbonh,nbona,0, &
         ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
         idumar,ntc,idumar,0,0,0, &
         idumar,rndfp,rndfs)
   
   ! RNDFP = # degrees of freedom for solute
   ! RNDFS = # degrees of freedom for solvent
   ! RNDF = total number of degrees of freedom.
   
   if (master) then
      write (6,'(a,f8.0)') &
            "|  # of SOLUTE  degrees of freedom (RNDFP): ",rndfp
      write (6,'(a,f8.0)') &
            "|  # of SOLVENT degrees of freedom (RNDFS): ",rndfs
   end if
   !    modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
   rndfp = rndfp - ndfmin + num_noshake
   rndf = rndfp+rndfs
   if (master) then
      write (6,'(a,f8.0,a,i6,a,f8.0)') &
            "|  NDFMIN = ",rndfp, "     NUM_NOSHAKE = ",num_noshake, "     CORRECTED RNDFP = ", rndfp
      write (6,'(a,f8.0)') &
            "|  TOTAL # of degrees of freedom (RNDF) = ", rndf
   end if

#endif

   call fix_degree_count(rndf) ! correct for extra points
   
#ifndef LES
   if(tsgld)then
     ! number of degrees of freedom in the SGLD part
     if(isgsta == 1 .and. isgend == nr)then
       rndfsgld=rndf
     else
       if(isgsta == 1)then
         rndfsgld=0
       else
         call degcnt(ibelly,nr,ix(ibellygp),isgsta-1,nbonh,nbona,0, &
            ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
            idumar,ntc,idumar,0,0,0,idumar,rndsg1,rndsg2)
         rndfsgld=rndsg1
       endif
       call degcnt(ibelly,nr,ix(ibellygp),isgend,nbonh,nbona,0, &
            ix(iibh),ix(ijbh),ix(iiba),ix(ijba),idumar, &
            idumar,ntc,idumar,0,0,0,idumar,rndsg1,rndsg2)
       rndfsgld=rndsg1-rndfsgld
     endif
   endif
#endif

#ifdef MPI /* SOFT CORE */
   if (ifsc /=0 ) call sc_degrees_o_freedom(ndfmin)
#endif

   ! End of degrees of freedom stuff
   !=======================================================================
   
   boltz2 = 8.31441d-3 * 0.5d0
   pconv = 1.6604345d+04  ! factor to convert the pressure kcal/mole to bar
   
   !     ---convert to kcal/mol units
   
   boltz2 = boltz2/4.184d0   ! k-sub-B/2
   dtx = dt*20.455d+00
   dtxinv = 1.0d0 / dtx
   dt5 = dtx * 0.5d0
   pconv = pconv*4.184d0
   
   ! FAC() are #deg freedom * kboltz / 2
   ! multiply by T to get expected kinetic energy
   ! FAC(1) is for total system
   
   fac(1) = boltz2*rndf
   fac(2) = boltz2*rndfp

   if(rndfp < 0.1d0) fac(2) = 1.d-6

#ifdef LES
   ! replaced solvent variables with LES ones
   ! since separate solvent coupling no longer used
   ! ASSUME SAME COUPLING CONSTANT FOR BOTH BATHS, just different target T
   
   ! will also have to accumulate LES and non-LES kinetic energies separately
   
   if (temp0les < 0.d0) then
      fac(3) = boltz2*rndfs
      if(rndfs < 0.1d0) fac(3) = 1.d-6
   else
      fac(3) = boltz2*rndfles
      if(rndfles < 0.1d0) fac(3) = 1.d-6
   end if
#else
   fac(3) = boltz2*rndfs
   if(rndfs < 0.1d0) fac(3) = 1.d-6
#endif
   if ( ipimd==CMD ) then
      if ( eq_cmd ) then
         fac(1) = boltz2 * dble( 3*natomCL )
      else
         fac(1) = boltz2 * dble( 3*(natomCL-1) )
      endif
   endif
   onefac(1) = 1.0d0/fac(1)
   onefac(2) = 1.0d0/fac(2)
   onefac(3) = 1.0d0/fac(3)
   factt = rndf/(rndf+ndfmin)
   
   ! these are "desired" kinetic energies based on
   ! # degrees freedom and target temperature
   ! they will be used for calculating the velocity scaling factor
   
   ekinp0 = fac(2)*temp0
#ifdef LES
   
   ! modified for LES temperature
   
   ekins0=0.d0
   ekinles0=0.d0
   if (temp0les < 0.d0) then
      ekins0 = fac(3) * temp0
      ekin0  = fac(1) * temp0
      if (master) &
            write (6,*) "Single temperature bath for LES and non-LES"
   else
      ekinles0 = fac(3)*temp0les
      ekin0  = ekinp0 + ekinles0
      if (master) then
         write (6,*) "LES particles coupled to separate bath"
         write (6,'(a,f8.2)')"    LES target temperature:    ",temp0les
         write (6,'(a,f8.2)')"    LES target kinetic energy: ",ekinles0
         write (6,'(a,f8.2)')"non-LES target temperature:    ",temp0
         write (6,'(a,f8.2)')"non-LES target kinetic energy: ",ekinp0
      end if
   end if
#else
   ekins0 = fac(3)*temp0
   ekin0  = fac(1)*temp0
#endif

#ifdef LES
   if ( ntt==4 ) call nose_hoover_init_LES(amass,v)
#else
   if ( ntt==4 ) call nose_hoover_init(amass,v)
#endif
   !     LN setup:
   gammai = gamma_ln/20.455d0
   c_implic = 1.d0/(1.d0+gammai*dt5)
   c_explic = 1.d0 - gammai*dt5
   c_ave    = 1.d0+gammai*dt5
   sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
#ifdef LES
   if( temp0les < 0.d0 ) then
      sdfacles = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
   else
      sdfacles = sqrt( 4.d0*gammai*boltz2*temp0les/dtx )
   endif
#endif
    if(tlangv .and. ifbox==0) then
       call get_position(nr,x,sysx,sysy,sysz,sysrange,0)
#ifdef MPI /* SOFT CORE */
       if (ifsc == 1) call sc_mix_position(sysx,sysy,sysz,clambda)
#endif
    end if
   !     Constant pH setup 
   !
   if (icnstph /= 0) then
      ! (separate stream of random numbers so choices stay sync'ed between MPI nodes)
      call amrset_gen(cnstph_rand_gen, ig)
      ! Initialize data for multi-site MC moves
      call cnstphinit(ix(icpstinf),ix(icptrsct),  &
                xx(lcpcrg),x,ix(icphidx),ix(icptpair))
   end if
   
   if (ntt == 1) dttp = dt/tautp
   if (ntp > 0) dtcp = comp * 1.0d-06 * dt / taup
   
   nrek = 4
   nrep = 15
   
   nvalid = 0
   nvalidi = 0
   nstep = 0
   fit = 0.d0
   fiti = 0.d0
   fit2 = 0.d0

   ener = 0.0d0 !Zeros all these arrays of size nren
   enert = 0.0d0
   enert2 = 0.0d0
   enert_old = 0.d0
   enert2_old = 0.d0
   edvdl = 0.d0
   edvdl_r = 0.d0
   ! for PIMD/NMPIMD/CMD/RPMD:
   totenert = 0.d0
   totenert2 = 0.d0

   ener(5) = 1.d0
   ener(6) = 1.d0
   ener(7) = box(1)
   ener(8) = box(2)
   ener(9) = box(3)
   
   ekcmt(1:4) = 0.d0
   nitp = 0
   nits = 0
   
   !=======================================================================
   !     ----- MAKE A FIRST DYNAMICS STEP -----
   !=======================================================================
   !  init = 3:  general startup if not continuing a previous run

   if( ipimd.eq.NMPIMD .or. ipimd.eq.CMD) then
      call trans_pos_cart_to_nmode( x )
   end if
  
   if( init == 3 ) then
      if (ntp > 0 .and. iamoeba==0 .and. ipimd==0) then
         xr(1:nr3) = x(1:nr3)
         
         ! ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
         !       OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
         !       MASS -----
         call ekcmr(nspm,nsp,tma,ekcmt,xr,v,amass,1,nr)
      end if
      
      ! ----- CALCULATE THE FORCE -----

      npbstep = nstep

      !   ---   set irespa to get full energies calculated on step "0":
      irespa = 0
      iprint = 1

      if(ipimd==NMPIMD .or. ipimd==CMD) then
         call trans_pos_nmode_to_cart(x,cartpos)
         call force(xx,ix,ih,ipairs,cartpos,f,ener(23),vir, &
                  xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                  do_list_update)

#if defined(MPI) && defined(LES)
         if ( ievb == 1 .and. i_qi > 0) then
            call evb_umb ( f, cartpos, real_mass, natom, istart3, iend3 )
            if( i_qi == 2 ) call qi_corrf_les ( cartpos, amass )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

         call trans_frc_cart_to_nmode(f)
         i3 = 3*(istart-1)

#if defined(MPI) && defined(LES)
         if ( ievb /= 0 .and. i_qi == 0 ) then
            call evb_umb ( f, x, real_mass, natom, istart3, iend3 )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      else if ( ilscivr == 1 )then
         ! prepare the Hessian Matrix of the potential for the LSC-IVR
         ! at this point, x is the position a bead at equilibrium
         ! initialize the LSC-IVR variables
         natom_lsc = natom
         ndof_lsc = natom * 3

         call lsc_init
         do ilsc = 1, natom_lsc
            mass_lsc(3*ilsc-2) = amass(ilsc)
            mass_lsc(3*ilsc-1) = amass(ilsc)
            mass_lsc(3*ilsc  ) = amass(ilsc)
         end do
         v2_lsc = 0.0d0
         do ilsc = 1, ndof_lsc
            ! ith vector of the Hesian matrix
            x_lsc = 0.0d0
            x_lsc(1:ndof_lsc) = x(1:ndof_lsc)
            x_lsc(ilsc) = x(ilsc) + dx_lsc
            call force(xx,ix,ih,ipairs,x_lsc,f_lsc,ener(23),vir, &
                       xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                       do_list_update)

#ifdef MPI
            call xdist( f_lsc )
#endif
            v2_lsc(1:ndof_lsc,ilsc) = f_lsc(1:ndof_lsc)
         enddo

         call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update)


#ifdef MPI
         call xdist(f)
#endif
         ! 2nd derivative of the potential:
         do ilsc = 1, ndof_lsc
            v2_lsc(1:ndof_lsc,ilsc) = &
                ( f(1:ndof_lsc) - v2_lsc(1:ndof_lsc,ilsc) )/dx_lsc
         end do

         ! get the iniital position of the momentum:
         call lsc_xp(x,v)

      else

         ! -- ti decomp
         if(idecomp > 0) then
            decpr = .false.
            if(mod(nstep+1,ntpr) == 0) decpr = .true.
         end if
         call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
                  xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
                  do_list_update)
 
#ifdef MPI
         if ( ievb /= 0 ) then
#ifdef LES
            call evb_umb_primitive ( f, x, real_mass, natom, istart, iend )
#else
            call evb_umb_primitive ( f, x, amass, natom, istart, iend )
#endif
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      endif

      if (icnstph /= 0 ) call cnstphwrite(ix(icpresst),0,ix(icptrsct))



      ! This FORCE call does not count as a "step". CALL NMRDCP to decrement
      ! local NMR step counter
      call nmrdcp

#ifdef MPI /* SOFT CORE */
      ! If softcore potentials are used, collect their dvdl contributions:
      if ( ifsc /= 0 ) then
         call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, commsander, ierr)
         sc_dvdl=0.0d0 ! zero for next step
         call mpi_reduce(sc_ener, sc_ener_tmp, 8, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, commsander, ierr)
         sc_ener(1:8) = sc_ener_tmp(1:8)
      end if
      if ( ifsc == 2 ) then
         ! If this is a perturb to nothing run, scale forces and calculate dvdl
         call sc_nomix_frc(f,nr3,ener,clambda)
         if( numtasks>1 ) then
            call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            call mpi_bcast(ener,51,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         end if
      end if

      if( icfe /= 0 )then
         ! ---free energies using thermodynamic integration (icfe /= 0)

         if( master ) then
            !  --- first, send the forces and energy to your partner:
            partner = ieor(masterrank,1)
            call mpi_sendrecv( f, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                               frcti, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                               partner, 5, commmaster, ist, ierr )
            call mpi_sendrecv( ener, 51, MPI_DOUBLE_PRECISION, partner, 5, &
                               ecopy, 51, MPI_DOUBLE_PRECISION, partner, 5, &
                               commmaster, ist, ierr)
            ! exchange sc-dvdl contributions between masters
            call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, &
                               5, sc_tot_dvdl_partner, 1, &
                               MPI_DOUBLE_PRECISION, partner, 5, &
                               commmaster, ist, ierr )
            if( masterrank==0 ) then
               call mix_frcti(frcti,ecopy,f,ener,nr3,clambda,klambda)
            else
               call mix_frcti(f,ener,frcti,ecopy,nr3,clambda,klambda)
            end if
         end if

         if( numtasks>1 ) then
            call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            call mpi_bcast(ener,51,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         end if

      end if
#endif /* MPI SOFT CORE */

      irespa = 1
      
      ! Reset quantities depending on TEMP0 and TAUTP (which may have been
      ! changed by MODWT during FORCE call).
      ! Recalculate target kinetic energies.
      
      ekinp0 = fac(2) * temp0

#ifdef LES
      
      ! modified for LES temperature, not solvent
      
      ekins0 = 0.d0
      ekinles0 = 0.d0
      if (temp0les < 0.d0) then
         ekins0 = fac(3) * temp0
         ekin0 = fac(1) * temp0
      else
         ekinles0 = fac(3) * temp0les
         ekin0 = ekinp0 + ekinles0
      end if
#else
      ekins0 = fac(3) * temp0
      ekin0 = fac(1) * temp0
#endif

      if (ntt == 1) dttp = dt / tautp
      
      if (ntp > 0) then
         ener(10) = volume
         ener(42) = tmass / (0.602204d0*volume)
         if( iamoeba == 0 ) then
            ekcmt(4) = 0.d0
            vir(4) = 0.d0
            pres(4) = 0.d0
            do m = 1,3
               ekcmt(m) = ekcmt(m) * 0.5d0
               ekcmt(4) = ekcmt(4) + ekcmt(m)
               vir(4) = vir(4) + vir(m)
               pres(m) = (pconv+pconv) * (ekcmt(m)-vir(m)) / volume
               pres(4) = pres(4) + pres(m)
            end do
            pres(4) = pres(4) / 3.d0
         end if
      end if
      
      ntnb = 0
      i3 = 0
      tempsu = 0.0d0
      
#ifdef LES
      ! added LES tempsu (actual LES sum of m*v**2 )
      tempsules = 0.0d0
#endif
      eke_cmd = 0.d0
      do j = 1,nrp
         winf = winv(j) * dt5
         aamass = amass(j)
         do m = 1,3
            i3 = i3+1
            rterm = v(i3)*v(i3) * aamass
#ifdef LES
            if (temp0les < 0.d0) then
               tempsu = tempsu + rterm
               if (ipimd.eq.CMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
                  eke_cmd = eke_cmd + aamass*v(i3)*v(i3) 
               endif
            else
               if (cnum(j) == 0) then
                  tempsu = tempsu + rterm
               else
                  tempsules = tempsules + rterm
               end if
            end if
#else
            if(ipimd.eq.CMD.and.mybeadid==1) then
               eke_cmd = eke_cmd + aamass*v(i3)*v(i3)
            end if
            tempsu = tempsu + rterm
#endif
            if(ipimd.ne.NMPIMD.and.ipimd.ne.CMD) v(i3) = v(i3) - f(i3) * winf
            if (vlim) v(i3) = sign(min(abs(v(i3)),vlimit),v(i3))
         end do
      end do

#ifdef MPI /* SOFT CORE */
      if ( ifsc /= 0 ) call calc_softcore_ekin(nrp,amass,v,v)
#endif
      
      do im=1,iscale
         v(nr3+im) = v(nr3+im) - f(nr3+im) * dt5 / scalm
         tempsu = tempsu + scalm * v(nr3+im)*v(nr3+im)
      end do
      ener(3) = tempsu * 0.5d0

#ifdef LES
      
      ! added for LES temperature using old solvent variable for ener(4)
      
      if (temp0les < 0.d0) then
         ener(4)=0.d0
         ener(2) = ener(3)
         ! for CMD:
         if( ipimd > 0 ) then
            ener(4) = equal_part + Epot_deriv ! "virial" estimate of KE
            ener(1) = ener(4)+ener(23)
         else
            ener(1) = ener(2)+ener(23)
         endif
         if (ipimd.eq.CMD) then
            ener(2) = eke_cmd*0.5d0
            ener(4) = ener(2)
         endif
      else
         ener(4) = tempsules * 0.5d0
         ener(2) = ener(3) + ener(4)
      end if
#else
      ! for better output for parallel PIMD/NMPIM/CMD/RPMD
      if (ipimd>0) then
         ener(1) = 0.
         ener(2) = 0.
         ener(3) = 0.
         ener(4) = 0.
         ener(10) = 0.
      endif
      ener(2) = ener(3)
      ener(1) = ener(2)+ener(23)
#endif

      if(ntt == 1) then
#ifdef LES
         if (temp0les >= 0.d0) then
            ekmh = max(ener(3),fac(2)*10.d0)
            ekmhles = max(ener(4),fac(3)*10.d0)
         else
            ekmh = max(ener(3),fac(1)*10.d0)
         end if
#else
         ekmh = max(ener(3),fac(1)*10.d0)
#endif
      end if
      
      ! endif for init=3 check
      
   end if  ! ( init == 3 )
   
   !-------------------------------------------------------------------------
   ! init = 4:  continuation of a previous trajectory
   !            this code also done for init=3
   !
   ! Note: if the last printed energy from the previous trajectory was
   !       at time "t", then the restrt file has velocities at time
   !       t + 0.5dt, and coordinates at time t + dt
   !-------------------------------------------------------------------------
   
   ekmh = 0.0d0
#ifdef LES
   ekmhles = 0.0d0
#endif

   i3 = 0
   do j = 1,nrp
      aamass = amass(j)
      do m = 1,3
         i3 = i3+1
         rterm = v(i3)*v(i3) * aamass
#  ifdef LES
         ! use copy number, not solute/solvent
         if (temp0les < 0.d0) then
            ! 1 bath
            ekmh = ekmh + rterm
         else
            if (cnum(j) == 0) then
               ekmh = ekmh + rterm
            else
               ekmhles = ekmhles + rterm
            end if
         end if
#  else
         ekmh = ekmh + rterm
#  endif
      end do
   end do

#ifdef MPI /* SOFT CORE */
   if ( ifsc /= 0 ) call calc_softcore_ekin(nrp,amass,v,v)
#endif
   
   do im=1,iscale
      ekmh = ekmh + scalm*v(nr3+im)*v(nr3+im)
   end do
   ekmh = ekmh * 0.5d0
#ifdef LES
   ekmhles = ekmhles * 0.5d0
#endif

   do i=1,nr3+iscale
      vold(i) = v(i)
   end do
   
   if (init == 4) then

   else
      
      !-------------------------------------------------------------------
      !           PRINT THE INITIAL ENERGIES AND TEMPERATURES
      !-------------------------------------------------------------------
      
      if (nstep <= 0 .and. master .and. facc /= 'A') then
         if(tsgld)then
           ener(48)=sgft
           ener(49)=tempsg
         endif
         rewind(7)
         call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)
#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener)
         if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif

         !--- BEGIN DIPOLE PRINTING CODE ---
         ! See code further on for comments-explanations
         call nmlsrc('dipoles',5,prndipfind)
         if(prndipfind /= 0 ) then
            write(6,*) '------------------------------- DIPOLE INFO ----------------------------------'
            write(6,9018) nstep,t
            read (5,'(a)') prndiptest
            call rgroup(natom,natc,nres,prndipngrp,ix(i02),ih(m02), &
                 ih(m04),ih(m06),ih(m08),ix(icnstrgp), &
                 jgroup,indx,irespw,npdec, &
                 xx(l60),xx(lcrdr),0,0,0,idecomp,5,.false.)
            rewind(5)
            if(prndipngrp > 0) then
               call printdip(prndipngrp,ix(icnstrgp),xx(lcrd), &
                    xx(l15),xx(linddip),xx(Lmass), natom)
            end if
            write(6,*) '----------------------------- END DIPOLE INFO --------------------------------'
         end if
         !--- END DIPOLE PRINTING CODE ---

         if (nmropt > 0) then
            call nmrptx(6)
            call nmrptx(7)
         end if
         call amflsh(7)
      end if
      if (nstlim == 0) return
      init = 4
   end if


   if(ntp > 0 .and. ipimd > 0 ) then
      REQUIRE(ipimd.eq.NMPIMD)
#ifdef LES
      call part_setup_cnst_press_pimd(Nkt,tau_vol)
#else     
      call full_setup_cnst_press_pimd(Nkt,tau_vol)
#endif
      e2 = 1.0/ (2.0*3.0)
      e4 = e2 / (4.0*5.0)
      e6 = e4 / (6.0*7.0)
      e8 = e6 / (8.0*9.0)
      x_lnv = log( box(1)*box(2)*box(3) ) / 3 
   end if

   ! For CMD.
   if ( ipimd==CMD ) then

      if ( .not.eq_cmd ) then

         ! De-activate thermostat for path-centroid.
#ifdef LES
         do iatom = 1, natom
         do idim  = 1, 3
            if ( cnum(iatom)==0 .or. cnum(iatom)==1 )  then
               activate = .false.
            else
               activate = .true.
            end if
            call Thermostat_switch(thermo(idim,iatom),activate)
         enddo
         enddo
         if ( .not.restart_cmd ) then
            ! Scale path-centroid velocity and set total momentum equal to zero.
            call part_scale_vel_centroid(v,amass,istart,iend)
            nstep_cmd = 0
            t_cmd = 0.d0
         else
            t_cmd = t
            nstep_cmd = int( t / dt )
         end if
#else
         if ( mybeadid.eq.1 ) then
            activate = .false.
         else
            activate = .true.
         end if
         do iatom = 1, natom
         do idim  = 1, 3
            call Thermostat_switch(thermo(idim,iatom),activate)
         enddo
         enddo
         if ( .not.restart_cmd ) then
            ! Scale path-centroid velocity and set total momentum equal to zero.
            call full_scale_vel_centroid(v,amass,istart,iend)
            nstep_cmd = 0
            t_cmd = 0.d0
         else
            nstep_cmd = nstep
            t_cmd = t
         end if
#endif /* LES */

      else

         nstep_cmd = nstep
         t_cmd = t

      end if

   end if  ! ipimd.eq.CMD and adiab_param<1.d0

#ifdef MPI
   ! If this is a replica run and we are on exchange > 1, restore the 
   !    old ekmh value since it was reset after we left runmd last time.
   !    DAN ROE: Only for ntt==1??
   if (rem>0.and.mdloop>1) then
      ! if (master) write(6,'(a,f16.6)') "Resetting ekmh to ",remd_ekmh
      ! write(50+worldrank,'(a,f16.6)') "Resetting ekmh to ",remd_ekmh
      ekmh=remd_ekmh
   endif
#endif

  
   !=======================================================================
   !     ----- MAIN LOOP FOR PERFORMING THE DYNAMICS STEP -----
   !           (at this point, the coordinates are a half-step "ahead"
   !           of the velocities; the variable EKMH holds the kinetic
   !           energy at these "-1/2" velocities, which are stored in
   !           the array VOLD.)
   !=======================================================================
   
   260 continue
   onstep = mod(irespa,nrespa) == 0
  
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  |  EVB reactive flux                                            |
!  +:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::+
!  |  Driver for coordinating backward and forward propagation as  |
!  |  well as for enforcing stopping criteria                      |
!  +---------------------------------------------------------------+

#if defined(MPI)
   if( ievb /= 0 .and. trim( adjustl( evb_dyn) ) == "react_flux" ) then
      REQUIRE( ipimd.eq.0 .or. ipimd.eq.NMPIMD )
      call react_flux ( x, v, f, winv, tempi * factt, dt5, dtx &
                      , nr, nstep, nstlim )
   endif
#endif
 
   !---------------------------------------------------------------
   !  ---Step 1a: do some setup for pressure calculations:
   !---------------------------------------------------------------

   if (ntp > 0 .and. iamoeba == 0 .and. ipimd==0) then
      ekcmt(1) = 0.d0
      ekcmt(2) = 0.d0
      ekcmt(3) = 0.d0
      xr(1:nr3) = x(1:nr3)
      
      ! ----- CALCULATE THE CENTER OF MASS ENERGY AND THE COORDINATES
      !       OF THE SUB-MOLECULES WITH RESPECT TO ITS OWN CENTER OF
      !       MASS -----
      
      call timer_start(TIME_EKCMR)
      call ekcmr(nspm,nsp,tma,ekcmt,xr,v,amass,istart,iend)
#ifdef MPI
      call trace_mpi('mpi_allreduce', &
            3,'MPI_DOUBLE_PRECISION',mpi_sum)
# ifdef USE_MPI_IN_PLACE
      call mpi_allreduce(MPI_IN_PLACE,ekcmt,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
# else
      call mpi_allreduce(ekcmt,mpitmp,3, &
            MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
      ekcmt(1) = mpitmp(1)
      ekcmt(2) = mpitmp(2)
      ekcmt(3) = mpitmp(3)
# endif
#endif
      call timer_stop(TIME_EKCMR)
   end if

   ! Constant pH setup
   if (icnstph /= 0) then
      if (ntnb.eq.1) then
         call cnstphupdatepairs(x,ix(icptrsct),ix(icphidx),ix(icptpair))
      end if
      if (mod(irespa,ntcnstph) == 0) then
         call cnstphbeginstep(ix(icpstinf),ix(icpresst), &
               ix(icptrsct), xx(lcpene),xx(lcpcrg),xx(l190),icpselres, &
               icpselstat, cnstph_rand_gen,ix(icptpair))
      end if
   end if

   !--------------------------------------------------------------
   !  ---Step 1b: Get the forces for the current coordinates:
   !--------------------------------------------------------------

   npbstep = nstep
   iprint = 0
   if( nstep == 0 .or. nstep+1 == nstlim ) iprint = 1

   if ( ipimd==NMPIMD .or. ipimd==CMD) then
      call trans_pos_nmode_to_cart(x,cartpos)
      call force(xx,ix,ih,ipairs,cartpos,f,ener(23),vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update)

#if defined(MPI) && defined(LES)
         if ( ievb == 1 .and. i_qi > 0) then
            call evb_umb ( f, cartpos, real_mass, natom, istart3, iend3 )
            if( i_qi == 2 ) call qi_corrf_les ( cartpos, amass )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

      call trans_frc_cart_to_nmode(f)

#if defined(MPI) && defined(LES)
         if ( ievb /= 0 .and. i_qi == 0 ) then
            call evb_umb ( f, x, real_mass, natom, istart3, iend3 )
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif

#endif

   else
      ! -- ti decomp
      if(idecomp > 0) then
         decpr = .false.
         if(mod(nstep+1,ntpr) == 0) decpr = .true.
      end if
      call force(xx,ix,ih,ipairs,x,f,ener(23),vir, &
               xx(l96),xx(l97),xx(l98),xx(l99), qsetup, &
               do_list_update)
!Antonios added:
!        ENER(23)=ENER(23)+ECHI
!        write(72,*)ECHI   This is the one that works. Prints one value for every step.
! It was used for the serial version before we decided to add the echi at evdw in ew_directe.h 
!Antonios end

#if defined(MPI)
         if ( ievb /= 0 ) then
#ifdef LES
            call evb_umb_primitive ( f, x, real_mass, natom, istart, iend )
#else
            call evb_umb_primitive ( f, x, amass, natom, istart, iend )
#endif
            evb_nrg(1) = evb_frc%evb_nrg
            evb_nrg(2) = evb_vel0%evb_nrg
            if( nbias > 0 ) evb_nrg(3) = sum( evb_bias%nrg_bias(:) )
         endif
#endif

   endif

   ! Constant pH transition evaluation
   if ((icnstph /= 0) .and. (mod(irespa,ntcnstph) == 0)) then
      call cnstphendstep(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
            ix(icptrsct), xx(lcpene),xx(lcpcrg),xx(l190),xx(l15),ener(39), &
            icpselres,icpselstat,cnstph_rand_gen)
      call cnstphwrite(ix(icpresst),icpselres,ix(icptrsct))
   end if

#ifdef MPI
   ! If softcore potentials are used, collect their dvdl contributions:
   if ( ifsc /= 0 ) then
      call mpi_reduce(sc_dvdl, sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_dvdl=0.0d0 ! zero for next step
      call mpi_reduce(sc_ener, sc_ener_tmp, 8, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commsander, ierr)
      sc_ener(1:8) = sc_ener_tmp(1:8)
   end if
   if ( ifsc == 2 ) then
      ! If this is a perturb to nothing run, scale forces and calculate dvdl
      call sc_nomix_frc(f,nr3,ener,clambda)

      if( numtasks>1 ) then
         call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(ener,51,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      end if
   end if

   if ( icfe /= 0 )then

      ! --- free energies using thermodynamic integration (icfe /= 0)

      !  --- first, send the forces, energy, and virial to your partner:

      if( master ) then
         partner = ieor(masterrank,1)
         call mpi_sendrecv( f, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                  frcti, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, partner, 5, &
                  commmaster, ist, ierr )
         call mpi_sendrecv( ener, 51, MPI_DOUBLE_PRECISION, partner, 5, &
                            ecopy, 51, MPI_DOUBLE_PRECISION, partner, 5, &
                            commmaster, ist, ierr)
         ! exchange sc-dvdl contributions between masters:
         call mpi_sendrecv( sc_tot_dvdl, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    sc_tot_dvdl_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                    commmaster, ist, ierr )

         ! ---- collect statistics for free energy calculations:

         if( onstep ) then
            if( masterrank==0 ) then
               if( klambda == 1 ) then
                  do i=1,nren
                     edvdl(i) = edvdl(i) - ener(i) + ecopy(i)
                     edvdl_r(i) = edvdl_r(i) - ener(i) + ecopy(i)
                  end do
               else
                  clfac = klambda*(1.d0 - clambda)**(klambda-1)
                  do i=1,nren
                     edvdl(i) = edvdl(i) - (ener(i) - ecopy(i))*clfac
                     edvdl_r(i) = edvdl_r(i) - (ener(i) - ecopy(i))*clfac
                  end do
               end if
            else
               if( klambda == 1 ) then
                  do i=1,nren
                     edvdl(i) = edvdl(i) + ener(i) - ecopy(i)
                     edvdl_r(i) = edvdl_r(i) + ener(i) - ecopy(i)
                  end do
               else
                  clfac = klambda*(1.d0 - clambda)**(klambda-1)
                  do i=1,nren
                     edvdl(i) = edvdl(i) + (ener(i) - ecopy(i))*clfac
                     edvdl_r(i) = edvdl_r(i) + (ener(i) - ecopy(i))*clfac
                  end do
               end if
            end if
            ! This includes the sc-dvdl contribution into the vdw-part 
            !    and potential energy parts of the dvdl-statistics
            if (ifsc == 1) then
               call adj_dvdl_stat(edvdl(24), edvdl_r(24), clambda)
               call adj_dvdl_stat(edvdl(23), edvdl_r(23), clambda)
            end if
         end if

         if( masterrank==0 ) then
            call mix_frcti(frcti,ecopy,f,ener,nr3,clambda,klambda)
         else
            call mix_frcti(f,ener,frcti,ecopy,nr3,clambda,klambda)
         endif
      endif

      if( numtasks>1 ) then
         call mpi_bcast(f,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(ener,51,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      end if
   end if  ! ( icfe /= 0 )

#endif /* MPI */

   ! Reset quantities depending on TEMP0 and TAUTP (which may have been
   ! changed by MODWT during FORCE call).
   ekinp0 = fac(2)*temp0

#ifdef LES
   ! TEMP0LES may have changed too
   
   ekinles0=0.d0
   ekins0=0.d0
   if (temp0les >= 0.d0) then
      ekinles0 = fac(3)*temp0les
      ekin0 = ekinp0 + ekinles0
   else
      ekins0 = fac(3)*temp0
      ekin0 = fac(1)*temp0
   end if
#else
   ekins0 = fac(3)*temp0
   ekin0 = fac(1)*temp0
#endif

   if (ntt == 1) dttp = dt/tautp
   
   !  Pressure coupling:
   if (ntp > 0.and.ipimd>0) then
      REQUIRE(ipimd.eq.NMPIMD)
      centvir=0.0

#ifdef LES
      do iatom=istart,iend
         if(cnum(iatom).eq.0.or.cnum(iatom).eq.1) then
            centvir=centvir-x(3*iatom-2)*f(3*iatom-2)
            centvir=centvir-x(3*iatom-1)*f(3*iatom-1)
            centvir=centvir-x(3*iatom  )*f(3*iatom)
         end if
      end do
#else
      if(mybeadid.eq.1) then
         do iatom=istart,iend
            centvir=centvir-x(3*iatom-2)*f(3*iatom-2)
            centvir=centvir-x(3*iatom-1)*f(3*iatom-1)
            centvir=centvir-x(3*iatom  )*f(3*iatom)
         end do         
      end if
#endif /* LES */

      if(iamoeba.eq.1) then
         atomvir=vir(1)+vir(2)+vir(3)
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,centvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,atomvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#  else
         call mpi_allreduce(centvir,mpitmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         centvir=mpitmp(1)
         tmp=0.0
         call mpi_allreduce(atomvir,tmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         atomvir=tmp
#  endif
#endif
      else
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,centvir,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,bnd_vir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         call mpi_allreduce(MPI_IN_PLACE,e14vir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#    ifndef LES
         call mpi_allreduce(MPI_IN_PLACE,atvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#    endif
#  else  
         call mpi_allreduce(centvir,tmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         centvir=tmp
         tmpvir=0.0
         call mpi_allreduce(bnd_vir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         bnd_vir=tmpvir
         tmpvir=0.0
         call mpi_allreduce(e14vir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         e14vir=tmpvir
#    ifndef LES
         tmpvir=0.0
         call mpi_allreduce(atvir,tmpvir,9,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         atvir=tmpvir
#    endif
#  endif
#endif
         atomvir=0.0
         atomvir=atomvir+atvir(1,1)+bnd_vir(1,1)+e14vir(1,1)
         atomvir=atomvir+atvir(2,2)+bnd_vir(2,2)+e14vir(2,2)
         atomvir=atomvir+atvir(3,3)+bnd_vir(3,3)+e14vir(3,3)
      end if
      pressure = (Nkt*3.0-centvir-(atomvir-Eimp_virial))/(3.0*volume)
      f_lnv_p = (pressure-pres0/pconv)*volume*3.0
   end if


   if (ntp > 0) then
      ener(10) = volume
      ener(42) = tmass / (0.602204d0*volume)
      if( iamoeba == 0 .and. ipimd==0 ) then
         ekcmt(4) = 0.d0
         vir(4) = 0.d0
         pres(4) = 0.d0
         do m = 1,3
            ekcmt(m) = ekcmt(m)*0.5d0
            ekcmt(4) = ekcmt(4)+ekcmt(m)
            vir(4) = vir(4)+vir(m)
            pres(m) = (pconv+pconv)*(ekcmt(m)-vir(m))/volume
            pres(4) = pres(4)+pres(m)
         end do
         pres(4) = pres(4)/3.d0
      end if
   end if

#ifdef MPI
! ------====== REMD ======------ 
! If rem>0 and mdloop==0, this is the first sander call and we don't want to 
!   actually do any MD or change the initial coordinates.
! Exit here since we only wanted to get the potential energy for the first 
!  subrem exchange probability calc.
   if (rem > 0 .and. mdloop.eq.0) then
      if (master) write (6,'(a,i3)') &
         'REMD: Exiting runmd after getting initial energies for replica',repnum
      goto 480 ! Go to the end of the runmd loop.
   endif  ! (rem>0 and mdloop==0)
#endif

   !----------------------------------------------------------------
   !  ---Step 1c: do randomization of velocities, if needed:
   !----------------------------------------------------------------
   ! ---Assign new random velocities every Vrand steps, if ntt=2

   resetvelo=.false.
   if (vrand /= 0 .and. ntt == 2) then
      if (mod((nstep+1),vrand) == 0) resetvelo=.true.
   end if

#ifdef MMTSB
   if ( mmtsb_switch == mmtsb_temp_rex .and. mmtsb_is_exchanged )  &
      resetvelo = .true.
#endif
   
   if (resetvelo) then
      ! DAN ROE: Why are only the masters doing this? Even if the velocities 
      !  are broadcast to the child processes, the wont the different # of random
      !  calls put the randomg num generators out of sync, or do we not care?
      
      if (master) then
         write (6,'(a,i8)') 'Setting new random velocities at step ', &
               nstep + 1
         call setvel(nr,v,winv,temp0*factt,init,iscale,scalm)

#ifdef MPI /* SOFT CORE */
         ! Make sure all common atoms have the same v (that of V0) in TI runs:
         if (icfe /=0 .and. ifsc /=0) call sc_sync_x(v,nr3)
#endif

#ifdef LES
     
         ! newvel call is fixed for the dual target temperatures

         if (temp0les >= 0.d0.and.temp0 /= temp0les) then
            vscalt = sqrt (temp0les/temp0) 
            do j=1,natom
              if(cnum(j) > 0) then
                i3 = 3*(j-1)
                v(i3+1) = v(i3+1) * vscalt
                v(i3+2) = v(i3+2) * vscalt
                v(i3+3) = v(i3+3) * vscalt
              endif
            end do
         end if
#endif
         if (ibelly > 0) call bellyf(nr,ix(ibellygp),v)
      end if
# ifdef MPI
      call trace_mpi('mpi_bcast',3*natom,'MPI_DOUBLE_PRECISION',0)
      call mpi_bcast(v, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
# endif

      ! At this point in the code, the velocities lag the positions
      ! by half a timestep.  If we intend for the velocities to be drawn 
      ! from a Maxwell distribution at the timepoint where the positions and 
      ! velocities are synchronized, we have to correct these newly 
      ! redrawn velocities by backing them up half a step using the 
      ! current force.
      ! Note that this fix only works for Newtonian dynamics.
      if( gammai==0.d0.and.(ipimd.ne.NMPIMD.or.ipimd.ne.CMD)) then                         
         i3 = 3*(istart-1)
         do j=istart,iend
            wfac = winv(j) * dt5
            v(i3+1) = v(i3+1) - f(i3+1)*wfac
            v(i3+2) = v(i3+2) - f(i3+2)*wfac
            v(i3+3) = v(i3+3) - f(i3+3)*wfac
            i3 = i3+3
         end do
      end if
      
   end if  ! (resetvelo)

   call timer_start(TIME_VERLET)

   !-----------------------------------------------------
   !  ---Step 2: Do the velocity update:
   !-----------------------------------------------------
   
   !step 2a: apply quenched MD if needed.  This is useful in NEB>0
   if (vv==1) call quench(f,v)
   
   !  Car-Parrinello on dipoles:  note that the (small?) kinetic energy
   !       of the dipoles is included in the epol energy
   if ( induced == 1 .and. indmeth == 3 ) call cp_dips(natom,xx(lpol),xx,dt)
   
   ! Nose'-Hoover thermostat (1st step).
   if ( ntt == 4 ) then

      Ekin2_tot = 0.d0
      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = dtx/amass(j)
         do idim = 1, 3
#ifdef LES
            if( ntp>0.and.ipimd.eq.NMPIMD .and. &
               (cnum(j).eq.0.or.cnum(j).eq.1) ) then
#else
            if(ntp>0.and.ipimd.eq.NMPIMD.and.mybeadid.eq.1) then
#endif
               exp1 = exp(-dt5*thermo(idim,j)%v(1)-dt5*v_lnv*c2_lnv)
               Ekin2_tot = Ekin2_tot + amass(j)*v(i3+idim)*v(i3+idim)
            else
               exp1 = exp( -dt5 * thermo(idim,j)%v(1) )
            end if
            exp2 = exp1*exp1
            vold(i3+idim)=v(i3+idim)
            v(i3+idim) = v(i3+idim) * exp2 + f(i3+idim) * wfac * exp1
         end do
         i3 = i3+3
      end do

      if(ntp>0.and.ipimd.eq.NMPIMD) then
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,Ekin2_tot,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
#  else
         call mpi_allreduce(Ekin2_tot,mpitmp,1,MPI_DOUBLE_PRECISION, &
                            mpi_sum,commworld,ierr)
         Ekin2_tot=mpitmp(1)
#  endif
#endif
         f_lnv_v = Ekin2_tot*(c2_lnv-1)
         tmp = exp(-dt5*thermo_lnv%v(1))
         v_lnv = tmp*(tmp*v_lnv+dtx*(f_lnv_v+f_lnv_p)/mass_lnv)
      end if

      call Thermostat_integrate_1(nchain,thermo,nthermo,dtx,ntp)

   else if( gammai == 0.d0 ) then

      !       ---Newtonian dynamics:
      
      ! Applying guiding force effect:
      if(tsgld)then
         call sgmdw(natom,istart,iend,dtx,rndfsgld,amass,winv,f,v,xx(lvsg))
         ener(48)=sgft
         ener(49)=tempsgi
      endif  

      i3 = 3*(istart-1)
      do j=istart,iend
         wfac = winv(j) * dtx
         v(i3+1) = v(i3+1) + f(i3+1)*wfac
         v(i3+2) = v(i3+2) + f(i3+2)*wfac
         v(i3+3) = v(i3+3) + f(i3+3)*wfac
         i3 = i3+3
      end do
      
   else if (tsgld) then

      !  Using SGLD algorithm:
      call sgldw(natom,istart,iend,dtx,temp0,rndfsgld,amass,winv, &
                 f,v,xx(lvsg))
      ener(48)=sgft
      ener(49)=tempsgi

   else  !  gamma_ln .ne. 0, which also implies ntt=3 (see mdread.f)

      !       ---simple model for Langevin dynamics, basically taken from
      !          Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      !          Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      !          should not be there.)

      !  Update Langevin parameters, since temp0 might have changed:
         sdfac = sqrt( 4.d0*gammai*boltz2*temp0/dtx )
! Qian add, save temperature for debye-huckel potential
  tempdebye = temp0
! Qian add end

#  ifdef LES
         sdfacles = sqrt( 4.d0*gammai*boltz2*temp0les/dtx )
#  endif


#ifdef MPI /* SOFT CORE */
      if (ifsc == 1) then
         call sc_lngdyn(winv,amass,v,f,sdfac,c_explic,c_implic, &
                        istart, iend, nr, dtx,worldrank)
      else
#endif

#ifndef LES
      if( ipimd.eq.1.or.ineb.eq.1) then
         do j=1,(mybeadid-1)*nr
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
         end do
      end if
#endif

      i3 = 3*(istart-1)
      do j=1,nr
         if( j<istart .or. j>iend ) then
         ! In order to generate the same sequence of pseudorandom numbers that
         ! you would using a single processor you have to go through the atoms 
         ! in order.  The unused results are thrown away
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            cycle
         end if
         wfac = winv(j) * dtx
         aamass = amass(j)         
#  ifdef LES
         if (temp0les >= 0 .and. temp0 /= temp0les) then
            if( cnum(j) == 0 ) then
               rsd = sdfac*sqrt(aamass)
               call gauss( 0.d0, rsd, fln )
               v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
               call gauss( 0.d0, rsd, fln )
               v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
               call gauss( 0.d0, rsd, fln )
               v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
            else
               rsdles = sdfacles*sqrt(aamass)
               call gauss( 0.d0, rsdles, fln )
               v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
               call gauss( 0.d0, rsdles, fln )
               v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
               call gauss( 0.d0, rsdles, fln )
               v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
            endif
         else
            rsd = sdfac*sqrt(aamass)
            call gauss( 0.d0, rsd, fln )
            v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
            call gauss( 0.d0, rsd, fln )
            v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
            call gauss( 0.d0, rsd, fln )
            v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
         endif
#  else
         rsd = sdfac*sqrt(aamass)
         call gauss( 0.d0, rsd, fln )
         v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
         call gauss( 0.d0, rsd, fln )
         v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
         call gauss( 0.d0, rsd, fln )
         v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic
#  endif
      
         i3 = i3+3
      end do

#ifndef LES
      if( ipimd>0) then
         do j=mybeadid*nr+1,nbead*nr
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
         end do
      end if
#endif

#ifdef MPI /* SOFT CORE */
      end if ! for (ifsc==1) call sc_lngdyn 
#endif
   end if  ! ( gammai == 0.d0 )

   !     --- consider vlimit
   
   if (vlim.and.ipimd==0) then
      vmax = 0.0d0
      do i=istart3,iend3
         vmax = max(vmax,abs(v(i)))
         v(i) = sign(min(abs(v(i)),vlimit),v(i))
      end do

      ! Only violations on the master node are actually reported
      ! to avoid both MPI communication and non-master writes.
      if (vmax > vlimit) then
         if (master) then
            write(6,'(a,i6,a,f10.4)') 'vlimit exceeded for step ',nstep, &
              '; vmax = ',vmax
         end if
      end if
   end if
   
   do im=1,iscale
      v(nr3+im) = (v(nr3+im) + f(nr3+im)*dtx/scalm)
   end do

   !-------------------------------------------------------------------
   !   Step 3: update the positions, putting the "old" positions into F:
   !-------------------------------------------------------------------


# ifdef LES
   if(ntp>0.and.ipimd.eq.NMPIMD) then
      aa = exp(dt5*v_lnv)
      arg2 = v_lnv*dt5*v_lnv*dt5
      poly = 1.0d0+arg2*(e2+arg2*(e4+arg2*(e6+arg2*e8)))
   endif

   i3 = 3*(istart-1)
   do j=istart,iend
      if(ntp>0.and.ipimd.eq.NMPIMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
         do idim = 1, 3
            f(i3+idim)=x(i3+idim)
            x(i3+idim)=aa*(x(i3+idim)*aa+v(i3+idim)*poly*dtx)
         enddo
      else
         do idim = 1, 3
            f(i3+idim) = x(i3+idim)
            x(i3+idim) = x(i3+idim)+v(i3+idim)*dtx
         enddo
      endif
      i3 = i3 + 3
   enddo

# else

   if(ntp>0.and.ipimd.eq.NMPIMD.and.mybeadid==1) then
      aa = exp(dt5*v_lnv)
      arg2 = v_lnv*dt5*v_lnv*dt5
      poly = 1.0d0+arg2*(e2+arg2*(e4+arg2*(e6+arg2*e8)))
      do i3=istart3,iend3
         f(i3)=x(i3)
         x(i3)=aa*(x(i3)*aa+v(i3)*poly*dtx)
      end do
   else
      do i3 = istart3,iend3
         f(i3) = x(i3)
         x(i3) = x(i3)+v(i3)*dtx
      end do
   end if

# endif /* LES */

   !Nose'-Hoover thermostat (2nd step).
   if ( ntt==4 ) then
      call Thermostat_integrate_2(nchain,thermo,nthermo,dtx,ntp)
      E_nhc = Thermostat_hamiltonian(nchain,thermo,nthermo)
   end if

   do i = 1,iscale
      f(nr3+i) = x(nr3+i)
      x(nr3+i) = x(nr3+i)+v(nr3+i)*dtx
   end do

   call timer_stop(TIME_VERLET)
   
   if (ntc /= 1) then

      !-------------------------------------------------------------------
      !   Step 4a: if shake is being used, update the new positions to fix
      !      the bond lengths.
      !-------------------------------------------------------------------
   
      call timer_start(TIME_SHAKE)
      qspatial=.false.
      call shake(nrp,nbonh,nbona,0,ix(iibh),ix(ijbh),ix(ibellygp), &
      winv,conp,skip,f,x,nitp,belly,ix(iifstwt),ix(noshake), &
      shkh,qspatial)
      call quick3(f,x,ix(iifstwr),natom,nres,ix(i02))
      if(nitp == 0) then
         erstop = .true.
         goto 480
      end if

      !-----------------------------------------------------------------
      !   Step 4b:   Now fix the velocities and calculate KE
      !-----------------------------------------------------------------
      
      !  ---re-estimate the velocities from differences in positions:

      if( .not.(ipimd==NMPIMD.and.ipimd==CMD.and.mybeadid.ne.1) ) then
         v(istart3:iend3) = (x(istart3:iend3)-f(istart3:iend3)) * dtxinv
      end if

      call timer_stop(TIME_SHAKE)
   end if
   call timer_start(TIME_VERLET)

   if(ineb>0.and.(mybeadid==1.or.mybeadid==nbead) ) then
      x(1:3*natom)=f(1:3*natom)
   end if

   if( ntt == 1 .or. onstep ) then

      !-----------------------------------------------------------------
      !   Step 4c: get the KE, either for averaging or for Berendsen:
      !-----------------------------------------------------------------
      
      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0
#ifdef LES
      ekeles = 0.d0
      ekphles = 0.d0 
#endif
      eke_cmd = 0.d0

      if (gammai == 0.0d0) then
         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
#ifdef LES
               if (temp0les < 0.d0) then
                  eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2
                  ekph = ekph + aamass*v(i3)**2
                  if(ipimd.eq.CMD.and.(cnum(j).eq.0.or.cnum(j).eq.1)) then
                    eke_cmd = eke_cmd + aamass*0.25d0*(v(i3)+vold(i3))**2
                  endif
               else
                  if (cnum(j) == 0) then
                     eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2
                     ekph = ekph + aamass*v(i3)**2
                  else
                     ekeles = ekeles + aamass*0.25d0*(v(i3)+vold(i3))**2
                     ekphles = ekphles + aamass*v(i3)**2
                  end if
               end if

#else
               eke = eke + aamass*0.25d0*(v(i3)+vold(i3))**2

               if(mybeadid==1) then
                  eke_cmd = eke_cmd + aamass*0.25d0*(v(i3)+vold(i3))**2
               end if
               ! try pseudo KE from Eq. 4.7b of Pastor, Brooks & Szabo,
               ! Mol. Phys. 65, 1409-1419 (1988):

               ekpbs = ekpbs + aamass*v(i3)*vold(i3)
               ekph = ekph + aamass*v(i3)**2
#endif
            end do
         end do

      else

         i3 = 3*(istart-1)
         do j=istart,iend
            aamass = amass(j)
            do m = 1,3
               i3 = i3+1
#ifdef LES
               if (temp0les < 0.d0) then
                  eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
               else
                  if (cnum(j) == 0) then
                     eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
                  else
                     ekeles = ekeles + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
                  end if
               end if
#else
               eke = eke + aamass*0.25d0*c_ave*(v(i3)+vold(i3))**2
#endif
            end do
         end do
      end if ! (if gammai == 0.0d0)

#ifdef MPI
      
      !  ---   sum up the partial kinetic energies:
      
#  ifdef LES
      if ( ipimd.eq.CMD ) then
         call mpi_reduce(eke_cmd,tmp_eke_cmd,1,MPI_DOUBLE_PRECISION, &
              mpi_sum,0,commsander,ierr)
         eke_cmd = tmp_eke_cmd
      endif
      if ( .not. mpi_orig .and. numtasks > 1 ) then
        if ( temp0les < 0 ) then
          mpitmp(1) = eke
          mpitmp(2) = ekph
#    ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,mpitmp,2, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(1)
          ekph = mpitmp(2)
#    else
          call mpi_allreduce(mpitmp,mpitmp(3),2, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(3)
          ekph = mpitmp(4)
#    endif
        else
          mpitmp(1) = eke
          mpitmp(2) = ekph
          mpitmp(3) = ekeles
          mpitmp(4) = ekphles
#    ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,mpitmp,4, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(1)
          ekph = mpitmp(2)
          ekeles = mpitmp(3)
          ekphles = mpitmp(4)
#    else
          call mpi_allreduce(mpitmp,mpitmp(5),4, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
          eke = mpitmp(5)
          ekph = mpitmp(6)
          ekeles = mpitmp(7)
          ekphles = mpitmp(8)
#    endif
        endif
      end if
#  else
      if ( .not. mpi_orig .and. numtasks > 1 ) then
         call trace_mpi('mpi_allreduce', &
               1,'MPI_DOUBLE_PRECISION',mpi_sum)
         mpitmp(1) = eke
         mpitmp(2) = ekph
         mpitmp(3) = ekpbs
#    ifdef USE_MPI_IN_PLACE
         call mpi_allreduce(MPI_IN_PLACE,mpitmp,3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(1)
         ekph = mpitmp(2)
         ekpbs = mpitmp(3)

#    else
         call mpi_allreduce(mpitmp,mpitmp(4),3, &
               MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
         eke = mpitmp(4)
         ekph = mpitmp(5)
         ekpbs = mpitmp(6)
#    endif
      end if
#  endif
      
      ! Calculate Ekin of the softcore part of the system
      if (ifsc /= 0 ) call calc_softcore_ekin(nrp,amass,v,vold)
      
#endif
      
      !         --- all processors handle the "extra" variables:

      do im=1,iscale
         eke = eke + scalm*0.25d0*(v(nr3+im)+vold(nr3+im))**2
         ekpbs = ekpbs + scalm*v(nr3+im)*vold(nr3+im)
         ekph = ekph + scalm*v(nr3+im)**2
      end do
      
      eke = eke * 0.5d0
      ekph = ekph * 0.5d0
      ekpbs = ekpbs * 0.5d0
#ifdef LES
      ekeles = ekeles * 0.5d0
      ekphles = ekphles * 0.5d0
#endif
      
      if( ntt == 1 ) then
#ifdef LES
         
         if (temp0les < 0.d0) then
            scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))
         else
            scaltp = sqrt(1.d0+2.d0*dttp*(ekinp0-eke)/(ekmh+ekph))
            scaltles = sqrt(1.d0+2.d0*dttp*(ekinles0-ekeles)/(ekmhles+ekphles))
         end if
#else
         
         !    --- following is from T.E. Cheatham, III and B.R. Brooks,
         !        Theor. Chem. Acc. 99:279, 1998.
         
         scaltp = sqrt(1.d0 + 2.d0*dttp*(ekin0-eke)/(ekmh+ekph))

         !    --- following is the "old" (amber7 and before) method:

         !  scaltpo = sqrt(1.d0 + dttp*(ekin0/ekph - 1.d0))
         !  write(6,*) 'scaltp: ',2.d0*dttp*(ekin0-eke)/(ekmh+ekph), &
         !            dttp*(ekin0/ekmh - 1.d0)

         !  following line reverts to the "old" behavior:
         !  scaltp = scaltpo

#endif

#ifdef MPI /* SOFT CORE */
         if (icfe /= 0) then
            if (ifsc == 1) then
               if (master) then
                  ! Linearly combine the scaling factors from both processes
                  ! the combined factor is broadcast to all nodes
                  ! the subroutine also correctly scales the softcore atom v's
                  call mix_temp_scaling(scaltp,clambda,v)
               end if
               call mpi_bcast(scaltp,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            end if
         end if
#endif

         do j = istart,iend
            i3=(j-1)*3+1
#ifdef LES
            if (temp0les > 0.d0 .and. cnum(j) /= 0 ) then
               v(i3  ) = v(i3  )*scaltles
               v(i3+1) = v(i3+1)*scaltles
               v(i3+2) = v(i3+2)*scaltles
            else
               v(i3  ) = v(i3  ) *scaltp
               v(i3+1) = v(i3+1) *scaltp
               v(i3+2) = v(i3+2) *scaltp
            end if
#else
            v(i3  ) = v(i3  ) *scaltp
            v(i3+1) = v(i3+1) *scaltp
            v(i3+2) = v(i3+2) *scaltp
#endif
         end do
         do im=1,iscale
            v(nr3+im) = v(nr3+im)*scaltp
         end do
      end if  ! (ntt == 1 )
      
   end if  ! ( ntt == 1 .or. onstep; end of step 4c )

   !-----------------------------------------------------------------
   !   Step 5:  several tasks related to dumping of trajectory information
   !-----------------------------------------------------------------
   
   itdump = .false.             ! Write coordinates this step?
   ivdump = .false.             ! Write velocities this step?
   ixdump = .false.             ! Write restart this step?
   ivscm  = .false.             ! Do com removal this step?

   !  --- Determine if trajectory, velocity, or restart
   !      writing is imminent, or if the center of mass
   !      motion will be removed.
   !      These require xdist of velocities or dipoles in parallel runs:
   !
   ! Modified so that when running REMD, writing can occur less often
   !  than exchanges (e.g. ntwx > nstlim)
   ! DAN ROE: Added two new variables, total_nstep and total_nstlim.
   !          For non-REMD runs, total_nstep=nstep+1 and total_nstlim=nstlim 
   !           just like before.
   !          For REMD runs, total_nstep=(mdloop-1)*nstlim+nstep+1, where
   !           mdloop is the current exchange - this is the current
   !           replica exchange MD step. total_nstlim=numexchg*nstlim, which is
   !           the maximum number of REMD steps.
   total_nstep=nstep+1
   total_nstlim=nstlim

#ifdef MPI
   if (rem>0) then
      total_nstep=(mdloop-1)*nstlim+nstep+1
      total_nstlim=nstlim*numexchg
   endif
#endif
   if (ntwx>0) itdump = mod(total_nstep,ntwx) == 0 ! Trajectory coords
   if (ntwv>0) ivdump = mod(total_nstep,ntwv) == 0 ! Velocity
   if( ntwr /= 0 ) then
      if ( mod(total_nstep, ntwr ) == 0 ) ixdump = .true. ! Restart
   endif
   if( total_nstep >= total_nstlim ) ixdump = .true. ! Final restart
   if( mod(total_nstep,nscm) == 0 ) ivscm =.true. ! C.o.M. removal
   if (ntwv == -1 .and. itdump) ivdump = .true. !Combined crdvel file


#ifdef MPI

   !-----------------------------------------------------------------
   !  --- now distribute the coordinates, and if necessary, dipoles and vel:
   !-----------------------------------------------------------------

   call timer_barrier( commsander )
   call timer_stop_start(TIME_VERLET,TIME_DISTCRD)
   if ( .not. mpi_orig .and. numtasks > 1 ) then
      call xdist(x)
   end if
   ! dac/knut change: force the coordinates to be the same on both masters.
   ! For certain compilers, addition may not be strictly commutative, so
   ! the forces on group 0 may be different by roundoff from the forces on 
   ! group 1.  This can lead to divergent trajectories.  The interval at
   ! which they are resynchronized is hard-wired here to 20, which seems to
   ! work fine in our tests.
   if( icfe /= 0 .and. mod(nstep+1,20) == 0 ) then

      ! In dual-topology this is done within softcore.f
      if (ifsc /= 1) then
         if( master ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                     0,commmaster,ierr)
      else
         if( master ) call sc_sync_x(x,nr3)
      end if
      if( numtasks>1 ) call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION, &
                                      0,commsander,ierr)
   end if
   call timer_stop(TIME_DISTCRD)

#endif  /* MPI */

   !           ----fix lone pair positions:
   if( numextra > 0 )call local_to_global(x,xx,ix)

#ifdef MPI
   if ( .not. mpi_orig .and. numtasks > 1 ) then
      call timer_start(TIME_DISTCRD)
      
      !  ---Here we provide every processor a full copy of the velocities
      !     for removal of center of mass motion, or for archiving.
      !     (Note: this is actually over-kill: for example, only the master
      !     node really needs the velocities for archiving.  But the extra
      !     overhead of doing it this way is probably small in most cases.)
      
      if( ivdump .or. ivscm .or. ixdump) then
         call xdist(v)
      endif
  
      if( ixdump .and. (induced == 1 .and. indmeth == 3 ) )then
         call xdist(xx(ldipvel))
         call xdist(xx(linddip))
      end if
      call timer_stop(TIME_DISTCRD)
   end if
   call timer_start(TIME_VERLET)

   !     ========================= END AMBER/MPI =========================
#endif  /* MPI */
   
   !-------------------------------------------------------------------
   !   Step 6: zero COM velocity if requested; used for preventing
   !   ewald "block of ice flying thru space" phenomenon, or accumulation
   !   of rotational momentum in vacuum simulations
   !-------------------------------------------------------------------

   if (ivscm) then
      if (mod(nstep,nsnb) == 0) ntnb = 1
      if( ifbox == 0 ) then
        if(tlangv)then
           ! Get current center of the system 
           call get_position(nr,x,vcmx,vcmy,vcmz,sysrange,0)

#ifdef MPI /* SOFT CORE */
           if (ifsc == 1) call sc_mix_position(vcmx,vcmy,vcmz,clambda)
#endif
           ! Center the system to the original center
           call re_position(nr,0,x,x, &
                            vcmx,vcmy,vcmz,sysx,sysy,sysz,sysrange,mv_flag,0)
        else
           !  ---Non-periodic simulation: remove both translation and rotation.
           !     Back the coords up 1/2 step, so that the correspond to the
           !     velocities; temporarily store in the F() array:
           f(1:nr3) = x(1:nr3) - v(1:nr3)*dt5
           !     --- now compute the com motion, remove it, and recompute (just
           !         to check that it is really gone.....)
           call cenmas(nr,f,v,amass,ekcm,xcm,vcm,acm,ekrot,ocm,4)
           call stopcm(nr,f,v,xcm,vcm,ocm, .true.)
           call cenmas(nr,f,v,amass,ekcm,xcm,vcm,acm,ekrot,ocm,4)
        endif
      else
        if(.not.tlangv)then        
         !    ---Periodic simulation: just remove the translational velocity:
         vcmx = 0.d0
         vcmy = 0.d0
         vcmz = 0.d0
         j = 1
         do i = 1, 3*natom,3
            aamass = amass(j)
            vcmx = vcmx + aamass * v(i)
            vcmy = vcmy + aamass * v(i+1)
            vcmz = vcmz + aamass * v(i+2)
            j = j + 1
         end do
         vcmx = vcmx * tmassinv
         vcmy = vcmy * tmassinv
         vcmz = vcmz * tmassinv
         vel2 = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz
         atempdrop = 0.5d0 * tmass * vel2 * onefac(1) !onefac(1) = 1.0d0/fac(1)
         vel = sqrt(vel2)
         if ( master ) write (6,'(a,f15.6,f9.2,a)') &
            'check COM velocity, temp: ',vel,atempdrop, '(Removed)'
         do i = 1, 3*natom, 3
            v(i)   = v(i)   - vcmx
            v(i+1) = v(i+1) - vcmy
            v(i+2) = v(i+2) - vcmz
         end do

#ifdef MPI /* SOFT CORE */
         if (icfe==1) then
            if (ifsc==1) then
               if (master) then
                  call sc_mix_velocities(v,nr3,clambda)
               end if
               call mpi_bcast(v,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
            end if
         end if
#endif

       end if  ! ( tlangv  )
      end if  ! ( ifbox == 0 )
   end if  ! (ivscm)
  
   !  Also zero out the non-moving velocities if a belly is active:
   if (belly) call bellyf(nr,ix(ibellygp),v)

   !-----------------------------------------------------------------
   !  --- put current velocities into VOLD
   !-----------------------------------------------------------------
   
   vold(istart3:iend3) = v(istart3:iend3)
   do im=1,iscale
      vold(nr3+im) = v(nr3+im)
   end do

   !-------------------------------------------------------------------
   !  Step 7: scale coordinates if constant pressure run:
   !-------------------------------------------------------------------
   if( ntp > 0 .and. ipimd > 0) then
      x_lnv_old = x_lnv
      x_lnv = x_lnv_old + v_lnv * dtx
      rmu(1:3) = exp( ( x_lnv - x_lnv_old ) )
      box(1:3) = box(1:3) * rmu(1:3)
      volume = box(1) * box(2) * box(3)
      ener(7:9) = box(1:3)
      ! only for NMPIMD in sander.LES
      ! (in sander.MPI volume, pressure and density printed in pimdout)
#ifdef LES
      ener(10) = volume
#else
      ener(10) = 0.
      totener(10) = volume
#endif
      call redo_ucell(rmu)
      call fill_tranvec()
      call ew_pscale(natom,x,amass,nspm,nsp,2)
   end if

   if( iamoeba == 0 ) then
      if (ntp == 1) then
         rmu(1) = (1.d0-dtcp*(pres0-pres(4)))**third
         rmu(2) = rmu(1)
         rmu(3) = rmu(1)
      else if (ntp > 1) then
         rmu(1) = (1.d0-dtcp*(pres0-pres(1)))**third
         rmu(2) = (1.d0-dtcp*(pres0-pres(2)))**third
         rmu(3) = (1.d0-dtcp*(pres0-pres(3)))**third
      end if
      if (ntp > 0) then
         box(1:3) = box(1:3)*rmu(1:3)
         ener(7:9) = box(1:3)
         
         !    WARNING!!   This is not correct for non-orthogonal boxes if
         !    NTP > 1 (i.e. non-isotropic scaling).  Currently general cell
         !    updates which allow cell angles to change are not implemented.
         !    The viral tensor computed for ewald is the general Nose Klein,
         !    however the cell response needs a more general treatment.
      
         call redo_ucell(rmu)
         ! keep tranvec up to date, rather than recomputing each MD step.
         call fill_tranvec()  ! tranvec is dependent on only ucell

#ifdef MPI /* SOFT CORE */
         ! if softcore potentials and the dual topology approach are used
         ! C.O.M. scaling has to be changed to account for different masses 
         ! of the same molecule in V0 and V1. This is quite inefficient and is
         ! therefore done in a separate routine in softcore.f
         ! only both masters actually do the computation for ifsc==1
         ! the scaled coordinates are then broadcast to the nodes
         if (icfe /= 0 .and. ifsc == 1) then
            if (master) then
               call sc_pscale(natom,x,amass,nspm,nsp,oldrecip,ucell)
            end if
            call mpi_bcast(x,nr3,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         else
#endif
            call ew_pscale(natom,x,amass,nspm,nsp,npscal)
#ifdef MPI /* SOFT CORE */
         end if
#endif
         if (ntr > 0 .and. nrc > 0) &
            call ew_pscale(natom,xc,amass,nspm,nsp,npscal)
      endif
      if (ipimd==NMPIMD.and.ntp>0) then
         ekcmt(4) = 0.d0
         vir(4) = 0.d0
         pres(4) = pressure*pconv
       endif
   else
      if (ntp>0) then
         if (ipimd==0) then     ! for classical AMOEBA
            ekcmt(4) = eke      !  for printing in prntmd()
            vir(4) = vir(1) + vir(2) + vir(3)
            pres(4) = (pressure_constant/volume)*(2.d0*eke - vir(4)) / 3.d0
         elseif (ipimd==NMPIMD) then     ! for NMPIMD AMOEBA
            ekcmt(4) = 0.d0
            vir(4) = 0.d0
            pres(4) = pressure*pconv
         endif
         call AM_RUNMD_scale_cell(natom,pres(4),dt,pres0,taup,x)
         call fill_tranvec()
      end if
   end if
   
#ifdef LES
   ener(3) = eke
   ener(4) = ekeles
   ener(2) = ener(3) + ener(4)
   if (ntt == 1 .and. onstep) then
      if ( temp0les < 0 ) then
        ekmh = max(ekph,fac(1)*10.d0)
      else
        ekmh = max(ekph,fac(2)*10.d0)
        ekmhles = max(ekphles,fac(3)*10.d0)
      endif
   end if

   if( ipimd > 0 ) then
      ener(4) = equal_part + Epot_deriv  ! "virial" estimate of KE
      ener(1) = ener(4)+ener(23)
   endif
#else
   if( ipimd > 0 ) then
      ! use a "virial" estimator for the KE, rather than one derived from the 
      ! bead velocities:
      totener(4) = equal_part + Epot_deriv
   else
      ener(4) = ekpbs + ener(23)  ! Pastor, Brooks, Szabo conserved quantity
                                  ! for harmonic oscillator: Eq. 4.7b of Mol.
                                  ! Phys. 65:1409-1419, 1988
   endif
   ener(3) = eke
   ener(2) = ener(3)
   if (ntt == 1 .and. onstep) then
      ekmh = max(ekph,fac(1)*10.d0)
   end if
#endif
   
   !     ---if velocities were reset, the KE is not accurate; fudge it
   !        here to keep the same total energy as on the previous step.
   !        Note that this only affects printout and averages for Etot
   !        and KE -- it has no effect on the trajectory, or on any averages
   !        of potential energy terms.
   
   if( resetvelo ) ener(2) = etot_save - ener(23)
   
   !     --- total energy is sum of KE + PE:
   
   if( ipimd >  0 ) then
      totener(1) = totener(4)+totener(23)
      etot_save = totener(2)+totener(23)
      if (ipimd==CMD) then
         etot_cmd = eke_cmd*0.5 + ener(23)
         totener(1)= etot_cmd
         ener(1) = etot_cmd
         ener(2) = eke_cmd*0.5
         ener(4) = ener(2)
      endif
   else
      ener(1) = ener(2)+ener(23)
      etot_save = ener(1)
   end if

   !-------------------------------------------------------------------
   !  Step 8:  update the step counter and the integration time:
   !-------------------------------------------------------------------
   
   nstep = nstep+1
   t = t+dt

   !For CMD
   if ( ipimd==CMD ) then
      nstep_cmd = nstep_cmd + 1
      t_cmd = t_cmd + dt
   end if
   
   !     ---full energies are only calculated every nrespa steps
   !     nvalid is the number of steps where all energies are calculated
   
   if ( onstep )then
      nvalid = nvalid + 1
      enert(1:nren) = enert(1:nren)+ener(1:nren)
      enert2(1:nren) = enert2(1:nren) + ener(1:nren)**2
#ifdef MPI
      if( ievb /= 0 ) then
         evb_nrg_ave(:) = evb_nrg_ave(:) + evb_nrg(:)
         evb_nrg_rms(:) = evb_nrg_rms(:) + evb_nrg(:)**2
      endif
      if ( ifsc /= 0 ) then
         sc_ener_ave(1:8) = sc_ener_ave(1:8) + sc_ener(1:8)
         sc_ener_rms(1:8) = sc_ener_rms(1:8) + sc_ener(1:8)**2
      end if
#endif
      if( nvalid == 1 ) etot_start = ener(1)

#ifndef LES
      if ( ipimd>0 .or. ineb>0 ) then
#  ifdef MPI
         if(master) call mpi_reduce(ener(2),totener(2),1,MPI_DOUBLE_PRECISION, &
                                    MPI_SUM,0,commmaster,ierr)
#  endif 
      endif

      ! Passing of dvdl=dV/dl for TI w.r.t. mass
      ! Note that ener(39) (in runmd and mix_frcti) = 
      !       = ener(17) = ene(21) (in force). All denote dvdl.
      if (ipimd>0 .and. itimass>0) totener(39) = ener(39)

      if(ipimd.eq.NMPIMD.and.ntp>0) then
         totener(14) = pressure * pconv
         totener(42) = tmass / (0.602204d0*volume)
      endif
      if(ipimd.eq.CMD) then
         totener(2)=eke_cmd*0.5d0
         totener(4)=totener(2)
         totener(1) = totener(2) + totener(23)
      endif
      totenert(1:nren) = totenert(1:nren)+totener(1:nren)
      totenert2(1:nren) = totenert2(1:nren) + totener(1:nren)**2
#endif /* LES */

   end if

   ! added for rbornstat
   if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1

   ntnb = 0
   if (mod(nstep,nsnb) == 0) ntnb = 1

   ! Since nstep has been incremented, total_nstep is now equal to
   !    (mdloop-1)*nstlim+nstep for REMD and nstep for MD.
   lout = mod(total_nstep,ntpr) == 0 .and. onstep

   irespa = irespa + 1
     
   ! reset pb-related flags
#ifdef MPI
   if(mytaskid == 0)then
#endif   
      if ( igb == 10 ) then
         if ( mod(nstep,npbgrid) == 0 .and. nstep /= nstlim ) pbgrid = .true.
         if ( mod(nstep,ntpr) == 0 .or. nstep == nstlim ) pbprint = .true.
         if ( mod(nstep,nsnbr) == 0 .and. nstep /= nstlim ) ntnbr = 1
         if ( mod(nstep,nsnba) == 0 .and. nstep /= nstlim ) ntnba = 1
      end if
#ifdef MPI
   endif
#endif   

   !-------------------------------------------------------------------
   !  Step 9:  output from this step if required:
   !-------------------------------------------------------------------
   
   !     ...only the master needs to do the output

   if (master) then
      
      !        -- restrt:
      
      if (ixdump) then

         if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
             call trans_pos_nmode_to_cart(x,cartpos)
             call trans_vel_nmode_to_cart(v,cartvel)
         end if

         ! NOTE - This assumes that if numextra > 0, then velocities are
         !        found in the array v...
         if (numextra > 0) call zero_extra_pnts_vec(v,ix)

         if( iwrap == 0 ) then
            nr = nrp
#ifdef LES
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    cartpos,cartvel,xx(lcrdr),box,t,temp0)
            else
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    x,v,xx(lcrdr),box,t,temp0les)
            endif
#else
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    cartpos,cartvel,xx(lcrdr),box,t,temp0)
            else
               call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                    x,v,xx(lcrdr),box,t,temp0)
            endif
#endif
         else
            
            ! --- use temp. array to hold coords. so that the master's values
            !     are always identical to those on all other nodes:
            
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1)=cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            end if

            call wrap_molecules(nspm,nsp,r_stack(l_temp))
            if(ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
            nr = nrp
#ifdef LES
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,temp0les)
#else
            call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, &
                  r_stack(l_temp),v,xx(lcrdr),box,t,temp0)
#endif
            call free_stack(l_temp,routine)
         end if  ! ( iwrap == 0 )

         if( igb == 0 .and. induced == 1 .and. indmeth == 3) &
            call wrt_dips(xx(linddip),xx(ldipvel),nr,t,title)

         if (icnstph /= 0) then
            call cnstphwriterestart(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
                  ix(icptrsct), xx(lcpene),xx(lcpcrg))
         end if

      end if  ! (ixdump)
      
      !     -- Coordinate archive:
      ! For formatted writes and replica exchange, write out a header line.

      if (itdump) then
#ifdef MPI
         ! Write out current replica#, exchange#, step#, and mytargettemp
         ! If mdloop==0 this is a normal md run (since REMD never calls corpac
         !  when mdloop==0) and we don't want the REMD header.
         ! total_nstep is set in step 5. 
         if (mdloop>0.and.loutfm) &
            write (MDCRD_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                              total_nstep, mytargettemp
#endif

         if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
             call trans_pos_nmode_to_cart(x,cartpos)
         end if
 
         if( iwrap == 0 ) then
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               call corpac(cartpos,1,nrx,MDCRD_UNIT,loutfm)
            else
               call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
            endif
            if(ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
         else
            call get_stack(l_temp,nr3,routine)
            if(.not. rstack_ok)then
               deallocate(r_stack)
               allocate(r_stack(1:lastrst),stat=alloc_ier)
               call reassign_rstack(routine)
            endif
            REQUIRE(rstack_ok)
            if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
               do iatom=1,natom
               do m=1,3
                  r_stack(l_temp+3*(iatom-1)+m-1) = cartpos(m,iatom)
               end do
               end do
            else
               do m=1,nr3
                  r_stack(l_temp+m-1) = x(m)
               end do
            endif
            
            call wrap_molecules(nspm,nsp,r_stack(l_temp))
            if (ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
            
            call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
            call corpac(box,1,3,MDCRD_UNIT,loutfm)
            call free_stack(l_temp,routine)
         end if
      end if  ! (itdump)

      !     Velocity archive:
      
      if (ivdump) then

         ! NOTE - This assumes that if numextra > 0, then velocities are
         !        found in the array v...
         if (numextra > 0) call zero_extra_pnts_vec(v,ix)

#ifdef MPI
         ! Write out current replica#, exchange#, step#, and mytargettemp
         ! If mdloop==0 this is a normal md run (since REMD never calls corpac
         !  when mdloop==0) and we don't want the REMD header.
         if (mdloop>0.and.loutfm) &
            write (MDVEL_UNIT,'(a,3(1x,i8),1x,f8.3)') "REMD ", repnum, mdloop, &
                              total_nstep, mytargettemp
#endif

          if(ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
             call trans_vel_nmode_to_cart ( v, cartvel )
             call corpac(cartvel,1,nrx,MDVEL_UNIT,loutfm)
          else
             call corpac(v,1,nrx,MDVEL_UNIT,loutfm)
          endif
      end if
      !     Energy archive:
      !      (total_nstep set in Step 5.)
      if (ntwe > 0) then
!Antonios modified:
         if (mod(total_nstep,ntwe) == 0.and.onstep) then 
               call mdeng(15,nstep,t,ener,onefac,ntp)
              write(37,*)EHBA
              write(36,*)EHBV
              write(38,*)ECHI
              write(39,*)ENPC
              write(40,*)ENCC
              ener(42) = ECHI
!Antonios end
         end if
      end if
      if (ioutfm > 0) then
         if (itdump) call end_binary_frame(MDCRD_UNIT)
         if (ivdump .and. ntwv>0 ) call end_binary_frame(MDVEL_UNIT)
      end if

#ifdef MPI
      if( ievb /= 0 ) call out_evb ( nstep )
#endif
      
      !     General printed output:
      
      if (lout) then
         if (facc /= 'A') rewind(7)

#ifdef LES
         ! Conserved quantity for Nose'-Hoover thermostat.
         if (ipimd>0.and.ntt==4) then
            Econserved = ener(2) + ener(23) + E_nhc
            Econserved = Econserved + Epot_spring
            if( ntp>0 ) Econserved = Econserved + pres0 / pconv * volume
            write(file_nhc,'(I10,F14.4)') nstep, Econserved
         endif
         if ( ipimd.eq.CMD ) then
            ener(2) = eke_cmd*0.5d0
            ener(4) = ener(2)
            ener(1) = ener(2) + ener(23)
         end if
#else
         if ( ipimd>0 .or. ineb>0 ) then
            if (ipimd>0) then
               ener(1) = 0.d0
               ener(2) = 0.d0
            endif
            ! Conserved quantity for Nose'-Hoover thermostat.
            if ( ipimd>0 .and. ntt==4 ) then
               Econserved = totener(2) + totener(23) + E_nhc
               Econserved = Econserved + Epot_spring
               if ( ntp>0 ) Econserved=Econserved+pres0/pconv*volume
#  ifdef MPI
               if ( worldrank.eq.0 ) &
#  endif
                  write(file_nhc,'(I10,F14.4)') nstep, Econserved
            endif
#  ifdef MPI
            if(worldrank.eq.0) &
#  endif
               call pimd_report(nstep,t,pimd_unit,totener,onefac)
         end if
#endif /* LES */
         call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)

#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener)
         if (ifsc /= 0) call sc_print_energies(7, sc_ener)
#endif

         ! Output for CMD.
#ifdef LES
         if (ipimd.eq.CMD) then

            ncmd = 0
            do iatom = 1, natom
               if ( cnum(iatom)==0 .or. cnum(iatom)==1 )  then
                  xcmd(ncmd+1) = x(3*iatom-2)
                  xcmd(ncmd+2) = x(3*iatom-1)
                  xcmd(ncmd+3) = x(3*iatom)
                  vcmd(ncmd+1) = v(3*iatom-2)
                  vcmd(ncmd+2) = v(3*iatom-1)
                  vcmd(ncmd+3) = v(3*iatom)
                  ncmd = ncmd+3
               endif
            enddo
            write(file_pos_cmd,'(10f8.3)') xcmd(1:ncmd)
            write(file_vel_cmd,'(10f8.3)') vcmd(1:ncmd)
            write(file_pos_cmd,'(10f8.3)') box(1:3)

            eke_cmd = eke_cmd * 0.5d0
            etot_cmd = eke_cmd + ener(23)

            if (eq_cmd) then
               temp_cmd = eke_cmd/boltz2/dble(3*natomCL)
            else
               temp_cmd = eke_cmd/boltz2/dble(3*(natomCL-1))
            endif

         endif
#else
         if (ipimd.eq.CMD.and.mybeadid.eq.1) then
            write(file_pos_cmd,'(10f8.3)') x(1:3*natom)
            write(file_vel_cmd,'(10f8.3)') v(1:3*natom)
            write(file_pos_cmd,'(10f8.3)') box(1:3)

            eke_cmd = eke_cmd * 0.5d0
            etot_cmd = eke_cmd + totener(23)

            if (eq_cmd) then
               temp_cmd = eke_cmd/boltz2/dble(3*natom)
            else
               temp_cmd = eke_cmd/boltz2/dble(3*(natom-1))
            endif
         end if
#endif /* LES */

          !--- Print QMMM Muliken Charges if needed ---
          if (qmmm_nml%ifqnt) then
            if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
              total_mulliken_charge=0.0d0
              total_cm3_chg=0.d0
              write(6,'("    Atomic Charges for Step",i8," :")') nstep
              if (qmmm_nml%dftb_chg == 1) then
                write(6,'("  Atom    Element       Mulliken Charge       CM3 Charge")')
                call qm2_dftb_cm3(qm2_struct%scf_mchg)
              else
                 write(6,'("  Atom    Element       Mulliken Charge")')
              end if
              do i=1,qmmm_struct%nquant_nlink
                 !Mulliken charges have already been calculated and stored.
                 mulliken_charge = qm2_struct%scf_mchg(i)
                 total_mulliken_charge=total_mulliken_charge+mulliken_charge
                 if (qmmm_nml%dftb_chg == 1) then
                   cm3_chg = cm3%qcm3(i)
                   total_cm3_chg =  total_cm3_chg + cm3_chg
                   write(6,'(" ",i5,"      ",A2,"        ",F14.3,"    ",F14.3)') i, &
                   element_sym(qmmm_struct%iqm_atomic_numbers(i)), &
                   mulliken_charge, cm3_chg
                 else
                   write(6,'(" ",i5,"      ",A2,"        ",F14.3)') i, &
                   element_sym(qmmm_struct%iqm_atomic_numbers(i)), &
                   mulliken_charge
                 endif
              end do
              if (qmmm_nml%dftb_chg == 1) then
                write(6,'(" Total Charges: ",8X,F12.3,6X,F12.3)')  &
                     total_mulliken_charge, total_cm3_chg
              else
                write(6,'(" Total Mulliken Charge =",F12.3)') &
                     total_mulliken_charge
              endif
            end if
          end if
         
         !--- BEGIN DIPOLE PRINTING CODE ---

         ! RCW 2nd Dec 2003 - also output dipole information if
         ! the dipoles namelist has been specified and corresponding
         ! groups defined.

         ! Check input unit 5 for namelist dipoles
         ! We expect to find &dipoles followed by a group
         ! specification of the dipoles to output.
         call nmlsrc('dipoles',5,prndipfind)

        if(prndipfind /= 0 ) then
           !We calculate the dipoles
           write(6,*) '------------------------------- DIPOLE INFO ----------------------------------'
           write(6,9018) nstep,t
           9018 format(/1x, 'NSTEP =',i7,1x,'TIME(PS) =',f10.3)

           !Get the groups for the dipoles - Ideally we only really want
           !to call this the once but for the time being I will call it
           !every time
         
           read (5,'(a)') prndiptest

           call rgroup(natom,natc,nres,prndipngrp,ix(i02),ih(m02), &
                ih(m04),ih(m06),ih(m08),ix(icnstrgp), &
                jgroup,indx,irespw,npdec, &
                xx(l60),xx(lcrdr),0,0,0,idecomp,5,.false.)

          ! Need to rewind input file after rgroup so it is available
          ! when we next loop through

           rewind(5)

           if(prndipngrp > 0) then
             !prndipngrp - holds number of groups specified + 1
             !ix(icnstrgp) - holds map of group membership for each atom
             !x(lcrd) - X,Y,Z coords of atoms - (3,*)
             !x(l15) - Partial Charges
             !x(linddip) - induced dipoles X,Y,Z for each atom (3,*)
             !x(Lmass) - Mass of each atom
             call printdip(prndipngrp,ix(icnstrgp),xx(lcrd), &
                  xx(l15),xx(linddip),xx(Lmass), natom)
           end if
           write(6,*) '----------------------------- END DIPOLE INFO --------------------------------'
        end if
        !--- END DIPOLE PRINTING CODE ---

         if (nmropt > 0) then
            call nmrptx(6)
            call nmrptx(7)
         end if
         call amflsh(7)
      end if
      
      !  Output running averages:
      ! DAN ROE: total_nstep==Total nstep REMD/MD, set in step 5
      if ( ntave > 0 )then
         if ( mod(total_nstep,ntave) == 0 .and. onstep )then
            write(6,542)
            tspan = ntave/nrespa
            do m = 1,nren
               enert_tmp(m) = enert(m)-enert_old(m)
               enert2_tmp(m) = enert2(m)-enert2_old(m)
               enert_old(m) = enert(m)
               enert2_old(m) = enert2(m)
               enert_tmp(m) = enert_tmp(m)/tspan
               enert2_tmp(m) = enert2_tmp(m)/tspan - &
                     enert_tmp(m)*enert_tmp(m)
               if ( enert2_tmp(m) < 0.d0 )enert2_tmp(m) = 0.d0
               enert2_tmp(m) = sqrt(enert2_tmp(m))
            end do
#ifdef MPI
            if( ievb /= 0 ) then
               evb_nrg_tmp (:) = evb_nrg_ave(:) - evb_nrg_old (:)
               evb_nrg_tmp2(:) = evb_nrg_rms(:) - evb_nrg_old2(:)
               evb_nrg_old (:) = evb_nrg_ave(:)
               evb_nrg_old2(:) = evb_nrg_rms(:)
               evb_nrg_tmp (:) = evb_nrg_tmp (:) / tspan
               evb_nrg_tmp2(:) = evb_nrg_tmp2(:) / tspan - evb_nrg_tmp(:)**2
               evb_nrg_tmp2(:) = max( evb_nrg_tmp2(:), 0.0d0 )
               evb_nrg_tmp2(:) = sqrt( evb_nrg_tmp2(:) )
            endif
            if ( ifsc /= 0 ) then
               do m = 1,8
                  sc_ener_tmp(m) = sc_ener_ave(m)-sc_ener_old(m)
                  sc_ener_tmp2(m) = sc_ener_rms(m)-sc_ener_old2(m)
                  sc_ener_old(m) = sc_ener_ave(m)
                  sc_ener_old2(m) = sc_ener_rms(m)
                  sc_ener_tmp(m) = sc_ener_tmp(m)/tspan
                  sc_ener_tmp2(m) = sc_ener_tmp2(m)/tspan - sc_ener_tmp(m)**2
                  if (sc_ener_tmp2(m) < 0.0d0) sc_ener_tmp2(m) = 0.0d0
                  sc_ener_tmp2(m) = sqrt(sc_ener_tmp2(m))
               end do
            end if         
            if( ievb /= 0 ) evb_frc%evb_ave = .true.
#endif
            write(6,540) ntave/nrespa
            call prntmd(nstep,izero,izero,t,enert_tmp,onefac,0,.false.)
#ifdef MPI
            if (ifsc /= 0) call sc_print_energies(6, sc_ener_tmp)
            if( ievb /= 0 ) evb_frc%evb_rms = .true.
#endif
            write(6,550)
            call prntmd(nstep,izero,izero,t,enert2_tmp,onefac,0,.true.)
#ifdef MPI /* SOFT CORE */
            if (ifsc /= 0) call sc_print_energies(6, sc_ener_tmp2)
#endif
            if( icfe > 0 ) then
               write(6,541) ntave/nrespa
               edvdl_r(1:51) = edvdl_r(1:51)/tspan
               edvdl_r(39) = enert_tmp(39)  ! fix for DV/DL output
               call prntmd(nstep,izero,izero,t,edvdl_r,onefac,0,.false.)
               edvdl_r(1:51) = 0.d0
            end if
            write(6,542)
         end if
      end if  ! ( ntave > 0 )
      
      !     --- end masters output ---
      
   end if  ! (master)

#ifdef MPI /* SOFT CORE */
   if (ntave > 0 .and. icfe > 0 .and. dynlmb > 0) then
      if ( mod(nstep,ntave) == 0 .and. onstep ) then
         ! For runs with dynamically changing lambda, raise lambda here
         ! and flush all buffers for the next averages
         clambda = clambda + dynlmb
         call sc_change_clambda(clambda)
         if (master) then
            sc_ener(1:8) = 0.0d0
            sc_ener_ave(1:8) = 0.0d0
            sc_ener_rms(1:8) = 0.0d0
            sc_ener_old(1:8) = 0.0d0
            sc_ener_old2(1:8) = 0.0d0
            enert = 0.0d0
            enert2 = 0.0d0
            enert_old = 0.d0
            enert2_old = 0.d0
            write (6,*)
            write (6,'(a,f12.4,a,f12.4)') &
                'Dynamically changing lambda: Increased clambda by ', &
                dynlmb, ' to ', clambda
            write (6,*)
         end if
      end if
   end if
#endif
   
   !=======================================================================
   
   !  ---major cycle back to new step unless we have reached our limit:
   

#ifdef MMTSB
   if ( mmtsb_switch /= mmtsb_off ) then
      if ( mod( nstep, mmtsb_iterations ) == 0 ) then
         write(6,'(a,i8)') &
               'MMTSB Replica Exchange iterations completed at NSTEP = ', &
               nstep
         ! apparently 23 is the magic number for potential energy.
         write(6,'(a,f12.4)') &
               'MMTSB Replica Exchange potential energy = ', ener(23)
         ! write coordinates; preferred format is pdb, but can't do that
         ! so write a restart file; server will post process with ambpdb.
#  ifdef LES
         call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, x,v, &
               xx(lcrdr),box,t,temp0les)
#  else
         call mdwrit(nstep,nrp,nr,nres,ntxo,ntr,ntb, x,v, &
               xx(lcrdr),box,t,temp0)
#  endif
         ! Chop up trajectory files for later continuous temp splicing.
         call close_dump_files
         ! contact server
         if ( mmtsb_switch == mmtsb_temp_rex ) then
            call mmtsb_newtemp( ener(23), temp_mmtsb, is_done_mmtsb )
         else if ( mmtsb_switch == mmtsb_lambda_rex ) then
            ! currently temp_mmtsb is ignored, but multidimensional soon
            call mmtsb_newlambda( unpert_pe_mmtsb, pert_pe_mmtsb, &
                  lambda_mmtsb, temp_mmtsb, is_done_mmtsb )
         end if
         call open_dump_files
         if ( is_done_mmtsb ) then
            goto 480
         end if

         ! in the future we may want amber based tracking of exchanges
         ! perhaps we can use the Simmerling group's code ?
         if ( mmtsb_switch == mmtsb_temp_rex ) then
            if ( abs( temp_mmtsb - temp0 ) <= TEN_TO_MINUS3 ) then
               ! no exchange, continue at the same reference temp.
               mmtsb_is_exchanged = .false.
               write(6,'(a,i8,a,f12.4)') &
                     'MMTSB Replica Exchange temperature unchanged'
            else
               ! exchange temp via changing the reference temp.
               ! the velocities will be randomly reset at the new temp via
               ! the resetvelo variable.
               mmtsb_is_exchanged = .true.
               write(6,'(a,f8.2,a,f8.2)') &
                     'MMTSB Replica Exchange temperature change from ', &
                     temp0, ' to ', temp_mmtsb
               temp0 = temp_mmtsb
            end if
         else if ( mmtsb_switch == mmtsb_lambda_rex ) then
            if ( abs( lambda_mmtsb - clambda ) <= TEN_TO_MINUS3 ) then
               ! no exchange, continue at the same lambda
               mmtsb_is_exchanged = .false.
               write(6,'(a,i8,a,f12.4)') &
                     'MMTSB Replica Exchange lambda unchanged'
            else
               ! exchange lambda
               ! the velocities will be randomly reset via
               ! the resetvelo variable.
               mmtsb_is_exchanged = .true.
               write(6,'(a,f8.2,a,f8.2)') &
                     'MMTSB Replica Exchange lambda change from ', &
                     clambda, ' to ', lambda_mmtsb
               clambda = lambda_mmtsb
            end if
         end if  ! ( mmtsb_switch == mmtsb_temp_rex )
      else
         ! not a replica exchange update iteration.
         mmtsb_is_exchanged = .false.
      end if  ! ( mod( nstep, mmtsb_iterations ) == 0 )
   end if  ! ( mmtsb_switch /= mmtsb_off )
#endif


   call trace_integer( 'end of step', nstep )
   call trace_output_mpi_tally( )
   call timer_stop(TIME_VERLET)
#if !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD)
   call ncsu_on_mdstep(ener(23), v, ekmh)
#endif /* !defined(DISABLE_NCSU) && defined(NCSU_ENABLE_BBMD) */
   if (nstep < nstlim) goto 260
   480 continue


#ifdef MPI     
! ------====== REMD Post-Dynamics ======------
   if(rem == 1) then
      remd_ekmh=ekmh
               
      ! ---=== HYBRID REMD ===---
      if (numwatkeep>=0) then 
         ! This is a hybrid REMD run. Get energy of stripped system for next
         !  exchange.
         call hybrid_remd_ene(xx,ix,ih,ipairs,qsetup,      &
                              numwatkeep,hybridgb,igb,nspm,t,temp0, &
                              ntb,cut,                              &
                              ener,vir,do_list_update,nstep,        &
                              nitp,nits,onefac,loutfm )
      endif ! numwatkeep>=0

      ! Set myeptot, mytemp, and mytargettemp
      if (mdloop>0) mytemp = ener(2) * onefac(1)
      myeptot = ener(23)
      mytargettemp = temp0
      if (master) write(6,'(a,f15.4,2(a,f6.2))') &
         "REMD: myEptot= ",myeptot," myTargetTemp= ", &
         mytargettemp," mytemp= ",mytemp
#  ifdef LES
   else if(rem == 2 ) then
      mytemp = ener(4) * onefac(3)
      myeptot = ener(40)
      mytargettemp = temp0les
#  endif
   endif ! rem==1
! ------====== END REMD Post-Dynamics ======------
#endif /* MPI */

   !=======================================================================
   !     ----- PRINT AVERAGES -----
   !=======================================================================
   
#  ifdef MPI
   ! -- ti decomp
   if (icfe /= 0 .and. idecomp /= 0) then
      if( idecomp == 1 .or. idecomp == 2 ) then
         call collect_dec(nrs)
      !else if( idecomp == 3 .or. idecomp == 4 ) then
      !   call collect_dec(npdec*npdec)
      end if
   end if

   ! Turn off avg. for REMD. 
   if (master.and.rem==0) then
#  else
   if (master) then
#  endif /*MPI*/
      tspan = nvalid
      if (nvalid > 0) then
         do m = 1,nren
            enert(m) = enert(m)/tspan
            enert2(m) = enert2(m)/tspan - enert(m)*enert(m)
            if(enert2(m) < 0.d0) enert2(m) = 0.d0
            enert2(m) =  sqrt(enert2(m))
            edvdl(m) = edvdl(m)/tspan

            ! for PIMD/NMPIMD/CMD/RPMD averages
            if (ipimd>0) then
               totenert(m) = totenert(m)/tspan
               totenert2(m) = totenert2(m)/tspan - totenert(m)*totenert(m)
               if(totenert2(m) < 0.d0) totenert2(m) = 0.d0
               totenert2(m) =  sqrt(totenert2(m))
            endif
         end do
#ifdef MPI
         if( ievb /= 0 ) then
            evb_nrg_ave(:) = evb_nrg_ave(:) / tspan
            evb_nrg_rms(:) = evb_nrg_rms(:) / tspan - evb_nrg_ave(:)**2
            evb_nrg_rms(:) = max( evb_nrg_rms(:), 0.0d0 ) 
            evb_nrg_rms(:) = sqrt( evb_nrg_rms(:) ) 
         endif
         if ( ifsc /= 0 ) then
            do m = 1,8
               sc_ener_ave(m) = sc_ener_ave(m)/tspan
               sc_ener_rms(m) = sc_ener_rms(m)/tspan - sc_ener_ave(m)**2
               if(sc_ener_rms(m) < 0.0d0) sc_ener_rms(m) = 0.0d0
               sc_ener_rms(m) = sqrt(sc_ener_rms(m))
            end do
         end if         
         if( ievb /= 0 ) evb_frc%evb_ave = .true.
#endif
         write(6,540) nvalid
         call prntmd(nstep,izero,izero,t,enert,onefac,0,.false.)
#ifdef MPI /* SOFT CORE */
         if (ifsc /= 0) call sc_print_energies(6, sc_ener_ave)
         if( ievb /= 0 ) evb_frc%evb_rms = .true.
         if ( ipimd > 0 .and. worldrank==0 ) then
            write(pimd_unit,540) nvalid
            call pimd_report(nstep,t,pimd_unit,totenert,onefac)
            write(pimd_unit,550)
            call pimd_report(nstep,t,pimd_unit,totenert2,onefac)
         endif
#endif
         if (nmropt > 0) call nmrptx(6)
         write(6,550)
         call prntmd(nstep,izero,izero,t,enert2,onefac,0,.true.)

#ifdef MPI
         if (ifsc /= 0) call sc_print_energies(6, sc_ener_rms)
         if (ifsc /= 0) call sc_print_dvdl_values()
         
         if( icfe > 0 ) then
            write(6,541) nvalid
            edvdl(39) = enert(39)  ! fix for DV/DL output
            call prntmd(nstep,izero,izero,t,edvdl,onefac,0,.false.)
            ! -- ti decomp
            if(worldrank == 0 .and. idecomp /= 0) then
               call checkdec(idecomp)
               if(idecomp == 1 .or. idecomp == 2) call printdec(ix)
            end if
         end if
#endif

         if (nmropt >= 1) then
            write(6,500)
            if (iredir(7) /= 0) call pcshift(-1,x,f)
            call ndvptx(x,f,ih(m04),ih(m02),ix(i02),nres,xx(l95), &
                  natom, xx(lwinv),ntb,xx(lnmr01),ix(inmr02),6)
         end if
         
         ! Print Born radii statistics
         
         if ((rbornstat == 1).and.(igb /= 0)) then

            ! Born radii stats collected every nrespai step not nrespa step
            tspan = nvalidi

            write(6,580) nstep
            write(6,590)
            do m = 1,natom
               xx(l188-1+m) = xx(l188-1+m)/tspan
               xx(l189-1+m) = xx(l189-1+m)/tspan - &
                     xx(l188-1+m)*xx(l188-1+m)
               xx(l189-1+m) = sqrt(xx(l189-1+m))
               write(6,600) m, xx(l186-1+m), xx(l187-1+m), &
                     xx(l188-1+m), xx(l189-1+m)
            end do
         end if
         
         do m = 2,nrek
            enert(m) = enert(m)*onefac(m-1)
            enert2(m) = enert2(m)*onefac(m-1)
         end do
         temp = enert(2)
      end if  ! (nvalid > 0)
      
   end if  ! (master)

#ifdef MPI
   if( ievb /= 0 ) then
      call evb_dealloc
#if defined(LES)
      if( master ) call evb_pimd_dealloc
#endif
   endif
#endif
   
   if( icfe /= 0 ) then
      deallocate( frcti, stat = ier )
      REQUIRE( ier == 0 )
   end if

   500 format(/,' NMR restraints on final step:'/)
   540 format(/5x,' A V E R A G E S   O V E R ',i7,' S T E P S',/)
   541 format(/5x,' DV/DL, AVERAGES OVER ',i7,' STEPS',/)
   542 format('|',79('='))
   550 format(/5x,' R M S  F L U C T U A T I O N S',/)
   580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
   590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
   600 format(i4,2x,4f12.4)
   call trace_exit( 'runmd' )
   return
end subroutine runmd 

subroutine quench(f,v) 
   implicit none
   
#include "md.h" 
!need access to vv - temp verlet scaling
#include "memory.h" 
!need access to natom

   _REAL_ f(*),v(*),dotproduct,force
   !f is the forces and v is the velocity
 
   integer index
   dotproduct = 0.d0
   force = 0.d0
 
   do index=1,3*natom
      force = force + f(index)**2
      dotproduct = dotproduct + v(index)*f(index)
   enddo
  
   if (force/=0.0d0) then 
     force = 1.0d0/sqrt(force)
     dotproduct = dotproduct*force
   end if
   
   if (dotproduct>0.0d0) then
      v(1:3*natom) = dotproduct*f(1:3*natom)*force
   else 
      !v(1:3*natom) = 0.0d0
      v(1:3*natom) = vfac*dotproduct*f(1:3*natom)*force
   end if
   
end subroutine quench

