#include "copyright.h"
#include "dprec.h"
#include "ncsu-config.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open input files and read cntrl namelist.
subroutine mdread1()

#if !defined(DISABLE_NCSU) && defined(MPI)
   use ncsu_sander_hooks, only : ncsu_on_mdread1 => on_mdread1
#endif
   
   use lmod_driver, only : read_lmod_namelist
   use qmmm_module, only : qmmm_nml,qmmm_struct, qm_gb
   use constants, only : RETIRED_INPUT_OPTION, zero, one, two, three, seven, eight
   use amoeba_mdin, only: AMOEBA_read_mdin, iamoeba
#ifdef DSSP
   use dssp, only: idssp
#endif
#ifdef MPI /* SOFT CORE */
   use softcore, only : scalpha, ifsc, scmask, logdvdl, dvdl_norest, dynlmb
#endif

   use nose_hoover_vars, only: nchain
   use lscivr_vars, only: ilscivr, icorf_lsc
   use pimd_vars, only: ipimd,ineb,itimass
   use cmd_vars, only: restart_cmd, eq_cmd, adiab_param
   use stack, only: lastist,lastrst
   use nmr, only: echoin
#ifdef APBS
   use apbs
#endif /* APBS */
   implicit none
#  include "box.h"
#  include "def_time.h"
#  include "ew_cntrl.h"
#  include "files.h"
#  include "md.h"
#  include "memory.h"
#  include "mmtsb.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "sgld.h"
#  include "ew_erfc_spline.h"
#ifdef LES
#  include "les.h"
#else
   _REAL_      temp0les
#endif

   character(len=4) watdef(4),watnam,owtnm,hwtnm1,hwtnm2

   _REAL_      dele
   integer     ierr
   integer     ifind
   integer     imcdo
   integer     itotst
   integer     jn
   logical     mdin_cntrl, mdin_lmod  ! true if these namelists exist in mdin
   integer :: ifqnt    ! local only here - after reading in the namelist gets copied into qmmm_nml%ifqnt
   integer     mxgrp
   integer     isgld
   character(len=8) date
   character(len=10) time
   _REAL_      dtemp  ! retired 
   _REAL_      dxm  ! retired 
   _REAL_      heat  ! retired 
   _REAL_      timlim ! retired
   
   namelist /cntrl/ irest,ibelly, &
         ntx,ntxo,ntcx,ig,tempi, &
         ntb,ntt,nchain,temp0,tautp, &
         ntp,pres0,comp,taup, &
         nscm,nstlim,t,dt, &
         ntc,ntcc,nconp,tol,ntf,ntn,nsnb, &
         cut,scnb,scee,dielc, &
         ntpr,ntwx,ntwv,ntwe,ntave,ntpp,ioutfm, &
         ntr,nrc,ntrx,taur,nmropt, &
         ivcap,cutcap,xcap,ycap,zcap,fcap, &
         xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &
         imin,drms,dele,dx0, &
         pencut,ipnlty,iscale,scalm,noeskp, &
         maxcyc,ncyc,ntmin,vlimit, &
         mxsub,ipol,jfastw,watnam,owtnm,hwtnm1,hwtnm2, iesp, &
         skmin, skmax, vv,vfac, tmode, ips, &
         isgld,isgsta,isgend,tsgavg,sgft,tempsg,&
         jar, iamoeba, &
         numexchg, repcrd, numwatkeep, hybridgb, &
         ntwprt,tausw, &
         ntwr,iyammp,imcdo, &
         igb,alpb,Arad,rgbmax,saltcon,offset,gbsa,vrand, &
         surften,iwrap,nrespa,nrespai,gamma_ln,extdiel,intdiel, &
         cut_inner,icfe,clambda,klambda, rbornstat,lastrst,lastist,  &
         itgtmd,tgtrmsd,tgtmdfrc,tgtfitmask,tgtrmsmask, &
         idecomp,temp0les,restraintmask,restraint_wt,bellymask, &
         noshakemask,crgmask, &
         mmtsb_switch, mmtsb_iterations,rdt,icnstph,solvph,ntcnstph, &
         ifqnt,ievb, ipimd, itimass, ineb,profile_mpi, ilscivr, icorf_lsc, &
#ifdef MPI /* SOFT CORE */
         scalpha, ifsc, scmask, logdvdl, dvdl_norest, dynlmb, &
#endif
#ifdef DSSP
         idssp, &
#endif
         restart_cmd, eq_cmd, adiab_param,  &
         dtemp, dxm, heat, timlim  !all retired 

   ! Define default water residue name and the names of water oxygen & hydrogens
   
   data watdef/'WAT ','O   ','H1  ','H2  '/
   
   !     ----- READ THE CONTROL DATA AND OPEN DIFFERENT FILES -----
   
   if (mdout /= "stdout" ) &
         call amopen(6,mdout,owrite,'F','W')
   call amopen(5,mdin,'O','F','R')
   write(6,9308)
   call date_and_time( DATE=date, TIME=time )
   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)
   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'
   
   ! Echo the file assignments to the user:
   
   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
         'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)  , &
         'RESTRT',restrt(1:70) , 'REFC'  ,refc(1:70)  , &
         'MDVEL' ,mdvel(1:70)  , 'MDEN'   ,mden(1:70) , &
         'MDCRD' ,mdcrd(1:70)  , 'MDINFO' ,mdinfo(1:70), &
         'INPDIP', inpdip(1:70), 'RSTDIP', rstdip(1:70), &
         'INPTRAJ', inptraj(1:70)
   
   ! Echo the input file to the user:
   call echoin(5,6)
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   read(5,'(a80)') title
   !       ----read input in namelist format, first setting up defaults
   
   dtemp = RETIRED_INPUT_OPTION
   dxm   = RETIRED_INPUT_OPTION
   heat  = RETIRED_INPUT_OPTION
   timlim = RETIRED_INPUT_OPTION
   irest = 0
   ibelly = 0
   ipol = 0
   iesp = 0
   ntx = 1
   ntxo = 1
   ig = 71277
   tempi = ZERO
   ntb = 1
   ntt = 0
   nchain = 1
   temp0 = 300.0d0
#ifdef LES
   ! alternate temp for LES copies, if negative then use single bath
   ! single bath not the same as 2 baths with same target T
   temp0les = -ONE
   rdt = ZERO
#endif
   ipimd =0
   itimass = 0   ! Default = no TI w.r.t. mass.
   ineb  =0

   tautp = ONE
   ntp = 0
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntf = 1
   nsnb = 25
   cut =  EIGHT
   scnb = TWO
   scee = 1.2d0
   dielc = ONE
   ntpr = 50
   ntwr = 500
   ntwx = 0
   ntwv = 0
   ntwe = 0
   ntave = 0
   ioutfm = 0
   ntr = 0
   ntrx = 1
   ivcap = 0
   natcap = 0
   fcap = 1.5d0
   cutcap = 0.0d0
   xcap = 0.0d0
   ycap = 0.0d0
   zcap = 0.0d0
   forth = 1.5d0
   xlorth = -1.0d0
   ylorth = -1.0d0
   zlorth = -1.0d0
   xorth = 47114711.0d0
   yorth = 47114711.0d0
   zorth = 47114711.0d0
   numexchg = 0
   repcrd   = 1

   profile_mpi = 0 !whether to write profile_mpi timing file - default = 0 (NO).

   ! number of waters to keep for hybrid model,
   ! numwatkeep: the number of closest
   ! waters to keep. close is defined as close to non-water.
   ! for simulations with ions, ions should be stripped too
   ! or at least ignored in the "closest" calculation. this
   ! is not currently done.

   ! if it stays at -1 then we keep all waters
   ! 0 would mean to strip them all

    numwatkeep=-1

   ! hybridgb: gb model to use with hybrid REMD.
   hybridgb=0
 
   ! carlos targeted MD, like ntr
   
   itgtmd=0
   tgtrmsd=0.
   tgtmdfrc=0.
   tgtfitmask=''
   tgtrmsmask=''

   pencut = 0.1d0
   taumet = 0.0001d0
   omega = 500.0d0
   ipnlty = 1
   scalm = 100.0d0
   iscale = 0
   noeskp = 1
   nmropt = 0
   jar = 0
   tausw = 0.1d0
   imin = 0
   isftrp = 0
   rwell = ONE
   maxcyc = 1
   ncyc = 10
   ntmin = 1
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   jfastw = 0
   watnam = '    '
   owtnm =  '    '
   hwtnm1 = '    '
   hwtnm2 = '    '
   ntwprt = 0
   igb = 0
   alpb = 0
   Arad = 15.0d0
   rgbmax = 25.d0
   saltcon = ZERO
   offset = 0.09d0
   iyammp = 0
   imcdo = -1
   gbsa = 0
   vrand=1000
   surften = 0.005d0
   iwrap = 0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   extdiel = 78.5d0
   intdiel = ONE
   gbgamma = ZERO
   gbbeta = ZERO
   gbalpha = ONE
   gbneckscale = ONE / THREE
   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   ievb = 0
   rbornstat = 0
   idecomp = 0
   lastrst = 1
   lastist = 1
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   noshakemask=''
   crgmask=''
   mmtsb_switch = mmtsb_off ! MMTSB Replica Exchange Off by Default
   mmtsb_iterations = 100   ! MMTSB Replica Exchange Frequency in Iterations

   icnstph = 0
   solvph = SEVEN
   ntcnstph = 10
   skmin = 50 !used by neb calculation
   skmax = 100 !used by neb calculation
   vv = 0 !velocity verlet -- off if vv/=1
   vfac = 0 !velocity verlet scaling factor, 0 by default
   tmode = 1 !default tangent mode for NEB calculation

   ifqnt = 0  ! no QMMM
   ips = 0    ! no isotropic periodic sum

   isgld = 0   ! no self-guiding
   isgsta=1    ! Begining index of SGLD range
   isgend=0    ! Ending index of SGLD range
   tsgavg=0.2d0    !  Local averaging time of SGLD simulation
   sgft=0.0d0      !  Guiding factor of SGLD simulation
   tempsg=1.0d0    !  Guiding temperature of SGLD simulation

   !     Check to see if "cntrl" namelist has been defined.
   mdin_cntrl=.false.
   mdin_ewald=.false.
   mdin_pb=.false.
#ifdef APBS
   mdin_apbs = .false.
#endif /* APBS */
   mdin_lmod=.false.
   mdin_amoeba=.false.
   iamoeba = 0
#ifdef MPI /* SOFT CORE */
   scalpha=0.5
   ifsc=0
   logdvdl=0
   dvdl_norest=0
   dynlmb=0
#endif
#ifdef DSSP
   idssp = 0
#endif

   call nmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.

   call nmlsrc('ewald',5,ifind)
   if (ifind /= 0) mdin_ewald=.true.

   call nmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.

#ifdef APBS
   call nmlsrc('apbs',5,ifind)
   if (ifind /= 0) mdin_apbs=.true.
#endif /* APBS */

   call nmlsrc('lmod',5,ifind)
   if (ifind /= 0) mdin_lmod=.true.

   call nmlsrc('amoeba',5,ifind)
   if (ifind /= 0) mdin_amoeba=.true.

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl,err=999)
   else
      write(6, '(1x,a,/)') 'Could not find cntrl namelist'
      call mexit(6,1)
   end if

   if (ifqnt>0) then
     qmmm_nml%ifqnt = .true.
     if (saltcon /= 0.0d0) then
       qm_gb%saltcon_on = .true.
     else
       qm_gb%saltcon_on = .false.
     end if
     if (alpb == 1) then
       qm_gb%alpb_on = .true.
     else
       qm_gb%alpb_on = .false.
     end if
     if (igb == 10) then
       write(6, '(1x,a,/)') 'QMMM is not compatible with Poisson Boltzmann (igb=10).'
       call mexit(6,1)
     end if
   else
     qmmm_nml%ifqnt = .false.
   end if

   if ( mdin_lmod ) then
      rewind 5
      call read_lmod_namelist()
   end if
   
   !--------------------------------------------------------------------
   !     --- vars have been read ---
   !--------------------------------------------------------------------
   
   write(6,9309)
   
   ! emit warnings for retired cntrl namelist variables

   if ( dtemp /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: dtemp has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( dxm /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: dxm has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
            ! '  The step length will be unlimited.'
   end if
   if ( heat /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: heat has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   
   if ( timlim /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: timlim has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if

   call printflags()

   !--------------------------------------------------------------------
   ! If user has requested ewald electrostatics, read some more input
   !--------------------------------------------------------------------
   
   if( igb == 0 ) call load_ewald_info(parm,inpcrd,ntp,ipol)

   !--------------------------------------------------------------------
   ! parameters for IPS and for SGLD:
   !--------------------------------------------------------------------
   
   tvips = .false.
   teips = .false.
   if( ips == 1 .or. ips == 2 ) then
      use_pme = 0
      eedmeth = 6
      teips = .true.
   end if
   if( ips == 1 .or. ips == 3 ) then
      vdwmeth = 2
      tvips = .true.
   end if
   tlangv=gamma_ln > 0.0d0
   tsgld= isgld > 0
   
   ishake = 0
   if (ntc > 1) ishake = 1
   
   !--------------------------------------------------------------------
   ! Set up some parameters for GB simulations:
   !--------------------------------------------------------------------
   
   if( igb == 2 .or. hybridgb == 2 ) then
      
      !       --- use our best guesses for Onufriev/Case GB  (GB^OBC I)
      
      gbgamma = 2.90912499999d0  ! (the "99999" to force roundoff on print)
      gbbeta = ZERO
      gbalpha = 0.8d0
   end if

   if( igb == 5 .or. hybridgb == 5 ) then
      
      !       --- use our second best guesses for Onufriev/Case GB (GB^OBC II)
      
      gbgamma = 4.851d0
      gbbeta = 0.8d0
      gbalpha = ONE
   end if
   
   if( igb == 7 ) then
      
      !       --- use parameters for Mongan et al. CFA GBNECK
      
      gbgamma = 2.50798245d0
      gbbeta = 1.90792938d0
      gbalpha = 1.09511284d0
      gbneckscale = 0.361825d0
   end if
   
   !--------------------------------------------------------------------
   ! If user has requested PB electrostatics, read some more input
   !--------------------------------------------------------------------

   if ( igb == 10 ) then
      call pb_read
   end if

#ifdef APBS
   if ( mdin_apbs ) then
      call apbs_read
   end if
#endif /* APBS */

   if( iamoeba == 1 ) then
      if( mdin_amoeba ) then
         call AMOEBA_read_mdin(5)
      else
        write(6,*) ' iamoeba is set but the &amoeba namelist was not found'
        call mexit(6,1)
      end if
   end if
   
   ! -------------------------------------------------------------------
   ! If the user has requested NMR restraints, do a cursory read of the
   ! restraints file(s) now to determine the amount of memory necessary
   ! for these restraints:
   ! -------------------------------------------------------------------
   
   if (jar == 1 ) nmropt = 1
   intreq = 0
   irlreq = 0
   if (nmropt > 0) then
      mxgrp = 0
      itotst = 1
      
      ! Set ITOTST to 0 if IMIN equals 1 (i.e. if minimization, not dynamics)
      ! This will cause any "time-averaged" requests to be over-ridden.
      
      if (imin == 1) then 
         itotst = 0
      end if 
      !         CALL AMOPEN(31,NMR,'O','F','R')
      call restlx(5,itotst,mxgrp,dt,6,ierr)
      !         CLOSE(31)
   end if
   
   ! Set the definition of the water molecule. The default definition is in
   ! WATDEF(4).
   
   read(watdef(1),'(A4)') iwtnm
   read(watdef(2),'(A4)') iowtnm
   read(watdef(3),'(A4)') ihwtnm(1)
   read(watdef(4),'(A4)') ihwtnm(2)
   if (watnam /= '    ') read(watnam,'(A4)') iwtnm
   if (owtnm /= '    ') read(owtnm, '(A4)') iowtnm
   if (hwtnm1 /= '    ') read(hwtnm1,'(A4)') ihwtnm(1)
   if (hwtnm2 /= '    ') read(hwtnm2,'(A4)') ihwtnm(2)

#if !defined(DISABLE_NCSU) && defined(MPI)
   call ncsu_on_mdread1()
#endif
   
   return
   999 continue   ! bad cntrl read
   write(6,*) 'error in reading namelist cntrl'
   call mexit(6,1)
   
   ! --- input file polar opts read err trapping:
   
   9308 format(/10x,55('-'),/10x, &
         'Amber 10 SANDER                              2008', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,12('|',a6,': ',a,/))
end subroutine mdread1 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read constant pH file and initialize.
subroutine cnstphread(stateinf,resstate,protcnt,trescnt,&
      statene,chrgdat,charge)
   use constants, only : AMBER_ELECTROSTATIC
   implicit none
#  include "dynph.h"
#  include "files.h"
#  include "random.h"

   type (const_ph_info) :: stateinf(0:TITR_RES_C-1)
   integer :: resstate(0:TITR_RES_C-1), protcnt(0:TITR_STATES_C-1)
   integer, intent(out) :: trescnt
   _REAL_, intent(out) ::  statene(0:TITR_STATES_C-1), chrgdat(0:ATOM_CHRG_C-1),charge(1:*)
   integer itres, iselres, istat, iatom
   integer icumstat, icumchrg
   character(len=40) :: resname(0:TITR_RES_C-1)
   
   common /cnstphresname/ resname
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene,trescnt,resname

   icumstat = -1
   icumchrg = 0

   write(6,'(a,a)') 'reading charge increments from file: ',cpin
   call amopen(CNSTPH_UNIT,cpin,'O','F','R')
   read(18, nml=cnstph)
   do iatom=0, ATOM_CHRG_C-1
      chrgdat(iatom) = chrgdat(iatom) * AMBER_ELECTROSTATIC
   end do
   close(CNSTPH_UNIT)

   !     Alter charges to match specified initial states
   do itres=0,trescnt-1
      do iatom=0, stateinf(itres)%num_atoms - 1 ! For each atom in selected residue
         charge(iatom+stateinf(itres)%first_atom) & !set atom charge to
               = chrgdat(stateinf(itres)%first_charge + iatom & !corresponding atom charge value
               + resstate(itres) * stateinf(itres)%num_atoms & !from selected new state
               )
      end do
   end do
   !     Error (overflow) checking
   if (trescnt > TITR_RES_C) then
      write(6,*) 'Too many titrating residues', &
            'alter dynph.h, recompile'
      stop
   end if
   do itres=0,trescnt-1         !Find res with reference to highest # state
      if (stateinf(itres)%first_state > icumstat) then
         icumstat = stateinf(itres)%first_state
         iselres = itres
      end if
   end do
   icumstat = stateinf(iselres)%first_state + stateinf(iselres)%num_states
   icumchrg = stateinf(iselres)%first_charge + &
         stateinf(iselres)%num_atoms * stateinf(iselres)%num_states
   
   if (icumstat > TITR_STATES_C) then
      write(6,*) 'Too many titrating states', &
            ' alter dynph.h, recompile'
      stop
   end if
   if (icumchrg > ATOM_CHRG_C) then
      write(6,*) 'Too much charge data', &
            'alter dynph.h, recompile'
      stop
   end if
end subroutine cnstphread 



!======================================================================
!          MDREAD2
!======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to defaults and print the inputable variables.
subroutine mdread2(x,ix,ih,ipairs)

   use lmod_driver, only : LMOD_NTMIN_LMOD, LMOD_NTMIN_XMIN, write_lmod_namelist
   use decomp, only : jgroup, indx, irespw
   use findmask
   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_rij_eqns, qm_gb, &
                           PM3, AM1, MNDO, PDDGPM3, PDDGMNDO, PM3CARB1, DFTB, RM1
   use constants, only : ZERO, ONE, TWO
   use parms, only: req,fmn
   use nblist, only: a,b,c,alpha,beta,gamma,nbflag,skinnb,sphere,nbtell, &
                     cutoffnb
   use amoeba_mdin, only : iamoeba,beeman_integrator
   use amoeba_runmd, only : AM_RUNMD_get_coords
#ifdef APBS
   use apbs
#endif /* APBS */

   use nose_hoover_vars, only: nchain
   use pimd_vars, only: ipimd, itimass, ineb
   use cmd_vars, only: restart_cmd, eq_cmd, adiab_param
#if defined(LES) && defined(MPI)
   use evb_pimd, only: evb_pimd_init
#endif
#ifdef MPI /* SOFT CORE */
   use softcore, only : ifsc, scmask, scalpha, dvdl_norest
! REMD
   use remd, only : rem, rremd
#endif
   implicit none
   _REAL_ x(*)
#  include "memory.h"
   integer ix(lasti),ipairs(*)
   character(len=4) ih(*)
   integer nbond
   integer atom1,atom2
   integer ntmp
   logical belly,konst
   character(len=1) atsymb,atsymb2
   character(len=2) atype
   integer ngrp,inerr,nr,iaci,ir,i,mxresat,j
   integer noshakegp( natom ), natnos
   integer crggp( natom )
   _REAL_ dummy,rvdw,dcharge
   logical newstyle
#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#ifdef MPI_DOUBLE_PRECISION
#  undef MPI_DOUBLE_PRECISION
#endif
#  include "mpif.h"
#  include "parallel.h"
   integer ist(MPI_STATUS_SIZE), partner, ierr, nbonh_c, num_noshake_c
#ifdef CRAY_PVP
#  define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   !     ========================= END AMBER/MPI =========================
#endif
#  include "files.h"
#  include "md.h"
#  include "box.h"
#  include "mmtsb.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "ew_mpole.h"
#  include "ew_legal.h"
#  include "def_time.h"
#  include "tgtmd.h"
#  include "sgld.h"
#ifdef LES
#  include "les.h"
#endif
   
   ! -------------------------------------------------------------------
   !      --- set up resat array, containing string identifying
   !          residue for each atom
   ! -------------------------------------------------------------------
   
   mxresat = min( natom, matom )
   ir = 0
   do 5 i=1,mxresat
      if (i >= ix(ir+i02)) ir=ir+1
      write(resat(i),'(a4,1x,a4,i4)') ih(m04+i-1), &
            ih(ir+m02-1),ir
      !                                    ---null terminator:
      resat(i)(14:14) = char(0)
   5 continue
   close(unit=8)
   
   ! -------------------------------------------------------------------
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   ! -------------------------------------------------------------------
   
   nrp = natom
   
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
   
   ! For 0< ipimd <= 3, no removal of COM motion
   if (ipimd>0.and.ipimd<=3) then
      nscm = nstlim + 1
      ndfmin = 0
   endif

   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 999999999
      ndfmin = 0
   end if
   if(nscm <= 0) nscm = 999999999
   if(gamma_ln > 0.0d0)ndfmin=0  ! No COM motion removal for LD simulation
   if(ntt == 4)ndfmin=0  ! No COM motion removal for Nose'-Hoover simulation
   init = 3
   if (irest > 0) init = 4
   if (scnb == ZERO ) scnb = TWO
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0
   
   !     ----- RESET THE CAP IF NEEDED -----
   
   ! ivcap == 0: Cap will be in effect if it is in the prmtop file (ifcap = 1)

   if(ivcap == 1) then
      ! Cap will be in effect even if not in prmtop file
      !   requires additional information in sander.in file as in the case of ivcap == 3, 4, or 5
      ifcap = 2
   else if(ivcap == 2) then
      ! Inactivate cap
      ifcap = 0
   else if(ivcap == 3) then
      ! Sphere -> not yet implemented
      ifcap = 3
   else if(ivcap == 4) then
      ! Orthorhombus
      ifcap = 4
   else if(ivcap == 5) then
      ! Shell of waters around solute
      ifcap = 5
   end if

   if( ig==-1 ) call microsec(ig)
   
   ! -------------------------------------------------------------------
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   ! -------------------------------------------------------------------
   
   nr = nrp
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt

   ! Error Checking for REMD
#ifdef MPI
   if (rem>0) then
      ! Make sure that the number of replicas is even so
      !  that they all have partners in the exchange step
      if (mod(numgroup,2).ne.0) then
         write (6,'(a)') "==================================="
         write (6,'(a)') "REMD requires an even # of replicas"
         write (6,'(a)') "==================================="
         call mexit (6,1)
      endif

      write(6,'(/a)') 'Replica exchange'
      write(6,'(5x,4(a,i8))') 'numexchg=',numexchg,', rem=',rem

      ! REPCRD option temporarily disabled
      if (repcrd==0) write(6,'(a)') &
         "REMD WARNING: repcrd temporarily disabled. Only replica &
         &trajectories/output can be written."

      ! Check for correct number of exchanges
      if (numexchg<0) then
         write(6,'(a)') "REMD ERROR: numexchg must be >= 0, "
         call mexit(6,1)
      endif

      ! Hybrid GB
      if (numwatkeep >= 0) then
         write(6,'(5x,4(a,i8))') 'numwatkeep=',numwatkeep,', hybridgb=',hybridgb
         ! Check that user specified GB model for hybrid REMD
         if (hybridgb/=2.and.hybridgb/=5.and.hybridgb/=1) then
            write(6,'(a)') "HYBRID REMD ERROR: hybridgb must be 1, 2, or 5."
            call mexit(6,1)
         endif
      else
         !Check that user did not specify GB model if no hybrid run.
         if (hybridgb/=0) then
            write(6,'(a)') &
            "HYBRID REMD ERROR: numwatkeep must be >= 0 when hybridgb is set."
            call mexit(6,1)
         endif
      endif
      ! RREMD
      if (rremd>0) then
         write(6,'(5x,4(a,i8))') "rremd=",rremd
      endif

      ! ntb>1 not allowed for remd 
      if (ntb>1) then
         write(6,'(a,i1)') "ERROR: ntb > 1 not allowed for rem > 0, ntb=",ntb
         call mexit(6,1)
      endif

#  ifdef LES
      ! DAN ROE: Temporarily disable LES REMD until it is verified with new
      !           REMD code
      if (rem==2) then
         write (6,*) "******* LES REM (rem==2) temporarily disabled. Stop. *******"
         call mexit(6,1)
      endif

      if (rem==2 .and. igb/=1) then
         write (6,*) ' partial REM (rem=2) only works with igb=1'
         call mexit(6,1)
      endif
#  else
      if (rem==2) then
         write(6,*) '*******  For rem ==2, partial REM'
         write(6,*) 'use sander.LES with topology created by addles'
         call mexit(6,1)
      endif
#  endif

   endif ! rem>0

#else /* MPI */
   ! Check if user set numexchg with no MPI
   if (numexchg>0) write(6,'(a)') &
      "WARNING: numexchg > 0 - for REMD run please recompile sander for &
      &parallel runs."

   ! Check if user set numwatkeep or hybridgb with no MPI - not sensible.
   if (numwatkeep>=0) write(6,'(a)') &
      "WARNING: numwatkeep >= 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."

   if (hybridgb>0) write(6,'(a)') &
      "WARNING: hybridgb > 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."
#endif /* MPI */
   ! End error checking for REMD

   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,4(a,i8))') 'iwrap   =',iwrap,', ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   write(6,'(5x,3(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', idecomp =',idecomp,', rbornstat=',rbornstat

   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', igb     =',igb,', nsnb    =',nsnb
   write(6,'(5x,4(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel

   if (( igb /= 0 .and. igb /= 10).or.hybridgb>0) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', gbalpha= ',gbalpha
      write(6,'(5x,3(a,f10.5))') 'gbbeta  =',gbbeta, &
            ', gbgamma =',gbgamma,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
   end if
   
   if( alpb /= 0 ) then
      write(6,'(5x,3(a,f10.5))') 'Arad =', Arad
   end if    

   write(6,'(5x,3(a,f10.5))') 'scnb    =',scnb, &
         ', scee    =',scee
   
   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr

   if( imin /= 0 ) then
      if( ipimd > 0 ) then
         write(6,'(/a)') 'pimd cannot be used in energy minimization'
         stop
      end if

      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( 0, 1, 2 )
         ! no specific output
      case ( LMOD_NTMIN_XMIN, LMOD_NTMIN_LMOD )
         call write_lmod_namelist( )
      case default
         ! invalid ntmin
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         stop
      end select
   else
      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i10))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit
      
      if( ntt == 1 ) then
         write(6,'(/a)') 'Berendsen (weak-coupling) temperature regulation:'
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', tautp   =', tautp
#ifdef LES
         write(6,'(5x,3(a,f10.5))') 'temp0LES   =',temp0les
#endif
      else if( ntt == 2 ) then
         write(6,'(/a)') 'Anderson (strong collision) temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig, ', vrand   =',vrand
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, ', tempi   =',tempi
      else if( ntt == 3 ) then
         write(6,'(/a)') 'Langevin dynamics temperature regulation:'
         write(6,'(5x,4(a,i8))') 'ig      =',ig
         write(6,'(5x,3(a,f10.5))') 'temp0   =',temp0, &
               ', tempi   =',tempi,', gamma_ln=', gamma_ln
      else if( ntt == 4 ) then
         write(6,'(/a)') 'Nose-Hoover chains'
         write(6,'(5x,(a,f10.5))') 'gamma_ln=', gamma_ln
         write(6,'(5x,(a,i8))') 'number of oscillators=', nchain
      end if

      if( ntp /= 0 ) then
         write(6,'(/a)') 'Pressure regulation:'
         write(6,'(5x,4(a,i8))') 'ntp     =',ntp
         write(6,'(5x,3(a,f10.5))') 'pres0   =',pres0, &
               ', comp    =',comp,', taup    =',taup
      end if

   end if

   if( ntc /= 1 ) then
      write(6,'(/a)') 'SHAKE:'
      write(6,'(5x,4(a,i8))') 'ntc     =',ntc,', jfastw  =',jfastw
      write(6,'(5x,3(a,f10.5))') 'tol     =',tol
   end if

   if( ifcap == 1 .or. ifcap == 2 .or. ifcap == 3 ) then
      write(6,'(/a)') 'Water cap:'
      write(6,'(5x,2(a,i8))') 'ivcap   =',ivcap,', natcap  =',natcap
      write(6,'(5x,2(a,f10.5))') 'fcap    =',fcap, ', cutcap  =',cutcap
      write(6,'(5x,3(a,f10.5))') 'xcap    =',xcap, ', ycap    =',ycap,    &
                                 ', zcap    =',zcap
   else if( ifcap == 4 ) then
      write(6,'(/a)') 'Orthorhombus:'
      write(6,'(5x,1(a,i8))')    'ivcap   =',ivcap
      write(6,'(5x,1(a,f10.5))') 'forth   =',forth
      write(6,'(5x,3(a,f10.5))') 'xlorth  =',xlorth,', ylorth  =',ylorth, &
                                 ', zlorth  =',zlorth
      write(6,'(5x,3(a,f10.5))') 'xorth   =',xorth, ', yorth   =',yorth,  &
                                 ', zorth   =',zorth
   else if( ifcap == 5 ) then
      write(6,'(/a)') 'Water shell:'
      write(6,'(5x,(a,i8,a,f10.5))') 'ivcap   =',ivcap,', cutcap  =',cutcap
   endif

   if( nmropt > 0 ) then
      write(6,'(/a)') 'NMR refinement options:'
      write(6,'(5x,4(a,i8))')'iscale  =',iscale,', noeskp  =',noeskp, &
            ', ipnlty  =',ipnlty,', mxsub   =',mxsub
      write(6,'(5x,3(a,f10.5))') 'scalm   =',scalm, &
            ', pencut  =',pencut,', tausw   =',tausw
   end if

   if( numextra > 0 ) then
      write(6,'(/a)') 'Extra-points options:'
      write(6,'(5x,4(a,i8))') 'frameon =',frameon, &
            ', chngmask=',chngmask
   end if

   if( ipol /= 0 ) then
      write(6,'(/a)') 'Polarizable options:'
      write(6,'(5x,4(a,i8))') 'indmeth =',indmeth, &
            ', maxiter =',maxiter,', irstdip =',irstdip, &
            ', scaldip =',scaldip
      write(6,'(5x,3(a,f10.5))') &
            'diptau  =',diptau,', dipmass =',dipmass
   end if

#ifdef MPI /* SOFT CORE */
   if( icfe /= 0 .or. ifsc/=0) then
      write(6,'(/a)') 'Free energy options:'
      write(6,'(5x,4(a,i8,a,i8))') 'icfe =',icfe,',  ifsc =',ifsc      
      write(6,'(5x,4(a,i8))') 'klambda =',klambda
      write(6,'(5x,3(a,f10.5))') 'clambda =',clambda
   end if
#endif

   ! Options for TI w.r.t. mass.
   select case (itimass)
   case (0)     ! Default: no TI wrt. mass.
   case (1,2)   ! 1 = use virial est., 2 = use thermodynamic est.
      write(6,'(/a)') 'Isotope effects (thermodynamic integration w.r.t. mass):'
      write(6,'(5x,4(a,i8))') 'itimass =',itimass      
      write(6,'(5x,3(a,f10.5))') 'clambda =',clambda
      if (icfe /= 0) then 
         write(6,'(/2x,a,i2,a,i2,a)') 'Error: Cannot do TI w.r.t. both potential (icfe =', &
            icfe, ') and mass (itimass =', itimass, ').'
         stop
      endif
      if (ipimd == 0) then 
         write(6,'(/2x,a)') 'Error (IPIMD=0): TI w.r.t. mass requires a PIMD run.'
         stop       
      endif
   case default ! Invalid itimass
      write(6,'(/2x,a,i2,a)') 'Error: Invalid ITIMASS (', itimass, ' ).'
      stop
   end select

!KFW
!   call mpi_bcast ( ievb, 1, MPI_INTEGER, 0, commworld, ierr )
!   call mpi_barrier ( commworld, ierr )

   if( ievb == 1 ) then
!KFW  write(6,'(/a)') 'EVB options:'
!KFW  write(6,'(5x,3(a,f10.5))') 'V11     =',v11,', V22    =', v22, &
!KFW         ', V12     =', v12
!kfw  write(6,'(5x,3(a,f10.5))') 'kevb    =',kevb,', evbt   =', evbt
   end if

   if( itgtmd /= 0 ) then
      write(6,'(/a)') 'Targeted molecular dynamics:'
      write(6,'(5x,3(a,f10.5))') 'tgtrmsd =',tgtrmsd, &
            ', tgtmdfrc=',tgtmdfrc
   end if

   if( ntb > 0 ) then
      write(6,'(/a)') 'Ewald parameters:'
      write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
            ', ew_type =',ew_type,', nbflag  =',nbflag, &
            ', use_pme =',use_pme
      write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
            ', eedmeth =',eedmeth,', netfrc  =',netfrc
      write(6, 9002) a, b, c
      write(6, 9003) alpha, beta, gamma
      write(6, 9004) nfft1, nfft2, nfft3
      write(6, 9006) cutoffnb, dsum_tol
      write(6, 9007) ew_coeff
      write(6, 9005) order
      9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
      9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
      9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
      9005 format (5x,'Interpolation order =',i5)
      9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
      9007 format (5x,'Ewald Coefficient =',f9.5)
   end if

   if( mmtsb_switch /= mmtsb_off ) then
      call mmtsb_print_banner()
      call mmtsb_init( temp0, clambda )
   end if

   if( icnstph /= 0) then
      write(6, '(/a)') 'Constant pH options:'
      write(6, '(5x,a,i8)') 'ntcnstph =', ntcnstph
      write(6, '(5x,a,f10.5)') 'solvph =', solvph
   end if
 
!---- QMMM Options ----

   if( qmmm_nml%ifqnt ) then
      write(6, '(/a)') 'QMMM options:'
      write(6, '(5x,"        ifqnt = True       nquant = ",i8)') &
                 qmmm_struct%nquant
      write(6, '(5x,"         qmgb = ",i8,"  qmcharge = ",i8,"   adjust_q = ",i8)') &
                 qmmm_nml%qmgb, qmmm_nml%qmcharge, qmmm_nml%adjust_q
      write(6, '(5x,"         spin = ",i8,"     qmcut = ",f8.4, "    qmshake = ",i8)') qmmm_nml%spin, &
                 qmmm_nml%qmcut, qmmm_nml%qmshake
      write(6, '(5x,"lnk_atomic_no = ",i8,"   lnk_dis = ",f8.4,"   qmmm_int = ",i8)') &
                 qmmm_nml%lnk_atomic_no,qmmm_nml%lnk_dis, qmmm_nml%qmmm_int
      if ( qmmm_nml%qmtheory == PM3 ) then 
         write(6, '(5x,"     qm_theory =     PM3")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == AM1 ) then
         write(6, '(5x,"     qm_theory =     AM1")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == MNDO ) then
         write(6, '(5x,"     qm_theory =    MNDO")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == PDDGPM3 ) then
         write(6, '(5x,"     qm_theory = PDDGPM3")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == PDDGMNDO ) then
         write(6, '(5x,"     qm_theory =PDDGMNDO")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == PM3CARB1 ) then
         write(6, '(5x,"     qm_theory =PM3CARB1")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == DFTB ) then
         write(6, '(5x,"     qm_theory =    DFTB")',ADVANCE='NO') 
      else if ( qmmm_nml%qmtheory == RM1 ) then
         write(6, '(5x,"     qm_theory =     RM1")',ADVANCE='NO') 
      else
         write(6, '(5x,"     qm_theory = UNKNOWN!")',ADVANCE='NO') 
      end if
      write (6, '(" verbosity = ",i8)') qmmm_nml%verbosity
      if (qmmm_nml%qmqm_analyt) then
        write(6, '(5x,"       qmqmdx = Analytical")')
      else
        write(6, '(5x,"       qmqmdx = Numerical")')
      end if
      if (qmmm_nml%tight_p_conv) then
        write(6, '(5x," tight_p_conv = True (converge density to SCFCRT)")')
      else
        write(6, '(5x," tight_p_conv = False (converge density to 0.05xSqrt[SCFCRT])")')
      end if
      write(6, '(5x,"      scfconv = ",e9.3,"  itrmax = ",i8)') qmmm_nml%scfconv, qmmm_nml%itrmax
      if (qmmm_nml%printcharges) then
        write(6, '(5x," printcharges = True ")',ADVANCE='NO')
      else
        write(6, '(5x," printcharges = False")',ADVANCE='NO')
      end if
      if (qmmm_nml%peptide_corr) then
        write(6, '(5x," peptide_corr = True")')
      else
        write(6, '(5x," peptide_corr = False")')
      end if
      if (qmmm_nml%qmqmrij_incore) then
        write(6, '(4x,"qmqmrij_incore = True ")',ADVANCE='NO')
      else
        write(6, '(4x,"qmqmrij_incore = False")',ADVANCE='NO')
      end if 
      if (qmmm_nml%qmmmrij_incore) then
        write(6, '(4x,"qmmmrij_incore = True ")')
      else
        write(6, '(4x,"qmmmrij_incore = False")')
      end if
      if (qmmm_nml%qmqm_erep_incore) then
        write(6, '(2x,"qmqm_erep_incore = True ")')
      else
        write(6, '(2x,"qmqm_erep_incore = False")')
      end if
      if (qmmm_nml%allow_pseudo_diag) then
        write(6, '(7x,"pseudo_diag = True ")',ADVANCE='NO')
        write(6, '("pseudo_diag_criteria = ",f8.4)') qmmm_nml%pseudo_diag_criteria
      else
        write(6, '(7x,"pseudo_diag = False")')
      end if
      write(6,   '(6x,"diag_routine = ",i8)') qmmm_nml%diag_routine
      !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
      !have put in the namelist and set the value to false.
      if (qmmm_nml%qm_ewald>0) then
        if (qmmm_nml%qm_pme) then
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = True ")') qmmm_nml%qm_ewald
        else
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
        end if
          write(6, '(10x,"  kmaxqx = ",i4," kmaxqy = ",i4," kmaxqz = ",i4," ksqmaxq = ",i4)') &
                  qmmm_nml%kmaxqx, qmmm_nml%kmaxqy, qmmm_nml%kmaxqz, qmmm_nml%ksqmaxq
      else
        write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
      end if
   end if

! ---------------------

#ifdef MPI
! --- MPI TIMING OPTIONS ---
      write(6, '(/a)') '| MPI Timing options:'
      write(6, '("|",5x," profile_mpi = ",i8)') profile_mpi
! Sanity check for profile_mpi
      call int_legal_range('profile_mpi',profile_mpi,0,1)
! --------------------------
#endif
 
   cut = cut*cut
   cut_inner = cut_inner*cut_inner
   
   
   !------------------------------------------------------------------------
   ! If user has requested generalized born electrostatics, set up variables
   !------------------------------------------------------------------------
   
   if( igb == 0 .and. gbsa > 0 ) then
      write(0,*) 'GB/SA calculation is performed only when igb>0'
      call mexit( 6,1 )
   end if
   if( gbsa == 2 .and. &
       ((imin == 0 .and. nstlim > 1) .or. &
        (imin == 1 .and. maxcyc > 1)) ) then
      write(0,*) 'GBSA=2 only works for single point energy calc'
      call mexit( 6,1 )
   end if
#ifdef APBS
   if( igb /= 0 .and. igb /= 10 .and. .not. mdin_apbs) then
#else
   if (( igb /= 0 .and. igb /= 10).or.hybridgb>0) then
#endif /* APBS */
#ifdef LES
      write(6,*) 'igb=1,5,7 are working with LES, no SA term included'
#endif
      ! igb7 uses special S_x screening params.
      ! overwrite the tinker values read from the prmtop
      if (igb == 7) then
         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if (atype(1:1) == 'C' .or. atype(1:1) == 'c') then
               x(l96+i-1) = 4.84353823306d-1
            else if (atype(1:1) == 'H' .or. atype(1:1) == 'h') then
               x(l96+i-1) = 1.09085413633d0
            else if (atype(1:1) == 'N' .or. atype(1:1) == 'n') then
               x(l96+i-1) = 7.00147318409d-1
            else if (atype(1:1) == 'O' .or. atype(1:1) == 'o') then
               x(l96+i-1) = 1.06557401132d0
            else if (atype(1:1) == 'S' .or. atype(1:1) == 's') then
               x(l96+i-1) = 6.02256336067d-1
            else if (atype(1:1) == 'P' .or. atype(1:1) == 'p') then
               x(l96+i-1) = 5d-1
            else
               x(l96+i-1) = 5d-1
            end if
         end do
      end if
      
      !       put fs(i)*(rborn(i) - offset) into the "fs" array
      
      fsmax = 0.d0
      do i=1,natom
         x(l96-1+i) = x(l96-1+i)*( x(l97-1+i) - offset )
         fsmax = max( fsmax, x(l96-1+i) )
         if (rbornstat == 1) then
            x(l186-1+i) = 0.d0
            x(l187-1+i) = 999.d0
            x(l188-1+i) = 0.d0
            x(l189-1+i) = 0.d0
         end if
      end do
      
      !     ---------------------------------------------------------------------
      !       ---get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !         T = 298.15, epsext=78.5,
      
      kappa = sqrt( 0.10806d0 * saltcon )
      
      !       ---scale kappa by 0.73 to account(?) for lack of ion exlcusions:
      
      kappa = 0.73d0* kappa

      !Set kappa for qmmm if needed
      qm_gb%kappa = kappa
      !     ---------------------------------------------------------------------
      
      if ( gbsa == 1 ) then
         
         !     --- assign parameters for calculating SASA according to the
         !         LCPO method ---
         
         do i=1,natom
            ix(i80+i-1)=0
         end do
         
         !         --- get the number of bonded neighbors for each atom:
         
         do i=1,nbona
            atom1=ix(iiba-1+i)/3+1
            atom2=ix(ijba-1+i)/3+1
            ix(i80+atom1-1)=ix(i80+atom1-1)+1
            ix(i80+atom2-1)=ix(i80+atom2-1)+1
         end do
         
         !         --- construct parameters for SA calculation; note that the
         !             radii stored in L165 are augmented by 1.4 Ang.
         
         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            nbond=ix(i80+i-1)
            if (atype == 'CT') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.77887d0
                  x(l175-1+i) = -0.28063d0
                  x(l180-1+i) = -0.0012968d0
                  x(l185-1+i) = 0.00039328d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.56482d0
                  x(l175-1+i) = -0.19608d0
                  x(l180-1+i) = -0.0010219d0
                  x(l185-1+i) = 0.0002658d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.23348d0
                  x(l175-1+i) = -0.072627d0
                  x(l180-1+i) = -0.00020079d0
                  x(l185-1+i) = 0.00007967d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.00000d0
                  x(l175-1+i) = 0.00000d0
                  x(l180-1+i) = 0.00000d0
                  x(l185-1+i) = 0.00000d0
               else
                  write(6,*) 'Unusual nbond for CT:', i, nbond, &
                     ' Using default carbon LCPO parameters'
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.77887d0
                  x(l175-1+i) = -0.28063d0
                  x(l180-1+i) = -0.0012968d0
                  x(l185-1+i) = 0.00039328d0
               end if
            else if (atype(1:1) == 'C' .or. atype(1:1) == 'c') then
               if (nbond == 2) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.51245d0
                  x(l175-1+i) = -0.15966d0
                  x(l180-1+i) = -0.00019781d0
                  x(l185-1+i) = 0.00016392d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.070344d0
                  x(l175-1+i) = -0.019015d0
                  x(l180-1+i) = -0.000022009d0
                  x(l185-1+i) = 0.000016875d0
               else
                  write(6,*) 'Unusual nbond for C :', i, nbond, &
                     ' Using default carbon LCPO parameters'
                  x(l165-1+i) = 1.70d0 + 1.4d0
                  x(l170-1+i) = 0.77887d0
                  x(l175-1+i) = -0.28063d0
                  x(l180-1+i) = -0.0012968d0
                  x(l185-1+i) = 0.00039328d0
               end if
            else if (atype == 'O ') then
               x(l165-1+i) = 1.60d0 + 1.4d0
               x(l170-1+i) = 0.68563d0
               x(l175-1+i) = -0.1868d0
               x(l180-1+i) = -0.00135573d0
               x(l185-1+i) = 0.00023743d0
            else if (atype == 'O2') then
               x(l165-1+i) = 1.60d0 + 1.4d0
               x(l170-1+i) = 0.88857d0
               x(l175-1+i) = -0.33421d0
               x(l180-1+i) = -0.0018683d0
               x(l185-1+i) = 0.00049372d0
            else if (atype(1:1) == 'O' .or. atype(1:1) == 'o') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.77914d0
                  x(l175-1+i) = -0.25262d0
                  x(l180-1+i) = -0.0016056d0
                  x(l185-1+i) = 0.00035071d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.49392d0
                  x(l175-1+i) = -0.16038d0
                  x(l180-1+i) = -0.00015512d0
                  x(l185-1+i) = 0.00016453d0
               else
                  write(6,*) 'Unusual nbond for O:', i, nbond, &
                     ' Using default oxygen LCPO parameters'
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.77914d0
                  x(l175-1+i) = -0.25262d0
                  x(l180-1+i) = -0.0016056d0
                  x(l185-1+i) = 0.00035071d0
               end if
            else if (atype == 'N3') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.078602d0
                  x(l175-1+i) = -0.29198d0
                  x(l180-1+i) = -0.0006537d0
                  x(l185-1+i) = 0.00036247d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.22599d0
                  x(l175-1+i) = -0.036648d0
                  x(l180-1+i) = -0.0012297d0
                  x(l185-1+i) = 0.000080038d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.051481d0
                  x(l175-1+i) = -0.012603d0
                  x(l180-1+i) = -0.00032006d0
                  x(l185-1+i) = 0.000024774d0
               else
                  write(6,*) 'Unusual nbond for N3:', i, nbond, &
                     ' Using default nitrogen LCPO parameters'
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.078602d0
                  x(l175-1+i) = -0.29198d0
                  x(l180-1+i) = -0.0006537d0
                  x(l185-1+i) = 0.00036247d0
               end if
            else if (atype(1:1) == 'N' .or. atype(1:1) == 'n') then
               if (nbond == 1) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.73511d0
                  x(l175-1+i) = -0.22116d0
                  x(l180-1+i) = -0.00089148d0
                  x(l185-1+i) = 0.0002523d0
               else if (nbond == 2) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.41102d0
                  x(l175-1+i) = -0.12254d0
                  x(l180-1+i) = -0.000075448d0
                  x(l185-1+i) = 0.00011804d0
               else if (nbond == 3) then
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.062577d0
                  x(l175-1+i) = -0.017874d0
                  x(l180-1+i) = -0.00008312d0
                  x(l185-1+i) = 0.000019849d0
               else
                  write(6,*) 'Unusual nbond for N:', i, nbond, &
                     ' Using default nitrogen LCPO parameters'
                  x(l165-1+i) = 1.65d0 + 1.4d0
                  x(l170-1+i) = 0.078602d0
                  x(l175-1+i) = -0.29198d0
                  x(l180-1+i) = -0.0006537d0
                  x(l185-1+i) = 0.00036247d0
               end if
            else if (atype == 'SH') then
               x(l165-1+i) = 1.90d0 + 1.4d0
               x(l170-1+i) = 0.7722d0
               x(l175-1+i) = -0.26393d0
               x(l180-1+i) = 0.0010629d0
               x(l185-1+i) = 0.0002179d0
            else if (atype(1:1) == 'S' .or. atype(1:1) == 's') then
               x(l165-1+i) = 1.90d0 + 1.4d0
               x(l170-1+i) = 0.54581d0
               x(l175-1+i) = -0.19477d0
               x(l180-1+i) = -0.0012873d0
               x(l185-1+i) = 0.00029247d0
            else if (atype(1:1) == 'P' .or. atype(1:1) == 'p') then
               if (nbond == 3) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.03873d0
                  x(l175-1+i) = -0.0089339d0
                  x(l180-1+i) = 0.0000083582d0
                  x(l185-1+i) = 0.0000030381d0
               else
                  write(6,*) 'Unusual nbond for P:', i, nbond, &
                     ' Using default phosphorus LCPO parameters'
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               end if
            else if (atype(1:1) == 'Z') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype(1:1) == 'H' .or. atype(1:1) == 'h') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype == 'MG') then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(l165-1+i) = 1.18d0 + 1.4d0
               !  The following values were taken from O.sp3 with two bonded 
               !  neighbors -> O has the smallest van der Waals radius 
               ! compared to all other elements which had been parametrized
               x(l170-1+i) = 0.49392d0
               x(l175-1+i) = -0.16038d0
               x(l180-1+i) = -0.00015512d0
               x(l185-1+i) = 0.00016453d0
            else
               ! write( 0,* ) 'bad atom type: ',atype
               ! call mexit( 6,1 )
               x(l165-1+i) = 1.70 + 1.4;
               x(l170-1+i) = 0.51245;
               x(l175-1+i) = -0.15966;
               x(l180-1+i) = -0.00019781;
               x(l185-1+i) = 0.00016392;
               write(6,'(a,a)') 'Using carbon SA parms for atom type', atype 
            end if
         end do  !  i=1,natom
         !
      else if ( gbsa == 2 ) then

         !     --- assign parameters for calculating SASA according to the
         !         ICOSA method; the radii are augmented by 1.4 A ---

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if (atype(1:1) == 'N' .or. atype(1:1) == 'n') then
               x(L165-1+i) = 1.55d0 + 1.4d0
            else if (atype(1:1) == 'C' .or. atype(1:1) == 'c') then
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atype(1:1) == 'H' .or. atype(1:1) == 'h' .or. &
                     atype == '1H' .or. &
                     atype == '2H' .or. &
                     atype == '3H') then
               x(L165-1+i) = 1.20d0 + 1.4d0
            else if (atype(1:1) == 'O' .or. atype(1:1) == 'o') then
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atype(1:1) == 'P' .or. atype(1:1) == 'p') then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atype(1:1) == 'S' .or. atype(1:1) == 's') then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atype == 'MG' .or. atype == 'mg') then
               !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !             Mg radius = 1.45A: Aqvist 1992
               x(L165-1+i) = 1.18d0 + 1.4d0
            else
               write( 0,* ) 'bad atom type: ',atype
               call mexit( 6,1 )
            end if  !  atype(1:1) == 'N'
            x(L170-1+i) = 0.0d0
            x(L175-1+i) = 0.0d0
            x(L180-1+i) = 0.0d0
            x(L185-1+i) = 0.0d0
            !  write(6,*) i,' ',atype,x(L165-1+i)
         end do  !  i=1,natom

      end if ! ( gbsa == 1 )
      
   end if  ! ( igb /= 0 .and. igb /= 10)
   
   !------------------------------------------------------------------------
   ! If user has requested Poisson-Boltzmann electrostatics, set up variables
   !------------------------------------------------------------------------
  
   if ( igb == 10 ) then
      call pb_init(ifcap,natom,nres,ntypes,nbonh,nbona,ix(i02),ix(i04),ix(i06),ix(i08),ix(i10),&
                   ix(iibh),ix(ijbh),ix(iiba),ix(ijba),ix(ibellygp),ih(m02),ih(m04),ih(m06),x(l15),x(l97))
   end if  ! ( igb == 10 ) 

   if (icnstph /= 0) then
      !     Read charge data and alter current charges accordingly
      call cnstphread(ix(icpstinf),ix(icpresst),ix(icpptcnt), &
            ix(icptrsct),x(lcpene),x(lcpcrg),x(l15))

      !     Fill proposed charges array from current charges
      do i=1,natom
         x(l190-1+i) = x(l15-1+i)
      end do
   end if

!  +---------------------------------------------------------------+
!  |  Read EVB input file                                          |
!  +---------------------------------------------------------------+

   if( ievb /= 0 ) then
#ifdef MPI
      call evb_input
      call evb_init
#  if defined(LES)
!KFW  call evb_pimd_init
#  endif
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      call mexit(6,1)
#endif
   endif

   if( iyammp /= 0 ) write( 6, '(a)' ) '  Using yammp non-bonded potential'
   
   ! -------------------------------------------------------------------
   !
   ! -- add check to see if the space in nmr.h is likely to be
   !     too small for this run:
   !     [Note: this check does *not* indicate if MXTAU and MXP are
   !      too small.  There is no easy way to ensure this, since
   !      the experimental intensities are read in a namelist
   !      command: if too many intensities are input, the read
   !      statment may cause a coredump before returning control
   !      to the main program.  Be careful.  sigh....]
   
   if (natom > matom .and. nmropt > 1) then
      write(6,*) 'WARNING: MATOM in nmr.h is smaller than the ', &
            natom,' atoms in this molecule.'
      write(6,*) 'Printout of NMR violations may be compromised.'
   end if
   
   ! -------------------------------------------------------------------
   !     --- checks on bogus data ---
   ! -------------------------------------------------------------------
   
   inerr = 0
      
   if( icfe < 0 .or. icfe > 1 ) then
      write(6,*) 'icfe must be 0 or 1 (icfe=2 is no longer supported)'
      inerr = 1
   end if
   if( icfe /= 0 .and. numgroup /= 2 ) then
      write(6,*) 'numgroup must be 2 if icfe is set'
      inerr = 1
   end if
   if (ievb>0) then
#ifdef MPI
!KFW  if( numgroup /= 2 ) then
!KFW     write(6,*) 'numgroup must be 2 if ievb is set'
!KFW     inerr = 1
!KFW  end if
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      inerr = 1
#endif
   end if
   if( igb > 0 .and. numextra > 0) then
      write(6,'(a)') 'Cannot use igb>0 with extra-point force fields'
      inerr = 1
   end if
   if (ips < 0 .or. ips > 3) then
      write(6,'(/2x,a,i3,a)') 'IPS (',ips,') must be 0,1,2, or 3'
      inerr = 1
   end if
   if (ips /= 0 .and. ipol /= 0 ) then
      write(6,'(/2x,a)') 'IPS and IPOL are inconsistent options'
      inerr = 1
   endif

   if (jar < 0 .or. jar > 1) then
      write(6,'(/2x,a,i3,a)') 'JAR (',jar,') must be 0 or 1'
      inerr = 1
   end if

   if (igb > 0 .and. ips > 0 ) then
      write(6,'(/2x,a,i3,a,i3,a)') 'IGB (',igb,') and ips (',ips, &
          ') cannot both be turned on'
      inerr = 1
   end if
   if (igb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 &
         .and. igb /=6 .and. igb /=7 .and. igb /=10) then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 0,1,2,5,6,7 or 10.'
      inerr = 1
   end if
   if (alpb /= 0 .and. alpb /= 1 )  then
      write(6,'(/2x,a,i3,a)') 'ALPB (',alpb,') must be 0 or 1.'
      inerr = 1
   end if
   if (alpb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 .and. igb /=7 )  then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 1,2,5, or 7 if ALPB > 0.'
      inerr = 1
   end if
#ifdef LES
   if( igb /= 0 .and. igb /= 1 .and. igb /= 5 .and. igb /=7 ) then
      write(6,'(/,a)') 'Error: LES is only compatible with IGB > 0,1,5,7'
      inerr = 1
   end if
   if( alpb /= 0) then
      write(6,'(/,a)') 'Error: LES is not compatible with ALPB'
      inerr = 1
   end if
   if( gbsa > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with GBSA > 0'
      inerr = 1
   end if
   if( qmmm_nml%ifqnt ) then
      write(6,'(/,a)') 'Error: LES is not compatible with QM/MM'
      inerr = 1
   end if
   if( ipol /= 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with IPOL > 0'
      inerr = 1
   end if
   if (temp0les >= 0.d0 .and. iscale > 0 ) then
      write (6,'(/,a)') 'Error: iscale cannot be used with temp0les'
      inerr = 1
   end if
#endif
   if (irest /= 0 .and. irest /= 1) then
      write(6,'(/2x,a,i3,a)') 'IREST (',irest,') must be 0 or 1.'
      inerr = 1
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6,'(/2x,a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1.'
      inerr = 1
   end if
   if (imin < 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be >= 0.'
      inerr = 1
   end if
   if (imin == 5 .and. ioutfm /= 0 ) then
      write(6,'(/2x,a,i3,a)') 'IMIN=5 currently requires formatted mdcrd/inptraj (ioutfm=0).'
      inerr = 1
   end if
   if (imin == 5 .and. ifbox /= 0 .and. ntb /= 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN=5 does not support periodic boundaries (ifbox>0, ntb>0).'
      inerr = 1
   end if


   if (iscale > mxvar) then
      write(6,9501) iscale,mxvar
      9501 format('ERROR: ISCALE (',i5,') exceeds MXVAR (',i5, &
            '). See nmr.h')
      inerr = 1
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      inerr = 1
   end if
   if (ntxo /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTXO (',ntxo,') must be 1.'
      write(6,'(/2x,a)') '  (ntx0=0 is no longer supported)'
      inerr = 1
   end if
   
   if (ntb /= 0 .and. ntb /= 1 .and. ntb /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTB (',ntb,') must be 0, 1 or 2.'
      inerr = 1
   end if
   if (ntb == 0 .and. iwrap == 1) then
      write(6,'(/2x,a)') 'Error: IWRAP=1 cannot be used without a periodic box.'
      inerr = 1
   end if
   
   if (ntt < 0 .or. ntt > 4) then
      write(6,'(/2x,a,i3,a)') 'NTT (',ntt,') must be between 0 and 4.'
      inerr = 1
   end if
   if (ntt == 1 .and. tautp < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUTP (',tautp,') < DT (step size)'
      inerr = 1
   end if
   if( ntt < 3 .or. ntt > 4 ) then
      if( gamma_ln > 0.d0 ) then
         write(6,'(a)') 'ntt must be 3 or 4 if gamma_ln > 0'
         inerr = 1
      end if
   end if
   
   if (ntp /= 0 .and. ntp /= 1 .and. ntp /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTP (',ntp,') must be 0, 1 or 2.'
      inerr = 1
   end if
   if (ntp > 0 .and. taup < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUP (',taup,') < DT (step size)'
      inerr = 1
   end if
   if (npscal < 0 .or. npscal > 1) then
      write(6,'(/2x,a,i3,a)') 'NPSCAL (',npscal,') must be 0 or 1.'
      inerr = 1
   end if
   
   if (ntc < 1 .or. ntc > 4) then
      write(6,'(/2x,a,i3,a)') 'NTC (',ntc,') must be 1,2,3 or 4.'
      inerr = 1
   end if
   if (jfastw < 0 .or. jfastw > 4) then
      write(6,'(/2x,a,i3,a)') 'JFASTW (',jfastw,') must be 0->4.'
      inerr = 1
   end if
   
   if (ntf < 1 .or. ntf > 8) then
      write(6,'(/2x,a,i3,a)') 'NTF (',ntf,') must be in 1..8.'
      inerr = 1
   end if
   
   if (scee == 0.0d0) then
      write(6,'(/2x,a)') 'SCEE must be set explicitly'
      inerr = 1
   end if
   
   if (ioutfm /= 0 .and. ioutfm /= 1) then
      write(6,'(/2x,a,i3,a)') 'IOUTFM (',ioutfm,') must be 0 or 1.'
      inerr = 1
   end if
   
   if (ntpr < 0) then
      write(6,'(/2x,a,i3,a)') 'NTPR (',ntpr,') must be >= 0.'
      inerr = 1
   end if
   if (ntwx < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWX (',ntwx,') must be >= 0.'
      inerr = 1
   end if
   if (ntwv < -1) then
      write(6,'(/2x,a,i3,a)') 'NTWV (',ntwv,') must be >= -1.'
      inerr = 1
   end if
   if (ntwv == -1 .and. ioutfm /= 1) then
      write (6, '(/2x,a)') 'IOUTFM must be 1 for NTWV == -1.'
      inerr = 1
   end if
   if (ntwv == -1 .and. ntwx == 0) then
      write (6, '(/2x,a)') 'NTWX must be > 0 for NTWV == -1.'
      inerr = 1
   end if
   if (ntwe < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWE (',ntwe,') must be >= 0.'
      inerr = 1
   end if
   if (ntave < 0) then
      write(6,'(/2x,a,i3,a)') 'NTAVE (',ntave,') must be >= 0.'
      inerr = 1
   end if
   if (ntr /= 0 .and. ntr /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTR (',ntr,') must be 0 or 1.'
      inerr = 1
   end if
   if (ntrx /= 0 .and. ntrx /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTRX (',ntrx,') must be 1 or 0.'
      inerr = 1
   end if
   if (nmropt < 0 .or. nmropt > 2) then
      write(6,'(/2x,a,i3,a)') 'NMROPT (',nmropt,') must be in 0..2.'
      inerr = 1
   end if
   
   if (idecomp < 0 .or. idecomp > 4) then
      write(6,'(/2x,a)') 'IDECOMP must be 0..4'
      inerr = 1
   end if
      
   ! check settings related to ivcap

   if(ivcap == 3 .or. ivcap == 4) then
      write(6,'(/2x,a)') 'IVCAP == 3 and IVCAP == 4 currently not implemented'
      inerr = 1
   endif
   if (ivcap < 0 .and. ivcap > 5) then
      write(6,'(/2x,a)') 'IVCAP must be 0 ... 5'
      inerr = 1
   end if
   if ((ivcap == 1 .or. ivcap == 5) .and. igb /= 10) then
      write(6,'(/2x,a)') 'IVCAP == 1,5 only works with igb == 10'
      inerr = 1
   end if
   if((ivcap == 1 .or. ivcap == 3 .or. ivcap == 5 ) .and. cutcap <= 0.0d0) then
      write(6,'(/2x,a)') 'For IVCAP == 1,3, or 5, cutcap must be > 0.0'
      inerr = 1
   endif
   if (ivcap == 4 .and. &
         (xlorth < ZERO .or. ylorth < ZERO .or. zlorth < ZERO .or. &
!  give magic numbers a name  srb aug 2007 !
          xorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          yorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          zorth > 47114710.0d0)) then
      write(6,'(/2x,a)') &
      'For IVCAP == 4, xlorth, ylorth, zlorth, xorth, yorth, zorth must be set'
      inerr = 1
   end if
   if ((ivcap == 3 .or. ivcap == 4) .and. ibelly == 0) then
      write(6,'(/2x,a,a)') &
         'For IVCAP == 3 or 4, ibelly must be 1 and all atoms', &
         '  not in the spherical or orthorhombic region must be set NOT moving'
      inerr = 1
   end if
   if (ivcap == 5 .and. (imin /= 1 .or. maxcyc > 1)) then
      write(6,'(/2x,a,a)') &
         'IVCAP == 5 only works for single-point energy calculation'
      inerr = 1
   end if

   ! check if ifbox variable from prmtop file matches actual angles:

   if ( igb == 0 .and. ntb /= 0 ) then
      if ( ifbox == 1 ) then
         if ( abs(alpha - 90.0d0) > 1.d-5 .or. &
           abs(beta  - 90.0d0) > 1.d-5 .or. &
           abs(gamma - 90.0d0) > 1.d-5 ) then
           ifbox =3
           write(6,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
         end if
      end if

      if ( ifbox == 2 ) then
         if ( abs(alpha - 109.4712190d0) > 1.d-5 .or. &
              abs(beta  - 109.4712190d0) > 1.d-5 .or. &
              abs(gamma - 109.4712190d0) > 1.d-5 ) then
              write(6,'(/2x,a)') &
              'Error: ifbox=2 in prmtop but angles are not correct'
              inerr = 1
         end if
      end if
   end if

   ! checks for targeted MD
   if (itgtmd /= 0 .and. itgtmd /= 1) then
      write(6,'(/2x,a,i3,a)') 'ITGTMD (',itgtmd,') must be 0 or 1.'
      inerr = 1
   end if
   if (itgtmd == 1 .and. ntr == 1) then
      if (len_trim(tgtfitmask) > 0 .or. len_trim(tgtrmsmask) <= 0) then
         write(6,'(/2x,a)') 'ITGTMD: tgtrmsmask (and not tgtfitmask) ' //  &
                            'should be specified if NTR=1'
         inerr = 1
      end if
   end if
   ! skip this test until fallback to rgroup() is supported
   !if (itgtmd == 1 .and. ntr == 0) then
   !   if (len_trim(tgtfitmask) == 0 .and. len_trim(tgtrmsmask) == 0) then
   !      write(6,'(/2x,a)')  &
   !        'ITGTMD: both tgtfitmask and tgtrmsmask should be specified if NTR=0'
   !      inerr = 1
   !   end if
   !end if
   
   !     -- consistency checking
   
   if (imin > 0.and.nrespa > 1)  then
      write(6,'(/2x,a)') 'For minimization, set nrespa,nrespai=1'
      inerr = 1
   end if
   if (ntp > 0 .and. nrespa > 1) then
      write(6,'(/2x,a)') 'nrespa must be 1 if ntp>0'
      inerr = 1
   end if
   if  (ntx < 4.and.init /= 3)  then
      write(6,'(/2x,a)') 'NTX / IREST inconsistency'
      inerr = 1
   end if
   if (ntb == 2 .and. ntp == 0) then
      write(6,'(/2x,a)') 'NTB set but no NTP option (must be 1 or 2)'
      inerr = 1
   end if
   if (ntp /= 0 .and. ntb /= 2) then
      write(6,'(/,a,a)')' NTP > 0 but not constant pressure P.B.C.', &
            ' (NTB = 2) must be used'
      inerr = 1
   end if
   if (ntb /= 0 .and. ifbox == 0 .and. ntp /= 0) then
      write(6,'(/,a)') ' (NTB /= 0 && NTP /= 0) but IFBOX == 0'
      write(6,'(/,a)') ' This combination is not supported'
      inerr = 1
   end if
   if (ntb /= 0 .and. &
         ( box(1) < 1.d0  .or. &
         box(2) < 1.d0  .or. &
         box(3) < 1.d0 ) ) then
      write(6,'(/,a,3f10.3)') ' BOX is too small: ',box(1),box(2),box(3)
      inerr = 1
   else if (ntb /= 0 .and. &
         (sqrt(cut) >= box(1)*0.5d0 .or. &
         sqrt(cut) >= box(2)*0.5d0 .or. &
         sqrt(cut) >= box(3)*0.5d0) ) then
      write(6,'(/,a)') ' CUT must be < half smallest box dimension'
      inerr = 1
   end if
   if (ntb /= 0 .and. igb > 0 ) then
      write(6,'(/,a)') ' igb>0 is only compatible with ntb=0'
      inerr = 1
   end if
#ifdef APBS
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 .and. &
      .not. mdin_apbs) then
#else
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 ) then
#endif /* APBS */
      write(6,'(/,a,f8.2)') ' unreasonably small cut for non-periodic run: ', &
         sqrt(cut)
      inerr = 1
   end if
   if ( rgbmax < 5.d0*fsmax ) then
      write(6,'(/,a,f8.2)') ' rgbmax must be at least ', 5.d0*fsmax
      inerr = 1
   end if
   if (icfe /= 0 .and. indmeth == 3 ) then
      write(6,'(/,a)') ' indmeth=3 cannot be used with icfe>0'
      inerr = 1
   end if
   if (icfe /= 0 .and. ibelly /= 0 ) then
      write(6,'(/,a)') ' ibelly cannot be used with icfe'
      inerr = 1
   end if

   ! Modification done by Ilyas Yildirim
   if (icfe == 1 .and. (klambda < 1 .or. klambda > 6)) then
     write(6,'(/,a)') ' klambda must be between 1 and 6'
     inerr = 1
   end if
   ! End of modification done by Ilyas Yildirim                                       

   if (clambda < 0.d0 .or. clambda > 1.d0 ) then
      write(6,'(/,a)') ' clambda must be between 0 and 1'
      inerr = 1
   end if

   if (icfe /= 0 .and. (idecomp == 3 .or. idecomp == 4)) then
      write(6,'(/,a)') ' Pairwise decomposition for thermodynamic integration not implemented'
      inerr = 1
   end if
   if (icfe /= 0 .and. idecomp /= 0 .and. ipol /= 0) then
      write(6,'(/,a)') ' IPOL is incompatible with IDECOMP and ICFE'
      inerr = 1
   end if
 
#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      if (icfe /= 1 .and. ifsc==1) then
         write (6,'(/,a)') ' Softcore potential requires a standard TI run, set icfe to 1'
         inerr = 1
      end if
      if ( igb > 0 ) then
         write (6,'(/,a)') ' Softcore potential is incompatible with GB (for now)'
         inerr = 1
      end if
      if ( ntf > 1 ) then
         write (6,'(/,a)') ' Softcore potentials require ntf=1 because SHAKE constraints on some bonds might be removed'
         inerr = 1
      end if
      if (clambda > 0.99 .or. clambda < 0.01) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with clambda < 0.01 or > 0.99'
         inerr = 1
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' Softcore potential requires linear mixing, set klambda to 1'
         inerr = 1
      end if
      if (imin == 1 .and. ntmin /= 2) then
         write (6,'(/,a)') ' Minimizations with ifsc=1 require the steepest descent algorithm.'
         write (6,'(/,a)') ' Set ntmin to 2 and restart'
         inerr = 1
      end if
   end if
#endif
  
   !check for neb
   if (idecomp > 0 .and. (ntr > 0 .or. ibelly > 0)) then
      write(6,'(/,a)') 'IDECOMP is not compatible with NTR or IBELLY'
      inerr = 1
   end if
   if (icnstph /= 0) then
      if (igb == 0) then
         write(6, '(/,a)') 'Constant pH requires GB implicit solvent'
         inerr = 1
      end if
      if (icfe /= 0) then
         write(6, '(/,a)') &
         'Constant pH and thermodynamic integration are incompatable'
         inerr = 1
      end if
   end if

#ifdef noVIRIAL
   if( ntp > 0 ) then
      write(6,'(/,a)') 'Error: constant pressure is incompatible with noVIRIAL'
      inerr = 1
   end if
#endif
   
   !-----------------------------------------------------
   !     ----sanity checks for Ewald
   !-----------------------------------------------------
   
   if( igb == 0 ) then
      call float_legal_range('skinnb: (nonbond list skin) ', &
            skinnb,skinlo,skinhi)
      
      !  --- Will check on sanity of settings after coords are read in
      !      and the extent of the system is determined.
      
      if(periodic == 1)then
         call float_legal_range('skinnb+cutoffnb: (nonbond list cut) ', &
               skinnb+cutoffnb,zero,sphere)
      end if
      if (ntb==0 .and. use_pme/=0) then
         write(6,'(/,a)') &
         'Using PME with a non-periodic simulation does not make sense. Set either ntb>0 of use_pme=0.'
         inerr = 1
      end if
      call float_legal_range('a: (unit cell size) ',a,boxlo,boxhi)
      call float_legal_range('b: (unit cell size) ',b,boxlo,boxhi)
      call float_legal_range('c: (unit cell size) ',c,boxlo,boxhi)
      call float_legal_range('alpha: (unit cell angle) ', &
            alpha,anglo,anghi)
      call float_legal_range('beta: (unit cell angle) ', &
            beta,anglo,anghi)
      call float_legal_range('gamma: (unit cell angle) ', &
            gamma,anglo,anghi)
      call int_legal_range('order: (interpolation order) ', &
            order,orderlo,orderhi)
      call opt_legal_range('verbose: ',verbose,0,4)
      call opt_legal_range('netfrc: ',netfrc,0,1)
      call opt_legal_range('nbflag: ',nbflag,0,1)
      call opt_legal_range('nbtell: ',nbtell,0,2)
      call opt_legal_range('ew_type: ',ew_type,0,1)
      call opt_legal_range('vdwmeth: ',vdwmeth,0,2)
      call opt_legal_range('eedmeth: ',eedmeth,1,6)
      call opt_legal_range('ee_type: ',ee_type,1,2)
      call opt_legal_range('maxiter: ',maxiter,1,50)
      call opt_legal_range('indmeth: ',indmeth,0,3)
      call opt_legal_range('fix_quad: ',fix_quad,0,1)
      call float_legal_range('eedtbdns: (erfc table density) ', &
            eedtbdns,denslo,denshi)
   end if  ! ( igb == 0 )

   if( ntb==2 .and. ipimd==1) then
      write(6,*) 'primitive PIMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ntb==2 .and. ipimd==3 ) then
      write(6,*) 'CMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ipimd==3 .and. adiab_param>=1.0 ) then
      write(6,*) 'For CMD adiab_param must be <=1'
      inerr=1
   endif

   if( ntb==2 .and. ipimd==4 ) then
      write(6,*) 'RPMD is incompatible with NTP ensemble'
      inerr=1
   endif

   if( ntt/=0 .and. ipimd==4 ) then
      write(6,*) 'RPMD is incompatible with NVT ensemble'
      inerr=1
   endif

   if( ntt/=4 .and. ipimd==2 ) then
      write(6,*) 'NMPIMD requires Nose-Hoover chains (ntt=4)'
      inerr=1
   endif

   if( ntt/=4 .and. ipimd==3 ) then
      write(6,*) 'CMD requires Nose-Hoover chains (ntt=4)'
      inerr=1
   endif

   if( ineb > 0 .and. ipimd > 0 ) then
      write(6,*) 'ineb>0 and ipimd>0 are incompatible options'
      inerr=1
   endif

   if( iamoeba == 1 )then
#ifdef LES
      write(6,*)'amoeba is incompatible with LES'
      inerr=1
#endif

      if( ntc > 1 ) then
         write(6,*) 'SHAKE (ntc>1) and amoeba are incompatible options'
         inerr=1
      end if
      if( ntp > 1 .and. beeman_integrator > 0 ) then
         write(6,*) 'ntp>1 is not consistent with the beeman integrator'
         inerr=1
      end if

   end if

    
   ! ---WARNINGS:
   
   if ( ibelly == 1 .and. igb == 0 .and. ntb /= 0 ) then
      write(6,'(/,a,/,a,/,a)') &
            'Warning: Although EWALD will work with belly', &
            '(for equilibration), it is not strictly correct!'
   end if
   
   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      call mexit(6,1)
   end if
   
   ! Load the restrained atoms (ntr=1) or the belly atoms (ibelly=1)
   ! or atoms for targeted md (itgtmd=1). Selections are read from
   ! &cntrl variables or, if these are not defined, it falls back to
   ! the old group input format.
   
   konst = ntr > 0
   dotgtmd = itgtmd > 0
   belly = .false.
   natc = 0
   ngrp = 0
   natbel = 0
   nattgtfit = 0  ! number of atoms for tgtmd fitting (=overlap)
   nattgtrms = 0  ! number of atoms for tgtmd rmsd calculation
   nrc = 0
   if(konst.or.dotgtmd) then
      if (ntrx <= 0) then
         call amopen(10,refc,'O','U','R')
      else
         call amopen(10,refc,'O','F','R')
      end if
      ! these messages should be written after "5.  REFERENCE..." ?
      if (konst) write(6,9408)
      if (dotgtmd) write(6,9409)

      call rdrest(natom,ntrx,x(lcrdr))
      close(10)

      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.)
#else
      call AMOEBA_check_newstyle_inpcrd(inpcrd,newstyle)
      if ( newstyle )then
         call AM_RUNMD_get_coords(natom,t,irest,ntb,x(lcrd),x(lvel))
      else
         if( irest == 1 .and. beeman_integrator > 0 ) then
            write(6,*) 'Cannot do a beeman_integrator restart with old-style coordinates'
            call mexit(6,1)
         end if
         call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.)
      endif
#endif
      
      ! VH - tgtmd change: preferably call atommask() instead of rgroup()
      if (konst) then
         if( len_trim(restraintmask) <= 0 ) then
            call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02), &
                  ih(m04),ih(m06),ih(m08),ix(icnstrgp),jgroup,indx,irespw,npdec, &
                  x(l60),x(lcrdr),konst,dotgtmd,belly,idecomp,5,.true.)
         else
            call atommask( natom, nres, 0, ih(m04), ih(m06), &
               ix(i02), ih(m02), x(lcrd), restraintmask, ix(icnstrgp) )

            ! for now, emulate the "GATHER ALL THE CONSTRAINED ATOMS TOGETHER"
            ! section of rgroup(); later, the various masks should be done
            ! differently, i.e. without the "gather", as in the following:
            !     x(l60:l60+natom-1) = restraint_wt
            !     natc = sum(ix(icnstrgp:icnstrgp+natom-1))

            natc = 0
            do i=1,natom
              if( ix(icnstrgp-1+i) <= 0 ) cycle
              natc = natc + 1
              ix(icnstrgp-1+natc) = i
              x(l60-1+natc) = restraint_wt
            end do
            write(6,'(a,a,a,i5,a)') '     Mask ', &
            restraintmask(1:len_trim(restraintmask)), ' matches ',natc,' atoms'
         end if
      end if
      nrc = natc
      
      if (dotgtmd) then
         if (len_trim(tgtfitmask) <= 0 .and. len_trim(tgtrmsmask) <= 0) then
            ! the following if-endif can be deleted when we stop
            ! supporting rgroup()
            if (konst) then
               ! cannot do both ntr and tgtmd together using old group format
               write(6,'(/2x,a)') 'NTR must be 0 for targeted MD (TGTMD=1)'
               call mexit(6,1)
            else  ! the following only for backward compatibility
               call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02), &
                     ih(m04),ih(m06),ih(m08),ix(icnstrgp), &
                     jgroup,indx,irespw,npdec, &
                     x(l60),x(lcrdr),konst,dotgtmd,belly,idecomp,5,.true.)
               ! tgtmd atoms are now stored in nattgt, igroup -> icnstrgp
               nattgtfit = natc
               nattgtrms = natc
               do i=1,nattgtfit
                  ix(itgtfitgp-1+i) = ix(icnstrgp-1+i)
                  ix(itgtrmsgp-1+i) = ix(icnstrgp-1+i)
               end do
            end if
         else
            if (ntr == 0) then  ! read tgtfitmask only if ntr=1 
               ! read in atom group for tgtmd fitting (=overlap region)
               call atommask( natom, nres, 0, ih(m04), ih(m06), &
                  ix(i02), ih(m02), x(lcrd), tgtfitmask, ix(itgtfitgp) )
               ! see comments above (for ntr) for the following reduction cycle
               nattgtfit = 0
               do i=1,natom
                 if( ix(itgtfitgp-1+i) <= 0 ) cycle
                 nattgtfit = nattgtfit + 1
                 ix(itgtfitgp-1+nattgtfit) = i
               end do
               write(6,'(a,a,a,i5,a)')  &
               '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
               '" matches ',nattgtfit,' atoms'
            end if
            ! read in atom group for tgtmd rmsd calculation
            call atommask( natom, nres, 0, ih(m04), ih(m06), &
               ix(i02), ih(m02), x(lcrd), tgtrmsmask, ix(itgtrmsgp) )
            nattgtrms = 0
            do i=1,natom
              if( ix(itgtrmsgp-1+i) <= 0 ) cycle
              nattgtrms = nattgtrms + 1
              ix(itgtrmsgp-1+nattgtrms) = i
            end do
            write(6,'(a,a,a,i5,a)')  &
            '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
            '" matches ',nattgtrms,' atoms'
         end if
      end if

   end if  ! (konst.or.dotgtmd)
   
   ! dotgtmd may be false here even if doing tgtmd
   ! this is so belly info is read properly? following existing KONST code
   
   dotgtmd=.false.
   konst = .false.
   belly = ibelly > 0
   ngrp = 0
   if(belly) then
      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.)
#else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.)
#endif
      write(6,9418)
      if( len_trim(bellymask) <= 0 ) then
         call rgroup(natom,natbel,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,indx,irespw,npdec, &
            x(l60),x(lcrdr),konst,dotgtmd,belly,idecomp,5,.true.)
      else
         call atommask( natom, nres, 0, ih(m04), ih(m06), &
            ix(i02), ih(m02), x(lcrd), bellymask, ix(ibellygp) )
         natbel = sum(ix(ibellygp:ibellygp+natom-1))
         write(6,'(a,a,a,i5,a)') '     Mask ', &
            bellymask(1:len_trim(bellymask)), ' matches ',natbel,' atoms'
      end if
   end if
   call setvar(ix,belly)

   !  see if the user has input a noshakemask string, and process it:
   natnos = 0
   if( len_trim(noshakemask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), noshakemask, noshakegp )
      natnos = sum(noshakegp(1:natom))
      write(6,*)
      write(6,'(a,a,a,i5,a)') 'Noshake mask ', &
         noshakemask(1:len_trim(noshakemask)), ' matches ',natnos,' atoms'
      call setnoshake(ix,noshakegp,ntc,num_noshake)
      if( ntf > 1 ) then
         write(6,'(a)') '   Setting ntf to 1'
         ntf = 1
      end if
   end if

#ifdef MPI /* SOFT CORE */
   ! lower charges if a crgmask is set
   if ( len_trim(crgmask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), crgmask, crggp )
      write(6,'(a,a,a,i5,a)') 'Zero-Charge Mask ',crgmask(1:len_trim(crgmask)), ' matches ',sum(crggp(1:natom)),' atoms'
      call remove_charges(crggp, natom, x(l15))
   end if
#endif

   konst = .false.
   belly = .false.
   if(idecomp > 0) then
      write(6,9428)
      call rgroup(natom,ntmp,nres,ngrp,ix(i02),ih(m02), &
            ih(m04),ih(m06),ih(m08),ix(ibellygp), &
            jgroup,indx,irespw,npdec, &
            x(l60),x(lcrdr),konst,dotgtmd,belly,idecomp,5,.true.)
   end if
   
   if( ibelly > 0 .and. igb > 0 ) then
      
      !          ---here, the only allowable belly has just the first
      !             NATBEL atoms in the moving part.  Check to see that this
      !             requirement is satisfied:
      
      do i=natbel+1,natom
         if( ix(ibellygp+i-1) /= 0 ) then
            write(6,*) 'When igb>0, the moving part must be at the'
            write(6,*) '   start of the molecule.  This does not seem'
            write(6,*) '   to be the case here.'
            write(6,*) 'natbel,i,igroup(i) = ' &
                  ,natbel,i,ix(ibellygp+i-1)
            call mexit(6,1)
         end if
      end do
   end if
   

   !     ----- CALCULATE THE SQUARE OF THE BOND PARAMETERS FOR SHAKE
   !           THE PARAMETERS ARE PUT SEQUENTIALLY IN THE ARRAY CONP -----
   
   do i=1,nbonh + nbona + nbper
      j = ix(iicbh+i-1)
      x(l50+i-1) = req(j)**2
   end do
   
#ifdef MPI
      if( icfe /= 0 ) then

      !  use the masses of the prmtop file for the first group for both groups:
      !  [only the master nodes communicate here, since non-master nodes
      !   have not yet allocated space]
      ! This leads to problems for dual topology runs, and is therefore skipped
      ! if ifsc is set to one, the masses from both prmtop files are used
         if (ifsc == 0) then
            call mpi_bcast(x(lmass),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
            call mpi_bcast(x(lwinv),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
         end if
         tmass = sum(x(lmass:lmass+natom-1))
         tmassinv = 1.d0/tmass

      !  next, do a minimal sanity check that the SHAKE parameters are
      !  consistent on the two processors:

         ! For Softcore this might be allowed
         ! Put a better check here later
         if( ntc == 2 .and. ifsc == 0) then
            partner = ieor(masterrank,1)
            call mpi_sendrecv( nbonh, 1, MPI_INTEGER, partner, 5, &
                               nbonh_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            call mpi_sendrecv( num_noshake, 1, MPI_INTEGER, partner, 5, &
                               num_noshake_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            if( nbonh - num_noshake /= nbonh_c - num_noshake_c ) then
               write(6,*) 'SHAKE lists are not compatible in the two groups!'
               call mexit(6,1)
            end if
         else if( ntc == 3 ) then
            write(6,*) 'ntc = 3 is not compatible with icfe>0'
            call mexit(6,1)
         end if

      end if
#endif

   if ( iamoeba /= 1 )then
      if( igb == 0 ) &
         call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10),fmn, &
         nspm,ix(i70),x(l75),tmass,tmassinv,x(lmass),x(lwinv),req)
   endif
  
   !  DEBUG input; force checking
   call load_debug(5)

   return
   ! -------------------------------------------------------------------------
   ! Standard format statements:
   
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9408 format(/4x,'LOADING THE CONSTRAINED ATOMS AS GROUPS',/)
   9409 format(/4x,'LOADING THE TARGETED MD ATOMS AS GROUPS',/)
   9418 format(/4x,'LOADING THE BELLY ATOMS AS GROUPS',/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)
end subroutine mdread2 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit defined preprocessor names, ie, flags.
subroutine printflags()

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(len=max_line_length) line  ! output string of active flags
   integer n                            ! len(line)
   
   line = '| Flags:'
   n = 8
   
#ifdef ISTAR2
   call printflags2(' ISTAR2',7,n,line,.false.)
#endif
#ifdef MPI
   call printflags2(' MPI',4,n,line,.false.)
# ifdef USE_MPI_IN_PLACE
   call printflags2(' USE_MPI_IN_PLACE',17,n,line,.false.)
# endif
#endif
#ifdef LES
   call printflags2(' LES',4,n,line,.false.)
#endif
#ifdef NMODE
   call printflags2(' NMODE',6,n,line,.false.)
#endif
#ifdef HAS_10_12
   call printflags2(' HAS_10_12',10,n,line,.false.)
#endif
#ifdef DNA_SHIFT
   call printflags2(' DNA_SHIFT',10,n,line,.false.)
#endif
#ifdef MMTSB
   call printflags2(' MMTSB',6,n,line,.false.)
#endif

#ifdef noVIRIAL
   call printflags2(' noVIRIAL',9,n,line,.false.)
#endif
   
   call printflags2(' ',1,n,line,.true.)
   return
end subroutine printflags 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Primitive pre-Fortran90 implementation of printflags.
subroutine printflags2(flag,flag_len,line_len,line,last)

   implicit none
   integer     max_line_length
   parameter ( max_line_length = 80 )

   character(*) flag                ! flag name with blank prefix, intent(in)
   integer flag_len                 ! len(flag), intent(in)
   integer line_len                 ! len(line), intent(inout)
   character(len=max_line_length) line ! intent(inout)
   logical last                     ! is this the last flag ?, intent(in)

   if (line_len + flag_len > max_line_length) then
      write( 6,'(a)') line
      ! begin another line
      line = '| Flags:'
      line_len=8
   end if
   line=line(1:line_len) // flag(1:flag_len)
   line_len=line_len+flag_len
   if(last)write( 6,'(a)') line
   return
end subroutine printflags2 

!-------------------------------------------------
!     --- FLOAT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range 

!-------------------------------------------------
!     --- INT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range 

!-------------------------------------------------
!     --- OPT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer option; abort on illegal values.
subroutine opt_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald OPTION CHECKING: ')
   60 format(1x,'option ',a,' has value ',i5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i5,' Upper limit: ',i5)
   63 format(1x,'Check the manual')
   return
end subroutine opt_legal_range 

!-------------------------------------------------
!     --- SANDER_BOMB ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
   implicit none
   character(len=*) routine,string1,string2

   write(6, '(1x,2a)') &
         'SANDER BOMB in subroutine ', routine
   write(6, '(1x,a)') string1
   write(6, '(1x,a)') string2
   call mexit(6,1)
end subroutine sander_bomb
!-------------------------------------------------

!-------------------------------------------------
!     --- remove_charges ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Zero charges on some atoms
subroutine remove_charges(crggp,natom,charge)
  use constants, only: INV_AMBER_ELECTROSTATIC
  implicit none
  integer natom, crggp(*),i
  _REAL_ charge(*), charge_removed
  
  charge_removed = 0.d0
  do i=1,natom
     if (crggp(i)==1) then
        charge_removed = charge_removed + charge(i) * INV_AMBER_ELECTROSTATIC
        write (6,'(a,f12.4,a,i5)') 'Removing charge of ', charge(i) * INV_AMBER_ELECTROSTATIC,' from atom ',i
        charge(i)=0
     end if
  end do
  write(6, '(a,f12.4,a)') 'Total charge of ',charge_removed,' removed'
  RETURN
end subroutine remove_charges
!-------------------------------------------------
