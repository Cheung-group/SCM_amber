#include "copyright.h"
#include "dprec.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine locmem here]
subroutine locmem()
   
   
   !     locmem:  partitions core array into storage for all
   !        the major arrays of the program.
   use nblist, only: cutoffnb,skinnb
   use amoeba_mdin, only: iamoeba,am_nbead
   implicit none
   
#  include "box.h"
#  include "nmr.h"
#  include "memory.h"
#  include "md.h"
#  include "tgtmd.h"
#  include "ew_cntrl.h"
#  include "extra_pts.h"
#  include "dynph.h"
#  include "sgld.h"
#ifdef MPI
#  include "parallel.h"
#endif
   integer none,ntbond,ntangl,ntdih,m7,istartr, &
         istarti, iendr,iendi, istomp,ida_max
   integer r_ptr,i_ptr,h_ptr,maxpr
   _REAL_ maxpr_float,natom_float,n2_float
   
   
   !     --- Identification of REAL arrays ---
   
   !     CG      ...  L15   ! PARTIAL CHARGES FOR ATOMS
   !     AMASS   ...  LWINV ! ATOMIC MASSES (inverted in rdparm - see Lmass)
   !     XCHRG   ...  Lpol  ! atomic polarizibilities
   !     C       ...  LCRD  ! COORDINATES
   !     F       ...  Lforce! FORCE
   !     V       ...  Lvel  ! VELOCITY for MD, work space for min
   !     VOLD    ...  Lvel2 ! OLD VELOCITY for MD
   !     XR      ...  L45   ! Coords rel. to COM of each molecule
   !     CONP    ...  L50   ! BOND PARAMETER FOR SHAKE
   !     XC      ...  LCRDR ! POSITION COORDINATE FOR CONSTRAINT
   !     WEIT    ...  L60   ! WEIGHT FOR POSITION CONSTRAINT
   !                  L65   ! polarization
   !                  Lmass ! masses
   !     TMA     ...  L75   ! SUB-MOLECULAR WEIGHT ARRAY IN RUNMD
   !                  L95   ! 3*Natom Real Scratch (for pol.) or Natom (for nmr)
   !                        ! also used for SKIP array in shake (2*ntbond)
   !                  L96   ! GB "fs" array
   !                  L97   ! GB "rborn" array
   !                  L98   ! GB "reff" array
   !                  L99   ! GB "onereff" array (1/reff)
   !                 Lfrctmp! 3*Natom + 40 Real Scratch( fdist + pol.), mpi only
   !                  L105  ! NMR "xstore" variable
   !                  L110  ! NMR "fnoe" variable
   !                  L115  ! NMR "ddep" variable
   !                  L120  ! NMR "dddep" variable
   !                  L125  ! NMR "dorat" variable
   !                  L130  ! NMR "ddrat" variable
   !                  L135  ! NMR "rate" variable
   !                  L140  ! NMR "trp" variable
   !                  L145  ! NMR "dint" variable
   !                  L165  ! GB/SA TDND "vdwrad" array
   !                  L170  ! GB/SA LCPO "P1" array
   !                  L175  ! GB/SA LCPO "P2" array
   !                  L180  ! GB/SA LCPO "P3" array
   !                  L185  ! GB/SA LCPO "P4" array
   !                  L186  ! GB max of rborn array
   !                  L187  ! GB min of rborn array
   !                  L188  ! GB ave of rborn array
   !                  L189  ! GB rms fluct of rborn array
   !                  L190  ! constant pH dcharge array
   !                  Lcpcrg! Constant pHstate charges
   !                  Lcpene! Constant pHstate energies
   
   !     --- Identification of Hollerith arrays
   
   !     LBRES   ...  m02  ! RESIDUE LABEL (nres)
   !     IGRAPH  ...  m04  ! ATOM NAMES (natom)
   !     ISYMBL  ...  m06  ! ATOM SYMBOL ARRAY (natom)
   !     ITREE   ...  m08  ! ATOM TREE STRUCTURE ARRAY (natom)
   !     n14     ...  m12  ! 1*natom
   !     ni14    ...  m14  ! 15*natoms 1-4 indicies
   !     iarx    ...  m16  ! 1*natom  scratch for nonbond+1-4
   
   
   !     --- Identification of Integer arrays ---
   
   !     IPRES   ...  I02      ITA     ...  I32
   !     IAC     ...  I04      JTA     ...  I34
   !     ICO     ...  I06      KTA     ...  I36
   !     IBLO    ...  I08      ICTA    ...  I38
   !     INB     ...  I10      IPH     ...  I40
   !     IBH     ...  Iibh     JPH     ...  I42
   !     JBH     ...  Ijbh     KPH     ...  I44
   !     ICBH    ...  Iicbh    LPH     ...  I46
   !     IBA     ...  Iiba     ICPH    ...  I48
   !     JBA     ...  Ijba     IPA     ...  I50
   !     ICBA    ...  Iicba    JPA     ...  I52
   !     ITH     ...  I24      KPA     ...  I54
   !     JTH     ...  I26      LPA     ...  I56
   !     KTH     ...  I28      ICPA    ...  I58
   !     ICTH    ...  I30

   !     Cnstr  IGROUP  ...  Icnstrgp
   !     Tgtfit IGROUP  ...  Itgtfitgp
   !     Tgtrms IGROUP  ...  Itgtrmsgp
   !     Belly  IGROUP  ...  Ibellygp
   !     JOIN    ...  I64
   !     ISTORE  ...  I65
   !     NSP     ...  I70 ! SUBMOLECULE INDEX ARRAY
   !     IAR1    ...  I78
   !     NUMBOND ...  I80
   !     NUMNEAR ...  I82
   !     Constant pH state metadata ... Icpstinf
   !     Constant pH residue states ... Icpresst
   !     Constant pH state protonation levels .. Icpptcnt
   !     Dynamic protonation num titrating residues ... Icptrsct

   
   ! The IVMxx pointers are used for the following:
   
   !           Iifstwt  ...  IFSTWT (fast shake bond array)
   !           Iifstwr  .... IFSTWR (fast shake residue array)
   
   !    NMR restraints/weight changes require two storage areas. These are:
   
   !           WORKN   ...  LNMR01 (X array)
   !           IWORKN  ...  INMR02 (IX array)
   !    SGLD arrays
   !     LVSG     ...   SYSTEMATIC  Velocity
   
   !     --- assign standard partition lengths ---
   
   maxdup = 2000
   none = 0
   ntbond = nbonh  + nbona + nbper
   ntangl = ntheth + ntheta + ngper
   ntdih  = nphih  + nphia + ndper + 2*maxdup
   m7     = nphih  + maxdup
   
   !-----------------------------------------------------------------------
   !     --- set pointers for real arrays ---
   !-----------------------------------------------------------------------

   r_ptr = 1
   call adj_mem_ptr( r_ptr, l15, natom )
   call adj_mem_ptr( r_ptr, lwinv, natom*am_nbead )
   if (ipol > 0) then
      call adj_mem_ptr( r_ptr, lpol, natom )
   else
      call adj_mem_ptr( r_ptr, lpol, 0 )
   end if
   call adj_mem_ptr( r_ptr, lcrd, 3*natom*am_nbead + mxvar )
   call adj_mem_ptr( r_ptr, lforce, 3*natom*am_nbead + mxvar + 40 )
   if (imin == 0) then
      call adj_mem_ptr( r_ptr, lvel,  3*natom*am_nbead + mxvar )
      call adj_mem_ptr( r_ptr, lvel2, 3*natom*am_nbead + mxvar )
   else
      call adj_mem_ptr( r_ptr, lvel, 6*(3*natom*am_nbead + mxvar) )
      call adj_mem_ptr( r_ptr, lvel2, 0 )
   end if
   call adj_mem_ptr( r_ptr, l45, 3*natom*am_nbead + mxvar )
   call adj_mem_ptr( r_ptr, l50, ntbond )
   
   ! positional restraints or carlos added targeted MD
   
   if (ntr > 0.or.itgtmd > 0) then
      call adj_mem_ptr( r_ptr, lcrdr, 3*natom + mxvar )
      call adj_mem_ptr( r_ptr, l60, natom )
   else
      call adj_mem_ptr( r_ptr, lcrdr, 0 )
      call adj_mem_ptr( r_ptr, l60, 0 )
   end if
   if (ipol > 0) then
      call adj_mem_ptr( r_ptr, l65, 3*natom )
   else
      call adj_mem_ptr( r_ptr, l65, 0 )
   end if
   
   !   SGLD pointers
   if(tsgld)then
      call adj_mem_ptr( r_ptr, lvsg, 3*natom )
   else
      call adj_mem_ptr( r_ptr, lvsg, 0 )
   endif
   !     --- real array NMR restraints/weight changes:
   
   call adj_mem_ptr( r_ptr, lmass, natom*am_nbead )
   call adj_mem_ptr( r_ptr, lnmr01, irlreq )
   
   call adj_mem_ptr( r_ptr, l75, natom )
   if (ipol > 0) then
      call adj_mem_ptr( r_ptr, l95, max(3*natom, 2*ntbond) )
   else if (nmropt > 0 ) then
      call adj_mem_ptr( r_ptr, l95, max(natom, 2*ntbond) )
   else
      call adj_mem_ptr( r_ptr, l95, 2*ntbond )
   end if
   if( igb /= 0 .or. hybridgb>0 ) then
      call adj_mem_ptr( r_ptr, l96, natom )
      call adj_mem_ptr( r_ptr, l97, natom )
#ifdef LES
      call adj_mem_ptr( r_ptr, l98, natom*ncopy )
      call adj_mem_ptr( r_ptr, l99, natom*ncopy )
#else
      call adj_mem_ptr( r_ptr, l98, natom )
      call adj_mem_ptr( r_ptr, l99, natom )
#endif
      if ( gbsa > 0 ) then
         call adj_mem_ptr( r_ptr, l165, natom )
         call adj_mem_ptr( r_ptr, l170, natom )
         call adj_mem_ptr( r_ptr, l175, natom )
         call adj_mem_ptr( r_ptr, l180, natom )
         call adj_mem_ptr( r_ptr, l185, natom )
      else
         call adj_mem_ptr( r_ptr, l165, 0 )
         call adj_mem_ptr( r_ptr, l170, 0 )
         call adj_mem_ptr( r_ptr, l175, 0 )
         call adj_mem_ptr( r_ptr, l180, 0 )
         call adj_mem_ptr( r_ptr, l185, 0 )
      end if
      if ( rbornstat == 1 ) then
         call adj_mem_ptr( r_ptr, l186, natom )
         call adj_mem_ptr( r_ptr, l187, natom )
         call adj_mem_ptr( r_ptr, l188, natom )
         call adj_mem_ptr( r_ptr, l189, natom )
      else
         call adj_mem_ptr( r_ptr, l186, 0 )
         call adj_mem_ptr( r_ptr, l187, 0 )
         call adj_mem_ptr( r_ptr, l188, 0 )
         call adj_mem_ptr( r_ptr, l189, 0 )
      end if
   else
      call adj_mem_ptr( r_ptr, l96,  0 )
      call adj_mem_ptr( r_ptr, l97,  0 )
      call adj_mem_ptr( r_ptr, l98,  0 )
      call adj_mem_ptr( r_ptr, l99,  0 )
      call adj_mem_ptr( r_ptr, l165,  0 )
      call adj_mem_ptr( r_ptr, l170,  0 )
      call adj_mem_ptr( r_ptr, l175,  0 )
      call adj_mem_ptr( r_ptr, l180,  0 )
      call adj_mem_ptr( r_ptr, l185,  0 )
      call adj_mem_ptr( r_ptr, l186,  0 )
      call adj_mem_ptr( r_ptr, l187,  0 )
      call adj_mem_ptr( r_ptr, l188,  0 )
      call adj_mem_ptr( r_ptr, l189,  0 )
   end if  ! ( igb /= 0 )

#ifdef MPI
   call adj_mem_ptr( r_ptr, lfrctmp, 3*natom*am_nbead + 40 )
#endif

   if (nmropt >= 2) then
      call adj_mem_ptr( r_ptr, l105, mxsub*isubr )
      call adj_mem_ptr( r_ptr, l110, 3*natom + mxvar )
      call adj_mem_ptr( r_ptr, l115, ma*ma )
      call adj_mem_ptr( r_ptr, l120, 3*ma*ma )
      call adj_mem_ptr( r_ptr, l125, 3*ma*ma )
      call adj_mem_ptr( r_ptr, l130, 3*ma*ma )
      call adj_mem_ptr( r_ptr, l135, ma*ma )
      call adj_mem_ptr( r_ptr, l140, ma*ma )
      call adj_mem_ptr( r_ptr, l145, 3*ma + mxvar )
   else
      call adj_mem_ptr( r_ptr, l110, 0 )
   end if
   call adj_mem_ptr( r_ptr, l150, 0)

   if ( icnstph /= 0 ) then
      call adj_mem_ptr( r_ptr, l190, natom)
      call adj_mem_ptr( r_ptr, lcpcrg, ATOM_CHRG_C)
      call adj_mem_ptr( r_ptr, lcpene, TITR_STATES_C)
   else
      call adj_mem_ptr( r_ptr, l190, 0)
      call adj_mem_ptr( r_ptr, lcpcrg, 0)
      call adj_mem_ptr( r_ptr, lcpene, 0)
   end if

   lastr = r_ptr
   
   !-----------------------------------------------------------------------
   !     --- Allocate Hollerith Space ---
   !-----------------------------------------------------------------------
   
   h_ptr = 1
   call adj_mem_ptr( h_ptr, m02, nres + 1 )
   call adj_mem_ptr( h_ptr, m04, natom )
   call adj_mem_ptr( h_ptr, m06, natom )
   call adj_mem_ptr( h_ptr, m08, natom )
   call adj_mem_ptr( h_ptr, m12, natom )
   if (ipol > 0) then
      call adj_mem_ptr( h_ptr, m14, 15*natom )
   end if
   call adj_mem_ptr( h_ptr, m16, 2*natom )
   lasth = h_ptr
   
   !-----------------------------------------------------------------------
   !     --- Static Integer Arrays ---
   !-----------------------------------------------------------------------
   
   i_ptr = 1
   call adj_mem_ptr( i_ptr, i02, nres + 1 )
   call adj_mem_ptr( i_ptr, i04, natom )
   call adj_mem_ptr( i_ptr, i06, ntypes*ntypes )
   call adj_mem_ptr( i_ptr, i08, natom )
   call adj_mem_ptr( i_ptr, i10, 2*nnb )
   iibh = i_ptr
   
   !     ----- BOND ARRAYS -----
   
   ijbh  = iibh  + ntbond
   iicbh = ijbh  + ntbond
   iiba  = iibh  + nbonh
   ijba  = ijbh  + nbonh
   iicba = iicbh + nbonh
   i24   = iicbh  + ntbond + nbper
   
   !     ----- ANGLE ARRAYS -----
   
   i26 = i24 + ntangl
   i28 = i26 + ntangl
   i30 = i28 + ntangl
   i32 = i24 + ntheth
   i34 = i26 + ntheth
   i36 = i28 + ntheth
   i38 = i30 + ntheth
   i40 = i30 + ntangl + ngper
   
   !     ----- DIHEDRAL ARRAYS -----
   
   i42 = i40 + ntdih
   i44 = i42 + ntdih
   i46 = i44 + ntdih
   i48 = i46 + ntdih
   i50 = i40 + m7
   i52 = i42 + m7
   i54 = i44 + m7
   i56 = i46 + m7
   i58 = i48 + m7
   icnstrgp = i48 + ntdih + ndper
   itgtfitgp = icnstrgp + natom   ! VH tgtmd - 2 groups
   itgtrmsgp = itgtfitgp + natom
   ibellygp = itgtrmsgp + natom
   noshake = ibellygp + natom
   i64 = noshake + nbonh + nbona
   i_ptr = i64 + natom
   if (nmropt >= 2) then
      call adj_mem_ptr( i_ptr, i65, mxsub*isubi )
   else
      call adj_mem_ptr( i_ptr, i65, 0 )
   end if
   call adj_mem_ptr( i_ptr, i70, natom + 1 )
   
   ! Allocate memory for NMR restraints/weight changes:
   
   call adj_mem_ptr( i_ptr, inmr02, intreq )
   
   ! Allocate the IVMxx array:
   
   call adj_mem_ptr( i_ptr, iifstwt, ntbond )
   call adj_mem_ptr( i_ptr, iifstwr, nres + 1 )
   
   !  --- right now, i78 (iar1) is not being created, needs no space:
   call adj_mem_ptr( i_ptr, i78, 0 )
   
   !  --- allocate array for numbond for surface area calculation:
   
   if( gbsa == 1 )then
      call adj_mem_ptr( i_ptr, i80, natom )
      call adj_mem_ptr( i_ptr, i82, 40*natom )
   else if( gbsa == 2 )then
      call adj_mem_ptr( i_ptr, i80, natom )
      call adj_mem_ptr( i_ptr, i82, 80*natom )
   else
   !  call adj_mem_ptr( i_ptr, i80, 0 )
      call adj_mem_ptr( i_ptr, i82, 0 )
   end if
   
   if(igb /= 0 .or.hybridgb>0 ) then
      call adj_mem_ptr( i_ptr, i86, natom )
   else
      call adj_mem_ptr( i_ptr, i86, 0 )
   end if

   if (icnstph /= 0) then
      call adj_mem_ptr( i_ptr, icpstinf, TITR_RES_C*STATEINF_FLD_C)
      call adj_mem_ptr( i_ptr, icpresst, TITR_RES_C)
      call adj_mem_ptr( i_ptr, icpptcnt, TITR_STATES_C)
      call adj_mem_ptr( i_ptr, icptrsct, 1)
      call adj_mem_ptr( i_ptr, icphidx,  TITR_RES_C)
      call adj_mem_ptr( i_ptr, icptpair, TITR_RES_C * 4)
   else
      call adj_mem_ptr( i_ptr, icpstinf, 0)
      call adj_mem_ptr( i_ptr, icpresst, 0)
      call adj_mem_ptr( i_ptr, icpptcnt, 0)
      call adj_mem_ptr( i_ptr, icptrsct, 0)
      call adj_mem_ptr( i_ptr, icphidx,  0)
      call adj_mem_ptr( i_ptr, icptpair, 0)
   end if
   lasti = i_ptr
   
   !     --- crude (but useful?) estimate for MAXPR:
   ! DAN ROE: Does this need to be changed for hybridgb
   if( igb /= 0 ) then
      maxpr = 1
   else
      if( numextra == 0 ) then
         maxpr_float = natom * (cutoffnb + skinnb)**3 / 3.0d0
      else   ! need more nonbon storage with extra points
         maxpr_float = natom * (cutoffnb + skinnb)**3 / 2.5d0
      end if
      
      !       --- cap at maximum possible number of pairs:
      
      natom_float = natom
      n2_float = natom_float*(natom_float-1.d0)/2.d0
      if( maxpr_float > n2_float ) maxpr_float = n2_float
      
      !       --- check that MAXPR fits into 32 bit integer:
      
      if( maxpr_float < 2.147d9 ) then
         maxpr = maxpr_float
      else
         write(6,'(a,e12.2)' ) &
               'Unreasonably large value for MAXPR: ',maxpr_float
         call mexit(6,1)
      end if
# ifdef MPI
      if(iamoeba.eq.0.and.periodic == 1) then
!Antonios changed here
!         if( numtasks <= 8 ) maxpr = 10*maxpr/numtasks
         if( numtasks <= 8 ) maxpr = maxpr/numtasks

         !  allow for some load imbalance in list at high processor number:
!         if( numtasks >  8 ) maxpr = 10*4*maxpr/(3*numtasks)
         if( numtasks >  8 ) maxpr = 4*maxpr/(3*numtasks)
      end if
!Antonios ends here
# endif
   end if

   lastpr = maxpr
   if( igb == 0 ) then
      istartr = lastr
      istarti = i_ptr
      iendr = lastr
      iendi = i_ptr
      call ewald_mem(maxpr,natom,nnb,istartr,iendr,istarti,iendi)
      lastr = iendr
      i_ptr = iendi
      lasti = i_ptr
      
      istartr = lastr
      istarti = i_ptr
      iendr = lastr
      iendi = i_ptr
      call debug_mem(natom,ntypes, &
            istartr,iendr,istarti,iendi)
      lastr = iendr
      i_ptr = iendi
      lasti = i_ptr
      
      call init_extra_pts1()
   else
      istartr = lastr
      istarti = i_ptr
      iendr = lastr
      iendi = i_ptr
      call debug_mem(natom,ntypes,istartr,iendr,istarti,iendi)
      lastr = iendr
      i_ptr = iendi
      lasti = i_ptr
   end if
   
   return
end subroutine locmem 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine adj_mem_ptr here]
subroutine adj_mem_ptr(mem_ptr,assign_ptr,size)
   implicit none
   integer mem_ptr,assign_ptr,size

   assign_ptr = mem_ptr
   mem_ptr = mem_ptr + size
   return
end subroutine adj_mem_ptr 
