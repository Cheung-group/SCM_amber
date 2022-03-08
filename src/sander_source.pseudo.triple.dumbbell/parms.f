 !<compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"

module parms


!    BC_PARMR is the number of reals in common RPARMS; BC_PARMI is the
!    number of ints in IPARMS.  (Change these if you change the sizes below).

! Antonios changed
#define BC_PARMR 2915000
#define BC_PARMI 3000

!RCW: WARNING - DO NOT DELETE THE COMMON BLOCKS HERE - THEY ARE A HACK TO FORCE
!               THE MEMORY LAYOUT TO BE LINEAR FOR DOING MPI BROADCASTS. ULTIMATELY
!               WE SHOULD GET RID OF THESE, MAKE EVERYTHING DYNAMIC AND DO THE BROADCASTS
!               IN A BETTER WAY.

integer, parameter :: num_bc_parmr = BC_PARMR, num_bc_parmi = BC_PARMI
integer, parameter :: MAX_BOND_TYPE = 5000 !NUMBND
integer, parameter :: MAX_ATOM_TYPE = 800  !NATYP

_REAL_ rk(MAX_BOND_TYPE),req(MAX_BOND_TYPE),tk(5000),teq(5000),pk(3000), &
      pn(3000),phase(3000),cn1(480000),cn2(480000),solty(480000), &
      gamc(2000),gams(2000),fmn(2000), &
      asol(480000),bsol(480000),hbcut(480000)
common/rparms/rk,req,tk,teq,pk, &
      pn,phase,cn1,cn2,solty, &
      gamc,gams,fmn, &
      asol,bsol,hbcut

integer ipn(3000)
! Antonios ends

common/iparms/ipn

logical :: charmm
integer, parameter :: MAXIMPRTYPES=1200
integer :: nimprtyp ! number of charmm improper types
_REAL_ cn114(1830),cn214(1830),rkub(900),rub(900)
_REAL_,dimension(1:MAXIMPRTYPES) :: pk_impr,phase_impr

! NPHB is the number of h-bond parameters. NIMPRP is the number of
! improper torsional parameters (NPTRA-NIMPRP is the number of regular
! torsional parameters).

integer       numbnd,numang,nptra,nphb,nimprp
common/prmlim/numbnd,numang,nptra,nphb,nimprp

end module parms
