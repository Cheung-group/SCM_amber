!+ Specification and control of Amber's Input/Output

!Antonios added
! File names
character(len=4096) groupbuffer
character(len=256) mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmrf, mincor, &
      vecs, radii, freqe,redir(9),rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file,pimdout, &
      inptraj,twhb,twvdw,twchi,fenpc,fencc

character owrite, facc
common /files/ groupbuffer, mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmrf, mincor, &
      vecs, radii, freqe, owrite, facc,rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file,pimdout, &
      inptraj,twhb,twvdw,twchi,fenpc,fencc
!Antonios end

! put this in a separate common block to stop the compiler from
! complaining about misalignment
integer numgroup, nslice
common/nmgrp/ numgroup, nslice

! File units
! An I/O Unit resource manager does not exist.
integer     MDCRD_UNIT
integer     INPTRAJ_UNIT
integer     MDEN_UNIT
integer     MDINFO_UNIT
integer     MDVEL_UNIT
!Antonios added
integer     TWHB_UNIT
integer     TWVDW_UNIT
integer     TWCHI_UNIT
integer     FENCC_UNIT
integer     FENPC_UNIT
!Antonios end
parameter ( MDINFO_UNIT =  7 )
parameter ( MDCRD_UNIT  = 12 )
parameter ( INPTRAJ_UNIT = 24 )
parameter ( MDEN_UNIT   = 15 )
parameter ( MDVEL_UNIT  = 13 )
!Antonios added
parameter ( TWHB_UNIT  = 37 )
parameter ( TWVDW_UNIT  = 36 )
parameter ( TWCHI_UNIT  = 38 )
parameter ( FENPC_UNIT  = 39 )
parameter ( FENCC_UNIT  = 40 )
!Antonios end
integer, parameter :: CNSTPH_UNIT = 18, CPOUT_UNIT = 19

! 18 was picked because CNSTPH uses it; conflicts are not expected.
integer     MMTSB_UNIT
parameter ( MMTSB_UNIT = 18 )

!!
!! EVB I/O unit
!!
   integer, parameter :: evb_unit = 75
   integer, parameter :: schlegel_unit = 80

!! FULL PIMD I/O unit
   integer, parameter :: pimd_unit = 277
! File related controls and options
character(len=80) title,title1
common/runhed/ title, title1

logical mdin_ewald,mdin_pb,mdin_amoeba

#ifdef APBS
logical mdin_apbs, sp_apbs
common/mdin_flags/mdin_ewald,mdin_pb,mdin_amoeba,mdin_apbs,sp_apbs
#else
common/mdin_flags/mdin_ewald,mdin_pb,mdin_amoeba
#endif /* APBS */

integer BC_HULP  ! size in integers of common HULP
parameter ( BC_HULP = 9 )

integer     ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave
common/hulp/ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave

!      NMRRDR : Contains information about input/output file redirection
!               REDIR and IREDIR contain information regarding
!               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
!               PCSHIFT and DIPOLE respectively. If IREDIR(I) > 0,
!               then that input/output has been redirected.

integer iredir(9)
common/nmrrdr/redir,iredir
