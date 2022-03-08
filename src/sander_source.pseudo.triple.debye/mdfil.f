#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdfil here]
subroutine mdfil
   
   !     Author: George Seibel; many subsequent modifications by many others.
#ifdef MPI
   use remd, only : rem, rremd, remlog, remtype, remstripcoord, saveenefile,&
                    clusterinfofile, reservoirname
#endif
   use cns_xref, only : is_xref_on

   implicit none
   
   !     Modified for multisander to allow reading command line from a string
   !     rather than the command line.  Use iargc_wrap and getarg_wrap instead
   !     of the intrinsics.
   
   !     OUTPUT: (to common)
   
#  include "files.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#  include "dprec.h"
#endif
   
   !     INTERNAL:
   
   character(len=80) arg
   !         temp for each of the whitespace delimited command line arguments
   integer iarg
   !         index of the current argument
   integer iargc_wrap
   !         wrapper to intrinsic that returns the index of the last argument
   !         from either the command line or a string
   integer last_arg_index
   !         index of the last argument
   
   !     --- default file names ---
   
   mdin   = 'mdin'
   mdout  = 'mdout'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   mdvel  = 'mdvel'
   mden   = 'mden'
   mdcrd  = 'mdcrd'
   inptraj = 'inptraj'
   mdinfo = 'mdinfo'
   vecs   = 'vecs'
   freqe  = 'dummy'
   rstdip = 'rstdip'
   inpdip = 'inpdip'
   mddip  = 'mddip'
   radii  = 'radii'
   cpin   = 'cpin'
   cpout  = 'cpout'
   cprestrt = 'cprestrt'
   evbin  = 'evbin'                                        ! EVB input file
   evbout = 'evbout'                                       ! EVB output file
   pimdout = 'pimdout'
!Antonios added. 29.1.10
   twhb = 'twhb'
   twvdw = 'twvdw'
   twchi = 'twchi'
!Antonios end
#ifdef MMTSB
   mmtsb_setup_file = 'mmtsb_setup.job'
#endif
   if (numgroup == 1) groups(:) = ' '

   !     --- default status of output: New
   
   owrite = 'N'

   facc = 'W'           ! default: overwriting
#ifdef MPI
   remlog = 'rem.log'   ! default log file name for REM
   remtype = 'rem.type' ! default log file name for REM simulation type info
   remstripcoord = ' ' ! default is no output for hybrid RREMD stripped 
                       !  trajectory unless a filename is specified.
   saveenefile = 'saveene' ! default reservoir structure energy file
   clusterinfofile = 'cluster.info' ! default dihedral cluster info file
   reservoirname = 'reserv/frame' ! default reservoir structure file name
#endif
   is_xref_on = .false.   !  default for cns_xref

   !     --- get command line arguments ---
   
   iarg = 0
   last_arg_index = iargc_wrap()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getarg_wrap(iarg,arg)

      if (arg == '-O') then
#ifdef ABSOFT_WINDOWS
         ! although Replace is standard Fortran 90 apparently Absoft f90.exe 
         ! cannot handle it
         owrite = 'U' !      status of output: unknown
#else
         owrite = 'R' !      status of output: Replace
#endif
      else if (arg == '-A') then
         owrite = 'U' !      status of output: unknown
         facc = 'A'
      else if (arg == '-i') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdin)
      else if (arg == '-o') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdout)
      else if (arg == '-p') then
         iarg = iarg + 1
         call getarg_wrap(iarg,parm)
      else if (arg == '-c') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpcrd)
      else if (arg == '-vecs') then
         iarg = iarg + 1
         call getarg_wrap(iarg,vecs)
      else if (arg == '-radii') then
         iarg = iarg + 1
         call getarg_wrap(iarg,radii)
      else if (arg == '-f') then
         iarg = iarg + 1
         call getarg_wrap(iarg,freqe)
      else if (arg == '-r') then
         iarg = iarg + 1
         call getarg_wrap(iarg,restrt)
      else if (arg == '-ref' .or. arg == '-z') then
         iarg = iarg + 1
         call getarg_wrap(iarg,refc)
      else if (arg == '-e') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mden)
      else if (arg == '-v') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdvel)
      else if (arg == '-x'.or.arg == '-t') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdcrd)
      else if (arg == '-y') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inptraj)
      else if (arg == '-inf') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdinfo)
      else if (arg == '-idip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpdip)
      else if (arg == '-rdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rstdip)
      else if (arg == '-mdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mddip)
      else if (arg == '-cpin') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cpin)
      else if (arg == '-cpout') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cpout)
      else if (arg == '-cprestrt') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cprestrt)
      else if (arg == '-evbin') then                       ! EVB input file
         iarg = iarg + 1
         call getarg_wrap(iarg,evbin)
      else if (arg == '-evbout') then                      ! EVB output file
         iarg = iarg + 1
         call getarg_wrap(iarg,evbout)
!Antonios added
      else if (arg == '-b') then
         iarg = iarg + 1
         call getarg_wrap(iarg,twhb)
      else if (arg == '-d') then
         iarg = iarg + 1
         call getarg_wrap(iarg,twvdw)
      else if (arg == '-chi') then
         iarg = iarg + 1
         call getarg_wrap(iarg,twchi)
!Antonios end

#ifdef MMTSB
      else if (arg == '-mmtsb') then
         iarg = iarg + 1
         call getarg_wrap(iarg, mmtsb_setup_file)
#endif
#ifdef MPI
      else if (arg(1:3) == '-p4') then
         iarg = iarg+1
      else if (arg == '-np') then
         iarg = iarg+1
      else if (arg == '-mpedbg') then
         continue
      else if (arg == '-dbx') then
         continue
      else if (arg == '-gdb') then
         continue

      else if (arg == '-nrecip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg,'(i5)',err=91) num_recip
         if (num_recip == numtasks) then
            num_direct=numtasks
         else
            num_direct=numtasks-num_recip
         end if
         
         !     Parse input options for multisander
         
      else if (arg == '-ng') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg,'(i5)',err=91) numgroup

      else if (arg == '-ng-nonsequential') then
         ng_sequential = .false.

      else if (arg == '-groupfile') then
         iarg = iarg + 1
         call getarg_wrap(iarg,groups)

      else if (arg == '-gpes') then
         iarg = iarg + 1
         call getarg_wrap(iarg,gpes)

      else if (arg == '-rem') then
         iarg = iarg + 1
         call getarg_wrap(iarg, arg)
         read(arg, '(i5)', err=91) rem

      else if (arg == '-nslice') then                      ! # of PIMD slices
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg, '(i5)', err=91) nslice

      else if (arg == '-remlog') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remlog)

      !     RREMD Options
      else if (arg == '-rremd') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg, '(i5)',err=91) rremd

      else if (arg == '-saveene') then
         iarg = iarg + 1
         call getarg_wrap(iarg, saveenefile)

      else if (arg == '-clusterinfo') then
         iarg = iarg + 1
         call getarg_wrap(iarg, clusterinfofile)

      else if (arg == '-reservoir') then
         iarg = iarg + 1
         call getarg_wrap(iarg, reservoirname)
      !     End RREMD options
      else if (arg == '-remtype') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remtype)

      else if (arg == '-hybridtraj') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remstripcoord)
      ! End REMD Options

#endif
#ifdef PUPIL_SUPPORT
      else if ((arg == '-ORBInitialPort') .or. &
               (arg == '-ORBInitialHost') .or. &
               (arg == '-OptPrint')) then
        iarg = iarg + 1
#endif /*PUPIL_SUPPORT*/
      else if (arg == '-pimdout') then
         iarg = iarg + 1
         call getarg_wrap(iarg,pimdout)

      else if (arg == '-cns') then
         is_xref_on = .true.
      else if (arg == ' ') then
         continue
      else
         write(6,'(/,5x,a,a)') 'mdfil: Error unknown flag: ',arg
         write(6,9000)
         call mexit(6, 1)
      end if 
   end do  !  while (iarg < last_arg_index)

   return
   
#ifdef MPI
   91 write(6,*) 'mdfil: Error "-nrecip", "-rem" and "-ng"', &
                 'require integer arguments'
   write(6,*)'                   '
   call mexit(6, 1)
#endif
   9000 format(/,5x, &
         'usage: sander  [-O|A] -i mdin -o mdout -p prmtop -c inpcrd ', &
         '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden ', &
         '-idip inpdip -rdip rstdip -mdip mddip ', &
#ifdef MPI
         '-ng numgroup -remlog rem.log -remtype rem.type -rem [0|1|2] ', &
         '-rremd [0|1|2|3] -saveene saveene -clusterinfo cluster.info ', &
         '-reservoir reserv/frame -hybridtraj hybrid.strip.crd' &
#endif
         '-ng numgroup -remlog remlog -rem [0|1|2] ', &
         '-inf mdinfo -radii radii -y inptraj]' &
         , /, 'Consult the manual for additional options.')
end subroutine mdfil 

!        '-O                Overwrite existing files.',
!        '-A                Append existing files.',
!        '-i MDIN           Namelist control input file',
!        '-o MDOUT          Output file.',
!        '-p PARM           ParmTop file.',
!        '-c INPCRD         Coordinate file.',
!        '-vecs VECS        ???',
!        '-radii RADII      ???',
!        '-f FREQE          ???',
!        '-r RESTRT         ???',
!        '-ref REFC         ???',
!        '-z REFC           alias for -ref.',
!        '-e MDEN           ???',
!        '-v MDVEL          ???',
!        '-x MDCRD          ???',
!        '-t MDCRD          alias for -x MDCRD',
!        '-inf MDINFO       ???',
!        '-y INPCRD         input trajectory for imin==5',
!        '-idip INPDIP      ???',
!        '-rdip RSTDIP      ???',
!        '-mdip MDDIP       ???',
!        '-ng NUMGROUP      number of separate sander groups',
!        '-rem [0|1|2]      type of REM simulation',
!        '-remlog REMLOG    the filename for REM log file',
!        '-rremd [0|1|2]    type of REMD reservoir',
!        '-remtype REMTYPE  the filename for REM simulation type output file',
!        '-hybridtraj TRAJOUT    the filename for hybrid stripped output traj',
!        '-saveene SAVEENE  the filename for energies of reservoir structures',
!        '-clusterinfo CLUSTERINFO    the filename with dihedral cluster info',
!        '-reservoir RESNAME the filename of reservoir structures',
!        '-cpin CPDAT       Constant pH state information ',
!        '-cprstrt CPDAT    Constant pH state restart information',
!        '-cpout CPOUT      Constant pH protonation output
!        '-evbin  EVBIN     EVB input file ',
!        '-evbout EVBOUT    EVB output file ',
!#ifdef MMTSB
!        '-mmtsb MMTSB      MMTSB Setup file; contents server generated'
!#endif
!#ifdef MPI
!        '-nrecip N     Set number of reciprocal tasks to N;'
!        '              if < numtasks, remaining tasks go to direct.'
!        '-p4 -np -mpedbg -dbx -gdb -- flags read by MPI'
!#endif



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for IARGC to support reading command line input from a string
integer function iargc_wrap()

   implicit none
   integer iargc
#ifdef MPI
#  include "files.h"
#  include "parallel.h"

   integer istart, iend
   integer ia, ie

   if (numgroup > 1 .and. groups(1:1) /= ' ') then
      ia = 0
      istart = 1
      iend = len(groupbuffer)

      do while (istart <= iend)
         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else
            do ie = istart, iend
              if (groupbuffer(ie:ie) == ' ') exit
            end do
            ia = ia + 1
            istart = ie
         end if
      end do
      iargc_wrap = ia
   else
#endif

      iargc_wrap = iargc()

#ifdef MPI
   end if
#endif
end function iargc_wrap 



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for GETARG to support grabbing command line arguments from a string
subroutine getarg_wrap(iarg, arg)

   implicit none
   integer iarg
   character(len=*) arg

   integer*4 which_argument

#ifdef MPI
#  include "files.h"
   integer istart, iend
   integer ia, ie

   if (groups(1:1) /= ' ') then

      ia = 0
      istart = 1
      iend = len(groupbuffer)

      do while (istart <= iend)
         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else
            do ie = istart, iend
              if (groupbuffer(ie:ie) == ' ') exit
            end do
            ia = ia + 1

            if (iarg == ia) then
               arg = groupbuffer(istart:ie)
               return
            end if
            istart = ie
         end if
      end do

   else
#endif

      ! Intrinsic getarg requires a 4 byte integer argument; 
      ! this guards the argument for builds with default 8 byte integers.
      which_argument = iarg
      call getarg(which_argument, arg)

#ifdef MPI
   end if
#endif
end subroutine getarg_wrap 

