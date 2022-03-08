! <compile=optimized>
#include "copyright.h"
#include "assert.h"
#include "dprec.h"

module remd

! MODULE: REMD
! ================== REPLICA EXCHANGE MOLECULAR DYNAMICS ====================
! Daniel R. Roe, 2007
! Based on original implementation in multisander.f by:
!    Guanglei Cui
!    Carlos Simmerling
!
! ---=== Subroutine list ===---
!   Called from outside remd.f:
!     remd_setup (from sander.f)
!     remd_cleanup (from sander.f)
!     remd_exchange (from sander.f)
!     hybrid_remd_ene (from runmd.f)
!   Internal:
!     subrem
!     remd_scale_velo
!     sorttemp
!     load_reservoir_structure
!     stripwat
!     calc_rremd_cluster
!     load_reservoir_files
!     templookup *This is a function.

implicit none

! DAN ROE: None of this should be available to sander if no MPI - REMD only
!          works in parallel.

#ifdef MPI

! ... Files for REMD output:
!remlog_unit        - records information for each exchange.
!remtype_unit       - records general information about replica runs.
!remstripcoord_unit - records the stripped coordinates during a hybrid run.
integer, parameter :: remlog_unit=99, remtype_unit=98, remstripcoord_unit=97

! ... Constants:
!   maxgroup     - largest number of replicas allowed in a REMD run.
!   maxreservoir - largest number of structures allowed in a reservoir.
!                  NOTE: if maxreservoir is changed, the file format in
!                  sander.f also needs to be changed.
!   maxdihclust  - maximum number of dihedral angles that can be binned
!                  (each unique combination of dihedral bins defines a cluster)
!                  for a RREMD run with defined weights (rremd==3).
!   maxclust     - maximum number of clusters that can be defined.
integer, parameter :: maxgroup=1024, maxreservoir=200000, maxdihclust=100, &
                      maxclust=100000

! ... logical variables:
!   DAN ROE: initmodwt, allocatepbmem removed
!   hybridwritetraj - if true a trajectory of stripped coords will be written
!                     during hybrid remd.
!   jumpright       - if true replica that controls the exchange will attempt
!                     to exchange with the replica at next highest temperature.
!   initmodwt       - If nmropt>0 modwt will reset temp0 to its initial value
!                     after each exchange. initmodwt will force a re-read of
!                     temp0 after each exchange.
logical, save :: jumpright, hybridwritetraj, initmodwt

! ... integer variables:
! mdloop      - Current exchange index. mdloop will be always 0 unless rem>0
! rem         - The type of replica exchange
!               0, regular MPI/MD
!               1, conventional REM
!               2, partial REM
! rremd       - The type of reservoir replica exchange
!               0, no reservoir
!               1, Boltzmann weighted reservoir
!               2, 1/N weighted reservoir
!               3, reservoir with weights defined by dihedral clustering
! DAN ROE: reserv_velo should probably be a logical
! reserv_velo - 0 if no velocities in reservoir restart files
!               1 if we have velo
! exchsuccess - The number of successful exchanges between each neighbor pair
!               for printing acceptance %
! repnum      - Replica number
! numreps     - Total # of replicas, set equal to numgroup
integer, save :: mdloop, rem, reserv_velo, rremd, repnum, numreps
integer, dimension(:), allocatable, save :: exchsuccess

! ... floats:
! mytemp        - collected temperature of current replica
! myeptot       - collected potential energy of current replica
! mytargettemp  - collected target temperature of current replica
! myscaling     - the velocity scaling factor that will be applied next
! newtargettemp - the updated target temperature for the next sander run
! remd_ekmh     - Store ekmh variable between runmd calls (exchanges). ekmh is
!                 the KE from old velocities.
! temptable()   - Store sorted list of replica temperatures.
! DAN ROE: common block not allowed to have save attribute
_REAL_ :: mytemp, myeptot, mytargettemp, myscaling, newtargettemp,remd_ekmh
! DAN ROE: Is it ok to make remd_ekmh part of the common block?
common/remf/mytemp, myeptot, mytargettemp, myscaling, newtargettemp, remd_ekmh
_REAL_, dimension(:), allocatable, save :: temptable

! ... strings:
! remlog          - replica exchange log filename
! remtype         - replica exchange simulation type info filename
! remstripcoord   - stripped coordinates during a hybrid rremd run.
! saveenefile     - contains energies of all structures in reservoir
! clusterinfofile - contains information for dihedral cluster binning (rremd==3)
! reservoirname   - base filename of reservoir files
! DAN ROE: Is allrgoupbuffer still needed?
! allgroupbuffer  - all of the group buffers so we can store and exchange the
!                   input file names
character (len=80), save :: remlog, remtype, remstripcoord, saveenefile, &
                            clusterinfofile, reservoirname

! ---=== RREMD ===---
! saveene()     - Energy of each structure in the reservoir. 
! restemp0      - reservoir temperature
! rremd_idx     - # of random structure from reservoir
! reservoirsize - total # of structures in reservoir
_REAL_, dimension(:), allocatable, save :: saveene
_REAL_, save :: restemp0
integer, save :: rremd_idx, reservoirsize

! ---=== Dihedral Clustering ===---
! clusternum() - cluster that reservoir structure belongs to
! nclust       - # clusters
! nclustdih    - # dihedrals used for clustering
! clusterid    - array of bin values for each dihedral, 1 set per cluster
! clustersize  - # structures from reservoir in this cluster (note 0 is unknown
!                cluster, currently set to 0 in multisander.f)
! dihclustat   - 2d, 4 atoms that define each of the dihedrals used for clustering
! dihclustnbin - #bins used for each of the dihedrals (each is dihclustnbin/360
!                in width)
! currdihid    - dih bin values for the current MD structure (to assign it to a
!                cluster)
! incluster    - cluster for the current MD structure
integer, dimension(maxreservoir), save :: clusternum
integer, dimension(maxclust,maxdihclust), save :: clusterid
integer, dimension(0:maxclust), save :: clustersize
integer, dimension(maxdihclust), save :: currdihid, dihclustnbin
integer, dimension(4,maxdihclust), save :: dihclustat
integer, save :: incluster,nclust,nclustdih

! Hybrid REMD Temp. Storage
_REAL_, dimension(:), allocatable, private :: hybrid_coord
_REAL_, dimension(:), allocatable, private :: hybrid_force

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! The following subroutines are used for REMD

!*********************************************************************
!               SUBROUTINE REMD_SETUP
!*********************************************************************
!  Set up REMD run - open log, set up temperature table, etc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine remd_setup(numexchg,hybridgb,numwatkeep,temp0,mxvar,natom)

   implicit none
#  include "parallel.h"
#  include "mpif.h"

   integer, intent(in) :: numexchg,hybridgb,numwatkeep,mxvar,natom
   _REAL_, intent(in) :: temp0

   integer i, j, ierror

!--------------------

   ! ---=== INITIALIZE VARIABLES ===--- 
   !mytemp = 0.0d0
   myeptot = 0.0d0
   mytargettemp = 0.0d0
   myscaling = -1.0d0
   newtargettemp = -99.0d0
   jumpright = .true.
   initmodwt = .true.
   rremd_idx=-1
   ! For storing ekmh between runmd calls
   ! DAN ROE: ntt==1 only?
   remd_ekmh=0.0d0

   ! Allocate memory for temptable and exchsuccess
   allocate(temptable(numreps), exchsuccess(numreps), stat=ierror)
   REQUIRE(ierror==0)

   ! initialize the exchange success array to zero (for printing acceptance
   ! rate exchsuccess(i) give the number of times that temperature i 
   ! has exch with i+1
   do i=1,numreps
       exchsuccess(i)=0
   enddo

   ! ---=== OPEN REMD FILES ===---
   !  Open remtype_unit to record initial REMD information
   if (worldrank == 0 .and. rem>0) then
      call amopen(remtype_unit,remtype,'U','F','W')
   endif

   ! R-REMD: load in the file containing energies for the reservoir.
   ! Only sander masters perform the read since they do the exchanges.
   ! Coordinates/velocities for reservoir are loaded as needed during the
   ! RREMD run.
   if (rremd > 0) call load_reservoir_files()

   ! Open and write out header to remlog file. Only overall master
   !  deals with the remlog.
   if (worldrank == 0) then
      write(6,'(a)') "REMD: Initializing remlog."
      call amopen(remlog_unit,remlog,'U','F','W')
      write(remlog_unit,'(a)') "# Replica Exchange log file"
      write(remlog_unit,'(a,i10)') "# numexchg is ",numexchg
      if (numwatkeep >= 0) then
         write(remlog_unit,'(a)') "# HYBRID REMD:"
         write(remlog_unit,'(a,i3,a)') "#   using igb= ",hybridgb,&
                                       " for exchange energies."
         write(remlog_unit,'(a,i10)') &
            "#   number of closest waters kept= ",numwatkeep
      endif
      if (rremd>0) write(remlog_unit,'(a,i3)') "# Reservoir REMD: ",rremd
      ! Write out REMD filenames
      write(remlog_unit,'(a)') "# REMD filenames:"
      write(remlog_unit,'(a,a)') "#   remlog= ",trim(remlog)
      write(remlog_unit,'(a,a)') "#   remtype= ",trim(remtype)
      if (rremd>0) then
         write(remlog_unit,'(a,a)') "#   saveene= ",trim(saveenefile)
         if (rremd==3) &
            write(remlog_unit,'(a,a)') "#   clusterinfo= ",trim(clusterinfofile)
         write(remlog_unit,'(a,a)') "#   reservoir= ",trim(reservoirname)
      endif
      write(remlog_unit,'(a)') &
         "# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0,&
         & Success rate (i,i+1), ResStruct#"
   end if ! worldrank==0, remlog header

   ! ---=== HYBRID REMD SETUP ===---
   hybridwritetraj=.false.
   if (numwatkeep>=0) then
      ! 1a- Allocate memory for temp. coord/force storage
      ! Note: This should be the same as in locmem.f except am_nbead is
      !  not known. amoeba may not work with hybrid remd.
      !   if (sanderrank==0) write(6,'(a)') "HYBRID REMD: Allocating memory."
      allocate( hybrid_coord(3*natom + mxvar), &
                hybrid_force(3*natom + mxvar + 40), stat=ierror)
      REQUIRE( ierror == 0 )
      
      ! If an output file was specified for hybrid REMD stripped
      !  coords trajectory (-hybridtraj FILE), open it. Only masters
      !  will writ to it.
      if (sanderrank==0) then
         if (remstripcoord(1:1).eq.' ') then
           write(6,'(a)') &
              "HYBRID REMD: Hybrid stripped traj file will not be written."
         else
            !write(remstripcoord,'(a16,i3.3)') "hybrid.stripped.", nodeid
            hybridwritetraj=.true.
            write(6,'(a,a)') &
               "HYBRID REMD: Opening hybrid stripped coord output: ",&
               remstripcoord
            write(remlog_unit,'(a,a)') "#   hybdridtraj= ",remstripcoord
            call amopen(remstripcoord_unit,remstripcoord,'U','F','W')
            write(remstripcoord_unit,'(a80)') "Hybrid stripped coords"
         endif
      endif
   endif
  
   call mpi_barrier(commworld, ierror)

   ! ---=== MAKE SORTED TEMPERATURE TABLE ===---
   ! Now make the sorted temp0 table, used for exchanges. Only masters
   !  need this since they do the exchanges.
   if (sanderrank==0) then
      ! Master processes gather all of the target temperatures
      ! DAN ROE: Since mytargettemp isn't set yet gather temp0;
      !  for rem==2 would need to gather temp0les
      ! Could probably just set mytargettemp=temp0 or
      !  mytargettemp=temp0les and then gather mytargettemp
      !write(6,*) "Sending temp0= ",temp0
      call mpi_allgather(temp0, 1, mpi_double_precision, &
                         temptable, 1, mpi_double_precision, &
                         commmaster, ierror)

      ! Sort temperatures
      call sorttemp(temptable)

      ! Check for duplicate temperatures 
      ! DAN ROE: Since table is already sorted, we can just check if 
      !          any two adjacent temperatures match.
      do i=1,numreps-1
         j=i+1
         ! DAN ROE: debug
         !write(6,'(a3,i3,a7,f6.2,a4,i3,a7,f6.2)') &
         !   "i= ",i," temp= ",temptable(i), &
         !   " j= ",j," temp= ",temptable(j)  
         if (temptable(i).eq.temptable(j)) then
            write (6,*) "================================"
            write (6,*) "TEMPERATURES OF 2 REPLICAS MATCH!"
            write (6,*) "T= ",temptable(i)
            write (6,*) "NOT ALLOWED FOR TEMPERATURE REMD"
            write (6,*) "================================"
            call mexit(6,1)
         endif
      enddo

   endif ! sanderrank==0, making sorted temptable

   return

end subroutine remd_setup


!*********************************************************************
!               SUBROUTINE REMD_CLEANUP
!*********************************************************************
!  Close REMD files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine remd_cleanup()

   implicit none
#  include "parallel.h"

   integer ier

!--------------------

   ! Close remlog and remtype files
   if (worldrank==0) then
      close(remlog_unit)
      close(remtype_unit)
   endif

   ! Deallocate temptable, exchsuccess. These are always allocated
   !  during a replica run.
   deallocate(temptable, stat=ier)
   REQUIRE(ier==0)
   deallocate(exchsuccess, stat=ier)
   REQUIRE(ier==0)

   ! Deallocate hybrid remd temp. storage if it was allocated
   if (allocated(hybrid_coord)) deallocate(hybrid_coord, stat=ier)
   REQUIRE(ier==0)
   if (allocated(hybrid_force)) deallocate(hybrid_force, stat=ier)
   REQUIRE(ier==0)

   if (allocated(saveene)) deallocate(saveene, stat=ier)
   REQUIRE(ier==0)

   ! Close the hybrid strip coordinates traj. file
   if (sanderrank==0.and.hybridwritetraj) close(remstripcoord_unit)

   return

end subroutine remd_cleanup


!*********************************************************************
!               SUBROUTINE REMD_EXCHANGE
!*********************************************************************
!  Perform actions necessary to exchange replicas in REMD
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine remd_exchange(x,v,nr3,natom,nr,temp0)

   implicit none

#  include "parallel.h"
#  include "mpif.h"

   _REAL_, dimension (*), intent(inout) :: x, v
   integer, intent(in) :: nr3, natom, nr
   _REAL_, intent(inout) :: temp0

   integer ierror

!-------------------

   ! ---=== RREMD DIHEDRAL CLUSTERING CALC.===---
   ! For rremd==3 calculate cluster of current structure.
   if (sanderrank==0 .and. rremd==3) then
      call calc_rremd_cluster(x)
      !write (6,*) "Current MD structure was assigned to cluster ",incluster
   endif

   ! ---=== REMD EXCHANGE CALCULATION ===---
   ! Attempt Exchange with subrem().
   ! The subrem call uses the energy from the previous runmd call.
   ! All procs call subrem: only masters do exchg but the others 
   !  keep their random # generator in sync.
   ! DAN ROE: Debug info
   !if (master) then
   !      write (6,*) "calling subrem for REMD exchange ",mdloop
   !      write (6,*) "current temp0 is ",temp0
   !      write (6,*) "myeptot is ",myeptot
   !      write (6,*) "mytargettemp is ",mytargettemp
   !endif
   call subrem()

   ! ---=== REMD COMMUNICATION ===---
   ! Broadcast reservoir structure index. All threads needs to know
   !  this to decide if they have to receive coords or not.
   if (rremd>0) call mpi_bcast(rremd_idx,1,mpi_integer,0,commsander,ierror)
   ! Broadcast mytemp common block (has mytemp, myeptot, mytargettemp,
   !  myscaling, and newtargettemp).
   call mpi_bcast(mytemp,5,mpi_double_precision,0,commsander,ierror)

   ! ---=== RREMD RESERVOIR LOADING ===---
   ! If we exchanged with the reservoir, swap the coordinates
   !  here (after the inpcrd have been read).
   if (rremd_idx>0) &
      call load_reservoir_structure(x,v,nr3,natom)

   ! ---=== REMD VELOCITY SCALING ===---
   ! REMD: If an attempt is accepted, set the target temperature to 
   !  the new one and rescale velocities.
   ! DAN ROE: Eventually take out debug info.
   call remd_scale_velo(v,temp0,nr,nr3)

   ! initmodwt forces modwt() in nmr.f to re-read values such as temp0 - 
   ! otherwise they will be reset to the initial value after each 
   ! exchange.
   initmodwt = .true.

   return

end subroutine remd_exchange


!*********************************************************************
!               SUBROUTINE SUBREM
!*********************************************************************
!  calculation of exchange probability for REMD
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine subrem()

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only: &
      ncsu_on_delta => on_delta, ncsu_on_exchange => on_exchange
#endif /* DISABLE_NCSU */


   implicit none

#  include "mpif.h"
#  include "parallel.h"

   _REAL_ delta, straw, o_scaling, metrop
   _REAL_ l_temp0, r_temp0, o_temp0
   _REAL_ l_eptot, r_eptot, o_eptot 
   _REAL_ l_temp,  r_temp,  o_temp
   _REAL_ alltempi(maxgroup), alltemp0(maxgroup), alleptot(maxgroup)
   ! For use in remlog output
   _REAL_ d_scaling(maxgroup), d_o_temp0(maxgroup)
   ! DAN ROE: for recording RREMD structure output in log
   integer d_rremd_idx(maxgroup) 
   ! DAN ROE: cluster sizes for rremd==3
   integer myclustersize, o_clustersize
   ! index is the position of the replica's T in the sorted T list
   integer myindex, l_index, r_index, o_index
   ! repnum is the actual replica number
   integer l_repnum,r_repnum,o_repnum

   integer i, ierror, istatus(mpi_status_size)

   ! temporary for exchange success, will allgather from all nodes then 
   ! add sum to actual exchsuccess array. also define an array for allgather sum
   integer texchsuccess(maxgroup),tempsuccess(maxgroup)
   _REAL_ exchfrac(maxgroup)

   logical exchange

#ifndef DISABLE_NCSU
   _REAL_ U_mm, U_mo, U_om, U_oo, beta_m, beta_o
#endif /* DISABLE_NCSU */

! ---------------------

   ! Initialize variables
   do i=1,numreps
      texchsuccess(i)=0
   enddo

   delta = 0.0d0
   straw = 0.0d0
   metrop= 0.0d0
   o_scaling = -0.1d0
   i = 0
   myindex = 0
   l_index = 0
   r_index = 0
   o_index = masterrank
   myclustersize=0
   o_clustersize=0

   ! Call amrand to get a random# which can eventually be used as a structure
   ! index in rremd. Do it here so the random # generators on all the
   ! threads stay in sync. Only the master of the highest T replica during
   ! jumpright will ever use this.
   ! Calling this even when we have no rremd will allow random # gen to stay
   ! in sync between remd and rremd runs.
   call amrand(straw)
   ! DAN ROE: FileDebug
   !   write(50+worldrank,'(i6,a,E16.6)') mdloop," straw= ",straw

   ! Only the sander masters do the exchanges.
   if (sanderrank == 0) then

      ! Gather current temperature, target temperature, and PE
      call mpi_allgather(mytemp, 1, mpi_double_precision, &
                         alltempi, 1, mpi_double_precision, &
                         commmaster, ierror)
      call mpi_allgather(mytargettemp, 1, mpi_double_precision, &
                         alltemp0, 1, mpi_double_precision, &
                         commmaster, ierror)
      call mpi_allgather(myEptot, 1, mpi_double_precision, &
                         alleptot, 1, mpi_double_precision, &
                         commmaster, ierror)

      ! alltemp0 and alleptot are indexed by replica #
      ! temptable is SORTED

      ! Find our position in the sorted temperature list
      myindex = templookup(mytargettemp, temptable)

      ! DAN ROE: Why do we look for both partners?

      ! Find neighbor Ts for exchange using Temperature list indexed
      !  by TEMPERATURE.
      ! Wrap around so replicas exchange in a circle, 
      !  i.e. highest and lowest T may attempt exchange.
      l_index = myindex - 1
      if (l_index < 1) l_index = numreps
      r_index = myindex + 1
      if (r_index > numreps) r_index = 1
      l_temp0 = temptable(l_index)
      r_temp0 = temptable(r_index)

      ! Now find replicas that have these Ts using Temperature list
      !  indexed by REPLICA
      l_repnum = templookup(l_temp0, alltemp0) 
      r_repnum = templookup(r_temp0, alltemp0) 

      ! Get the energies of these neighbor replicas
      l_eptot = alleptot(l_repnum)
      r_eptot = alleptot(r_repnum)

      ! Set up partner replica information.
      ! By definition, only even replica #s will initiate the exchange.
      if(mod(myindex, 2) == 0) then
         ! this replica will calculate delta
         if(jumpright) then
            ! Partner is to the right, i.e. at higher T
            ! DAN ROE: Integrate check for reservoir here?
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         else
            ! Partner is to the left, i.e. at lower T
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         end if
      else 
         ! This replica will not calculate delta, right and left are
         !  switched compared to controlling replica.
         ! Only needs o_eptot for writing debug info
         if(jumpright) then
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         else
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         end if
      end if  ! (even/odd rank check)

#ifndef DISABLE_NCSU
      call ncsu_on_delta(o_repnum - 1, &
         mod(myindex, 2).eq.0, U_mm, U_mo, U_om, U_oo)
#endif /* DISABLE_NCSU */

      ! RREMD: If jumpright and this replica has highest T, then attempt an
      !  exchange with a random structure in the reservoir.
      ! Read coordinates from the file when we are in sander.
      ! Initialize rremd_idx, so that it is -1 for any replica other than the 
      !  one that will use the reservoir coordinates.
      ! Change the partners for RREMD: lowest does not exchange, and highest 
      !  uses the reservoir.
      if (rremd > 0) then
         rremd_idx=-1
         if (jumpright) then
            if (myindex.eq.1) then
               ! Lowest T, will not do exchange or even wait for results from
               !  highest T
               o_index=-999
            elseif (myindex.eq.numreps) then
               ! Highest T, will exchange with reservoir
               o_temp0 = restemp0
               o_index = -900
               ! Get random structure index
               ! amrand has already been called at the start of the routine so
               !  that all threads stay in sync.
               !call amrand(straw)
               rremd_idx = straw*reservoirsize+1
               ! DAN ROE: Random number should not be outside reservoir size!
               ! By definition it wont be unless something breaks in amrand()
               if (rremd_idx.gt.reservoirsize) rremd_idx = reservoirsize
               ! Get PE, read in previously from saveene
               o_eptot=saveene(rremd_idx)
            endif 
         endif ! jumpright
      endif ! rremd check
      
      ! Calculate the exchange probablity if even # replica.
      if(mod(myindex, 2) == 0) then
         ! RREMD: Here we calculate delta beta term for normal exchanges but
         !  only beta for exchanges with non-Boltzmann weighting (rremd>1)
         ! Reservoir exchanges are for top replica and jumpright only!
         ! Note that 503.01 is 1/kB in internal amber units.
         if (rremd.eq.2.and.myindex.eq.numreps.and.jumpright) then
            ! No delta beta, assign weight of 1/N to each structure in
            !  reservoir. The exchange criterion becomes:
            !  metrop = exp[-beta_replica * (E_reservoir-E_replica)]
            ! NB 1/N delta:
            delta = (o_eptot - myEptot) * 503.01d0 / mytargettemp
            !metrop = exp(-delta)
         elseif (rremd.eq.3.and.myindex.eq.numreps.and.jumpright) then
            ! No delta beta, weight of each structure in the reservoir is
            !  user defined and has been read previously.
            ! runmd has determined which cluster the MD structure is in
            !  (incluster).
            ! Weights of each structure are obtained by clustersize(clusternum(i))
            !  or for the MD structure clustersize(incluster)
            ! Note that MD structures not present in the reservoir are given a
            !  weight of zero currently. In principle we could add to the reservoir,
            !  but that is for a future release.
            ! This makes the exchange criterion:
            !  metrop =  w_rep / w_res * exp[-beta_rep * (E_res-E_rep)]
            ! NB M/N delta:
            delta = (o_eptot - myEptot) * 503.01d0 / mytargettemp
            !metrop = exp(-delta)
            myclustersize=clustersize(incluster)
            o_clustersize=clustersize(clusternum(rremd_idx))
            !write (6,*) "myclustersize is ",myclustersize
            !write (6,*) "o_clustersize is ",o_clustersize
            !write (6,*) "metrop weight factor is ", &
            !   myclustersize/o_clustersize
            !metrop = metrop * myclustersize / o_clustersize
         else
            ! Std REMD, or RREMD with Boltzmann weighted reservoir, or NB RREMD
            !  but not exchanging with the reservoir this time.
            delta = (myEptot - o_eptot) * (mytargettemp - o_temp0) * 503.01d0 &
                  / (mytargettemp * o_temp0)
            !metrop = exp(-delta)
         endif ! REMD exchange calculation

#ifndef DISABLE_NCSU
         ! from Y.Sugita at al (JCP v=113 p=6042 year=2000)
         beta_m = 503.01D0/myTargetTemp
         beta_o = 503.01D0/o_temp0
         delta = delta + beta_m*(U_mo - U_mm) - beta_o*(U_oo - U_om)
#endif /* DISABLE_NCSU */

         metrop = exp(-delta)
         if (rremd.eq.3.and.myindex.eq.numreps.and.jumpright) &
            metrop = metrop * myclustersize / o_clustersize

         ! Get random number between 0 and 1
         call amrand(straw)
         ! DAN ROE: FileDebug
!         write(50+worldrank,'(i6,a,E16.6)') mdloop," straw= ",straw

!         if(delta < 0.0d0 ) then
!            if (debugremd) write (41+worldrank,*) &
!               "energy of higher T structure is lower, &
!               &exchange is successful"
!            write (6,*) "energy of higher T structure is lower, &
!                        &normal metropolis would accept"
!         endif 

         ! Check for exchange
         exchange=.false.
         if (straw < metrop) then
            ! Accept
            exchange=.true.
         else
            ! Reject
            exchange=.false.
         endif

         ! If exchanged, set scaling, otherwise dont
         if (exchange) then
            ! Set velocity scaling   
            myscaling = sqrt(o_temp0 / mytargettemp)
            o_scaling = 1.0d0 / myscaling
            ! If RREMD and exchanging with the reservoir we are changing 
            !  structure, not temperature, so myscaling is inverted.
            ! Dont scale if we aren't reading velocities in the reservoir.
            if (rremd.gt.0.and.myindex.eq.numreps.and.jumpright) then
               if (reserv_velo==1) then
                  myscaling = o_scaling
               else
                  myscaling=1.d0
               endif
            endif
            ! Increment exchsuccess for lower replica #
            ! Use index, not repnum, since this is a property of temperatures
            if (jumpright) then
               ! this is the lower temperature since this is the controlling rep
               ! NOTE THAT RANKS ARE 0 TO NUMREPS-1 BUT USE 1 TO NUMREPS FOR 
               ! SUCCESS	
               texchsuccess(myindex) = 1
            else
               ! other is the lower temperature
               texchsuccess(myindex-1) = 1
            endif
         else
            ! No exchange
            myscaling = -1.0d0
            o_scaling = -1.0d0
         end if ! exchange

         !write (6,*) "done with exchange calc"

! CARLOS: DO WE NEED BARRIER?

         call mpi_barrier(commmaster, ierror)

         ! send the results to the partner
         ! NOTE THAT WE NEED TO SUBTRACT 1 FROM REPNUM SINCE SUBREM USES
         ! INDEX OF 1 TO NUMREPS WHILE THE MPI PROCESSES ARE 0 TO NUMREPS-1

         ! RREMD: Dont communicate if we exchanged with reservoir
         if (rremd.eq.0.or.o_index.gt.0) then
            !write (6,*) "sending o_scaling of ",o_scaling,&
            !            "to replica ",o_repnum
            call mpi_send(o_scaling, 1, mpi_double_precision, &
                          o_repnum-1, 0, commmaster, ierror)
         endif
      else 
         ! Not the replica controlling the exchange
         !write (6,*) "NOT controlling exchange",o_index
         ! call rand to keep in sync with replicas that calculated Metropolis
         call amrand(straw)
         ! DAN ROE: FileDebug
!         write(50+worldrank,'(i6,a,E16.6)') mdloop," straw= ",straw

! CARLOS: NEED BARRIER?
         call mpi_barrier(commmaster, ierror)
         ! receive the scaling data; if it is >0 we know exchange succeeded
         ! RREMD: Don't communicate if our partner is exchanging with
         !  the reservoir
         if (rremd.eq.0.or.o_index.gt.0) then 
            call mpi_recv(myscaling, 1, mpi_double_precision, &
                          o_repnum-1, 0, commmaster, istatus, ierror)
            !write (6,*) "received scaling data ",myscaling,&
            !            "from replica ",o_repnum
         else
            ! RREMD and lowest T replica, partner attempted exchange with reservoir
            ! DAN ROE: Why the period?
            myscaling=-1.
         endif
      end if ! replica controlling the exchange (even #)

      ! toggle exchange direction
      
      jumpright = .not. jumpright
      
      if(myscaling < 0.0d0) then
         newtargettemp = mytargettemp
         exchange=.false.
         ! If RREMD, didnt exchange with the bath so get rid of the rremd_idx 
         !  that we set. Ok to do for all replicas.
         ! DAN ROE: DO this later
         !rremd_idx=-1
      else 
         ! RREMD: don't change temp0 if we exchanged with the reservoir.
         !  Instead, change coordinates to those of the new structure.
         ! HOW? Write over restart?
         ! OK to leave my_scaling since technically the new structure needs
         !  scaling.
         if (rremd.eq.0.or.o_index.gt.0) then
            newtargettemp=o_temp0
            exchange=.true.
         ! DAN ROE: No NCSU for RREMD?
#ifndef DISABLE_NCSU
            call ncsu_on_exchange(o_repnum - 1)
#endif /* DISABLE_NCSU */
         else ! RREMD highest replica
            newtargettemp = mytargettemp
            exchange=.true.
            !write (6,*) "exchanged with the structure reservoir!",rremd_idx
         endif
      endif

! CARLOS: NEED BARRIER?
      call mpi_barrier(commmaster, ierror)
 
      ! REM: gather exchange log data

      call mpi_gather(myscaling, 1, mpi_double_precision, &
                      d_scaling, 1, mpi_double_precision, &
                      0, commmaster, ierror)
      call mpi_gather(newtargettemp, 1, mpi_double_precision, &
                      d_o_temp0, 1, mpi_double_precision, &
                      0, commmaster, ierror)

      ! DAN ROE: RREMD: gather reservoir structure data
      call mpi_gather(rremd_idx,1,MPI_INTEGER, &
                      d_rremd_idx, 1, MPI_INTEGER, &
                      0, commmaster, ierror)

      call mpi_allreduce(texchsuccess,tempsuccess,numreps, &
               MPI_INTEGER,mpi_sum,commmaster,ierror)

      ! add the current successes to overall array

      do i=1,numreps
         exchsuccess(i)=exchsuccess(i)+tempsuccess(i)

         ! multiple fraction by 2.0 since we attempt this pair every OTHER 
         ! exchange attempt (alternating directions of exchange)

         exchfrac(i) = float(exchsuccess(i)) / float(mdloop) * 2.0
      enddo
            

      if(worldrank == 0) then

      ! only overall master writes the log
      ! DAN ROE: added d_rremd_idx write

         write(unit = remlog_unit, fmt = '(a,i8)') '# exchange ', mdloop
         do i = 1, numreps
            write(remlog_unit, '(i2, 6f10.2, i8)') &
               i, & ! Replica #
               d_scaling(i), & ! scaling factor
               alltempi(i), & ! current temperature
               alleptot(i), & ! current potential energy
               alltemp0(i), & ! current target temperature
               d_o_temp0(i),& ! next target temperature
               exchfrac(templookup(alltemp0(i),temptable)),& ! current exchange success fraction
                d_rremd_idx(i) ! structure# attempted from  reservoir
         end do
         ! Write out reservoir information.
         !if (rremd>0) then
         !
         !endif
      end if

      ! DAN ROE: All sander masters write REMD exchange info
      !if (rremd>0) then
         write(6,'(26("="),a,26("="))') "REMD EXCHANGE CALCULATION"
         write(6,'(a6,i10,a8,i1)') "Exch= ",mdloop," RREMD= ",rremd
         ! This Replica Information
         write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
           "Replica         "," Temp= ",mytargettemp," Indx= ",myindex, &
           " Rep#= ",repnum," EPot= ",myeptot
         ! Partner Information
         ! use .not.jumpright since jumpright has already been toggled
         if (rremd>0.and..not.jumpright.and.myindex.eq.numreps) then
            ! Partner is Reservoir
            write(6,'(a16,a7,f6.2,a10,i8,a7,f10.2)') &
               "Reservoir       "," Temp= ",restemp0," Struct#= ",rremd_idx, &
               " EPot= ",o_eptot
            if (exchange) then
            write(6,'(a20)') "ReservoirExchange= T"
            else
               write(6,'(a20)') "ReservoirExchange= F"
            endif
            if (rremd==3) then
               ! Reservoir has weights
               write(6,'(a11,i10)') "mycluster= ",incluster
               write(6,'(a15,i10,a16,i10)') &
                  "myclustersize= ",myclustersize,&
                  " o_clustersize= ",o_clustersize
            endif
         else if (rremd>0.and..not.jumpright.and.myindex.eq.1) then
            ! Lowest T, partner would be highest T but that is exchanging
            ! with the reservoir so lowest T has no partner.
            write(6,'(a)') &
               "No partner, highest T exchanging w/ Reservoir."
         else
            ! Partner is normal replica 
            write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
               "Partner         "," Temp= ",o_temp0," Indx= ",o_index, &
               " Rep#= ",o_repnum," EPot= ",o_eptot
         endif ! not jumpright and myindex.eq.numreps
         ! Exchange calculation information
         if (mod(myindex,2)==0) then
            ! This replica controlled exchange and calculated metrop
            write(6,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "Metrop= ",metrop," delta= ",delta," o_scaling= ",o_scaling
         else
            ! This replica did not control exchange
            write(6,'(a)') "Not controlling exchange."
         endif ! mod(myindex,2)
         ! Write random #, scaling, and success
         write(6,'(a8,E16.6,a12,f10.2,a10,L1)' ) &
            "Rand=   ",straw," MyScaling= ",myscaling," Success= ",exchange
         write(6,'(24("="),a,24("="))') "END REMD EXCHANGE CALCULATION"
      !endif ! rremd>0, RREMD information writeout
      
      ! If no exchange occured reset rremd_idx
      if (.not.exchange) rremd_idx=-1

   else  ! not part of commmaster

! call rand to keep rand generator in sync with masters 
! (they called it for exch)
         call amrand(straw)
         ! DAN ROE: FileDebug
!         write(50+worldrank,'(i6,a,E16.6)') mdloop," straw= ",straw
   end if

   return
end subroutine subrem


!*********************************************************************
!               SUBROUTINE REMD_SCALE_VELO
!*********************************************************************
! Scale velocities based on new temps after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine remd_scale_velo(v,temp0,nr,nr3)

   implicit none
#  include "parallel.h"
#  ifdef LES
#     include "les.h"
#  endif

   _REAL_, dimension(*), intent(inout) :: v
   _REAL_, intent(inout) :: temp0
   integer, intent(in) :: nr,nr3

   integer i

!--------------------

   ! REMD: If an attempt is accepted, set the target temperature to 
   !  the new one and rescale velocities.
   ! DAN ROE: Eventually take out debug info.
   if (rem==1) then
      if (sanderrank==0) &
         write (6,'(2(a,f6.2))') &
            "REMD: checking to see if bath T has changed: ",&
            temp0,"->",newtargettemp
      ! All processes set temperature. newtargettemp is set in subrem
      if (newtargettemp>0.0) temp0=newtargettemp
      if (myscaling > 0.0) then
         ! All processes scale velocities.
         ! DAN ROE: This could potentially be divided up as in runmd
         !  since when there are mutiple threads per group each thread 
         !  only ever knows about its own subset of velocities anyway.
         if (sanderrank==0) then
            write (6,'(a,f8.3,a,f8.3)') &
               "REMD: scaling velocities by ",myscaling,&
               " to match new bath T ",temp0
         endif
         do i = 1,nr3
            v(i) = v(i) * myscaling
         enddo
      endif
#  ifdef LES
   elseif (rem==2) then
      if (newtargettemp>0.0) temp0les=newtargettemp
      if (myscaling>0.0) then
         do i=1,nr
            if (cnum(i) > 0.0) then
               v(3*i-2) = v(3*i-2) * myscaling
               v(3*i-1) = v(3*i-1) * myscaling
               v(3*i)   = v(3*i)   * myscaling
            endif
         enddo
      endif
#  endif  /* LES */
   endif ! rem==1, velocity scaling and bath T

   return

end subroutine remd_scale_velo


!*********************************************************************
!               FUNCTION TEMPLOOKUP
!*********************************************************************
! lookup temp in templist and return its index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function templookup(temp, templist)

   implicit none
#  include "parallel.h"

   _REAL_, intent(in) :: temp
   _REAL_, dimension(maxgroup), intent(in) :: templist

   integer i
   
   do i = 1, numreps
      if(abs(temp-templist(i)) < 1.0d-6) then
         templookup = i
         return
      end if
   end do
end function templookup


!*********************************************************************
!               SUBROUTINE SORTTEMP
!*********************************************************************
! sort temp ascendingly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sorttemp(temp)

   implicit none

#  include "parallel.h"

   _REAL_, dimension(maxgroup), intent(inout) :: temp

   _REAL_ tempt
   integer i, j

   do i = 1, numreps
      do j = i + 1, numreps
         if(temp(j) < temp(i)) then
            tempt = temp(i)
            temp(i) = temp(j)
            temp(j) = tempt
         end if
      end do
   end do
end subroutine sorttemp

!*********************************************************************
!               SUBROUTINE LOAD_RESERVOIR_STRUCTURE
!*********************************************************************
! Load the specified structure from the reservoir into coords
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine load_reservoir_structure(x,v,nr3,natom)

   implicit none
#  include "mpif.h"
#  include "parallel.h"

   _REAL_, dimension(*), intent(inout) :: x, v
   integer, intent(in) :: nr3, natom

   integer i, ierror
   character(len=80) line
   character(len=90) framename

!--------------------

   ! ----===== RREMD RESERVOIR LOADING =====----
   ! If we exchanged with the reservoir, swap the coordinates
   !  here (after the inpcrd have been read), Master process only.
   if (sanderrank==0) then
      ! Read restart file into coords. Need to set a filename for
      !  correct structure.
      write (6,'(a34)') "=========Reservoir Read==========="
      ! DAN ROE: integer size will need to correspond to maxreservoir
      ! Debug
      write (6,'(a,a)') "reservoirname=",trim(reservoirname)
      write (6,'(a,i6.6)') "rremd_idx=",rremd_idx
      ! DAN ROE: Should put a check in here to make sure framename
      !  doesn't get blown up
      write (framename,'(a,a1,i6.6)') &
         reservoirname(1:index(reservoirname," ")-1),".",rremd_idx
      write (6,'(a,a)') "Reservoir Filename= ",trim(framename)
      call amopen(39,framename,'O','F','R') 
      ! Read title, # atoms
      read (39,"(a)") line
      write (6,"(a7,a)") "Title= ", line
      ! DAN ROE: should this read have a format?
      read (39,*) i 
      if (i.ne.natom) then
         write (6,*) "restart file has wrong #atoms ",framename,i,natom
         backspace (39)
         read (39,"(a)") line
         write (6,"(a)") line
         call mexit(6,1)
      endif
      ! Read coordinates
      read (39,"(6(f12.7))") (x(i),i=1,nr3)
      ! DAN ROE: Debug
      write (6,'(a,i10)') "RREMD: coords read for frame ",rremd_idx
      write (6,*) (x(i),i=1,10)
      ! Read Velocities
      if (reserv_velo==1) then
         read (39,"(6(f12.7))") (v(i),i=1,nr3)
         write (6,'(a,i10)') "RREMD: velocities read for frame ",rremd_idx
         write (6,*) (v(i),i=1,10)
      endif
      close (unit=39)
      ! DAN ROE: Should be some error checking on reservoir read.
      ! CARLOS: IMPORTANT! Scale the reservoir velocities!
      !  Scaling is set in subrem
      ! If we have const P, need box change. Currently unsupported
      ! DAN ROE: checked for in mdread.
      !if (ifbox >= 1 .and. ntp > 0) then
      !   write (6,*) "const P not allowed for RREMD"
      !   call mexit (6,1)
      !endif
   endif ! sanderrank==0 

   ! Now master needs to broadcast coords and velo for reservoir
   !  to all processes in its group.
   call mpi_bcast(x,nr3,mpi_double_precision,0,commsander,ierror)
   if (reserv_velo==1) then
      call mpi_bcast(v,nr3,mpi_double_precision,0,commsander,ierror)
   endif
   ! DAN ROE: Debug
   if (sanderrank==0) then
      !write (6,*) "coordinates broadcast "
      !write (6,*) (x(lcrd+i-1),i=1,10)
      !if (reserv_velo==1) then
      !   write (6,*) "velocities broadcast "
      !   write (6,*) (x(lvel+i-1),i=1,10)
      !endif
      write (6,'(a38)') "==========End Reservoir Read=========="
   endif
   
   ! ----===== END RREMD RESERVOIR LOADING =====----

   return

end subroutine load_reservoir_structure


!*********************************************************************
!               SUBROUTINE HYBRID_REMD_ENE
!*********************************************************************
! Get energy of stripped structure in hybrid remd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine hybrid_remd_ene( &
              x,ix,ih,ipairs,qsetup,                        &
              numwatkeep,hybridgb,igb,nspm,t,temp0,                  &
              ntb,cut,                                               &
              ener,vir,do_list_update,nstep,nitp,nits,onefac,loutfm  ) 
   implicit none
#  include "parallel.h"
#  include "memory.h"

! sander.f
   _REAL_ x(*)
   integer ix(*), ipairs(*)
   character(len=4) ih(*)
   logical qsetup
! md.h
   integer numwatkeep, hybridgb, igb, nspm
   _REAL_ t,temp0
! memory.h
!   integer natom, nres
! box.h
   integer ntb
   _REAL_ cut
! runmd.f
   _REAL_ ener(*), vir(*), onefac(*)
   logical do_list_update,loutfm
   integer nstep, nitp, nits

! nstep, nitp, nits, and onefac only needed for printmd

! Temporary storage for natom, 
!  ntb, and cut during stripped coord call to force.
! DAN ROE: Is dynamic allocation the way to go?
! NOTE: This allocation currently won't work with amoeba.
   integer hybrid_natom, hybrid_ntb, ier, i
   _REAL_ hybrid_cut

!--------------------

! This is a hybrid REMD run. Get energy of stripped system for next
!  exchange.
! DAN ROE: Note: hybrid code is placed here since it needs access to force()
! 1- First the current coordinates and forces are placed in a temporary
!    array. The original ntb, cut, and natom values are saved. 
! 2- Then strip away all but numwatkeep waters from the temporary coords. 
!    natom is changed to reflect the new system size.
! 3- Set igb, ntb, and cut to hybrid values, and call force using the 
!    temporary coord (now stripped) and force arrays. The force call will use
!    GB based on the hybridgb variable. After we return from 
!    force we should have the correct PE for the stripped system. 
! DAN ROE: All threads need to do this because all need
!  to call force. All threads know about numwatkeep since it is a 
!  namelist variable.

   if (sanderrank==0) write(6,'(17("="),a,i10,17("="))') &
      "HYBRID REMD: energy calc for exch ",mdloop+1

! DAN ROE: Debug
!         if (master) then
!            write (6,*) "Pre-strip coordinates: "
!            write (6,*) (xx(lcrd+i-1),i=1,10)
!            write (6,*) "Pre-force velocities:"
!            write (6,*) (xx(lvel+i-1),i=1,10)
!            write (6,*) "Pre-force forces:"
!            write (6,*) (xx(lforce+i-1),i=1,10)
!            call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)
!         endif

   ! 1- Store coords and forces in temp arrays,
   !    backup natom, ntb, and cut
   do i=1, natom
      hybrid_coord(3*i - 2) = x(lcrd+3*i - 3)
      hybrid_coord(3*i - 1) = x(lcrd+3*i - 2)
      hybrid_coord(3*i)     = x(lcrd+3*i - 1)
      hybrid_force(3*i - 2) = x(lforce+3*i - 3)
      hybrid_force(3*i - 1) = x(lforce+3*i - 2)
      hybrid_force(3*i)     = x(lforce+3*i - 1)
   enddo
   hybrid_natom=natom
   hybrid_ntb=ntb
   hybrid_cut=cut

   ! 2- Strip waters. This will change the coords (hybrid_coord) and
   !    # of atoms (natom). The coords will be imaged.
   if (sanderrank==0) write(6,'(a)') "HYBRID REMD: Stripping waters"
   call stripwat(natom,hybrid_coord,ix(i02),ih(m02),x(lmass), &
                 ix(i70),nspm,nres,numwatkeep)
   if (sanderrank==0) then
      write(6,'(a,i8)') "HYBRID REMD: New natom= ",natom
      !write (6,*) "Post-strip coordinates: "
      !write (6,*) (xx(lcrd+i-1),i=1,10)
      ! DAN ROE: write to the stripped coordinate file
      !write (remstripcoord_unit,'(a,3(1x,i8),1x,f8.3)') &
      !                          "REMD ", repnum, mdloop, &
      !                          remstep, mytargettemp
      if (hybridwritetraj) &
         call corpac(hybrid_coord,1,natom*3,remstripcoord_unit,loutfm)
   endif
   call mpi_barrier(commworld, ier)

   ! 3- Call force to calculate energy using the stripped
   !    coordinates and GB (based on hybridgb).
   ! Make sure PBC off, and change the cutoff
   igb=hybridgb
   ntb=0
   cut=9801.0d0 !  = 99 * 99
   if (sanderrank==0) write(6,'(a)') "HYBRID REMD: Calling force."
   call force(x,ix,ih,ipairs,hybrid_coord,hybrid_force,ener(23),vir, &
              x(l96),x(l97),x(l98),x(l99),qsetup, &
              do_list_update)
! DAN ROE: Debug
   if (sanderrank==0) then
!            write (6,*) "Post-force coordinates: "
!            write (6,*) (xx(lcrd+i-1),i=1,10)
!            write (6,*) "Post-force velocities:"
!            write (6,*) (xx(lvel+i-1),i=1,10)
!            write (6,*) "Post-force forces:"
!            write (6,*) (xx(lforce+i-1),i=1,10)
      call prntmd(nstep,nitp,nits,t,ener,onefac,7,.false.)
   endif
   if (sanderrank==0) write(6,'(a,f12.6,a,f6.2)') &
      "HYBRID REMD: myEptot= ",ener(23)," myTargetTemp= ",temp0

   ! 4- Restore original natom, igb, cut
   if (sanderrank==0) write(6,'(a)') "HYBRID REMD: Restoring..."
   natom=hybrid_natom
   igb=0
   ntb=hybrid_ntb
   cut=hybrid_cut

! DAN ROE: Debug
!         if (master) then
!            write (6,*) "Restored coordinates: "
!            write (6,*) (xx(lcrd+i-1),i=1,10)
!            write (6,*) "Restored velocities:"
!            write (6,*) (xx(lvel+i-1),i=1,10)
!            write (6,*) "Restored forces:"
!            write (6,*) (xx(lforce+i-1),i=1,10)
!         endif
   if (sanderrank==0) write(6,'(25("="),a,25("="))') &
      "END HYBRID REMD energy calc."
   
   return

end subroutine hybrid_remd_ene


!*********************************************************************
!               SUBROUTINE STRIPWAT
!*********************************************************************
! strip water from structure for replica exchange with mixed solvent models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine stripwat(nrp,x,ipres,lbres,amass,nsp,nspm,nres,numwatkeep)

   implicit none
#  include "parallel.h"
! box.h needed for imaging water
#  include "box.h"

   integer, parameter :: size=10000

   integer nrp, ipres(*), nsp(*), nspm, nres, numwatkeep
   _REAL_  x(*), amass(*)
   character(len=4) lbres(*)

   _REAL_  x2(size), closedist(size), xcm(3)
   _REAL_  r2, xij, yij, zij, xi, yi, zi, tempr, aamass, tmassinv  

   integer watresnum(size) 
   integer firstwat, watpointer, newnrp, totwat, tempi 
   integer i, i1, j, j1, j2, k

   logical keepwat(size), templ

!----------------------------------------------------------
   ! find out where water starts

      do i=1,nres

         if (lbres(i).eq."WAT") then
            firstwat=i
            go to 4120
         endif
      enddo
4120  continue

      !if (sanderrank==0) write (6,*) "first water residue is # ",firstwat
      ! call amflsh(6)

      ! next we need to reimage the system so that waters surround the solute
      ! we could use minimum image type calculation to get the
      ! closest distance but we still need the water imaged properly
      ! for the GB calculation (which is not periodic)
      ! follow the code for imaging from runmd's iwrap=1 code
      ! first center the system on the CM of the solute

      xcm(1) = 0.d0
      xcm(2) = 0.d0
      xcm(3) = 0.d0

      ! here tmassinv is only for non-water

      tmassinv=0.d0
      i = 0
      do k=1,firstwat-1
         do j =ipres(k),ipres(k+1)-1
            aamass = amass(j)
            xcm(1) = xcm(1) + x(i+1)*aamass
            xcm(2) = xcm(2) + x(i+2)*aamass
            xcm(3) = xcm(3) + x(i+3)*aamass
            i = i + 3
            tmassinv=tmassinv+aamass
         enddo
      end do

      tmassinv=1.d0/tmassinv

      xcm(1) = xcm(1) * tmassinv
      xcm(2) = xcm(2) * tmassinv
      xcm(3) = xcm(3) * tmassinv

      ! center all atoms, not just solute

      do i=1,nrp
         x(3*i-2) = x(3*i-2)-xcm(1)
         x(3*i-1) = x(3*i-1)-xcm(2)
         x(3*i)   = x(3*i)-xcm(3)
      enddo

      ! now re-image the box

      call wrap_molecules(nspm,nsp,x)
      if(ifbox == 2) call wrap_to(nspm,nsp,x,box)

      ! now start setting closest distance between water and solute

      ! in fact we do not need the distance, the distance squared is
      ! fine since we only want to sort them

      !if (sanderrank==0) &
      !   write (6,*) "looping over the waters to check distance"

      watpointer=0

      do i= firstwat,nres

         ! CHECK TO MAKE SURE IT IS WATER

         if (lbres(i).ne."WAT ") then
            if (sanderrank==0) then
               write (6,*) "solvent molecule is not water: ",lbres(i)
               write (6,*) "stopping water search"
            endif
            call mexit(6,1)
         endif
         
         watpointer=watpointer+1

         ! closedist(i) is the distance from the water to the 
         ! closest solute atom

         closedist(watpointer)=9999
         watresnum(watpointer)=i
           
         !  do i1=ipres(i),ipres(i+1)-1

         ! for water, just take distance to oxygen (first atom in water)

         i1=ipres(i)
         xi = x(3*i1-2)
         yi = x(3*i1-1)
         zi = x(3*i1)

         ! loop over non-water residues/atoms
         
         do j=1,firstwat-1
            do j1=ipres(j),ipres(j+1)-1

               xij = xi - x(3*j1-2)
               yij = yi - x(3*j1-1)
               zij = zi - x(3*j1  )
               r2 = xij*xij + yij*yij + zij*zij
         
               if (r2.lt.closedist(watpointer)) closedist(watpointer)=r2
            enddo
         enddo
      enddo
         
      totwat=watpointer

      ! now we have a list of all waters - their residue number and distance
      ! to solute

      ! sort them by distance

      do i=1,totwat-1
         do j=i+1,totwat

           ! write (6,*) "working on ",i,j

            if (closedist(i).gt.closedist(j)) then
               tempr=closedist(i)
               tempi=watresnum(i)
               closedist(i)=closedist(j)
               watresnum(i)=watresnum(j)
               closedist(j)=tempr
               watresnum(j)=tempi

            endif
         enddo
      enddo

      ! now set save flags for closest numwatkeep

      !if (sanderrank==0) write (6,*) "numwatkeep is ",numwatkeep," of ", totwat

      do i=1,numwatkeep
         keepwat(i)=.true.
      enddo
      do i=numwatkeep+1,totwat
         keepwat(i)=.false.
      enddo

      ! now sort them back into the order by watresnum

      do 4140 i=firstwat,nres

         ! i is the residue number we are looking for

         ! i1 is the current water at this residue number
         ! this is in the 1 to totwat sense so we can pull those
         ! indices for swapping waters in list

         i1=i-firstwat+1

         ! look to see where this water is, and put it in place

         do j=1,totwat

            if (watresnum(j).eq.i) then

               ! found it, so swap them

               tempr=closedist(i1)
               tempi=watresnum(i1)
               templ=keepwat(i1)
               closedist(i1)=closedist(j)
               watresnum(i1)=watresnum(j)
               keepwat(i1)=keepwat(j)
               closedist(j)=tempr
               watresnum(j)=tempi
               keepwat(j)=templ

               ! get next i

               go to 4140 
            endif
         enddo
4140  enddo 
             
      ! now go through and write the restart file. we need to write to temp
      ! array since we can't remove atoms while writing the file itself since
      ! there is more than 1 molecule per line

      newnrp=0
      do i=1,nres
         if (i.lt.firstwat.or.keepwat(i-firstwat+1)) then
            do j=ipres(i),ipres(i+1)-1
               newnrp=newnrp+1
               x2(3*newnrp-2) = x(3*j-2)
               x2(3*newnrp-1) = x(3*j-1)
               x2(3*newnrp)   = x(3*j)
            enddo
         endif
      enddo

      ! now copy this array to old one, reset nrp

      nrp=newnrp
      do i=1,nrp
        x(3*i-2) = x2(3*i-2)
        x(3*i-1) = x2(3*i-1)
        x(3*i)   = x2(3*i)
      enddo
       
   return

end subroutine stripwat


!*********************************************************************
!               SUBROUTINE CALC_RREMD_CLUSTER 
!*********************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [calculate cluster of current coordinates based on dihedral binning, 
!   used with rremd==3. Dihedral calc. based on Jcoupling routine from nmr.f]
! Called from runmd(). Sets incluster. Only masters call
! this subroutine.
subroutine calc_rremd_cluster(x)
   implicit none

! Input:
!    X(I)  : Coordinate array.
   _REAL_, dimension(*), intent(in) :: x

! Internal constants
   _REAL_, parameter :: small=1.0d-14, small2=1.0d-5, pi=3.14159265358979323846d0
! Internal variables
   _REAL_ xij(3),xkj(3),xkl(3),t(6),dc(6)
   _REAL_ ajcoef(3)
   _REAL_ rij2, rkj2, rkl2, value 
   _REAL_ dx, dy, dz, gx, gy, gz
   _REAL_ bi, bk, ct, z1, z2, s
   integer m, i, j
!  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
!        I4  atom number.
   integer i1, i2, i3, i4

!--------------------

   do i=1,nclustdih
      i1=3*(dihclustat(1,i)-1)
      i2=3*(dihclustat(2,i)-1)
      i3=3*(dihclustat(3,i)-1)
      i4=3*(dihclustat(4,i)-1)

      ! ----- Dihedral Calculation -----
      !write (6,*) "calculating dih for ",i1,i2,i3,i4

      ! Calculate the underlying torsion:

      rij2 = 0.0d0
      rkj2 = 0.0d0
      rkl2 = 0.0d0

      do m=1,3
         xij(m) = x(i1+m) - x(i2+m)
         xkj(m) = x(i3+m) - x(i2+m)
         xkl(m) = x(i3+m) - x(i4+m)
         rij2 = rij2 + xij(m)**2
         rkj2 = rkj2 + xkj(m)**2
         rkl2 = rkl2 + xkl(m)**2
      enddo

      ! Calculate ij X jk AND kl X jk

      dx = xij(2)*xkj(3) - xij(3)*xkj(2)
      dy = xij(3)*xkj(1) - xij(1)*xkj(3)
      dz = xij(1)*xkj(2) - xij(2)*xkj(1)

      gx = xkj(3)*xkl(2) - xkj(2)*xkl(3)
      gy = xkj(1)*xkl(3) - xkj(3)*xkl(1)
      gz = xkj(2)*xkl(1) - xkj(1)*xkl(2)

      ! Calculate the magnitudes of above vectors, and their dot product:

      bi = dx*dx + dy*dy + dz*dz
      bk = gx*gx + gy*gy + gz*gz
      ct = dx*gx + dy*gy + dz*gz

      ! If this is a linear dihedral, we cannot calculate a value for it,
      ! so set value to -999 and return

      if (bk < 0.01d0 .or. bi < 0.01d0) then
         value=-999.
         return
      end if

      bi = sqrt(bi)
      bk = sqrt(bk)
      z1 = 1.0d0/bi
      z2 = 1.0d0/bk
      ct = ct*z1*z2

      !write (6,*) "ct is ",ct

      if (ct > 1.0d0-small) ct = 1.0d0-small
      if (ct < -1.0d0+small) ct = -1.0d0+small
      value = acos(ct)

      s = xkj(1)*(dz*gy-dy*gz) + xkj(2)*(dx*gz-dz*gx) + &
          xkj(3)*(dy*gx-dx*gy)

      if (s < 0.0d0) value = -value
      value = pi - value

      ! now convert to degrees!

      value = value * 180.d0/pi

      ! ----- End Dihedral Calculation -----
      
      ! dihedral is in 0 to 360 range. Convert to -180 to 180
      if (value>180.d0) value=value-360

      ! Shift all values up to 0 to facilitate binning
      currdihid(i)=(value+180.d0) / (360.d0 / float(dihclustnbin(i)))
      !write (6,"(a4,i5,a7,3(1x,f8.3))") &
      !       "dih:",i," atom1:",x(i1+1),x(i1+2),x(i1+3)
      !write (6,"(a4,i5,a7,3(1x,f8.3))") &
      !       "dih:",i," atom2:",x(i2+1),x(i2+2),x(i2+3)
      !write (6,"(a4,i5,a7,3(1x,f8.3))") &
      !       "dih:",i," atom3:",x(i3+1),x(i3+2),x(i3+3)
      !write (6,"(a4,i5,a7,3(1x,f8.3))") &
      !       "dih:",i," atom4:",x(i4+1),x(i4+2),x(i4+3)
      !write (6,*) "dih:",i," value:",value," bin:",currdihid(i)
   enddo ! 1, nclustdih

   ! Set the default cluster to #0, which means "no match". clustersize(0) was
   !  set to 0 in multisander.f. This means that "no match" structures will 
   !  automatically be rejected. 
   ! Eventually the reservoir could have "no match" structures added to it.
      incluster=0
      do i=1,nclust
         do j=1,nclustdih
         ! As soon as a bin doesn't match, exit
            if (clusterid(i,j).ne.currdihid(j)) go to 4600
         enddo
         incluster=i
         go to 4610 ! Exit as soon as a cluster matches
4600  enddo
4610  continue
      write (6,'(a,i10)') &
         "RREMD: Current MD structure was assigned to cluster ",incluster

   return

end subroutine calc_rremd_cluster 


!*********************************************************************
!               SUBROUTINE LOAD_RESERVOIR_FILES 
!*********************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [load energy/cluster information for doing reservoir REMD]
! All threads call this subroutine
subroutine load_reservoir_files()

   implicit none

#  include "parallel.h"
#  include "mpif.h"

   integer numatoms, iseed, ierror, i, nodeid
! i1-4 for reading dihedral atom nums
   integer j, i1, i2, i3, i4

! --------------------
   if (sanderrank==0) then
      ! Open saveene, which contains energies for each reservoir struct
      !write (41+worldrank,*) "loading ensemble energies"
      call amopen(32,saveenefile,'O','F','R')

      ! Read saveene file header line.
      ! Format:
      ! <# reservoir structures> <reservoir T> <#atoms> <random seed>
      !   <velocity flag [1 if reading velocity from reservoir]>
      ! Read # atoms - we dont seem to know natom yet.
      ! DAN ROE: Check if we need #atoms read here.
      ! Note: For cleanliness we could make the unit here the remlog
      !  unit instead of the more arbitrary 32 since remlog is not used
      !  yet.
      read (32,*,err=5000) reservoirsize,restemp0,numatoms,iseed,reserv_velo

      ! Write the reservoir type and other info to remtype 
      if (worldrank==0) then
         write (remtype_unit,'(a,a)') &
            "RREMD: Info from saveene file ",saveenefile
         write (remtype_unit,'(a,i5)') "  NumAtoms= ",numatoms
         write (remtype_unit,'(a,i5)') "  ReservoirSize= ",reservoirsize
         write (remtype_unit,'(a,f6.2)') "  ReservoirTemp(K)= ",restemp0
         write (remtype_unit,'(a,i10)') "  ReservoirRandomSeed= ",iseed
         if (reserv_velo.eq.1) then
            write (remtype_unit,'(a)') &
               "  Velocities will be read from&
               & reservoir restart files"
         else
            write (remtype_unit,'(a)') &
               "  Velocities will be assigned to&
               & structure after exchange"
         endif
        ! Print reservoir type information
         write (remtype_unit,'(a,i5)') "RREMD: Reservoir type ",rremd
         if (rremd.eq.1) then
            write (remtype_unit,'(a)') &
               "  Boltzmann weighted reservoir,&
               & exchange uses delta beta"
         elseif (rremd.eq.2) then
            write (remtype_unit,'(a)') &
               "  Non-Boltzmann 1/N weighted reservoir,&
               & exchange uses beta"
         elseif (rremd.eq.3) then
            write (remtype_unit,'(a)') &
               "  Non-Boltzmann weighted reservoir with&
               & defined weights"
            write (remtype_unit,'(a)') &
               "  (Currently only works via dihedral clustering)"
            write (remtype_unit,'(a)') &
               "  Exchange uses beta and weights."
            write (remtype_unit,'(a)') &
               "  FORCING WEIGHT OF 0 FOR STRUCTURES&
               & NOT IN RESERVOIR!"
         else
            write (remtype_unit,'(a,i5)') &
               "Unknown reservoir type: rremd=", rremd
            !call mexit(6,1)
         endif
      endif ! worldrank==0 remtype file write

      ! Exit if unknown reservoir type
      ! DAN ROE: Should this exit just be called above?
      if (rremd.lt.0 .or. rremd.gt.3) call mexit(6,1)

      ! Check reservoir size limits
      ! Currently the format for reservoir files allows a max integer
      ! width of 6. To accomodate a larger reservoir the format string
      ! in load_reservoir_structure would have to be changed.
      if (reservoirsize.gt.999999) then
         write (6,*) "RREMD ERROR: Reservoir size limit is currently &
                     &999,999. To accomodate larger reservoir sizes &
                     &edit the reservoir file format string in &
                     &load_reservoir_structure (remd.f)"
         call mexit(6,1)
      endif

      ! Allocate memory
      allocate(saveene(reservoirsize), stat=ierror)
      REQUIRE(ierror==0)

      ! Read energy (and cluster number for rremd==3) for each structure
      ! All we really need to load are the energies; if the exchange
      !  is successful we can grab the coords (and velocity if necessary
      !  from the disk during the run.
      ! DAN ROE: For rremd<3, If cluster #s are present this will break.
      !  Should add a check. Maybe make cluster #s a separate file?
      do i=1,reservoirsize
         if (rremd==1.or.rremd==2) then
            ! 1 = Boltzmann weighted, read energies
            ! 2 = 1/N reservoir, read energies
            read(32,*,err=5000,iostat=ierror) saveene(i)
         elseif (rremd==3) then 
            ! Weighted reservoir, read energies and corresponding cluster #s
            read(32,*,err=5000,iostat=ierror) saveene(i), clusternum(i)
            ! DAN ROE: Since I started numbering clusters at 0, add 1
            !  Change later to be consistent.
            !clusternum(i)=clusternum(i)+1
            if (clusternum(i)<1) then
               write(6,*) &
                  "Error: Cluster# < 1 not allowed."
               goto 5000
            endif
         endif
         if (ierror<0.and.i<reservoirsize) then
            write(6,*) &
               "Error: EOF reached before all values read."
            goto 5000
         endif
         ! DAN ROE: Debug
         if (worldrank==0) then
            write (remtype_unit,*) "frame,energy ", i,saveene(i)
            if (rremd==3) write(remtype_unit,*) "  cluster ",clusternum(i)
         endif
      enddo
      go to 5010

5000  write (6,*) "RREMD: Error in reading saveene!"
      call mexit(6,1)

5010  continue
      close (unit=32)

      ! Read in the cluster info file for rremd==3
      ! DAN ROE: Put in error checking on reads
      if (rremd==3) then
         call amopen(32,clusterinfofile,'O','F','R')
         ! First read number of dihedrals to bin
         read (32,*,err=5020) nclustdih
         if (nclustdih > maxdihclust) then
            !write (6,*) "NumDihedrals > maxdihclust in remd.f: ",maxdihclust
            if (worldrank==0) then
               write (6,*) &
                  "NumDihedrals > maxdihclust in remd.f: ",maxdihclust
            endif
            call mexit(6,1)
         endif
         ! Now read atom #s and bins for each dihedral angle.
         ! Atom #s should start from 1
         do i=1,nclustdih
            ! Format: atom#1 atom#2 atom#3 atom#4 bins
            read (32,"(10(i10,1x))",err=5020) i1,i2,i3,i4,j
            dihclustat(1,i)=i1
            dihclustat(2,i)=i2
            dihclustat(3,i)=i3
            dihclustat(4,i)=i4
            dihclustnbin(i)=j
         enddo
         ! Read number of clusters
         read (32,*,err=5020) nclust
         if (nclust > maxclust) then
           !write (6,*) "NumClusters > maxclust in remd.f: ",maxclust
            if (worldrank==0) then
               write (6,*) &
                  "NumClusters > maxclust in remd.f: ",maxclust
            endif
            call mexit(6,1)
         endif
         ! Read cluster weight and ID (bin values) for each cluster
         ! Set the cluster size for cluster 0, which corresponds to an MD 
         !  structure not being present in the reservoir. The current code
         !  does not add this to a new cluster.
         ! This makes the reservoir exchange rigorous - structures not in
         !  the reservoir have a weight of zero, and therefore can't be
         !  exchanged.
         clustersize(0)=0
         do i=1,nclust
            ! Format: Cluster# Weight Bin1...BinNclustdih
            ! DAN ROE: This assumes clusters are in order!
            read (32,"(2(i10,1x),16(i3))",err=5020) &
               j, clustersize(i),(clusterid(i,j),j=1,nclustdih)
         enddo
         goto 5030

5020     write (6,*) "RREMD: Error in reading clusterinfo!"
         call mexit(6,1)

5030     continue
         close (unit=32)

         ! Overall master Write out information to remtype
         if (worldrank==0) then
            write(remtype_unit,'(a,a)') &
               "RREMD: clusterinfo file ",clusterinfofile
            write (remtype_unit,*) "  NumDihedrals= ", nclustdih
            do i=1,nclustdih
               write (remtype_unit,"(a14,i5,a7,4(1x,i5),a7,i5)") &
                  "    Dihedral #",i," atoms:",dihclustat(1,i), &
                  dihclustat(2,i),dihclustat(3,i),dihclustat(4,i), &
                  " bins: ",dihclustnbin(i)
            enddo
            write (remtype_unit,*) "  NumClusters= ",nclust
            do i=1,nclust
               write (remtype_unit,*) &
                  "    Cluster ",i," has ",clustersize(i)," members"
               write (remtype_unit,'(a12,16(i3))') &
                  "      Bins= ",(clusterid(i,j),j=1,nclustdih)
            enddo
         endif
      endif ! rremd == 3 
   endif ! sanderrank == 0

   ! Set random # generator based on the seed in saveene; this is so
   !  we can do independent runs with different sequence of structures
   !  from the reservoir. Here we have no concept of the ig seed in 
   !  md.in and we don't want to hardcode, so saveene is a convenient
   !  location. All threads do this.
   ! DAN ROE: FileDebug
   ! Send iseed to all threads since everyone needs to call amrset
   ! NOTE: repnum = nodeid + 1
   nodeid = repnum - 1
   call mpi_bcast(iseed,1,mpi_integer,0,commworld,ierror)
   !   write (50+worldrank,'(a,i4,a,i10)') &
   !      " amrset for Res. on", worldrank,&
   !      " with value ",iseed+17*nodeid
   call amrset(iseed + 17 * nodeid)
   ! Broadcast reserv_velo to all so we know it in sander() when we 
   !  read reservoir.
   ! DAN ROE: Does every process need this or just masters?
   call mpi_bcast(reserv_velo,1,mpi_integer,0,commworld,ierror)

   return

end subroutine load_reservoir_files

#endif /* MPI */

end module remd

