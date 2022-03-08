! <compile=optimized>
#include "dprec.h"
#include "assert.h"

!+ MMTSB Replica Exchange

!     Module MMTSB_Replica_Exchange

! Description:
! Module for replica exchange using Multiscale Modeling Tools in
! Structural Biology (MMTSB), mmtsb.scripps.edu.
! Currently, both temperature and lambda replica exchange are available.
! These calculations are enabled via compilation with the preprocessor
! name MMTSB defined.
! Documentation does not exist yet.
! The public mmtsb routines include:
!   Subroutine Mmtsb_init( )
!          Initialize a replica exchange calculation.
!   Subroutine Mmtsb_newtemp( pot_energy, temp, isdone )
!          Update temperature in replica exchange.
!   Subroutine Mmtsb_print_banner( )
!          Print details about this replica exchange calculation.

! History:
! $Id: mmtsb.f,v 9.0 2006/04/03 23:35:55 case Exp $

! Code Description:
!   Language:           Fortran 77.
!   Software Standards: Internal Amber Standards.


! **********************************************************************
!  Copyright 2003                                                      *
!                                                                      *
!   Modified BSD license                                               *
!                                                                      *
!   Redistribution and use in source and binary forms, with or without *
!   modification, are permitted provided that the following conditions *
!   are met:                                                           *
!                                                                      *
!    1.Redistributions of source code must retain the above copyright  *
!      notice, this list of conditions and the following disclaimer.   *
!    2.Redistributions in binary form must reproduce the above         *
!      copyright notice, this list of conditions and the following     *
!      disclaimer in the documentation and/or other materials provided *
!      with the distribution.                                          *
!    3.The name of the author may not be used to endorse or promote    *
!      products derived from this software without specific prior      *
!      written permission.                                             *
!                                                                      *
!   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS''                  *
!   AND ANY EXPRESS OR IMPLIED WARRANTIES,                             *
!   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                         *
!   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR                      *
!   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT                   *
!   SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,                         *
!   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                       *
!   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT                          *
!   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR                     *
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                        *
!   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON                       *
!   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,                      *
!   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE                    *
!   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE                    *
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE                           *
!   POSSIBILITY OF SUCH DAMAGE.                                        *
!                                                                      *
!  To report bugs, suggest enhancements, etc., contact:                *
!    MMTSB, mmtsb@scripps.edu                                          *
!                                                                      *
! **********************************************************************

!     implicit none

!     Private
!     Public :: Mmtsb_init
!     Public :: Mmtsb_newtemp
!     Public :: Mmtsb_print_banner


!     Contains
! public routines in alphabetical order

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mmtsb_init here]
subroutine mmtsb_init( initial_temp, initial_lambda )
   !
   ! Initialize a replica exchange calculation obtaining the initial
   ! value(s) of the independent variable(s) from the server.
   ! For temperature replica exchange, if the input temp, which is the
   ! reference temp from the amber input file, is different than the
   ! temp obtained from the server then an exchange is considered
   ! to have occurred (consequently, the velocities will be reset).
   ! This will be called from mdread.f.
   !
   !------------------------------------------------------------------
   use constants, only : TEN_TO_MINUS3, zero
   implicit none
#  include "files.h"
#  include "mmtsb.h"
#  include "md.h"

   _REAL_      initial_lambda  ! reference lambda., intent(inout)
   _REAL_      initial_temp    ! reference temp., intent(inout)

   integer     chainlength
   logical     is_done_mmtsb
   _REAL_      lambda_mmtsb
   _REAL_      temp_mmtsb

   if( mmtsb_switch /= mmtsb_off ) then
#ifdef MMTSB
      mmtsb_is_exchanged = .false.
      write(6,'(5x,a,a)') 'Reading MMTSB setup from file: ' , &
            mmtsb_setup_file(1:43)
      call amopen( mmtsb_unit, mmtsb_setup_file, 'O','F','R' )
      read( mmtsb_unit, * )     chainlength
      read( mmtsb_unit, * )     sendfiles
      read( mmtsb_unit, '(A)' ) jobid
      read( mmtsb_unit, '(A)' ) servername
      read( mmtsb_unit, '(A)' ) serverport
      read( mmtsb_unit, '(A)' ) serverid
      read( mmtsb_unit, '(A)' ) datadir
      close( mmtsb_unit)

      !write(6,'(5x,a,i8)')   'chainlength  = ', chainlength ! unused by Amber
      write(6,'(5x,a,a60)')  'sendfiles    = ', sendfiles
      write(6,'(5x,a,a60)')  'jobid        = ', jobid
      write(6,'(a,i8)') '|    servername   = ', servername
      write(6,'(5x,a,a60)')  'serverport   = ', serverport
      write(6,'(a,i8)') '|    serverid     = ', serverid
      write(6,'(5x,a,a60)')  'datadir      = ', datadir

      ! Does the energy need to be computed for the first
      ! call to the server - guess not.
      if ( mmtsb_switch == mmtsb_temp_rex ) then
         temp_mmtsb = initial_temp  ! this may not be necessary.
         call mmtsb_newtemp( zero, temp_mmtsb, is_done_mmtsb )
      else if ( mmtsb_switch == mmtsb_lambda_rex ) then
         ! the server does not need the initial lambda
         call mmtsb_newlambda( zero, zero, lambda_mmtsb, temp_mmtsb, &
               is_done_mmtsb )
      else
         write(6,'(a)') &
               'MMTSB Error: Unknown mmtsb_switch value !'
         call mexit( 6, 1 )
      end if
      if ( is_done_mmtsb ) then
         write(6,'(a)') &
               'MMTSB Error: Finished upon initialization call to server !'
         call mexit( 6, 1 )
      end if

      if ( mmtsb_switch == mmtsb_temp_rex ) then
         if ( abs( temp_mmtsb - initial_temp ) <= TEN_TO_MINUS3 ) then
            ! no exchange
            mmtsb_is_exchanged = .false.
         else
            ! temp exchange occurred
            ! the velocities will be randomly reset at the new temp
            mmtsb_is_exchanged = .true.
            initial_temp = temp_mmtsb
         end if
         write(6,'(5x,a,f10.5)') 'temp0    = ', initial_temp
      else if ( mmtsb_switch == mmtsb_lambda_rex ) then
         if( icfe == 0 ) then
            write(6,'(a)') &
                  'MMTSB Warning: Invalid icfe for lambda replica exchange!'
         end if
         initial_temp = temp_mmtsb  ! server initial temp has precedence
         if ( abs( lambda_mmtsb - initial_lambda )<=TEN_TO_MINUS3) then
            ! no exchange
            mmtsb_is_exchanged = .false.
         else
            ! lambda exchange occurred
            ! the velocities will be randomly reset
            mmtsb_is_exchanged = .true.
            initial_lambda = lambda_mmtsb
         end if
         write(6,'(5x,a,f10.5)') 'temp0    = ', initial_temp
         write(6,'(5x,a,f10.5)') 'clambda  = ', initial_lambda
      end if

#else
      write(6,'(a)') &
            'MMTSB not enabled but requested with nonzero mmtsb_switch !'
      write(6,'(5x,a)') 'Recompile with MMTSB defined, -DMMTSB.'
      call mexit( 6, 1 )
#endif
   end if  ! ( mmtsb_switch /= mmtsb_off )

end subroutine mmtsb_init


#ifdef MMTSB
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mmtsb_newlambda here]
subroutine mmtsb_newlambda( unpert_pe, pert_pe,lambda,temp,isdone)
   !
   ! Update lambda in replica exchange.
   ! Eventually multidimensional replica exchange, lambda and temp,
   ! will be supported, but for now only the initial temp may be
   ! obtained.
   ! This will be called from runmd.f
   ! Developer note: perhaps lambda should be intent(inout) and
   ! mmtsb_is_exchanged assigned here.
   !
   !------------------------------------------------------------------
   use constants, only : TEN_TO_PLUS3, one
   implicit none
#  include "mmtsb.h"

   _REAL_  unpert_pe  ! lambda=0 Potential En. of this cycle, intent(in)
   _REAL_  pert_pe    ! lambda=1 Potential En. of this cycle, intent(in)
   _REAL_  lambda     ! New lambda,                      intent(out)
   _REAL_  temp       ! initial temp,                    intent(out)
   logical isdone     ! replica exchange calculation finished,intent(out)

   external getbias
   _REAL_   newtemp
   external newtemp

   _REAL_   krg      ! unused in our implementation
   _REAL_   krho     ! unused in our implementation
   _REAL_   lastrho  ! unused in our implementation
   _REAL_   rho      ! unused in our implementation

   ! variables are required as arguments to these C functions
   ! a negative one value indicates inactivity
   krg     = -one
   krho    = -one
   rho     = -one
   ! this value indicates that the 7th argument is a potential energy
   lastrho = -TEN_TO_PLUS3
   temp = newtemp( servername, serverport, serverid, jobid, &
         datadir, unpert_pe, pert_pe, lastrho, sendfiles )
   isdone = temp < -9999.0
   call getbias( lambda, krg, rho, krho )

end subroutine mmtsb_newlambda


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mmtsb_newtemp here]
subroutine mmtsb_newtemp( pot_energy, temp, isdone )
   !
   ! Update temperature in replica exchange.
   ! This will be called from runmd.f
   ! Developer note: perhaps temp should be intent(inout) and
   ! mmtsb_is_exchanged assigned here.
   !
   !------------------------------------------------------------------
   use constants, only : one
   implicit none
#  include "mmtsb.h"

   _REAL_  pot_energy ! Final Potential Energy of this cycle, intent(in)
   _REAL_  temp       ! New temperature,                      intent(out)
   logical isdone     ! replica exchange calculation finished,intent(out)

   _REAL_   newtemp
   external newtemp
   _REAL_   lastrg   ! unused in temperature replica exchange
   _REAL_   lastrho  ! unused in temperature replica exchange

   ! variables are required as arguments to this C function
   ! a negative value indicates inactivity
   lastrg  = -one
   lastrho = -one
   temp = newtemp( servername, serverport, serverid, jobid, &
         datadir, pot_energy, lastrg, lastrho, sendfiles )
   isdone = temp < -9999.0

end subroutine mmtsb_newtemp
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mmtsb_print_banner here]
subroutine mmtsb_print_banner( )
   !
   ! Print details about this replica exchange calculation.
   ! This will be called from mdread.f.
   !
   !------------------------------------------------------------------
   implicit none
#  include "mmtsb.h"

   if( mmtsb_switch /= mmtsb_off ) then
      write(6,'(/a)') 'MMTSB Replica Exchange:'
      write(6,'(5x,a,i8)') 'mmtsb_iterations  = ', mmtsb_iterations
   end if

end subroutine mmtsb_print_banner


! private routines in alphabetical order


!     End module MMTSB_Replica_Exchange

