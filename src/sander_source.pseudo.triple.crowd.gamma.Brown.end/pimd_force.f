#include "dprec.h"
#include "assert.h"

!------------------------------------------------------------------------------
subroutine pimd_part_spring_force( x, f, mass, Epot_spring, Epot_deriv, dvdl )
! computes the harmonic-spring energy and forces arising from discretized
! path integral action.
! computes the kinetic term of the virial estimator for the internal energy
! Computes the estimators for the thermonynamic integration w.r.t. mass.
!------------------------------------------------------------------------------
  use constants, only: hbar, kB
  use pimd_vars, only: Eimp_virial, x_centroid, nbead, dmdlm, itimass
#ifdef MPI
#ifdef LES
  use evb_pimd,  only: natomCL, vel0_bead, bead_dcrypt
  use miller,    only: i_qi, gradRC, div_ndx
  use evb_parm, only: nbias
  use evb_data,  only: f_v, F_QI, G_QI
#endif 
#endif

  implicit none

#include "md.h"
#include "les.h"
#include "memory.h"

  _REAL_, intent(in) :: x( 3, natom ) !! coordinates
  _REAL_, intent(in) :: mass(natom)
  _REAL_ :: f( 3, natom )             !! forces
  _REAL_ :: Epot_spring, Epot_deriv, beta, kT, dvdl
  integer :: idim, iatom, icopy, istart, iend, ierr
  _REAL_ :: coeff, tmp, tmp_x, tmp_f, tmp_e, tmp_ti, sum_e, sum_d, sum_v

  _REAL_ :: coeff_F, coeff_G1, coeff_G2
  _REAL_ :: F_sum1, F_sum2, F_sum3, F_sum4, G_sum, Ftot, Gtot
  _REAL_ :: ftmp, fv

  integer :: m, n, mm, i3

  _REAL_ , intrinsic :: sqrt, dot_product

#ifdef MPI
# include  "mpif.h"
  _REAL_ :: ener(4), temp_ener(4) ! JVAN
#ifdef LES
  _REAL_ :: dx(3*natomCL)
#endif
#endif

#include "parallel.h"

#ifdef MPI
  istart = iparpt(mytaskid) + 1
  iend   = iparpt(mytaskid+1) 
#else
  istart = 1
  iend   = natom
#endif

  kT = kB * temp0
  beta = 1.0d0 / kT

  coeff = dble( ncopy ) / ( hbar * beta )**2 
  sum_e = 0.d0
  sum_d = 0.d0
  sum_v = 0.d0
  dvdl = 0.d0 ! JVAN

  !For QI.
  coeff_F = 2.d0 / dble( ncopy )
  coeff_G1 = 6.d0 * dble( ncopy ) / beta / beta
  coeff_G2 = 4.d0 * coeff / beta

  do iatom = istart, iend 
     if(cnum(iatom).eq.0) then
        x_centroid(1:3,iatom) = x(1:3,iatom)
     else if( cnum(iatom).eq.1) then
        x_centroid(1:3,iatom) = 0.d0
        do icopy=1,ncopy
           x_centroid(1:3,iatom)=x_centroid(1:3,iatom)+x(1:3,iatom+icopy-1)
        end do
        x_centroid(1:3,iatom) = x_centroid(1:3,iatom)/ncopy
        do icopy=2,ncopy
           x_centroid(1:3,iatom+icopy-1)=x_centroid(1:3,iatom)
        end do
     end if
  end do 

!  +----------------------------------------------------------------------+
!  |  Virial estimator for kinetic energy                                 |
!  |                                                                      |
!  |  K = f/(2B) + 1/P sum(s=1,P)[ 1/2 (x_s - x_c) . dV(x_s)/dx_s ]       |
!  +----------------------------------------------------------------------+
        
  do iatom = istart,iend
     if(cnum(iatom).eq.0) then
        ! Contribution from classical atoms to the estimator of dV/dl 
        ! for TI w.r.t. mass. (Same for thermodynamic and virial estimators.)
        if (itimass > 0) dvdl = dvdl - 1.5d0 * kT * dmdlm(iatom)
        cycle
     else
     do idim    = 1, 3
        tmp = x(idim, iatom) - x_centroid(idim,iatom)
        sum_d = sum_d - tmp * f(idim,iatom)
        sum_v = sum_v - f(idim,iatom)*x(idim,iatom) 
        if ( cnum(iatom) == 1 ) then
           tmp_e = ( x(idim,iatom) - x( idim,iatom+1) )**2
           tmp_f = 2.0d0 * x( idim, iatom ) - x( idim, iatom+1 ) &
                                            - x( idim, iatom+ncopy-1)
        else if ( cnum(iatom) == ncopy ) then
           tmp_e = ( x(idim,iatom) - x(idim,iatom-ncopy+1))**2
           tmp_f = 2.0d0 * x(idim,iatom) - x(idim,iatom-ncopy+1 ) &
                                         - x(idim,iatom-1)
        else
           tmp_e = ( x(idim,iatom) - x(idim,iatom+1) )**2
           tmp_f = 2.0d0 * x(idim,iatom) - x(idim,iatom+1)&
                                         - x(idim,iatom-1)
        endif
        sum_e = sum_e + mass(iatom)*tmp_e

        !--------------------------------------------------------!
        ! Calculation of estimators of dV/dl for TI w.r.t. mass. !
        ! dvdl is computed from Eq. 3.16 or Eq. 3.7  in          !
        ! Vanicek & Miller, JCP, 2007                            !
        !--------------------------------------------------------!
        if (itimass > 0) then
          select case (itimass)
            case (1)   ! virial-like estimator, Eq. 3.16
              tmp_ti = kT/ncopy - f(idim,iatom) * tmp
            case (2)   ! thermodynamic-like estimator, Eq. 3.7
              tmp_ti = kT - coeff * mass(iatom) * tmp_e 
          end select
          dvdl = dvdl - dmdlm(iatom) * 0.5d0 * tmp_ti
        end if
        
        f(idim,iatom) = f(idim,iatom) - coeff*mass(iatom)*tmp_f
     enddo
     end if
  enddo

#ifdef MPI
#ifdef LES

!  call qi_corrf_les ( x, mass )
   goto 888

!  +----------------------------------------------------------+
!  |  QI dynamical factors                                    |
!  +----------------------------------------------------------+

   if( i_qi == 2 ) then
      ftmp = 1.0d0
      F_sum1 = 0.d0
      F_sum2 = 0.d0
      F_sum3 = 0.d0
      F_sum4 = 0.d0
      G_sum = 0.d0
      do n = 1, nbias
         i3 = 1
         do m = 1, natomCL
 
            ! for "imaginary-time" velocity correlation function.
            ! eq.(2.24) JCP 120, 3086 (2004).
            mm = bead_dcrypt( m, div_ndx(n) )
            if ( cnum(mm) /= 0 ) then
               if ( div_ndx(n) == ncopy ) then
                  dx(i3)   = x(1,mm-ncopy+1) - x(1,mm-1)
                  dx(i3+1) = x(2,mm-ncopy+1) - x(2,mm-1)
                  dx(i3+2) = x(3,mm-ncopy+1) - x(3,mm-1)
               else
                  dx(i3)   = x(1,mm+1) - x(1,mm-1)
                  dx(i3+1) = x(2,mm+1) - x(2,mm-1)
                  dx(i3+2) = x(3,mm+1) - x(3,mm-1)
               endif
            else
               dx(i3) = 0.d0
               dx(i3+1) = 0.d0
               dx(i3+2) = 0.d0
            endif
            i3 = i3 + 3
 
            ! for F and G factors.
            ! eq.(2.29) JCP 120, 3086 (2004).
            mm = bead_dcrypt( m, 1 )
            if ( cnum(mm) == 1 ) then
               do idim = 1, 3
                  F_sum1 = F_sum1 + mass(mm) * ( x(idim,mm) - x(idim,mm+ncopy-1) )**2
                  F_sum1 = F_sum1 + mass(mm) * ( x(idim,mm+ncopy/2) - x(idim,mm+ncopy/2-1) )**2
                  G_sum = G_sum + mass(mm) * ( x(idim,mm) - x(idim,mm+ncopy-1) )**2
               enddo
               do icopy = 1, ncopy/2-1
                  do idim = 1, 3
                     F_sum1 = F_sum1  &
                        + mass(mm) * ( x(idim,mm+icopy) - x(idim,mm+icopy-1) )**2
                     F_sum2 = F_sum2  &
                        + mass(mm) * ( x(idim,mm+icopy+ncopy/2) - x(idim,mm+icopy+ncopy/2-1) )**2
                  enddo
               enddo
               do icopy = 1, ncopy-1
                  do idim = 1, 3
                     G_sum = G_sum + mass(mm) * ( x(idim,mm+icopy) - x(idim,mm+icopy-1) )**2
                  enddo
               enddo
               do icopy = 1, ncopy/2-1
                  F_sum3 = F_sum3 + vel0_bead(icopy)
                  F_sum4 = F_sum4 + vel0_bead(icopy+ncopy/2)
               enddo
            endif
            
         enddo
         ! for "imaginary-time" velocity correlation function.
         ftmp = ftmp * dot_product( gradRC(:,n), dx(:) )
      enddo
      ! for "imaginary-time" velocity correlation function.
      fv = - 0.25d0 * dble(ncopy) * coeff * ftmp 
 
      Ftot = - coeff * ( F_sum1 - F_sum2 ) + coeff_F * ( F_sum3 -F_sum4 )
      Gtot = coeff_G1 - coeff_G2 * G_sum
 
   endif

 888 continue

#endif
#endif

  Epot_spring = 0.5d0 * coeff * sum_e
  Epot_deriv = 0.5d0 * sum_d
  Eimp_virial = sum_v

#ifdef MPI
  ener(1) = Epot_spring
  ener(2) = Epot_deriv
  ener(3) = Eimp_virial
  ener(4) = dvdl
  call mpi_allreduce( ener, temp_ener, 4, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr )
  Epot_spring = temp_ener(1)
  Epot_deriv  = temp_ener(2)
  Eimp_virial = temp_ener(3)
  dvdl = temp_ener(4)
#endif 

end subroutine



!------------------------------------------------------------------------------
subroutine pimd_full_spring_force ( x, f, mass, Epot_spring, Epot_deriv, dvdl )
! computes the harmonic-spring energy and forces arising from discretized
! path integral action.
! computes the kinetic term of the virial estimator for the internal energy
! Computes the estimators for the thermonynamic integration w.r.t. mass.
!------------------------------------------------------------------------------
  use pimd_vars, only: nbead,Eimp_virial,x_centroid, dmdlm, itimass
  use constants, only: hbar, kB   
  use full_pimd_vars, only: xall,mybeadid
  implicit none
# include "md.h"
# include "memory.h"
  _REAL_, intent(in) :: x( 3, natom ) !! coordinates
  _REAL_, intent(in) :: mass(natom)
  _REAL_ :: f( 3, natom )             !! forces
  _REAL_ :: Epot_spring, Epot_deriv, beta, kT, dvdl
  integer :: idim, iatom, istart, iend, st, ierr,ibead,prev_bead,next_bead
  integer :: prev_node, next_node, tag_prev, tag_next
  _REAL_ :: coeff, tmp, tmp_x, tmp_f, tmp_e, tmp_ti, sum_e, sum_d, sum_v

#ifdef MPI
# include  "mpif.h"
  _REAL_ :: ener(4), temp_ener(4)
#endif
#include "files.h"
#include "parallel.h"

  REQUIRE(ng_sequential)
/*
  next_node = worldrank+sandersize
  if(next_node>=worldsize) next_node = next_node-worldsize

  prev_node = worldrank-sandersize
  if(prev_node<0) prev_node = prev_node + worldsize

  tag_prev = 111
  tag_next = 222

  call mpi_sendrecv(x,3*natom,MPI_DOUBLE_PRECISION,next_node,tag_prev, &
                    xprev,3*natom,MPI_DOUBLE_PRECISION,prev_node,tag_prev, &
                    commworld,st,ierr)
  
  call mpi_sendrecv(x,3*natom,MPI_DOUBLE_PRECISION,prev_node,tag_next, &
                    xnext,3*natom,MPI_DOUBLE_PRECISION,next_node,tag_next, &
                    commworld,st,ierr)
*/

#ifdef MPI

  REQUIRE( allocated(xall) )

  if( sanderrank.eq.0) then
      call mpi_allgather(x   ,3*natom,MPI_DOUBLE_PRECISION,&
                         xall,3*natom,MPI_DOUBLE_PRECISION,&
                         commmaster,ierr)
  endif

  call mpi_bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#endif



  istart = iparpt(mytaskid)+1
  iend   = iparpt(mytaskid+1)


  kT = kB * temp0
  beta = 1.0d0/kT
  coeff = dble( nbead ) / ( hbar * beta )**2 

  sum_e = 0.d0
  sum_d = 0.d0
  sum_v = 0.d0
  dvdl = 0.d0

  x_centroid(1:3,1:natom)=0.0
  do ibead = 1, nbead
  do iatom = istart,iend
     x_centroid(1:3,iatom)=x_centroid(1:3,iatom)+xall(1:3,iatom,ibead)
  end do
  end do

  x_centroid(1:3,1:natom)=x_centroid(1:3,1:natom)/nbead   

  prev_bead = mybeadid-1
  next_bead = mybeadid+1
  
  if(prev_bead.eq.0) prev_bead=nbead
  if(next_bead.eq.nbead+1) next_bead=1

  do iatom = istart,iend
  do idim    = 1, 3
     tmp   = ( x(idim,iatom) - x_centroid(idim,iatom) )
     tmp_e = ( x(idim,iatom) - xall(idim,iatom,prev_bead) )**2
     tmp_f = 2.0d0 * x(idim,iatom) - xall(idim,iatom,prev_bead) &
                                   - xall(idim,iatom,next_bead)
     sum_e = sum_e + mass(iatom) * tmp_e
     sum_d = sum_d - f(idim,iatom) * tmp
     sum_v = sum_v - f(idim,iatom) * x(idim,iatom)

     !--------------------------------------------------------!
     ! Calculation of estimators of dV/dl for TI w.r.t. mass. !
     ! dvdl is computed from Eq. 3.16 or Eq. 3.7  in          !
     ! Vanicek & Miller, JCP, 2007                            !
     !--------------------------------------------------------!
     if (itimass > 0) then
        select case (itimass)
          case (1)   ! virial-like estimator, Eq. 3.16
            tmp_ti = kT/nbead - f(idim,iatom) * tmp
          case (2)   ! thermodynamic-like estimator, Eq. 3.7
            tmp_ti = kT - coeff * mass(iatom) * tmp_e 
        end select
        dvdl = dvdl - dmdlm(iatom) * 0.5d0 * tmp_ti
     end if

     f(idim,iatom) = f(idim,iatom) - coeff * mass(iatom)*tmp_f
  enddo
  enddo

#ifdef MPI
  ener(1) = sum_e
  ener(2) = sum_d
  ener(3) = sum_v
  ener(4) = dvdl
  call mpi_allreduce(ener,temp_ener,4,MPI_DOUBLE_PRECISION,mpi_sum,commworld,ierr)
  sum_e = temp_ener(1)
  sum_d = temp_ener(2)
  sum_v = temp_ener(3)
  dvdl = temp_ener(4)
#endif

  Epot_spring = 0.5d0 * coeff * sum_e
  Epot_deriv = 0.5d0 * sum_d
  Eimp_virial = sum_v 


end subroutine pimd_full_spring_force


subroutine dispatch( ntime, begin_in, step_in, begin_out, step_out, input, output )
   implicit none
   integer ntime, begin_in, step_in, begin_out, step_out
   integer i, id_in, id_out
   _REAL_  input(*), output(*)

   do i=1,ntime
      id_in = begin_in  + (i-1) * step_in
      id_out= begin_out + (i-1) * step_out 
      output( id_out ) = input( id_in )
   end do
end subroutine

subroutine rotate( rotat, n, xc )
   implicit none
   integer i,i3,n
   _REAL_  rotat(3,3), xc(3*n), tmpx, tmpy, tmpz

   do i = 1, n
      i3 = 3*i-3
      tmpx = rotat(1,1)*xc(i3+1) + rotat(1,2)*xc(i3+2) + rotat(1,3)*xc(i3+3)
      tmpy = rotat(2,1)*xc(i3+1) + rotat(2,2)*xc(i3+2) + rotat(2,3)*xc(i3+3)
      tmpz = rotat(3,1)*xc(i3+1) + rotat(3,2)*xc(i3+2) + rotat(3,3)*xc(i3+3)
      xc(i3+1) = tmpx
      xc(i3+2) = tmpy
      xc(i3+3) = tmpz
   end do

end subroutine 

subroutine dispatch3( ntime, begin_in, step_in, begin_out, step_out, input, output )
   implicit none
   integer ntime, begin_in, step_in, begin_out, step_out
   _REAL_ input(*), output(*)

   call dispatch( ntime, 3*begin_in - 2, 3*step_in, 3*begin_out-2, 3*step_out, input, output )
   call dispatch( ntime, 3*begin_in - 1, 3*step_in, 3*begin_out-1, 3*step_out, input, output )
   call dispatch( ntime, 3*begin_in    , 3*step_in, 3*begin_out  , 3*step_out, input, output )
end subroutine dispatch3

subroutine normalize(vector,dim)
   !normalize a vector of components of a dim dimensional vector
   !Very small vectors need to be considered zero for simulated annealing.
   !These are trapped below and not normalized.
   implicit none
    _REAL_ vector(*),vlength
    integer dim, i
   
   vlength = 0.d0
   do i = 1, dim 
      vlength = vlength + vector(i)*vector(i)
   end do
   if (vlength>1.0d-6) then
      vlength = 1.0d0/sqrt(vlength)
      vector(1:dim) = vector(1:dim)*vlength 
   else 
      vector(1:dim) = 0.d0
   end if      

end subroutine normalize

subroutine part_neb_forces( x, f, v, energy)
   use pimd_vars, only: nbead, nebrms, nrg_all
   use neb_vars, only:  tangents, springforce
   implicit none
#include "md.h" 
!need access to skmin, skmax, and tanmode
#include "memory.h" 
!need access to natom
#include "les.h"
   logical rmsok
   integer i, rep, iatom, iatm3, index
   _REAL_ eplus,eminus,spring,spring2,ezero,emax
   _REAL_ dotproduct,dotproduct2
  !external function
  ! _REAL_ ddot
   
  !passed in variables:
     !x is the coordinates
     !f is the forces
     !energy is the total energy of the system
     !amass is the atomic mass array
     !v is the velocities
   _REAL_ x(*),f(*),energy,v(*)


   !ezero is the higher energy of the two end points
   ezero = max( nrg_all(1), nrg_all(ncopy) )
   
   !emax is the highest energy of all the points
   emax = nrg_all(1)
   do rep = 2, ncopy
      emax = max(emax, nrg_all(rep))
   end do
   
   !spring2 is the spring constant of the second point in path
   if (skmin.EQ.skmax) then
      spring2 = skmax   
   else if (nrg_all(2)>ezero) then
      spring2 = skmax - skmin*((emax-max(nrg_all(1),nrg_all(2)))/  &
          (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

   energy = 0.d0
   do rep = 2, ncopy-1 
      !calculate spring constant for rep and rep+1
 
      spring = spring2
      if (skmin.EQ.skmax) then
         spring2 = skmax
      else if (nrg_all(rep+1)>ezero.AND.emax/=ezero) then
         spring2 = skmax - skmin*((emax-max(nrg_all(rep),nrg_all(rep-1)))/ &
                   (emax-ezero))
      else
         spring2 = skmax - skmin
      end if

      tangents = 0.0  

      if (tmode.EQ.1) then
         !calculate the tangents (for all images except first and last)
         if (nrg_all(rep+1)>nrg_all(rep).AND.nrg_all(rep)>nrg_all(rep-1)) then
            do iatom = 1, natom
               if(cnum(iatom).eq.rep) then
                   iatm3 = 3*(iatom - 1)
                   tangents(iatm3+1) = x(iatm3 + 4) - x(iatm3 + 1)
                   tangents(iatm3+2) = x(iatm3 + 5) - x(iatm3 + 2)
                   tangents(iatm3+3) = x(iatm3 + 6) - x(iatm3 + 3)
                end if
             end do
         else if (nrg_all(rep+1)<nrg_all(rep).AND.nrg_all(rep)<nrg_all(rep-1)) then
            do iatom = 1, natom
               if(cnum(iatom).eq.rep) then
                  iatm3 = 3*(iatom - 1)
                  tangents(iatm3+1) = x(iatm3 + 1) - x(index - 2)
                  tangents(iatm3+2) = x(index + 2) - x(index - 1)
                  tangents(iatm3+3) = x(index + 3) - x(index )
               end if
            end do  
         else if (nrg_all(rep+1)>nrg_all(rep-1)) then
            eplus = max(abs(nrg_all(rep+1)-nrg_all(rep)), &
               abs(nrg_all(rep-1)-nrg_all(rep)))
            eminus = min(abs(nrg_all(rep+1)-nrg_all(rep)), &
               abs(nrg_all(rep-1)-nrg_all(rep)))
         
            do iatom = 1, natom
               if(cnum(iatom).eq.rep) then
                  iatm3 = 3*(iatom - 1)
                  tangents(iatm3+1) = (x(iatm3 + 4) - x(iatm3 + 1))*eplus
                  tangents(iatm3+2) = (x(iatm3 + 5) - x(iatm3 + 2))*eplus
                  tangents(iatm3+3) = (x(iatm3 + 6) - x(iatm3 + 3))*eplus
                  tangents(iatm3+1) = tangents(iatm3+1) + (x(iatm3 + 1) - x(iatm3 - 2))*eminus
                  tangents(iatm3+2) = tangents(iatm3+2) + (x(iatm3 + 2) - x(iatm3 - 1))*eminus
                  tangents(iatm3+3) = tangents(iatm3+3) + (x(iatm3 + 3) - x(iatm3 ))*eminus
               end if
            end do  
         else !nrg_all(rep+1)<=nrg_all(rep-1)
            eplus = max(abs(nrg_all(rep+1)-nrg_all(rep)), &
               abs(nrg_all(rep-1)-nrg_all(rep)))
            eminus = min(abs(nrg_all(rep+1)-nrg_all(rep)), &
               abs(nrg_all(rep-1)-nrg_all(rep)))

            do iatom = 1, natom
               if(cnum(iatom).eq.rep) then
                  iatm3 = 3*(iatom - 1)
                  tangents(iatm3+1) = (x(iatm3 + 4) - x(iatm3 + 1))*eminus
                  tangents(iatm3+2) = (x(iatm3 + 5) - x(iatm3 + 2))*eminus
                  tangents(iatm3+3) = (x(iatm3 + 6) - x(iatm3 + 3))*eminus
                  tangents(iatm3+1) = tangents(iatm3+1) + (x(iatm3 + 1) - x(iatm3 - 2))*eplus
                  tangents(iatm3+2) = tangents(iatm3+2) + (x(iatm3 + 2) - x(iatm3 - 1))*eplus
                  tangents(iatm3+3) = tangents(iatm3+3) + (x(iatm3 + 3) - x(iatm3 ))*eplus
               endif
            end do  
         end if 
      else !tmode.NE.1 so use strict tangents definition
         do iatom = 1, natom
            if(cnum(iatom).eq.rep) then
               iatm3 = 3*(iatom - 1)
               tangents(iatm3+1) = x(iatm3 + 4) - x(iatm3 + 1)
               tangents(iatm3+2) = x(iatm3 + 5) - x(iatm3 + 2)
               tangents(iatm3+3) = x(iatm3 + 6) - x(iatm3 + 3)
            end if
         end do
      end if
      call normalize(tangents, 3 * natom)
!2451 format('spring = ',e10.3)
!2452 format('spring2 = ',e10.3)

      dotproduct = 0.d0
      dotproduct2 = 0.d0
  
      do iatom = 1, natom
         if(cnum(iatom).eq.rep) then
            iatm3 = 3*(iatom - 1)
            dotproduct = dotproduct + f(index+1)*tangents(iatm3+1)
            dotproduct = dotproduct + f(index+2)*tangents(iatm3+2)
            dotproduct = dotproduct + f(index+3)*tangents(iatm3+3)

            springforce(iatm3+1) = (x(iatm3+4) - x(iatm3+1))*spring2 - &
                                   (x(iatm3+1) - x(iatm3-2))*spring
            springforce(iatm3+2) = (x(iatm3+5) - x(iatm3+2))*spring2 - &
                                   (x(iatm3+2) - x(iatm3-1))*spring
            springforce(iatm3+3) = (x(iatm3+6) - x(iatm3+3))*spring2 - &
                                   (x(iatm3+3) - x(iatm3))*spring

            energy = energy + 0.5*spring2*(x(iatm3+4)-x(iatm3+1))*(x(iatm3+4)-x(iatm3+1))
            energy = energy + 0.5*spring2*(x(iatm3+5)-x(iatm3+2))*(x(iatm3+5)-x(iatm3+2))
            energy = energy + 0.5*spring2*(x(iatm3+6)-x(iatm3+3))*(x(iatm3+6)-x(iatm3+3))

            energy = energy + 0.5*spring*(x(iatm3+1)-x(iatm3-2))*(x(iatm3+1)-x(iatm3-2))
            energy = energy + 0.5*spring*(x(iatm3+2)-x(iatm3-1))*(x(iatm3+2)-x(iatm3-1))
            energy = energy + 0.5*spring*(x(iatm3+3)-x(iatm3  ))*(x(iatm3+3)-x(iatm3  ))

            dotproduct2 = dotproduct2 + springforce(iatm3+1)*tangents(iatm3+1)
            dotproduct2 = dotproduct2 + springforce(iatm3+2)*tangents(iatm3+2)
            dotproduct2 = dotproduct2 + springforce(iatm3+3)*tangents(iatm3+3)
         end if
      end do

      do iatom = 1, natom
         if(cnum(iatom).eq.rep) then
            iatm3 = 3*(iatom - 1)
            f(iatm3+1) = f(iatm3+1) - dotproduct*tangents(iatm3+1)
            f(iatm3+2) = f(iatm3+2) - dotproduct*tangents(iatm3+2)         
            f(iatm3+3) = f(iatm3+3) - dotproduct*tangents(iatm3+3)

            f(iatm3+1) = f(iatm3+1) + dotproduct2*tangents(iatm3+1)
            f(iatm3+2) = f(iatm3+2) + dotproduct2*tangents(iatm3+2)         
            f(iatm3+3) = f(iatm3+3) + dotproduct2*tangents(iatm3+3)
         end if
      end do

   end do !rep
   
   dotproduct = 0.d0
   dotproduct2= 0.d0
   
   !note that there are fewer degrees of freedom for the neb case
   do index = 1,3*natom
      dotproduct = dotproduct + f(index)*f(index)
      dotproduct2= dotproduct2+ v(index)*v(index)
   enddo
end subroutine part_neb_forces


subroutine full_neb_forces(mass,x,f,v,epot,energy)
   use pimd_vars, only: nbead, natomCL, nebrms,nrg_all
   use full_pimd_vars, only: xall,mybeadid
   use neb_vars, only: fitgroup, rmsgroup, springforce, tangents,neb_rotat
   implicit none
#include "md.h" 
!need access to skmin, skmax, and tanmode
#include "memory.h" 
!need access to natom
#ifdef MPI
# include "mpif.h"
# include "parallel.h"
  integer ierr
#endif
   _REAL_ mass(*)
   logical rmsok
   integer i, rep, iatom, iatm3, index
   _REAL_ eplus,eminus,spring,spring2,ezero,emax
   _REAL_ dotproduct,dotproduct2
   _REAL_ rmsdvalue, rotat(3,3)
   _REAL_ dottmp(1)
 !external function
  ! _REAL_ ddot
   
  !passed in variables:
     !x is the coordinates
     !f is the forces
     !energy is the total energy of the system
     !amass is the atomic mass array
     !v is the velocities
   _REAL_ x(3,*),f(*),epot,v(*)
   _REAL_ energy


#ifdef MPI
  if( sanderrank.eq.0) then
      call mpi_allgather(x   ,3*natom,MPI_DOUBLE_PRECISION,&
                         xall,3*natom,MPI_DOUBLE_PRECISION,&
                         commmaster,ierr)
      call mpi_allgather(epot,1,MPI_DOUBLE_PRECISION, &
                         nrg_all,1, MPI_DOUBLE_PRECISION, &
                         commmaster,ierr)
  endif

  call mpi_bcast(f,3*natom,MPI_DOUBLE_PRECISION,0,commsander,ierr)
  call mpi_bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,0,commsander,ierr)
  call mpi_bcast(nrg_all,nbead,MPI_DOUBLE_PRECISION,0,commsander,ierr)

#endif


   nattgtfit = natomCL
   nattgtrms = natomCL

   do i = 1, natomCL
      fitgroup(i) = i
      rmsgroup(i) = i
   end do
   !start with rms fit


   do rep = 2, nbead
      call rmsfit( xall(1,1,rep), xall(1,1,rep-1), mass, fitgroup, rmsgroup, rmsdvalue, rmsok) 

      if( rep.eq.mybeadid) then
          !x(1:3,1:natom)=xall(1:3,1:natom,rep)
          do iatom=1,natom
             x(1,iatom) = xall(1,iatom,rep)
             x(2,iatom) = xall(2,iatom,rep)
             x(3,iatom) = xall(3,iatom,rep)
          end do

          call rotate( neb_rotat, natom, f )
          call rotate( neb_rotat, natom, v )
      end if
   end do

   !ezero is the higher energy of the two end points
   ezero = max( nrg_all(1), nrg_all(nbead) )
   
   !emax is the highest energy of all the points
   emax = nrg_all(1)
   do rep = 2, nbead
      emax = max(emax, nrg_all(rep))
   end do
   
   !spring2 is the spring constant of the second point in path
   if (skmin.EQ.skmax) then
      spring2 = skmax   
   else if (nrg_all(2)>ezero) then
      spring2 = skmax - skmin*((emax-max(nrg_all(1),nrg_all(2)))/  &
          (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

   energy = 0.d0
   rep = mybeadid

   nebrms = 0.d0

   if( mybeadid==1.or.mybeadid==nbead) return
   !calculate spring constant for rep and rep+1
 
   spring = spring2
   if (skmin.EQ.skmax) then
      spring2 = skmax
   else if (nrg_all(rep+1)>ezero.AND.emax/=ezero) then
      spring2 = skmax - skmin*((emax-max(nrg_all(rep),nrg_all(rep-1)))/ &
                (emax-ezero))
   else
      spring2 = skmax - skmin
   end if

   tangents = 0.0  


   if (tmode.EQ.1) then
      !calculate the tangents (for all images except first and last)
      if (nrg_all(rep+1)>nrg_all(rep).AND.nrg_all(rep)>nrg_all(rep-1)) then
         do iatom = 1, natom
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = xall(1,iatom,rep+1) - xall(1,iatom,rep)
            tangents(iatm3+2) = xall(2,iatom,rep+1) - xall(2,iatom,rep)
            tangents(iatm3+3) = xall(3,iatom,rep+1) - xall(3,iatom,rep)
         end do
      else if (nrg_all(rep+1)<nrg_all(rep).AND.nrg_all(rep)<nrg_all(rep-1)) then
         do iatom = 1, natomCL
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = xall(1,iatom,rep) - xall(1,iatom,rep-1)
            tangents(iatm3+2) = xall(2,iatom,rep) - xall(2,iatom,rep-1)
            tangents(iatm3+3) = xall(3,iatom,rep) - xall(3,iatom,rep-1)
         end do  
      else if (nrg_all(rep+1)>nrg_all(rep-1)) then
         eplus = max( abs(nrg_all(rep+1)-nrg_all(rep)), &
            abs(nrg_all(rep-1)-nrg_all(rep)) )
         eminus = min(abs(nrg_all(rep+1)-nrg_all(rep)), &
            abs(nrg_all(rep-1)-nrg_all(rep)) )
         
         do iatom = 1, natomCL
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = (xall(1,iatom,rep+1) - xall(1,iatom,rep))*eplus
            tangents(iatm3+2) = (xall(2,iatom,rep+1) - xall(2,iatom,rep))*eplus
            tangents(iatm3+3) = (xall(3,iatom,rep+1) - xall(3,iatom,rep))*eplus
            tangents(iatm3+1) = tangents(iatm3+1) + (xall(1,iatom,rep) - xall(1,iatom,rep-1))*eminus
            tangents(iatm3+2) = tangents(iatm3+2) + (xall(2,iatom,rep) - xall(2,iatom,rep-1))*eminus
            tangents(iatm3+3) = tangents(iatm3+3) + (xall(3,iatom,rep) - xall(3,iatom,rep-1))*eminus
         end do  
      else !nrg_all(rep+1)<=nrg_all(rep-1)
         eplus = max(abs(nrg_all(rep+1)-nrg_all(rep)), &
            abs(nrg_all(rep-1)-nrg_all(rep)))
         eminus = min(abs(nrg_all(rep+1)-nrg_all(rep)), &
            abs(nrg_all(rep-1)-nrg_all(rep)))
         do iatom = 1, natomCL
            iatm3 = 3*(iatom - 1)
            tangents(iatm3+1) = (xall(1,iatom,rep+1) - xall(1,iatom,rep))*eminus
            tangents(iatm3+2) = (xall(2,iatom,rep+1) - xall(2,iatom,rep))*eminus
            tangents(iatm3+3) = (xall(3,iatom,rep+1) - xall(3,iatom,rep))*eminus
            tangents(iatm3+1) = tangents(iatm3+1) + (xall(1,iatom,rep) - xall(1,iatom,rep-1))*eplus
            tangents(iatm3+2) = tangents(iatm3+2) + (xall(2,iatom,rep) - xall(2,iatom,rep-1))*eplus
            tangents(iatm3+3) = tangents(iatm3+3) + (xall(3,iatom,rep) - xall(3,iatom,rep-1))*eplus
         end do  
      end if 
   else !tmode.NE.1 so use strict tangents definition
      do iatom = 1, natomCL
         iatm3 = 3*(iatom - 1)
         tangents(iatm3+1) = xall(1,iatom,rep+1) - xall(1,iatom,rep)
         tangents(iatm3+2) = xall(2,iatom,rep+1) - xall(2,iatom,rep)
         tangents(iatm3+3) = xall(3,iatom,rep+1) - xall(3,iatom,rep)
      end do
   end if

   call normalize(tangents, 3 * natomCL)
!2451 format('spring = ',e10.3)
!2452 format('spring2 = ',e10.3)

   dotproduct = 0.d0
   dotproduct2 = 0.d0
  
   do iatom = 1, natomCL
      iatm3 = 3*(iatom - 1)
      dotproduct = dotproduct + f(iatm3+1)*tangents(iatm3+1)
      dotproduct = dotproduct + f(iatm3+2)*tangents(iatm3+2)
      dotproduct = dotproduct + f(iatm3+3)*tangents(iatm3+3)

      springforce(iatm3+1) = (xall(1,iatom,rep+1) - xall(1,iatom,rep))*spring2 - &
                             (xall(1,iatom,rep) - xall(1,iatom,rep-1))*spring
      springforce(iatm3+2) = (xall(2,iatom,rep+1) - xall(2,iatom,rep))*spring2 - &
                             (xall(2,iatom,rep) - xall(2,iatom,rep-1))*spring
      springforce(iatm3+3) = (xall(3,iatom,rep+1) - xall(3,iatom,rep))*spring2 - &
                             (xall(3,iatom,rep) - xall(3,iatom,rep-1))*spring
      energy = energy + 0.5*spring2*(xall(1,iatom,rep+1)-xall(1,iatom,rep))*(xall(1,iatom,rep+1)-xall(1,iatom,rep))
      energy = energy + 0.5*spring2*(xall(2,iatom,rep+1)-xall(2,iatom,rep))*(xall(2,iatom,rep+1)-xall(2,iatom,rep))
      energy = energy + 0.5*spring2*(xall(3,iatom,rep+1)-xall(3,iatom,rep))*(xall(3,iatom,rep+1)-xall(3,iatom,rep))

      energy = energy + 0.5*spring*(xall(1,iatom,rep)-xall(1,iatom,rep-1))*(xall(1,iatom,rep)-xall(1,iatom,rep-1))
      energy = energy + 0.5*spring*(xall(2,iatom,rep)-xall(2,iatom,rep-1))*(xall(2,iatom,rep)-xall(2,iatom,rep-1))
      energy = energy + 0.5*spring*(xall(3,iatom,rep)-xall(3,iatom,rep-1))*(xall(3,iatom,rep)-xall(3,iatom,rep-1))

      dotproduct2 = dotproduct2 + springforce(iatm3+1)*tangents(iatm3+1)
      dotproduct2 = dotproduct2 + springforce(iatm3+2)*tangents(iatm3+2)
      dotproduct2 = dotproduct2 + springforce(iatm3+3)*tangents(iatm3+3)
   end do

   do iatom = 1, natomCL
      iatm3 = 3*(iatom - 1)
      f(iatm3+1) = f(iatm3+1) - dotproduct*tangents(iatm3+1)
      f(iatm3+2) = f(iatm3+2) - dotproduct*tangents(iatm3+2)         
      f(iatm3+3) = f(iatm3+3) - dotproduct*tangents(iatm3+3)

      f(iatm3+1) = f(iatm3+1) + dotproduct2*tangents(iatm3+1)
      f(iatm3+2) = f(iatm3+2) + dotproduct2*tangents(iatm3+2)         
      f(iatm3+3) = f(iatm3+3) + dotproduct2*tangents(iatm3+3)
   end do

   dotproduct = 0.d0
   dotproduct2= 0.d0
   
   !note that there are fewer degrees of freedom for the neb case
   do iatm3 = 1,3*natom
      dotproduct = dotproduct + f(iatm3)*f(iatm3)
      dotproduct2= dotproduct2+ v(iatm3)*v(iatm3)
   enddo

   nebrms = dotproduct
end subroutine


subroutine pimd_neb_energy_report(filenumber)
   use pimd_vars, only: nebrms, nrg_all,nbead
   implicit none
   integer filenumber,ix
   _REAL_ sum
   
   sum = 0
   do ix = 1, nbead
      sum = sum + nrg_all(ix)
   enddo

   write(filenumber,'(a)') "NEB replicate breakdown:"
   do ix = 1, nbead
      write(filenumber,'(a,i3,a,f13.4)') "Energy for replicate ",ix," = ",nrg_all(ix)
   enddo
   write(filenumber,'(a,f13.4)') "Total Energy of replicates = ",sum
   write(filenumber,'(a,f13.6)') "NEB RMS = ",nebrms

end subroutine pimd_neb_energy_report





! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  QI correlation functions (LES-PIMD)                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine qi_corrf_les ( x, mass )

#ifdef MPI
#ifdef LES

   use constants, only: hbar, kB
   use evb_pimd,  only: natomCL, vel0_bead, bead_dcrypt
   use miller,    only: i_qi, gradRC, div_ndx
   use evb_parm,  only: nbias
   use evb_data,  only: f_v, F_QI, G_QI

   implicit none

#  include "md.h"
#  include "les.h"
#  include "memory.h"
#  include "mpif.h"
#  include "parallel.h"

  _REAL_, intent(in) :: x(3,natom) !! coordinates
  _REAL_, intent(in) :: mass(natom)

   !  ..........................................................................

  integer :: idim, icopy
  integer :: m, n, mm, i3

  _REAL_ :: dx(3*natomCL)
  _REAL_ :: beta, coeff 
  _REAL_ :: coeff_F, coeff_G1, coeff_G2
  _REAL_ :: F_sum1, F_sum2, F_sum3, F_sum4, G_sum, Ftot, Gtot
  _REAL_ :: ftmp, fv
  _REAL_ , intrinsic :: sqrt, dot_product

!  +----------------------------------------------------------+
!  |  QI dynamical factors                                    |
!  +----------------------------------------------------------+

   if( i_qi == 2 ) then

      beta = 1.0d0 / ( kB * temp0 )
      coeff = dble( ncopy ) / ( hbar * beta )**2

      coeff_F  = 2.0d0 / dble( ncopy )
      coeff_G1 = 6.0d0 * dble( ncopy ) / beta / beta
      coeff_G2 = 4.0d0 * coeff / beta

      ftmp = 1.0d0
      F_sum1 = 0.d0
      F_sum2 = 0.d0
      F_sum3 = 0.d0
      F_sum4 = 0.d0
      G_sum = 0.d0
      do n = 1, nbias
         i3 = 1
         do m = 1, natomCL
 
            ! for "imaginary-time" velocity correlation function.
            ! eq.(2.24) JCP 120, 3086 (2004).
!  +---------------------------------------------------------------------------+
!  |  f_v [ Eq. 2.24, JCP 120, 3086 (2004) ]                                   |
!  |                                                                           |
!  |  f_v(x_1,x_2,...,x_P) = m(iP/2hbarB)^2 grad_RC(x_P) . [x_1 - x_(P-1)]     |
!  |                         grad_RC(x_(P/2)) . [x_(P/2+1) - x_(P/2-1)]        |
!  +---------------------------------------------------------------------------+

            mm = bead_dcrypt( m, div_ndx(n) )
            if ( cnum(mm) /= 0 ) then
               if ( div_ndx(n) == ncopy ) then
                  dx(i3)   = x(1,mm-ncopy+1) - x(1,mm-1)
                  dx(i3+1) = x(2,mm-ncopy+1) - x(2,mm-1)
                  dx(i3+2) = x(3,mm-ncopy+1) - x(3,mm-1)
               else
                  dx(i3)   = x(1,mm+1) - x(1,mm-1)
                  dx(i3+1) = x(2,mm+1) - x(2,mm-1)
                  dx(i3+2) = x(3,mm+1) - x(3,mm-1)
               endif
            else
               dx(i3) = 0.d0
               dx(i3+1) = 0.d0
               dx(i3+2) = 0.d0
            endif

            i3 = i3 + 3
 
!  +---------------------------------------------------------------------------+
!  |  F and G factors [ Eq. 2.29, JCP 120, 3086 (2004) ]                       |
!  |                                                                           |
!  |  F(x_1,x_2,...,x_P) = -mP/hbar^2/B^2 {sum(1,P/2) - sum(P/2+1,P)}          |
!  |       (x_k - x_k-1)^2 + 2/P {sum(1,P/2) - sum(P/2+1,P)} V(x_k)            |
!  |                                                                           |
!  |  G(x_1,x_2,...,x_P) = 2dP/B^2 - 4mP/hbar^2/B^3 sum(1,P) [x_k - x_(k-1)]^2 |
!  +---------------------------------------------------------------------------+

            mm = bead_dcrypt( m, 1 )
            if ( cnum(mm) == 1 ) then
               do idim = 1, 3
                  F_sum1 = F_sum1 + mass(mm) * ( x(idim,mm) - x(idim,mm+ncopy-1) )**2
                  F_sum1 = F_sum1 + mass(mm) * ( x(idim,mm+ncopy/2) - x(idim,mm+ncopy/2-1) )**2
                  G_sum = G_sum + mass(mm) * ( x(idim,mm) - x(idim,mm+ncopy-1) )**2
               enddo
               do icopy = 1, ncopy/2-1
                  do idim = 1, 3
                     F_sum1 = F_sum1  &
                        + mass(mm) * ( x(idim,mm+icopy) - x(idim,mm+icopy-1) )**2
                     F_sum2 = F_sum2  &
                        + mass(mm) * ( x(idim,mm+icopy+ncopy/2) - x(idim,mm+icopy+ncopy/2-1) )**2
                  enddo
               enddo
               do icopy = 1, ncopy-1
                  do idim = 1, 3
                     G_sum = G_sum + mass(mm) * ( x(idim,mm+icopy) - x(idim,mm+icopy-1) )**2
                  enddo
               enddo
               do icopy = 1, ncopy/2-1
                  F_sum3 = F_sum3 + vel0_bead(icopy)
                  F_sum4 = F_sum4 + vel0_bead(icopy+ncopy/2)
               enddo
            endif
            
         enddo
         ! for "imaginary-time" velocity correlation function.
         ftmp = ftmp * dot_product( gradRC(:,n), dx(:) )

!        write(6,*) '>>> dx(:) = ', dx(:)
!        write(6,*) '>>> gradRC = ', gradRC(:,n), n

      enddo
      ! for "imaginary-time" velocity correlation function.
      fv = - 0.25d0 * dble(ncopy) * coeff * ftmp 
 
      Ftot = - coeff * ( F_sum1 - F_sum2 ) + coeff_F * ( F_sum3 -F_sum4 )
      Gtot = coeff_G1 - coeff_G2 * G_sum

      f_v  = fv
      F_QI = Ftot
      G_QI = Gtot 

!     write(6,*) '<<< coeff, F_sum1, F_sum2 = ', coeff, F_sum1, F_sum2
!     write(6,*) '<<< coeff_F, F_sum3, F_sum4 = ', coeff_F, F_sum3, F_sum4
!     write(6,*) '<<< vel0_bead = ', vel0_bead(:)

!   write(6,*) '>>> ncopy = ', ncopy

!   write(6,*) '>>> f_v = ', fv, coeff, ftmp
!   write(6,*) '>>>   F = ', F_QI
!   write(6,*) '>>>   G = ', Gtot
   

   endif

#endif
#endif

   end subroutine qi_corrf_les



