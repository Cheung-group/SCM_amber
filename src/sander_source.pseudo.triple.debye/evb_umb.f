! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_UMB                                                                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "dprec.h"

   subroutine evb_umb ( f, q, mass, natom, istart3, iend3 )

   use constants, only: kB
   use evb_parm,  only: k_umb, r0_umb, evb_dyn, nbias, dbonds_RC, bond_RC &
                      , out_RCdot
   use evb_data,  only: evb_frc, evb_bias, RCdot
   use evb_check, only: full_evb_debug, dbonds_debug, bond_debug 
#ifdef LES
   use evb_pimd,  only: bead_dcrypt, natomCL
   use miller,    only: gradRC, i_qi, div_ndx
#endif

   implicit none

#include "md.h"

   integer, intent(   in) :: natom, istart3, iend3
   _REAL_ , intent(   in) :: q(natom*3)
   _REAL_ , intent(   in) :: mass(natom)
   _REAL_ , intent(inout) :: f(natom*3)

   !  ..........................................................................

   integer :: m, n, mm, nn, m_dx, mm_dx
   integer :: idx, jdx, kdx
   _REAL_  :: fharm(natom*3), dr(3), fr(3), rij, rkj, rij_inv, rkj_inv
   _REAL_  :: evb_fbias(natom*3,nbias), RC, pi
   _REAL_   , intrinsic :: sqrt, acos
   character, intrinsic :: trim, adjustl


   select case( trim( adjustl( evb_dyn) ) )

!  +---------------------------------------------------------------------------+
!  |  Difference of 2 bonds RC harmonic umbrella sampling                      |
!  |                                                                           |
!  |  /\ = r_ij - r_kj                                                         |
!  |                                                                           |
!  |  V'(/\) = V_0 + 0.5 k_evb * ( /\ - /\_0 )^2                               |
!  |                                                                           |
!  |  dV' / dR = ( dV' / d/\ ) * ( d/\ / dR )                                  |
!  |           = dV_0 / dR + k_evb * ( /\ - /\_0 ) * d/dR ( r_ij - r_kj )      |
!  +---------------------------------------------------------------------------+

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )

         fharm(:) = 0.0d0

         do n = 1, nbias

            idx = ( dbonds_RC(n)%iatom - 1 ) * 3
            jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
            kdx = ( dbonds_RC(n)%katom - 1 ) * 3

            do nn = 1, 3
               dr(nn) = q(idx+nn) - q(jdx+nn)
            enddo 

            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(idx+nn) = fharm(idx+nn) - fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
            enddo 

            do nn = 1, 3 
               dr(nn)  = q(kdx+nn) - q(jdx+nn)
            enddo

            rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rkj_inv = 1.0d0 / rkj

            do nn = 1, 3
               fr(nn) = dr(nn) * rkj_inv
               fharm(kdx+nn) = fharm(kdx+nn) + fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) - fr(nn)
            enddo

#ifdef LES
            if( i_qi > 0 ) then
               do m = 1, natomCL
                  mm = bead_dcrypt( m, div_ndx(n) )
                   m_dx = (  m - 1 ) * 3
                  mm_dx = ( mm - 1 ) * 3
                  gradRC(m_dx+1,n) = fharm(mm_dx+1)
                  gradRC(m_dx+2,n) = fharm(mm_dx+2)
                  gradRC(m_dx+3,n) = fharm(mm_dx+3)
               enddo
            endif
#endif
            RC = rij - rkj

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
            evb_fbias(:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:)
#ifdef DEBUG_EVB
            if( dbonds_debug .or. full_evb_debug ) &
               call dbonds_anal2num ( q, evb_fbias(:,n), natom*3 )
#endif

         enddo

         do n = 1, nbias
            do m = istart3, iend3
               f(m) = f(m) + evb_fbias(m,n)
            enddo
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
         enddo

!  +---------------------------------------------------------------------------+
!  |  Bond RC harmonic umbrella sampling                                       |
!  |                                                                           |
!  |  /\ = r_ij                                                                |
!  |                                                                           |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                            |
!  |                                                                           |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                            |
!  |              = dV_0 / dR + k_evb * ( /\ - eta ) * d/dR ( r_ij )           |
!  +---------------------------------------------------------------------------+

      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )

         fharm(:) = 0.0d0

         do n = 1, nbias

            idx = ( bond_RC(n)%iatom - 1 ) * 3
            jdx = ( bond_RC(n)%jatom - 1 ) * 3

            do nn = 1, 3
               dr(nn) = q(idx+nn) - q(jdx+nn)
            enddo

            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(idx+nn) = fharm(idx+nn) - fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
            enddo

#ifdef LES
            if( i_qi > 0 ) then
               do m = 1, natomCL
                  mm = bead_dcrypt( m, div_ndx(n) )
                   m_dx = (  m - 1 ) * 3
                  mm_dx = ( mm - 1 ) * 3
                  gradRC(m_dx+1,n) = fharm(mm_dx+1)
                  gradRC(m_dx+2,n) = fharm(mm_dx+2)
                  gradRC(m_dx+3,n) = fharm(mm_dx+3)
               enddo
            endif
#endif

            RC = rij 

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
            evb_fbias(:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:)
#ifdef DEBUG_EVB
            if( bond_debug .or. full_evb_debug ) &
               call bond_anal2num ( q, evb_fbias(:,n), natom*3 )
#endif

         enddo

         do n = 1, nbias
            do m = istart3, iend3
               f(m) = f(m) + evb_fbias(m,n)
            enddo
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
         enddo

   end select

   if( out_RCdot ) then

      pi = acos( -1.0d0 )
      RCdot = 0.0d0
      do n = 1, natom
         idx = ( n - 1 ) * 3
         RCdot = RCdot + ( fharm(idx+1)**2 + fharm(idx+2)**2 &
                         + fharm(idx+3)**2 ) / mass(n)
      enddo

      RCdot = sqrt( 2.0d0 * kB * temp0 * RCdot / pi )

   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_umb


