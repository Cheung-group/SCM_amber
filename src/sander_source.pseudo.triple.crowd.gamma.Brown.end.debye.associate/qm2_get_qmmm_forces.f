! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_get_qmmm_forces(dxyzqm,qm_xcrd,dxyzmm)

! This routine calculates force on the MM atoms due to the QM atoms
! and vice versa. The forces are added to dxyzmm and dxyzqm respectively.

! Currently only analytical derivatives are supported.
                                                                                                                                                             
!
!     Variable Definitions:
!
!     qm_coords - Cartesian coordinates of QM atoms.
!     dxyzqm - Coupled potential energy derivatives with respect to
!              movement of QM atoms.
!     qm_xcrd - Cartesian coordinates of QM and MM atoms. In same order as cutoff list.
!               also contains scaled charge in 4th element.
!     dxyzmm - Coupled potential energy derivatives with respect to
!              movement of MM atoms.
!     Note this routine does not need iqmatoms since it gets mm
!     atoms from the pair list and puts QM forces in it's own
!     array which gets moved to the main force array later.

! Current Code by Ross Walker (TSRI, 2004)
! Derivative optimisations by Ross Walker and Mike Crowley (TSRI, 2004)

      use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, qmmm_mpi

      implicit none

!Passed in
      _REAL_ , intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)
      _REAL_ , intent(in) :: qm_xcrd(4,qmmm_struct%qm_mm_pairs)
      _REAL_ , intent(out) :: dxyzmm(*)

!Local
      _REAL_ psum(36), psum_light
      _REAL_ pair_force(3)
      _REAL_ qm_atom_coord(3), fqm(3)
      _REAL_ qm_atom_core, qm_atom_alpa
      integer jj, jf, jl, ij, i, k, ii, j
      integer n_atomic_orb !Number of atomic orbitals on qm atom
      integer loop_count !Keeps track of number of times through nquant * ni_mm loop
      integer inner_loop_count !xyz offset for loop over pairs.
      logical heavy_atom
      
      loop_count=0

      !do jj=1,qmmm_struct%nquant_nlink
      do jj=qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
         jf=qm2_params%orb_loc(1,jj)
         jl=qm2_params%orb_loc(2,jj)
         qm_atom_coord(1) = qmmm_struct%qm_coords(1,jj)
         qm_atom_coord(2) = qmmm_struct%qm_coords(2,jj)
         qm_atom_coord(3) = qmmm_struct%qm_coords(3,jj)
         qm_atom_core = qm2_params%core_chg(jj)
         qm_atom_alpa = qm2_params%cc_exp_params(jj)
         n_atomic_orb = qm2_params%natomic_orbs(jj)
         heavy_atom = (n_atomic_orb>1)
         fqm(1:3)=0.0d0
! Split into heavy and light atoms here - means code duplication but is good for speed.
         if (heavy_atom) then
!        Get elements of density matrix involving only orbitals centered
!        on current qm atom.
           ij = 0
           do i=jf,jl
              k = qm2_params%pascal_tri1(i) + jf - 1
              do j=jf,i
                ij = ij + 1
                k = k + 1
                psum(ij) = qm2_struct%den_matrix(k)
              end do
           end do
!RCW: We could maybe move this outside of the loop for optimisation purposes.
!        Loop over MM atoms that are in the interaction list for current
!        QM atom. We have the loop repeated twice here. Once for when we
!        have the qm_mm 1 elec repulsion integrals in memory and once for
!        when we need to calculate them on the fly.
!        This is excess code but it avoids the if(...) being inside the
!        inner loop.
!We don't have QM-MM 1 e repul in memory so calc on the fly.
           inner_loop_count=1
           do ii=1,qmmm_struct%qm_mm_pairs
              loop_count=loop_count+1
!             Get analytical derivatives with respect to x, y, and z
!             directions for the interaction of the current QM-MM pair.
              call qm2_deriv_qmmm_heavy(jj,loop_count, &
                                  psum,qm_atom_coord,qm_xcrd(1,ii),n_atomic_orb, &
                                  pair_force,qm_atom_core,qm_atom_alpa)
  
              dxyzmm(inner_loop_count) = dxyzmm(inner_loop_count) - pair_force(1)
              dxyzmm(inner_loop_count+1) = dxyzmm(inner_loop_count+1) - pair_force(2)
              dxyzmm(inner_loop_count+2) = dxyzmm(inner_loop_count+2) - pair_force(3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3
           end do
        else !If (heavy_atom)
          !Light atom - see above for comments
           k = qm2_params%pascal_tri1(jf) + jf 
           psum_light = qm2_struct%den_matrix(k)
           inner_loop_count=1
           do ii=1,qmmm_struct%qm_mm_pairs
              loop_count=loop_count+1
              call qm2_deriv_qmmm_light(jj,loop_count,psum_light,qm_atom_coord,qm_xcrd(1,ii), &
                                        pair_force,qm_atom_core,qm_atom_alpa)
              dxyzmm(inner_loop_count) = dxyzmm(inner_loop_count) - pair_force(1)
              dxyzmm(inner_loop_count+1) = dxyzmm(inner_loop_count+1) - pair_force(2)
              dxyzmm(inner_loop_count+2) = dxyzmm(inner_loop_count+2) - pair_force(3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3
           end do
        end if !if (heavy_atom)
        dxyzqm(1:3,jj) = dxyzqm(1:3,jj) + fqm(1:3) 
      end do !jj=1,qmmm_struct%nquant

      RETURN
end subroutine qm2_get_qmmm_forces

subroutine qm2_deriv_qmmm_light(iqm,loop_count,psum_light,xyz_qm,xyz_mm,pair_force, &
                         qm_atom_core,alpa)
!For light atoms
!See heavy version of routine for comments
!  Current Version: Ross Walker (TSRI, 2005)

      use qmmm_module, only :  qmmm_nml,qm2_params, qmmm_struct, qm2_struct, qm2_rij_eqns, &
                               alph_mm, EXPONENTIAL_CUTOFF
      use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, A2_TO_BOHRS2xAU_TO_EV, one, &
                            EV_TO_KCAL, AU_TO_KCAL, zero, two
      implicit none

! ON RETURN, pair_force HOLDS ANALYTICAL DERIVATIVES                                   

!Passed in
      _REAL_, intent(in) :: xyz_qm(3),xyz_mm(4),psum_light
      _REAL_, intent(out) :: pair_force(3)
      integer, intent(in) :: iqm, loop_count
      _REAL_, intent(in) ::  qm_atom_core, alpa

!Local
      _REAL_ FABX, FABY, FABZ ,FNUCX, FNUCY, FNUCZ
      _REAL_ r2, rij, onerij, rr2
      _REAL_ vec_qm_mm1, vec_qm_mm2, vec_qm_mm3
      _REAL_ c1, mm_charge
      _REAL_ ee, sqrtaee
      _REAL_ DGX, DGY, DGZ
      _REAL_ EXP1, EXP2, EXP3, EXP4
      _REAL_ qmmm_erep, anam1, oner2, temp_real, temp_real2
      integer :: i, qmitype

      vec_qm_mm1=xyz_qm(1) - xyz_mm(1)
      vec_qm_mm2=xyz_qm(2) - xyz_mm(2)
      vec_qm_mm3=xyz_qm(3) - xyz_mm(3)
      mm_charge = xyz_mm(4)
#include "qm2_array_locations.h"
      if (qmmm_nml%qmmmrij_incore) then
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RIJ=qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        EXP1=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)
        EXP2=qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        SQRTAEE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
      else
        r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
        RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
        oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
!        RIJ=one/oneRIJ
        RIJ=R2*oneRIJ
        EXP1 = exp(-alpa*rij)
        EXP2 = exp(-ALPH_MM*rij)
        SQRTAEE=one/sqrt(RR2+qm2_params%multip_2c_elec_params(3,iqm)**2)
      end if

      qmmm_erep = AU_TO_EV*SQRTAEE

      EE=-A2_TO_BOHRS2xAU_TO_EV*SQRTAEE*SQRTAEE*SQRTAEE
!        S-orbital of QM atom:
      DGX=vec_qm_mm1*EE
      DGY=vec_qm_mm2*EE
      DGZ=vec_qm_mm3*EE

! BEGIN --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
!           Note: This could technically be done at the same time we do the
!                 energy in hcore_qmmm.
      C1=qm_atom_core*mm_charge

!     ****   START OF THE AM1 and PM3 RM1 etc SPECIFIC DERIVATIVE CODE   ***
!     ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!     ANALYTICAL DERIVATIVES

      if ( qmmm_nml%qmmm_int==2 .and. qmmm_struct%AM1_OR_PM3 ) then
        qmitype = qmmm_struct%qm_atom_type(iqm)
        ANAM1=zero
        oner2=oneRIJ*oneRIJ
        do I=1,qm2_params%num_fn(qmitype)
          temp_real=RIJ-qm2_params%FN3(i,qmitype)
          temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
          if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
            ANAM1=ANAM1+qm2_params%FN1(i,qmitype)* &
                 (oner2+two*qm2_params%FN2(i,qmitype)*temp_real*oneRIJ)*EXP(-temp_real2)
          end if
        end do
        ANAM1=-ANAM1*c1*onerij
        FNUCX=ANAM1*vec_qm_mm1
        FNUCY=ANAM1*vec_qm_mm2
        FNUCZ=ANAM1*vec_qm_mm3
     else
        FNUCX=zero; FNUCY=zero; FNUCZ=zero
     endif

!     ****   END OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***

      EXP3 = (EXP1+EXP2)*abs(c1)
      EXP4 = qmmm_erep*onerij*(alpa*EXP1 + ALPH_MM*EXP2)*abs(c1)
      FNUCX = FNUCX+dgx*c1-vec_qm_mm1*EXP4+dgX*EXP3
      FNUCY = FNUCY+dgy*c1-vec_qm_mm2*EXP4+dgy*EXP3
      FNUCZ = FNUCZ+dgz*c1-vec_qm_mm3*EXP4+dgz*EXP3

! END --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM

!     MM CORE AFFECTING AO'S ON QM ATOM.
      mm_charge=-mm_charge*psum_light
      FABX=mm_charge*DGX
      FABY=mm_charge*DGY
      FABZ=mm_charge*DGZ                    

      pair_force(1) = (FABX+FNUCX)*EV_TO_KCAL                                          
      pair_force(2) = (FABY+FNUCY)*EV_TO_KCAL                                           
      pair_force(3) = (FABZ+FNUCZ)*EV_TO_KCAL                                           


      RETURN                                                                    
end subroutine qm2_deriv_qmmm_light

subroutine qm2_deriv_qmmm_heavy(iqm,loop_count,psum,xyz_qm,xyz_mm,n_atomic_orb,pair_force, &
                         qm_atom_core,alpa)

!     This routine computes the analytical energy derivatives for the QM-MM
!     interaction energy arising from a single QM-MM pair.  The contibutions
!     to the derivatives come from the electron-core and core-core interactions.

!     Variable Definitions:
!
!     psum   - Density matrix elements for orbitals centered on QM atom.
!     xyz_qm - Cartesian coordinates of QM atom.
!     xyz_mm - Cartesian coordinates of MM atom. and charge
! n_atomic_orb - Number of atomic orbitals
! pair_force - Energy derivatives in the x, y, and z directions for
!              the interaction of the QM-MM atom pair.  The algebraic
!              sign of pair_force corresponds to dE/dR for the MM atom and
!              -dE/dR for the QM atom.
!
!  Current Version: Ross Walker (TSRI, 2004)
      use qmmm_module, only :  qmmm_nml,qm2_params, qm2_struct, qm2_rij_eqns, qmmm_struct, &
                               MAX_VALENCE_DIMENSION, alph_mm, AXIS_TOL, EXPONENTIAL_CUTOFF
      use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, HALF_AU_TO_EV, FOURTH_AU_TO_EV, &
                            A2_TO_BOHRS2xAU_TO_EV, one, zero, four, two, half, EV_TO_KCAL
      implicit none

! ON RETURN, pair_force HOLDS ANALYTICAL DERIVATIVES                                   

!Passed in
      _REAL_, intent(in) :: xyz_qm(3),xyz_mm(4),psum(36)
      _REAL_, intent(out) :: pair_force(3)
      integer, intent(in) :: iqm, loop_count, n_atomic_orb
      _REAL_, intent(in) ::  qm_atom_core, alpa

!Local
      _REAL_ FABX, FABY, FABZ ,FNUCX, FNUCY, FNUCZ
      _REAL_ r2, rij, onerij, rr2, rr, one_rija0ev
      _REAL_ vec_qm_mm1, vec_qm_mm2, vec_qm_mm3
      _REAL_ c1, bb, mm_charge
      _REAL_ sqrtaee, dze, qzze, qxxe
      _REAL_ SQRTRRMDDADE, SQRTRRADDADE, SQRTRR2AQE
      _REAL_ SQRTRRAQQAQE, SQRTRRMQQAQE, SQRTRRAQQ2AQE
      _REAL_ Xtdx(3), Xtdy(3), Xtdz(3)
      _REAL_ Ytdx(3), Ytdy(3), Ytdz(3)
      _REAL_ Ztdx(3), Ztdy(3), Ztdz(3)
      _REAL_ TX(3),TY(3),TZ(3), TZ3i, TZ3i2
      _REAL_ RXY2, RYZ2, RZX2, oneRXY
      logical LRXY2, LRYZ2, LRZX2
      _REAL_ TERMX, TERMY, TERMZ
      _REAL_ DGX(4), DGY(4), DGZ(4)
      _REAL_ DRX(MAX_VALENCE_DIMENSION)
      _REAL_ DRY(MAX_VALENCE_DIMENSION)
      _REAL_ DRZ(MAX_VALENCE_DIMENSION)
      _REAL_ EXP1, EXP2, EXP3, EXP4
      _REAL_ DD, QQ
      _REAL_ qm_mm_e_repul(4)
      _REAL_ :: anam1, oner2, temp_real, temp_real2, tmp1, tmp2
      integer:: i, qmitype
      integer isp, m, n, mn
      integer k, kk, l, ll

!****************************************************************************
!*                                                                          *
!*             CALCULATION OF ANALYTICAL DERIVATIVES                        *
!*                                                                          *
!* Routine inlined and optimised by Ross Walker and Mike Crowley (TSRI 2004)*
!****************************************************************************
                                                                               

      vec_qm_mm1=xyz_qm(1) - xyz_mm(1)
      vec_qm_mm2=xyz_qm(2) - xyz_mm(2)
      vec_qm_mm3=xyz_qm(3) - xyz_mm(3)
      mm_charge = xyz_mm(4)
#include "qm2_array_locations.h"
      if (qmmm_nml%qmmmrij_incore) then
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RIJ=qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        one_rija0ev = oneRIJ*A_TO_BOHRS*AU_TO_EV
        EXP1=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)
        EXP2=qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        SQRTAEE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
!Heavy Atom Specific
          SQRTRRADDADE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count)
          SQRTRRMDDADE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count)
          SQRTRR2AQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count)
          SQRTRRAQQAQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count)
          SQRTRRMQQAQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE,loop_count)
          SQRTRRAQQ2AQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count)
!End heavy atom specific
      else
        r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
        RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
        oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
        RIJ=one/oneRIJ
        RR=RIJ*A_TO_BOHRS
        one_rija0ev = A_TO_BOHRS * oneRIJ*AU_TO_EV
        EXP1 = exp(-alpa*rij)
        EXP2 = exp(-ALPH_MM*rij)
        SQRTAEE=one/sqrt(RR2+qm2_params%multip_2c_elec_params(3,iqm)**2)
!Heavy atom specific
          DD=qm2_params%multip_2c_elec_params(1,iqm)
          QQ=qm2_params%multip_2c_elec_params(2,iqm)
          tmp1 = qm2_params%multip_2c_elec_params(4,iqm)**2
          tmp2 = qm2_params%multip_2c_elec_params(5,iqm)**2
          SQRTRRADDADE=one/SQRT((RR+DD)**2+tmp1)
          SQRTRRMDDADE=one/SQRT((RR-DD)**2+tmp1)
          SQRTRR2AQE=one/SQRT(RR2+tmp2)
          SQRTRRAQQAQE=one/SQRT((RR+two*QQ)**2+tmp2)
          SQRTRRMQQAQE=one/SQRT((RR-two*QQ)**2+tmp2)
          SQRTRRAQQ2AQE=one/SQRT(RR2+(four*(QQ**2)+tmp2))
!end heavy_atom specific
      end if

      qm_mm_e_repul(1) = AU_TO_EV*SQRTAEE
      qm_mm_e_repul(2) = HALF_AU_TO_EV*(SQRTRRADDADE - SQRTRRMDDADE)
      qm_mm_e_repul(3) = qm_mm_e_repul(1) + FOURTH_AU_TO_EV*(SQRTRRAQQAQE + SQRTRRMQQAQE)-HALF_AU_TO_EV*SQRTRR2AQE
      qm_mm_e_repul(4) = qm_mm_e_repul(1) + HALF_AU_TO_EV*(SQRTRRAQQ2AQE - SQRTRR2AQE)

!     Returns the
!     derivatives of the electron-core interaction energies in a local
!     diatomic frame.  At most only four terms are computed because there
!     are at most four unique electron-core interactions for the QM-MM
!     pair: (ss|  ), (so|  ), (oo|  ), (pp|  ).  Note that it is not
!     necessary to specify whether an x, y, or z derivative is being
!     evaluated because the formula is the same for all three directions.
!      one_rija0 = one_rija0

      SQRTAEE = -SQRTAEE*SQRTAEE*SQRTAEE*A2_TO_BOHRS2xAU_TO_EV
!        S-orbital of QM atom:
      DGX(1) = vec_qm_mm1*SQRTAEE
      DGY(1) = vec_qm_mm2*SQRTAEE
      DGZ(1) = vec_qm_mm3*SQRTAEE

      DD=qm2_params%multip_2c_elec_params(1,iqm)*one_rija0ev
      SQRTRRADDADE=SQRTRRADDADE*SQRTRRADDADE*SQRTRRADDADE
      SQRTRRMDDADE=SQRTRRMDDADE*SQRTRRMDDADE*SQRTRRMDDADE

      DZE = half*(A2_TO_BOHRS2xAU_TO_EV*(SQRTRRMDDADE-SQRTRRADDADE)-DD*(SQRTRRMDDADE+SQRTRRADDADE))

      SQRTRR2AQE=SQRTRR2AQE*SQRTRR2AQE*SQRTRR2AQE
      SQRTRRAQQ2AQE=SQRTRRAQQ2AQE*SQRTRRAQQ2AQE*SQRTRRAQQ2AQE
      QXXE  = SQRTAEE+half*A2_TO_BOHRS2xAU_TO_EV*(SQRTRR2AQE-SQRTRRAQQ2AQE)

      SQRTRRAQQAQE=SQRTRRAQQAQE*SQRTRRAQQAQE*SQRTRRAQQAQE
      SQRTRRMQQAQE=SQRTRRMQQAQE*SQRTRRMQQAQE*SQRTRRMQQAQE

      QZZE  = SQRTAEE+half*(one_rija0ev*qm2_params%multip_2c_elec_params(2,iqm)*(SQRTRRMQQAQE-SQRTRRAQQAQE) - &
                            half*A2_TO_BOHRS2xAU_TO_EV*(SQRTRRMQQAQE+SQRTRRAQQAQE-two*SQRTRR2AQE))

      DGX(2)=vec_qm_mm1*DZE
      DGX(3)=vec_qm_mm1*QZZE
      DGX(4)=vec_qm_mm1*QXXE

      DGY(2)=vec_qm_mm2*DZE
      DGY(3)=vec_qm_mm2*QZZE
      DGY(4)=vec_qm_mm2*QXXE

      DGZ(2)=vec_qm_mm3*DZE
      DGZ(3)=vec_qm_mm3*QZZE
      DGZ(4)=vec_qm_mm3*QXXE

!     This routine
!     takes the electron-core integral derivatives DG which have been
!     computed for a local frame and rotates them to the molecular frame.
!     The rotated derivatives are returned in DR

!     Determines the transformation (TX,TY,TZ) and derivatives of the transformation
!     (TDX,TDY,TDZ) involved in rotating from a local frame to the molecular frame 
!     for calculation of electron-core interactions.
      RXY2=vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2   
      RYZ2=vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
      RZX2=vec_qm_mm3*vec_qm_mm3+vec_qm_mm1*vec_qm_mm1
      LRXY2 = RXY2 < AXIS_TOL
      LRYZ2 = RYZ2 < AXIS_TOL
      LRZX2 = RZX2 < AXIS_TOL

!     Zeros entire array of 3
      XTDX=zero
      YTDX=zero
      ZTDX=zero
      XTDY=zero
      YTDY=zero
      ZTDY=zero
      XTDZ=zero
      YTDZ=zero
      ZTDZ=zero

      IF(.NOT.(LRXY2 .OR. LRYZ2 .OR. LRZX2)) then 
        oneRXY = one/sqrt(RXY2)
      !Ross Walker + Mike Crowley - rearranged order here slightly for speed.
        TZ(3)=oneRIJ/oneRXY
        TZ3i = one/TZ(3)  !Inverse of TZ(3) to avoid other divisions.
        TZ3i2 = TZ3i*TZ3i  !Square of 1/TZ(3)

        TX(1)=vec_qm_mm1*oneRIJ
        TX(2)=vec_qm_mm2*oneRIJ
        TX(3)=vec_qm_mm3*oneRIJ

        TY(1)=-TX(2)*SIGN(one,TX(1))*TZ3i
        TY(2)=ABS(TX(1)*TZ3i)
        TY(3)=zero

        TZ(1)=-TX(1)*TX(3)*TZ3i
        TZ(2)=-TX(2)*TX(3)*TZ3i

        TERMX = TX(1)*oneRIJ
        TERMY = TX(2)*oneRIJ
        TERMZ = TX(3)*oneRIJ

        XTDX(1)=oneRIJ-TX(1)*TERMX
        YTDX(1)=-TX(1)*TERMY
        ZTDX(1)=-TX(1)*TERMZ

        XTDX(2)=-TX(2)*TERMX
        YTDX(2)=oneRIJ-TX(2)*TERMY
        ZTDX(2)=-TX(2)*TERMZ

        XTDX(3)=-TX(3)*TERMX
        YTDX(3)=-TX(3)*TERMY
        ZTDX(3)=oneRIJ-TX(3)*TERMZ

        XTDZ(3)=TX(1)*oneRXY-TZ(3)*TERMX
        YTDZ(3)=TX(2)*oneRXY-TZ(3)*TERMY
        ZTDZ(3)=-TZ(3)*TERMZ

        XTDY(1)=-XTDX(2)*TZ3i+TX(2)*XTDZ(3)*TZ3i2
        YTDY(1)=-YTDX(2)*TZ3i+TX(2)*YTDZ(3)*TZ3i2
        ZTDY(1)=-ZTDX(2)*TZ3i+TX(2)*ZTDZ(3)*TZ3i2

        XTDY(1)=XTDY(1)*sign(one,TX(1))
        YTDY(1)=YTDY(1)*sign(one,TX(1))
        ZTDY(1)=ZTDY(1)*sign(one,TX(1))

        XTDY(2)=XTDX(1)*TZ3I-TX(1)*XTDZ(3)*TZ3I2
        YTDY(2)=YTDX(1)*TZ3I-TX(1)*YTDZ(3)*TZ3I2
        ZTDY(2)=ZTDX(1)*TZ3I-TX(1)*ZTDZ(3)*TZ3I2

        XTDY(2)=XTDY(2)*sign(one,TX(1))
        YTDY(2)=YTDY(2)*sign(one,TX(1))
        ZTDY(2)=ZTDY(2)*sign(one,TX(1))

! Don't need to zero these again as they were zeroed above
!        XTDY(3)=0.0D0
!        YTDY(3)=0.0D0
!        ZTDY(3)=0.0D0

!Note: Ross Walker and Mike Crowley, we could factor out TZ3I here or we could
!pre-compute -TX(3)*TZ3I etc etc. But for the moment we will leave it as is since
!this is really just doing the compiler's work.
        XTDZ(1)=-TX(3)*XTDX(1)*TZ3I-TX(1)*XTDX(3)*TZ3I &
               +TX(1)*TX(3)*XTDZ(3)*TZ3I2
        YTDZ(1)=-TX(3)*YTDX(1)*TZ3I-TX(1)*YTDX(3)*TZ3I &
               +TX(1)*TX(3)*YTDZ(3)*TZ3I2
        ZTDZ(1)=-TX(3)*ZTDX(1)*TZ3I-TX(1)*ZTDX(3)*TZ3I &
               +TX(1)*TX(3)*ZTDZ(3)*TZ3I2

        XTDZ(2)=-TX(3)*XTDX(2)*TZ3I-TX(2)*XTDX(3)*TZ3I &
               +TX(2)*TX(3)*XTDZ(3)*TZ3I2
        YTDZ(2)=-TX(3)*YTDX(2)*TZ3I-TX(2)*YTDX(3)*TZ3I &
               +TX(2)*TX(3)*YTDZ(3)*TZ3I2
        ZTDZ(2)=-TX(3)*ZTDX(2)*TZ3I-TX(2)*ZTDX(3)*TZ3I &
               +TX(2)*TX(3)*ZTDZ(3)*TZ3I2
      elseif (LRXY2) THEN
!     MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=zero
        TX(2)=zero
        TX(3)=sign(one,vec_qm_mm3)
        TY(1)=zero
        TY(2)=one
        TY(3)=zero
        TZ(1)=TX(3)
        TZ(2)=zero
        TZ(3)=zero
      !X Axis
        XTDX(1)=oneRIJ
        XTDZ(3)=-oneRIJ
      !Y Axis
        YTDX(2)=oneRIJ
        YTDY(3)=-TX(3)*oneRIJ
      !Z Axis
      elseif (LRYZ2) THEN
!     MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=sign(one,vec_qm_mm1)
        TX(2)=zero
        TX(3)=zero
        TY(1)=zero
        TY(2)=TX(1)
        TY(3)=zero
        TZ(1)=zero
        TZ(2)=zero
        TZ(3)=one
      !X Axis
      !Y Axis
        YTDX(2)=oneRIJ
        YTDY(1)=-oneRIJ
   !Z Axis
        ZTDX(3)=oneRIJ
        ZTDZ(1)=-TX(1)*oneRIJ
      else !if (LRZX2) THEN
!     MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=zero
        TX(2)=sign(one,vec_qm_mm2)
        TX(3)=zero
        TY(1)=-TX(2)
        TY(2)=zero
        TY(3)=zero
        TZ(1)=zero
        TZ(2)=zero
        TZ(3)=one
      !X Axis
        XTDX(1)=oneRIJ
        XTDY(2)=oneRIJ
      !Y Axis
      !Z Axis
        ZTDX(3)=oneRIJ
        ZTDZ(2)=-TX(2)*oneRIJ
      end if


      ISP=0
      do K=1,n_atomic_orb
         KK=K-1
         do L=K,n_atomic_orb
            LL=L-1
            ISP=ISP+1
            IF(LL == 0) THEN
               DRX(ISP)=DGX(1)                       !(SS/SS)
               DRY(ISP)=DGY(1)                       !(SS/SS)
               DRZ(ISP)=DGZ(1)                       !(SS/SS)
            ELSEIF(KK == 0) THEN
               DRX(ISP)=DGX(2)*TX(LL)+qm_mm_e_repul(2)*XTDX(LL)   !(SP/SS)
               DRY(ISP)=DGY(2)*TX(LL)+qm_mm_e_repul(2)*YTDX(LL)   !(SP/SS)
               DRZ(ISP)=DGZ(2)*TX(LL)+qm_mm_e_repul(2)*ZTDX(LL)   !(SP/SS)
            ELSE
               DRX(ISP)=DGX(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(XTDX(KK)*TX(LL)+TX(KK)*XTDX(LL))  &
              +DGX(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(XTDY(KK)*TY(LL)+TY(KK)*XTDY(LL)   &
              +XTDZ(KK)*TZ(LL)+TZ(KK)*XTDZ(LL))

               DRY(ISP)=DGY(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(YTDX(KK)*TX(LL)+TX(KK)*YTDX(LL))  &
              +DGY(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(YTDY(KK)*TY(LL)+TY(KK)*YTDY(LL)   &
              +YTDZ(KK)*TZ(LL)+TZ(KK)*YTDZ(LL))

               DRZ(ISP)=DGZ(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(ZTDX(KK)*TX(LL)+TX(KK)*ZTDX(LL))  &
              +DGZ(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(ZTDY(KK)*TY(LL)+TY(KK)*ZTDY(LL)   &
              +ZTDZ(KK)*TZ(LL)+TZ(KK)*ZTDZ(LL))
            ENDIF
         end do
      end do

! BEGIN --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
!           Note: This could technically be done at the same time we do the
!                 energy in hcore_qmmm.

      C1=qm_atom_core*mm_charge

!     ****   START OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***
!     ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!     ANALYTICAL DERIVATIVES

      if ( qmmm_nml%qmmm_int==2 .and. qmmm_struct%AM1_OR_PM3 ) then
        qmitype = qmmm_struct%qm_atom_type(iqm)
        ANAM1=zero
        oner2=oneRIJ*oneRIJ
        do I=1,qm2_params%num_fn(qmitype)
          temp_real=RIJ-qm2_params%FN3(i,qmitype)
          temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
          if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
            ANAM1=ANAM1+qm2_params%FN1(i,qmitype)* &
                 (oner2+two*qm2_params%FN2(i,qmitype)*temp_real*oneRIJ)*EXP(-temp_real2)
          end if
        end do
        ANAM1=-ANAM1*c1*onerij
        FNUCX=ANAM1*vec_qm_mm1
        FNUCY=ANAM1*vec_qm_mm2
        FNUCZ=ANAM1*vec_qm_mm3
     else
        FNUCX=zero; FNUCY=zero; FNUCZ=zero
     endif

!     ****   END OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***

      FNUCX = FNUCX+dgX(1)*c1
      FNUCY = FNUCY+dgy(1)*c1
      FNUCZ = FNUCZ+dgz(1)*c1
      EXP3 = (EXP1+EXP2)*abs(c1)
      EXP4 = qm_mm_e_repul(1)*onerij*(alpa*EXP1 + ALPH_MM*EXP2)*abs(c1)
      FNUCX = FNUCX-vec_qm_mm1*EXP4
      FNUCY = FNUCY-vec_qm_mm2*EXP4
      FNUCZ = FNUCZ-vec_qm_mm3*EXP4

      FNUCX = FNUCX+dgX(1)*EXP3
      FNUCY = FNUCY+dgy(1)*EXP3
      FNUCZ = FNUCZ+dgz(1)*EXP3

! END --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM

      FABX=zero
      FABY=zero
      FABZ=zero
!     MM CORE AFFECTING AO'S ON QM ATOM.
      ISP=0
      do M=1,n_atomic_orb
         BB=one
         do N=M,n_atomic_orb
            MN=M+qm2_params%pascal_tri1(N)
            ISP=ISP+1
            FABX=FABX-BB*mm_charge*PSUM(MN)*DRX(ISP)
            FABY=FABY-BB*mm_charge*PSUM(MN)*DRY(ISP)
            FABZ=FABZ-BB*mm_charge*PSUM(MN)*DRZ(ISP)                    
            BB=two
         end do
      end do

      pair_force(1) = (FABX+FNUCX)*EV_TO_KCAL                                          
      pair_force(2) = (FABY+FNUCY)*EV_TO_KCAL                                           
      pair_force(3) = (FABZ+FNUCZ)*EV_TO_KCAL                                           

      RETURN                                                                    
end subroutine qm2_deriv_qmmm_heavy


