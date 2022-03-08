! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_get_qm_forces(dxyzqm)
!Current code maintained by: Ross Walker (TSRI 2004)

!This routine calculates the derivatives of the energy for QM-QM
!interactions.

! qmmm_struct%qm_coords - QM Atom Coordinates
! dxyzqm  - Returned with the forces in for each QM atom.

      use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, qmmm_mpi, &
                              PM3, AM1, PDDGPM3, PM3CARB1, RM1
      use constants, only : EV_TO_KCAL
      implicit none

!Passed in
      _REAL_, intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)                                  

!Local
      _REAL_ e_repul(22) !Used when qmqm_erep_incore = false
      _REAL_ pair_force(3)
      integer loop_count !Keeps track of number of times through nquant * (nquant-1)/2 loop
      _REAL_ psum(36) !36 = max combinations with heavy and heavy = 4 orbs * 4 orbs (Note, no d orb support)
      _REAL_ xyz_qmi(3), xyz_qmj(3), vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
      integer natqmi, natqmj, qmitype, qmjtype
      integer ii, iif, iil, jj, jjf, jjl, ij
      integer i,j,k,l
      integer n_atomic_orbi, n_atomic_orbj
      integer jstart, jend
      _REAL_ aa,ee,deriv,angle,refh,heat,sum
      _REAL_ corei, corej
      _REAL_ alphai, alphaj, betasas, betasap, betapas, betapap
      _REAL_ bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j
      _REAL_ qqi, qqi2, qqj, qqj2, ddi,ddj
      _REAL_ htype, fqmii(3)

#define CHNGE 1.D-4
#define HALFCHNGE 5.D-5
!one/CHNGE = 10000
#define oneCHNGE 10000
#define DELADJ 1.0D-8
#define TWOONEDELADJ 50000000
   if (qmmm_nml%qmqm_analyt) then !We do analytical derivatives
   !RCW: Note: there is a lot of repeated code in the two options below
   !but this is necessary to factor this if statement out of the inner loop.
      loop_count = 0
#ifdef MPI
      do ii = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
         jstart =  qmmm_mpi%nquant_nlink_jrange(1,ii)
         jend = qmmm_mpi%nquant_nlink_jrange(2,ii)
#else
      do II=2,qmmm_struct%nquant_nlink
         jstart = 1
         jend = ii-1
#endif
         !Loop over all pairs of quantum atoms
!   GET FIRST ATOM INFO                                                             
         iif=qm2_params%orb_loc(1,II)                                                          
         iil=qm2_params%orb_loc(2,II) 
!             n_atomic_orbi = iil - iif + 1
         n_atomic_orbi = qm2_params%natomic_orbs(ii)
         natqmi=qmmm_struct%iqm_atomic_numbers(II)               
         corei=qm2_params%core_chg(ii)
         alphai=qm2_params%cc_exp_params(ii)
         ddi = qm2_params%multip_2c_elec_params(1,ii)
         qqi = qm2_params%multip_2c_elec_params(2,ii)
         qqi2 = qqi*qqi
         bdd1i = qm2_params%multip_2c_elec_params(3,ii)
         bdd2i = qm2_params%multip_2c_elec_params(4,ii)
         bdd3i = qm2_params%multip_2c_elec_params(5,ii)
         xyz_qmi(1:3)=qmmm_struct%qm_coords(1:3,II)        
         qmitype = qmmm_struct%qm_atom_type(ii)
         fqmii(1:3) = 0.0d0
         do JJ=jstart,jend  !jj=1,ii-1
!   GET SECOND ATOM INFO                                 
           jjf=qm2_params%orb_loc(1,JJ)                                                       
           jjl=qm2_params%orb_loc(2,JJ)                                                        
!           n_atomic_orbj = jjl - jjf + 1
           n_atomic_orbj = qm2_params%natomic_orbs(jj)
           natqmj=qmmm_struct%iqm_atomic_numbers(JJ)                                                      
           corej=qm2_params%core_chg(jj)
           alphaj=qm2_params%cc_exp_params(jj)
           ddj = qm2_params%multip_2c_elec_params(1,jj)
           qqj = qm2_params%multip_2c_elec_params(2,jj)
           qqj2 = qqj*qqj
           bdd1j = qm2_params%multip_2c_elec_params(3,jj)
           bdd2j = qm2_params%multip_2c_elec_params(4,jj)
           bdd3j = qm2_params%multip_2c_elec_params(5,jj)
           xyz_qmj(1:3)=qmmm_struct%qm_coords(1:3,JJ)    
           vec_qm_qm1=xyz_qmi(1) - xyz_qmj(1)
           vec_qm_qm2=xyz_qmi(2) - xyz_qmj(2)
           vec_qm_qm3=xyz_qmi(3) - xyz_qmj(3)
           qmjtype = qmmm_struct%qm_atom_type(jj)
           betasas = qm2_params%betasas(qmitype,qmjtype)
           betasap = qm2_params%betasap(qmitype,qmjtype)
           betapas = qm2_params%betasap(qmjtype,qmitype)
           betapap = qm2_params%betapap(qmitype,qmjtype)

           IJ=0                                                       
!  FORM DIATOMIC MATRICES                                                       
           do I=jjf,jjl                                              
             K=qm2_params%pascal_tri1(i)+jjf-1    
             do J=jjf,I                                            
                IJ=IJ+1                                              
                K=K+1                                                
                psum(IJ)=qm2_struct%den_matrix(K)
             end do
           end do

! GET SECOND ATOM FIRST ATOM INTERSECTION                                       
           do I=iif,iil                                              
             L=qm2_params%pascal_tri1(i)
             K=L+jjf-1                                                
             do J=jjf,jjl                                           
                IJ=IJ+1                                              
                K=K+1                                                
                psum(IJ)=qm2_struct%den_matrix(K)                                            
             end do
             K=L+iif-1                                                
             do L=iif,I                                            
                 K=K+1                                                
                 IJ=IJ+1                                              
                 psum(IJ)=qm2_struct%den_matrix(K)                                            
             end do
           end do
           loop_count=loop_count+1
           if (qmmm_nml%qmqm_erep_incore) then
             CALL qm2_deriv_qm_analyt(ii,jj,loop_count,qm2_struct%qm_qm_e_repul(1,loop_count), &
                       psum,natqmi,natqmj,n_atomic_orbj,n_atomic_orbi, &
                       corei,corej,betasas,betasap,betapas,betapap,vec_qm_qm1,vec_qm_qm2,  &
                       vec_qm_qm3,alphai, alphaj,pair_force, qqi, qqi2, qqj, qqj2, ddi, ddj, &
                       bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j, &
                       qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_ss_eqn_adb(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy2(1,1,qmitype,qmjtype), &
                       qmitype,qmjtype)
           else
!Same call as above, just qm2_struct%qm_qm_e_repul(1,loop_count) replaced with local e_repul
             CALL qm2_deriv_qm_analyt(ii,jj,loop_count,e_repul, &
                       psum,natqmi,natqmj,n_atomic_orbj,n_atomic_orbi, &
                       corei,corej,betasas,betasap,betapas,betapap,vec_qm_qm1,vec_qm_qm2,  &
                       vec_qm_qm3,alphai, alphaj,pair_force, qqi, qqi2, qqj, qqj2, ddi, ddj, &
                       bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j, &
                       qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_ss_eqn_adb(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy2(1,1,qmitype,qmjtype), &
                       qmitype,qmjtype)
           end if
           fqmii(1:3) = fqmii(1:3) + pair_force(1:3)
           dxyzqm(1:3,JJ)=dxyzqm(1:3,JJ)+pair_force(1:3)                          
         end do
         dxyzqm(1:3,II)=dxyzqm(1:3,II)-fqmii(1:3)
      end do

!************** END ANALYTICAL DERIVATIVES **************
   else !We will do (pseudo numerical derivatives)
!************** PSEUDO NUMERICAL DERIVATIVES **************
#ifdef MPI
      do ii = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
         jstart =  qmmm_mpi%nquant_nlink_jrange(1,ii)
         jend = qmmm_mpi%nquant_nlink_jrange(2,ii)
#else
         do II=2,qmmm_struct%nquant_nlink
           jstart = 1
           jend = ii-1
#endif
          !Loop over all pairs of quantum atoms
          iif=qm2_params%orb_loc(1,II)
          iil=qm2_params%orb_loc(2,II)
          qmitype = qmmm_struct%qm_atom_type(ii)
          natqmi=qmmm_struct%iqm_atomic_numbers(II)
          xyz_qmi(1)=qmmm_struct%qm_coords(1,II)
          xyz_qmi(2)=qmmm_struct%qm_coords(2,II)
          xyz_qmi(3)=qmmm_struct%qm_coords(3,II)
          do JJ=jstart,jend !jj=1,ii-1
!  FORM DIATOMIC MATRICES
            jjf=qm2_params%orb_loc(1,JJ)
            jjl=qm2_params%orb_loc(2,JJ)
!   GET FIRST ATOM
            qmjtype = qmmm_struct%qm_atom_type(jj)
            natqmj=qmmm_struct%iqm_atomic_numbers(JJ)
            xyz_qmj(1)=qmmm_struct%qm_coords(1,JJ)
            xyz_qmj(2)=qmmm_struct%qm_coords(2,JJ)
            xyz_qmj(3)=qmmm_struct%qm_coords(3,JJ)
            IJ=0
            do I=jjf,jjl
              K=qm2_params%pascal_tri1(i)+jjf-1
              do J=jjf,I
                IJ=IJ+1
                K=K+1
                psum(IJ)=qm2_struct%den_matrix(K)
              end do
            end do
! GET SECOND ATOM FIRST ATOM INTERSECTION
            do I=iif,iil
               L=qm2_params%pascal_tri1(i)
               K=L+jjf-1
               do J=jjf,jjl
                  IJ=IJ+1
                  K=K+1
                  psum(IJ)=qm2_struct%den_matrix(K)
               end do
               K=L+iif-1
               do L=iif,I
                   K=K+1
                   IJ=IJ+1
                   psum(IJ)=qm2_struct%den_matrix(K)
               end do
            end do
            do K=1,3
              xyz_qmi(K)=xyz_qmi(K)-HALFCHNGE
              CALL qm2_dhc(psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,jjf, &
                       jjl,iif,iil,AA)
              xyz_qmi(K)=xyz_qmi(K)+CHNGE
              CALL qm2_dhc(psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,jjf, &
                       jjl,iif,iil,EE)
              xyz_qmi(K)=xyz_qmi(K)-HALFCHNGE
              DERIV=(AA-EE)*EV_TO_KCAL*oneCHNGE
              dxyzqm(K,II)=dxyzqm(K,II)-DERIV
              dxyzqm(K,JJ)=dxyzqm(K,JJ)+DERIV
            end do
          end do
      end do
!************** END PSEUDO NUMERICAL DERIVATIVES **************
   end if

   if(qmmm_nml%peptide_corr) then
!  NOW ADD IN MOLECULAR-MECHANICS CORRECTION TO THE H-N-C=O TORSION            
     if (qmmm_nml%qmtheory==PM3 .OR. qmmm_nml%qmtheory==PDDGPM3 .OR. qmmm_nml%qmtheory == PM3CARB1 ) then
       htype = 7.1853D0                                                      
     elseif (qmmm_nml%qmtheory==AM1 .OR. qmmm_nml%qmtheory==RM1) then
       htype = 3.3191D0                                                      
     else !Assume MNDO
       htype = 6.1737D0                                                      
     end if
!Parallel
     do I=qmmm_mpi%mytaskid+1,qm2_struct%n_peptide_links,qmmm_mpi%numthreads !1,n_peptide_links
       do J=1,4                                    
         do K=1,3                                
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))-DELADJ
           CALL qm2_dihed(qmmm_struct%qm_coords,qm2_struct%peptide_links(1,I), &
                          qm2_struct%peptide_links(2,I),qm2_struct%peptide_links(3,I), &
                          qm2_struct%peptide_links(4,I),ANGLE)
           REFH=HTYPE*SIN(ANGLE)**2         
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))+DELADJ*2.D0
           CALL qm2_dihed(qmmm_struct%qm_coords,qm2_struct%peptide_links(1,I),qm2_struct%peptide_links(2,I), &
                          qm2_struct%peptide_links(3,I),qm2_struct%peptide_links(4,I),ANGLE)
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))-DELADJ
           HEAT=HTYPE*SIN(ANGLE)**2         
           SUM=(REFH-HEAT)*TWOONEDELADJ
           dxyzqm(K,qm2_struct%peptide_links(J,I))=dxyzqm(K,qm2_struct%peptide_links(J,I))-SUM 
         end do
       end do                                   
     end do                                    
   end if                                           
   RETURN                                         
end subroutine qm2_get_qm_forces

subroutine qm2_deriv_qm_analyt(iqm,jqm,loop_count,qm_qm_e_repul,PSUM,natqmi, &
           natqmj,n_atomic_orbj,n_atomic_orbi,corei,corej,betasas,betasap,betapas,betapap, &
           vec_qm_qm1,vec_qm_qm2, vec_qm_qm3,alphai, alphaj, pair_force, &
           qqi, qqi2, qqj, qqj2, ddi, ddj, bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, &
           bdd3j,zz_sxs_over_sas,ss_eqn_adb,zz_sxp_over_sapij, &
           zz_sxp_over_sapji,sp_eqn_xx1ij,sp_eqn_xx2ij,sp_eqn_xx1ji,sp_eqn_xx2ji, &
           sp_eqn_xyij, sp_eqn_xyji,zz_pxp_over_pap,pp_eqn_xxy1,pp_eqn_xxy2, &
           qmitype,qmjtype)

!************************************************************************        
!*                                                                      *        
!*         CALCULATION OF ANALYTICAL DERIVATIVES                        *      
!*                                                                      *
!* Current routine maintained by Ross Walker (TSRI, 2004)               *
!* Inlining and optimisation by Ross Walker (TSRI, 2004)                *  
!*                                                                      *        
!************************************************************************        


!*** Variable Definitions:
!
!    qm_qm_e_repul = QM-QM electron repulsion integrals for this QM-QM pair
!             psum = Density matrix elements for orbitals centered on this
!                    QM-QM pair. - For coulomb term
!           natqmi = atomic number of QM atom i
!           natqmj = atomic number of QM atom j
!    n_atomic_orbj = number of atomic orbitals on atom j
!    n_atomic_orbi = number of atomic orbitals on atom i
!       vec_qm_qmx = xyz_qmi(x) - xyz_qmj(x) 
!           alphai =
!           alphaj =
!       pair_force = Returned with cartesian force on this QM-QM pair

   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_rij_eqns, qm2_struct, qm2_params, &
                           AM1, PM3, MNDO, PDDGPM3, PDDGMNDO, PM3CARB1, RM1, &
                           AXIS_TOL, OVERLAP_CUTOFF, EXPONENTIAL_CUTOFF

   use constants, only : A_TO_BOHRS, A3_TO_BOHRS3, A2_TO_BOHRS2, AU_TO_EV, zero, &
                         half, one, two, four, fourth, eighth, sixteenth, thirtysecond, &
                         EV_TO_KCAL

   implicit none

!Passed in
   integer, intent(in) :: iqm,jqm,loop_count,qmitype,qmjtype
   _REAL_, intent(inout) :: qm_qm_e_repul(22) !for qmqm_erep_incore=true it is in, for =false it is out.
   _REAL_, intent(in) :: psum(36)
   integer, intent(in) :: natqmi, natqmj
   integer, intent(in) :: n_atomic_orbj, n_atomic_orbi
   _REAL_, intent(in) :: corei,corej,betasas,betasap,betapas,betapap
   _REAL_, intent(in) :: vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
   _REAL_, intent(in) :: alphai, alphaj, qqi, qqi2, qqj, qqj2, ddi, ddj
   _REAL_, intent(in) :: bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j
   _REAL_, intent(out):: pair_force(3)
   _REAL_, intent(in) :: zz_sxs_over_sas(6,6),ss_eqn_adb(6,6)
   _REAL_, intent(in) :: zz_sxp_over_sapij(6,6), zz_sxp_over_sapji(6,6)
   _REAL_, intent(in) :: sp_eqn_xx1ij(6,6), sp_eqn_xx2ij(6,6)
   _REAL_, intent(in) :: sp_eqn_xx1ji(6,6), sp_eqn_xx2ji(6,6)
   _REAL_, intent(in) :: sp_eqn_xyij(6,6), sp_eqn_xyji(6,6)
   _REAL_, intent(in) :: zz_pxp_over_pap(6,6)
   _REAL_, intent(in) :: pp_eqn_xxy1(6,6), pp_eqn_xxy2(6,6)

!Local
   _REAL_ FAAX, FAAY, FAAZ,FABX, FABY, FABZ,FNUCX, FNUCY, FNUCZ
   _REAL_ dsx,dsy,dsz,DSPX(3), DSPY(3), DSPZ3
   _REAL_ DPXPX(3), DPYPY(3), DPZPZ(3), DPCROSS
   _REAL_ dgx(22), dgy(22), dgz(22)
   _REAL_ drx(100), dry(100), drz(100)
   _REAL_ r2, oner2, rij, onerij, rr, rr2
   _REAL_ c1, f3, ddx, ddy, ddz
   _REAL_ part1x, part1y, part1z
   _REAL_ part2x, part2y, part2z
   _REAL_ part3x, part3y, part3z
   _REAL_ anam1
   _REAL_ bb, aa
   integer total_atomic_orb, qm_atomi_orb_start, orb_offset
   integer isp, k, l
   integer m, n, mn, kk, ll, kl, mk, nk, ml, nl

   _REAL_ vec_qm_qm_onerij1, vec_qm_qm_onerij2, vec_qm_qm_onerij3
   _REAL_ vec2_qm_qm1, vec2_qm_qm2, vec2_qm_qm3
   _REAL_ vec_qm_qm123, vec_qm_qm12, vec_qm_qm23, vec_qm_qm13
   integer i, j
   _REAL_ ABN, ADBR2, SS
   _REAL_ temp_real, temp_real2, temp_real3, temp_real4
   _REAL_ EXP1i, EXP1j, SQRTAEE, bdd1ij

!Specific for PDDG
   _REAL_ :: PDDG_EXP1, PDDG_EXP2, PDDG_EXP3, PDDG_EXP4, PDDG_CORR

!Originally in delri
   _REAL_ termx, termy, termz, ee, ade, aqe, dze, qzze, qxxe
   _REAL_ aed, aeq, edz, eqzz, eqxx
   _REAL_ add, adq, aqd, aqq, dxdx, dzdz,dzqxx,qxxdz,dzqzz,qzzdz
   _REAL_ qxxqxx, qxxqyy, qxxqzz, qzzqxx, qzzqzz, dxqxz, qxzdx, qxzqxz

!Originally in delmol
   integer mm, nn
   _REAL_ xTDX(3),xTDY(3),xTDZ(3)
   _REAL_ yTDX(3),yTDY(3),yTDZ(3)
   _REAL_ zTDX(3),zTDY(3),zTDZ(3)
   _REAL_ TX(3),TY(3),TZ(3), TZ3i, TZ3i2
   _REAL_ xtemp1, ytemp1, ztemp1, xtemp2, ytemp2, ztemp2

!Originally for rotat
   _REAL_ RXY2, RYZ2, RZX2, oneRXY
   logical LRXY2, LRYZ2, LRZX2

   FAAX=zero; FAAY=zero; FAAZ=zero
   FABX=zero; FABY=zero; FABZ=zero

   if (n_atomic_orbi > 1 .OR. n_atomic_orbj >1) then      
     vec_qm_qm123 = vec_qm_qm1*vec_qm_qm2*vec_qm_qm3
     vec_qm_qm12 = vec_qm_qm1*vec_qm_qm2
     vec_qm_qm23 = vec_qm_qm2*vec_qm_qm3
     vec_qm_qm13 = vec_qm_qm1*vec_qm_qm3
   end if

   vec2_qm_qm1 = vec_qm_qm1*vec_qm_qm1
   vec2_qm_qm2 = vec_qm_qm2*vec_qm_qm2
   vec2_qm_qm3 = vec_qm_qm3*vec_qm_qm3

#include "qm2_array_locations.h"   
   if (qmmm_nml%qmqmrij_incore) then
        rr2=qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS2,loop_count)
     oneRIJ=qm2_rij_eqns%qmqmrijdata(QMQMONERIJ,loop_count) 
        RIJ=qm2_rij_eqns%qmqmrijdata(QMQMRIJ,loop_count)
         RR=qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS,loop_count)
      EXP1i=qm2_rij_eqns%qmqmrijdata(QMQMEXP1I,loop_count)
      EXP1j=qm2_rij_eqns%qmqmrijdata(QMQMEXP1J,loop_count)
      SQRTAEE=qm2_rij_eqns%qmqmrijdata(QMQMSQRTAEE,loop_count)
      if (qmmm_struct%PDDG_IN_USE) then
!PDDG Specific terms
        PDDG_EXP1 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP1,loop_count)
        PDDG_EXP2 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP2,loop_count)
        PDDG_EXP3 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP3,loop_count)
        PDDG_EXP4 = qm2_rij_eqns%qmqmrijdata(QMQMPDDG_EXP4,loop_count)
      end if
   else
     r2=vec2_qm_qm1+vec2_qm_qm2+vec2_qm_qm3                                          
     rr2=r2 * A2_TO_BOHRS2
     oneRIJ = one/sqrt(r2) !1/sqrt is faster than doing sqrt.
     RIJ = r2*oneRIJ !one/oneRIJ
     RR=RIJ * A_TO_BOHRS                                                                
     EXP1i = EXP(-alphai*RIJ)
     EXP1j = EXP(-alphaj*RIJ)
     bdd1ij = bdd1i+bdd1j
     bdd1ij = bdd1ij*bdd1ij
     SQRTAEE=1.0d0/sqrt(RR2+bdd1ij)
     if (qmmm_struct%PDDG_IN_USE) then
!PDDG Specific terms
       PDDG_EXP1 = EXP(-10.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge1(jqm))**2)
       PDDG_EXP2 = EXP(-10.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge2(jqm))**2)
       PDDG_EXP3 = EXP(-10.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge1(jqm))**2)
       PDDG_EXP4 = EXP(-10.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge2(jqm))**2)
     end if
   end if

   vec_qm_qm_onerij1 = vec_qm_qm1 * oneRIJ
   vec_qm_qm_onerij2 = vec_qm_qm2 * oneRIJ
   vec_qm_qm_onerij3 = vec_qm_qm3 * oneRIJ

   total_atomic_orb=n_atomic_orbi+n_atomic_orbj

   qm_atomi_orb_start=n_atomic_orbj+1

!If we don't have the 1-e repul integrals for this QM-QM pair in memory we need
!to calculate them now.
   if (.NOT. qmmm_nml%qmqm_erep_incore) then
     call qm2_repp(iqm,jqm,rr,rr2,qm_qm_e_repul,SQRTAEE)
   end if

   C1=corei*corej                                                

!   ****   START OF THE AM1 and PM3 RM1 etc SPECIFIC DERIVATIVE CODE   ***
!   ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!    ANALYTICAL DERIVATIVES
   IF( qmmm_struct%AM1_OR_PM3 )THEN
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
      do i=1,qm2_params%num_fn(qmjtype)
        temp_real=RIJ-qm2_params%FN3(i,qmjtype)
        temp_real2=qm2_params%FN2(i,qmjtype)*temp_real*temp_real
        if (temp_real2 < EXPONENTIAL_CUTOFF) then !Skip doing the exponential if it is essentially zero
          ANAM1=ANAM1+qm2_params%FN1(i,qmjtype)* &
                      (oner2+two*qm2_params%FN2(i,qmjtype)*temp_real*oneRIJ)*EXP(-temp_real2)
        end if
      end do
      ANAM1=-ANAM1*c1
      FNUCX=ANAM1*vec_qm_qm_oneRIJ1
      FNUCY=ANAM1*vec_qm_qm_oneRIJ2
      FNUCZ=ANAM1*vec_qm_qm_oneRIJ3
   else
      FNUCX=zero; FNUCY=zero; FNUCZ=zero
   endif
!   ****   END OF THE AM1 and PM3 RM1 etc SPECIFIC DERIVATIVE CODE   ***

!   ****   PDDG SPECIFIC TERMS ****
   !PDDG Specific terms
   if (qmmm_struct%PDDG_IN_USE) then
         PDDG_EXP1 = -20.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge1(jqm)) * PDDG_EXP1
         PDDG_EXP2 = -20.0D0 * (RIJ - qm2_params%pddge1(iqm) - qm2_params%pddge2(jqm)) * PDDG_EXP2
         PDDG_EXP3 = -20.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge1(jqm)) * PDDG_EXP3
         PDDG_EXP4 = -20.0D0 * (RIJ - qm2_params%pddge2(iqm) - qm2_params%pddge2(jqm)) * PDDG_EXP4
         PDDG_CORR = qm2_params%PDDG_TERM1(qmitype,qmjtype)*PDDG_EXP1 + &
                     qm2_params%PDDG_TERM2(qmitype,qmjtype)*PDDG_EXP2 + &
                     qm2_params%PDDG_TERM3(qmitype,qmjtype)*PDDG_EXP3 + &
                     qm2_params%PDDG_TERM4(qmitype,qmjtype)*PDDG_EXP4

         FNUCX = FNUCX+PDDG_CORR*vec_qm_qm_oneRIJ1
         FNUCY = FNUCY+PDDG_CORR*vec_qm_qm_oneRIJ2
         FNUCZ = FNUCZ+PDDG_CORR*vec_qm_qm_oneRIJ3
    end if
!   ****   END PDDG SPECIFIC TERMS ****

!   THE FIRST DERIVATIVES OF OVERLAP INTEGRALS                                  
      !Loop over atomic orbitals for QM atom i
      !for PX, PY, PZ - So K=0 = S orbital, K=1 = PX, 2 = PY, 3 = PZ
      ! CALCULATE OVERLAP DERIVATIVES, STORE RESULTS IN DS                     

   if (RR2 < OVERLAP_CUTOFF) then  !10A cutoff on overlaps
!ALL Atom pairs must have at least an S-S interaction - Min 1 orbital per atom.
! (S/S) Terms
     dsx=zero ; dsy=zero; dsz = zero
     do I=1,6 !1 to NGAUSS
       do J=1,6
         ADBR2=zz_sxs_over_sas(i,j)*RR2
!Check overlap is non-zero before starting
         if(ADBR2 < EXPONENTIAL_CUTOFF) then
           SS=ss_eqn_adb(i,j)*EXP(-ADBR2)
           DSx=DSx+vec_qm_qm1*SS
           DSy=DSy+vec_qm_qm2*SS
           DSz=DSz+vec_qm_qm3*SS
         end if
       end do
     end do

!    COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
     orb_offset=1+qm2_params%pascal_tri1(qm_atomi_orb_start)
     temp_real=betasas*PSUM(orb_offset)
     FABx=FABx+temp_real*DSx
     FABy=FABy+temp_real*DSy
     FABz=FABz+temp_real*DSz

! END OF S/S terms, if both QM atoms have only 1 orbital then we don't do the loops below.

!NOW DO S with P orbitals if necessary (K=0)
     if (n_atomic_orbj>1) then
       DSPX=zero ; DSPY=zero !Zeros entire arrays of 3
       DSPZ3 = zero
       do j=1,6
         do i=1,6
           ADBR2=zz_sxp_over_sapij(i,j)*RR2
           !Check overlap is non-zero - otherwise we can skip this bit.
           if (ADBR2<EXPONENTIAL_CUTOFF) then
             ADBR2 = EXP(-ADBR2)
!   (S/PX-x) TERM 
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm1)
             DSPX(1)=DSPX(1)+ABN*ADBR2
!   (S/PY-y) TERM
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm2)
             DSPY(2)=DSPY(2)+ABN*ADBR2
!   (S/PZ-z) TERM
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm3)
             DSPZ3=DSPZ3+ABN*ADBR2

             ABN=ADBR2*sp_eqn_xyij(i,j)
!   (PX-y/S) TERM
             DSPX(2)=DSPX(2)+ABN*vec_qm_qm12
!   (PX-z/S) TERM
             DSPX(3)=DSPX(3)+ABN*vec_qm_qm13
!   (PY-z/S) TERM
             DSPY(3)=DSPY(3)+ABN*vec_qm_qm23
!   (PY-x/S) TERM
             !DSPY(1)=DSPY(1)-ABN*vec_qm_qm12 DSPY(1)=DSPX(2)
!   (PZ-x/S) TERM
             !DSPZ(1)=DSPZ(1)-ABN*vec_qm_qm13 DSPZ(1)=DSPX(3)
!   (PZ-y/S) TERM
             !DSPZ(2)=DSPZ(2)-ABN*vec_qm_qm23 DSPZ(2)=DSPY(3)
           end if !(ADBR2<EXPONENTIAL_CUTOFF)
         end do
       end do
       !   COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
       temp_real = betasap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DSPX(1)
       FABy=FABy+temp_real*DSPX(2)
       FABz=FABz+temp_real*DSPX(3)
       temp_real = betasap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DSPX(2) !DSPY(1)=DSPX(2)
       FABy=FABy+temp_real*DSPY(2)
       FABz=FABz+temp_real*DSPY(3)
       temp_real = betasap*PSUM(orb_offset+3)
       FABx=FABx+temp_real*DSPX(3) !DSPZ(1)=DSPX(3)
       FABy=FABy+temp_real*DSPY(3) !DSPZ(2)=DSPY(3)
       FABz=FABz+temp_real*DSPZ3
     end if ! if (n_atomic_orbj>1)

!Now do P orbitals of atom 1 with S of atom 2 (K>0 and L=0)
     if (n_atomic_orbi>1) then
       DSPX=zero ; DSPY=zero !Zeros entire arrays of 3
       DSPZ3=zero
       do I=1,6
         do J=1,6
           ADBR2=zz_sxp_over_sapji(j,i)*RR2
           if (ADBR2<EXPONENTIAL_CUTOFF) then !Only do the rest if the overlap is not zero
             ADBR2 = EXP(-ADBR2)
!   (PX-x/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm1)
             DSPX(1)=DSPX(1)+ABN*ADBR2
!   (PY-y/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm2)
             DSPY(2)=DSPY(2)+ABN*ADBR2
!   (PZ-z/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm3)
             DSPZ3=DSPZ3+ABN*ADBR2

             ABN=ADBR2*sp_eqn_xyji(j,i)
!   (PX-y/S) TERM
             DSPX(2)=DSPX(2)-ABN*vec_qm_qm12
!   (PX-z/S) TERM
             DSPX(3)=DSPX(3)-ABN*vec_qm_qm13
!   (PY-z/S) TERM
             DSPY(3)=DSPY(3)-ABN*vec_qm_qm23
!   (PY-x/S) TERM
             !DSPY(1)=DSPY(1)-ABN*vec_qm_qm12 DSPY(1)=DSPX(2)
!   (PZ-x/S) TERM
             !DSPZ(1)=DSPZ(1)-ABN*vec_qm_qm13 DSPZ(1)=DSPX(3)
!   (PZ-y/S) TERM
             !DSPZ(2)=DSPZ(2)-ABN*vec_qm_qm23 DSPZ(2)=DSPY(3)
           end if !(ADBR2<EXPONENTIAL_CUTOFF)
         end do
       end do
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+1))
       FABx=FABx+temp_real*DSPX(1)
       FABy=FABy+temp_real*DSPX(2)
       FABz=FABz+temp_real*DSPX(3)
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+2))
       FABx=FABx+temp_real*DSPX(2) !DSPY(1)=DSPX(2)
       FABy=FABy+temp_real*DSPY(2)
       FABz=FABz+temp_real*DSPY(3)
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+3))
       FABx=FABx+temp_real*DSPX(3) !DSPZ(1)=DSPX(3)
       FABy=FABy+temp_real*DSPY(3) !DSPZ(2)=DSPY(3)
       FABz=FABz+temp_real*DSPZ3
     end if !if (n_atomic_orbi>1)

!Ross Walker - PP combinations that are the same have been condensed
!to one variable for speed and to save memory.

!Finally do P's of atom 1 with P's of atom 2 (if necessary)
     if (n_atomic_orbi > 1 .and. n_atomic_orbj > 1) then
       DPXPX=zero; DPYPY=zero; DPZPZ=zero !Zeros entire arrays of 3
       DPCROSS=zero
       do I=1,6
          do J=1,6
            ADBR2=zz_pxp_over_pap(i,j)*RR2
            if (ADBR2<EXPONENTIAL_CUTOFF) then !Only do the next bit if the overlap is not zero
              ADBR2=EXP(-ADBR2)
!    (PX / PX) - x
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm1)*ADBR2
              DPXPX(1)=DPXPX(1)+ABN*vec_qm_qm1
!    (PY / PY) - y
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm2)*ADBR2
              DPYPY(2)=DPYPY(2)+ABN*vec_qm_qm2
!    (PZ / PZ) - z
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm3)*ADBR2
              DPZPZ(3)=DPZPZ(3)+ABN*vec_qm_qm3

              ABN = (pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm1)*ADBR2
!    (PX / PX) - y
              DPXPX(2)=DPXPX(2)+ABN*vec_qm_qm2
!    (PX / PX) - z
              DPXPX(3)=DPXPX(3)+ABN*vec_qm_qm3
!    (PX / PY) - x
              !DPXPY(1)=DPXPY(1)+ABN*vec_qm_qm2 PXPY(1)=PXPX(2)
!    (PX / PZ) - x
              !DPXPZ(1)=DPXPZ(1)+ABN*vec_qm_qm3 PXPZ(1)=PXPX(3)
!    (PY / PX) - x
              !DPYPX(1)=DPYPX(1)+ABN*vec_qm_qm2 PXPY=PYPX
!    (PZ / PX) - x
              !DPZPX(1)=DPZPX(1)+ABN*vec_qm_qm3 PXPZ=PZPX

              ABN=(pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm2)*ADBR2
!    (PY / PY) - x
              DPYPY(1)=DPYPY(1)+ABN*vec_qm_qm1
!    (PY / PY) - z
              DPYPY(3)=DPYPY(3)+ABN*vec_qm_qm3
!    (PX / PY) - y
!              DPXPY(2)=DPXPY(2)+ABN*vec_qm_qm1 PXPY(2)=PYPX(2)=PYPY(1)
!    (PY / PZ) - y
!              DPYPZ(2)=DPYPZ(2)+ABN*vec_qm_qm3 PYPZ(2)=PZPY(2)=PYPY(3)
!    (PY / PX) - y
              !DPYPX(2)=DPYPX(2)+ABN*vec_qm_qm1 PXPY=PYPX
!    (PZ / PY) - y
              !DPZPY(2)=DPZPY(2)+ABN*vec_qm_qm3 PYPZ=PZPY

              ABN=(pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm3)*ADBR2
!    (PZ / PZ) - x
              DPZPZ(1)=DPZPZ(1)+ABN*vec_qm_qm1
!    (PZ / PZ) - y
              DPZPZ(2)=DPZPZ(2)+ABN*vec_qm_qm2
!    (PX / PZ) - z
!              DPXPZ(3)=DPXPZ(3)+ABN*vec_qm_qm1 PXPZ(3)=PZPX(3)=PZPZ(1)
!    (PY / PZ) - z
!              DPYPZ(3)=DPYPZ(3)+ABN*vec_qm_qm2 PYPZ(3)=PZPY(3)=PZPZ(2)
!    (PZ / PY) - z
              !DPZPY(3)=DPZPY(3)+ABN*vec_qm_qm2 PYPZ=PZPY
!    (PZ / PX) - z
              !DPZPX(3)=DPZPX(3)+ABN*vec_qm_qm1 PXPZ=PZPX

              ABN=pp_eqn_xxy2(i,j) * ADBR2*vec_qm_qm123
              DPCROSS=DPCROSS+ABN
!    (PX / PY) - z
              !DPXPY(3)=DPXPY(3)+ABN !PXPY(3)=PYPX(3)=DPCROSS
!    (PX / PZ) - y
              !DPXPZ(2)=DPXPZ(2)+ABN !PXPZ(2)=PZPX(2)=DPCROSS
!    (PY / PZ) - x
              !DPYPZ(1)=DPYPZ(1)+ABN !PYPZ(1)=PZPY(1)=DPCROSS
!    (PZ / PY) - x
              !DPZPY(1)=DPZPY(1)+ABN PYPZ=PZPY
!    (PY / PX) - z
              !DPYPX(3)=DPYPX(3)+ABN PXPY=PYPX
!    (PZ / PX) - y
              !DPZPX(2)=DPZPX(2)+ABN PXPZ=PZPX
           end if !ADBR2>EXPONENTIAL_CUTOFF
         end do
       end do
!      COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+1)
       temp_real = betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(1)
       FABy=FABy+temp_real*DPXPX(2)
       FABz=FABz+temp_real*DPXPX(3)
       temp_real = betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPXPX(2) !PYPX(1)=PXPY(1)=PXPX(2)
       FABy=FABy+temp_real*DPYPY(1) !PXPY(2)=PYPX(2)=PYPY(1)
       FABz=FABz+temp_real*DPCROSS  !PXPY(3)=PYPX(3)=DPCROSS
       temp_real = betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPXPX(3) !PXPZ(1)=PZPX(1)=PXPX(3)
       FABy=FABy+temp_real*DPCROSS  !PXPZ(2)=PZPX(2)=DPCROSS
       FABz=FABz+temp_real*DPZPZ(1) !PXPZ(3)=PZPX(3)=PZPZ(1)
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+2)
       temp_real=betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(2) !PYPX(1)=PXPY(1)=PXPX(2)
       FABy=FABy+temp_real*DPYPY(1) !PXPY(2)=PYPX(2)=PYPY(1)
       FABz=FABz+temp_real*DPCROSS  !PXPY(3)=PYPX(3)=DPCROSS 
       temp_real=betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPYPY(1)
       FABy=FABy+temp_real*DPYPY(2)
       FABz=FABz+temp_real*DPYPY(3)
       temp_real=betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPCROSS  !PYPZ(1)=PZPY(1)=DPCROSS
       FABy=FABy+temp_real*DPYPY(3) !PYPZ(2)=PZPY(2)=PYPY(3)
       FABz=FABz+temp_real*DPZPZ(2) !PYPZ(3)=PZPY(3)=PZPZ(2)
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+3)
       temp_real = betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(3) !PXPZ(1)=PZPX(1)=PXPX(3)
       FABy=FABy+temp_real*DPCROSS  !PXPZ(2)=PZPX(2)=DPCROSS
       FABz=FABz+temp_real*DPZPZ(1) !PXPZ(3)=PZPX(3)=PZPZ(1)
       temp_real = betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPCROSS  !PYPZ(1)=PZPY(1)=DPCROSS
       FABy=FABy+temp_real*DPYPY(3) !PYPZ(2)=PZPY(2)=PYPY(3)
       FABz=FABz+temp_real*DPZPZ(2) !PYPZ(3)=PZPY(3)=PZPZ(2)
       temp_real = betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPZPZ(1)
       FABy=FABy+temp_real*DPZPZ(2)
       FABz=FABz+temp_real*DPZPZ(3)
     end if !n_atomic_orbi >1 .and. n_atomic_orbj > 1
   end if !if (RR2 < 100.D0 * A2_TO_BOHRS2) then  !10A cutoff on overlaps
!Code is linear from this point on for a given pair.

!NOTE: Ross Walker - In the equations below the only thing that varies for this pair
!during a simulation is the value of RR so we may want to consider pre-computing
!a lot of this. E.g. For a given pair AEE is constant.

!  S-orbital of QM-QM pairs - all QM pairs have s-s interactions.
   TERMX=vec_qm_qm_onerij1*AU_TO_EV*A_TO_BOHRS
   TERMY=vec_qm_qm_onerij2*AU_TO_EV*A_TO_BOHRS
   TERMZ=vec_qm_qm_onerij3*AU_TO_EV*A_TO_BOHRS
   EE    =-RR*(SQRTAEE)**3
   DGX(1)=TERMX*EE
   DGY(1)=TERMY*EE
   DGZ(1)=TERMZ*EE
   if(n_atomic_orbi > 2) then
!   HEAVY ATOM-HYDROGEN
      ADE=(bdd2i+bdd1j)**2

      AQE=(bdd3i+bdd1j)**2

      DZE   = (RR+ddi)/(SQRT((RR+ddi)**2+ADE))**3 &
              -(RR-ddi)/(SQRT((RR-ddi)**2+ADE))**3

      temp_real = (two*RR)/(SQRT(RR2+AQE))**3

      QZZE  =-(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQE))**3 &
             -(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQE))**3 &
             +temp_real

      QXXE  =-(two*RR)/(SQRT(RR2+four*qqi2+AQE))**3 &
             +temp_real

      temp_real = -DZE*half

      DGX(2)=TERMX*temp_real
      DGY(2)=TERMY*temp_real
      DGZ(2)=TERMZ*temp_real
      temp_real = EE+QZZE*fourth
      DGX(3)=TERMX*temp_real
      DGY(3)=TERMY*temp_real
      DGZ(3)=TERMZ*temp_real
      temp_real = EE+QXXE*fourth
      DGX(4)=TERMX*temp_real
      DGY(4)=TERMY*temp_real
      DGZ(4)=TERMZ*temp_real
   end if
   if (n_atomic_orbj > 2) then
!   HYDROGEN-HEAVY ATOM
      AED=(bdd1i+bdd2j)**2

      AEQ=(bdd1i+bdd3j)**2

      EDZ   = (RR-ddj)/(SQRT((RR-ddj)**2+AED))**3 &
             -(RR+ddj)/(SQRT((RR+ddj)**2+AED))**3

      temp_real = (two*RR)/(SQRT(RR2+AEQ))**3

      EQZZ  =-(RR-two*qqj)/(SQRT((RR-two*qqj)**2+AEQ))**3 &
             -(RR+two*qqj)/(SQRT((RR+two*qqj)**2+AEQ))**3 &
             +temp_real

      EQXX  =-(two*RR)/(SQRT(RR2+four*qqj2+AEQ))**3 &
             +temp_real

      temp_real = -EDZ*half
      DGX(5)=TERMX*temp_real
      DGY(5)=TERMY*temp_real
      DGZ(5)=TERMZ*temp_real
      temp_real = EE+EQZZ*fourth
      DGX(11)=TERMX*temp_real
      DGY(11)=TERMY*temp_real
      DGZ(11)=TERMZ*temp_real
      temp_real = EE+EQXX*fourth
      DGX(12)=TERMX*temp_real
      DGY(12)=TERMY*temp_real
      DGZ(12)=TERMZ*temp_real
   end if
   if (n_atomic_orbi > 2 .and. n_atomic_orbj > 2) then
!   HEAVY ATOM-HEAVY ATOM
      ADD=(bdd2i+bdd2j)**2
      ADQ=(bdd2i+bdd3j)**2
      AQD=(bdd3i+bdd2j)**2
      AQQ=(bdd3i+bdd3j)**2

      DXDX  =-(two*RR)/(SQRT(RR2+(ddi-ddj)**2+ADD))**3 &
             +(two*RR)/(SQRT(RR2+(ddi+ddj)**2+ADD))**3

      DZDZ  =-(RR+ddi-ddj)/(SQRT((RR+ddi-ddj)**2+ADD))**3 &
             -(RR-ddi+ddj)/(SQRT((RR-ddi+ddj)**2+ADD))**3 &
             +(RR-ddi-ddj)/(SQRT((RR-ddi-ddj)**2+ADD))**3 &
             +(RR+ddi+ddj)/(SQRT((RR+ddi+ddj)**2+ADD))**3

      DZQXX = two*(RR+ddi)/(SQRT((RR+ddi)**2+four*qqj2+ADQ))**3 &
             -two*(RR-ddi)/(SQRT((RR-ddi)**2+four*qqj2+ADQ))**3 &
             -two*(RR+ddi)/(SQRT((RR+ddi)**2+ADQ))**3 &
             +two*(RR-ddi)/(SQRT((RR-ddi)**2+ADQ))**3

      QXXDZ = two*(RR-ddj)/(SQRT((RR-ddj)**2+four*qqi2+AQD))**3 &
             -two*(RR+ddj)/(SQRT((RR+ddj)**2+four*qqi2+AQD))**3 &
             -two*(RR-ddj)/(SQRT((RR-ddj)**2+AQD))**3 &
             +two*(RR+ddj)/(SQRT((RR+ddj)**2+AQD))**3

      DZQZZ = (RR+ddi-two*qqj)/(SQRT((RR+ddi-two*qqj)**2+ADQ))**3 &
             -(RR-ddi-two*qqj)/(SQRT((RR-ddi-two*qqj)**2+ADQ))**3 &
             +(RR+ddi+two*qqj)/(SQRT((RR+ddi+two*qqj)**2+ADQ))**3 &
             -(RR-ddi+two*qqj)/(SQRT((RR-ddi+two*qqj)**2+ADQ))**3 &
             +two*(RR-ddi)/(SQRT((RR-ddi)**2+ADQ))**3 &
             -two*(RR+ddi)/(SQRT((RR+ddi)**2+ADQ))**3

      QZZDZ = (RR+two*qqi-ddj)/(SQRT((RR+two*qqi-ddj)**2+AQD))**3 &
             -(RR+two*qqi+ddj)/(SQRT((RR+two*qqi+ddj)**2+AQD))**3 &
             +(RR-two*qqi-ddj)/(SQRT((RR-two*qqi-ddj)**2+AQD))**3 &
             -(RR-two*qqi+ddj)/(SQRT((RR-two*qqi+ddj)**2+AQD))**3 &
             -two*(RR-ddj)/(SQRT((RR-ddj)**2+AQD))**3 &
             +two*(RR+ddj)/(SQRT((RR+ddj)**2+AQD))**3

      QXXQXX=-(two*RR)/(SQRT(RR2+four*(qqi-qqj)**2+AQQ))**3 &
             -(two*RR)/(SQRT(RR2+four*(qqi+qqj)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QXXQYY=-(four*RR)/(SQRT(RR2+four*qqi2+four*qqj2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QXXQZZ= &
             -two*(RR-two*qqj)/(SQRT((RR-two*qqj)**2+four*qqi2+AQQ))**3 &
             -two*(RR+two*qqj)/(SQRT((RR+two*qqj)**2+four*qqi2+AQQ))**3 &
             +two*(RR-two*qqj)/(SQRT((RR-two*qqj)**2+AQQ))**3 &
             +two*(RR+two*qqj)/(SQRT((RR+two*qqj)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QZZQXX= &
             -two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+four*qqj2+AQQ))**3 &
             -two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+four*qqj2+AQQ))**3 &
             +two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQQ))**3 &
             +two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QZZQZZ= &
           -(RR+two*qqi-two*qqj)/(SQRT((RR+two*qqi-two*qqj)**2+AQQ))**3 &
           -(RR+two*qqi+two*qqj)/(SQRT((RR+two*qqi+two*qqj)**2+AQQ))**3 &
           -(RR-two*qqi-two*qqj)/(SQRT((RR-two*qqi-two*qqj)**2+AQQ))**3 &
           -(RR-two*qqi+two*qqj)/(SQRT((RR-two*qqi+two*qqj)**2+AQQ))**3 &
             +two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQQ))**3 &
             +two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQQ))**3 &
             +two*(RR-2.D0*qqj)/(SQRT((RR-2.D0*qqj)**2+AQQ))**3 &
             +two*(RR+2.D0*qqj)/(SQRT((RR+2.D0*qqj)**2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      DXQXZ = two*(RR-qqj)/(SQRT((RR-qqj)**2+(ddi-qqj)**2+ADQ))**3 &
             -two*(RR+qqj)/(SQRT((RR+qqj)**2+(ddi-qqj)**2+ADQ))**3 &
             -two*(RR-qqj)/(SQRT((RR-qqj)**2+(ddi+qqj)**2+ADQ))**3 &
             +two*(RR+qqj)/(SQRT((RR+qqj)**2+(ddi+qqj)**2+ADQ))**3

      QXZDX = two*(RR+qqi)/(SQRT((RR+qqi)**2+(qqi-ddj)**2+AQD))**3 &
             -two*(RR-qqi)/(SQRT((RR-qqi)**2+(qqi-ddj)**2+AQD))**3 &
             -two*(RR+qqi)/(SQRT((RR+qqi)**2+(qqi+ddj)**2+AQD))**3 &
             +two*(RR-qqi)/(SQRT((RR-qqi)**2+(qqi+ddj)**2+AQD))**3

      QXZQXZ=-two*(RR+qqi-qqj)/(SQRT((RR+qqi-qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR+qqi+qqj)/(SQRT((RR+qqi+qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR-qqi-qqj)/(SQRT((RR-qqi-qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             -two*(RR-qqi+qqj)/(SQRT((RR-qqi+qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR+qqi-qqj)/(SQRT((RR+qqi-qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             -two*(RR+qqi+qqj)/(SQRT((RR+qqi+qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             -two*(RR-qqi-qqj)/(SQRT((RR-qqi-qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             +two*(RR-qqi+qqj)/(SQRT((RR-qqi+qqj)**2+(qqi+qqj)**2+AQQ))**3

      temp_real = DZDZ*fourth
      DGX(6)=TERMX*temp_real
      DGY(6)=TERMY*temp_real
      DGZ(6)=TERMZ*temp_real
      temp_real = DXDX*fourth
      DGX(7)=TERMX*temp_real
      DGY(7)=TERMY*temp_real
      DGZ(7)=TERMZ*temp_real
      temp_real = -(EDZ*half+QZZDZ*eighth)
      DGX(8)=TERMX*temp_real
      DGY(8)=TERMY*temp_real
      DGZ(8)=TERMZ*temp_real
      temp_real = -(EDZ*half+QXXDZ*eighth)
      DGX(9)=TERMX*temp_real
      DGY(9)=TERMY*temp_real
      DGZ(9)=TERMZ*temp_real
      temp_real = -QXZDX*eighth
      DGX(10)=TERMX*temp_real
      DGY(10)=TERMY*temp_real
      DGZ(10)=TERMZ*temp_real
      temp_real = -(DZE*half+DZQZZ*eighth) 
      DGX(13)=TERMX*temp_real
      DGY(13)=TERMY*temp_real
      DGZ(13)=TERMZ*temp_real
      temp_real = -(DZE*half+DZQXX*eighth)
      DGX(14)=TERMX*temp_real
      DGY(14)=TERMY*temp_real
      DGZ(14)=TERMZ*temp_real
      temp_real = -DXQXZ*eighth
      DGX(15)=TERMX*temp_real
      DGY(15)=TERMY*temp_real
      DGZ(15)=TERMZ*temp_real
      temp_real2 = EE+EQZZ*fourth
      temp_real = temp_real2+QZZE*fourth+QZZQZZ*sixteenth
      DGX(16)=TERMX*temp_real
      DGY(16)=TERMY*temp_real
      DGZ(16)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQZZ*sixteenth
      DGX(17)=TERMX*temp_real
      DGY(17)=TERMY*temp_real
      DGZ(17)=TERMZ*temp_real
      temp_real2 = EE+EQXX*fourth
      temp_real = temp_real2+QZZE*fourth+QZZQXX*sixteenth
      DGX(18)=TERMX*temp_real
      DGY(18)=TERMY*temp_real
      DGZ(18)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQXX*sixteenth
      DGX(19)=TERMX*temp_real
      DGY(19)=TERMY*temp_real
      DGZ(19)=TERMZ*temp_real
      temp_real = QXZQXZ*sixteenth
      DGX(20)=TERMX*temp_real
      DGY(20)=TERMY*temp_real
      DGZ(20)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQYY*sixteenth
      DGX(21)=TERMX*temp_real
      DGY(21)=TERMY*temp_real
      DGZ(21)=TERMZ*temp_real
      temp_real = (QXXQXX-QXXQYY)*thirtysecond
      DGX(22)=TERMX*temp_real
      DGY(22)=TERMY*temp_real
      DGZ(22)=TERMZ*temp_real
   end if

   IF(n_atomic_orbi>1.OR.n_atomic_orbj>1) then
      RXY2 = vec2_qm_qm1 + vec2_qm_qm2
      RYZ2 = vec2_qm_qm2 + vec2_qm_qm3
      RZX2 = vec2_qm_qm3 + vec2_qm_qm1
      LRXY2 = RXY2 < AXIS_TOL
      LRYZ2 = RYZ2 < AXIS_TOL
      LRZX2 = RZX2 < AXIS_TOL

      XTDX=zero; YTDX=zero; ZTDX=zero
      XTDY=zero; YTDY=zero; ZTDY=zero
      XTDZ=zero; YTDZ=zero; ZTDZ=zero

      if (.NOT.(LRXY2 .OR. LRYZ2 .OR. LRZX2)) then
        TERMX = vec_qm_qm_onerij1 * oneRIJ !vec_qm_qm1/(RIJ*RIJ)
        TERMY = vec_qm_qm_onerij2 * oneRIJ
        TERMZ = vec_qm_qm_onerij3 * oneRIJ
        oneRXY  = one/sqrt(RXY2)
        !Ross Walker - rearranged order here slightly for speed.
        TZ(3) = oneRIJ/oneRXY
        TZ3i = RIJ * oneRXY !Inverse of TZ(3) to avoid other divisions.
        TZ3i2 = TZ3i*TZ3i   !Square of 1/TZ(3)

        TX(1) = vec_qm_qm_onerij1
        TX(2) = vec_qm_qm_onerij2
        TX(3) = vec_qm_qm_onerij3

        TY(1) = -TX(2)*SIGN(one,TX(1)) * TZ3i
        TY(2) = ABS(TX(1) * TZ3i)
        TY(3) = zero

        TZ(1) = -TX(1)*TX(3)*TZ3i
        TZ(2) = -TX(2)*TX(3)*TZ3i
        !TZ(3) = see above

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
!       XTDY(3)=0.0D0
!       YTDY(3)=0.0D0
!       ZTDY(3)=0.0D0
                                                                                                                                                             
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
      elseif (LRXY2) then
!      MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=zero
         TX(2)=zero
         TX(3)=sign(one,vec_qm_qm3)
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
      elseif (LRYZ2) then
!      MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=sign(one,vec_qm_qm1)
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
      else !if (LRZX2) then
!      MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=zero
         TX(2)=sign(one,vec_qm_qm2)
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
      !At least one atom has more than one atomic orbital so need to consider S and P interactions.
      isp = 0
      do K=qm_atomi_orb_start,total_atomic_orb
         KK=K-qm_atomi_orb_start
         do L=K,total_atomic_orb
            LL=L-qm_atomi_orb_start
            do M=1,n_atomic_orbj
               MM=M-1
               do N=M,n_atomic_orbj
                  NN=N-1
                  ISP=ISP+1
                  IF(NN == 0)THEN
                     IF(LL == 0) THEN
!   (SS/SS)
                        DRX(ISP)=DGX(1)
                        DRY(ISP)=DGY(1)
                        DRZ(ISP)=DGZ(1)
                     ELSEIF(KK == 0) THEN
!   (SP/SS)
                        DRX(ISP)=DGX(2)*TX(LL)+qm_qm_e_repul(2)*xTDX(LL)
                        DRY(ISP)=DGY(2)*TX(LL)+qm_qm_e_repul(2)*yTDX(LL)
                        DRZ(ISP)=DGZ(2)*TX(LL)+qm_qm_e_repul(2)*zTDX(LL)
                     ELSE
!   (PP/SS)
                        DRX(ISP)=DGX(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(xTDX(KK)    &
                                *TX(LL)+TX(KK)*xTDX(LL))+DGX(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(xTDY(KK)*TY(LL) &
                                +TY(KK)*xTDY(LL)+xTDZ(KK)*TZ(LL)+TZ(KK)*xTDZ(LL))
                        DRY(ISP)=DGY(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(yTDX(KK)    &
                                *TX(LL)+TX(KK)*yTDX(LL))+DGY(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(yTDY(KK)*TY(LL) &
                                +TY(KK)*yTDY(LL)+yTDZ(KK)*TZ(LL)+TZ(KK)*yTDZ(LL))
                        DRZ(ISP)=DGZ(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(zTDX(KK)    &
                                *TX(LL)+TX(KK)*zTDX(LL))+DGZ(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(zTDY(KK)*TY(LL) &
                                +TY(KK)*zTDY(LL)+zTDZ(KK)*TZ(LL)+TZ(KK)*zTDZ(LL))
                     ENDIF
                  ELSEIF(MM == 0) THEN
                     IF(LL == 0) THEN
!   (SS/SP)
                        DRX(ISP)=DGX(5)*TX(NN)+qm_qm_e_repul(5)*xTDX(NN)
                        DRY(ISP)=DGY(5)*TX(NN)+qm_qm_e_repul(5)*yTDX(NN)
                        DRZ(ISP)=DGZ(5)*TX(NN)+qm_qm_e_repul(5)*zTDX(NN)
                     ELSEIF(KK == 0) THEN
!   (SP/SP)
                        DRX(ISP)=DGX(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(xTDX(LL)*TX(NN) &
                                +TX(LL)*xTDX(NN))+DGX(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)     &
                                +xTDZ(LL)*TZ(NN)+TZ(LL)*xTDZ(NN))
                        DRY(ISP)=DGY(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(yTDX(LL)*TX(NN) &
                                +TX(LL)*yTDX(NN))+DGY(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)     &
                                +yTDZ(LL)*TZ(NN)+TZ(LL)*yTDZ(NN))
                        DRZ(ISP)=DGZ(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(zTDX(LL)*TX(NN) &
                                +TX(LL)*zTDX(NN))+DGZ(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)     &
                                +zTDZ(LL)*TZ(NN)+TZ(LL)*zTDZ(NN))
                     ELSE
!   (PP/SP)
                        DRX(ISP)=DGX(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(xTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*xTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *xTDX(NN))+DGX(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((xTDY(KK)*TY(LL)+TY(KK)*xTDY(LL)     &
                                +xTDZ(KK)*TZ(LL)+TZ(KK)*xTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*xTDX(NN))+DGX(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(xTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+xTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)+xTDZ(LL)*TZ(NN)+TZ(LL) &
                                *xTDZ(NN))+TX(LL)*(xTDY(KK)*TY(NN)+TY(KK)*xTDY(NN)       &
                                +xTDZ(KK)*TZ(NN)+TZ(KK)*xTDZ(NN)))
                        DRY(ISP)=DGY(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(yTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*yTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *yTDX(NN))+DGY(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((yTDY(KK)*TY(LL)+TY(KK)*yTDY(LL)     &
                                +yTDZ(KK)*TZ(LL)+TZ(KK)*yTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*yTDX(NN))+DGY(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(yTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+yTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)+yTDZ(LL)*TZ(NN)+TZ(LL) &
                                *yTDZ(NN))+TX(LL)*(yTDY(KK)*TY(NN)+TY(KK)*yTDY(NN)       &
                                +yTDZ(KK)*TZ(NN)+TZ(KK)*yTDZ(NN)))
                        DRZ(ISP)=DGZ(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(zTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*zTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *zTDX(NN))+DGZ(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((zTDY(KK)*TY(LL)+TY(KK)*zTDY(LL)     &
                                +zTDZ(KK)*TZ(LL)+TZ(KK)*zTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*zTDX(NN))+DGZ(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(zTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+zTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)+zTDZ(LL)*TZ(NN)+TZ(LL) &
                                *zTDZ(NN))+TX(LL)*(zTDY(KK)*TY(NN)+TY(KK)*zTDY(NN)       &
                                +zTDZ(KK)*TZ(NN)+TZ(KK)*zTDZ(NN)))
                     ENDIF
                  ELSEIF(LL == 0) THEN
!   (SS/PP)
                     DRX(ISP)=DGX(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(xTDX(MM)*TX(NN)+TX(MM) &
                             *xTDX(NN))+DGX(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(xTDY(MM)*TY(NN)+TY(MM)*xTDY(NN)+xTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*xTDZ(NN))
                     DRY(ISP)=DGY(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(yTDX(MM)*TX(NN)+TX(MM) &
                             *yTDX(NN))+DGY(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(yTDY(MM)*TY(NN)+TY(MM)*yTDY(NN)+yTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*yTDZ(NN))
                     DRZ(ISP)=DGZ(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(zTDX(MM)*TX(NN)+TX(MM) &
                             *zTDX(NN))+DGZ(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(zTDY(MM)*TY(NN)+TY(MM)*zTDY(NN)+zTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*zTDZ(NN))
                  ELSEIF(KK == 0) THEN
!   (SP/PP)
                     DRX(ISP)=DGX(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(xTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*xTDX(MM)*TX(NN)+TX(LL)*TX(MM)*xTDX(NN))+DGX(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(xTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(xTDY(MM)*TY(NN)&
                             +TY(MM)*xTDY(NN)+xTDZ(MM)*TZ(NN)+TZ(MM)*xTDZ(NN)))+DGX(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(xTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+xTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(xTDY(MM)&
                             *TX(NN)+TY(MM)*xTDX(NN)+xTDY(NN)*TX(MM)+TY(NN)*xTDX(MM))+TZ(LL)  &
                             *(xTDZ(MM)*TX(NN)+TZ(MM)*xTDX(NN)+xTDZ(NN)*TX(MM)+TZ(NN)*xTDX(MM)))
                     DRY(ISP)=DGY(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(yTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*yTDX(MM)*TX(NN)+TX(LL)*TX(MM)*yTDX(NN))+DGY(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(yTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(yTDY(MM)*TY(NN)&
                             +TY(MM)*yTDY(NN)+yTDZ(MM)*TZ(NN)+TZ(MM)*yTDZ(NN)))+DGY(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(yTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+yTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(yTDY(MM)&
                             *TX(NN)+TY(MM)*yTDX(NN)+yTDY(NN)*TX(MM)+TY(NN)*yTDX(MM))+TZ(LL)  &
                             *(yTDZ(MM)*TX(NN)+TZ(MM)*yTDX(NN)+yTDZ(NN)*TX(MM)+TZ(NN)*yTDX(MM)))
                     DRZ(ISP)=DGZ(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(zTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*zTDX(MM)*TX(NN)+TX(LL)*TX(MM)*zTDX(NN))+DGZ(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(zTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(zTDY(MM)*TY(NN)&
                             +TY(MM)*zTDY(NN)+zTDZ(MM)*TZ(NN)+TZ(MM)*zTDZ(NN)))+DGZ(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(zTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+zTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(zTDY(MM)&
                             *TX(NN)+TY(MM)*zTDX(NN)+zTDY(NN)*TX(MM)+TY(NN)*zTDX(MM))+TZ(LL)  &
                             *(zTDZ(MM)*TX(NN)+TZ(MM)*zTDX(NN)+zTDZ(NN)*TX(MM)+TZ(NN)*zTDX(MM)))
                  ELSE
!   (PP/PP)
                     DRX(ISP)=DGX(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(xTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*xTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*xTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*xTDX(NN))+DGX(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((xTDY(KK)*TY(LL)+TY(KK)*xTDY(LL)+xTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *xTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(xTDX(MM)*TX(NN)+TX(MM)    &
                             *xTDX(NN)))+DGX(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((xTDX(KK)*TX(LL)+TX(KK)*xTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(xTDY(MM)*TY(NN)+TY(MM)*xTDY(NN)+xTDZ(MM)*TZ(NN)+TZ(MM)*xTDZ(NN)))
                     DRY(ISP)=DGY(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(yTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*yTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*yTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*yTDX(NN))+DGY(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((yTDY(KK)*TY(LL)+TY(KK)*yTDY(LL)+yTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *yTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(yTDX(MM)*TX(NN)+TX(MM)    &
                             *yTDX(NN)))+DGY(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((yTDX(KK)*TX(LL)+TX(KK)*yTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(yTDY(MM)*TY(NN)+TY(MM)*yTDY(NN)+yTDZ(MM)*TZ(NN)+TZ(MM)*yTDZ(NN)))
                     DRZ(ISP)=DGZ(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(zTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*zTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*zTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*zTDX(NN))+DGZ(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((zTDY(KK)*TY(LL)+TY(KK)*zTDY(LL)+zTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *zTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(zTDX(MM)*TX(NN)+TX(MM)    &
                             *zTDX(NN)))+DGZ(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((zTDX(KK)*TX(LL)+TX(KK)*zTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(zTDY(MM)*TY(NN)+TY(MM)*zTDY(NN)+zTDZ(MM)*TZ(NN)+TZ(MM)*zTDZ(NN)))
                     DRX(ISP)=DRX(ISP)+DGX(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(xTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*xTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*xTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*xTDY(NN)+xTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*xTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*xTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*xTDZ(NN))+DGX(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     DRY(ISP)=DRY(ISP)+DGY(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(yTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*yTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*yTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*yTDY(NN)+yTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*yTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*yTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*yTDZ(NN))+DGY(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     DRZ(ISP)=DRZ(ISP)+DGZ(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(zTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*zTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*zTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*zTDY(NN)+zTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*zTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*zTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*zTDZ(NN))+DGZ(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
!      TO AVOID COMPILER DIFFICULTIES THIS IS DIVIDED
                     xTEMP1=  xTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+xTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(xTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+xTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(xTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+xTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     yTEMP1=  yTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+yTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(yTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+yTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(yTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+yTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     zTEMP1=  zTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+zTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(zTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+zTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(zTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+zTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                                                                                                                                                             
                     xTEMP2=  TX(KK)*(TX(MM)*(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)+xTDZ(LL)*TZ(NN)+TZ(LL)*xTDZ(NN))    &
                             +TX(NN)*(xTDY(LL)*TY(MM)+TY(LL)*xTDY(MM)+xTDZ(LL)*TZ(MM)+TZ(LL)*xTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(xTDY(KK)*TY(NN)+TY(KK)*xTDY(NN)+xTDZ(KK)*TZ(NN)+TZ(KK)*xTDZ(NN))+TX(NN)   &
                             *(xTDY(KK)*TY(MM)+TY(KK)*xTDY(MM)+xTDZ(KK)*TZ(MM)+TZ(KK)*xTDZ(MM)))
                     yTEMP2=  TX(KK)*(TX(MM)*(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)+yTDZ(LL)*TZ(NN)+TZ(LL)*yTDZ(NN))    &
                             +TX(NN)*(yTDY(LL)*TY(MM)+TY(LL)*yTDY(MM)+yTDZ(LL)*TZ(MM)+TZ(LL)*yTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(yTDY(KK)*TY(NN)+TY(KK)*yTDY(NN)+yTDZ(KK)*TZ(NN)+TZ(KK)*yTDZ(NN))+TX(NN)   &
                             *(yTDY(KK)*TY(MM)+TY(KK)*yTDY(MM)+yTDZ(KK)*TZ(MM)+TZ(KK)*yTDZ(MM)))
                     zTEMP2=  TX(KK)*(TX(MM)*(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)+zTDZ(LL)*TZ(NN)+TZ(LL)*zTDZ(NN))    &
                             +TX(NN)*(zTDY(LL)*TY(MM)+TY(LL)*zTDY(MM)+zTDZ(LL)*TZ(MM)+TZ(LL)*zTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(zTDY(KK)*TY(NN)+TY(KK)*zTDY(NN)+zTDZ(KK)*TZ(NN)+TZ(KK)*zTDZ(NN))+TX(NN)   &
                             *(zTDY(KK)*TY(MM)+TY(KK)*zTDY(MM)+zTDZ(KK)*TZ(MM)+TZ(KK)*zTDZ(MM)))
                                                                                                                                                             
                     DRX(ISP)=DRX(ISP)+qm_qm_e_repul(20)*(xTEMP1+xTEMP2)
                     DRY(ISP)=DRY(ISP)+qm_qm_e_repul(20)*(yTEMP1+yTEMP2)
                     DRZ(ISP)=DRZ(ISP)+qm_qm_e_repul(20)*(zTEMP1+zTEMP2)
                     DRX(ISP)=DRX(ISP)+DGX(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(xTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*xTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*xTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*xTDZ(NN)+xTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*xTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*xTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*xTDY(NN))
                     DRY(ISP)=DRY(ISP)+DGY(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(yTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*yTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*yTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*yTDZ(NN)+yTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*yTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*yTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*yTDY(NN))
                     DRZ(ISP)=DRZ(ISP)+DGZ(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(zTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*zTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*zTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*zTDZ(NN)+zTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*zTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*zTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*zTDY(NN))
                                                                                                                                                             
                     DRX(ISP)=DRX(ISP)+DGX(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((xTDY(KK)*TZ(LL)+TY(KK)*xTDZ(LL)+xTDZ(KK)*TY(LL)+TZ(KK)        &
                             *xTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(xTDY(MM)  &
                             *TZ(NN)+TY(MM)*xTDZ(NN)+xTDZ(MM)*TY(NN)+TZ(MM)*xTDY(NN)))
                     DRY(ISP)=DRY(ISP)+DGY(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((yTDY(KK)*TZ(LL)+TY(KK)*yTDZ(LL)+yTDZ(KK)*TY(LL)+TZ(KK)        &
                             *yTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(yTDY(MM)  &
                             *TZ(NN)+TY(MM)*yTDZ(NN)+yTDZ(MM)*TY(NN)+TZ(MM)*yTDY(NN)))
                     DRZ(ISP)=DRZ(ISP)+DGZ(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((zTDY(KK)*TZ(LL)+TY(KK)*zTDZ(LL)+zTDZ(KK)*TY(LL)+TZ(KK)        &
                             *zTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(zTDY(MM)  &
                             *TZ(NN)+TY(MM)*zTDZ(NN)+zTDZ(MM)*TY(NN)+TZ(MM)*zTDY(NN)))
                  ENDIF
               end do !N=M,n_atomic_orbj
            end do !M=1,n_atomic_orbj
         end do !L=K,total_atomic_orb
      end do !K=qm_atomi_orb_start,total_atomic_orb
   else
   ! Just have 1 atomic orbital interacting with 1 atomic orbital so only need
   ! to do (SS/SS)
      DRX(1)=DGX(1)
      DRY(1)=DGY(1)
      DRZ(1)=DGZ(1)
   end if !IF(n_atomic_orbi>1.OR.n_atomic_orbj>1)

!   THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM                              
                                                                              
!      CORE-CORE TERMS, MNDO, AM1, RM1 and PM3
                                                                               
!  SPECIAL TREATMENT FOR N-H AND O-H TERMS                                      
   IF(natqmi == 1 .AND. (natqmj == 7 .OR. natqmj == 8)) THEN
      F3=one+EXP1i+RIJ*EXP1j
      temp_real3 = alphai*EXP1i+(alphaj*RIJ-one)*EXP1j
      temp_real3 = temp_real3 * qm_qm_e_repul(1)
      DDX=(dgx(1)*F3-vec_qm_qm_oneRIJ1*temp_real3)*C1
      DDY=(dgy(1)*F3-vec_qm_qm_oneRIJ2*temp_real3)*C1
      DDZ=(dgz(1)*F3-vec_qm_qm_oneRIJ3*temp_real3)*C1
   ELSEIF((natqmi == 7 .OR. natqmi == 8).AND. natqmj == 1) THEN
      F3=one+EXP1j+RIJ*EXP1i
      temp_real3 = alphaj*EXP1j+(alphai*RIJ-one)*EXP1i
      temp_real3 = temp_real3 * qm_qm_e_repul(1)
      DDX=(dgx(1)*F3-vec_qm_qm_oneRIJ1*temp_real3)*C1
      DDY=(dgy(1)*F3-vec_qm_qm_oneRIJ2*temp_real3)*C1
      DDZ=(dgz(1)*F3-vec_qm_qm_oneRIJ3*temp_real3)*C1
   ELSE
      PART1x=dgx(1)*C1
      PART1y=dgy(1)*C1
      PART1z=dgz(1)*C1
      temp_real3 = (EXP1i+EXP1j)*ABS(C1)
      PART3x=dgx(1)*temp_real3
      PART3y=dgy(1)*temp_real3
      PART3z=dgz(1)*temp_real3
      temp_real4 = -qm_qm_e_repul(1)*(alphai*EXP1i+alphaj*EXP1j)*ABS(C1)
      PART2x=temp_real4*vec_qm_qm_oneRIJ1
      PART2y=temp_real4*vec_qm_qm_oneRIJ2
      PART2z=temp_real4*vec_qm_qm_oneRIJ3
      DDx=PART1x+PART2x+PART3x
      DDy=PART1y+PART2y+PART3y
      DDz=PART1z+PART2z+PART3z
   ENDIF
   FNUCX=DDX+FNUCX
   FNUCY=DDY+FNUCY
   FNUCZ=DDZ+FNUCZ


!  FIRST, CORE-ELECTRON ATTRACTION DERIVATIVES (MNDO, AM1 and PM3)             
!  ATOM CORE I AFFECTING A.O.S ON J                                     
   ISP=0                                                               
   do M=1,n_atomic_orbj                                                      
      BB=one                                                        
      do N=M,n_atomic_orbj                                                    
         MN=M+qm2_params%pascal_tri1(N)
         ISP=ISP+1      
         temp_real = BB*corei*PSUM(MN)                                               
         FABx=FABx-temp_real*DRX(ISP)                    
         FABy=FABy-temp_real*DRY(ISP)                    
         FABz=FABz-temp_real*DRZ(ISP)                    
         BB=two                                                       
      end do
   end do
!  ATOM CORE J AFFECTING A.O.S ON I                                     
   K=n_atomic_orbj
   K=qm2_params%pascal_tri2(n_atomic_orbj)
   ISP=-K+1                                                            
   do M=qm_atomi_orb_start,total_atomic_orb                                                      
      BB=one
      do N=M,total_atomic_orb                                                    
         MN=M+qm2_params%pascal_tri1(N)
         ISP=ISP+K              
         temp_real = BB*corej*PSUM(MN)                                       
         FABx=FABx-temp_real*DRX(ISP)                    
         FABy=FABy-temp_real*DRY(ISP)                    
         FABz=FABz-temp_real*DRZ(ISP)                    
         BB=two
      end do
   end do

!   NOW FOR COULOMB AND EXCHANGE TERMS (MNDO, AM1 and PM3)                           
   ISP=0                                                               
   do K=qm_atomi_orb_start,total_atomic_orb                                                      
      AA=one                                                     
      KK=qm2_params%pascal_tri1(K)
      do L=K,total_atomic_orb                                                    
         LL=qm2_params%pascal_tri1(L)
         KL=K+LL                                                 
         do M=1,n_atomic_orbj                                                
            BB=one                                                   
            MK=M+KK                                                 
            ML=M+LL                                                 
            do N=M,n_atomic_orbj                                              
               ISP=ISP+1                                               
               MN=M+qm2_params%pascal_tri1(N)
!    COULOMB TERM
               temp_real=AA*BB*PSUM(KL)*PSUM(MN)                                                              
               FAAx=FAAx+temp_real*DRX(ISP)           
               FAAy=FAAy+temp_real*DRY(ISP)           
               FAAz=FAAz+temp_real*DRZ(ISP)           
!  EXCHANGE TERM                                                              
               NK=N+KK                                                 
               NL=N+LL 
               temp_real = AA*BB*0.25D0*(PSUM(MK)*PSUM(NL)+PSUM(NK)*PSUM(ML))                                 
               FAAx = FAAx-temp_real*DRX(ISP) 
               FAAy = FAAy-temp_real*DRY(ISP) 
               FAAz = FAAz-temp_real*DRZ(ISP) 
               BB=two
            end do !N=M,n_atomic_orbj
         end do !M=1,n_atomic_orbj
         AA=two
      end do ! L=K,total_atomic_orb
   end do ! K=qm_atomi_orb_start,total_atomic_orb

   pair_force(1)=FAAx+FABx+FNUCX
   pair_force(1) = -pair_force(1)*EV_TO_KCAL                                 
   pair_force(2)=FAAy+FABy+FNUCY                                           
   pair_force(2) = -pair_force(2)*EV_TO_KCAL                                 
   pair_force(3)=FAAz+FABz+FNUCZ
   pair_force(3) = -pair_force(3)*EV_TO_KCAL                                 

   RETURN                                                                    
end subroutine qm2_deriv_qm_analyt

subroutine qm2_dhc(P,iqm, jqm,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi, &
                  natqmj,iif,iil,jjf,jjl,DENER)
!***********************************************************************
!
! Ross Walker (SDSC, 2006) : Do 'Pseudo' Numerical Derivatives for QM
!
!***********************************************************************
      use qmmm_module, only : qm2_params, OVERLAP_CUTOFF, qmmm_nml, qm2_struct
      use constants, only : A2_TO_BOHRS2, EV_TO_KCAL

      implicit none

!Passed in
      _REAL_ P(*)
      _REAL_, intent(in)  :: xyz_qmi(3),xyz_qmj(3)
      integer, intent(in) :: iqm, jqm, natqmi, natqmj, qmitype, qmjtype
      integer, intent(in) :: iif,iil,jjf,jjl
      _REAL_, intent(out) :: DENER

! Local
      integer n_atomic_orbi, n_atomic_orbj, ja, jb, ia,ib,i,j
      integer j1,jj,i1, linear, i2, ii
! dimension 36 = max of 8*(8+1)/2 = 4 orbs with 4 orbs - S,3P with S,3P
      _REAL_ H(36), F(36)
      _REAL_ SHMAT(4,4),E1B(10), E2A(10), W(100)
      _REAL_ enuclr, ee, vec_qm_qm1, vec_qm_qm2, vec_qm_qm3, R2
      integer orb_loc(2,2),KR

      !qm2_Helect is a function
      _REAL_ qm2_helect

      orb_loc(1,1)=1
      orb_loc(2,1)=iil-iif+1
      orb_loc(1,2)=orb_loc(2,1)+1
      orb_loc(2,2)=orb_loc(1,2)+jjl-jjf
      n_atomic_orbi = orb_loc(2,2)-orb_loc(1,2)+1
      n_atomic_orbj = orb_loc(2,1)-orb_loc(1,1)+1
      LINEAR=qm2_params%pascal_tri2(orb_loc(2,2))
      do I=1,LINEAR
         F(I)=0.0D0
         H(I)=0.0D0
      end do
      IA=orb_loc(1,1)
      IB=orb_loc(2,1)
      JA=orb_loc(1,2)
      JB=orb_loc(2,2)
      J=2
      I=1
! RCW: Caution, i and j reversed here.

      vec_qm_qm1 = (xyz_qmj(1)-xyz_qmi(1))
      vec_qm_qm2 = (xyz_qmj(2)-xyz_qmi(2))
      vec_qm_qm3 = (xyz_qmj(3)-xyz_qmi(3))
      R2 = (vec_qm_qm1*vec_qm_qm1+vec_qm_qm2*vec_qm_qm2+vec_qm_qm3*vec_qm_qm3)*A2_TO_BOHRS2
      if (R2 < OVERLAP_CUTOFF) then
        SHMAT=0.0d0
        CALL qm2_h1elec(R2,xyz_qmj(1),                                  &
               xyz_qmi(1),n_atomic_orbj, n_atomic_orbi, SHMAT,           &
               qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmjtype,qmitype), &
               qm2_params%atom_orb_ss_eqn(1,1,qmjtype,qmitype),          &
               qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
               qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
               qm2_params%atom_orb_sp_ovlp(1,1,qmjtype,qmitype),         &
               qm2_params%atom_orb_sp_ovlp(1,1,qmitype,qmjtype),         &
               qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmjtype,qmitype), &
               qm2_params%atom_orb_pp_ovlp_ieqj1(1,1,qmjtype,qmitype),   &
               qm2_params%atom_orb_pp_ovlp_ieqj2(1,1,qmjtype,qmitype),   &
               qm2_params%atom_orb_pp_ovlp_inj(1,1,qmjtype,qmitype),     &
               qm2_params%betasas(qmjtype,qmitype),                      &
               qm2_params%betasap(qmjtype,qmitype),                      &
               qm2_params%betasap(qmitype,qmjtype),                      &
               qm2_params%betapap(qmjtype,qmitype))
        J1=0
        do J=JA,JB
           JJ=qm2_params%pascal_tri1(j)
           J1=J1+1
           I1=0
           do I=IA,IB
              JJ=JJ+1
              I1=I1+1
              H(JJ)=SHMAT(I1,J1)
              F(JJ)=SHMAT(I1,J1)
           end do
        end do
      end if !(R2 < OVERLAP_CUTOFF)

      KR=1
      CALL qm2_rotate_qmqm(-1,iqm,jqm,natqmi,natqmj,xyz_qmi,xyz_qmj,            &
                  W(KR),KR,E2A,E1B,ENUCLR,qmitype,qmjtype)

!    * ENUCLR IS SUMMED OVER CORE-CORE REPULSION INTEGRALS.
      I2=0
      do I1=IA,IB
         II=qm2_params%pascal_tri1(i1)+IA-1
         do J1=IA,I1
            II=II+1
            I2=I2+1
            H(II)=H(II)+E1B(I2)
            F(II)=F(II)+E1B(I2)
         end do
      end do
      I2=0
      do I1=JA,JB
         II=qm2_params%pascal_tri1(i1)+JA-1
         do J1=JA,I1
            II=II+1
            I2=I2+1
            H(II)=H(II)+E2A(I2)
            F(II)=F(II)+E2A(I2)
         end do
      end do
      CALL qm2_fock2_2atm(F,P,W,orb_loc)
      EE=qm2_HELECT(orb_loc(2,2)-1,P,H,F)
      DENER=EE+ENUCLR
      RETURN
end subroutine qm2_dhc




