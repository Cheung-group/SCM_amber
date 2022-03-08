! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_hcore_qmqm(COORD,H,W,ENUCLR)
!***********************************************************************
! Current code, optimisation and inlining by: Ross Walker (TSRI, 2005)
!
! This routine is responsible for generating the one-electron matrix and
! the two electron integrals via calls to qm2_h1elec and qm2_rotate_qmqm.
! qm2_h1elec has been inlined in this code for speed.
!
! Current Version: Ross Walker (TSRI, 2005)
!
!IN -
! COORD = QM coordinates
!
!OUT-
! H = One electron matix
! W = Two electron integrals
! ENUCLR = Nuclear energy 
!***********************************************************************

      use qmmm_module, only : qmmm_nml, qmmm_struct, qm2_struct, qm2_params, qm2_rij_eqns, &
                              qmmm_mpi, OVERLAP_CUTOFF
      use constants, only : A2_TO_BOHRS2, A_TO_BOHRS
      implicit none

#include "qm2_array_locations.h"
!Passed in
      _REAL_, intent(in) :: COORD(3,qmmm_struct%nquant_nlink)
      _REAL_, intent(out) :: W(qm2_struct%n2el)
      _REAL_, intent(out) ::  ENUCLR
      _REAL_, intent(out) :: H(qm2_struct%matsize)

!Local
      _REAL_ E1B(10),E2A(10),SHMAT(4,4)
      integer i, first_si, last_pi, first_pi, ni, i1, i2, j
      integer kr, j1, first_sj, last_pj, ii, j2, jj
      integer loop_count, jstart, jend
      integer n_atomic_orbi, n_atomic_orbj, qmitype, qmjtype
      _REAL_ enuc, elec_ke_p, vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
      _REAL_ half_num, R2

! FILL THE DIAGONALS as we don't do them in the loop below.
      do i = qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
        first_si=qm2_params%orb_loc(1,i)
        first_pi=first_si+1
        last_pi=qm2_params%orb_loc(2,i)

        I2=qm2_params%pascal_tri2(first_si)
        H(i2)=qm2_params%orb_elec_ke(1,i)
        elec_ke_p=qm2_params%orb_elec_ke(2,i)
        do I1=first_pi,last_pi
          I2=qm2_params%pascal_tri2(I1)
          H(I2) = elec_ke_p
        end do
      end do
      loop_count=0
#ifdef MPI
      KR = qmmm_mpi%two_e_offset+1
      do i = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
         jstart =  qmmm_mpi%nquant_nlink_jrange(1,i)
         jend = qmmm_mpi%nquant_nlink_jrange(2,i)
#else
      KR=1                                                                      
      do I=2,qmmm_struct%nquant_nlink
         jstart = 1
         jend = i-1
#endif
         first_si=qm2_params%orb_loc(1,I)
         last_pi=qm2_params%orb_loc(2,I)                
         n_atomic_orbi = qm2_params%natomic_orbs(i)
         qmitype = qmmm_struct%qm_atom_type(i) 
         NI=qmmm_struct%iqm_atomic_numbers(I)
!   FILL THE ATOM-OTHER ATOM ONE-ELECTRON MATRIX<PSI(LAMBDA)|PSI(SIGMA)>        
         do J=jstart, jend
            loop_count=loop_count+1
            first_sj=qm2_params%orb_loc(1,J)
            last_pj=qm2_params%orb_loc(2,J) 
            qmjtype = qmmm_struct%qm_atom_type(j)

            !Calculate Overlap Integrals using a Gaussian Expansion
            !STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970
            !Fill SHMAT with a 4x4 array of overlaps, in order S,PX,PY,PZ
            !     R2   =  INTERATOMIC DISTANCE^2 IN BOHRS2
            if (qmmm_nml%qmqmrij_incore) then
              R2=qm2_rij_eqns%qmqmrijdata(QMQMRIJBOHRS2,loop_count)
            else
              vec_qm_qm1 = (coord(1,i)-coord(1,j))
              vec_qm_qm2 = (coord(2,i)-coord(2,j))
              vec_qm_qm3 = (coord(3,i)-coord(3,j))
              R2 = (vec_qm_qm1*vec_qm_qm1+vec_qm_qm2*vec_qm_qm2+vec_qm_qm3*vec_qm_qm3)*A2_TO_BOHRS2
            end if
            if (R2 < OVERLAP_CUTOFF) then
              n_atomic_orbj = qm2_params%natomic_orbs(j)
              SHMAT=0.0d0

              CALL qm2_h1elec(R2,COORD(1,I),COORD(1,J),n_atomic_orbi,n_atomic_orbj,SHMAT, &
                            qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_ss_eqn(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
                            qm2_params%atom_orb_sp_ovlp(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_sp_ovlp(1,1,qmjtype,qmitype), &
                            qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_pp_ovlp_ieqj1(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_pp_ovlp_ieqj2(1,1,qmitype,qmjtype), &
                            qm2_params%atom_orb_pp_ovlp_inj(1,1,qmitype,qmjtype), &
                            qm2_params%betasas(qmitype,qmjtype),qm2_params%betasap(qmitype,qmjtype), &
                            qm2_params%betasap(qmjtype,qmitype), qm2_params%betapap(qmitype,qmjtype))

              I2=0
              do I1=first_si,last_pi
                 II=qm2_params%pascal_tri1(i1)+first_sj-1
                 I2=I2+1
                 J2=0
                 JJ=MIN(I1,last_pj)
                 do J1=first_sj,JJ                                                   
                    II=II+1                                                       
                    J2=J2+1                                                       
                    H(II)=H(II)+SHMAT(I2,J2) 
                 end do
              end do
            end if !(R2 < OVERLAP_CUTOFF)
                                                                   
!   CALCULATE THE TWO-ELECTRON INTEGRALS, W; THE ELECTRON NUCLEAR TERMS         
!   E1B AND E2A; AND THE NUCLEAR-NUCLEAR TERM ENUC.                             
            CALL qm2_rotate_qmqm(loop_count,i,j,NI,qmmm_struct%iqm_atomic_numbers(J),COORD(1,I),COORD(1,J), &
                       W(KR),KR,E1B,E2A,ENUC,qmitype,qmjtype)

            ENUCLR = ENUCLR + ENUC                    
                                                                               
!   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM I.                     
           if(I == J) then
             half_num=0.5D0
           else
             half_num=1.0D0
           end if
                                                                               
           I2=0                                                                
           do I1=first_si,last_pi                                                      
              II=qm2_params%pascal_tri1(i1)+first_si-1
              do J1=first_si,I1
                 II=II+1                                                       
                 I2=I2+1                                                       
                 H(II)=H(II)+E1B(I2)*HALF_num
              end do
           end do
                                                                               
!   ADD ON THE ELECTRON-NUCLEAR ATTRACTION TERM FOR ATOM J.                     
                                                                               
           I2=0                                                                
           do I1=first_sj,last_pj
              II=qm2_params%pascal_tri1(i1)+first_sj-1                                              
              do J1=first_sj,I1                                                   
                 II=II+1                                                       
                 I2=I2+1                                                       
                 H(II)=H(II)+E2A(I2)*HALF_num
              end do
           end do
         end do  ! J=1,iminus
      end do !  I=1,qmmm_struct%nquant_nlink
      RETURN                                                                    
end subroutine qm2_hcore_qmqm

                                                                       
