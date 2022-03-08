! <compile=optimized>
#include "copyright.h"
#include "dprec.h"
subroutine qm2_h1elec(R2,XI,XJ, &
                      n_atomic_orbi,n_atomic_orbj,SHMAT,sxs_over_sas,ss_eqn, &
                      sxp_over_sap,pxs_over_pas,sp_ovlp,ps_ovlp,pxp_over_pap, &
                      pp_ovlp_ieqj1, pp_ovlp_ieqj2,pp_ovlp_inj,betasas,betasap, &
                      betapas,betapap)
!***********************************************************************
!
!  qm2_h1elec forms the one-electron matrix between two atoms and
!  calculates the overlaps.
!
!  Current code optimised by Ross Walker (TSRI, 2005)
!
!   ON INPUT
!               XI   = COORDINATES OF FIRST ATOM.
!               XJ   = COORDINATES OF SECOND ATOM.
!  n_atomic_orbi,j   = number of atomic orbitals on i and j.
!                                                 
!   ON OUTPUT   SHMAT = MATRIX OF ONE-ELECTRON INTERACTIONS.
!                                                           
!***********************************************************************
      use qmmm_module, only : EXPONENTIAL_CUTOFF
      use constants, only : A2_TO_BOHRS2, A_TO_BOHRS, half
      implicit none

!Passed In
      _REAL_, intent(in) :: R2,XI(3),XJ(3)
      integer, intent(in) :: n_atomic_orbi, n_atomic_orbj
      _REAL_, intent(out) :: SHMAT(4,4)
      _REAL_, intent(in) :: betasas,betasap,betapas,betapap
      !These are passed in here rather than used from the module to avoid
      !having to do a look up on a 4 dimensional array within a loop.
      _REAL_, intent(in) :: sxs_over_sas(6,6), ss_eqn(6,6), sxp_over_sap(6,6)
      _REAL_, intent(in) :: pxs_over_pas(6,6), sp_ovlp(6,6), ps_ovlp(6,6)
      _REAL_, intent(in) :: pxp_over_pap(6,6), pp_ovlp_ieqj1(6,6), pp_ovlp_ieqj2(6,6)
      _REAL_, intent(in) :: pp_ovlp_inj(6,6)
      
!Local
      _REAL_ BIJSP, BIJPS, BIJPP
      _REAL_ SH, vec_qm_qm(3)
      _REAL_ ADBR2, TOMB
      integer i,j,k,l, ii,jj


!***********************************************************************        
!   CALCULATE THE OVERLAP INTEGRALS USING A GAUSSIAN EXPANSION         *        
!         STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970     *        
!                                                                      *
!         FILL SHMAT=  4X4 ARRAY OF OVERLAPS, IN ORDER S,PX,PY,PZ      *        
!***********************************************************************        

!------------------------------
!Current Code and optimisation:
!       Ross Walker (TSRI, 2004)
!------------------------------

!     R2   =  INTERATOMIC DISTANCE^2 IN BOHRS2

      vec_qm_qm(1) = (XI(1)-XJ(1))*A_TO_BOHRS
      vec_qm_qm(2) = (XI(2)-XJ(2))*A_TO_BOHRS
      vec_qm_qm(3) = (XI(3)-XJ(3))*A_TO_BOHRS

!All atoms have S-orbitals
!    S-S
      do K=1,6 !1 to NGAUSS
         do L=1,6
           ADBR2=sxs_over_sas(k,l)*R2
!          CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
           IF(ADBR2 < EXPONENTIAL_CUTOFF) THEN
              SHMAT(1,1)=SHMAT(1,1)+ss_eqn(k,l)*EXP(-ADBR2)
           ENDIF
         end do
      end do
!Multiply by S-S beta factor
      SHMAT(1,1)=SHMAT(1,1)*half*betasas

      !Does atom J have p orbitals?
      if (n_atomic_orbj>1) then
         do K=1,6 !1 to NGAUSS
            do L=1,6
              ADBR2=sxp_over_sap(k,l)*R2
!             CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
              IF(ADBR2 < EXPONENTIAL_CUTOFF) THEN
                SH=sp_ovlp(k,l)*EXP(-ADBR2)
                SHMAT(1,2)=SHMAT(1,2)+SH*vec_qm_qm(1)
                SHMAT(1,3)=SHMAT(1,3)+SH*vec_qm_qm(2)
                SHMAT(1,4)=SHMAT(1,4)+SH*vec_qm_qm(3)
              ENDIF
            end do
         end do
!Multiply by S-P beta factor
         BIJSP=half*betasap
         SHMAT(1,2)=SHMAT(1,2)*BIJSP
         SHMAT(1,3)=SHMAT(1,3)*BIJSP
         SHMAT(1,4)=SHMAT(1,4)*BIJSP
      end if

      !Does atom I have p orbitals?
      if (n_atomic_orbi>1) then
         do K=1,6 !1 to NGAUSS
            do L=1,6
              ADBR2=pxs_over_pas(l,k)*R2
!             CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
              IF(ADBR2 < EXPONENTIAL_CUTOFF) THEN
                SH=-ps_ovlp(l,k)*EXP(-ADBR2)
                SHMAT(2,1)=SHMAT(2,1)+SH*vec_qm_qm(1)
                SHMAT(3,1)=SHMAT(3,1)+SH*vec_qm_qm(2)
                SHMAT(4,1)=SHMAT(4,1)+SH*vec_qm_qm(3)
              ENDIF
            end do
         end do
!Multiply by P-S beta factor
         BIJPS=half*betapas
         SHMAT(2,1)=SHMAT(2,1)*BIJPS
         SHMAT(3,1)=SHMAT(3,1)*BIJPS
         SHMAT(4,1)=SHMAT(4,1)*BIJPS
      end if

      !Now do both have p orbitals
      if (n_atomic_orbi>1 .AND. n_atomic_orbj>1) then
        BIJPP=half*betapap
        do I=1,n_atomic_orbi-1
           ii=i+1
           do J=1,n_atomic_orbj-1
           jj=j+1
!    P-P
            TOMB=vec_qm_qm(i)*vec_qm_qm(j)
            do K=1,6 !1 to NGAUSS
               do L=1,6
                 ADBR2=pxp_over_pap(k,l)*R2
!                CHECK OF OVERLAP IS NON-ZERO BEFORE DOING THE EXPONENTIAL
                 IF(ADBR2 < EXPONENTIAL_CUTOFF) THEN
                   IF(ii.EQ.jj) then
                     SH=EXP(-ADBR2)*(pp_ovlp_ieqj1(k,l)*TOMB+pp_ovlp_ieqj2(k,l))
                   else
                     SH=EXP(-ADBR2)*TOMB*pp_ovlp_inj(k,l)
                   end if
                   SHMAT(ii,jj)=SHMAT(ii,jj)+SH
                 ENDIF
               end do !I=1,6
            end do !K=1,6
!Multiply by P-P beta factor
            SHMAT(ii,jj)=SHMAT(ii,jj)*BIJPP
          end do !i=1,n_atomic_orbj-1
        end do !j=1,n_atomic_orbi-1
      end if

      RETURN
end subroutine qm2_h1elec


