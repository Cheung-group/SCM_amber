#include "dprec.h"
#include "assert.h"

module amoeba_mdin
  implicit none
  private
  integer,save :: iamoeba=0,do_amoeba_valence=1,do_amoeba_nonbond=1, &
                  do_vdw_taper=1,do_vdw_longrange=1,am_nbead=1
  _REAL_,save :: sor_coefficient = 0.75d0
  _REAL_,save :: dipole_scf_tol = 0.01d0
  integer,save :: dipole_scf_iter_max = 50
#ifndef DISABLE_AMOEBA_CG
  integer,save :: dipole_scf_use_cg = 0
  integer,save :: dipole_scf_cg_niter = 4
#endif /* DISABLE_AMOEBA_CG */
  _REAL_,save :: ee_dsum_cut=7.d0
  _REAL_,save :: ee_damped_cut=4.5d0
  _REAL_,save :: vdw_taper = 0.9d0
  _REAL_,save :: thole_expon_coeff=0.39d0
  _REAL_,save :: compress = 0.000046d0
  integer, save :: amoeba_verbose = 0
  integer,save :: beeman_integrator = 0
  public AMOEBA_read_mdin,iamoeba,do_amoeba_valence, &
         do_amoeba_nonbond,amoeba_verbose,beeman_integrator, &
         sor_coefficient,dipole_scf_tol,dipole_scf_iter_max, &
#ifndef DISABLE_AMOEBA_CG
         dipole_scf_use_cg, dipole_scf_cg_niter, &
#endif /* DISABLE_AMOEBA_CG */
         ee_dsum_cut,ee_damped_cut,thole_expon_coeff,vdw_taper, &
         compress,do_vdw_taper,do_vdw_longrange,am_nbead
  contains
!-------------------------------------------------------------------------------
subroutine AMOEBA_read_mdin(nf)
  implicit none
  integer,intent(in) :: nf

  integer            :: do_bond=1,do_ureyb=1,do_reg_angle=1,  &
                        do_trig_angle=1,do_opbend=1,do_torsion=1, &
                        do_pi_torsion=1,do_strbend=1,do_torsion_torsion=1, &
                        do_str_torsion=1, &
                        do_recip=1,do_adjust=1,do_direct=1,do_self=1, &
                        do_vdw=1,do_induced=1
  integer ifind
  namelist/amoeba/do_amoeba_valence,do_amoeba_nonbond, &
                  do_bond,do_ureyb,do_reg_angle,  &
                  do_trig_angle,do_opbend,do_torsion,do_str_torsion, &
                  do_pi_torsion,do_strbend,do_torsion_torsion, &
                  do_recip,do_adjust,do_direct,do_self, &
                  do_vdw,do_induced,amoeba_verbose,beeman_integrator, & 
                  sor_coefficient,dipole_scf_tol, &
                  dipole_scf_iter_max,&
#ifndef DISABLE_AMOEBA_CG
                  dipole_scf_use_cg, dipole_scf_cg_niter, &
#endif /* DISABLE_AMOEBA_CG */
                  ee_dsum_cut,ee_damped_cut, &
                  thole_expon_coeff,vdw_taper,do_vdw_taper,do_vdw_longrange, &
                  compress,am_nbead

# include "files.h"
  
  read(nf,nml=amoeba)
  call AM_VAL_set_user_bit(do_bond,do_ureyb,do_reg_angle,do_trig_angle, &
                          do_opbend,do_torsion,do_str_torsion, &
                          do_pi_torsion,do_strbend, &
                          do_torsion_torsion)
  call AM_NONBOND_set_user_bit(do_recip,do_adjust,do_direct,do_self, &
                               do_vdw,do_induced)
  
end subroutine AMOEBA_read_mdin
!----------------------------------------------------------
end module amoeba_mdin
!-------------------------------------------------------------------------------
