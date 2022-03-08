#include "dprec.h"

!------------------------------------------------------------------------------
module pimd_vars

! internal variables used in PIMD
!------------------------------------------------------------------------------
  implicit none
  save

  integer :: PRIPIMD = 1   ! Primitive PIMD
  integer :: NMPIMD = 2    ! Normal-mode PIMD
  integer :: CMD = 3       ! Centroid MD
  integer :: RPMD = 4      ! Ring-polymer MD
  integer :: ipimd
  integer :: ineb
  integer :: nbead
  integer :: natomCL
  integer :: nthermo
  ! itimass = cntrl option for thermodynamic integration w.r.t. mass.
  ! itimass = 0: no TI(default), 1: virial estimator, 2: thermodynamic estimator
  integer :: itimass
  _REAL_  :: equal_part
  _REAL_  :: Epot_external, Epot_spring, Epot_deriv,Eimp_virial
  _REAL_  :: nebrms
  _REAL_  :: nbead_inv
  _REAL_  :: bnd_vir(3,3)
  _REAL_  :: tau_vol=20.4550d0
  ! dmdlm = d mass / d lambda / mass (for TI w.r.t. mass)
  _REAL_, allocatable :: dmdlm(:)
  _REAL_, allocatable :: x_centroid(:,:), real_mass(:), nrg_all(:) &
                       , tempbuf(:,:), cartpos(:,:), cartvel(:,:)
end module

module full_pimd_vars

   integer            :: mybeadid
   _REAL_             :: totener(60)
   _REAL_             :: totenert(60), totenert2(60)
   _REAL_,allocatable :: xall(:,:,:) 
end module 


module part_pimd_vars
   
  _REAL_, allocatable :: massCL( : )
  _REAL_, allocatable :: frcx_copy(:,:), frcx_temp(:,:)
  _REAL_, allocatable :: xrep(:), frep(:), vrep(:), ftmp(:,:)
  _REAL_, allocatable :: pimd_mmchg(:)

end module 

module neb_vars
  _REAL_ :: neb_rotat(3,3)
  _REAL_, allocatable :: tangents(:), springforce(:)
  integer, allocatable :: fitgroup(:), rmsgroup(:)   
end module  
