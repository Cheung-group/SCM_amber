subroutine nose_hoover_init_LES ( mass, v )
! initialize Nose-Hoover chain thermostats
  use pimd_vars, only: nbead,NMPIMD,CMD,ipimd
  use constants, only: pi, hbar, kB
  use full_pimd_vars, only : mybeadid
  use nose_hoover_module, only: nose_hoover_module_init,  &
                                Thermostat_init,  &
                                Thermostat_link,  &
                                Thermostat_type
  use nose_hoover_vars, only: file_nhc, nchain, thermo, nthermo, tau

  use cmd_vars, only: restart_cmd, adiab_param,  &
                      omega_nmode
  
  implicit none
# include "les.h"
# include "memory.h"
# include "md.h"
# include "parallel.h"
  _REAL_, intent(in) :: mass( natom )
  _REAL_ :: v( 3, natom )
  _REAL_ :: kT, beta
  _REAL_ :: sigma, dummy, tmp
  integer :: iseed, ibead,j
  integer :: idim, iatom
  type(Thermostat_type) :: thermo_dummy

  if ( worldrank.eq.0 ) then
     open( file_nhc, file = "NHC.dat" )
  endif

  kT = kB * temp0
  beta = 1.d0 / kT

  call nose_hoover_module_init
  nthermo = 3 * natom

  allocate(thermo(3,natom))

  if (ipimd==NMPIMD.or.ipimd==CMD) then
     if (irest==1) then
        call trans_vel_from_cart_to_nmode(v)
     endif
  endif

  if (ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
     do iatom = 1, natom
        if (cnum(iatom).eq.0) then
           tau = 1.d0/omega_nmode(1)
        else
           tau = 1.d0/omega_nmode(cnum(iatom))
        endif
        do idim  = 1, 3
           call Thermostat_init(nchain,thermo(idim,iatom),1,0,kT,tau)
           call Thermostat_link(thermo(idim,iatom),mass(iatom),v(idim,iatom))
        enddo
     enddo
  else 
     do iatom = 1, natom
        if (cnum(iatom).eq.0) then
           tau = hbar * beta
        else
           tau = hbar*beta/sqrt(dble(nbead))
        endif
        do idim  = 1, 3
           call Thermostat_init(nchain,thermo(idim,iatom),1,0,kT,tau)
           call Thermostat_link(thermo(idim,iatom),mass(iatom),v(idim,iatom))
        enddo
     enddo
  endif

end subroutine nose_hoover_init_LES


