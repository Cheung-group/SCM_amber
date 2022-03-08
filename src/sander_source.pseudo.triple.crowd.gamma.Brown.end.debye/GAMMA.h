#include "dprec.h"
!!11.15.2010
!!Antonios added 11.15.10 for size dependence friction

_REAL_ GAMMA_BD,RRADIUS,gamma_array,gammai_as,sdfac_as,rsd_as, &
      c_implic_as,c_explic_as,fln_as
integer nrsize,gamma_index,rsd_index,natom_as


COMMON/gamaparms/GAMMA_BD,RRADIUS,nrsize,natom_as
COMMON/xfriction/gamma_array(10000),gammai_as(10000),sdfac_as(10000), &
      rsd_as(10000), c_implic_as(10000),c_explic_as(10000),fln_as(10000)
